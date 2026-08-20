package service

import (
	"context"
	"encoding/json"
	"fmt"

	"task136-seqannot/internal/store"
)

// ReconcileAll is the restart-recovery entry point. It walks every job that is
// not in a terminal state (done/failed) and resumes it: completed steps are
// skipped (already committed atomically) and the remaining steps are
// re-executed. After ReconcileAll returns, every previously-pending job is
// either done or failed with a recorded error.
func (svc *Service) ReconcileAll(ctx context.Context) ([]store.JobRow, error) {
	pending, err := svc.store.ListPendingJobs(ctx)
	if err != nil {
		return nil, err
	}
	resumed := make([]store.JobRow, 0, len(pending))
	for _, j := range pending {
		// Reset the error field; running will overwrite on completion.
		j.Error = ""
		_ = svc.store.SaveJob(ctx, j)
		if err := svc.runJob(ctx, j.ID); err != nil {
			// failure is recorded inside runJob; keep going on other jobs
			_ = err
		}
		final, err := svc.store.GetJob(ctx, j.ID)
		if err == nil {
			resumed = append(resumed, final)
		}
	}
	return resumed, nil
}

// ReplayJob re-derives every step's result from the stored sequence + params
// and compares it against the persisted results. It returns equal=true when
// every step's recomputed hash matches the stored hash, which is the
// idempotency guarantee: a completed job replayed after a restart must produce
// byte-identical results.
type ReplayReport struct {
	JobID  string         `json:"job_id"`
	Equal  bool           `json:"equal"`
	Steps  []ReplayStep   `json:"steps"`
}

type ReplayStep struct {
	Name      string `json:"name"`
	Stored    string `json:"stored_hash"`
	Recomputed string `json:"recomputed_hash"`
	Equal     bool   `json:"equal"`
}

func (svc *Service) ReplayJob(ctx context.Context, jobID string) (ReplayReport, error) {
	job, err := svc.store.GetJob(ctx, jobID)
	if err != nil {
		return ReplayReport{}, err
	}
	var steps []StepSpec
	if err := json.Unmarshal([]byte(job.Params), &steps); err != nil {
		return ReplayReport{}, fmt.Errorf("service: bad job params: %w", err)
	}
	storedSteps, err := svc.store.ListSteps(ctx, jobID)
	if err != nil {
		return ReplayReport{}, err
	}
	storedHash := map[string]string{}
	for _, s := range storedSteps {
		storedHash[s.StepName] = s.ResultHash
	}
	report := ReplayReport{JobID: jobID, Equal: true}
	for i, sp := range steps {
		key := stepKey(i, sp)
		res, err := svc.executeStep(ctx, job.SequenceID, sp)
		if err != nil {
			report.Equal = false
			report.Steps = append(report.Steps, ReplayStep{
				Name: key, Stored: storedHash[key], Recomputed: "",
				Equal: false,
			})
			continue
		}
		h := resultHash(res)
		eq := h == storedHash[key]
		if !eq {
			report.Equal = false
		}
		report.Steps = append(report.Steps, ReplayStep{
			Name: key, Stored: storedHash[key], Recomputed: h, Equal: eq,
		})
	}
	return report, nil
}

// ResumeJob is the per-job resume entry point exposed via the API. It is the
// same as the reconcile path but scoped to one job.
func (svc *Service) ResumeJob(ctx context.Context, jobID string) (store.JobRow, error) {
	j, err := svc.store.GetJob(ctx, jobID)
	if err != nil {
		return store.JobRow{}, err
	}
	if j.Status == StatusDone {
		return j, nil
	}
	if err := svc.runJob(ctx, jobID); err != nil {
		return store.JobRow{}, err
	}
	return svc.store.GetJob(ctx, jobID)
}
