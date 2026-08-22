package service

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"

	"task136-seqannot/internal/store"
)

// JobSpec describes a submitted analysis job. Each Step is one analysis item
// to run over the target sequence; the Steps execute in order and each step's
// result is committed atomically with its progress record.
type JobSpec struct {
	SequenceID string    `json:"sequence_id"`
	Steps      []StepSpec `json:"steps"`
}

// StepSpec is one analysis item. Name is the analysis kind; Params carries
// the item-specific arguments (e.g. min_aa_len for ORF, pattern for motif).
type StepSpec struct {
	Name   string          `json:"name"`
	Params json.RawMessage `json:"params"`
}

// JobStatus values for the persisted analysis_jobs.status column.
const (
	StatusPending = "pending"
	StatusRunning  = "running"
	StatusDone     = "done"
	StatusFailed   = "failed"
)

// SubmitJob validates and persists a job in pending state, then runs it
// synchronously to completion. The synchronous execution is fine because the
// store serializes all access on a single connection; the persisted step
// records make the job crash-resumable regardless of whether execution is
// sync or async.
func (svc *Service) SubmitJob(ctx context.Context, spec JobSpec) (store.JobRow, error) {
	if _, err := svc.store.GetSequence(ctx, spec.SequenceID); err != nil {
		return store.JobRow{}, err
	}
	if len(spec.Steps) == 0 {
		return store.JobRow{}, fmt.Errorf("service: job has no steps")
	}
	params, _ := json.Marshal(spec.Steps)
	jobID := newID("job")
	row := store.JobRow{
		ID: jobID, SequenceID: spec.SequenceID, Params: string(params),
		Status: StatusPending, Results: "{}", CreatedAt: stamp(), UpdatedAt: stamp(),
	}
	if err := svc.store.SaveJob(ctx, row); err != nil {
		return store.JobRow{}, err
	}
	if err := svc.runJob(ctx, jobID); err != nil {
		return store.JobRow{}, err
	}
	return svc.store.GetJob(ctx, jobID)
}

// runJob executes the steps of a job, resuming from the last committed step.
// Each step's result is written in its own transaction together with the
// job_steps row and the progress bump, so a crash between steps leaves the job
// in a consistent "partially done, resumable" state.
func (svc *Service) runJob(ctx context.Context, jobID string) error {
	job, err := svc.store.GetJob(ctx, jobID)
	if err != nil {
		return err
	}
	var steps []StepSpec
	if err := json.Unmarshal([]byte(job.Params), &steps); err != nil {
		return fmt.Errorf("service: bad job params: %w", err)
	}
	// load already-completed steps so we can skip them on resume
	done, err := svc.store.ListSteps(ctx, jobID)
	if err != nil {
		return err
	}
	doneSet := map[string]bool{}
	for _, d := range done {
		if d.Status == StatusDone {
			doneSet[d.StepName] = true
		}
	}
	results := map[string]any{}
	if job.Results != "" && job.Results != "{}" {
		_ = json.Unmarshal([]byte(job.Results), &results)
	}
	// mark running. Clear any error left by a prior failed run so that, while
	// the job is executing, a caller never observes the stale failure message.
	// The success path clears it again on completion; the failure path
	// overwrites it with the new error.
	job.Status = StatusRunning
	job.Error = ""
	job.UpdatedAt = stamp()
	if err := svc.store.SaveJob(ctx, job); err != nil {
		return err
	}
	for i, sp := range steps {
		stepKey := stepKey(i, sp)
		if doneSet[stepKey] {
			continue
		}
		res, err := svc.executeStep(ctx, job.SequenceID, sp)
		if err != nil {
			// record failure and stop; already-done steps remain
			job, _ = svc.store.GetJob(ctx, jobID)
			job.Status = StatusFailed
			job.Error = err.Error()
			job.UpdatedAt = stamp()
			_ = svc.store.SaveJob(ctx, job)
			return err
		}
		hash := resultHash(res)
		// commit step result + progress atomically
		results[stepKey] = res
		resultsJSON, _ := json.Marshal(results)
		err = svc.store.InTx(ctx, func(tx store.DBTX) error {
			if err := store.UpsertStep(tx, ctx, store.StepRow{
				JobID: jobID, StepName: stepKey, Status: StatusDone,
				ResultHash: hash, UpdatedAt: stamp(),
			}); err != nil {
				return err
			}
			_, err := tx.ExecContext(ctx, `UPDATE analysis_jobs SET results=?, progress=?, updated_at=? WHERE id=?`,
				string(resultsJSON), i+1, stamp(), jobID)
			return err
		})
		if err != nil {
			return err
		}
	}
	// all steps done
	job, err = svc.store.GetJob(ctx, jobID)
	if err != nil {
		return err
	}
	// Clear any error left by a previous failed run. A job that failed, was
	// resumed, and now completed must present a clean done state: status=done
	// with an empty error, never the stale message from the prior failure.
	job.Status = StatusDone
	job.Progress = len(steps)
	job.Error = ""
	job.UpdatedAt = stamp()
	return svc.store.SaveJob(ctx, job)
}

// stepKey produces a stable key for a step within a job. Two steps of the same
// name with the same params are allowed (and produce distinct keys by index).
func stepKey(i int, sp StepSpec) string {
	return fmt.Sprintf("%d:%s", i, sp.Name)
}

// executeStep runs one analysis item and returns its result. The result is a
// pure function of the stored sequence + step params, which makes the job
// idempotent and replay-verifiable.
func (svc *Service) executeStep(ctx context.Context, seqID string, sp StepSpec) (any, error) {
	switch sp.Name {
	case "composition":
		return svc.Composition(ctx, seqID)
	case "reverse_complement":
		return map[string]any{"sequence": mustRevComp(svc, ctx, seqID)}, nil
	case "translate":
		return svc.Translate(ctx, seqID)
	case "orf":
		var p struct {
			MinAALen int `json:"min_aa_len"`
		}
		if len(sp.Params) > 0 {
			if err := json.Unmarshal(sp.Params, &p); err != nil {
				return nil, err
			}
		}
		return svc.FindORFs(ctx, seqID, p.MinAALen)
	case "motif":
		var p struct {
			Pattern string `json:"pattern"`
			MotifID string `json:"motif_id"`
		}
		if len(sp.Params) > 0 {
			if err := json.Unmarshal(sp.Params, &p); err != nil {
				return nil, err
			}
		}
		return svc.SearchMotif(ctx, seqID, p.Pattern, p.MotifID)
	case "restriction":
		var p struct {
			EnzymeIDs []string `json:"enzyme_ids"`
		}
		if len(sp.Params) > 0 {
			if err := json.Unmarshal(sp.Params, &p); err != nil {
				return nil, err
			}
		}
		return svc.RestrictionMap(ctx, seqID, p.EnzymeIDs)
	case "digest":
		var p struct {
			EnzymeIDs []string `json:"enzyme_ids"`
		}
		if len(sp.Params) > 0 {
			if err := json.Unmarshal(sp.Params, &p); err != nil {
				return nil, err
			}
		}
		return svc.Digest(ctx, seqID, p.EnzymeIDs)
	default:
		return nil, fmt.Errorf("service: unknown step %q", sp.Name)
	}
}

func mustRevComp(svc *Service, ctx context.Context, id string) string {
	r, err := svc.ReverseComplement(ctx, id)
	if err != nil {
		return ""
	}
	return r
}

// resultHash returns a stable hex digest of a step result for idempotency
// verification. The hash is computed over the canonical JSON of the result.
func resultHash(v any) string {
	b, err := json.Marshal(v)
	if err != nil {
		return ""
	}
	sum := sha256.Sum256(b)
	return hex.EncodeToString(sum[:])
}
