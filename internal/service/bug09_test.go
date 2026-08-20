package service

import (
	"context"
	"path/filepath"
	"testing"

	"task136-seqannot/internal/seq"
	"task136-seqannot/internal/store"
)

func TestBug09_ResumeProgressCountsCommittedSteps(t *testing.T) {
	st, err := store.Open(filepath.Join(t.TempDir(), "test.db"))
	if err != nil { t.Fatal(err) }
	defer st.Close()
	ResetCounters()
	svc := New(st)
	ctx := context.Background()
	seqRow, err := svc.CreateSequence(ctx, SequenceRequest{Name: "s", Residues: "ATGC", Type: string(seq.TypeLinear)})
	if err != nil { t.Fatal(err) }
	job := store.JobRow{ID: "job-y", SequenceID: seqRow.ID, Params: `[{"name":"composition"},{"name":"reverse_complement"},{"name":"unknown"}]`, Status: StatusRunning, Progress: 0, Results: "{}", CreatedAt: "1", UpdatedAt: "1"}
	if err := st.SaveJob(ctx, job); err != nil { t.Fatal(err) }
	if err := st.InTx(ctx, func(tx store.DBTX) error {
		return store.UpsertStep(tx, ctx, store.StepRow{JobID: job.ID, StepName: "0:composition", Status: StatusDone, ResultHash: "x", UpdatedAt: "1"})
	}); err != nil { t.Fatal(err) }
	if err := st.InTx(ctx, func(tx store.DBTX) error {
		return store.UpsertStep(tx, ctx, store.StepRow{JobID: job.ID, StepName: "1:reverse_complement", Status: StatusDone, ResultHash: "x", UpdatedAt: "1"})
	}); err != nil { t.Fatal(err) }
	if err := st.InTx(ctx, func(tx store.DBTX) error {
		_, err := tx.ExecContext(ctx, `UPDATE analysis_jobs SET results=? WHERE id=?`, `{"0:composition":{}}`, job.ID)
		return err
	}); err != nil { t.Fatal(err) }
	if _, err := svc.ResumeJob(ctx, job.ID); err == nil { t.Fatal("unknown step unexpectedly succeeded") }
	got, err := st.GetJob(ctx, job.ID)
	if err != nil { t.Fatal(err) }
	if got.Progress != 2 || got.Status != StatusFailed {
		t.Fatalf("resumed progress=%d status=%s", got.Progress, got.Status)
	}
}
