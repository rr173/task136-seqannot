package service

import (
	"context"
	"errors"
	"path/filepath"
	"testing"

	"task136-seqannot/internal/seq"
	"task136-seqannot/internal/store"
)

func TestBug10_CancelledStepPersistsFailure(t *testing.T) {
	st, err := store.Open(filepath.Join(t.TempDir(), "test.db"))
	if err != nil { t.Fatal(err) }
	defer st.Close()
	ResetCounters()
	svc := New(st)
	ctx := context.Background()
	seqRow, err := svc.CreateSequence(ctx, SequenceRequest{Name: "s", Residues: "ATGC", Type: string(seq.TypeLinear)})
	if err != nil { t.Fatal(err) }
	job := store.JobRow{ID: "job-z", SequenceID: seqRow.ID, Params: `[{"name":"composition"}]`, Status: StatusPending, Results: "{}", CreatedAt: "1", UpdatedAt: "1"}
	if err := st.SaveJob(ctx, job); err != nil { t.Fatal(err) }
	cancelled, cancel := context.WithCancel(ctx)
	cancel()
	if err := svc.runJob(cancelled, job.ID); !errors.Is(err, context.Canceled) {
		t.Fatalf("runJob error=%v want context canceled", err)
	}
	got, err := st.GetJob(ctx, job.ID)
	if err != nil { t.Fatal(err) }
	if got.Status != StatusFailed || got.Error == "" {
		t.Fatalf("cancelled job was not durably failed: %+v", got)
	}
}
