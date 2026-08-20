package service

import (
	"context"
	"path/filepath"
	"testing"

	"task136-seqannot/internal/seq"
	"task136-seqannot/internal/store"
)

func TestBug08_ResumeClearsPreviousError(t *testing.T) {
	st, err := store.Open(filepath.Join(t.TempDir(), "test.db"))
	if err != nil { t.Fatal(err) }
	defer st.Close()
	ResetCounters()
	svc := New(st)
	ctx := context.Background()
	seqRow, err := svc.CreateSequence(ctx, SequenceRequest{Name: "s", Residues: "ATGC", Type: string(seq.TypeLinear)})
	if err != nil { t.Fatal(err) }
	job := store.JobRow{ID: "job-x", SequenceID: seqRow.ID, Params: `[{"name":"composition"}]`, Status: StatusFailed, Error: "temporary failure", Results: "{}", CreatedAt: "1", UpdatedAt: "1"}
	if err := st.SaveJob(ctx, job); err != nil { t.Fatal(err) }
	if _, err := svc.ResumeJob(ctx, job.ID); err != nil { t.Fatal(err) }
	got, err := st.GetJob(ctx, job.ID)
	if err != nil { t.Fatal(err) }
	if got.Status != StatusDone || got.Error != "" {
		t.Fatalf("resumed job=%+v", got)
	}
}
