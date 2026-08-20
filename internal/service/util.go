package service

import (
	"strconv"

	"task136-seqannot/internal/annot"
	"task136-seqannot/internal/store"
)

// newID generates a deterministic-enough unique id for an entity of the given
// kind. The engine does not require global uniqueness across processes (the
// SQLite primary key is the source of truth), but within a run ids must not
// collide so that selfcheck scenarios (which reset the store between cases) can
// assert deterministic row counts. We therefore use a monotonic counter per
// kind, seeded to 0 at process start; the counter is process-local and the
// id is only ever used as a primary key, never as an idempotency key for
// business results.
var counters = map[string]int{}

// ResetCounters zeroes every id counter. The selfcheck harness calls it before
// each scenario so persisted ids are deterministic from 1. It is exported only
// for the selfcheck package; production code never needs it.
func ResetCounters() {
	for k := range counters {
		delete(counters, k)
	}
}

func newID(kind string) string {
	counters[kind]++
	return kind + "-" + padLeft(strconv.Itoa(counters[kind]), 6)
}

func padLeft(s string, width int) string {
	for len(s) < width {
		s = "0" + s
	}
	return s
}

// toFeature converts a stored feature row into the annot.Feature domain type.
func toFeature(r store.FeatureRow) annot.Feature {
	strand := annot.Plus
	if r.Strand == "-" {
		strand = annot.Minus
	}
	return annot.Feature{
		ID: r.ID, SequenceID: r.SequenceID, Type: annot.FeatureType(r.Type),
		Start: r.Start, End: r.End, Strand: strand, Label: r.Label,
	}
}
