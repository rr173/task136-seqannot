package annot

import "testing"

func TestBug04_CircularQueryKeepsWrappedInterval(t *testing.T) {
	features := []Feature{
		{ID: "wrap", SequenceID: "s", Type: "gene", Start: 9, End: 2},
		{ID: "middle", SequenceID: "s", Type: "gene", Start: 4, End: 6},
	}
	got := QueryOverlap(features, 10, true, 8, 2)
	if len(got) != 1 || got[0].ID != "wrap" {
		t.Fatalf("wrapped query returned %+v", got)
	}
}
