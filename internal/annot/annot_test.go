package annot

import "testing"

func TestValidateLinear(t *testing.T) {
	f := Feature{Start: 3, End: 5, Strand: Plus}
	if err := f.Validate(10, false); err != nil {
		t.Fatal(err)
	}
	bad := Feature{Start: 5, End: 3, Strand: Plus}
	if err := bad.Validate(10, false); err == nil {
		t.Fatal("expected start>end error on linear")
	}
	if err := bad.Validate(10, true); err != nil {
		t.Fatalf("circular should allow wrap: %v", err)
	}
}

func TestSpansWrap(t *testing.T) {
	f := Feature{Start: 8, End: 2, Strand: Plus} // wraps on n=10
	if !f.Spans(9, 10) || !f.Spans(10, 10) || !f.Spans(1, 10) || !f.Spans(2, 10) {
		t.Fatal("wrap spans wrong")
	}
	if f.Spans(5, 10) {
		t.Fatal("5 should not be covered by wrap [8,2]")
	}
}

func TestQueryOverlap(t *testing.T) {
	feats := []Feature{
		{SequenceID: "s", Start: 1, End: 3, Strand: Plus, Label: "A"},
		{SequenceID: "s", Start: 5, End: 7, Strand: Plus, Label: "B"},
	}
	out := QueryOverlap(feats, 10, false, 3, 5)
	if len(out) != 2 {
		t.Fatalf("want 2 overlapping got %d", len(out))
	}
}

func TestQueryOverlapCircular(t *testing.T) {
	feats := []Feature{
		{SequenceID: "s", Start: 9, End: 2, Strand: Plus, Label: "wrap"},
	}
	out := QueryOverlap(feats, 10, true, 1, 2)
	if len(out) != 1 {
		t.Fatalf("circular overlap missed: %+v", out)
	}
}

// TestQueryOverlapCircularWrapQuery covers a wrap-around query on a circular
// molecule: querying [9,2] must return only features touching the arc
// [9,10] U [1,2] and must not pick up features sitting in the uncovered
// middle (positions 3..8). A linear query is normalized, so the same
// [9,2] on a linear molecule should match the [2,9] span instead.
func TestQueryOverlapCircularWrapQuery(t *testing.T) {
	feats := []Feature{
		{SequenceID: "s", Start: 1, End: 2, Strand: Plus, Label: "head"},   // in wrap arc
		{SequenceID: "s", Start: 9, End: 10, Strand: Plus, Label: "tail"}, // in wrap arc
		{SequenceID: "s", Start: 4, End: 6, Strand: Plus, Label: "mid"},   // uncovered middle
		{SequenceID: "s", Start: 8, End: 4, Strand: Plus, Label: "fwrap"}, // wrapping feature, touches tail+mid+head
	}

	out := QueryOverlap(feats, 10, true, 9, 2)
	got := make(map[string]bool, len(out))
	for _, f := range out {
		got[f.Label] = true
	}
	// Arc [9,10] U [1,2]: head and tail touch it; mid (positions 4..6) does
	// not. fwrap ([8,10] U [1,4]) touches the arc at the tail and head, so
	// it must be included.
	for _, want := range []string{"head", "tail", "fwrap"} {
		if !got[want] {
			t.Errorf("circular wrap query: missing %q in %+v", want, got)
		}
	}
	if got["mid"] {
		t.Errorf("circular wrap query: mid must not be matched (uncovered middle), got %+v", got)
	}
}

func TestQueryOverlapLinearNormalizesReversed(t *testing.T) {
	feats := []Feature{
		{SequenceID: "s", Start: 4, End: 6, Strand: Plus, Label: "mid"},
	}
	// On a linear molecule a reversed query [9,2] is normalized to [2,9],
	// which does cover position 4..6.
	out := QueryOverlap(feats, 10, false, 9, 2)
	if len(out) != 1 {
		t.Fatalf("linear reversed query should normalize and match: %+v", out)
	}
}
