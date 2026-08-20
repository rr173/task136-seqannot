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
