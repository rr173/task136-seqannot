package align

import "testing"

func TestGlobalIdentical(t *testing.T) {
	r := Global("ACGT", "ACGT", Params{Match: 1, Mismatch: -1, Gap: -1})
	if r.Score != 4 {
		t.Fatalf("score=%d want 4", r.Score)
	}
	if r.AlignedA != "ACGT" || r.AlignedB != "ACGT" {
		t.Fatalf("aligned=%q/%q", r.AlignedA, r.AlignedB)
	}
}

func TestGlobalWithGap(t *testing.T) {
	// ACGT vs ACGGT: one insertion, optimal global alignment
	r := Global("ACGT", "ACGGT", Params{Match: 1, Mismatch: -1, Gap: -1})
	// score should be 4 matches - 1 gap = 3
	if r.Score != 3 {
		t.Fatalf("score=%d want 3", r.Score)
	}
}

func TestLocalFindsBestRegion(t *testing.T) {
	// xxxxACGTxxxx vs yyyyyACGTyyyy -> local best is ACGT/ACGT score 4
	r := Local("AAAAACGTAAAA", "TTTTACGTTTTT", Params{Match: 2, Mismatch: -1, Gap: -2})
	if r.Score != 8 {
		t.Fatalf("local score=%d want 8", r.Score)
	}
	if r.AlignedA != "ACGT" || r.AlignedB != "ACGT" {
		t.Fatalf("local aligned=%q/%q", r.AlignedA, r.AlignedB)
	}
}

func TestGlobalMismatch(t *testing.T) {
	r := Global("AAAA", "TTTT", Params{Match: 1, Mismatch: -1, Gap: -1})
	// all mismatches: score -4
	if r.Score != -4 {
		t.Fatalf("score=%d want -4", r.Score)
	}
}

func TestDeterminism(t *testing.T) {
	p := Params{Match: 2, Mismatch: -1, Gap: -1}
	r1 := Global("GATTACA", "GCATGCU", p)
	r2 := Global("GATTACA", "GCATGCU", p)
	if r1 != r2 {
		t.Fatalf("non-deterministic: %+v != %+v", r1, r2)
	}
}
