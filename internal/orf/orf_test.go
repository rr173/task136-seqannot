package orf

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestFindSimpleORF(t *testing.T) {
	// ATG AAA TAA -> M K *  in frame 0
	s, _ := seq.NewSequence("t", "ATGAAATAA", seq.TypeLinear, "", "")
	orfs := Find(s, Options{MinAALen: 1})
	if len(orfs) == 0 {
		t.Fatal("expected at least 1 ORF")
	}
	var f0 *ORF
	for i := range orfs {
		if orfs[i].Frame == 0 && orfs[i].Strand == '+' {
			f0 = &orfs[i]
			break
		}
	}
	if f0 == nil {
		t.Fatalf("no frame0 + ORF among %v", orfs)
	}
	if f0.Protein != "MK*" {
		t.Fatalf("protein=%q want MK*", f0.Protein)
	}
	if f0.Start != 1 || f0.End != 9 {
		t.Fatalf("coords=[%d,%d] want [1,9]", f0.Start, f0.End)
	}
	if f0.Length != 3 {
		t.Fatalf("length=%d want 3", f0.Length)
	}
}

func TestFindMinAALen(t *testing.T) {
	// ATG TAA -> M * (aaLen excluding stop = 1)
	s, _ := seq.NewSequence("t", "ATGTAA", seq.TypeLinear, "", "")
	short := Find(s, Options{MinAALen: 5})
	for _, o := range short {
		if o.StopCodon != "" && o.Start == 1 {
			t.Fatalf("min filter should drop the 1-aa ORF, got %+v", o)
		}
	}
	long := Find(s, Options{MinAALen: 1})
	found := false
	for _, o := range long {
		if o.StopCodon != "" && o.Start == 1 {
			found = true
		}
	}
	if !found {
		t.Fatal("min=1 should keep the ORF")
	}
}

func TestFindOverlappingORFsAllReported(t *testing.T) {
	// Two forward frames with overlapping ORFs:
	// frame0: ATG AAA TAA  (1..9)
	// frame1: TGA ... is a stop in frame1 -> need careful construction.
	// Use a sequence with ORFs in frame 0 and frame 1 that overlap.
	// "ATG ATG AAA TAA TGA" — frame0 ATG@1, second ATG@4, stop TAA@13 (M M K *)
	// frame1 starts at offset 1: TGA@2 is a stop immediately.
	s, _ := seq.NewSequence("t", "ATGATGAAATAATGA", seq.TypeLinear, "", "")
	orfs := Find(s, Options{MinAALen: 1})
	// We must have at least the frame0 ORF M MK* and not lose it.
	hasFrame0 := false
	for _, o := range orfs {
		if o.Frame == 0 && o.Protein == "MMK*" {
			hasFrame0 = true
		}
	}
	if !hasFrame0 {
		t.Fatalf("missing frame0 MMK* among %+v", orfs)
	}
}

func TestFindReverseStrandORF(t *testing.T) {
	// Build a sequence whose reverse complement contains a clean ORF.
	// We want rc = "ATG CCC TAA". rc of that is "TTAGGGCAT".
	// So original = "TTAGGGCAT"; reverse-complement back is ATGCCCTAA -> M P *
	s, _ := seq.NewSequence("t", "TTAGGGCAT", seq.TypeLinear, "", "")
	orfs := Find(s, Options{MinAALen: 1})
	hasMP := false
	for _, o := range orfs {
		if o.Protein == "MP*" {
			hasMP = true
		}
	}
	if !hasMP {
		t.Fatalf("reverse ORF MP* not found among %+v", orfs)
	}
}

func TestFindCircularWrap(t *testing.T) {
	// A circular sequence where the start codon is the last codon of the frame
	// (it spans the origin) and the only in-frame stop sits across the origin,
	// before the start in linear codon order. Forward frame 0 of "TAAATGAT"
	// (n=8, circular): codons ATG (positions 3-5, start) then wrap-completed
	// TTA->ATT (positions 6-8,1,2) ... actually frame 0 offset 0 yields codons
	// ATG(1-3)? No: we use a known-good case. GTAACCCAT forward frame 1 emits
	// ATG(8,9,1) start and closes on TAA(2-4) stop across the origin.
	s, _ := seq.NewSequence("t", "GTAACCCAT", seq.TypeCircular, "", "")
	orfs := Find(s, Options{MinAALen: 1})
	var got *ORF
	for i := range orfs {
		if orfs[i].Frame == 1 {
			got = &orfs[i]
			break
		}
	}
	if got == nil {
		t.Fatalf("frame 1 ORF missing among %+v", orfs)
	}
	// start at the origin-spanning ATG (position 8), stop across the origin
	// (TAA at 2-4): End wraps below Start, protein is M*.
	if got.Start != 8 || got.End != 4 {
		t.Fatalf("coords=[%d,%d] want [8,4] (origin-spanning)", got.Start, got.End)
	}
	if got.StartCodon != "ATG" || got.StopCodon != "TAA" {
		t.Fatalf("codons start=%q stop=%q want ATG/TAA", got.StartCodon, got.StopCodon)
	}
	if got.Protein != "M*" {
		t.Fatalf("protein=%q want M*", got.Protein)
	}
	if got.Length != 2 {
		t.Fatalf("length=%d want 2 (M + stop)", got.Length)
	}
}

func TestFindNoStopPartial(t *testing.T) {
	// ATG AAA AAA ... no stop -> partial ORF
	s, _ := seq.NewSequence("t", "ATGAAAAAA", seq.TypeLinear, "", "")
	orfs := Find(s, Options{MinAALen: 1})
	hasPartial := false
	for _, o := range orfs {
		if o.StopCodon == "" && o.StartCodon == "ATG" {
			hasPartial = true
		}
	}
	if !hasPartial {
		t.Fatalf("expected a partial ORF among %+v", orfs)
	}
}
