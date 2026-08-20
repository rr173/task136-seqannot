package seq

import "testing"

func TestBasesMatch(t *testing.T) {
	cases := []struct {
		a, b byte
		want  bool
	}{
		{'A', 'A', true}, {'A', 'T', false}, {'A', 'R', true}, {'G', 'R', true}, {'C', 'R', false},
		{'T', 'Y', true}, {'C', 'Y', true}, {'G', 'B', true}, {'A', 'B', false}, {'N', 'A', true},
		{'a', 'r', true}, {'a', 'g', false},
	}
	for _, c := range cases {
		if got := BasesMatch(c.a, c.b); got != c.want {
			t.Errorf("BasesMatch(%q,%q)=%v want %v", c.a, c.b, got, c.want)
		}
	}
}

func TestReverseComplement(t *testing.T) {
	got := ReverseComplement("ATGC")
	want := "GCAT"
	if got != want {
		t.Fatalf("ReverseComplement(ATGC)=%q want %q", got, want)
	}
	// palindrome complements to itself
	if got := ReverseComplement("GAATTC"); got != "GAATTC" {
		t.Fatalf("palindrome RC=%q", got)
	}
	// IUPAC degenerate: R<->Y
	if got := ReverseComplement("AGR"); got != "YCT" {
		t.Fatalf("RC(AGR)=%q want YCT", got)
	}
}

func TestIsValidSequence(t *testing.T) {
	if !IsValidSequence("ACGTN") {
		t.Fatal("ACGTN should be valid")
	}
	if IsValidSequence("ACGTX") {
		t.Fatal("X should be invalid")
	}
	if IsValidSequence("") {
		t.Fatal("empty should be invalid")
	}
}

func TestComposeGC(t *testing.T) {
	s := Sequence{Residues: "ATGC", Type: TypeLinear}
	c := Compose(s)
	// GC fraction = (G+C)/4 = 2/4 = 0.5 -> 5000bp
	if c.GCBp != 5000 {
		t.Fatalf("GCBp=%d want 5000", c.GCBp)
	}
	// counts: A1 C1 G1 T1
	if c.Counts['G'] != 1 || c.Counts['C'] != 1 {
		t.Fatalf("counts=%v", c.Counts)
	}
}

func TestComposeMWCircular(t *testing.T) {
	lin := Compose(Sequence{Residues: "ATGC", Type: TypeLinear})
	cir := Compose(Sequence{Residues: "ATGC", Type: TypeCircular})
	// linear terminal -61960, circular +79000 -> difference 140960
	if cir.MWBp-lin.MWBp != 140960 {
		t.Fatalf("MW delta=%d want 140960", cir.MWBp-lin.MWBp)
	}
}

func TestCodonAmino(t *testing.T) {
	a, ok := CodonAmino('A', 'T', 'G')
	if !ok || a != Metionine {
		t.Fatalf("ATG=%q ok=%v want M", a, ok)
	}
	a, ok = CodonAmino('T', 'A', 'A')
	if !ok || a != Stop {
		t.Fatalf("TAA=%q ok=%v want *", a, ok)
	}
	// ambiguous: NNN could be anything
	_, unamb := CodonAmino('N', 'N', 'N')
	if unamb {
		t.Fatal("NNN should be ambiguous")
	}
}

func TestTranslate6Frames(t *testing.T) {
	// A short sequence with a start and stop in frame 0:
	// ATG AAA TAA  -> M K *
	s := Sequence{Residues: "ATGAAATAA", Type: TypeLinear}
	frames := Translate6(s)
	if len(frames) != 6 {
		t.Fatalf("want 6 frames got %d", len(frames))
	}
	f0, _ := FrameFor(frames, 0)
	if f0.Aminos != "MK*" {
		t.Fatalf("frame0 aminos=%q want MK*", f0.Aminos)
	}
	// frame 0 first codon is start, last is stop
	if !f0.Codons[0].IsStart || !f0.Codons[2].IsStop {
		t.Fatalf("codon flags wrong: %+v", f0.Codons)
	}
}

func TestParseFasta(t *testing.T) {
	in := ">seq1 a test\nATGC\nGGGG\n>seq2\nTTTT\n"
	recs, err := ParseFasta(in)
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 2 {
		t.Fatalf("want 2 records got %d", len(recs))
	}
	if recs[0].ID != "seq1" || recs[0].Desc != "a test" || recs[0].Residues != "ATGCGGGG" {
		t.Fatalf("rec0=%+v", recs[0])
	}
	if recs[1].Residues != "TTTT" {
		t.Fatalf("rec1=%+v", recs[1])
	}
}

func TestParseFastaBareSequence(t *testing.T) {
	recs, err := ParseFasta("ATGCATGC")
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 1 || recs[0].Residues != "ATGCATGC" {
		t.Fatalf("bare seq parse=%+v", recs)
	}
}

func TestWrapIndex(t *testing.T) {
	cases := []struct {
		pos, n, want int
	}{
		{5, 10, 5}, {10, 10, 0}, {12, 10, 2}, {-1, 10, 9}, {0, 10, 0},
	}
	for _, c := range cases {
		if got := WrapIndex(c.pos, c.n); got != c.want {
			t.Errorf("WrapIndex(%d,%d)=%d want %d", c.pos, c.n, got, c.want)
		}
	}
}

func TestFormatFastaRoundTrip(t *testing.T) {
	in := []FastaRecord{{Header: "x", ID: "x", Residues: "ATGCATGC"}}
	out := FormatFasta(in, 4)
	parsed, err := ParseFasta(out)
	if err != nil {
		t.Fatal(err)
	}
	if parsed[0].Residues != "ATGCATGC" {
		t.Fatalf("round trip failed: %q", parsed[0].Residues)
	}
}
