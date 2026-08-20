package motif

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestBug02_ReverseStrandMotifUsesReverseWindow(t *testing.T) {
	s, err := seq.NewSequence("reverse", "AACG", seq.TypeLinear, "", "")
	if err != nil {
		t.Fatal(err)
	}
	m, err := Compile("m1", "target", "CGT", "")
	if err != nil {
		t.Fatal(err)
	}
	hits := Search(s, m)
	for _, h := range hits {
		if h.Strand == '-' && h.Start == 2 && h.End == 4 && h.Matched == "CGT" {
			return
		}
	}
	t.Fatalf("reverse-complement hit missing or projected incorrectly: %+v", hits)
}
