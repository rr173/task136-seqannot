package seq

import "testing"

func TestBug05_FastaStripsAllSequenceWhitespace(t *testing.T) {
	recs, err := ParseFasta(">x\nATG\tC\nA C\n")
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 1 || recs[0].Residues != "ATGCAC" {
		t.Fatalf("unexpected parser result: %+v", recs)
	}
}
