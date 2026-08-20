package restrict

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestBug03_CircularReverseCutCoordinates(t *testing.T) {
	s, err := seq.NewSequence("circular", "TGAAC", seq.TypeCircular, "", "")
	if err != nil {
		t.Fatal(err)
	}
	e, err := Compile("e1", "enzyme", "TCA", 1)
	if err != nil {
		t.Fatal(err)
	}
	for _, site := range FindSites(s, e) {
		if site.Strand == '-' && site.Start == 1 && site.End == 3 && site.CutPos == 2 {
			return
		}
	}
	t.Fatalf("reverse circular cut projection missing: %+v", FindSites(s, e))
}
