package orf

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestBug01_CircularORFWrapsOrigin(t *testing.T) {
	s, err := seq.NewSequence("circular", "GTAACCCAT", seq.TypeCircular, "", "")
	if err != nil {
		t.Fatal(err)
	}
	got := Find(s, Options{MinAALen: 0})
	for _, o := range got {
		if o.Frame == 1 && o.Start == 8 && o.End == 4 && o.Protein == "M*" {
			return
		}
	}
	t.Fatalf("circular ORF crossing origin missing: %+v", got)
}
