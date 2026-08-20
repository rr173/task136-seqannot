package motif

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestCompileAndSearch(t *testing.T) {
	m, err := Compile("m1", "RBS", "RBS", "")
	if err != nil {
		t.Fatal(err)
	}
	// R = A|G, B = C|G|T, S = G|C
	// sequence "AGC" -> R matches A, B matches G, S matches C -> hit at 1
	s, _ := seq.NewSequence("t", "AGC", seq.TypeLinear, "", "")
	hits := Search(s, m)
	if len(hits) != 1 {
		t.Fatalf("want 1 hit got %d: %+v", len(hits), hits)
	}
	if hits[0].Start != 1 || hits[0].End != 3 || hits[0].Strand != '+' {
		t.Fatalf("hit=%+v", hits[0])
	}
}

func TestSearchBothStrands(t *testing.T) {
	// pattern GG; GG on + strand (pos 1-2), and CC on + (pos 3-4) pairs as GG
	// on the reverse strand, so we expect at least one '+' and one '-' hit.
	m, _ := Compile("m", "gg", "GG", "")
	s, _ := seq.NewSequence("t", "GGCC", seq.TypeLinear, "", "")
	hits := Search(s, m)
	plus, minus := 0, 0
	for _, h := range hits {
		if h.Strand == '+' {
			plus++
		} else {
			minus++
		}
	}
	if plus == 0 {
		t.Fatal("expected a + strand hit")
	}
	if minus == 0 {
		t.Fatal("expected a - strand hit")
	}
}

func TestCompileInvalid(t *testing.T) {
	if _, err := Compile("m", "X", "ZQ", ""); err == nil {
		t.Fatal("expected error for invalid code")
	}
}

func TestSearchCircularWrap(t *testing.T) {
	// circular: pattern GGG; sequence "GGAGG" circular should find a wrap hit
	m, _ := Compile("m", "g3", "GGG", "")
	s, _ := seq.NewSequence("t", "GGAGG", seq.TypeCircular, "", "")
	hits := Search(s, m)
	if len(hits) == 0 {
		t.Fatal("circular should find wrap hits")
	}
}
