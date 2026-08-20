package restrict

import (
	"testing"

	"task136-seqannot/internal/seq"
)

func TestFindSitesLinear(t *testing.T) {
	// EcoRI: GAATTC, cut offset 1 (G|AATTC)
	e, _ := Compile("e1", "EcoRI", "GAATTC", 1)
	s, _ := seq.NewSequence("t", "GAATTCGAATTC", seq.TypeLinear, "", "")
	sites := FindSites(s, e)
	if len(sites) != 2 {
		t.Fatalf("want 2 sites got %d: %+v", len(sites), sites)
	}
	if sites[0].Start != 1 || sites[0].End != 6 {
		t.Fatalf("site0=%+v", sites[0])
	}
	if sites[1].Start != 7 || sites[1].End != 12 {
		t.Fatalf("site1=%+v", sites[1])
	}
}

func TestFindSitesDegenerate(t *testing.T) {
	// site "GGNCC" (BsrI-like) - N matches anything
	e, _ := Compile("e", "BsrI", "GGNCC", 1)
	s, _ := seq.NewSequence("t", "GGACC", seq.TypeLinear, "", "")
	sites := FindSites(s, e)
	if len(sites) != 1 {
		t.Fatalf("want 1 degenerate site got %d: %+v", len(sites), sites)
	}
}

func TestFindSitesReverseStrand(t *testing.T) {
	// site GG; sequence GGCC has GG on + (2-3) and CC on + (3-4) which pairs
	// as GG on the reverse strand -> a '-' strand site must be present.
	e, _ := Compile("e", "GG", "GG", 0)
	s, _ := seq.NewSequence("t", "GGCC", seq.TypeLinear, "", "")
	sites := FindSites(s, e)
	hasMinus := false
	for _, st := range sites {
		if st.Strand == '-' {
			hasMinus = true
		}
	}
	if !hasMinus {
		t.Fatalf("expected a reverse-strand site among %+v", sites)
	}
}

func TestDigestLinear(t *testing.T) {
	e, _ := Compile("e", "EcoRI", "GAATTC", 1)
	s, _ := seq.NewSequence("t", "GAATTCGAATTC", seq.TypeLinear, "", "")
	frags := Digest(s, []Enzyme{e})
	// cuts at position 2 and 8 (1-based, bond after). fragments: 2, 6, 4
	want := []int{2, 6, 4}
	if len(frags) != len(want) {
		t.Fatalf("frags=%v want %v", frags, want)
	}
	for i := range want {
		if frags[i] != want[i] {
			t.Fatalf("frags=%v want %v", frags, want)
		}
	}
}

func TestDigestCircular(t *testing.T) {
	e, _ := Compile("e", "EcoRI", "GAATTC", 1)
	// circular with two sites -> two arcs
	s, _ := seq.NewSequence("t", "GAATTCGAATTC", seq.TypeCircular, "", "")
	frags := Digest(s, []Enzyme{e})
	if len(frags) != 2 {
		t.Fatalf("circular frags=%v want 2", frags)
	}
	sum := 0
	for _, f := range frags {
		sum += f
	}
	if sum != s.Length() {
		t.Fatalf("fragment sum %d != length %d", sum, s.Length())
	}
}

func TestDigestNoCuts(t *testing.T) {
	e, _ := Compile("e", "EcoRI", "GAATTC", 1)
	s, _ := seq.NewSequence("t", "ACGTACGT", seq.TypeLinear, "", "")
	frags := Digest(s, []Enzyme{e})
	if len(frags) != 1 || frags[0] != 8 {
		t.Fatalf("no-cut frags=%v want [8]", frags)
	}
}
