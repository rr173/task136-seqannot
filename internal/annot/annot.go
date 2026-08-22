package annot

import (
	"fmt"

	"task136-seqannot/internal/seq"
)

// FeatureType is a label like "gene", "CDS", "promoter", "exon". The engine
// does not constrain the vocabulary; callers may use any non-empty string.
type FeatureType string

// Strand is '+' or '-'; empty defaults to '+' on persistence.
type Strand byte

const (
	Plus  Strand = '+'
	Minus Strand = '-'
)

// Feature annotates a region of a sequence. Start/End are 1-based inclusive on
// the forward strand. For linear sequences Start<=End; for circular sequences
// Start may be greater than End, meaning the feature wraps around the origin
// (the region is [Start, n] U [1, End]).
type Feature struct {
	ID         string
	SequenceID string
	Type       FeatureType
	Start      int
	End        int
	Strand     Strand
	Label      string
}

// Validate checks the feature against a sequence of length n and circularity.
// Returns an error if coordinates are out of range, start==0, or (for linear
// molecules) start > end.
func (f Feature) Validate(n int, circular bool) error {
	if f.Start < 1 || f.End < 1 {
		return fmt.Errorf("annot: start/end must be >= 1")
	}
	if f.Start > n || f.End > n {
		return fmt.Errorf("annot: coordinate out of range [1,%d]", n)
	}
	if !circular && f.Start > f.End {
		return fmt.Errorf("annot: start>end on linear sequence")
	}
	if f.Strand != 0 && f.Strand != Plus && f.Strand != Minus {
		return fmt.Errorf("annot: bad strand %q", f.Strand)
	}
	return nil
}

// NormalizedStrand returns '+' for a zero strand (treated as forward).
func (f Feature) NormalizedStrand() Strand {
	if f.Strand == 0 {
		return Plus
	}
	return f.Strand
}

// Spans reports whether a 1-based query position q is covered by the feature
// on a molecule of length n. For circular features that wrap (Start>End), the
// covered set is [Start, n] U [1, End].
func (f Feature) Spans(q, n int) bool {
	if q < 1 || q > n {
		return false
	}
	if f.Start <= f.End {
		return q >= f.Start && q <= f.End
	}
	// wrap
	return q >= f.Start || q <= f.End
}

// Overlaps reports whether two features on the same molecule overlap. Two
// features overlap if any 1-based position is shared between them, including
// the circular-wrap case.
func Overlaps(a, b Feature, n int, circular bool) bool {
	for q := 1; q <= n; q++ {
		if a.Spans(q, n) && b.Spans(q, n) {
			return true
		}
	}
	return false
}

// IntervalOverlap is the non-circular O(1) overlap used by the query layer when
// the molecule is linear; it is the standard half-open comparison lifted to
// inclusive coordinates.
func IntervalOverlap(a, b Feature) bool {
	if a.SequenceID != b.SequenceID {
		return false
	}
	return a.Start <= b.End && b.Start <= a.End
}

// QueryOverlap returns the features whose region overlaps the query interval
// (1-based inclusive). For circular molecules a feature wraps when its
// Start>End, and the query interval likewise wraps when qstart>qend (the
// region is [qstart,n] U [1,qend]); only features intersecting that arc are
// returned, not those sitting in the uncovered middle. On a linear molecule a
// reversed query (qstart>qend) is normalized to [qend,qstart]. The check is
// position-scan based for correctness across the wrap boundary.
func QueryOverlap(features []Feature, n int, circular bool, qstart, qend int) []Feature {
	// On a circular molecule qstart>qend is a wrap-around query and must be
	// preserved; swapping it would invert the arc and include the uncovered
	// middle. Only linear queries are normalized.
	if !circular && qstart > qend {
		qstart, qend = qend, qstart
	}
	var out []Feature
	q := Feature{Start: qstart, End: qend}
	for _, f := range features {
		if Overlaps(f, q, n, circular) {
			out = append(out, f)
		}
	}
	return out
}

// Compile-time check that seq is referenced (keeps the import honest even when
// only the sequence length semantics are used).
var _ = seq.TypeLinear
