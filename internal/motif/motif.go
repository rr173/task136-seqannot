package motif

import (
	"fmt"

	"task136-seqannot/internal/seq"
)

// Motif is a named IUPAC pattern (e.g. "RGSCW"). The pattern is compiled into
// a base-set per position for O(n*m) scanning. A motif matches a window of the
// sequence when, at every position, the sequence base's set and the motif
// code's set intersect.
type Motif struct {
	ID          string `json:"id"`
	Name        string `json:"name"`
	Pattern     string `json:"pattern"`
	Description string `json:"description"`
	Sets        []uint8 `json:"-"` // compiled base-set per pattern position
}

// Hit is one match of a motif on a sequence. Coordinates are 1-based inclusive
// on the original forward strand; Strand is '+' for a hit on the forward
// strand and '-' for a hit found on the reverse complement (reported back in
// forward coordinates).
type Hit struct {
	MotifID string `json:"motif_id"`
	Pattern string `json:"pattern"`
	Start   int    `json:"start"`   // 1-based inclusive, forward strand
	End     int    `json:"end"`     // 1-based inclusive, forward strand
	Strand  byte   `json:"strand"`
	Matched string `json:"matched"` // the actual sequence window that matched
}

// ErrBadMotif is returned when a pattern contains a non-IUPAC code.
var ErrBadMotif = fmt.Errorf("motif: pattern contains invalid IUPAC code")

// Compile validates and compiles a motif pattern. The pattern must be non-empty
// and contain only IUPAC nucleotide codes (case-insensitive).
func Compile(id, name, pattern, desc string) (Motif, error) {
	up := seq.Normalize(pattern)
	if up == "" {
		return Motif{}, fmt.Errorf("motif: empty pattern")
	}
	sets := make([]uint8, len(up))
	for i := 0; i < len(up); i++ {
		s, ok := seq.SetFor(up[i])
		if !ok {
			return Motif{}, fmt.Errorf("%w: %q at %d", ErrBadMotif, up[i], i)
		}
		sets[i] = s
	}
	return Motif{ID: id, Name: name, Pattern: up, Description: desc, Sets: sets}, nil
}

// MatchAt reports whether m matches the window of s starting at 0-based offset
// start on strand strand ('+' forward, '-' reverse). On the '+' strand `start`
// is a forward index and the window reads forward positions start, start+1, ...
// On the '-' strand `start` is an index into the reverse complement: the window
// reads forward positions n-1-start, n-1-(start+1), ... (right-to-left, i.e. in
// RC/motif order), taking the complement of each forward base so that the
// matched bytes are the actual bases paired with the motif in its own
// orientation. Forward-strand 1-based inclusive coordinates are projected from
// the RC window by Search. Circular sequences wrap.
func (m Motif) MatchAt(s seq.Sequence, start int, strand byte) (Hit, bool) {
	mlen := len(m.Sets)
	n := s.Length()
	if mlen == 0 || n == 0 {
		return Hit{}, false
	}
	var matched []byte
	for i := 0; i < mlen; i++ {
		// Forward strand: window advances left-to-right over forward
		// positions. Reverse strand: window is the RC slice beginning at RC
		// index `start`, so it walks forward positions right-to-left
		// (n-1-start, n-1-(start+1), ...) to read the bases in motif order.
		var pos int
		if strand == '-' {
			pos = n - 1 - (start + i)
		} else {
			pos = start + i
		}
		if s.IsCircular() {
			pos = seq.WrapIndex(pos, n)
		} else if pos < 0 || pos >= n {
			return Hit{}, false
		}
		var b byte
		if strand == '-' {
			b = seq.Complement(s.Residues[pos])
		} else {
			b = s.Residues[pos]
		}
		bs, ok := seq.SetFor(b)
		if !ok {
			return Hit{}, false
		}
		if bs&m.Sets[i] == 0 {
			return Hit{}, false
		}
		matched = append(matched, b)
	}
	// Compute forward-strand 1-based inclusive coordinates.
	var fwdStart, fwdEnd int
	if strand == '-' {
		// window on rc corresponds to forward positions [n-1-(start+mlen-1) .. n-1-start]
		fwdStart = n - (start + mlen) // 0-based
		fwdEnd = n - 1 - start
		if s.IsCircular() {
			fwdStart = seq.WrapIndex(fwdStart, n)
		}
		fwdStart++
		fwdEnd++
	} else {
		fwdStart = start + 1
		fwdEnd = start + mlen
		if s.IsCircular() {
			fwdEnd = seq.WrapIndex(fwdEnd-1, n) + 1
		}
	}
	return Hit{
		MotifID: m.ID,
		Pattern: m.Pattern,
		Start:   fwdStart,
		End:     fwdEnd,
		Strand:  strand,
		Matched: string(matched),
	}, true
}

// Search scans both strands of s for all non-overlapping-by-start hits of m.
// Every start position (0-based) on the forward strand and on the reverse
// complement strand is tested; on a circular sequence the scan wraps up to one
// full revolution so motifs spanning the origin are found. Hits are returned
// in order of forward-strand start, then '+' before '-' at the same start.
func Search(s seq.Sequence, m Motif) []Hit {
	n := s.Length()
	if n == 0 || len(m.Sets) == 0 {
		return nil
	}
	var hits []Hit
	limit := n
	if s.IsCircular() {
		// allow start positions across a full revolution so wrap-around
		// windows are tested; dedup by (start, strand) below.
		limit = n
	}
	for strandIdx := 0; strandIdx < 2; strandIdx++ {
		strand := byte('+')
		if strandIdx == 1 {
			strand = '-'
		}
		for start := 0; start < limit; start++ {
			if h, ok := m.MatchAt(s, start, strand); ok {
				hits = append(hits, h)
			}
		}
	}
	return hits
}
