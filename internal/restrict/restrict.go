package restrict

import (
	"fmt"

	"task136-seqannot/internal/seq"
)

// Enzyme is a restriction enzyme: a recognition site (IUPAC, possibly
// degenerate) and a cut offset relative to the start of the site. The cut
// offset is 0-based from the 5' end of the recognition site on the forward
// strand; negative offsets mean the enzyme cuts upstream of its site. The
// engine treats every enzyme as a simple "cut at site_start + cut_offset"
// cutter (Type II with a fixed offset); ambiguous multi-cut enzymes are not
// modeled, which keeps digestion deterministic.
type Enzyme struct {
	ID        string `json:"id"`
	Name      string `json:"name"`
	Site      string `json:"site"` // upper-cased IUPAC
	CutOffset int    `json:"cut_offset"` // 0-based offset from site start to the phosphodiester cut
	Sets      []uint8 `json:"-"`
}

// CutSite is one occurrence of an enzyme's recognition site on a sequence,
// with the 1-based inclusive forward-strand coordinates of the site and the
// cut position expressed as a 1-based position: the bond between (cut) and
// (cut+1) is cleaved. For circular sequences coordinates wrap.
type CutSite struct {
	EnzymeID string `json:"enzyme_id"`
	Enzyme   string `json:"enzyme"`
	Start    int    `json:"start"`    // 1-based inclusive, forward strand
	End      int    `json:"end"`      // 1-based inclusive, forward strand
	Strand   byte   `json:"strand"`
	CutPos   int    `json:"cut_pos"` // 1-based; bond (CutPos, CutPos+1) is cut (wraps on circular)
}

// ErrBadSite is returned for a recognition site containing a non-IUPAC code.
var ErrBadSite = fmt.Errorf("restrict: invalid site code")

// Compile validates a recognition site into an Enzyme. The site must be
// non-empty and IUPAC-legal; cut_offset has no range restriction (callers may
// pass negatives).
func Compile(id, name, site string, cut int) (Enzyme, error) {
	up := seq.Normalize(site)
	if up == "" {
		return Enzyme{}, fmt.Errorf("restrict: empty site")
	}
	sets := make([]uint8, len(up))
	for i := 0; i < len(up); i++ {
		s, ok := seq.SetFor(up[i])
		if !ok {
			return Enzyme{}, fmt.Errorf("%w: %q at %d", ErrBadSite, up[i], i)
		}
		sets[i] = s
	}
	return Enzyme{ID: id, Name: name, Site: up, CutOffset: cut, Sets: sets}, nil
}

// FindSites scans both strands of s for every occurrence of e's recognition
// site. Strand '+' hits are on the forward strand; '-' hits are found on the
// reverse complement but reported in forward coordinates. For circular
// sequences the scan wraps one full revolution so sites spanning the origin
// are detected. The order is forward-strand start, then '+' before '-'.
func FindSites(s seq.Sequence, e Enzyme) []CutSite {
	n := s.Length()
	if n == 0 || len(e.Sets) == 0 {
		return nil
	}
	slen := len(e.Sets)
	var sites []CutSite
	for strandIdx := 0; strandIdx < 2; strandIdx++ {
		strand := byte('+')
		if strandIdx == 1 {
			strand = '-'
		}
		for start := 0; start < n; start++ {
			if matchSite(s, start, e.Sets, strand) {
				sites = append(sites, buildCutSite(e, s, start, strand, n))
			}
		}
	}
	// On circular sequences sites whose window wraps past the end are also
	// caught by the start scan because matchSite uses WrapIndex internally.
	_ = slen
	return sites
}

// matchSite reports whether e.Sets matches the window of s starting at 0-based
// `start` on the given strand, wrapping on circular sequences.
func matchSite(s seq.Sequence, start int, sets []uint8, strand byte) bool {
	n := s.Length()
	for i := 0; i < len(sets); i++ {
		var pos int
		if s.IsCircular() {
			pos = seq.WrapIndex(start+i, n)
		} else {
			pos = start + i
			if pos >= n {
				return false
			}
		}
		var b byte
		if strand == '-' {
			b = seq.Complement(s.Residues[pos])
		} else {
			b = s.Residues[pos]
		}
		bs, ok := seq.SetFor(b)
		if !ok {
			return false
		}
		if bs&sets[i] == 0 {
			return false
		}
	}
	return true
}

// buildCutSite projects a 0-based start into the public 1-based CutSite. For
// forward-strand hits the cut is at start + cutOffset (0-based) -> bond between
// that base and the next. For reverse-strand hits the cut mirrors across the
// site.
func buildCutSite(e Enzyme, s seq.Sequence, start int, strand byte, n int) CutSite {
	slen := len(e.Sets)
	var fwdStart, fwdEnd int
	if strand == '-' {
		fwdStart = n - (start + slen)
		fwdEnd = n - 1 - start
		if s.IsCircular() {
			fwdStart = seq.WrapIndex(fwdStart, n)
		}
	} else {
		fwdStart = start
		fwdEnd = start + slen - 1
		if s.IsCircular() {
			fwdEnd = seq.WrapIndex(fwdEnd, n)
		}
	}
	// cut position (1-based): on + strand, site_start_1based + cutOffset gives
	// the base after which the bond is cut. On - strand the cut is mirrored.
	cut1 := fwdStart + 1 + e.CutOffset
	if strand == '-' {
		cut1 = fwdEnd + 1 - e.CutOffset
	}
	if s.IsCircular() {
		cut1 = seq.WrapIndex(cut1-1, n) + 1
	}
	return CutSite{
		EnzymeID: e.ID,
		Enzyme:   e.Name,
		Start:    fwdStart + 1,
		End:      fwdEnd + 1,
		Strand:   strand,
		CutPos:   cut1,
	}
}

// Digest simulates cutting s with a set of enzymes and returns the fragment
// lengths. The fragments are the spans between consecutive cut positions
// (sorted, with circular wraparound handled by connecting the last cut back to
// the first). For a linear sequence with cuts at positions c1<c2<..<ck the
// fragments are [1,c1], (c1,c2], ..., (ck,n]. For a circular sequence the
// fragments are the arcs between consecutive cuts around the ring.
//
// Digestion is deterministic: the fragment list is a pure function of (sequence
// length, sorted unique cut positions, molecule type).
func Digest(s seq.Sequence, enzymes []Enzyme) []int {
	n := s.Length()
	if n == 0 {
		return nil
	}
	var cuts []int
	for _, e := range enzymes {
		for _, cs := range FindSites(s, e) {
			cuts = append(cuts, cs.CutPos)
		}
	}
	cuts = uniqueSorted(cuts)
	if len(cuts) == 0 {
		// no cuts: one whole molecule
		return []int{n}
	}
	var fragments []int
	if s.IsCircular() {
		// arcs between consecutive cuts, wrapping
		for i := 0; i < len(cuts); i++ {
			next := cuts[(i+1)%len(cuts)]
			cur := cuts[i]
			d := next - cur
			if d <= 0 {
				d += n
			}
			fragments = append(fragments, d)
		}
	} else {
		// linear: head fragment [1..cuts[0]], middles, tail
		fragments = append(fragments, cuts[0])
		for i := 1; i < len(cuts); i++ {
			fragments = append(fragments, cuts[i]-cuts[i-1])
		}
		fragments = append(fragments, n-cuts[len(cuts)-1])
	}
	return fragments
}

func uniqueSorted(xs []int) []int {
	if len(xs) == 0 {
		return xs
	}
	// insertion sort (small slices) + dedup
	out := make([]int, 0, len(xs))
	for _, v := range xs {
		inserted := false
		for i, w := range out {
			if v < w {
				out = append(out, 0)
				copy(out[i+1:], out[i:])
				out[i] = v
				inserted = true
				break
			}
			if v == w {
				inserted = true
				break
			}
		}
		if !inserted {
			out = append(out, v)
		}
	}
	return out
}
