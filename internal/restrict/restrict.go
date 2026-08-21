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
	ID         string `json:"id"`
	Name       string `json:"name"`
	Site       string `json:"site"`        // upper-cased IUPAC
	CutOffset  int    `json:"cut_offset"`  // 0-based offset from site start to the phosphodiester cut
	Sets       []uint8 `json:"-"`
	palindrome bool   // caches whether Site is its own reverse complement (IUPAC)
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
	return Enzyme{ID: id, Name: name, Site: up, CutOffset: cut, Sets: sets, palindrome: seq.ReverseComplement(up) == up}, nil
}

// FindSites scans both strands of s for every occurrence of e's recognition
// site. Strand '+' hits are on the forward strand; '-' hits are found on the
// reverse complement but reported in forward coordinates. For circular
// sequences the scan wraps one full revolution so sites spanning the origin
// are detected. The order is forward-strand start, then '+' before '-'.
//
// When the recognition site is its own reverse complement (a palindromic site
// such as GAATTC), reading the reverse strand reproduces the very same window
// as the forward strand, so a '-' match at any start is the same physical cut
// as the '+' match there — not a second cut. Scanning '-' in that case would
// invent a non-existent pseudo-hit, so the '-' strand is skipped for
// palindromic enzymes. Non-palindromic sites (e.g. TCA, GG) have distinct
// reverse-strand occurrences that are real, separate cuts and are reported.
func FindSites(s seq.Sequence, e Enzyme) []CutSite {
	n := s.Length()
	if n == 0 || len(e.Sets) == 0 {
		return nil
	}
	var sites []CutSite
	for strandIdx := 0; strandIdx < 2; strandIdx++ {
		strand := byte('+')
		if strandIdx == 1 {
			strand = '-'
		}
		if strand == '-' && e.palindrome {
			// reverse-strand scan would duplicate forward hits
			continue
		}
		for start := 0; start < n; start++ {
			if matchSite(s, start, e.Sets, strand) {
				sites = append(sites, buildCutSite(e, s, start, strand, n))
			}
		}
	}
	// On circular sequences sites whose window wraps past the end are caught by
	// the start scan because matchSite uses WrapIndex internally, so no extra
	// wrap-around pass is needed.
	return sites
}

// matchSite reports whether e.Sets matches the window of s starting at 0-based
// `start` on the given strand, wrapping on circular sequences.
//
// On the '+' strand the recognition site is read 5'->3' along the forward
// window [start, start+slen). On the '-' strand the site is read 5'->3' along
// the reverse complement of that same window, i.e. the i-th site base is
// compared against Complement(forward[start+slen-1-i]) — the window is matched
// in reverse-complement order, not merely complemented. Matching only the
// complement (without reversal) misses real reverse-strand sites and invents
// pseudo-hits that have no biological reverse-strand site, so this is the
// correct projection of a reverse-strand recognition event.
func matchSite(s seq.Sequence, start int, sets []uint8, strand byte) bool {
	n := s.Length()
	slen := len(sets)
	for i := 0; i < slen; i++ {
		// Window offset of the forward base compared against site position i.
		// Forward strand reads the window left-to-right (offset i); reverse
		// strand reads it right-to-left (offset slen-1-i), complemented.
		off := i
		if strand == '-' {
			off = slen - 1 - i
		}
		var pos int
		if s.IsCircular() {
			pos = seq.WrapIndex(start+off, n)
		} else {
			pos = start + off
			if pos >= n {
				return false
			}
		}
		b := s.Residues[pos]
		if strand == '-' {
			b = seq.Complement(b)
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

// buildCutSite projects a 0-based window start into the public 1-based CutSite.
// The matched window always spans forward positions [start, start+slen) (wrapped
// on circular sequences). For a forward-strand hit the recognition site reads
// left-to-right over that window, so the cut offset is measured from fwdStart.
// For a reverse-strand hit the site reads right-to-left over the same window
// (its reverse complement), so the cut offset mirrors to be measured from
// fwdEnd: the phosphodiester cut lies at fwdEnd - cutOffset, consistent with
// the recognition direction. Both branches yield a 1-based CutPos with the bond
// (CutPos, CutPos+1) cleaved (wrapping on circular sequences).
func buildCutSite(e Enzyme, s seq.Sequence, start int, strand byte, n int) CutSite {
	slen := len(e.Sets)
	fwdStart := start
	fwdEnd := start + slen - 1
	if s.IsCircular() {
		fwdEnd = seq.WrapIndex(fwdEnd, n)
	}
	// cut position (1-based): on '+' the offset runs forward from fwdStart; on
	// '-' it mirrors and runs back from fwdEnd.
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
