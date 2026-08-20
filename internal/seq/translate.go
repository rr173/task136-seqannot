package seq

import "fmt"

// Frame identifies one of the six reading frames. Frames 0..2 are the forward
// strands starting at offset 0,1,2; frames 3..5 are the reverse-complement
// strand starting at offset 0,1,2 of ReverseComplement(s). The engine reports
// every frame's amino-acid translation and ORFs are derived per frame.
type Frame int

// FrameInfo describes a single translated frame.
type FrameInfo struct {
	Frame  int     `json:"frame"`  // 0..5
	Strand byte    `json:"strand"`  // '+' or '-'
	Offset int     `json:"offset"`  // 0,1,2
	Aminos string  `json:"aminos"`  // amino acid string; trailing partial codon is dropped
	Codons []Codon `json:"codons"` // one per full codon in frame order
}

// Codon annotates a single translated codon with its amino acid and the
// 1-based genomic coordinates it covers on the original (forward) strand.
type Codon struct {
	Frame    int   `json:"frame"`
	Strand   byte  `json:"strand"`   // '+' or '-'
	Start    int   `json:"start"`     // 1-based inclusive on the original strand
	End      int   `json:"end"`       // 1-based inclusive
	Amino    Amino `json:"amino"`
	IsStart  bool  `json:"is_start"`
	IsStop   bool  `json:"is_stop"`
	Bases    string `json:"bases"`
	Resolved bool  `json:"resolved"` // false when ambiguity prevented a single amino
}

// Translate6 computes all six reading frames for a sequence. For reverse
// frames the translation is performed on ReverseComplement(s) but every
// codon's Start/End coordinates are projected back onto the original forward
// strand, so callers always address the genome in a single coordinate space.
func Translate6(s Sequence) []FrameInfo {
	n := len(s.Residues)
	if n < 3 {
		return []FrameInfo{}
	}
	frames := make([]FrameInfo, 0, 6)
	for off := 0; off < 3; off++ {
		frames = append(frames, translateForward(s.Residues, n, off))
	}
	rc := ReverseComplement(s.Residues)
	for off := 0; off < 3; off++ {
		frames = append(frames, translateReverse(rc, n, off))
	}
	return frames
}

// translateForward builds a '+' frame starting at offset off on the forward
// strand. Codon coordinates are simply [off+3k+1, off+3k+3].
func translateForward(res string, n, off int) FrameInfo {
	aminos := make([]byte, 0, (n-off)/3)
	codons := make([]Codon, 0, (n-off)/3)
	for i := off; i+2 < n; i += 3 {
		c0, c1, c2 := res[i], res[i+1], res[i+2]
		am, resolv := CodonAmino(c0, c1, c2)
		three := [3]byte{upperByte(c0), upperByte(c1), upperByte(c2)}
		codons = append(codons, Codon{
			Frame: off,
			Start: i + 1,
			End:   i + 3,
			Amino: am,
			IsStart: IsStartCodon(three),
			IsStop:  IsStopCodon(three),
			Bases:   string(three[:]),
			Resolved: resolv,
		})
		aminos = append(aminos, byte(am))
	}
	return FrameInfo{Frame: off, Strand: '+', Offset: off, Aminos: string(aminos), Codons: codons}
}

// translateReverse builds a '-' frame. Translation happens on the reverse
// complement, but each codon's coordinates must point back to the original
// forward strand. For an rc codon k (0-based) at rc-offset `off`, the codon
// occupies rc positions [off+3k .. off+3k+2] on the rc string, which map to
// forward positions [n-1-(off+3k+2) .. n-1-(off+3k)] = [n-off-3k-3 ..
// n-off-3k-1] (0-based). We report these as 1-based inclusive.
func translateReverse(rc string, n, off int) FrameInfo {
	aminos := make([]byte, 0, (n-off)/3)
	codons := make([]Codon, 0, (n-off)/3)
	for i := off; i+2 < len(rc); i += 3 {
		c0, c1, c2 := rc[i], rc[i+1], rc[i+2]
		am, resolv := CodonAmino(c0, c1, c2)
		three := [3]byte{upperByte(c0), upperByte(c1), upperByte(c2)}
		// forward 0-based positions for this rc codon:
		fwdEnd := n - 1 - i       // last (3') base on forward strand
		fwdStart := n - 1 - (i + 2) // first (5') base on forward strand
		if fwdStart < 0 {
			break
		}
		codons = append(codons, Codon{
			Frame: off + 3,
			Strand: '-',
			Start: fwdStart + 1,
			End:   fwdEnd + 1,
			Amino: am,
			IsStart: IsStartCodon(three),
			IsStop:  IsStopCodon(three),
			Bases:   string(three[:]),
			Resolved: resolv,
		})
		aminos = append(aminos, byte(am))
	}
	return FrameInfo{Frame: off + 3, Strand: '-', Offset: off, Aminos: string(aminos), Codons: codons}
}

// FrameFor returns the FrameInfo for the given 0..5 frame index, or an error
// if the index is out of range.
func FrameFor(frames []FrameInfo, frame int) (FrameInfo, error) {
	for _, f := range frames {
		if f.Frame == frame {
			return f, nil
		}
	}
	return FrameInfo{}, fmt.Errorf("seq: frame %d not found", frame)
}
