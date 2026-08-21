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
//
// For a circular molecule, the frame is a closed loop: if the offset leaves a
// partial codon at the 3' end of the linear pass, that codon is completed by
// wrapping around the origin (e.g. for GTAACCCAT in forward offset 1 the last
// codon is positions 8,9,1 = ATG). Without this completion an ORF whose start
// sits at the very end of the molecule and whose stop lands across the origin
// is never emitted - the start codon would not exist in any frame's codon
// list. Linear sequences are left unchanged: their trailing partial codon is
// not translatable.
func Translate6(s Sequence) []FrameInfo {
	n := len(s.Residues)
	if n < 3 {
		return []FrameInfo{}
	}
	frames := make([]FrameInfo, 0, 6)
	for off := 0; off < 3; off++ {
		frames = append(frames, translateForward(s, off))
	}
	rc := ReverseComplement(s.Residues)
	for off := 0; off < 3; off++ {
		frames = append(frames, translateReverse(s, rc, off))
	}
	return frames
}

// translateForward builds a '+' frame starting at offset off on the forward
// strand. Linear codon coordinates are [off+3k+1, off+3k+3]. For a circular
// molecule the trailing partial codon (when (n-off)%3 != 0) is completed from
// the origin and becomes an origin-spanning codon whose End wraps below Start.
func translateForward(s Sequence, off int) FrameInfo {
	res := s.Residues
	n := len(res)
	aminos := make([]byte, 0, (n-off)/3+1)
	codons := make([]Codon, 0, (n-off)/3+1)
	i := off
	for ; i+2 < n; i += 3 {
		c0, c1, c2 := res[i], res[i+1], res[i+2]
		am, resolv := CodonAmino(c0, c1, c2)
		three := [3]byte{upperByte(c0), upperByte(c1), upperByte(c2)}
		codons = append(codons, Codon{
			Frame:    off,
			Strand:   '+',
			Start:    i + 1,
			End:      i + 3,
			Amino:    am,
			IsStart:  IsStartCodon(three),
			IsStop:   IsStopCodon(three),
			Bases:    string(three[:]),
			Resolved: resolv,
		})
		aminos = append(aminos, byte(am))
	}
	// circular: complete the origin-spanning codon from the partial tail.
	if s.Type == TypeCircular {
		rem := (n - off) % 3
		if rem != 0 {
			tail := n - rem // 0-based index of the first trailing base
			b0 := res[WrapIndex(tail, n)]
			b1 := res[WrapIndex(tail+1, n)]
			b2 := res[WrapIndex(tail+2, n)]
			am, resolv := CodonAmino(b0, b1, b2)
			three := [3]byte{upperByte(b0), upperByte(b1), upperByte(b2)}
			codons = append(codons, Codon{
				Frame:    off,
				Strand:   '+',
				Start:    tail + 1, // 1-based start of the trailing base run
				End:      WrapIndex(tail+2, n) + 1,
				Amino:    am,
				IsStart:  IsStartCodon(three),
				IsStop:   IsStopCodon(three),
				Bases:    string(three[:]),
				Resolved: resolv,
			})
			aminos = append(aminos, byte(am))
		}
	}
	return FrameInfo{Frame: off, Strand: '+', Offset: off, Aminos: string(aminos), Codons: codons}
}

// translateReverse builds a '-' frame. Translation happens on the reverse
// complement, but each codon's coordinates must point back to the original
// forward strand. For an rc codon k (0-based) at rc-offset `off`, the codon
// occupies rc positions [off+3k .. off+3k+2] on the rc string, which map to
// forward positions [n-1-(off+3k+2) .. n-1-(off+3k)] = [n-off-3k-3 ..
// n-off-3k-1] (0-based). We report these as 1-based inclusive.
//
// For a circular molecule the trailing partial rc codon is completed by
// wrapping the rc string at its origin; the resulting codon spans the forward
// origin (End may wrap below Start, as on the forward strand).
func translateReverse(s Sequence, rc string, off int) FrameInfo {
	n := len(rc)
	aminos := make([]byte, 0, (n-off)/3+1)
	codons := make([]Codon, 0, (n-off)/3+1)
	i := off
	for ; i+2 < n; i += 3 {
		c0, c1, c2 := rc[i], rc[i+1], rc[i+2]
		am, resolv := CodonAmino(c0, c1, c2)
		three := [3]byte{upperByte(c0), upperByte(c1), upperByte(c2)}
		// forward 0-based positions for this rc codon:
		fwdEnd := n - 1 - i // last (3') base on forward strand
		fwdStart := n - 1 - (i + 2) // first (5') base on forward strand
		if fwdStart < 0 {
			break
		}
		codons = append(codons, Codon{
			Frame:    off + 3,
			Strand:   '-',
			Start:    fwdStart + 1,
			End:      fwdEnd + 1,
			Amino:    am,
			IsStart:  IsStartCodon(three),
			IsStop:   IsStopCodon(three),
			Bases:    string(three[:]),
			Resolved: resolv,
		})
		aminos = append(aminos, byte(am))
	}
	// circular: complete the origin-spanning rc codon from the partial tail.
	if s.Type == TypeCircular {
		rem := (n - off) % 3
		if rem != 0 {
			tail := n - rem // 0-based rc index of the first trailing rc base
			b0 := rc[WrapIndex(tail, n)]
			b1 := rc[WrapIndex(tail+1, n)]
			b2 := rc[WrapIndex(tail+2, n)]
			am, resolv := CodonAmino(b0, b1, b2)
			three := [3]byte{upperByte(b0), upperByte(b1), upperByte(b2)}
			// rc positions tail..tail+2 wrap to forward positions, 5'..3'
			// on the reverse strand == smallest forward idx (after wrap) ..
			// largest forward idx.
			fwdStart := n - 1 - WrapIndex(tail+2, n) // 5' base (first translated)
			fwdEnd := n - 1 - WrapIndex(tail, n)     // 3' base (last translated)
			codons = append(codons, Codon{
				Frame:    off + 3,
				Strand:   '-',
				Start:    fwdStart + 1,
				End:      fwdEnd + 1,
				Amino:    am,
				IsStart:  IsStartCodon(three),
				IsStop:   IsStopCodon(three),
				Bases:    string(three[:]),
				Resolved: resolv,
			})
			aminos = append(aminos, byte(am))
		}
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
