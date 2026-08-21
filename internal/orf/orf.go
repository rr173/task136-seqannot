package orf

import (
	"task136-seqannot/internal/seq"
)

// ORF is an Open Reading Frame: a maximal run of codons starting at a start
// codon (ATG) and ending at the first in-frame stop codon (inclusive of the
// stop). Coordinates are 1-based inclusive on the original forward strand;
// the strand is '+' for forward frames and '-' for reverse frames.
type ORF struct {
	Frame      int    `json:"frame"`       // 0..5
	Strand     byte   `json:"strand"`      // '+' or '-'
	Start      int    `json:"start"`       // 1-based inclusive, original forward strand
	End        int    `json:"end"`         // 1-based inclusive, original forward strand
	Length     int    `json:"length"`      // amino acid length INCLUDING the stop residue (*)
	StartCodon string `json:"start_codon"`
	StopCodon  string `json:"stop_codon"`
	Protein    string `json:"protein"`    // amino string including trailing '*'
}

// Options configures ORF detection.
type Options struct {
	MinAALen int // minimum amino-acid length excluding the stop codon; 0 means no floor
}

// Find scans a sequence for ORFs in all six frames and returns every ORF that
// is at least opts.MinAALen amino acids long (the stop does not count toward
// the length). Overlapping ORFs in different (or the same) frames are all
// reported; the engine never merges or drops overlapping hits. For circular
// sequences, an ORF whose stop codon sits past the origin wraps around and its
// End coordinate is reported modulo length.
//
// Detection rule per frame: walk codons left-to-right; on a start codon open
// an ORF; on the next in-frame stop codon close it. Nested starts (a second
// start before a stop) are NOT reported separately - only the outermost
// (earliest) start is kept, matching the standard "longest ORF from first
// ATG" convention. If a frame has no stop before the sequence ends (or, for
// circular genomes, wraps a full lap without finding a stop), the run is
// reported as a partial ORF with an empty StopCodon and End set to the last
// codon position; partials still respect MinAALen.
func Find(s seq.Sequence, opts Options) []ORF {
	frames := seq.Translate6(s)
	var orfs []ORF
	for _, f := range frames {
		orfs = append(orfs, findInFrame(f, s, opts)...)
	}
	return orfs
}

// findInFrame walks one frame's codons. It tracks an open ORF (start codon
// index) and closes it at the first in-frame stop. For a circular molecule
// the codon list already includes the origin-spanning codon completed by the
// translator, so a start sitting at the very end of the molecule with its stop
// across the origin is visible here: the forward pass opens it, and the
// wrap-around pass closes it by continuing from codon 0. The wrap stops if it
// finds a stop or returns to the open start (a full lap with no stop), in
// which case the run is not a bounded ORF and is dropped.
func findInFrame(f seq.FrameInfo, s seq.Sequence, opts Options) []ORF {
	if len(f.Codons) == 0 {
		return nil
	}
	circular := s.IsCircular()
	codons := f.Codons
	var out []ORF
	openIdx := -1

	// emit closes the open ORF at stopIdx. When wrap is true the stop sits
	// across the origin: the protein runs codons[openIdx:] then codons[:stopIdx+1].
	emit := func(stopIdx int, wrap bool) {
		if openIdx < 0 {
			return
		}
		startCodon := codons[openIdx]
		stopCodon := codons[stopIdx]
		prot := buildProteinWrap(codons, openIdx, stopIdx, wrap)
		aaLen := codonSpan(openIdx, stopIdx, len(codons), wrap) - 1 // exclude stop
		if opts.MinAALen > 0 && aaLen < opts.MinAALen {
			openIdx = -1
			return
		}
		out = append(out, ORF{
			Frame:      f.Frame,
			Strand:     f.Strand,
			Start:      startCodon.Start,
			End:        stopCodon.End,
			Length:     aaLen + 1, // include stop
			StartCodon: startCodon.Bases,
			StopCodon:  stopCodon.Bases,
			Protein:    prot,
		})
		openIdx = -1
	}

	// forward pass
	for i := 0; i < len(codons); i++ {
		c := codons[i]
		if c.IsStart && openIdx < 0 {
			openIdx = i
		}
		if c.IsStop && openIdx >= 0 {
			emit(i, false)
		}
	}

	// still open at the end of the linear pass
	if openIdx >= 0 {
		if circular {
			// wrap across the origin: scan codon (openIdx+1) mod L, then onward,
			// until a stop closes the ORF or we circle back to openIdx (no stop).
			for steps := 1; steps < len(codons); steps++ {
				j := (openIdx + steps) % len(codons)
				c := codons[j]
				if c.IsStop {
					emit(j, true)
					break
				}
			}
			openIdx = -1 // full lap with no stop: not a bounded ORF
		} else {
			// partial ORF (no stop found): report if it clears MinAALen
			lastIdx := len(codons) - 1
			aaLen := lastIdx - openIdx
			if opts.MinAALen == 0 || aaLen >= opts.MinAALen {
				startCodon := codons[openIdx]
				out = append(out, ORF{
					Frame:      f.Frame,
					Strand:     f.Strand,
					Start:      startCodon.Start,
					End:        codons[lastIdx].End,
					Length:     aaLen + 1,
					StartCodon: startCodon.Bases,
					StopCodon:  "",
					Protein:    buildProteinWrap(codons, openIdx, lastIdx, false),
				})
			}
			openIdx = -1
		}
	}
	return out
}

// codonSpan counts the codons in an ORF from startIdx to stopIdx inclusive.
// Without wrap this is stopIdx-startIdx+1; with wrap the span crosses the end
// of the list: codons[startIdx:] then codons[:stopIdx+1].
func codonSpan(startIdx, stopIdx, length int, wrap bool) int {
	if !wrap {
		return stopIdx - startIdx + 1
	}
	return (length - startIdx) + (stopIdx + 1)
}

// buildProtein concatenates the amino letters of a codon slice, preserving
// unresolved codons as 'X'.
func buildProtein(codons []seq.Codon) string {
	b := make([]byte, 0, len(codons))
	for _, c := range codons {
		b = appendCodonAA(b, c)
	}
	return string(b)
}

// buildProteinWrap builds the amino string for an ORF whose codon run may wrap
// around the end of the codon list (an origin-spanning ORF). For a non-wrapping
// run it is equivalent to buildProtein(codons[startIdx:stopIdx+1]).
func buildProteinWrap(codons []seq.Codon, startIdx, stopIdx int, wrap bool) string {
	if !wrap {
		return buildProtein(codons[startIdx : stopIdx+1])
	}
	b := make([]byte, 0, len(codons))
	for i := startIdx; i < len(codons); i++ {
		b = appendCodonAA(b, codons[i])
	}
	for i := 0; i <= stopIdx; i++ {
		b = appendCodonAA(b, codons[i])
	}
	return string(b)
}

func appendCodonAA(b []byte, c seq.Codon) []byte {
	if c.Resolved {
		return append(b, byte(c.Amino))
	}
	return append(b, 'X')
}
