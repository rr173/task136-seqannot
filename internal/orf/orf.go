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
	n := s.Length()
	var orfs []ORF
	for _, f := range frames {
		orfs = append(orfs, findInFrame(f, s, n, opts)...)
	}
	return orfs
}

// findInFrame walks one frame's codons. It tracks an open ORF (start codon
// index) and closes it at the first stop. For circular molecules, after the
// linear pass the frame may continue from the beginning until it reaches its
// own start or completes a full revolution.
func findInFrame(f seq.FrameInfo, s seq.Sequence, n int, opts Options) []ORF {
	if len(f.Codons) == 0 {
		return nil
	}
	circular := s.IsCircular()
	var out []ORF
	openIdx := -1
	codons := f.Codons

	closeORF := func(stopIdx int) {
		if openIdx < 0 {
			return
		}
		startCodon := codons[openIdx]
		stopCodon := codons[stopIdx]
		prot := buildProtein(codons[openIdx : stopIdx+1])
		aaLen := stopIdx - openIdx // excludes stop
		if opts.MinAALen > 0 && aaLen < opts.MinAALen {
			openIdx = -1
			return
		}
		out = append(out, ORF{
			Frame:      f.Frame,
			Strand:     '+',
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
			closeORF(i)
		}
	}

	// partial ORF (no stop found): report if it clears MinAALen
	if openIdx >= 0 {
		lastIdx := len(codons) - 1
		aaLen := lastIdx - openIdx
		if opts.MinAALen == 0 || aaLen >= opts.MinAALen {
			startCodon := codons[openIdx]
			prot := buildProtein(codons[openIdx:])
			out = append(out, ORF{
				Frame:      f.Frame,
				Strand:     '+',
				Start:      startCodon.Start,
				End:        codons[lastIdx].End,
				Length:     aaLen + 1,
				StartCodon: startCodon.Bases,
				StopCodon:  "",
				Protein:    prot,
			})
		}
		openIdx = -1
	}

	// circular wrap: continue the frame across the origin until we either
	// find a stop or return to the openIdx position (full revolution).
	if circular && openIdx >= 0 {
		// reopen at the wrapped start: the partial we just emitted already
		// covered the linear tail; a circular ORF that wraps needs the stop
		// to be found before completing a full lap. We re-scan from codon 0.
		wrapped := scanCircularWrap(codons, n, f.Frame, opts)
		out = append(out, wrapped...)
	}
	return out
}

// scanCircularWrap handles the rare case of an ORF that opens near the 3' end
// of a circular genome and closes only after wrapping past the origin. Because
// reverse frames already project onto forward coordinates and a full-lap scan
// would re-report the same linear ORFs, we restrict this to ORFs that are
// still open at the end of the linear pass (no stop found at all).
func scanCircularWrap(codons []seq.Codon, n, frame int, opts Options) []ORF {
	// find the last start with no stop after it
	if len(codons) == 0 {
		return nil
	}
	var openIdx = -1
	for i := 0; i < len(codons); i++ {
		if codons[i].IsStart {
			openIdx = i
		}
		if codons[i].IsStop {
			openIdx = -1
		}
	}
	if openIdx < 0 {
		return nil
	}
	// wrap: continue from codon 0 until a stop or full revolution
	startCodon := codons[openIdx]
	var prot []byte
	prot = append(prot, startCodon.Bases...)
	total := len(codons) - openIdx
	i := 0
	for steps := 0; steps < len(codons); steps++ {
		c := codons[i]
		if c.IsStop {
			aaLen := total - 1 // exclude stop
			if opts.MinAALen > 0 && aaLen < opts.MinAALen {
				return nil
			}
			// End wraps: stop codon's forward End is < startCodon.Start
			return []ORF{{
				Frame:      frame,
				Strand:     '+',
				Start:      startCodon.Start,
				End:        c.End,
				Length:     aaLen + 1,
				StartCodon: startCodon.Bases,
				StopCodon:  c.Bases,
				Protein:    string(prot),
			}}
		}
		total++
		prot = append(prot, c.Bases...)
		i = (i + 1) % len(codons)
	}
	return nil // full revolution, no stop: not a valid ORF
}

// buildProtein concatenates the amino letters of a codon slice, preserving
// unresolved codons as 'X'.
func buildProtein(codons []seq.Codon) string {
	b := make([]byte, 0, len(codons))
	for _, c := range codons {
		if c.Resolved {
			b = append(b, byte(c.Amino))
		} else {
			b = append(b, 'X')
		}
	}
	return string(b)
}
