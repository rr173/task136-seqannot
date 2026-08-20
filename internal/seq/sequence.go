package seq

import (
	"errors"
	"fmt"
	"time"
)

// SequenceType distinguishes linear from circular molecules. Circular
// molecules permit features, ORFs, motif hits and restriction sites to wrap
// across the origin; linear molecules reject wraparound coordinates.
type SequenceType string

const (
	TypeLinear   SequenceType = "linear"
	TypeCircular SequenceType = "circular"
)

// Sequence is the engine's canonical representation of a stored DNA molecule.
// Residues are always normalized to upper case and validated against the IUPAC
// alphabet at construction. Coordinates in the public API are 1-based and
// inclusive unless a method documents otherwise; circular wrap is handled by
// WrapIndex under the hood.
type Sequence struct {
	ID         string
	Name       string
	Type       SequenceType
	Residues   string
	Description string
	FastaHeader string
	CreatedAt  time.Time
}

// Length reports the number of residues.
func (s Sequence) Length() int { return len(s.Residues) }

// BaseAt returns the residue at 1-based position p (1..Length). Out-of-range
// positions on a circular sequence wrap; on a linear sequence they error.
func (s Sequence) BaseAt(p int) (byte, error) {
	n := s.Length()
	if n == 0 {
		return 0, errors.New("sequence: empty")
	}
	if s.Type == TypeCircular {
		idx := WrapIndex(p-1, n)
		return s.Residues[idx], nil
	}
	if p < 1 || p > n {
		return 0, fmt.Errorf("sequence: position %d out of range [1,%d]", p, n)
	}
	return s.Residues[p-1], nil
}

// IsCircular reports whether the molecule wraps around the origin.
func (s Sequence) IsCircular() bool { return s.Type == TypeCircular }

// Composition computes the nucleotide counts and the GC fraction. Counts use
// IUPAC semantics: a degenerate base contributes fractionally to each base it
// can represent. GC content is (G+C)/(fully-resolved-bases); ambiguous bases
// contribute their partial weight to G and C. The result is expressed in
// integer basis points (1bp = 0.01%) so callers can compare deterministically
// without floating point; the human-facing percentage is bp/100.
type Composition struct {
	Length   int           `json:"length"`
	Counts   map[byte]int  `json:"counts"` // canonical A/C/G/T counts (fractional contributions rounded)
	GCBp     int           `json:"gc_bp"`   // GC fraction in basis points (0..10000)
	MWBp     int           `json:"mw_bp"`   // molecular weight in milli-Daltons (milli-amu) per strand
	Sequence string        `json:"sequence"`
}

// Compose computes composition statistics for a validated sequence. It never
// mutates s. The molecular weight uses the standard single-stranded DNA
// approximate formula:
//
//	MW(Da) = (#A*313.21) + (#T*304.2) + (#G*329.21) + (#C*289.18) - 61.96
//
// (terminal correction per phosphodiester backbone) and is returned as
// milli-Daltons (Da*1000) rounded half-up to keep integers exact. For a
// circular sequence the terminal correction is replaced by +79.0 (one extra
// phosphodiester bond completes the ring) but the engine reports the linear
// formula per strand to keep the value a pure function of residues + type.
func Compose(s Sequence) Composition {
	n := len(s.Residues)
	counts := map[byte]int{'A': 0, 'C': 0, 'G': 0, 'T': 0}
	gcNum, gcDen := 0, 0
	// Weight accumulators in milli-Da (Da*1000): 313210, 304200, 329210, 289180.
	mwMilli := 0
	for i := 0; i < n; i++ {
		b := s.Residues[i]
		set, ok := iupacSet[b]
		if !ok {
			continue
		}
		members := setMembersCount(set)
		gcDen += members
		var weight int
		if set&bitsA != 0 {
			counts['A']++
			weight += 313210
		}
		if set&bitsC != 0 {
			counts['C']++
			gcNum++
			weight += 289180
		}
		if set&bitsG != 0 {
			counts['G']++
			gcNum++
			weight += 329210
		}
		if set&bitsT != 0 {
			counts['T']++
			weight += 304200
		}
		// each degenerate base contributes its per-member weight once
		mwMilli += weight / members
	}
	mwMilli += terminalCorrectionMilli(s.Type)
	gc := 0
	if gcDen > 0 {
		// round half-up to basis points
		gc = roundDiv(gcNum*10000, gcDen)
	}
	return Composition{
		Length:   n,
		Counts:   counts,
		GCBp:     gc,
		MWBp:     mwMilli,
		Sequence: s.Residues,
	}
}

func setMembersCount(s uint8) int {
	c := 0
	if s&bitsA != 0 {
		c++
	}
	if s&bitsC != 0 {
		c++
	}
	if s&bitsG != 0 {
		c++
	}
	if s&bitsT != 0 {
		c++
	}
	return c
}

// terminalCorrectionMilli returns the per-strand terminal correction for the
// molecular weight formula in milli-Da. Linear single strands lose 61.96 Da
// (one water-equivalent terminal group) -> -61960 milli-Da; circular strands
// add one phosphodiester (~79.0 Da) -> +79000 milli-Da.
func terminalCorrectionMilli(t SequenceType) int {
	if t == TypeCircular {
		return 79000
	}
	return -61960
}

// roundDiv returns round(a/b) with half-up rounding for non-negative a and b>0.
func roundDiv(a, b int) int {
	if b == 0 {
		return 0
	}
	return (a*2 + b) / (b * 2)
}

// NewSequence validates and normalizes a raw residue string into a Sequence
// value (without persisting; persistence is the store's job). It is the single
// construction chokepoint so every Sequence in the system carries validated
// upper-case residues. type defaults to linear when empty.
func NewSequence(name, residues string, typ SequenceType, desc, header string) (Sequence, error) {
	res := Normalize(residues)
	if !IsValidSequence(res) {
		return Sequence{}, fmt.Errorf("seq: invalid base in %q", residues)
	}
	if name == "" {
		return Sequence{}, errors.New("seq: empty name")
	}
	if typ == "" {
		typ = TypeLinear
	}
	if typ != TypeLinear && typ != TypeCircular {
		return Sequence{}, fmt.Errorf("seq: bad type %q", typ)
	}
	return Sequence{
		Name:        name,
		Type:        typ,
		Residues:    res,
		Description: desc,
		FastaHeader: header,
		CreatedAt:   time.Now(), //not used for idempotency; results are residue+param pure
	}, nil
}
