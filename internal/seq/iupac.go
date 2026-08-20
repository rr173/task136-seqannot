package seq

import "strings"

// Base bitmask over the four canonical DNA nucleotides.
const (
	bitsA = 1 << iota
	bitsC
	bitsG
	bitsT
)

// IUPAC nucleotide code -> set of matching canonical bases (as a bitmask).
//
// A=adenine, C=cytosine, G=guanine, T=thymine, R=puRine(A|G), Y=pYrimidine(C|T),
// S=Strong(G|C), W=Weak(A|T), K=Keto(G|T), M=aMino(A|C), B=not A(C|G|T),
// D=not C(A|G|T), H=not G(A|C|T), V=not T(A|C|G), N=aNy(A|C|G|T).
var iupacSet = map[byte]uint8{
	'A': bitsA,
	'C': bitsC,
	'G': bitsG,
	'T': bitsT,
	'R': bitsA | bitsG,
	'Y': bitsC | bitsT,
	'S': bitsG | bitsC,
	'W': bitsA | bitsT,
	'K': bitsG | bitsT,
	'M': bitsA | bitsC,
	'B': bitsC | bitsG | bitsT,
	'D': bitsA | bitsG | bitsT,
	'H': bitsA | bitsC | bitsT,
	'V': bitsA | bitsC | bitsG,
	'N': bitsA | bitsC | bitsG | bitsT,
}

// complementBits swaps A<->T and C<->G within a base-set bitmask, preserving
// the IUPAC semantics of the complement (e.g. R<->Y, S<->S, W<->W, K<->M,
// B<->V, D<->H, N<->N).
func complementBits(s uint8) uint8 {
	var out uint8
	if s&bitsA != 0 {
		out |= bitsT
	}
	if s&bitsT != 0 {
		out |= bitsA
	}
	if s&bitsC != 0 {
		out |= bitsG
	}
	if s&bitsG != 0 {
		out |= bitsC
	}
	return out
}

// IsValidBase reports whether b is a recognized single IUPAC nucleotide code.
func IsValidBase(b byte) bool {
	_, ok := iupacSet[b]
	return ok
}

// SetFor returns the canonical-base bitmask for an IUPAC code and whether the
// code is recognized. It is the exported accessor used by the motif and
// restriction packages to compile patterns without re-implementing the table.
func SetFor(b byte) (uint8, bool) {
	s, ok := iupacSet[upperByte(b)]
	return s, ok
}

// IsValidSequence reports whether s contains only valid IUPAC codes. It is
// case-insensitive (uppercased before lookup) and rejects empty input only when
// the caller asks for it.
func IsValidSequence(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		if !IsValidBase(upperByte(s[i])) {
			return false
		}
	}
	return true
}

// Normalize upper-cases every byte of s in place (ASCII only; IUPAC letters are
// all ASCII so this is sufficient and allocation-free).
func Normalize(s string) string {
	b := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		b[i] = upperByte(s[i])
	}
	return string(b)
}

func upperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - ('a' - 'A')
	}
	return b
}

// BasesMatch reports whether the sequence base sb can match the query code qb
// under IUPAC degeneracy rules. Both operands may be degenerate; a match
// exists if the two base-sets share at least one canonical nucleotide. Either
// byte being an unknown code means no match.
func BasesMatch(sb, qb byte) bool {
	ss, ok1 := iupacSet[upperByte(sb)]
	qs, ok2 := iupacSet[upperByte(qb)]
	if !ok1 || !ok2 {
		return false
	}
	return ss&qs != 0
}

// Complement returns the IUPAC complement of a single base. Returns 0 for an
// unknown code (callers must validate first).
func Complement(b byte) byte {
	s, ok := iupacSet[upperByte(b)]
	if !ok {
		return 0
	}
	return setToCode(complementBits(s))
}

// setToCode maps a base-set bitmask back to the canonical IUPAC letter. When
// several letters share the same set (there are none for the standard table),
// the first registered wins; the table is bijective by construction.
func setToCode(s uint8) byte {
	for c, set := range iupacSet {
		if set == s {
			return c
		}
	}
	return 'N'
}

// ReverseComplement returns the reverse complement of s. It complements each
// base under IUPAC rules and reverses the order, so the result reads 3'->5'
// on the original strand. Unknown codes are rejected upstream.
func ReverseComplement(s string) string {
	n := len(s)
	out := make([]byte, n)
	for i := 0; i < n; i++ {
		out[n-1-i] = Complement(s[i])
	}
	return string(out)
}

// WrapIndex maps a possibly out-of-range position to a valid index on a
// sequence of length n by taking it modulo n. It is the coordinate primitive
// used by every circular-aware feature (ORFs, motifs, restriction sites,
// annotations). n must be > 0.
func WrapIndex(pos, n int) int {
	pos %= n
	if pos < 0 {
		pos += n
	}
	return pos
}

// WrapLen reports how many distinct positions a half-open interval [start, end)
// spans on a circular sequence of length n, accounting for wraparound. start
// and end are absolute (possibly >= n or < 0) positions; the result is the
// forward distance modulo n in [0, n].
func WrapLen(start, end, n int) int {
	if n <= 0 {
		return 0
	}
	d := (end - start) % n
	if d < 0 {
		d += n
	}
	return d
}

// NormalizeRange clamps a [start, end) half-open interval onto a linear
// sequence of length n. It returns the 1-based inclusive [lo, hi] pair used
// throughout the public API, or false if the interval is empty or fully out
// of range. For circular sequences callers should use WrapIndex instead.
func NormalizeRange(start, end, n int) (int, int, bool) {
	if start < 0 || end < 0 || start >= end || start >= n {
		return 0, 0, false
	}
	if end > n {
		end = n
	}
	return start + 1, end, true
}

// joinStrands concatenates two strings; a tiny helper kept in this package so
// motif/restriction callers share one implementation.
func joinStrands(a, b string) string {
	var sb strings.Builder
	sb.WriteString(a)
	sb.WriteString(b)
	return sb.String()
}
