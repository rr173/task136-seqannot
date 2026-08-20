package seq

// Standard genetic code: 64 sense codons -> single-letter amino acid.
// Start codons (the canonical ATG and the alternative near-cognates that the
// engine reports as "start-capable") are tracked separately so ORF detection
// can mark the initiating methionine. Stop codons map to '*'.

// Amino is the one-letter amino acid code for a codon; '*' denotes a stop.
type Amino byte

const (
	Stop      Amino = '*'
	Metionine Amino = 'M' // canonical start residue
)

// codonTable maps a 3-letter uppercase codon (A/C/G/T only; IUPAC ambiguity is
// resolved by CodonAmino with set expansion) to its amino acid under the
// standard bacterial/archaeal/plastid "translation table 1" (NCBI).
var codonTable = map[[3]byte]Amino{
	{'T', 'T', 'T'}: 'F', {'T', 'T', 'C'}: 'F', {'T', 'T', 'A'}: 'L', {'T', 'T', 'G'}: 'L',
	{'C', 'T', 'T'}: 'L', {'C', 'T', 'C'}: 'L', {'C', 'T', 'A'}: 'L', {'C', 'T', 'G'}: 'L',
	{'A', 'T', 'T'}: 'I', {'A', 'T', 'C'}: 'I', {'A', 'T', 'A'}: 'I', {'A', 'T', 'G'}: 'M',
	{'G', 'T', 'T'}: 'V', {'G', 'T', 'C'}: 'V', {'G', 'T', 'A'}: 'V', {'G', 'T', 'G'}: 'V',
	{'T', 'C', 'T'}: 'S', {'T', 'C', 'C'}: 'S', {'T', 'C', 'A'}: 'S', {'T', 'C', 'G'}: 'S',
	{'C', 'C', 'T'}: 'P', {'C', 'C', 'C'}: 'P', {'C', 'C', 'A'}: 'P', {'C', 'C', 'G'}: 'P',
	{'A', 'C', 'T'}: 'T', {'A', 'C', 'C'}: 'T', {'A', 'C', 'A'}: 'T', {'A', 'C', 'G'}: 'T',
	{'G', 'C', 'T'}: 'A', {'G', 'C', 'C'}: 'A', {'G', 'C', 'A'}: 'A', {'G', 'C', 'G'}: 'A',
	{'T', 'A', 'T'}: 'Y', {'T', 'A', 'C'}: 'Y', {'T', 'A', 'A'}: Stop, {'T', 'A', 'G'}: Stop,
	{'C', 'A', 'T'}: 'H', {'C', 'A', 'C'}: 'H', {'C', 'A', 'A'}: 'Q', {'C', 'A', 'G'}: 'Q',
	{'A', 'A', 'T'}: 'N', {'A', 'A', 'C'}: 'N', {'A', 'A', 'A'}: 'K', {'A', 'A', 'G'}: 'K',
	{'G', 'A', 'T'}: 'D', {'G', 'A', 'C'}: 'D', {'G', 'A', 'A'}: 'E', {'G', 'A', 'G'}: 'E',
	{'T', 'G', 'T'}: 'C', {'T', 'G', 'C'}: 'C', {'T', 'G', 'A'}: Stop, {'T', 'G', 'G'}: 'W',
	{'C', 'G', 'T'}: 'R', {'C', 'G', 'C'}: 'R', {'C', 'G', 'A'}: 'R', {'C', 'G', 'G'}: 'R',
	{'A', 'G', 'T'}: 'S', {'A', 'G', 'C'}: 'S', {'A', 'G', 'A'}: 'R', {'A', 'G', 'G'}: 'R',
	{'G', 'G', 'T'}: 'G', {'G', 'G', 'C'}: 'G', {'G', 'G', 'A'}: 'G', {'G', 'G', 'G'}: 'G',
}

// startCodons are the codons treated as translation starts. ATG is canonical;
// the near-cognates are NOT treated as starts by this engine (only ATG
// initiates an ORF), which keeps ORF detection deterministic and matches the
// common bioinformatic convention for prokaryotic-style scanning.
var startCodons = map[[3]byte]bool{
	{'A', 'T', 'G'}: true,
}

// IsStartCodon reports whether the 3 canonical bases form an initiator.
func IsStartCodon(c [3]byte) bool { return startCodons[c] }

// IsStopCodon reports whether the 3 canonical bases terminate translation.
func IsStopCodon(c [3]byte) bool {
	a, ok := codonTable[c]
	return ok && a == Stop
}

// CodonAmino resolves an amino acid for a codon that may contain IUPAC
// ambiguity. If every possible expansion of the codon yields the same amino
// acid, that amino is returned with unambiguous=true. If expansions disagree
// (e.g. the codon spans a stop and a sense) the codon is reported as
// ambiguous ('X') with unambiguous=false. Purely unknown bases yield 'X'.
// Codons whose length is not 3 are treated as 'X' and unambiguous=false.
func CodonAmino(c0, c1, c2 byte) (Amino, bool) {
	if !IsValidBase(c0) || !IsValidBase(c1) || !IsValidBase(c2) {
		return 'X', false
	}
	expansions := expandCodon([3]byte{upperByte(c0), upperByte(c1), upperByte(c2)})
	if len(expansions) == 0 {
		return 'X', false
	}
	first, ok := codonTable[expansions[0]]
	if !ok {
		return 'X', false
	}
	for _, e := range expansions[1:] {
		a, ok := codonTable[e]
		if !ok || a != first {
			return 'X', false
		}
	}
	return first, true
}

// expandCodon enumerates every canonical (A/C/G/T) codon consistent with an
// IUPAC-degenerate triplet. The order is lexicographic by base (A<C<G<T) so
// the expansion is deterministic.
func expandCodon(c [3]byte) [][3]byte {
	bases := make([][]byte, 3)
	for i, b := range c {
		bases[i] = codeToCanonical(b)
	}
	var out [][3]byte
	var rec func(int, [3]byte)
	rec = func(i int, cur [3]byte) {
		if i == 3 {
			cp := cur
			out = append(out, cp)
			return
		}
		for _, b := range bases[i] {
			cur[i] = b
			rec(i+1, cur)
		}
	}
	rec(0, [3]byte{})
	return out
}

// codeToCanonical returns the canonical bases a single IUPAC code can stand
// for, in fixed A<C<G<T order.
func codeToCanonical(b byte) []byte {
	s, ok := iupacSet[upperByte(b)]
	if !ok {
		return nil
	}
	out := make([]byte, 0, 4)
	if s&bitsA != 0 {
		out = append(out, 'A')
	}
	if s&bitsC != 0 {
		out = append(out, 'C')
	}
	if s&bitsG != 0 {
		out = append(out, 'G')
	}
	if s&bitsT != 0 {
		out = append(out, 'T')
	}
	return out
}
