package align

import "fmt"

// SubstitutionMatrix scores the alignment of two amino acids (or nucleotides,
// represented as bytes). The default nucleotide matrix scores +match for an
// exact base identity and -mismatch otherwise; ambiguity is treated as a
// mismatch unless the two codes share a canonical base.
type SubstitutionMatrix struct {
	Score func(a, b byte) int
}

// NucleotideMatrix returns a simple match/mismatch scoring function for DNA
// bases. match is the reward for identical bases; mismatch is the penalty for
// a difference. IUPAC-degenerate bases match if their base-sets intersect.
func NucleotideMatrix(match, mismatch int) SubstitutionMatrix {
	return SubstitutionMatrix{Score: func(a, b byte) int {
		if baseIntersect(a, b) {
			return match
		}
		return mismatch
	}}
}

// baseIntersect reports whether two IUPAC codes share at least one canonical
// nucleotide. It mirrors seq.BasesMatch but the align package is deliberately
// independent of the seq package to keep the DP core dependency-free.
func baseIntersect(a, b byte) bool {
	return setOf(a)&setOf(b) != 0
}

func setOf(b byte) uint8 {
	switch b {
	case 'A':
		return 1
	case 'C':
		return 2
	case 'G':
		return 4
	case 'T':
		return 8
	case 'U':
		return 8 // treat U as T
	case 'R':
		return 1 | 4
	case 'Y':
		return 2 | 8
	case 'S':
		return 4 | 2
	case 'W':
		return 1 | 8
	case 'K':
		return 4 | 8
	case 'M':
		return 1 | 2
	case 'B':
		return 2 | 4 | 8
	case 'D':
		return 1 | 4 | 8
	case 'H':
		return 1 | 2 | 8
	case 'V':
		return 1 | 2 | 4
	case 'N':
		return 1 | 2 | 4 | 8
	}
	return 0
}

// Result is an alignment of two sequences. The AlignedA/AlignedB strings are
// the same length and use '-' to denote a gap. Score is the total alignment
// score. Mode is "global" (Needleman-Wunsch) or "local" (Smith-Waterman).
type Result struct {
	Mode     string `json:"mode"`
	Score    int    `json:"score"`
	AlignedA string `json:"aligned_a"`
	AlignedB string `json:"aligned_b"`
}

// Params configures an alignment.
type Params struct {
	Match    int
	Mismatch int
	Gap      int // gap penalty (applied per gap position); must be <= 0
}

// Global aligns two sequences end-to-end with the Needleman-Wunsch algorithm
// using affine-ish linear gap penalties (one penalty per gap position). The
// traceback is reconstructed with a tie-break that prefers diagonal (match)
// over up/left, making the alignment a deterministic pure function of the two
// input strings and params.
func Global(a, b string, p Params) Result {
	// pad gap penalty to be non-positive
	gap := p.Gap
	if gap > 0 {
		gap = -gap
	}
	mtx := NucleotideMatrix(p.Match, p.Mismatch)
	score, sa, sb := dp(a, b, mtx.Score, gap, true)
	return Result{Mode: "global", Score: score, AlignedA: sa, AlignedB: sb}
}

// Local aligns two sequences with Smith-Waterman, returning the highest-scoring
// local alignment. The traceback starts from the single highest-scoring cell;
// ties are broken by (row, then column) so the result is deterministic.
func Local(a, b string, p Params) Result {
	gap := p.Gap
	if gap > 0 {
		gap = -gap
	}
	mtx := NucleotideMatrix(p.Match, p.Mismatch)
	score, sa, sb := dp(a, b, mtx.Score, gap, false)
	return Result{Mode: "local", Score: score, AlignedA: sa, AlignedB: sb}
}

// dp runs the dynamic programming. For global it uses Needleman-Wunsch with
// full traceback from (m,n); for local it uses Smith-Waterman with zero-floor
// and traceback from the max cell until a zero is reached. The tie-break order
// in fill is diag > up > left, and the traceback follows the recorded pointer.
func dp(a, b string, score func(byte, byte) int, gap int, global bool) (int, string, string) {
	m, n := len(a), len(b)
	if m == 0 && n == 0 {
		return 0, "", ""
	}
	// dp[i][j] = best score aligning a[:i] with b[:j]
	dp := make([][]int, m+1)
	ptr := make([][]byte, m+1)
	for i := range dp {
		dp[i] = make([]int, n+1)
		ptr[i] = make([]byte, n+1)
	}
	if global {
		for i := 1; i <= m; i++ {
			dp[i][0] = dp[i-1][0] + gap
			ptr[i][0] = 'u'
		}
		for j := 1; j <= n; j++ {
			dp[0][j] = dp[0][j-1] + gap
			ptr[0][j] = 'l'
		}
	}
	maxI, maxJ, maxVal := 0, 0, 0
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			diag := dp[i-1][j-1] + score(a[i-1], b[j-1])
			up := dp[i-1][j] + gap
			left := dp[i][j-1] + gap
			best := diag
			dir := byte('d')
			if up > best {
				best = up
				dir = 'u'
			}
			if left > best {
				best = left
				dir = 'l'
			}
			if !global && best < 0 {
				best = 0
				dir = 0
			}
			dp[i][j] = best
			ptr[i][j] = dir
			if !global && best > maxVal {
				maxVal = best
				maxI, maxJ = i, j
			}
		}
	}
	if global {
		maxVal = dp[m][n]
		maxI, maxJ = m, n
	}
	// traceback
	var ra, rb []byte
	i, j := maxI, maxJ
	for {
		if global {
			if i == 0 && j == 0 {
				break
			}
		} else {
			if dp[i][j] == 0 {
				break
			}
		}
		d := ptr[i][j]
		switch {
		case d == 'd' || (i > 0 && j > 0 && (d == 0 && false)):
			ra = append(ra, a[i-1])
			rb = append(rb, b[j-1])
			i--
			j--
		case d == 'u':
			ra = append(ra, a[i-1])
			rb = append(rb, '-')
			i--
		case d == 'l':
			ra = append(ra, '-')
			rb = append(rb, b[j-1])
			j--
		default:
			// for local, d==0 with dp==0 already broke; for global we only
			// reach here on the boundaries which are handled by the loop above
			if i > 0 && j == 0 {
				ra = append(ra, a[i-1])
				rb = append(rb, '-')
				i--
			} else if j > 0 && i == 0 {
				ra = append(ra, '-')
				rb = append(rb, b[j-1])
				j--
			} else {
				break
			}
		}
	}
	reverse(ra)
	reverse(rb)
	return maxVal, string(ra), string(rb)
}

func reverse(xs []byte) {
	for i, j := 0, len(xs)-1; i < j; i, j = i+1, j-1 {
		xs[i], xs[j] = xs[j], xs[i]
	}
}

// Describe returns a compact human-readable summary of an alignment result,
// used by the API layer.
func (r Result) Describe() string {
	return fmt.Sprintf("%s score=%d\n%s\n%s", r.Mode, r.Score, r.AlignedA, r.AlignedB)
}
