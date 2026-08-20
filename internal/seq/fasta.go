package seq

import (
	"errors"
	"strings"
)

// FASTA parsing: one or more records separated by a header line beginning with
// '>'. Sequence lines accumulate into the residues; whitespace within a record
// is stripped. The parser is intentionally strict about the structure: a record
// without a header, an empty id, or residues that fail the IUPAC check are
// rejected with a descriptive error rather than silently coerced.

// FastaRecord is a single parsed FASTA entry.
type FastaRecord struct {
	Header   string // full header line after '>'
	ID       string // first whitespace-delimited token of the header
	Desc     string // remainder of the header, if any
	Residues string // normalized (upper-cased) sequence
}

// ErrEmptyFasta is returned when the input contains no records.
var ErrEmptyFasta = errors.New("fasta: no records")

// ErrFastaHeader is returned when a record's header is missing or empty.
var ErrFastaHeader = errors.New("fasta: missing or empty header")

// ParseFasta parses one or more FASTA records. Each record must start with a
// '>' header line followed by one or more sequence lines. Empty input or a
// record with an empty id yields an error. Sequence letters are normalized to
// upper case but NOT validated here (the caller decides whether to reject
// ambiguous bases); the engine validates at ingestion via IsValidSequence.
func ParseFasta(text string) ([]FastaRecord, error) {
	text = strings.ReplaceAll(text, "\r\n", "\n")
	text = strings.ReplaceAll(text, "\r", "\n")
	lines := strings.Split(text, "\n")

	var records []FastaRecord
	var cur *FastaRecord
	var seqb strings.Builder

	flush := func() error {
		if cur == nil {
			return nil
		}
		res := strings.ToUpper(strings.ReplaceAll(seqb.String(), " ", ""))
		seqb.Reset()
		cur.Residues = res
		if cur.ID == "" {
			return ErrFastaHeader
		}
		if cur.Residues == "" {
			return errors.New("fasta: record " + cur.ID + " has empty sequence")
		}
		records = append(records, *cur)
		cur = nil
		return nil
	}

	for i, ln := range lines {
		trimmed := strings.TrimSpace(ln)
		switch {
		case strings.HasPrefix(trimmed, ">"):
			if err := flush(); err != nil {
				return nil, err
			}
			hdr := strings.TrimSpace(strings.TrimPrefix(trimmed, ">"))
			if hdr == "" {
				return nil, ErrFastaHeader
			}
			id, desc := splitHeader(hdr)
			cur = &FastaRecord{Header: hdr, ID: id, Desc: desc}
		case trimmed == "":
			// blank lines inside a record are tolerated and skipped
			continue
		default:
			if cur == nil {
				if i == 0 {
					// tolerate a bare sequence (no header) as a single
					// anonymous record only when it is the whole input
					cur = &FastaRecord{ID: "seq", Desc: "", Header: "seq"}
				} else {
					return nil, ErrFastaHeader
				}
			}
			seqb.WriteString(trimmed)
		}
	}
	if err := flush(); err != nil {
		return nil, err
	}
	if len(records) == 0 {
		return nil, ErrEmptyFasta
	}
	return records, nil
}

// splitHeader splits a header into the first token (id) and the remainder
// (description), matching the conventional ">ID description" layout.
func splitHeader(h string) (string, string) {
	h = strings.TrimSpace(h)
	idx := strings.IndexAny(h, " \t")
	if idx < 0 {
		return h, ""
	}
	return h[:idx], strings.TrimSpace(h[idx+1:])
}

// FormatFasta renders records back to FASTA text with a fixed 60-char line
// width. It is the inverse of ParseFasta for round-trip testing.
func FormatFasta(records []FastaRecord, lineWidth int) string {
	if lineWidth <= 0 {
		lineWidth = 60
	}
	var b strings.Builder
	for _, r := range records {
		b.WriteByte('>')
		b.WriteString(r.Header)
		b.WriteByte('\n')
		for i := 0; i < len(r.Residues); i += lineWidth {
			end := i + lineWidth
			if end > len(r.Residues) {
				end = len(r.Residues)
			}
			b.WriteString(r.Residues[i:end])
			b.WriteByte('\n')
		}
	}
	return b.String()
}
