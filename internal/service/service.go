package service

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"time"

	"task136-seqannot/internal/align"
	"task136-seqannot/internal/annot"
	"task136-seqannot/internal/motif"
	"task136-seqannot/internal/orf"
	"task136-seqannot/internal/restrict"
	"task136-seqannot/internal/seq"
	"task136-seqannot/internal/store"
)

// Service is the orchestration layer between the HTTP handlers and the domain
// packages + store. It owns no mutable in-process state beyond the store, so
// it is safe to share across goroutines; all shared state lives in SQLite.
type Service struct {
	store *store.Store
}

// New constructs a Service over an opened store.
func New(s *store.Store) *Service {
	return &Service{store: s}
}

// SequenceRequest is the body for creating a sequence.
type SequenceRequest struct {
	Name        string `json:"name"`
	Residues    string `json:"residues"`
	Type        string `json:"type"`
	Description string `json:"description"`
	FastaHeader string `json:"fasta_header"`
}

// CreateSequence validates and persists a new sequence, assigning an id.
func (svc *Service) CreateSequence(ctx context.Context, req SequenceRequest) (store.SeqRow, error) {
	typ := seq.SequenceType(req.Type)
	s, err := seq.NewSequence(req.Name, req.Residues, typ, req.Description, req.FastaHeader)
	if err != nil {
		return store.SeqRow{}, err
	}
	id := newID("seq")
	row := store.SeqRow{
		ID: id, Name: s.Name, Type: s.Type, Residues: s.Residues,
		Description: s.Description, FastaHeader: s.FastaHeader, Length: s.Length(),
		CreatedAt: stamp(),
	}
	if err := svc.store.SaveSequence(ctx, row); err != nil {
		return store.SeqRow{}, err
	}
	return row, nil
}

// ImportFasta parses FASTA text and persists the first record (the engine
// treats the first record as the canonical one; multi-record import is not
// supported to keep ids deterministic per request).
func (svc *Service) ImportFasta(ctx context.Context, fastaText string) (store.SeqRow, error) {
	records, err := seq.ParseFasta(fastaText)
	if err != nil {
		return store.SeqRow{}, err
	}
	r := records[0]
	s, err := seq.NewSequence(r.ID, r.Residues, seq.TypeLinear, r.Desc, r.Header)
	if err != nil {
		return store.SeqRow{}, err
	}
	row := store.SeqRow{
		ID: newID("seq"), Name: s.Name, Type: s.Type, Residues: s.Residues,
		Description: s.Description, FastaHeader: s.FastaHeader, Length: s.Length(),
		CreatedAt: stamp(),
	}
	if err := svc.store.SaveSequence(ctx, row); err != nil {
		return store.SeqRow{}, err
	}
	return row, nil
}

// loadSeq rebuilds a seq.Sequence from a stored row, validating on the way out
// so a corrupted row surfaces as an error rather than a silent zero value.
func loadSeq(r store.SeqRow) (seq.Sequence, error) {
	s, err := seq.NewSequence(r.Name, r.Residues, r.Type, r.Description, r.FastaHeader)
	if err != nil {
		return seq.Sequence{}, err
	}
	s.ID = r.ID
	return s, nil
}

// Composition returns GC content and molecular weight for a stored sequence.
func (svc *Service) Composition(ctx context.Context, id string) (seq.Composition, error) {
	r, err := svc.store.GetSequence(ctx, id)
	if err != nil {
		return seq.Composition{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return seq.Composition{}, err
	}
	return seq.Compose(s), nil
}

// ReverseComplement returns the reverse complement string of a stored sequence.
func (svc *Service) ReverseComplement(ctx context.Context, id string) (string, error) {
	r, err := svc.store.GetSequence(ctx, id)
	if err != nil {
		return "", err
	}
	s, err := loadSeq(r)
	if err != nil {
		return "", err
	}
	return seq.ReverseComplement(s.Residues), nil
}

// TranslateResult is the JSON shape returned by the translate API.
type TranslateResult struct {
	Frames []seq.FrameInfo `json:"frames"`
}

// Translate computes the six reading frames for a stored sequence.
func (svc *Service) Translate(ctx context.Context, id string) (TranslateResult, error) {
	r, err := svc.store.GetSequence(ctx, id)
	if err != nil {
		return TranslateResult{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return TranslateResult{}, err
	}
	return TranslateResult{Frames: seq.Translate6(s)}, nil
}

// ORFResult is the JSON shape for ORF detection.
type ORFResult struct {
	ORFs []orf.ORF `json:"orfs"`
}

// FindORFs detects ORFs in a stored sequence with a minimum amino-acid length.
func (svc *Service) FindORFs(ctx context.Context, id string, minAALen int) (ORFResult, error) {
	r, err := svc.store.GetSequence(ctx, id)
	if err != nil {
		return ORFResult{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return ORFResult{}, err
	}
	return ORFResult{ORFs: orf.Find(s, orf.Options{MinAALen: minAALen})}, nil
}

// MotifSearchResult is the JSON shape for motif hits.
type MotifSearchResult struct {
	Hits []motif.Hit `json:"hits"`
}

// SearchMotif searches a stored sequence for a motif identified by either an
// inline pattern or a persisted motif id.
func (svc *Service) SearchMotif(ctx context.Context, seqID, pattern, motifID string) (MotifSearchResult, error) {
	r, err := svc.store.GetSequence(ctx, seqID)
	if err != nil {
		return MotifSearchResult{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return MotifSearchResult{}, err
	}
	var m motif.Motif
	if motifID != "" {
		mr, err := svc.store.GetMotif(ctx, motifID)
		if err != nil {
			return MotifSearchResult{}, err
		}
		m, err = motif.Compile(mr.ID, mr.Name, mr.Pattern, mr.Description)
		if err != nil {
			return MotifSearchResult{}, err
		}
	} else {
		m, err = motif.Compile(newID("motif"), "ad-hoc", pattern, "")
		if err != nil {
			return MotifSearchResult{}, err
		}
	}
	return MotifSearchResult{Hits: motif.Search(s, m)}, nil
}

// RestrictionMap finds cut sites for a set of enzymes (by id or inline).
func (svc *Service) RestrictionMap(ctx context.Context, seqID string, enzymeIDs []string) (restrictSitesResult, error) {
	r, err := svc.store.GetSequence(ctx, seqID)
	if err != nil {
		return restrictSitesResult{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return restrictSitesResult{}, err
	}
	enzymes, err := svc.resolveEnzymes(ctx, enzymeIDs)
	if err != nil {
		return restrictSitesResult{}, err
	}
	var sites []restrict.CutSite
	for _, e := range enzymes {
		sites = append(sites, restrict.FindSites(s, e)...)
	}
	return restrictSitesResult{Sites: sites}, nil
}

// Digest simulates restriction digestion and returns fragment lengths.
func (svc *Service) Digest(ctx context.Context, seqID string, enzymeIDs []string) (digestResult, error) {
	r, err := svc.store.GetSequence(ctx, seqID)
	if err != nil {
		return digestResult{}, err
	}
	s, err := loadSeq(r)
	if err != nil {
		return digestResult{}, err
	}
	enzymes, err := svc.resolveEnzymes(ctx, enzymeIDs)
	if err != nil {
		return digestResult{}, err
	}
	return digestResult{Fragments: restrict.Digest(s, enzymes)}, nil
}

type restrictSitesResult struct {
	Sites []restrict.CutSite `json:"sites"`
}

type digestResult struct {
	Fragments []int `json:"fragments"`
}

func (svc *Service) resolveEnzymes(ctx context.Context, ids []string) ([]restrict.Enzyme, error) {
	var out []restrict.Enzyme
	for _, id := range ids {
		er, err := svc.store.GetEnzyme(ctx, id)
		if err != nil {
			return nil, err
		}
		e, err := restrict.Compile(er.ID, er.Name, er.Site, er.CutOffset)
		if err != nil {
			return nil, err
		}
		out = append(out, e)
	}
	return out, nil
}

// AlignmentRequest configures a pairwise alignment.
type AlignmentRequest struct {
	SeqAID  string `json:"seq_a_id"`
	SeqBID  string `json:"seq_b_id"`
	Mode    string `json:"mode"` // "global" or "local"
	Match   int    `json:"match"`
	Mismatch int   `json:"mismatch"`
	Gap     int    `json:"gap"`
}

// Align two stored sequences and persist the result.
func (svc *Service) Align(ctx context.Context, req AlignmentRequest) (store.AlignmentRow, error) {
	a, err := svc.store.GetSequence(ctx, req.SeqAID)
	if err != nil {
		return store.AlignmentRow{}, err
	}
	b, err := svc.store.GetSequence(ctx, req.SeqBID)
	if err != nil {
		return store.AlignmentRow{}, err
	}
	params := align.Params{Match: req.Match, Mismatch: req.Mismatch, Gap: req.Gap}
	var res align.Result
	switch req.Mode {
	case "global", "":
		res = align.Global(a.Residues, b.Residues, params)
	case "local":
		res = align.Local(a.Residues, b.Residues, params)
	default:
		return store.AlignmentRow{}, fmt.Errorf("service: bad align mode %q", req.Mode)
	}
	paramsJSON, _ := json.Marshal(params)
	row := store.AlignmentRow{
		ID: newID("aln"), SeqAID: req.SeqAID, SeqBID: req.SeqBID, Mode: res.Mode,
		Score: res.Score, AlignedA: res.AlignedA, AlignedB: res.AlignedB,
		Params: string(paramsJSON), CreatedAt: stamp(),
	}
	if err := svc.store.SaveAlignment(ctx, row); err != nil {
		return store.AlignmentRow{}, err
	}
	return row, nil
}

// Motif CRUD helpers ----------------------------------------------------------

// CreateMotif persists a compiled motif after validation.
func (svc *Service) CreateMotif(ctx context.Context, name, pattern, desc string) (store.MotifRow, error) {
	if _, err := motif.Compile(newID("m"), name, pattern, desc); err != nil {
		return store.MotifRow{}, err
	}
	row := store.MotifRow{ID: newID("motif"), Name: name, Pattern: seq.Normalize(pattern), Description: desc, CreatedAt: stamp()}
	if err := svc.store.SaveMotif(ctx, row); err != nil {
		return store.MotifRow{}, err
	}
	return row, nil
}

// CreateEnzyme persists a restriction enzyme after validation.
func (svc *Service) CreateEnzyme(ctx context.Context, name, site string, cut int) (store.EnzymeRow, error) {
	if _, err := restrict.Compile(newID("e"), name, site, cut); err != nil {
		return store.EnzymeRow{}, err
	}
	row := store.EnzymeRow{ID: newID("enzyme"), Name: name, Site: seq.Normalize(site), CutOffset: cut, CreatedAt: stamp()}
	if err := svc.store.SaveEnzyme(ctx, row); err != nil {
		return store.EnzymeRow{}, err
	}
	return row, nil
}

// Feature CRUD --------------------------------------------------------------

// CreateFeature validates and persists an annotation feature.
func (svc *Service) CreateFeature(ctx context.Context, f store.FeatureRow) (store.FeatureRow, error) {
	r, err := svc.store.GetSequence(ctx, f.SequenceID)
	if err != nil {
		return store.FeatureRow{}, err
	}
	feat := toFeature(f)
	if err := feat.Validate(r.Length, r.Type == seq.TypeCircular); err != nil {
		return store.FeatureRow{}, err
	}
	f.ID = newID("feat")
	f.CreatedAt = stamp()
	if f.Strand == "" {
		f.Strand = "+"
	}
	if err := svc.store.SaveFeature(ctx, f); err != nil {
		return store.FeatureRow{}, err
	}
	return f, nil
}

// QueryFeatures returns features on a sequence, optionally filtered to those
// overlapping a [start,end] interval.
func (svc *Service) QueryFeatures(ctx context.Context, seqID string, start, end int) ([]store.FeatureRow, error) {
	rows, err := svc.store.ListFeatures(ctx, seqID)
	if err != nil {
		return nil, err
	}
	r, err := svc.store.GetSequence(ctx, seqID)
	if err != nil {
		return nil, err
	}
	if start == 0 && end == 0 {
		return rows, nil
	}
	feats := make([]annot.Feature, 0, len(rows))
	for _, fr := range rows {
		feats = append(feats, toFeature(fr))
	}
	matched := annot.QueryOverlap(feats, r.Length, r.Type == seq.TypeCircular, start, end)
	out := make([]store.FeatureRow, 0, len(matched))
	for _, mf := range matched {
		for _, fr := range rows {
			if fr.ID == mf.ID {
				out = append(out, fr)
			}
		}
	}
	return out, nil
}

// DeleteFeature removes a feature by id.
func (svc *Service) DeleteFeature(ctx context.Context, id string) error {
	return svc.store.DeleteFeature(ctx, id)
}

// errInvalid is returned for caller-validation failures that should map to 400.
var errInvalid = errors.New("invalid request")

// stamp returns an RFC3339Nano UTC timestamp. It is the single source of
// wall-clock time in the engine and is only used for created_at/updated_at
// bookkeeping; analysis results never depend on it, preserving idempotency.
func stamp() string {
	return time.Now().UTC().Format(time.RFC3339Nano)
}
