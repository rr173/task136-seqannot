package store

import (
	"context"
	"time"

	"task136-seqannot/internal/seq"
)

// SeqRow is the persisted form of a sequence.
type SeqRow struct {
	ID          string             `json:"id"`
	Name        string             `json:"name"`
	Type        seq.SequenceType   `json:"type"`
	Residues    string             `json:"residues"`
	Description string             `json:"description"`
	FastaHeader string             `json:"fasta_header"`
	Length      int                `json:"length"`
	CreatedAt   string             `json:"created_at"`
}

// SaveSequence persists a sequence. It is idempotent on the primary key:
// re-inserting the same id replaces the row (useful for the replay path).
func (s *Store) SaveSequence(ctx context.Context, r SeqRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO sequences(id,name,type,residues,description,fasta_header,length,created_at)
			VALUES(?,?,?,?,?,?,?,?)
			ON CONFLICT(id) DO UPDATE SET name=excluded.name,type=excluded.type,residues=excluded.residues,
				description=excluded.description,fasta_header=excluded.fasta_header,length=excluded.length`,
			r.ID, r.Name, string(r.Type), r.Residues, r.Description, r.FastaHeader, r.Length, r.CreatedAt)
		return err
	})
}

// GetSequence loads a sequence by id.
func (s *Store) GetSequence(ctx context.Context, id string) (SeqRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx,
		`SELECT id,name,type,residues,description,fasta_header,length,created_at FROM sequences WHERE id=?`, id)
	var r SeqRow
	var typ string
	if err := row.Scan(&r.ID, &r.Name, &typ, &r.Residues, &r.Description, &r.FastaHeader, &r.Length, &r.CreatedAt); err != nil {
		if isNoRows(err) {
			return SeqRow{}, ErrNotFound
		}
		return SeqRow{}, err
	}
	r.Type = seq.SequenceType(typ)
	return r, nil
}

// ListSequences returns all sequences ordered by id.
func (s *Store) ListSequences(ctx context.Context) ([]SeqRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx,
		`SELECT id,name,type,residues,description,fasta_header,length,created_at FROM sequences ORDER BY id`)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []SeqRow
	for rows.Next() {
		var r SeqRow
		var typ string
		if err := rows.Scan(&r.ID, &r.Name, &typ, &r.Residues, &r.Description, &r.FastaHeader, &r.Length, &r.CreatedAt); err != nil {
			return nil, err
		}
		r.Type = seq.SequenceType(typ)
		out = append(out, r)
	}
	return out, rows.Err()
}

// DeleteSequence removes a sequence (cascade removes its jobs/features).
func (s *Store) DeleteSequence(ctx context.Context, id string) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `DELETE FROM sequences WHERE id=?`, id)
		return err
	})
}

// JobRow is the persisted form of an analysis job.
type JobRow struct {
	ID         string `json:"id"`
	SequenceID string `json:"sequence_id"`
	Params     string `json:"params"`
	Status     string `json:"status"`
	Progress   int    `json:"progress"`
	Results    string `json:"results"`
	Error      string `json:"error"`
	CreatedAt  string `json:"created_at"`
	UpdatedAt  string `json:"updated_at"`
}

// SaveJob upserts an analysis job row.
func (s *Store) SaveJob(ctx context.Context, r JobRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO analysis_jobs(id,sequence_id,params,status,progress,results,error,created_at,updated_at)
			VALUES(?,?,?,?,?,?,?,?,?)
			ON CONFLICT(id) DO UPDATE SET sequence_id=excluded.sequence_id,params=excluded.params,status=excluded.status,
				progress=excluded.progress,results=excluded.results,error=excluded.error,updated_at=excluded.updated_at`,
			r.ID, r.SequenceID, r.Params, r.Status, r.Progress, r.Results, r.Error, r.CreatedAt, r.UpdatedAt)
		return err
	})
}

// GetJob loads a job by id.
func (s *Store) GetJob(ctx context.Context, id string) (JobRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx,
		`SELECT id,sequence_id,params,status,progress,results,error,created_at,updated_at FROM analysis_jobs WHERE id=?`, id)
	var r JobRow
	if err := row.Scan(&r.ID, &r.SequenceID, &r.Params, &r.Status, &r.Progress, &r.Results, &r.Error, &r.CreatedAt, &r.UpdatedAt); err != nil {
		if isNoRows(err) {
			return JobRow{}, ErrNotFound
		}
		return JobRow{}, err
	}
	return r, nil
}

// ListJobs returns all jobs ordered by created_at then id.
func (s *Store) ListJobs(ctx context.Context) ([]JobRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx,
		`SELECT id,sequence_id,params,status,progress,results,error,created_at,updated_at FROM analysis_jobs ORDER BY created_at, id`)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []JobRow
	for rows.Next() {
		var r JobRow
		if err := rows.Scan(&r.ID, &r.SequenceID, &r.Params, &r.Status, &r.Progress, &r.Results, &r.Error, &r.CreatedAt, &r.UpdatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

// ListPendingJobs returns jobs not in a terminal state (status pending/running).
func (s *Store) ListPendingJobs(ctx context.Context) ([]JobRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx,
		`SELECT id,sequence_id,params,status,progress,results,error,created_at,updated_at FROM analysis_jobs WHERE status IN ('pending','running') ORDER BY created_at, id`)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []JobRow
	for rows.Next() {
		var r JobRow
		if err := rows.Scan(&r.ID, &r.SequenceID, &r.Params, &r.Status, &r.Progress, &r.Results, &r.Error, &r.CreatedAt, &r.UpdatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

// StepRow is a single analysis step's progress record.
type StepRow struct {
	JobID      string `json:"job_id"`
	StepName   string `json:"step_name"`
	Status     string `json:"status"`
	ResultHash string `json:"result_hash"`
	UpdatedAt  string `json:"updated_at"`
}

// UpsertStep records a step's completion within a transaction (caller passes the
// tx so the step update and the job results update commit atomically).
func UpsertStep(tx DBTX, ctx context.Context, r StepRow) error {
	_, err := tx.ExecContext(ctx, `INSERT INTO job_steps(job_id,step_name,status,result_hash,updated_at)
		VALUES(?,?,?,?,?)
		ON CONFLICT(job_id,step_name) DO UPDATE SET status=excluded.status,result_hash=excluded.result_hash,updated_at=excluded.updated_at`,
		r.JobID, r.StepName, r.Status, r.ResultHash, r.UpdatedAt)
	return err
}

// ListSteps returns the step records for a job.
func (s *Store) ListSteps(ctx context.Context, jobID string) ([]StepRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx,
		`SELECT job_id,step_name,status,result_hash,updated_at FROM job_steps WHERE job_id=? ORDER BY step_name`, jobID)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []StepRow
	for rows.Next() {
		var r StepRow
		if err := rows.Scan(&r.JobID, &r.StepName, &r.Status, &r.ResultHash, &r.UpdatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

// MotifRow is a persisted motif.
type MotifRow struct {
	ID          string `json:"id"`
	Name        string `json:"name"`
	Pattern     string `json:"pattern"`
	Description string `json:"description"`
	CreatedAt   string `json:"created_at"`
}

func (s *Store) SaveMotif(ctx context.Context, r MotifRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO motifs(id,name,pattern,description,created_at) VALUES(?,?,?,?,?)
			ON CONFLICT(id) DO UPDATE SET name=excluded.name,pattern=excluded.pattern,description=excluded.description`,
			r.ID, r.Name, r.Pattern, r.Description, r.CreatedAt)
		return err
	})
}

func (s *Store) ListMotifs(ctx context.Context) ([]MotifRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx, `SELECT id,name,pattern,description,created_at FROM motifs ORDER BY id`)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []MotifRow
	for rows.Next() {
		var r MotifRow
		if err := rows.Scan(&r.ID, &r.Name, &r.Pattern, &r.Description, &r.CreatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

func (s *Store) GetMotif(ctx context.Context, id string) (MotifRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx, `SELECT id,name,pattern,description,created_at FROM motifs WHERE id=?`, id)
	var r MotifRow
	if err := row.Scan(&r.ID, &r.Name, &r.Pattern, &r.Description, &r.CreatedAt); err != nil {
		if isNoRows(err) {
			return MotifRow{}, ErrNotFound
		}
		return MotifRow{}, err
	}
	return r, nil
}

func (s *Store) DeleteMotif(ctx context.Context, id string) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `DELETE FROM motifs WHERE id=?`, id)
		return err
	})
}

// EnzymeRow is a persisted restriction enzyme.
type EnzymeRow struct {
	ID        string `json:"id"`
	Name      string `json:"name"`
	Site      string `json:"site"`
	CutOffset int    `json:"cut_offset"`
	CreatedAt string `json:"created_at"`
}

func (s *Store) SaveEnzyme(ctx context.Context, r EnzymeRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO enzymes(id,name,site,cut_offset,created_at) VALUES(?,?,?,?,?)
			ON CONFLICT(id) DO UPDATE SET name=excluded.name,site=excluded.site,cut_offset=excluded.cut_offset`,
			r.ID, r.Name, r.Site, r.CutOffset, r.CreatedAt)
		return err
	})
}

func (s *Store) ListEnzymes(ctx context.Context) ([]EnzymeRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx, `SELECT id,name,site,cut_offset,created_at FROM enzymes ORDER BY id`)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []EnzymeRow
	for rows.Next() {
		var r EnzymeRow
		if err := rows.Scan(&r.ID, &r.Name, &r.Site, &r.CutOffset, &r.CreatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

func (s *Store) GetEnzyme(ctx context.Context, id string) (EnzymeRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx, `SELECT id,name,site,cut_offset,created_at FROM enzymes WHERE id=?`, id)
	var r EnzymeRow
	if err := row.Scan(&r.ID, &r.Name, &r.Site, &r.CutOffset, &r.CreatedAt); err != nil {
		if isNoRows(err) {
			return EnzymeRow{}, ErrNotFound
		}
		return EnzymeRow{}, err
	}
	return r, nil
}

func (s *Store) DeleteEnzyme(ctx context.Context, id string) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `DELETE FROM enzymes WHERE id=?`, id)
		return err
	})
}

// FeatureRow is a persisted annotation feature.
type FeatureRow struct {
	ID         string `json:"id"`
	SequenceID string `json:"sequence_id"`
	Type       string `json:"type"`
	Start      int    `json:"start"`
	End        int    `json:"end"`
	Strand     string `json:"strand"`
	Label      string `json:"label"`
	CreatedAt  string `json:"created_at"`
}

func (s *Store) SaveFeature(ctx context.Context, r FeatureRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO features(id,sequence_id,type,start,end,strand,label,created_at) VALUES(?,?,?,?,?,?,?,?)`,
			r.ID, r.SequenceID, r.Type, r.Start, r.End, r.Strand, r.Label, r.CreatedAt)
		return err
	})
}

func (s *Store) ListFeatures(ctx context.Context, seqID string) ([]FeatureRow, error) {
	rows, err := s.WithCtx(ctx).QueryContext(ctx, `SELECT id,sequence_id,type,start,end,strand,label,created_at FROM features WHERE sequence_id=? ORDER BY start`, seqID)
	if err != nil {
		return nil, err
	}
	defer rows.Close()
	var out []FeatureRow
	for rows.Next() {
		var r FeatureRow
		if err := rows.Scan(&r.ID, &r.SequenceID, &r.Type, &r.Start, &r.End, &r.Strand, &r.Label, &r.CreatedAt); err != nil {
			return nil, err
		}
		out = append(out, r)
	}
	return out, rows.Err()
}

func (s *Store) GetFeature(ctx context.Context, id string) (FeatureRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx, `SELECT id,sequence_id,type,start,end,strand,label,created_at FROM features WHERE id=?`, id)
	var r FeatureRow
	if err := row.Scan(&r.ID, &r.SequenceID, &r.Type, &r.Start, &r.End, &r.Strand, &r.Label, &r.CreatedAt); err != nil {
		if isNoRows(err) {
			return FeatureRow{}, ErrNotFound
		}
		return FeatureRow{}, err
	}
	return r, nil
}

func (s *Store) DeleteFeature(ctx context.Context, id string) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `DELETE FROM features WHERE id=?`, id)
		return err
	})
}

// AlignmentRow is a persisted pairwise alignment result.
type AlignmentRow struct {
	ID        string `json:"id"`
	SeqAID    string `json:"seq_a_id"`
	SeqBID    string `json:"seq_b_id"`
	Mode      string `json:"mode"`
	Score     int    `json:"score"`
	AlignedA  string `json:"aligned_a"`
	AlignedB  string `json:"aligned_b"`
	Params    string `json:"params"`
	CreatedAt string `json:"created_at"`
}

func (s *Store) SaveAlignment(ctx context.Context, r AlignmentRow) error {
	return s.InTx(ctx, func(tx DBTX) error {
		_, err := tx.ExecContext(ctx, `INSERT INTO alignments(id,seq_a_id,seq_b_id,mode,score,aligned_a,aligned_b,params,created_at) VALUES(?,?,?,?,?,?,?,?,?)`,
			r.ID, r.SeqAID, r.SeqBID, r.Mode, r.Score, r.AlignedA, r.AlignedB, r.Params, r.CreatedAt)
		return err
	})
}

func (s *Store) GetAlignment(ctx context.Context, id string) (AlignmentRow, error) {
	row := s.WithCtx(ctx).QueryRowContext(ctx, `SELECT id,seq_a_id,seq_b_id,mode,score,aligned_a,aligned_b,params,created_at FROM alignments WHERE id=?`, id)
	var r AlignmentRow
	if err := row.Scan(&r.ID, &r.SeqAID, &r.SeqBID, &r.Mode, &r.Score, &r.AlignedA, &r.AlignedB, &r.Params, &r.CreatedAt); err != nil {
		if isNoRows(err) {
			return AlignmentRow{}, ErrNotFound
		}
		return AlignmentRow{}, err
	}
	return r, nil
}

// nowStamp returns an RFC3339 timestamp; kept here so queries don't depend on
// a global clock (the service layer injects timestamps for idempotency).
func nowStamp() string { return time.Now().UTC().Format(time.RFC3339Nano) }
