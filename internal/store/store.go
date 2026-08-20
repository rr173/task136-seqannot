package store

import (
	"context"
	"database/sql"
	"errors"
	"fmt"

	_ "modernc.org/sqlite"
)

// DBTX is the common interface over *sql.DB and *sql.Tx so that callers can
// run queries either at the top level or inside a transaction with the same
// code. This is required because the engine runs with SetMaxOpenConns(1); any
// read performed inside a transaction must go through the tx handle, not a
// second connection from the pool, or it will deadlock (see intx-writes-must-use-tx).
type DBTX interface {
	ExecContext(ctx context.Context, query string, args ...any) (sql.Result, error)
	QueryContext(ctx context.Context, query string, args ...any) (*sql.Rows, error)
	QueryRowContext(ctx context.Context, query string, args ...any) *sql.Row
}

// Store wraps a SQLite database connection with the engine's schema and the
// CRUD queries. It is concurrency-safe only via the single-connection pool;
// callers must wrap multi-step operations in InTx.
type Store struct {
	db *sql.DB
}

// Open creates or opens a SQLite database at path and applies the schema
// migrations. The connection pool is pinned to a single open connection so
// that writers and transactional readers never block each other across
// separate connections.
func Open(path string) (*Store, error) {
	dsn := path + "?_pragma=busy_timeout(5000)&_pragma=journal_mode(WAL)&_pragma=foreign_keys(on)"
	db, err := sql.Open("sqlite", dsn)
	if err != nil {
		return nil, fmt.Errorf("store: open: %w", err)
	}
	db.SetMaxOpenConns(1)
	if _, err := db.Exec(migrations); err != nil {
		db.Close()
		return nil, fmt.Errorf("store: migrate: %w", err)
	}
	return &Store{db: db}, nil
}

// Close releases the database connection.
func (s *Store) Close() error { return s.db.Close() }

// DB exposes the underlying connection for the rare caller that needs direct
// access (e.g. the httpapi health check). Transactional code must use InTx.
func (s *Store) DB() *sql.DB { return s.db }

// InTx runs fn inside a single transaction. fn receives a DBTX bound to the
// transaction; every read and write inside fn MUST go through that DBTX (not
// s.db), otherwise the single-connection pool deadlocks waiting for a second
// connection. The transaction is committed on a nil return and rolled back on
// any error.
func (s *Store) InTx(ctx context.Context, fn func(DBTX) error) error {
	tx, err := s.db.BeginTx(ctx, nil)
	if err != nil {
		return fmt.Errorf("store: begin: %w", err)
	}
	defer func() {
		if p := recover(); p != nil {
			_ = tx.Rollback()
			panic(p)
		}
	}()
	if err := fn(tx); err != nil {
		_ = tx.Rollback()
		return err
	}
	if err := tx.Commit(); err != nil {
		return fmt.Errorf("store: commit: %w", err)
	}
	return nil
}

// WithCtx returns a DBTX bound to the top-level connection (NOT a transaction).
// Use it for single-statement reads. Multi-statement sequences must use InTx.
func (s *Store) WithCtx(ctx context.Context) DBTX {
	return ctxStmt{db: s.db, ctx: ctx}
}

// ctxStmt adapts *sql.DB + a context to the DBTX interface for single-statement
// reads outside a transaction.
type ctxStmt struct {
	db  *sql.DB
	ctx context.Context
}

func (c ctxStmt) ExecContext(_ context.Context, q string, args ...any) (sql.Result, error) {
	return c.db.ExecContext(c.ctx, q, args...)
}
func (c ctxStmt) QueryContext(_ context.Context, q string, args ...any) (*sql.Rows, error) {
	return c.db.QueryContext(c.ctx, q, args...)
}
func (c ctxStmt) QueryRowContext(_ context.Context, q string, args ...any) *sql.Row {
	return c.db.QueryRowContext(c.ctx, q, args...)
}

// ErrNotFound is returned by single-row getters when no row matches.
var ErrNotFound = errors.New("store: not found")

// migrations is the schema applied on every Open. It is idempotent (CREATE
// TABLE IF NOT EXISTS). The schema is the authoritative data model; see the
// design doc for the entity descriptions.
const migrations = `
CREATE TABLE IF NOT EXISTS sequences (
	id            TEXT PRIMARY KEY,
	name          TEXT NOT NULL,
	type          TEXT NOT NULL,
	residues       TEXT NOT NULL,
	description   TEXT NOT NULL DEFAULT '',
	fasta_header  TEXT NOT NULL DEFAULT '',
	length        INTEGER NOT NULL,
	created_at    TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS analysis_jobs (
	id          TEXT PRIMARY KEY,
	sequence_id TEXT NOT NULL,
	params      TEXT NOT NULL,
	status      TEXT NOT NULL,
	progress    INTEGER NOT NULL DEFAULT 0,
	results     TEXT NOT NULL DEFAULT '{}',
	error       TEXT NOT NULL DEFAULT '',
	created_at  TEXT NOT NULL,
	updated_at  TEXT NOT NULL,
	FOREIGN KEY (sequence_id) REFERENCES sequences(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS job_steps (
	job_id      TEXT NOT NULL,
	step_name   TEXT NOT NULL,
	status      TEXT NOT NULL DEFAULT 'pending',
	result_hash TEXT NOT NULL DEFAULT '',
	updated_at  TEXT NOT NULL,
	PRIMARY KEY (job_id, step_name),
	FOREIGN KEY (job_id) REFERENCES analysis_jobs(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS motifs (
	id          TEXT PRIMARY KEY,
	name        TEXT NOT NULL,
	pattern     TEXT NOT NULL,
	description TEXT NOT NULL DEFAULT '',
	created_at  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS enzymes (
	id        TEXT PRIMARY KEY,
	name      TEXT NOT NULL,
	site      TEXT NOT NULL,
	cut_offset INTEGER NOT NULL,
	created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS features (
	id          TEXT PRIMARY KEY,
	sequence_id TEXT NOT NULL,
	type        TEXT NOT NULL,
	start       INTEGER NOT NULL,
	end         INTEGER NOT NULL,
	strand      TEXT NOT NULL DEFAULT '+',
	label       TEXT NOT NULL DEFAULT '',
	created_at  TEXT NOT NULL,
	FOREIGN KEY (sequence_id) REFERENCES sequences(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_features_seq ON features(sequence_id);

CREATE TABLE IF NOT EXISTS alignments (
	id       TEXT PRIMARY KEY,
	seq_a_id TEXT NOT NULL,
	seq_b_id TEXT NOT NULL,
	mode     TEXT NOT NULL,
	score    INTEGER NOT NULL,
	aligned_a TEXT NOT NULL,
	aligned_b TEXT NOT NULL,
	params   TEXT NOT NULL,
	created_at TEXT NOT NULL
);
`
