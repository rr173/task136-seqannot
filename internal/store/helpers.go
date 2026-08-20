package store

import (
	"database/sql"
	"errors"
)

// isNoRows reports whether err is sql.ErrNoRows. Wrapped to centralize the
// dependency on the database/sql sentinel for callers that compare with
// ErrNotFound instead.
func isNoRows(err error) bool {
	return errors.Is(err, sql.ErrNoRows)
}
