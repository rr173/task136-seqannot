package selfcheck

import "io"

// readAll is a thin wrapper over io.ReadAll, kept here so scenarios.go does
// not need to import io directly.
func readAll(r io.Reader) ([]byte, error) { return io.ReadAll(r) }
