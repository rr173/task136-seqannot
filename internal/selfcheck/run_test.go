package selfcheck

import "testing"

// TestRun exposes Run() to `go test`. The smoke-test entrypoint in main calls
// Run() directly; this test gives CI the same coverage.
func TestRun(t *testing.T) {
	if err := Run(); err != nil {
		t.Fatal(err)
	}
}
