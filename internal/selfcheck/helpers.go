package selfcheck

import (
	"net/http"
	"os"
)

// mkTempDir creates a fresh temp directory for a scenario and returns its
// path. We use os.MkdirTemp so it works both inside go test and from the
// smoke-test entrypoint (which is not a test context). The directory is left
// for the OS to reap; scenarios use small DBs.
func mkTempDir() (string, error) {
	return os.MkdirTemp("", "seqannot-selfcheck-*")
}

// contains reports whether s contains substr.
func contains(s, substr string) bool {
	return len(s) >= len(substr) && indexOf(s, substr) >= 0
}

func indexOf(s, sub string) int {
	for i := 0; i+len(sub) <= len(s); i++ {
		if s[i:i+len(sub)] == sub {
			return i
		}
	}
	return -1
}

// asFloat coerces a JSON-decoded number (always float64) to a float64,
// returning 0 when the value is not a number. The selfcheck asserts against
// numeric fields that arrive as float64 (see any-numeric-int-vs-float64 which
// warns about float64-only comparisons for `any` config; here the bodies are
// JSON so float64 is correct).
func asFloat(v any) float64 {
	if f, ok := v.(float64); ok {
		return f
	}
	return 0
}

// asString coerces a JSON-decoded value to a string.
func asString(v any) string {
	if s, ok := v.(string); ok {
		return s
	}
	return ""
}

// statusOK reports whether a status is in the 2xx range.
func statusOK(s int) bool { return s >= 200 && s < 300 }

// noopSilence avoids "imported and not used" when a helper is conditionally
// referenced; http is used by the frontend scenario's Get.
var _ = http.Get
