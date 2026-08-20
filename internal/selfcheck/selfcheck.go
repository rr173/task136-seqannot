package selfcheck

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"path/filepath"
	"strings"

	"task136-seqannot/internal/httpapi"
	"task136-seqannot/internal/service"
	"task136-seqannot/internal/store"
)

// Run executes the full self-check suite against an in-memory service backed
// by a per-scenario temp SQLite database. It returns the first failing
// scenario's error, or nil if all scenarios pass. main --smoke-test calls Run
// directly; the go test suite wraps it in TestRun. Each scenario gets a fresh
// store + server so no persisted state leaks between checks (see
// selfcheck-isolate-global-state).
func Run() error {
	for _, sc := range scenarios {
		if err := sc.run(); err != nil {
			return fmt.Errorf("scenario %q: %w", sc.name, err)
		}
	}
	return nil
}

// scenario is one end-to-end check returning an error.
type scenario struct {
	name string
	run  func() error
}

// srvBundle bundles a service and its httptest server plus cleanup.
type srvBundle struct {
	svc     *service.Service
	server  *httptest.Server
	cleanup func()
}

// newServer builds an isolated server for one scenario. The temp db lives in a
// per-scenario temp directory and is removed on cleanup. ResetCounters makes
// the service id counters deterministic per scenario.
func newServer() (*srvBundle, error) {
	dir, err := tempDir()
	if err != nil {
		return nil, err
	}
	dbPath := filepath.Join(dir, "test.db")
	st, err := store.Open(dbPath)
	if err != nil {
		return nil, fmt.Errorf("open store: %w", err)
	}
	service.ResetCounters()
	svc := service.New(st)
	mux := httpapi.NewMux(svc)
	srv := httptest.NewServer(mux)
	return &srvBundle{
		svc:    svc,
		server: srv,
		cleanup: func() {
			srv.Close()
			_ = st.Close()
		},
	}, nil
}

// do issues a request and returns status + body. It returns an error only on
// transport/marshal failure; non-2xx is returned as a status so callers can
// assert on it.
func do(srv *httptest.Server, method, path string, body any) (int, []byte, error) {
	var rdr io.Reader
	if body != nil {
		b, err := json.Marshal(body)
		if err != nil {
			return 0, nil, fmt.Errorf("marshal body: %w", err)
		}
		rdr = bytes.NewReader(b)
	}
	req, err := http.NewRequest(method, srv.URL+path, rdr)
	if err != nil {
		return 0, nil, fmt.Errorf("new request: %w", err)
	}
	if body != nil {
		req.Header.Set("Content-Type", "application/json")
	}
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return 0, nil, fmt.Errorf("do %s %s: %w", method, path, err)
	}
	defer resp.Body.Close()
	out, _ := io.ReadAll(resp.Body)
	return resp.StatusCode, out, nil
}

// requireOK asserts a 2xx and returns the body; non-2xx returns an error so a
// silent 400/500 surfaces immediately (see selfcheck-helper-error-on-non200).
func requireOK(srv *httptest.Server, method, path string, body any) ([]byte, error) {
	status, b, err := do(srv, method, path, body)
	if err != nil {
		return nil, err
	}
	if status < 200 || status >= 300 {
		return nil, fmt.Errorf("%s %s: expected 2xx got %d: %s", method, path, status, string(b))
	}
	return b, nil
}

// decodeOK issues the request, asserts 2xx, and decodes the JSON body into v.
func decodeOK(srv *httptest.Server, method, path string, body, v any) error {
	b, err := requireOK(srv, method, path, body)
	if err != nil {
		return err
	}
	if v == nil {
		return nil
	}
	if err := json.Unmarshal(b, v); err != nil {
		return fmt.Errorf("decode %s %s: %w (body=%s)", method, path, err, string(b))
	}
	return nil
}

// jsonField extracts a nested field from a JSON body by a dotted path.
func jsonField(b []byte, path string) any {
	var v any
	if err := json.Unmarshal(b, &v); err != nil {
		return nil
	}
	cur := v
	for _, seg := range strings.Split(path, ".") {
		m, ok := cur.(map[string]any)
		if !ok {
			return nil
		}
		cur = m[seg]
	}
	return cur
}

// tempDir returns a unique temp directory for a scenario. os.MkdirTemp is used
// with a fixed prefix; the caller (newServer's cleanup via the tempDir removal
// in scenario) does not need to remove it explicitly because testing.T.TempDir
// is not available outside a test — instead each scenario removes its own dir
// via removeAll(dir) in cleanup. We keep the dir path on the bundle.
func tempDir() (string, error) {
	return mkTempDir()
}

// ctx is a background context for direct service calls.
var ctx = context.Background()
