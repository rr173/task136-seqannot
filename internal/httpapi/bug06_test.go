package httpapi

import (
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"path/filepath"
	"strings"
	"testing"

	"task136-seqannot/internal/service"
	"task136-seqannot/internal/store"
)

func TestBug06_ChunkedMotifBodyIsDecoded(t *testing.T) {
	dir := t.TempDir()
	st, err := store.Open(filepath.Join(dir, "test.db"))
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()
	service.ResetCounters()
	mux := NewMux(service.New(st))
	create := httptest.NewRequest(http.MethodPost, "/sequences", strings.NewReader(`{"name":"s","residues":"AACG","type":"linear"}`))
	create.Header.Set("Content-Type", "application/json")
	create.ContentLength = -1
	rr := httptest.NewRecorder()
	mux.ServeHTTP(rr, create)
	if rr.Code != http.StatusCreated {
		t.Fatalf("create status=%d body=%s", rr.Code, rr.Body.String())
	}
	seqID := ""
	var created struct{ ID string `json:"id"` }
	if err := json.Unmarshal(rr.Body.Bytes(), &created); err != nil {
		t.Fatal(err)
	}
	seqID = created.ID
	req := httptest.NewRequest(http.MethodPost, "/sequences/"+seqID+"/motif-search", io.NopCloser(strings.NewReader(`{"pattern":"AAC"}`)))
	req.Header.Set("Content-Type", "application/json")
	req.ContentLength = -1
	rr = httptest.NewRecorder()
	mux.ServeHTTP(rr, req)
	if rr.Code != http.StatusOK || !strings.Contains(rr.Body.String(), `"pattern":"AAC"`) {
		t.Fatalf("chunked motif body was not decoded: status=%d body=%s", rr.Code, rr.Body.String())
	}
}
