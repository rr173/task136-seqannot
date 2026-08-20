package httpapi

import (
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestBug07_JSONTrailingValuesRejected(t *testing.T) {
	req := httptest.NewRequest(http.MethodPost, "/", strings.NewReader(`{"x":1}{"x":2}`))
	if err := readJSON(req, &struct{ X int `json:"x"` }{}); err == nil {
		t.Fatal("trailing JSON value was accepted")
	}
}
