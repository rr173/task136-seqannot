package httpapi

import (
	"net/http"
)

func (h handlers) reconcile(w http.ResponseWriter, r *http.Request) {
	rows, err := h.svc.ReconcileAll(r.Context())
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) healthz(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(http.StatusOK)
	_, _ = w.Write([]byte(`{"status":"ok"}`))
}
