package httpapi

import (
	"encoding/json"
	"errors"
	"net/http"

	"task136-seqannot/internal/service"
	"task136-seqannot/internal/store"
)

// writeJSON marshals v as JSON and writes it with the given status. It also
// sets the content type; on marshal failure it writes a 500.
func writeJSON(w http.ResponseWriter, status int, v any) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(status)
	if err := json.NewEncoder(w).Encode(v); err != nil {
		// best-effort; header already sent
		_, _ = w.Write([]byte(`{"error":"encode failed"}`))
	}
}

// readJSON decodes the request body into v. It limits the body to 8 MiB to
// avoid unbounded reads of large sequence payloads.
func readJSON(r *http.Request, v any) error {
	r.Body = http.MaxBytesReader(nil, r.Body, 8<<20)
	dec := json.NewDecoder(r.Body)
	dec.DisallowUnknownFields()
	if err := dec.Decode(v); err != nil {
		return err
	}
	return nil
}

// mapError translates a service/store error into an HTTP status + message. The
// conventions are: ErrNotFound -> 404; validation failures (sentinel prefix
// "seq:", "orf:", "motif:", "restrict:", "annot:", "service:") -> 400;
// everything else -> 500.
func mapError(err error) (int, string) {
	if err == nil {
		return http.StatusOK, ""
	}
	if errors.Is(err, store.ErrNotFound) {
		return http.StatusNotFound, err.Error()
	}
	msg := err.Error()
	switch {
	case hasPrefix(msg, "seq:"), hasPrefix(msg, "orf:"), hasPrefix(msg, "motif:"),
		hasPrefix(msg, "restrict:"), hasPrefix(msg, "annot:"), hasPrefix(msg, "service:"),
		hasPrefix(msg, "fasta:"):
		return http.StatusBadRequest, msg
	}
	return http.StatusInternalServerError, msg
}

func hasPrefix(s, p string) bool {
	return len(s) >= len(p) && s[:len(p)] == p
}

// writeErr writes a mapped error as a JSON {"error":...} body.
func writeErr(w http.ResponseWriter, err error) {
	status, msg := mapError(err)
	if status == http.StatusOK {
		status = http.StatusInternalServerError
	}
	writeJSON(w, status, map[string]string{"error": msg})
}

// ok writes a 200 with the JSON body.
func ok(w http.ResponseWriter, v any) { writeJSON(w, http.StatusOK, v) }
func created(w http.ResponseWriter, v any) { writeJSON(w, http.StatusCreated, v) }

// --- Sequence handlers ---

func (h handlers) createSequence(w http.ResponseWriter, r *http.Request) {
	var req service.SequenceRequest
	if err := readJSON(r, &req); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.CreateSequence(r.Context(), req)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) importFasta(w http.ResponseWriter, r *http.Request) {
	var body struct {
		Fasta string `json:"fasta"`
	}
	if err := readJSON(r, &body); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.ImportFasta(r.Context(), body.Fasta)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) listSequences(w http.ResponseWriter, r *http.Request) {
	rows, err := h.svc.ListSequences(r.Context())
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) getSequence(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	row, err := h.svc.GetSequence(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, row)
}

func (h handlers) deleteSequence(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	if err := h.svc.DeleteSequence(r.Context(), id); err != nil {
		writeErr(w, err)
		return
	}
	w.WriteHeader(http.StatusNoContent)
}

func (h handlers) composition(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	c, err := h.svc.Composition(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, c)
}

func (h handlers) reverseComplement(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	rc, err := h.svc.ReverseComplement(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, map[string]string{"reverse_complement": rc})
}
