package httpapi

import (
	"net/http"
	"strconv"

	"task136-seqannot/internal/service"
)

func (h handlers) translate(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	res, err := h.svc.Translate(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, res)
}

func (h handlers) findORFs(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	minAA, _ := strconv.Atoi(r.URL.Query().Get("min_aa_len"))
	res, err := h.svc.FindORFs(r.Context(), id, minAA)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, res)
}

func (h handlers) motifSearch(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	var body struct {
		Pattern  string `json:"pattern"`
		MotifID  string `json:"motif_id"`
	}
	// motif-search accepts either a JSON body or query params for convenience
	if r.ContentLength > 0 {
		if err := readJSON(r, &body); err != nil {
			writeErr(w, err)
			return
		}
	} else {
		body.Pattern = r.URL.Query().Get("pattern")
		body.MotifID = r.URL.Query().Get("motif_id")
	}
	res, err := h.svc.SearchMotif(r.Context(), id, body.Pattern, body.MotifID)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, res)
}

func (h handlers) restrictionMap(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	ids, err := readIDs(r)
	if err != nil {
		writeErr(w, err)
		return
	}
	res, err := h.svc.RestrictionMap(r.Context(), id, ids)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, res)
}

func (h handlers) digest(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	ids, err := readIDs(r)
	if err != nil {
		writeErr(w, err)
		return
	}
	res, err := h.svc.Digest(r.Context(), id, ids)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, res)
}

func (h handlers) align(w http.ResponseWriter, r *http.Request) {
	var req service.AlignmentRequest
	if err := readJSON(r, &req); err != nil {
		writeErr(w, err)
		return
	}
	if req.Match == 0 && req.Mismatch == 0 && req.Gap == 0 {
		// apply sensible defaults: +1/-1/-1
		req.Match, req.Mismatch, req.Gap = 1, -1, -1
	}
	row, err := h.svc.Align(r.Context(), req)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) getAlignment(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	row, err := h.svc.GetAlignment(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, row)
}

// readIDs parses a JSON body {enzyme_ids:[...]} or a repeated query param
// enzyme_ids=a&enzyme_ids=b into a string slice.
func readIDs(r *http.Request) ([]string, error) {
	if r.ContentLength > 0 {
		var body struct {
			EnzymeIDs []string `json:"enzyme_ids"`
		}
		if err := readJSON(r, &body); err != nil {
			return nil, err
		}
		return body.EnzymeIDs, nil
	}
	q := r.URL.Query()["enzyme_ids"]
	return q, nil
}
