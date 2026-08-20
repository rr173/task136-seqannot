package httpapi

import (
	"net/http"
	"strconv"

	"task136-seqannot/internal/service"
	"task136-seqannot/internal/store"
)

func (h handlers) submitJob(w http.ResponseWriter, r *http.Request) {
	var spec service.JobSpec
	if err := readJSON(r, &spec); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.SubmitJob(r.Context(), spec)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) listJobs(w http.ResponseWriter, r *http.Request) {
	rows, err := h.svc.ListJobs(r.Context())
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) getJob(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	row, err := h.svc.GetJob(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, row)
}

func (h handlers) resumeJob(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	row, err := h.svc.ResumeJob(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, row)
}

func (h handlers) replayJob(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	rep, err := h.svc.ReplayJob(r.Context(), id)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rep)
}

// Motif/Enzyme/Feature handlers ----------------------------------------------

func (h handlers) createMotif(w http.ResponseWriter, r *http.Request) {
	var body struct {
		Name    string `json:"name"`
		Pattern string `json:"pattern"`
		Description string `json:"description"`
	}
	if err := readJSON(r, &body); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.CreateMotif(r.Context(), body.Name, body.Pattern, body.Description)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) listMotifs(w http.ResponseWriter, r *http.Request) {
	rows, err := h.svc.ListMotifs(r.Context())
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) deleteMotif(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	if err := h.svc.DeleteMotif(r.Context(), id); err != nil {
		writeErr(w, err)
		return
	}
	w.WriteHeader(http.StatusNoContent)
}

func (h handlers) createEnzyme(w http.ResponseWriter, r *http.Request) {
	var body struct {
		Name string `json:"name"`
		Site string `json:"site"`
		Cut  int    `json:"cut_offset"`
	}
	if err := readJSON(r, &body); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.CreateEnzyme(r.Context(), body.Name, body.Site, body.Cut)
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) listEnzymes(w http.ResponseWriter, r *http.Request) {
	rows, err := h.svc.ListEnzymes(r.Context())
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) deleteEnzyme(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	if err := h.svc.DeleteEnzyme(r.Context(), id); err != nil {
		writeErr(w, err)
		return
	}
	w.WriteHeader(http.StatusNoContent)
}

func (h handlers) createFeature(w http.ResponseWriter, r *http.Request) {
	seqID := r.PathValue("id")
	var body struct {
		Type   string `json:"type"`
		Start  int    `json:"start"`
		End    int    `json:"end"`
		Strand string `json:"strand"`
		Label  string `json:"label"`
	}
	if err := readJSON(r, &body); err != nil {
		writeErr(w, err)
		return
	}
	row, err := h.svc.CreateFeature(r.Context(), store.FeatureRow{
		SequenceID: seqID, Type: body.Type, Start: body.Start, End: body.End,
		Strand: body.Strand, Label: body.Label,
	})
	if err != nil {
		writeErr(w, err)
		return
	}
	created(w, row)
}

func (h handlers) listFeatures(w http.ResponseWriter, r *http.Request) {
	seqID := r.PathValue("id")
	start, _ := strconv.Atoi(r.URL.Query().Get("start"))
	end, _ := strconv.Atoi(r.URL.Query().Get("end"))
	rows, err := h.svc.QueryFeatures(r.Context(), seqID, start, end)
	if err != nil {
		writeErr(w, err)
		return
	}
	ok(w, rows)
}

func (h handlers) deleteFeature(w http.ResponseWriter, r *http.Request) {
	id := r.PathValue("id")
	if err := h.svc.DeleteFeature(r.Context(), id); err != nil {
		writeErr(w, err)
		return
	}
	w.WriteHeader(http.StatusNoContent)
}
