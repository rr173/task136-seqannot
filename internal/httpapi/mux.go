package httpapi

import (
	"net/http"

	"task136-seqannot/internal/service"
	"task136-seqannot/internal/webfs"
)

// NewMux wires every route onto a fresh mux and returns it. Both main and
// selfcheck call this so they share the exact same handlers and route table.
// The service is injected (no package-global state) so each httptest scenario
// gets isolated storage.
func NewMux(svc *service.Service) http.Handler {
	mux := http.NewServeMux()
	h := handlers{svc: svc}

	// Sequences
	mux.HandleFunc("POST /sequences", h.createSequence)
	mux.HandleFunc("GET /sequences", h.listSequences)
	mux.HandleFunc("GET /sequences/{id}", h.getSequence)
	mux.HandleFunc("DELETE /sequences/{id}", h.deleteSequence)
	mux.HandleFunc("POST /sequences/{id}/composition", h.composition)
	mux.HandleFunc("POST /sequences/{id}/reverse-complement", h.reverseComplement)
	mux.HandleFunc("POST /sequences/{id}/translate", h.translate)
	mux.HandleFunc("POST /sequences/{id}/orfs", h.findORFs)
	mux.HandleFunc("POST /sequences/{id}/motif-search", h.motifSearch)
	mux.HandleFunc("POST /sequences/{id}/restriction-map", h.restrictionMap)
	mux.HandleFunc("POST /sequences/{id}/digest", h.digest)
	mux.HandleFunc("POST /sequences/import-fasta", h.importFasta)

	// Motifs
	mux.HandleFunc("POST /motifs", h.createMotif)
	mux.HandleFunc("GET /motifs", h.listMotifs)
	mux.HandleFunc("DELETE /motifs/{id}", h.deleteMotif)

	// Enzymes
	mux.HandleFunc("POST /enzymes", h.createEnzyme)
	mux.HandleFunc("GET /enzymes", h.listEnzymes)
	mux.HandleFunc("DELETE /enzymes/{id}", h.deleteEnzyme)

	// Alignments
	mux.HandleFunc("POST /alignments", h.align)
	mux.HandleFunc("GET /alignments/{id}", h.getAlignment)

	// Jobs
	mux.HandleFunc("POST /jobs", h.submitJob)
	mux.HandleFunc("GET /jobs", h.listJobs)
	mux.HandleFunc("GET /jobs/{id}", h.getJob)
	mux.HandleFunc("POST /jobs/{id}/resume", h.resumeJob)
	mux.HandleFunc("POST /jobs/{id}/replay", h.replayJob)

	// Features
	mux.HandleFunc("POST /sequences/{id}/features", h.createFeature)
	mux.HandleFunc("GET /sequences/{id}/features", h.listFeatures)
	mux.HandleFunc("DELETE /features/{id}", h.deleteFeature)

	// Admin
	mux.HandleFunc("POST /admin/reconcile", h.reconcile)

	// Health + frontend
	mux.HandleFunc("GET /healthz", h.healthz)
	mux.Handle("GET /", http.FileServerFS(webfs.WebFS()))

	return mux
}

// handlers holds the injected service for the route closures.
type handlers struct {
	svc *service.Service
}
