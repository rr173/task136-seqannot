package selfcheck

import (
	"fmt"
	"net/http"
)

// scenarios is the ordered list of end-to-end checks Run() executes.
var scenarios = []scenario{
	{"sequence-crud", scenarioSequenceCRUD},
	{"composition-gc", scenarioComposition},
	{"translate-6frames", scenarioTranslate},
	{"orf-detection", scenarioORF},
	{"motif-search", scenarioMotif},
	{"restriction-digest", scenarioRestriction},
	{"alignment", scenarioAlignment},
	{"job-submit-replay", scenarioJobReplay},
	{"restart-recovery", scenarioRestart},
	{"features", scenarioFeatures},
	{"frontend-served", scenarioFrontend},
	{"error-handling", scenarioErrors},
}

// makeSeq creates a sequence and returns its id (or an error).
func makeSeq(b *srvBundle, residues, typ string) (string, error) {
	body := map[string]any{"name": "s1", "residues": residues, "type": typ}
	var resp map[string]any
	if err := decodeOK(b.server, "POST", "/sequences", body, &resp); err != nil {
		return "", err
	}
	id := asString(resp["id"])
	if id == "" {
		return "", fmt.Errorf("no id in create response: %+v", resp)
	}
	return id, nil
}

func scenarioSequenceCRUD() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGCAAATG", "linear")
	if err != nil {
		return err
	}
	var got map[string]any
	if err := decodeOK(b.server, "GET", "/sequences/"+id, nil, &got); err != nil {
		return err
	}
	if asString(got["residues"]) != "ATGCAAATG" {
		return fmt.Errorf("get residues=%v want ATGCAAATG", got["residues"])
	}
	var list []map[string]any
	if err := decodeOK(b.server, "GET", "/sequences", nil, &list); err != nil {
		return err
	}
	if len(list) != 1 {
		return fmt.Errorf("list len=%d want 1", len(list))
	}
	if _, err := requireOK(b.server, "DELETE", "/sequences/"+id, nil); err != nil {
		return err
	}
	status, _, err := do(b.server, "GET", "/sequences/"+id, nil)
	if err != nil {
		return err
	}
	if status != http.StatusNotFound {
		return fmt.Errorf("after delete get status=%d want 404", status)
	}
	return nil
}

func scenarioComposition() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGCATGC", "linear")
	if err != nil {
		return err
	}
	var c map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+id+"/composition", nil, &c); err != nil {
		return err
	}
	if asFloat(c["length"]) != 8 {
		return fmt.Errorf("length=%v want 8", c["length"])
	}
	if asFloat(c["gc_bp"]) != 5000 {
		return fmt.Errorf("gc_bp=%v want 5000", c["gc_bp"])
	}
	return nil
}

func scenarioTranslate() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGAAATAA", "linear")
	if err != nil {
		return err
	}
	var res map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+id+"/translate", nil, &res); err != nil {
		return err
	}
	frames, ok := res["frames"].([]any)
	if !ok || len(frames) != 6 {
		return fmt.Errorf("want 6 frames got %v", res["frames"])
	}
	f0 := frames[0].(map[string]any)
	if asString(f0["aminos"]) != "MK*" {
		return fmt.Errorf("frame0 aminos=%v want MK*", f0["aminos"])
	}
	return nil
}

func scenarioORF() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGAAATAA", "linear")
	if err != nil {
		return err
	}
	var res map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+id+"/orfs?min_aa_len=1", nil, &res); err != nil {
		return err
	}
	orfs, ok := res["orfs"].([]any)
	if !ok || len(orfs) == 0 {
		return fmt.Errorf("want >=1 orf got %v", res["orfs"])
	}
	// verify the frame0 ORF protein is MK*
	var foundMK bool
	for _, o := range orfs {
		om := o.(map[string]any)
		if asString(om["protein"]) == "MK*" {
			foundMK = true
		}
	}
	if !foundMK {
		return fmt.Errorf("frame0 MK* ORF not found among %+v", orfs)
	}
	return nil
}

func scenarioMotif() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "AGCAGCAGC", "linear")
	if err != nil {
		return err
	}
	var res map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+id+"/motif-search",
		map[string]any{"pattern": "RGC"}, &res); err != nil {
		return err
	}
	hits, ok := res["hits"].([]any)
	if !ok || len(hits) == 0 {
		return fmt.Errorf("want >=1 motif hit got %v", res["hits"])
	}
	return nil
}

func scenarioRestriction() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	var enz map[string]any
	if err := decodeOK(b.server, "POST", "/enzymes",
		map[string]any{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1}, &enz); err != nil {
		return err
	}
	eid := asString(enz["id"])
	seqID, err := makeSeq(b, "GAATTCGAATTC", "linear")
	if err != nil {
		return err
	}
	var rmap map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+seqID+"/restriction-map",
		map[string]any{"enzyme_ids": []string{eid}}, &rmap); err != nil {
		return err
	}
	sites, ok := rmap["sites"].([]any)
	if !ok || len(sites) != 2 {
		return fmt.Errorf("want 2 sites got %v", rmap["sites"])
	}
	var dig map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+seqID+"/digest",
		map[string]any{"enzyme_ids": []string{eid}}, &dig); err != nil {
		return err
	}
	frags, ok := dig["fragments"].([]any)
	if !ok {
		return fmt.Errorf("no fragments: %v", dig)
	}
	if len(frags) != 3 {
		return fmt.Errorf("want 3 fragments got %d", len(frags))
	}
	return nil
}

func scenarioAlignment() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	aID, err := makeSeq(b, "ACGT", "linear")
	if err != nil {
		return err
	}
	bID, err := makeSeq(b, "ACGGT", "linear")
	if err != nil {
		return err
	}
	var aln map[string]any
	if err := decodeOK(b.server, "POST", "/alignments",
		map[string]any{"seq_a_id": aID, "seq_b_id": bID, "mode": "global", "match": 1, "mismatch": -1, "gap": -1}, &aln); err != nil {
		return err
	}
	if asFloat(aln["score"]) != 3 {
		return fmt.Errorf("score=%v want 3", aln["score"])
	}
	return nil
}

func scenarioJobReplay() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGAAATAAATGAAATAA", "linear")
	if err != nil {
		return err
	}
	var job map[string]any
	if err := decodeOK(b.server, "POST", "/jobs", map[string]any{
		"sequence_id": id,
		"steps": []map[string]any{
			{"name": "composition"},
			{"name": "translate"},
			{"name": "orf", "params": map[string]any{"min_aa_len": 1}},
		},
	}, &job); err != nil {
		return err
	}
	if asString(job["status"]) != "done" {
		return fmt.Errorf("job status=%v want done", job["status"])
	}
	jid := asString(job["id"])
	var rep map[string]any
	if err := decodeOK(b.server, "POST", "/jobs/"+jid+"/replay", nil, &rep); err != nil {
		return err
	}
	if rep["equal"] != true {
		return fmt.Errorf("replay equal=%v want true: %+v", rep["equal"], rep)
	}
	return nil
}

func scenarioRestart() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGAAATAA", "linear")
	if err != nil {
		return err
	}
	var job map[string]any
	if err := decodeOK(b.server, "POST", "/jobs", map[string]any{
		"sequence_id": id,
		"steps": []map[string]any{
			{"name": "composition"},
			{"name": "translate"},
		},
	}, &job); err != nil {
		return err
	}
	if asString(job["status"]) != "done" {
		return fmt.Errorf("status=%v want done", job["status"])
	}
	if _, err := requireOK(b.server, "POST", "/admin/reconcile", nil); err != nil {
		return err
	}
	var resumed map[string]any
	if err := decodeOK(b.server, "POST", "/jobs/"+asString(job["id"])+"/resume", nil, &resumed); err != nil {
		return err
	}
	if asString(resumed["status"]) != "done" {
		return fmt.Errorf("resume status=%v want done", resumed["status"])
	}
	return nil
}

func scenarioFeatures() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	id, err := makeSeq(b, "ATGCAAATGCAATG", "linear")
	if err != nil {
		return err
	}
	var feat map[string]any
	if err := decodeOK(b.server, "POST", "/sequences/"+id+"/features",
		map[string]any{"type": "gene", "start": 1, "end": 6, "strand": "+", "label": "g1"}, &feat); err != nil {
		return err
	}
	var feats []map[string]any
	if err := decodeOK(b.server, "GET", "/sequences/"+id+"/features?start=3&end=5", nil, &feats); err != nil {
		return err
	}
	if len(feats) != 1 {
		return fmt.Errorf("overlap query len=%d want 1", len(feats))
	}
	status, _, err := do(b.server, "POST", "/sequences/"+id+"/features",
		map[string]any{"type": "x", "start": 6, "end": 3, "strand": "+"})
	if err != nil {
		return err
	}
	if status != http.StatusBadRequest {
		return fmt.Errorf("bad coords status=%d want 400", status)
	}
	return nil
}

func scenarioFrontend() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	// the index page must be served
	resp, err := http.Get(b.server.URL + "/")
	if err != nil {
		return fmt.Errorf("GET /: %w", err)
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("GET / status=%d want 200", resp.StatusCode)
	}
	// the frontend's app.js must be served and must reference a real business
	// API endpoint, proving the UI is wired to the Go computation engine.
	resp2, err := http.Get(b.server.URL + "/app.js")
	if err != nil {
		return fmt.Errorf("GET /app.js: %w", err)
	}
	defer resp2.Body.Close()
	if resp2.StatusCode != http.StatusOK {
		return fmt.Errorf("GET /app.js status=%d want 200", resp2.StatusCode)
	}
	body, err := readAll(resp2.Body)
	if err != nil {
		return err
	}
	if !contains(string(body), "/sequences") {
		return fmt.Errorf("frontend app.js does not reference /sequences")
	}
	return nil
}

func scenarioErrors() error {
	b, err := newServer()
	if err != nil {
		return err
	}
	defer b.cleanup()
	status, _, err := do(b.server, "GET", "/sequences/missing", nil)
	if err != nil {
		return err
	}
	if status != http.StatusNotFound {
		return fmt.Errorf("missing seq status=%d want 404", status)
	}
	status, _, err = do(b.server, "POST", "/sequences", map[string]any{"name": "x", "residues": "XYZ", "type": "linear"})
	if err != nil {
		return err
	}
	if status != http.StatusBadRequest {
		return fmt.Errorf("invalid base status=%d want 400", status)
	}
	status, _, err = do(b.server, "POST", "/motifs", map[string]any{"name": "x", "pattern": "ZQ"})
	if err != nil {
		return err
	}
	if status != http.StatusBadRequest {
		return fmt.Errorf("invalid motif status=%d want 400", status)
	}
	return nil
}
