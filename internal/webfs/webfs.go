package webfs

import (
	"embed"
	"io/fs"
	"net/http"
)

// The web/ directory holds the native (no-build) frontend. It is embedded so
// the single Go binary serves both the API and the UI; main and selfcheck
// share the same embedded filesystem via this package.
//
//go:embed web/*
var webFS embed.FS

// WebFS returns the embedded filesystem rooted at web/ as an fs.FS, suitable
// for http.FileServerFS (Go 1.22+).
func WebFS() fs.FS {
	sub, err := fs.Sub(webFS, "web")
	if err != nil {
		// should never happen with a compile-time embed
		panic(err)
	}
	return sub
}

// WebFSHTTP returns the embedded filesystem as an http.FileSystem for callers
// that prefer the older http.FileServer (kept for compatibility).
func WebFSHTTP() http.FileSystem {
	return http.FS(WebFS())
}

// WebIndex returns the raw bytes of index.html, used by the smoke test to
// assert the frontend is served.
func WebIndex() ([]byte, error) {
	return webFS.ReadFile("web/index.html")
}

