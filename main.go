package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"net/http"
	"os"
	"os/signal"
	"syscall"
	"time"

	"task136-seqannot/internal/httpapi"
	"task136-seqannot/internal/selfcheck"
	"task136-seqannot/internal/service"
	"task136-seqannot/internal/store"
)

func main() {
	// Recognize the bare --smoke-test / -smoke-test flag before the flag
	// package parses, so `go run . --smoke-test` works as the skill mandates.
	var smoke bool
	args := os.Args[1:]
	filtered := make([]string, 0, len(args))
	for _, a := range args {
		switch a {
		case "--smoke-test", "-smoke-test":
			smoke = true
		default:
			filtered = append(filtered, a)
		}
	}
	os.Args = append([]string{os.Args[0]}, filtered...)

	addr := flag.String("addr", ":8080", "listen address")
	dbPath := flag.String("db", "seqannot.db", "sqlite database path")
	flag.Parse()

	if smoke {
		if err := selfcheck.Run(); err != nil {
			fmt.Fprintf(os.Stderr, "smoke-test failed: %v\n", err)
			os.Exit(1)
		}
		fmt.Println("smoke-test OK")
		return
	}

	st, err := store.Open(*dbPath)
	if err != nil {
		log.Fatalf("open store: %v", err)
	}
	defer st.Close()

	svc := service.New(st)
	// On startup, resume any jobs left non-terminal by a previous run.
	if rows, err := svc.ReconcileAll(context.Background()); err != nil {
		log.Printf("reconcile on startup: %v", err)
	} else {
		log.Printf("startup reconcile: %d job(s) reviewed", len(rows))
	}

	mux := httpapi.NewMux(svc)
	srv := &http.Server{
		Addr:              *addr,
		Handler:           mux,
		ReadHeaderTimeout: 10 * time.Second,
	}

	go func() {
		log.Printf("seqannot listening on %s", *addr)
		if err := srv.ListenAndServe(); err != nil && err != http.ErrServerClosed {
			log.Fatalf("listen: %v", err)
		}
	}()

	stop := make(chan os.Signal, 1)
	signal.Notify(stop, os.Interrupt, syscall.SIGTERM)
	<-stop
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()
	_ = srv.Shutdown(ctx)
}
