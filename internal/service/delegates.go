package service

import (
	"context"

	"task136-seqannot/internal/store"
)

// The following methods are thin proxies to the store. They exist on the
// service layer (rather than having httpapi call the store directly) so that
// the service remains the single dependency of the HTTP layer and so any
// future authorization, caching, or validation hook has one place to live.

func (svc *Service) ListSequences(ctx context.Context) ([]store.SeqRow, error) {
	return svc.store.ListSequences(ctx)
}
func (svc *Service) GetSequence(ctx context.Context, id string) (store.SeqRow, error) {
	return svc.store.GetSequence(ctx, id)
}
func (svc *Service) DeleteSequence(ctx context.Context, id string) error {
	return svc.store.DeleteSequence(ctx, id)
}

func (svc *Service) ListJobs(ctx context.Context) ([]store.JobRow, error) {
	return svc.store.ListJobs(ctx)
}
func (svc *Service) GetJob(ctx context.Context, id string) (store.JobRow, error) {
	return svc.store.GetJob(ctx, id)
}

func (svc *Service) ListMotifs(ctx context.Context) ([]store.MotifRow, error) {
	return svc.store.ListMotifs(ctx)
}
func (svc *Service) DeleteMotif(ctx context.Context, id string) error {
	return svc.store.DeleteMotif(ctx, id)
}

func (svc *Service) ListEnzymes(ctx context.Context) ([]store.EnzymeRow, error) {
	return svc.store.ListEnzymes(ctx)
}
func (svc *Service) DeleteEnzyme(ctx context.Context, id string) error {
	return svc.store.DeleteEnzyme(ctx, id)
}

func (svc *Service) GetAlignment(ctx context.Context, id string) (store.AlignmentRow, error) {
	return svc.store.GetAlignment(ctx, id)
}
