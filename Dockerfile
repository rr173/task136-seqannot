# Build stage: Go 1.26.3 with Go module mode, CGO off.
FROM docker.m.daocloud.io/library/golang:1.26.3-bookworm AS builder

WORKDIR /src

# Download module dependencies first to leverage cache.
COPY go.mod go.sum ./
RUN go mod download

# Copy the rest of the source (including the embedded web/ frontend).
COPY . .

ENV CGO_ENABLED=0 \
    GOTOOLCHAIN=local \
    GOPROXY=https://goproxy.cn,direct \
    GOSUMDB=sum.golang.google.cn

RUN go build -o /out/seqannot .

# Runtime stage: small alpine image with the binary and a writable data dir.
FROM docker.m.daocloud.io/library/alpine:3.20

RUN apk add --no-cache ca-certificates tzdata && \
    addgroup -S app && adduser -S -G app app

WORKDIR /data
COPY --from=builder /out/seqannot /usr/local/bin/seqannot

USER app
ENV SEQANNOT_DB=/data/seqannot.db
EXPOSE 8080

ENTRYPOINT ["/usr/local/bin/seqannot"]
CMD ["--addr", ":8080"]
