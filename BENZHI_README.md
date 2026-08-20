# seqannot — DNA 序列分析与注释引擎

本项目为生物信息分析平台提供 DNA 序列结构化解析与特征注释的后端引擎。主要输入是 DNA 序列（裸序列或 FASTA），主要输出是六读码框翻译、开放阅读框（ORF）、IUPAC 基序检索命中、限制酶切点与酶切片段、双序列比对结果，以及可恢复的分析任务与持久化特征注释。所有数据持久化到本地 SQLite，进程重启后能恢复全部状态并把未完成的分析任务继续跑完。

## 本地命令

```bash
go build ./...              # 编译
go run . --addr :8080       # 启动服务（浏览器访问 http://localhost:8080）
go test ./...               # 全部测试
go run . --smoke-test       # 健康自检（端到端 12 个场景）
```

## 标准本地验证链

```bash
GOTOOLCHAIN=local GOPROXY=https://goproxy.cn,direct GOSUMDB=sum.golang.google.cn go mod vendor
CGO_ENABLED=0 GOTOOLCHAIN=local go test -mod=vendor ./...
CGO_ENABLED=0 GOTOOLCHAIN=local go vet -mod=vendor ./...
CGO_ENABLED=0 GOTOOLCHAIN=local go build -mod=vendor ./...
CGO_ENABLED=0 GOTOOLCHAIN=local go run -mod=vendor . --smoke-test
```

## 启动与页面

```bash
go run . --addr :8080 --db seqannot.db
```

浏览器打开 `http://localhost:8080/` 即可访问前端操作页面：提交序列（裸序列或 FASTA）→ 选择分析项（组成/翻译/ORF/motif/限制酶/酶切）→ 查看结果；可管理 motif 库与限制酶目录；可提交分析任务并执行重启重算（reconcile）。

同时验证页面与业务 API 的 smoke-test：

```bash
go run . --smoke-test
```

该 smoke-test 覆盖：序列 CRUD、GC 组成、六读码框翻译、ORF 检测、motif 检索、限制酶切点与酶切、双序列比对、分析任务提交与重放校验（幂等）、重启恢复、特征注释区间查询、前端页面与业务 API 的联动、错误处理（404/400）。

## 依赖

- Go `1.26.3`（`GOTOOLCHAIN=local`，不写 toolchain 行）。
- 持久化：SQLite，经 `modernc.org/sqlite v1.52.0`（对应 SQLite `3.46.1`），`CGO_ENABLED=0`。无其他第三方依赖。
- 前端：原生 HTML/CSS/JavaScript，无构建依赖（不计入后端代码行数）。

## benzhi 评测 Docker

`build_benzhi_docker.sh <镜像名> <平台>` 通过 `benzhi.Dockerfile` 构建评测镜像：

```bash
bash ./build_benzhi_docker.sh go-task-benzhi:amd64 linux/amd64
docker run -it go-task-benzhi:amd64          # 进入容器
bash ./build_benzhi_docker.sh go-task-benzhi:arm64 linux/arm64
docker run -it go-task-benzhi:arm64
```

`benzhi.Dockerfile` 基于 `golang:1.26.3`，先 `go mod download` 再 `go build ./...`，始终保留 Go 工具链与 `CMD ["bash"]`。

## 双架构交付

主 `Dockerfile` 使用 `docker.m.daocloud.io/library/golang:1.26.3-bookworm` 构建、`alpine:3.20` 运行，`CGO_ENABLED=0`，同时支持 `linux/amd64` 与 `linux/arm64`：

```bash
docker buildx build --platform linux/amd64 --load -t task136-seqannot:amd64 .
docker buildx build --platform linux/arm64 --load -t task136-seqannot:arm64 .
docker run --rm task136-seqannot:amd64 --smoke-test
docker run --rm task136-seqannot:arm64 --smoke-test
```
