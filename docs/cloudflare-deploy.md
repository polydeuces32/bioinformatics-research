# Deploying the Dashboard to Cloudflare

## Why Containers, not Pages

Cloudflare Pages only serves static files — it cannot run Streamlit, which
requires a persistent Python process with WebSocket connections for its
reactive UI. [Cloudflare Containers](https://developers.cloudflare.com/containers/)
runs the actual Streamlit server inside a container, fronted by a Cloudflare
Worker that proxies requests (including the WebSocket upgrade) to it.

## What's here

```
cloudflare/
  wrangler.jsonc   — Worker + Container config
  src/index.ts     — Worker: routes every request to one Container instance
  package.json     — wrangler, @cloudflare/containers, typescript
  tsconfig.json
dashboard/
  Dockerfile       — builds the Streamlit image (build context = repo root)
  requirements.txt — lean deps (streamlit, plotly, pandas, numpy) — no
                      pydeseq2/gseapy/lifelines/GEOparse, since the
                      dashboard only reads results/tables/*.tsv
```

The Dockerfile's build context is the **repository root** (set via
`image_build_context` in `wrangler.jsonc`), so it can `COPY results/tables`
in addition to `dashboard/`. This means **the deployed dashboard is a
snapshot of whatever is in `results/tables/` at build time** — since
`results/` is gitignored, run the pipeline locally first if you want real
data baked in, or leave it empty to deploy the "run the pipeline first"
placeholder state (an empty `results/tables/.gitkeep` is committed so the
`COPY` never fails on a fresh clone).

The Worker routes every request — including WebSocket upgrades — to a
**single** Container instance (`getContainer()` with no name), because
Streamlit keeps server-side session state per browser tab. This is a
single-instance research dashboard, not a horizontally-scaled multi-tenant
service. The container sleeps after 10 minutes of inactivity and cold-starts
on the next request.

## Prerequisites

- Docker (or a Docker-compatible CLI + daemon) running locally — Wrangler
  builds the image with it.
- Node.js 18+.
- A Cloudflare account with Workers + Containers enabled (Containers is a
  paid-plan feature; check current pricing/availability on your account).
- A Cloudflare API token with Workers/Containers deploy permissions (or run
  `wrangler login` interactively instead, if deploying from a machine with a
  browser).

## Deploy

```bash
# 1. (optional) populate results/tables/ with real output first
python scripts/qc_counts.py
python scripts/deseq2_analysis.py
python scripts/annotate_results.py
python scripts/pathway_enrichment.py
python scripts/expand_drug_targets.py
python scripts/repurposing_scorecard.py

# 2. install Worker tooling
cd cloudflare
npm install

# 3. authenticate (pick one)
npx wrangler login                       # interactive browser login
# — or —
export CLOUDFLARE_API_TOKEN=...          # for CI / non-interactive environments

# 4. deploy
npx wrangler deploy
```

First deploy takes a few minutes (image build + push + container
provisioning). Subsequent deploys are faster (cached image layers). Check
status with:

```bash
npx wrangler containers list
```

## Redeploying with fresh data

The dashboard doesn't run the pipeline itself — it only serves whatever was
baked into the image. To publish updated results: re-run the pipeline
scripts locally to regenerate `results/tables/`, then `npx wrangler deploy`
again to rebuild and push a new image snapshot.

## Local dev / verification without a Cloudflare account

These don't require Cloudflare credentials:

```bash
# Run the dashboard directly (fastest path)
streamlit run dashboard/app.py

# Or build/run the exact container that would be deployed
docker build -f dashboard/Dockerfile -t crc-dashboard .
docker run -p 8501:8501 crc-dashboard

# Type-check the Worker
cd cloudflare && npm install && npm run typecheck

# Validate the wrangler config + Dockerfile build (stops before pushing)
cd cloudflare && npx wrangler deploy --dry-run
```
