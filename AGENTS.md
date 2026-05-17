# Agents

See `CLAUDE.md` for project goals, stack, and rules.

## Cursor Cloud specific instructions

### Environment

- Python 3.12 with pip (no conda). Dependencies are installed via pip by the update script.
- `~/.local/bin` must be on PATH for CLI tools (`snakemake`, `jupyter`, `gseapy`). Run `export PATH="$HOME/.local/bin:$PATH"` if needed.

### Running the pipeline

The pipeline is a sequence of Python scripts in `scripts/`. See `README.md` "Run the Pipeline" for the full ordered list. Key points:

- **Pre-computed data** exists in `data/processed/` (`counts_filtered.tsv`, `metadata.tsv`), so `deseq2_analysis.py` can run without downloading raw data.
- Raw data download (`download_geo.py` and the README `curl` commands) requires internet and can take a while. The `data/raw/` directory must be created first.
- Output directories `results/tables/`, `results/figures/`, and `logs/` must exist before running scripts. Create them with `mkdir -p results/tables results/figures logs`.
- Several scripts call external APIs (NCBI eUtils, Enrichr, OpenTargets) — these are unauthenticated public endpoints requiring no API keys.
- `plotext` renders terminal-based plots inline; `matplotlib` saves to `results/figures/`. Both work headless.

### Linting and testing

- No formal test suite or linter configuration exists in this repo. Use `python3 -c "import <module>"` to verify library availability.
- The `scripts/test_api.py` script tests OpenTargets GraphQL API connectivity.
- To validate the core pipeline, run `python3 scripts/deseq2_analysis.py` (uses pre-computed data, takes ~15s).
