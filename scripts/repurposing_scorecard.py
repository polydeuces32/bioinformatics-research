"""
Composite drug-repurposing candidate scorecard (Phase A2 of
docs/drug-repurposing-roadmap.md).

Combines three independent lines of in silico evidence for each druggable
gene into one ranked table, so a candidate like MMP7 can be judged against
the full DEG set rather than in isolation:

  1. Differential expression strength  (|log2FoldChange|, padj)
  2. Druggability                       (max clinical phase, approved drugs)
  3. Pathway centrality                 (member of a significant Hallmark pathway)

Requires results/tables/drug_targets_expanded.tsv, which is produced by
scripts/expand_drug_targets.py, and the pathway enrichment tables produced by
scripts/pathway_enrichment.py.

Every run also appends a row to results/tables/scorecard_run_history.tsv
(timestamp, git commit, top candidate, go/no-go count). This is the
iteration-loop mechanism described in docs/self-improvement-loop.md: as
backlog evidence sources (literature mining, connectivity mapping) are
implemented, re-running this script lets you diff successive rows to see
whether the ranking is stable or drifting, instead of comparing runs by hand.

Usage:
    python scripts/repurposing_scorecard.py
"""
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

DRUG_TARGETS_PATH = "results/tables/drug_targets_expanded.tsv"
PATHWAY_UP_PATH = "results/tables/pathways_upregulated.tsv"
PATHWAY_DOWN_PATH = "results/tables/pathways_downregulated.tsv"
OUT_PATH = "results/tables/repurposing_scorecard.tsv"
HISTORY_PATH = "results/tables/scorecard_run_history.tsv"

# Weights are equal by design — no evidence line is assumed more reliable
# than another at the in silico stage. Adjust only with a documented
# rationale (see docs/drug-repurposing-roadmap.md, Phase A go/no-go).
WEIGHT_DE = 1.0
WEIGHT_DRUGGABILITY = 1.0
WEIGHT_PATHWAY = 1.0


def zscore(series: pd.Series) -> pd.Series:
    std = series.std()
    if std == 0 or pd.isna(std):
        return pd.Series(0.0, index=series.index)
    return (series - series.mean()) / std


def git_commit_short() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def log_run_history(out_df: pd.DataFrame, n_go: int) -> None:
    """Append this run's summary to an immutable history log for drift tracking."""
    top = out_df.iloc[0] if len(out_df) else None
    row = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit_short(),
        "n_candidates_scanned": len(out_df),
        "top_gene": top["gene"] if top is not None else "none",
        "top_score": round(float(top["composite_score"]), 4) if top is not None else float("nan"),
        "n_go_no_go_pass": n_go,
    }
    history_path = Path(HISTORY_PATH)
    history_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not history_path.exists()
    pd.DataFrame([row]).to_csv(history_path, sep="\t", mode="a", header=write_header, index=False)
    print(f"Logged run to: {HISTORY_PATH}")


def load_pathway_genes(path: str) -> set:
    try:
        df = pd.read_csv(path, sep="\t")
    except FileNotFoundError:
        print(f"  warning: {path} not found — run pathway_enrichment.py first")
        return set()
    genes = set()
    for gene_str in df.get("Genes", pd.Series(dtype=str)).dropna():
        genes.update(g.strip().upper() for g in gene_str.split(";"))
    return genes


def main():
    drug_df = pd.read_csv(DRUG_TARGETS_PATH, sep="\t")
    drug_df = drug_df[drug_df["total_drugs"] > 0].copy()
    print(f"Loaded {len(drug_df)} druggable candidates from {DRUG_TARGETS_PATH}")

    pathway_genes = load_pathway_genes(PATHWAY_UP_PATH) | load_pathway_genes(PATHWAY_DOWN_PATH)
    print(f"Loaded {len(pathway_genes)} genes across significant Hallmark pathways")

    drug_df["de_score"] = zscore(drug_df["log2FoldChange"].abs()) + zscore(-np.log10(drug_df["padj"].clip(lower=1e-300)))
    drug_df["druggability_score"] = zscore(drug_df["max_clinical_phase"]) + zscore(np.log1p(drug_df["total_drugs"]))
    drug_df["in_significant_pathway"] = drug_df["gene"].str.upper().isin(pathway_genes)
    drug_df["pathway_score"] = drug_df["in_significant_pathway"].astype(float)

    drug_df["composite_score"] = (
        WEIGHT_DE * zscore(drug_df["de_score"])
        + WEIGHT_DRUGGABILITY * zscore(drug_df["druggability_score"])
        + WEIGHT_PATHWAY * zscore(drug_df["pathway_score"])
    )

    out_cols = [
        "gene", "log2FoldChange", "padj", "max_clinical_phase",
        "approved_drugs", "in_significant_pathway", "composite_score",
    ]
    out_df = drug_df[out_cols].sort_values("composite_score", ascending=False)

    top_decile_cutoff = out_df["composite_score"].quantile(0.9)
    out_df["top_decile"] = out_df["composite_score"] >= top_decile_cutoff

    print(f"\n=== Repurposing Candidate Scorecard (top 15) ===")
    print(out_df.head(15).to_string(index=False))

    n_go = int((out_df["top_decile"] & (out_df["max_clinical_phase"] >= 3) & out_df["in_significant_pathway"]).sum())
    print(
        f"\n{n_go} candidate(s) meet all three Phase A go/no-go criteria "
        f"(top decile score, phase>=3 drug, significant pathway membership)."
    )

    out_df.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nSaved: {OUT_PATH}")

    log_run_history(out_df, n_go)


if __name__ == "__main__":
    main()
