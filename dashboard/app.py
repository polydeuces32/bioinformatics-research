"""
Read-only visualization frontend for the pipeline outputs in results/tables/.

This dashboard does not run any analysis itself — it only renders whatever
has already been produced by the scripts in scripts/. If a results file is
missing, the corresponding section shows instructions instead of failing.

Usage:
    streamlit run dashboard/app.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

RESULTS_DIR = Path("results/tables")

st.set_page_config(
    page_title="CRC Transcriptomics — Results Dashboard",
    layout="wide",
)


def load_tsv(name: str) -> pd.DataFrame | None:
    path = RESULTS_DIR / name
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t")


def missing(name: str, script: str) -> None:
    st.info(f"`{name}` not found. Run `python scripts/{script}` first.")


st.title("Transcriptomic Signature of Colorectal Cancer — Results Dashboard")
st.caption(
    "Dry-lab, hypothesis-generating results only. See the repository README "
    "and docs/drug-repurposing-roadmap.md for scope and limitations."
)

tab_deg, tab_pathway, tab_drug, tab_history = st.tabs(
    ["Differential Expression", "Pathway Enrichment", "Drug Repurposing", "Run History"]
)

# --- Differential Expression ---
with tab_deg:
    all_res = load_tsv("deseq2_all_results.tsv")
    sig = load_tsv("deseq2_significant_annotated.tsv")

    if all_res is None:
        missing("deseq2_all_results.tsv", "deseq2_analysis.py")
    else:
        n_sig = len(sig) if sig is not None else int(
            ((all_res["padj"] < 0.05) & (all_res["log2FoldChange"].abs() > 1)).sum()
        )
        c1, c2, c3 = st.columns(3)
        c1.metric("Genes tested", f"{len(all_res):,}")
        c2.metric("Significant DEGs (FDR<0.05, |log2FC|>1)", f"{n_sig:,}")
        c3.metric(
            "Upregulated / Downregulated",
            f"{(sig['log2FoldChange'] > 0).sum() if sig is not None else '—'} / "
            f"{(sig['log2FoldChange'] < 0).sum() if sig is not None else '—'}",
        )

        plot_df = all_res.dropna(subset=["padj", "log2FoldChange"]).copy()
        plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
        plot_df["significant"] = (plot_df["padj"] < 0.05) & (plot_df["log2FoldChange"].abs() > 1)

        fig = px.scatter(
            plot_df, x="log2FoldChange", y="neg_log10_padj", color="significant",
            color_discrete_map={True: "crimson", False: "lightgray"},
            hover_data=["GeneID"], title="Volcano plot — tumor vs. normal",
            labels={"neg_log10_padj": "-log10(padj)"},
        )
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
        fig.add_vline(x=1, line_dash="dash", line_color="gray")
        fig.add_vline(x=-1, line_dash="dash", line_color="gray")
        st.plotly_chart(fig, width="stretch")

        if sig is not None:
            st.subheader("Top 20 DEGs by adjusted p-value")
            top = sig.sort_values("padj").head(20)[
                ["gene_symbol", "baseMean", "log2FoldChange", "padj"]
            ]
            st.dataframe(top, width="stretch", hide_index=True)
        else:
            missing("deseq2_significant_annotated.tsv", "annotate_results.py")

# --- Pathway Enrichment ---
with tab_pathway:
    up = load_tsv("pathways_upregulated.tsv")
    down = load_tsv("pathways_downregulated.tsv")

    if up is None and down is None:
        missing("pathways_upregulated.tsv / pathways_downregulated.tsv", "pathway_enrichment.py")
    else:
        col_up, col_down = st.columns(2)
        for col, df, label, color in [
            (col_up, up, "Upregulated", "crimson"),
            (col_down, down, "Downregulated", "steelblue"),
        ]:
            with col:
                st.subheader(f"{label} Hallmark pathways (FDR<0.05)")
                if df is None or df.empty:
                    st.write("No significant pathways.")
                    continue
                plot_df = df.head(10).copy()
                plot_df["neg_log10_fdr"] = -np.log10(plot_df["Adjusted P-value"].clip(lower=1e-300))
                fig = px.bar(
                    plot_df.sort_values("neg_log10_fdr"),
                    x="neg_log10_fdr", y="Term", orientation="h",
                    color_discrete_sequence=[color],
                    labels={"neg_log10_fdr": "-log10(FDR)", "Term": ""},
                )
                st.plotly_chart(fig, width="stretch")
                st.dataframe(
                    df[["Term", "Overlap", "Adjusted P-value"]],
                    width="stretch", hide_index=True,
                )

# --- Drug Repurposing ---
with tab_drug:
    scorecard = load_tsv("repurposing_scorecard.tsv")
    expanded = load_tsv("drug_targets_expanded.tsv")

    if scorecard is None:
        if expanded is None:
            missing("drug_targets_expanded.tsv", "expand_drug_targets.py")
        else:
            missing("repurposing_scorecard.tsv", "repurposing_scorecard.py")
    else:
        n_pass = int(
            (scorecard["top_decile"] & (scorecard["max_clinical_phase"] >= 3)
             & scorecard["in_significant_pathway"]).sum()
        )
        st.metric(
            "Candidates passing all 3 go/no-go criteria",
            f"{n_pass} / {len(scorecard)}",
            help="Top decile composite score + drug phase>=3 + significant pathway membership. "
                 "See docs/drug-repurposing-roadmap.md.",
        )

        fig = px.bar(
            scorecard.sort_values("composite_score"),
            x="composite_score", y="gene", orientation="h",
            color="in_significant_pathway",
            hover_data=["max_clinical_phase", "approved_drugs"],
            title="Repurposing candidate scorecard",
        )
        st.plotly_chart(fig, width="stretch")

        st.dataframe(
            scorecard.sort_values("composite_score", ascending=False),
            width="stretch", hide_index=True,
        )
        st.caption(
            "In silico evidence only — no wet-lab or clinical validation. "
            "See docs/drug-repurposing-roadmap.md for the Phase A/B/C plan."
        )

# --- Run History (drift tracking) ---
with tab_history:
    history = load_tsv("scorecard_run_history.tsv")
    if history is None:
        missing("scorecard_run_history.tsv", "repurposing_scorecard.py")
    else:
        st.write(
            "Each row is one run of `repurposing_scorecard.py`. See "
            "docs/self-improvement-loop.md — this is a versioned research log, "
            "not an autonomous learning process."
        )
        st.dataframe(history, width="stretch", hide_index=True)
        if history["top_gene"].nunique() > 1:
            st.warning(
                "Top candidate has changed across runs — ranking drift detected. "
                "Investigate before citing a specific gene externally."
            )
        if len(history) > 1:
            fig = px.line(
                history, x="timestamp_utc", y="top_score", color="top_gene", markers=True,
                title="Top candidate score over successive runs",
            )
            st.plotly_chart(fig, width="stretch")
