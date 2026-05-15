import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import plotext as plt

np.random.seed(42)

# --- Load data ---
expr = pd.read_csv("data/raw/TCGA_COAD_expression.gz", sep="\t", index_col=0)
clin = pd.read_csv("data/raw/TCGA_COAD_clinical.tsv", sep="\t", index_col=0)

# --- Genes to test ---
genes = ["MMP7", "CLDN2", "CLDN1", "CPNE7", "SIM2"]

# --- Use tumor samples only ---
tumor_samples = [s for s in expr.columns if s.endswith("-01")]
expr_tumor = expr[tumor_samples]
common = expr_tumor.columns.intersection(clin.index)
expr_tumor = expr_tumor[common]
clin_sub = clin.loc[common].copy()

# --- Build survival columns ---
clin_sub["days_to_death"] = pd.to_numeric(clin_sub["days_to_death"], errors="coerce")
clin_sub["days_to_last_followup"] = pd.to_numeric(clin_sub["days_to_last_followup"], errors="coerce")
clin_sub["OS_days"] = clin_sub["days_to_death"].fillna(clin_sub["days_to_last_followup"])
clin_sub["OS_event"] = (clin_sub["vital_status"] == "DECEASED").astype(int)

results = []

for gene in genes:
    if gene not in expr_tumor.index:
        print(f"{gene}: not found in expression matrix")
        continue

    clin_sub["expr"] = expr_tumor.loc[gene]
    surv = clin_sub[["OS_days", "OS_event", "expr"]].dropna()
    surv = surv[surv["OS_days"] > 0].copy()

    median_expr = surv["expr"].median()
    surv["group"] = surv["expr"].apply(lambda x: "High" if x >= median_expr else "Low")

    high = surv[surv["group"] == "High"]
    low  = surv[surv["group"] == "Low"]

    lr = logrank_test(
        high["OS_days"], low["OS_days"],
        event_observed_A=high["OS_event"],
        event_observed_B=low["OS_event"]
    )

    results.append({
        "gene": gene,
        "n": len(surv),
        "deaths": surv["OS_event"].sum(),
        "logrank_p": lr.p_value,
        "high_median_survival": high["OS_days"].median(),
        "low_median_survival": low["OS_days"].median()
    })

    print(f"\n{'='*40}")
    print(f"Gene: {gene} | n={len(surv)} | deaths={surv['OS_event'].sum()}")
    print(f"Log-rank p = {lr.p_value:.4f}")
    print(f"High expr median survival: {high['OS_days'].median():.0f} days")
    print(f"Low  expr median survival: {low['OS_days'].median():.0f} days")

# --- Summary table ---
res_df = pd.DataFrame(results).sort_values("logrank_p")
print(f"\n{'='*40}")
print("=== Per-Gene Survival Summary ===")
print(res_df.to_string(index=False))

# --- Terminal plot: p-values ---
genes_plot = res_df["gene"].tolist()
pvals = (-np.log10(res_df["logrank_p"])).tolist()
colors = ["red+" if p >= -np.log10(0.05) else "orange" for p in pvals]

plt.bar(genes_plot, pvals, color=colors)
plt.hline(-np.log10(0.05), color="white")
plt.title("Per-gene survival — -log10(log-rank p) | red = significant")
plt.xlabel("Gene")
plt.ylabel("-log10(p)")
plt.plotsize(118, 43)
plt.theme("dark")
plt.show()

# --- Save ---
res_df.to_csv("results/tables/survival_per_gene.tsv", sep="\t", index=False)
print("\nSaved: results/tables/survival_per_gene.tsv")
