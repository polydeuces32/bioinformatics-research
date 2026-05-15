import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import plotext as plt

np.random.seed(42)

# --- Load data ---
print("Loading TCGA expression...")
expr = pd.read_csv("data/raw/TCGA_COAD_expression.gz", sep="\t", index_col=0)
print(f"Expression matrix: {expr.shape}")

print("Loading clinical data...")
clin = pd.read_csv("data/raw/TCGA_COAD_clinical.tsv", sep="\t", index_col=0)
print(f"Clinical matrix: {clin.shape}")

# --- Load your validated signature ---
sig = pd.read_csv("results/tables/deseq2_significant_annotated.tsv", sep="\t", index_col=0)
top_up = sig[sig["log2FoldChange"] > 1].sort_values("padj").head(20)["gene_symbol"].dropna().tolist()

# Find signature genes present in TCGA
available_genes = [g for g in top_up if g in expr.index]
print(f"\nSignature genes found in TCGA: {len(available_genes)}/{len(top_up)}")
print(f"Genes: {available_genes}")

# --- Use tumor samples only ---
tumor_samples = [s for s in expr.columns if s.endswith("-01")]
expr_tumor = expr[tumor_samples]

# Match directly — both use full 15-char IDs
common = expr_tumor.columns.intersection(clin.index)
expr_tumor = expr_tumor[common]
clin_sub = clin.loc[common].copy()
print(f"Matched tumor samples: {len(common)}")

# --- Compute signature score (mean z-score of top upregulated genes) ---
gene_expr = expr_tumor.loc[available_genes].T
gene_expr_z = (gene_expr - gene_expr.mean()) / gene_expr.std()
clin_sub["signature_score"] = gene_expr_z.mean(axis=1)

# --- Build survival columns ---
clin_sub["days_to_death"] = pd.to_numeric(clin_sub["days_to_death"], errors="coerce")
clin_sub["days_to_last_followup"] = pd.to_numeric(clin_sub["days_to_last_followup"], errors="coerce")
clin_sub["OS_days"] = clin_sub["days_to_death"].fillna(clin_sub["days_to_last_followup"])
clin_sub["OS_event"] = (clin_sub["vital_status"] == "DECEASED").astype(int)

# Drop missing
surv = clin_sub[["OS_days", "OS_event", "signature_score"]].dropna()
surv = surv[surv["OS_days"] > 0].copy()
print(f"Samples for survival analysis: {len(surv)}")
print(f"Events (deaths): {surv['OS_event'].sum()}")

# --- Cox PH model ---
print("\nFitting Cox PH model...")
cph = CoxPHFitter()
cph.fit(surv, duration_col="OS_days", event_col="OS_event")
print("\n=== Cox Proportional Hazards Model ===")
cph.print_summary()

# --- Kaplan-Meier: high vs low signature ---
median_score = surv["signature_score"].median()
surv["group"] = surv["signature_score"].apply(lambda x: "High" if x >= median_score else "Low")

high = surv[surv["group"] == "High"]
low  = surv[surv["group"] == "Low"]

lr = logrank_test(
    high["OS_days"], low["OS_days"],
    event_observed_A=high["OS_event"],
    event_observed_B=low["OS_event"]
)

print(f"\n=== Kaplan-Meier Log-rank Test ===")
print(f"High signature n={len(high)}, deaths={high['OS_event'].sum()}")
print(f"Low  signature n={len(low)},  deaths={low['OS_event'].sum()}")
print(f"Log-rank p-value: {lr.p_value:.4f}")

# --- Terminal KM plot ---
kmf_high = KaplanMeierFitter()
kmf_low  = KaplanMeierFitter()
kmf_high.fit(high["OS_days"], high["OS_event"], label="High")
kmf_low.fit(low["OS_days"],  low["OS_event"],  label="Low")

t_high = kmf_high.survival_function_.index.tolist()
s_high = kmf_high.survival_function_["High"].tolist()
t_low  = kmf_low.survival_function_.index.tolist()
s_low  = kmf_low.survival_function_["Low"].tolist()

plt.plot(t_high, s_high, color="red+", label="High signature")
plt.plot(t_low,  s_low,  color="green+", label="Low signature")
plt.title(f"Kaplan-Meier — TCGA-COAD | Log-rank p={lr.p_value:.4f}")
plt.xlabel("Days")
plt.ylabel("Survival probability")
plt.plotsize(118, 43)
plt.theme("dark")
plt.show()

# --- Save ---
surv.to_csv("results/tables/survival_results.tsv", sep="\t")
print("\nSaved: results/tables/survival_results.tsv")
