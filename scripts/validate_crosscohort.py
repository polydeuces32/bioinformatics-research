import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import plotext as plt

np.random.seed(42)

# --- Load your DEG signature from GSE156451 ---
sig = pd.read_csv("results/tables/deseq2_significant_annotated.tsv", sep="\t", index_col=0)
top_up = sig[sig["log2FoldChange"] > 1].sort_values("padj").head(50)
top_dn = sig[sig["log2FoldChange"] < -1].sort_values("padj").head(50)
signature = pd.concat([top_up, top_dn])
sig_genes = signature["gene_symbol"].tolist()
print(f"Signature genes: {len(sig_genes)} ({len(top_up)} up, {len(top_dn)} down)")

# --- Load GSE50760 FPKM files ---
data_dir = Path("data/raw/GSE50760")
files = sorted(data_dir.glob("*.txt.gz"))

frames = []
for f in files:
    name = f.stem  # e.g. GSM1228184_AMC_2.1_FPKM
    parts = name.split("_")
    # parts: ['GSM1228184', 'AMC', '2.1', 'FPKM']
    sample_id = parts[0]
    suffix = parts[2].split(".")[-1]  # '1', '2', or '3'
    if suffix not in ["1", "2"]:
        continue
    condition = "tumor" if suffix == "1" else "normal"
    df = pd.read_csv(f, sep="\t", index_col=0, header=0)
    df.columns = ["fpkm"]
    df["sample"] = sample_id
    df["condition"] = condition
    frames.append(df)

combined = pd.concat(frames)
matrix = combined.pivot_table(index=combined.index, columns="sample", values="fpkm")
meta50760 = combined[["sample","condition"]].drop_duplicates().set_index("sample")

print(f"GSE50760 matrix: {matrix.shape}")
print(f"Conditions: {meta50760['condition'].value_counts().to_dict()}")

# --- Test directional concordance ---
results = []
for _, row in signature.iterrows():
    gene = row["gene_symbol"]
    direction = "up" if row["log2FoldChange"] > 0 else "down"

    if gene not in matrix.index:
        continue

    tumor_samples = meta50760[meta50760["condition"] == "tumor"].index
    normal_samples = meta50760[meta50760["condition"] == "normal"].index

    tumor_expr = matrix.loc[gene, tumor_samples].dropna()
    normal_expr = matrix.loc[gene, normal_samples].dropna()

    if len(tumor_expr) < 5 or len(normal_expr) < 5:
        continue

    median_diff = tumor_expr.median() - normal_expr.median()
    validated_direction = "up" if median_diff > 0 else "down"
    concordant = direction == validated_direction

    stat, pval = stats.mannwhitneyu(tumor_expr, normal_expr, alternative="two-sided")

    results.append({
        "gene": gene,
        "discovery_lfc": row["log2FoldChange"],
        "discovery_direction": direction,
        "validation_median_diff": median_diff,
        "validation_direction": validated_direction,
        "concordant": concordant,
        "mannwhitney_p": pval
    })

res_df = pd.DataFrame(results)
concordance_rate = res_df["concordant"].mean()
n_concordant = res_df["concordant"].sum()

print(f"\n=== Cross-cohort Validation ===")
print(f"Genes tested:     {len(res_df)}")
print(f"Concordant:       {n_concordant} ({concordance_rate:.1%})")
print(f"Discordant:       {len(res_df) - n_concordant}")

# Binomial test: is concordance > 50% (random)?
binom = stats.binomtest(n_concordant, len(res_df), p=0.5, alternative="greater")
print(f"Binomial p-value: {binom.pvalue:.2e}  (H0: concordance = random)")

# --- Terminal plot ---
plt.bar(
    ["Concordant", "Discordant"],
    [n_concordant, len(res_df) - n_concordant],
    color=["green+", "red+"]
)
plt.title(f"Cross-cohort Validation — GSE156451 → GSE50760 ({concordance_rate:.1%} concordance)")
plt.plotsize(60, 20)
plt.theme("dark")
plt.show()

res_df.to_csv("results/tables/crosscohort_validation.tsv", sep="\t", index=False)
print("\nSaved: results/tables/crosscohort_validation.tsv")
