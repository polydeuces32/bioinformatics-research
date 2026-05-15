import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(42)

COUNTS_PATH = "data/raw/GSE156451_counts.tsv.gz"

counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)

cols = counts.columns.tolist()
labels = ["tumor" if int(c.replace("GSM","")) <= 4731744 else "normal" for c in cols]
metadata = pd.DataFrame({"sample": cols, "condition": labels}).set_index("sample")

# 1. Library sizes
lib_sizes = counts.sum(axis=0)
print("=== Library Sizes ===")
print(f"Min:    {lib_sizes.min():,}")
print(f"Max:    {lib_sizes.max():,}")
print(f"Median: {lib_sizes.median():,.0f}")

# 2. Zero-count genes
zero_genes = (counts.sum(axis=1) == 0).sum()
print(f"\n=== Zero-count genes: {zero_genes} / {len(counts)} ===")

# 3. Filter low-count genes (standard: keep genes with >10 counts in at least 10 samples)
keep = (counts > 10).sum(axis=1) >= 10
counts_filtered = counts[keep]
print(f"\n=== After filtering: {counts_filtered.shape[0]} genes remain ===")

# 4. Plot library sizes
fig, ax = plt.subplots(figsize=(14, 4))
colors = ["tomato" if l == "tumor" else "steelblue" for l in metadata["condition"]]
ax.bar(range(len(lib_sizes)), lib_sizes.values / 1e6, color=colors)
ax.set_xlabel("Sample")
ax.set_ylabel("Library size (millions)")
ax.set_title("GSE156451 Library Sizes — red=tumor, blue=normal")
plt.tight_layout()
plt.savefig("results/figures/library_sizes.png", dpi=150)
print("\nSaved: results/figures/library_sizes.png")

# Save filtered counts and metadata
counts_filtered.to_csv("data/processed/counts_filtered.tsv", sep="\t")
metadata.to_csv("data/processed/metadata.tsv", sep="\t")
print("Saved: data/processed/counts_filtered.tsv")
print("Saved: data/processed/metadata.tsv")
