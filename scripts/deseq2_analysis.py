import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

np.random.seed(42)

counts = pd.read_csv("data/processed/counts_filtered.tsv", sep="\t", index_col=0).T
metadata = pd.read_csv("data/processed/metadata.tsv", sep="\t", index_col=0)

# Align
counts = counts.loc[metadata.index]
print(f"Input: {counts.shape[0]} samples x {counts.shape[1]} genes")

# Run DESeq2
dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition")
dds.deseq2()

# Extract results (tumor vs normal)
stat_res = DeseqStats(dds, contrast=["condition", "tumor", "normal"])
stat_res.summary()
results = stat_res.results_df

# Filter significant DEGs
sig = results[(results["padj"] < 0.05) & (results["log2FoldChange"].abs() > 1)]
print(f"\nSignificant DEGs (FDR<0.05, |log2FC|>1): {len(sig)}")
print(sig.sort_values("padj").head(10))

# Save
results.to_csv("results/tables/deseq2_all_results.tsv", sep="\t")
sig.sort_values("padj").to_csv("results/tables/deseq2_significant.tsv", sep="\t")
print("\nSaved results.")
