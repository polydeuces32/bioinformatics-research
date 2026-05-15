import numpy as np
import pandas as pd
import gseapy as gp
import plotext as plt

np.random.seed(42)

# Load full DESeq2 results — ranked by stat for GSEA
results = pd.read_csv("results/tables/deseq2_all_results.tsv", sep="\t", index_col=0)
sig = pd.read_csv("results/tables/deseq2_significant_annotated.tsv", sep="\t", index_col=0)

# Build gene symbol list for enrichment (upregulated only)
up_genes = sig[sig["log2FoldChange"] > 1].sort_values("padj")["gene_symbol"].dropna().tolist()
dn_genes = sig[sig["log2FoldChange"] < -1].sort_values("padj")["gene_symbol"].dropna().tolist()

print(f"Upregulated genes: {len(up_genes)}")
print(f"Downregulated genes: {len(dn_genes)}")

# Run enrichment against Hallmark gene sets
print("\nRunning enrichment — upregulated genes...")
enr_up = gp.enrichr(
    gene_list=up_genes,
    gene_sets=["MSigDB_Hallmark_2020"],
    organism="human",
    outdir=None,
    verbose=False
)

print("Running enrichment — downregulated genes...")
enr_dn = gp.enrichr(
    gene_list=dn_genes,
    gene_sets=["MSigDB_Hallmark_2020"],
    organism="human",
    outdir=None,
    verbose=False
)

# Filter significant pathways
up_res = enr_up.results[enr_up.results["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
dn_res = enr_dn.results[enr_dn.results["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")

print(f"\n=== Upregulated Pathways (FDR<0.05): {len(up_res)} ===")
print(up_res[["Term","Overlap","Adjusted P-value"]].head(10).to_string())

print(f"\n=== Downregulated Pathways (FDR<0.05): {len(dn_res)} ===")
print(dn_res[["Term","Overlap","Adjusted P-value"]].head(10).to_string())

# Terminal plot — top 10 up pathways
if len(up_res) > 0:
    top_paths = up_res.head(10)
    names = [t.replace("HALLMARK_","").replace("_"," ")[:30] for t in top_paths["Term"]]
    pvals = (-np.log10(top_paths["Adjusted P-value"])).tolist()

    plt.bar(names, pvals, orientation="horizontal", color="red+")
    plt.title("Upregulated Hallmark Pathways — -log10(FDR)")
    plt.plotsize(100, 25)
    plt.theme("dark")
    plt.show()

# Save
up_res.to_csv("results/tables/pathways_upregulated.tsv", sep="\t", index=False)
dn_res.to_csv("results/tables/pathways_downregulated.tsv", sep="\t", index=False)
print("\nSaved pathway results.")
