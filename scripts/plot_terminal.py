import pandas as pd
import plotext as plt

sig = pd.read_csv("results/tables/deseq2_significant_annotated.tsv", sep="\t", index_col=0)

# Top 10 upregulated + top 10 downregulated
top_up = sig[sig["log2FoldChange"] > 1].sort_values("padj").head(10)
top_dn = sig[sig["log2FoldChange"] < -1].sort_values("padj").head(10)
combined = pd.concat([top_up, top_dn]).sort_values("log2FoldChange", ascending=False)

genes = combined["gene_symbol"].tolist()
lfc = combined["log2FoldChange"].tolist()

colors = ["red+" if v > 2 else "red" if v > 1 else "green+" if v < -2 else "green" for v in lfc]

plt.bar(genes, lfc, orientation="horizontal", color=colors)
plt.title("Top DEGs | RED = upregulated in tumor | GREEN = downregulated in tumor")
plt.xlabel("log2FC (tumor vs normal)")
plt.plotsize(100, 35)
plt.theme("dark")
plt.show()
