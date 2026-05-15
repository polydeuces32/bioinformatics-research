import pandas as pd
import requests

sig = pd.read_csv("results/tables/deseq2_significant.tsv", sep="\t", index_col=0)

# Query Ensembl REST API to convert Entrez IDs to gene symbols
def entrez_to_symbol(entrez_ids, chunk_size=200):
    symbols = {}
    ids = [str(i) for i in entrez_ids]
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {"db": "gene", "id": ",".join(chunk), "retmode": "json"}
        r = requests.get(url, params=params, timeout=30)
        data = r.json().get("result", {})
        for eid in chunk:
            info = data.get(eid, {})
            symbols[int(eid)] = info.get("name", "unknown")
        print(f"Annotated {min(i+chunk_size, len(ids))}/{len(ids)}")
    return symbols

print("Fetching gene symbols...")
symbol_map = entrez_to_symbol(sig.index.tolist())
sig["gene_symbol"] = sig.index.map(symbol_map)

# Show top 20
top = sig.sort_values("padj").head(20)[["gene_symbol","baseMean","log2FoldChange","padj"]]
print("\nTop 20 DEGs:")
print(top.to_string())

sig.to_csv("results/tables/deseq2_significant_annotated.tsv", sep="\t")
print("\nSaved annotated results.")
