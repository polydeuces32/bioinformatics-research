"""
Systematic druggability scan across all significant DEGs (Phase A1 of
docs/drug-repurposing-roadmap.md).

Generalizes scripts/drug_targets.py, which only checks a hand-picked list of
14 genes, to the full significant DEG set produced by deseq2_analysis.py +
annotate_results.py. This avoids cherry-picking: every significant gene gets
the same druggability check, ranked by DE evidence.

Usage:
    python scripts/expand_drug_targets.py --top-n 100 --direction both
"""
import argparse
import time

import pandas as pd
import requests

SIG_PATH = "results/tables/deseq2_significant_annotated.tsv"
OUT_PATH = "results/tables/drug_targets_expanded.tsv"

ENSEMBL_XREF_URL = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

OT_QUERY = """
query targetDrugs($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    approvedName
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          name
          maximumClinicalStage
        }
        diseases {
          disease { name }
        }
      }
    }
  }
}
"""

PHASE_MAP = {
    "PHASE_1": 1,
    "PHASE_2": 2,
    "PHASE_3": 3,
    "PHASE_4": 4,
    "APPROVAL": 4,
    "PRECLINICAL": 0,
}

REQUEST_TIMEOUT = 15
MAX_RETRIES = 3
RETRY_BACKOFF_S = 2


def request_with_retry(method, url, **kwargs):
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.request(method, url, timeout=REQUEST_TIMEOUT, **kwargs)
            r.raise_for_status()
            return r
        except (requests.RequestException,) as e:
            last_err = e
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF_S * attempt)
    raise last_err


def symbol_to_ensembl(symbol: str) -> str | None:
    try:
        r = request_with_retry(
            "GET",
            ENSEMBL_XREF_URL.format(symbol=symbol),
            headers={"Content-Type": "application/json"},
            params={"object_type": "gene"},
        )
        rows = r.json()
        for row in rows:
            if row.get("id", "").startswith("ENSG"):
                return row["id"]
    except Exception as e:
        print(f"  {symbol}: Ensembl lookup failed — {e}")
    return None


def query_opentargets(ensembl_id: str) -> dict:
    r = request_with_retry(
        "POST",
        OT_URL,
        json={"query": OT_QUERY, "variables": {"ensemblId": ensembl_id}},
    )
    return r.json()


def main():
    parser = argparse.ArgumentParser(
        description="Scan all significant DEGs for druggability via OpenTargets."
    )
    parser.add_argument(
        "--top-n", type=int, default=100,
        help="Max number of DEGs to scan, ranked by padj (default: 100)."
    )
    parser.add_argument(
        "--direction", choices=["up", "down", "both"], default="both",
        help="Restrict to upregulated, downregulated, or both directions."
    )
    args = parser.parse_args()

    sig = pd.read_csv(SIG_PATH, sep="\t", index_col=0)
    sig = sig.dropna(subset=["gene_symbol"])
    sig = sig[sig["gene_symbol"] != "unknown"]

    if args.direction == "up":
        sig = sig[sig["log2FoldChange"] > 0]
    elif args.direction == "down":
        sig = sig[sig["log2FoldChange"] < 0]

    sig = sig.sort_values("padj").drop_duplicates(subset="gene_symbol").head(args.top_n)
    print(f"Scanning {len(sig)} significant DEGs (direction={args.direction})...\n")

    results = []
    for _, row in sig.iterrows():
        symbol = row["gene_symbol"]
        ensembl_id = symbol_to_ensembl(symbol)
        if ensembl_id is None:
            print(f"{symbol}: no Ensembl ID found, skipping")
            continue

        try:
            data = query_opentargets(ensembl_id)
        except Exception as e:
            print(f"{symbol}: OpenTargets query failed — {e}")
            continue

        if "errors" in data:
            print(f"{symbol}: GraphQL error — {data['errors'][0]['message'][:100]}")
            continue

        target = data.get("data", {}).get("target") or {}
        if not target:
            continue

        drug_data = target.get("drugAndClinicalCandidates") or {}
        drug_rows = drug_data.get("rows") or []

        max_phase = 0
        approved_drugs = []
        for r in drug_rows:
            phase = PHASE_MAP.get(r.get("maxClinicalStage") or "", 0)
            drug = r.get("drug") or {}
            is_approved = str(drug.get("maximumClinicalStage") or "") == "APPROVAL"
            if phase > max_phase:
                max_phase = phase
            if is_approved:
                diseases = r.get("diseases") or []
                disease = (diseases[0].get("disease") or {}).get("name", "") if diseases else ""
                approved_drugs.append(f"{drug.get('name', '')} ({disease})")

        results.append({
            "gene": symbol,
            "log2FoldChange": row["log2FoldChange"],
            "padj": row["padj"],
            "ensembl_id": ensembl_id,
            "total_drugs": drug_data.get("count", 0),
            "max_clinical_phase": max_phase,
            "approved_drugs": "; ".join(approved_drugs[:3]) if approved_drugs else "none",
        })
        print(f"{symbol}: {drug_data.get('count', 0)} drugs, max_phase={max_phase}")

    if not results:
        print("\nNo druggable targets found in scanned set.")
        return

    out_df = pd.DataFrame(results).sort_values(
        ["max_clinical_phase", "padj"], ascending=[False, True]
    )
    print(f"\n=== Druggable targets found: {(out_df['total_drugs'] > 0).sum()} / {len(out_df)} ===")
    print(out_df.head(20).to_string(index=False))

    out_df.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nSaved: {OUT_PATH}")


if __name__ == "__main__":
    main()
