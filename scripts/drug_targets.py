import requests
import pandas as pd
import plotext as plt

genes = {
    "MMP7":   "ENSG00000137673",
    "CLDN1":  "ENSG00000163347",
    "CLDN2":  "ENSG00000165376",
    "CDH3":   "ENSG00000062038",
    "TRIP13": "ENSG00000071539",
    "ANLN":   "ENSG00000011426",
    "INHBA":  "ENSG00000122641",
    "ATAD2":  "ENSG00000156802",
    "ETV4":   "ENSG00000175832",
    "SIM2":   "ENSG00000159263",
    "GRIN2D": "ENSG00000105464",
    "KRT80":  "ENSG00000167767",
    "ESM1":   "ENSG00000164283",
    "CPNE7":  "ENSG00000178773",
}

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY = """
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
          drugType
        }
        diseases {
          disease {
            name
          }
        }
      }
    }
  }
}
"""

PHASE_MAP = {
    "Phase I": 1,
    "Phase II": 2,
    "Phase III": 3,
    "Phase IV": 4,
    "Approved": 4,
    "Preclinical": 0,
}

results = []
print(f"Querying OpenTargets for {len(genes)} genes...\n")

for symbol, ensembl_id in genes.items():
    try:
        r = requests.post(
            OT_URL,
            json={"query": QUERY, "variables": {"ensemblId": ensembl_id}},
            timeout=15
        )
        data = r.json()

        if "errors" in data:
            print(f"{symbol}: GraphQL error — {data['errors'][0]['message'][:100]}")
            continue

        target = data.get("data", {}).get("target", {})
        if not target:
            print(f"{symbol}: no data returned")
            continue

        approved_symbol = target.get("approvedSymbol", symbol)
        full_name = target.get("approvedName", "")
        drug_data = target.get("drugAndClinicalCandidates") or {}
        drug_count = drug_data.get("count", 0)
        drug_rows = drug_data.get("rows", []) or []

        max_phase = 0
        approved_drugs = []
        for row in drug_rows:
            phase_raw = row.get("maxClinicalStage") or ""
            phase = PHASE_MAP.get(phase_raw, 0)

            drug = row.get("drug") or {}
            drug_name = drug.get("name", "")
            drug_max_stage = drug.get("maximumClinicalStage") or ""
            is_approved = "IV" in str(drug_max_stage) or str(drug_max_stage).lower() == "approved"

            diseases = row.get("diseases") or []
            disease = ""
            if diseases:
                d = diseases[0].get("disease") or {}
                disease = d.get("name", "")

            if phase > max_phase:
                max_phase = phase
            if is_approved:
                approved_drugs.append(f"{drug_name} ({disease})")

        results.append({
            "gene": approved_symbol,
            "full_name": full_name,
            "total_drugs": drug_count,
            "max_clinical_phase": max_phase,
            "approved_drugs": "; ".join(approved_drugs[:3]) if approved_drugs else "none"
        })
        print(f"{approved_symbol}: {drug_count} drugs, max phase={max_phase}, approved={len(approved_drugs)}")

    except Exception as e:
        print(f"{symbol}: error — {e}")

if not results:
    print("No results returned.")
else:
    res_df = pd.DataFrame(results).sort_values("max_clinical_phase", ascending=False)
    print(f"\n=== Drug Target Summary ===")
    print(res_df.to_string(index=False))

    plot_df = res_df[res_df["total_drugs"] > 0]
    if len(plot_df) > 0:
        gene_names = plot_df["gene"].tolist()
        phases = plot_df["max_clinical_phase"].tolist()
        colors = []
        for p in phases:
            if p >= 4:
                colors.append("green+")
            elif p == 3:
                colors.append("green")
            elif p >= 1:
                colors.append("orange")
            else:
                colors.append("red")

        plt.bar(gene_names, phases, color=colors)
        plt.title("Drug Targetability | green=approved, orange=clinical trial, red=preclinical")
        plt.ylabel("Max Clinical Phase")
        plt.plotsize(118, 43)
        plt.theme("dark")
        plt.show()

    res_df.to_csv("results/tables/drug_targets.tsv", sep="\t", index=False)
    print("\nSaved: results/tables/drug_targets.tsv")
