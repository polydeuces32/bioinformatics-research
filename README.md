# Transcriptomic Signature of Colorectal Cancer

A reproducible, replication-first analysis of colorectal cancer transcriptomics using only public data and open-source tools. Identifies a validated 90-gene signature and a drug repurposing hypothesis (MMP7-doxycycline).

**Author:** Giancarlo Vizhnay (Independent Researcher, NYC)
**Preprint:** [`results/preprint_draft.md`](results/preprint_draft.md)

---

## Key Findings

- **5,955 differentially expressed genes** between tumor and normal colon (GSE156451, n=95)
- **100% directional concordance** in independent validation cohort (GSE50760, binomial p = 8.08e-28)
- **Canonical CRC pathways activated:** E2F, MYC, G2-M, Wnt/β-catenin, EMT
- **MMP7 identified as drug repurposing target** — inhibited by doxycycline (WHO Essential Medicines List)
- **Signature is diagnostic, not prognostic** (null TCGA-COAD survival association — honest reporting of negative results)

---

## Repository Structure

---

## Reproducing the Analysis

### Setup

```bash
conda env create -f environment.yml
conda activate bioinformatics-research
```

### Download Data

```bash
# Discovery cohort (raw counts from NCBI pipeline)
curl -L "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE156451&format=file&file=GSE156451_raw_counts_GRCh38.p13_NCBI.tsv.gz" \
  -o data/raw/GSE156451_counts.tsv.gz

# Validation cohort (FPKM)
curl -L "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE50760&format=file" \
  -o data/raw/GSE50760_RAW.tar
mkdir -p data/raw/GSE50760
tar -xf data/raw/GSE50760_RAW.tar -C data/raw/GSE50760/

# TCGA-COAD survival data
curl -L "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.COAD.sampleMap%2FHiSeqV2.gz" \
  -o data/raw/TCGA_COAD_expression.gz
curl -L "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix" \
  -o data/raw/TCGA_COAD_clinical.tsv
```

### Run the Pipeline

```bash
python scripts/qc_counts.py
python scripts/deseq2_analysis.py
python scripts/annotate_results.py
python scripts/validate_crosscohort.py
python scripts/pathway_enrichment.py
python scripts/survival_analysis.py
python scripts/survival_per_gene.py
python scripts/drug_targets.py
```

All intermediate results are written to `results/tables/`.

---

## Data Sources

| Dataset | Source | Access |
|---|---|---|
| GSE156451 (discovery) | NCBI GEO | Open |
| GSE50760 (validation) | NCBI GEO | Open |
| TCGA-COAD (survival) | UCSC Xena | Open |
| OpenTargets v4 (drug targets) | GraphQL API | Open |

No restricted-access data is used.

---

## Reproducibility

- Random seed fixed globally (`seed=42`)
- All package versions pinned in `environment.yml`
- All data downloaded via scripted commands (no manual GUI steps)
- Terminal plots via `plotext` (no figure file dependencies)

---

## Limitations

This is a dry-lab (in silico) analysis. Findings are hypothesis-generating only. Causal inference, mechanism, and clinical utility require functional experiments.

---

## License

MIT — code free to use, modify, redistribute. Data citations belong to original submitters.

---

## Citation

If you use this work, please cite the preprint (bioRxiv submission pending).
