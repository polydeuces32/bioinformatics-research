# Drug Repurposing Roadmap — MMP7 / Doxycycline (Option 2)

Status: **in silico hypothesis only**. Nothing in this document or the linked
scripts constitutes clinical or preclinical evidence. All steps below are
either (a) implemented in this repo as additional computational evidence, or
(b) explicitly out of repo scope and require external wet-lab/clinical
partners.

## Current state

- `scripts/drug_targets.py` queries OpenTargets for a **curated list of 14**
  signature genes and finds MMP7 → doxycycline (WHO Essential Medicines List,
  max clinical phase = approved, but for non-oncology indications).
- This is a single-target, single-source lookup. It has not been
  cross-checked against expression magnitude, pathway centrality, or
  independent literature evidence.

## Phase A — Strengthen the in silico case (this repo, no wet lab)

Goal: turn "one gene we happened to check has an approved drug" into a
ranked, reproducible, falsifiable candidate list.

| # | Feature | Script | Status |
|---|---|---|---|
| A1 | Scan **all** significant DEGs for druggability (not just 14 curated genes) | `scripts/expand_drug_targets.py` | Added |
| A2 | Composite repurposing score: DE magnitude + FDR + druggability + pathway membership | `scripts/repurposing_scorecard.py` | Added |
| A3 | Literature co-mention evidence (PubMed E-utilities: gene + drug + "colorectal cancer") | Not implemented | Backlog |
| A4 | Connectivity-mapping check (LINCS L1000 signature reversal) against public CMap-style API | Not implemented | Backlog |
| A5 | Sensitivity check: does the candidate ranking survive removing the top 1–2 outlier genes? | Not implemented | Backlog |
| A6 | Version each scorecard run and track ranking stability/drift as A3–A5 land | `scripts/repurposing_scorecard.py` (history log) + [`self-improvement-loop.md`](self-improvement-loop.md) | Added |

Go/no-go for Phase A → Phase B: candidate must (1) rank in the top decile of
the composite score, (2) have an approved or Phase 3+ drug, (3) sit in a
pathway independently flagged as significant in `pathway_enrichment.py`.
MMP7/doxycycline currently satisfies (2); (1) and (3) are now checked
automatically by `repurposing_scorecard.py`.

## Phase B — Preclinical validation (out of repo scope)

Requires a wet-lab partner. Not executable from this codebase.

- CRC cell-line dose-response (MMP7 activity assay + proliferation/invasion)
  under doxycycline treatment vs. vehicle control.
- Organoid or xenograft validation if cell-line data is positive.
- Mechanistic confirmation that doxycycline's effect is MMP7-mediated
  (e.g. MMP7 knockdown/rescue) rather than an off-target antibiotic effect.

## Phase C — Translational / clinical (out of repo scope)

- Retrospective cohort analysis: MMP7 expression vs. outcomes in patients
  incidentally exposed to tetracyclines (requires clinical data access +
  IRB).
- If Phase B is positive: investigator-initiated Phase 2 repurposing trial.
- Regulatory: repurposing an approved, off-patent drug follows a 505(b)(2)-style
  pathway (US) rather than a full new-drug IND, which lowers cost/time
  relative to a novel compound — but still requires a sponsor and IRB-approved
  trial.

## Explicit non-goals

- This repo will not fabricate wet-lab or clinical data.
- No step here claims doxycycline is effective against CRC in patients.
- Any script producing a "score" or "rank" is a prioritization aid for human
  researchers, not a clinical claim.
