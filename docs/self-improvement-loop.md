# Iteration Loop ("Self-Improvement") for the Drug Repurposing Scorecard

## What this is not

This is **not** an autonomous machine-learning system that retrains or
changes its own logic. No model in this repo learns from its own output.
"Self-improvement" here means something narrower and honest: a **repeatable,
versioned research cycle** where each run's output is logged, so a human can
see whether adding new evidence changes the ranking — and decide whether to
trust it more or less over time.

Framing it any other way would violate this project's own rules (no
fabricated results, no unverified claims of capability).

## The loop

```
 ┌─────────────┐   ┌──────────────┐   ┌────────────┐   ┌─────────┐   ┌────────┐
 │ 1. Evidence │ → │ 2. Score all │ → │ 3. Go/No-Go│ → │ 4. Log  │ → │ 5. Diff│
 │   sources   │   │  candidates  │   │   check    │   │  run    │   │ vs prev│
 └─────────────┘   └──────────────┘   └────────────┘   └─────────┘   └───┬────┘
        ▲                                                                 │
        └─────────────────── implement next backlog item ◄────────────────┘
```

1. **Evidence sources** — currently: DE strength (`deseq2_analysis.py`),
   druggability (`expand_drug_targets.py`), pathway membership
   (`pathway_enrichment.py`). Backlog items A3–A5 in
   [`drug-repurposing-roadmap.md`](drug-repurposing-roadmap.md) (literature
   co-mention, connectivity mapping, outlier sensitivity) are additional
   evidence sources not yet implemented.
2. **Score** — `scripts/repurposing_scorecard.py` combines whatever evidence
   sources currently exist into one composite score per gene.
3. **Go/No-Go** — the same script checks the fixed rule from the roadmap doc
   (top decile + drug phase ≥3 + significant pathway) and counts how many
   candidates pass.
4. **Log** — every run appends one row to
   `results/tables/scorecard_run_history.tsv`:
   `timestamp_utc, git_commit, n_candidates_scanned, top_gene, top_score, n_go_no_go_pass`.
   This file is append-only — never edit past rows by hand.
5. **Diff** — compare the new row to prior rows. Two possible outcomes:
   - **Stable**: `top_gene` and `n_go_no_go_pass` don't change as new evidence
     is added → the candidate ranking is robust to the evidence sources
     tried so far.
   - **Drift**: the top candidate changes, or previously-passing candidates
     stop passing → treat this as a signal to investigate *why*, not as a
     failure. Record the reason in the roadmap doc before trusting either
     ranking.

## When to run the loop

- After implementing any backlog item (A3/A4/A5) in the roadmap doc.
- After re-running the discovery pipeline against a new or updated cohort.
- Before citing a specific gene/drug pair as a "top candidate" anywhere
  (preprint, README, external communication) — check the latest history row
  first.

## What this loop deliberately does not do

- It does not auto-tune `WEIGHT_DE` / `WEIGHT_DRUGGABILITY` / `WEIGHT_PATHWAY`
  in `repurposing_scorecard.py`. Weight changes require a documented
  rationale and a new commit, per the comment in that file — never a silent,
  automatic adjustment.
- It does not delete or overwrite history rows. If a run was made in error,
  add a new corrective row rather than editing history.
- It does not promote a candidate to Phase B (preclinical) automatically.
  That decision stays with a human, per the roadmap doc's explicit non-goals.
