# `predict/` — running every tool against the whole CFTR2 list

`notebooks/00–12` introduce one tool each. This folder does the cross-tool step: take the
**full public CFTR2 variant list as the benchmark**, run every predictor in `toolkit.py`
against it, and report both **coverage** (does a prediction exist at all?) and
**performance** (is it right?).

| Notebook | What it does |
|---|---|
| [`13_cftr2_benchmark.ipynb`](13_cftr2_benchmark.ipynb) | Splits all 2,097 CFTR2 variants into the three parts below, benchmarks every tool at its published threshold, writes `outputs/`. |

## The three parts

The categories are mutually exclusive and ordered — missense first, because SpliceAI scores
essentially every SNV in the region and would otherwise swallow the missense set:

| Part | Rule | n | Share |
|---|---|---|---|
| **1 · missense tools** | ≥1 of AlphaMissense / EVE / ESM1b / REVEL / PrimateAI scores it | 803 | 38.3% |
| **2 · splice tools only** | no missense score, but SpliceAI / Pangolin scores it | 510 | 24.3% |
| **3 · no tool at all** | nothing scores it | 784 | 37.4% |

Weighted by CFTR2 allele frequency instead of variant count, part 3 is **~73%** — F508del
alone is ~65% of CF alleles and is a deletion, invisible to every predictor here.

## Outputs

Written to `../outputs/` (gitignored — regenerated, never committed):

| File | Contents |
|---|---|
| `predict_cftr2_benchmark_ALL.csv` | master table: 2,097 CFTR2 variants × every tool score + part |
| `predict_part1_missense_tools.csv` | part 1 rows |
| `predict_part2_splice_tools_only.csv` | part 2 rows |
| `predict_part3_no_tool.csv` | part 3 rows — the variants nothing can score |
| `predict_tool_coverage.csv` | per-tool: join key, variants scored, % of CFTR2, learning type |
| `predict_tool_performance.csv` | per-tool TP/FP/FN/TN + sensitivity/specificity, naive **and** domain-restricted |
| `predict_category_summary.csv` | counts by part × variant class |
| `predict_cadd_rescue.csv` | live CADD on the part-3 variants CADD can reach (6 of 784) |
| `predict_coverage.png` | coverage bar chart + PrimateAI threshold diagnostic |

## Requirements

Everything here is **REAL data except Pangolin** (`source='DEMO'`: `build_pangolin.py` runs
the real model but only over 5 curated alleles, so its coverage is not benchmark-scale).
The notebook needs these extracts built first — see [`../data/README.md`](../data/README.md):

- `data/cftr2_2026-01-30.csv` (`build_cftr2.py`) — the benchmark itself, required
- `data/eve_cftr_2021-08.csv`, `esm1b_cftr.csv`, `revel_cftr_v1.3.csv`, `primateai_cftr.csv`,
  `spliceai_cftr_2021_v1.3.csv` — one `build_*.py` each
- `_tmp_fetch/alphamissense_cftr.tsv` — cache-backed
- network access for the live CADD v1.7 API (§9; 6 calls, degrades gracefully)

A missing extract leaves that tool's column empty — the notebook still runs, but its
coverage numbers will be wrong. Check the `source` column, and rebuild before treating
anything here as a finding.

```bash
jupyter nbconvert --to notebook --execute --inplace predict/13_cftr2_benchmark.ipynb
```

## Caveats that bound every number here

- **CFTR2 is not independent of ClinVar** — they cross-cite, and REVEL is trained on
  ClinVar-lineage labels, so REVEL's agreement is partly circular (notebook 12).
- **42 negatives.** CFTR2 has ~1,245 CF-causing vs ~42 Non CF-causing variants, so
  specificity is noisy and accuracy is close to meaningless. Rank by sensitivity.
- **Coverage ≠ competence** — REVEL emits scores for 12 nonsense variants; SpliceAI blankets
  missense SNVs it should not be judged on. That is why performance is reported twice, once
  restricted to each tool's actual domain.
