# `predict/` — running every tool against the whole CFTR2 list

`notebooks/00–12` introduce one tool each. This folder does the cross-tool step: take the
**full public CFTR2 variant list as the benchmark**, run every predictor in `toolkit.py`
against it, and report both **coverage** (does a prediction exist at all?) and
**performance** (is it right?).

| Notebook | What it does |
|---|---|
| [`13_cftr2_benchmark.ipynb`](13_cftr2_benchmark.ipynb) | Splits all 2,097 CFTR2 variants three ways, re-scores that split by *mechanism*, benchmarks every tool at its published threshold, writes `outputs/`. |

## Two answers to "is this variant predictable?"

**Does a score exist?** — the three parts, mutually exclusive and ordered (missense first,
because the splice tools between them blanket nearly everything with coordinates):

| Part | Rule | n | Share |
|---|---|---|---|
| **1 · missense tools** | ≥1 of AlphaMissense / EVE / ESM1b / REVEL / PrimateAI scores it | 803 | 38.3% |
| **2 · splice tools only** | no missense score, but SpliceAI / Pangolin scores it | 1,089 | 51.9% |
| **3 · no tool at all** | nothing scores it | 205 | 9.8% |

**Does anything address the mechanism?** — i.e. does some tool speak to *why this variant
breaks CFTR*, rather than merely returning a number for it:

| | no score exists | mechanism unaddressed |
|---|---|---|
| by variant | 205 (9.8%) | **951 (45.4%)** |
| by CF allele frequency | 2.3% | **79.9%** |
| of CF-causing variants | 111 / 1,245 | **788 / 1,245** |

The second column is the honest one. **746 of those 951 carry a score anyway** — a splice
verdict on a frameshift, or a missense score on a stop-gain. Nonsense and frameshift
variants need loss-of-function / NMD reasoning and F508del needs protein folding; this
toolkit models none of that.

**F508del is the whole lesson in one row.** Both splice tools reach it — SpliceAI **0.01**,
Pangolin **0.05** — and both are correct (it does not disrupt splicing) and useless (it
misfolds, and nothing here measures folding). Two independent models agreeing on the right
answer to the wrong question is still no answer.

## Both splice tools now score indels — and what that bought

- **SpliceAI**: `build_spliceai.py` reads Illumina's precomputed **`raw.indel`** release
  alongside the masked SNVs → 1.51 M CFTR indel records, 417 of the 713 CFTR2 indels.
  (`masked.indel` is commonly a 0-byte failed download, so the extract is **mixed
  masked/raw** — rows carry `score_type`; see notebook 09 for the size of that seam.)
- **Pangolin**: no precomputed release exists, so `build_pangolin.py` runs the model over
  every CFTR2 variant with GRCh38 coordinates (~1,892 of 2,097, ~4 min on a GPU) — 575 of
  the 713 indels.

That reach is mostly an illusion. It moved ~580 indels out of "no score" and only **52** into
"mechanism addressed"; on indels the two score **0.04** and **0.05** sensitivity, missing 365
of 380 and 478 of 505 CF-causing ones respectively — correctly, as splice biology. The real
payoff is different: two **independently-built** splice models now measured at scale, and they
agree (*r* ≈ 0.98; of 311 non-coding variants both score, they disagree on 12 at the 0.5 cut).
Agreement between two independent models is stronger evidence than either alone.

## Performance, in domain and out

Each tool at its own published cut-point, against CFTR2 CF-causing vs Non CF-causing:

| Tool | domain | n | sensitivity | specificity |
|---|---|---|---|---|
| AlphaMissense | missense | 292 | 0.938 | 0.818 |
| EVE | missense | 291 | 0.981 | 0.394 |
| ESM1b | missense | 292 | 0.923 | 0.727 |
| REVEL | missense | 292 | 0.907 | 0.667 |
| PrimateAI | missense | 245 | 0.202 | 1.000 |
| SpliceAI | non-coding | 158 | **0.968** | 1.000 |
| Pangolin | non-coding | 158 | **0.948** | 1.000 |
| SpliceAI | *indels (out of domain)* | 381 | 0.039 | 1.000 |
| Pangolin | *indels (out of domain)* | 507 | 0.053 | 1.000 |

Judged on everything they happened to score instead, SpliceAI reads **0.18** and Pangolin
**0.17** — a 20-fold swing from the same models on the same day, driven entirely by which
variants you let into the denominator. **Restrict a tool to its domain before judging it.**

## Outputs

Written to `../outputs/` (gitignored — regenerated, never committed):

| File | Contents |
|---|---|
| `predict_cftr2_benchmark_ALL.csv` | master table: 2,097 CFTR2 variants × every tool score + part + mechanism |
| `predict_part1_missense_tools.csv` | part 1 rows |
| `predict_part2_splice_tools_only.csv` | part 2 rows |
| `predict_part3_no_tool.csv` | part 3 rows — the variants nothing scores |
| `predict_mechanism_unaddressed.csv` | the 951 whose mechanism no tool here models |
| `predict_tool_coverage.csv` | per-tool: join key, variants scored, % of CFTR2, learning type |
| `predict_tool_performance.csv` | per-tool TP/FP/FN/TN + sensitivity/specificity, naive **and** domain-restricted |
| `predict_category_summary.csv` | counts by part × variant class × mechanism |
| `predict_coverage.png` | coverage bars, SpliceAI-vs-Pangolin agreement, PrimateAI threshold diagnostic |

## Requirements

All **REAL data** — including Pangolin, now that it is run at CFTR2 scale. The notebook needs
these extracts built first; see [`../data/README.md`](../data/README.md):

- `data/cftr2_2026-01-30.csv` (`build_cftr2.py`) — the benchmark itself, required
- `data/eve_cftr_2021-08.csv`, `esm1b_cftr.csv`, `revel_cftr_v1.3.csv`, `primateai_cftr.csv`,
  `spliceai_cftr_2021_v1.3.csv` — one `build_*.py` each
- `data/pangolin_cftr.csv` — `python build_pangolin.py` (runs the model; ~4 min on a GPU)
- `_tmp_fetch/alphamissense_cftr.tsv` — cache-backed
- network access for the live CADD v1.7 API (§10; it reaches none of the remaining gap, which
  is itself the finding)

A missing extract leaves that tool's column empty — the notebook still runs, but its coverage
numbers will be wrong. Check the `source` column, and rebuild before treating anything here as
a finding.

```bash
jupyter nbconvert --to notebook --execute --inplace predict/13_cftr2_benchmark.ipynb
```

## Caveats that bound every number here

- **CFTR2 is not independent of ClinVar** — they cross-cite, and REVEL is trained on
  ClinVar-lineage labels, so REVEL's agreement is partly circular (notebook 12).
- **42 negatives.** CFTR2 has ~1,245 CF-causing vs ~42 Non CF-causing variants, so specificity
  is noisy and accuracy is close to meaningless. Rank by sensitivity.
- **Coverage ≠ competence** — the splice tools call ~95% of CF-causing indels "no splice
  impact"; REVEL emits scores for 12 nonsense variants. That is why performance is reported
  twice, once restricted to each tool's actual domain.
- **The SpliceAI extract mixes masked SNVs with raw indels**, because Illumina's
  `masked.indel` release is usually a broken download. Masking is an all-or-nothing
  annotation filter (it zeroes gains at annotated splice sites and losses at unannotated
  positions), so `raw >= masked` always. Measured: it changes **nothing** in this benchmark —
  identical sensitivity/specificity, and 0 of 1,311 CFTR2 SNVs cross the 0.5 cut either way.
  Still check `score_type` before comparing scores across rows. Full working in notebook 09.
- **The mechanism split is a judgement call**, not a measurement: it assumes a splice model
  says nothing useful about a frameshift unless it actively flags disruption (≥ 0.5). The rule
  is in the notebook (§7) and one cell away from being changed.
