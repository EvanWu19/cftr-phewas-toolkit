# CFTR Variant Toolkit

A beginner-friendly, **provenance-honest** walkthrough of the computational tools
used to interpret CFTR variants — the missense pathogenicity predictors, the
splice predictors, the population-frequency reference, and the clinical/functional
truth sets. Built as a teaching companion to the project's A1 (missense triage)
and A2 (splice discordance) analyses.

Each tool gets its own Jupyter notebook explaining **what it is → what the score
means → the threshold and why → how to get the real data**. A final pair of
notebooks reproduces the combined worklists (07) and covers the methodology traps
— circular reasoning and training-data leakage (08).

---

## ⚠️ Read this first: REAL vs DEMO

The single most important thing about this toolkit — and about the analyses it
documents — is which numbers are *real* and which are small *hand-curated demo
tables* shipped so the notebooks run anywhere.

| Source | Status | Coverage | Notebook |
|---|---|---|---|
| gnomAD v4 (allele frequency) | **REAL** | ~2,466 missense / ~1,085 non-coding | 01 |
| **AlphaMissense** | **REAL** | genome-wide (all possible CFTR missense) | 02 |
| ClinVar | **REAL** | genome-wide | 05 |
| CADD | **REAL** (live API) | per-variant | 06 |
| EVE, ESM1b | **DEMO** | ~13 curated variants | 03 |
| REVEL, PrimateAI | **DEMO** | ~13 curated variants | 04 |
| CFTR2 class | **DEMO** | ~13 curated variants | 05 |
| SpliceAI, Pangolin | **DEMO** | 9 curated splice variants | 06 |

**The only genome-wide real predictor here is AlphaMissense.** Every DataFrame a
loader returns carries a `source` column (`REAL` / `DEMO`) so the two are never
confused. Never quote a DEMO value as a finding.

---

## What the headline numbers in the summary report actually mean

The earlier one-page summary reported `2496 / 413 / 403 / 10 / 4 / 1094`.
Notebook **07** reproduces each from the real data:

- **2,496** = ~2,466 real gnomAD missense variants + ~30 hand-curated variants.
- **413 = 403 upgrade + 10 downgrade**, and it is a **two-source comparison —
  AlphaMissense vs ClinVar** (AM pathogenic while ClinVar is uncertain, or AM
  benign while ClinVar is pathogenic). It is **not** a five-tool vote.
- **4** ("VUS but ≥3/5 tools pathogenic") lives **entirely in the ~13-variant
  DEMO set** — four of the five predictors are demo.
- **1,094** splice variants "scored" really means **9 scored** (with DEMO
  SpliceAI/Pangolin values); the other ~1,085 have no splice score at all.

See notebook **08** for why "predictor disagrees with ClinVar" is only evidence
when the predictor never trained on ClinVar-lineage labels (REVEL did; the
unsupervised tools did not).

---

## Notebooks

| # | File | Covers | Real data? |
|---|---|---|---|
| 00 | `notebooks/00_overview_and_setup.ipynb` | setup + the provenance map | — |
| 01 | `notebooks/01_gnomad.ipynb` | gnomAD — population allele frequency as orthogonal evidence | REAL |
| 02 | `notebooks/02_alphamissense.ipynb` | AlphaMissense — the one real genome-wide predictor | REAL |
| 03 | `notebooks/03_eve_esm1b.ipynb` | EVE (evolutionary model) + ESM1b (protein language model) | demo |
| 04 | `notebooks/04_revel_primateai.ipynb` | REVEL (supervised ensemble) + PrimateAI + **circularity** | demo |
| 05 | `notebooks/05_clinvar_cftr2.ipynb` | ClinVar + CFTR2 — clinical & functional truth sets | ClinVar REAL |
| 06 | `notebooks/06_splice_cadd.ipynb` | SpliceAI + Pangolin + CADD — splicing | CADD REAL |
| 07 | `notebooks/07_integration_A1_A2.ipynb` | **reproduce the A1/A2 worklists honestly** | mixed |
| 08 | `notebooks/08_decircularization_benchmark.ipynb` | **circular reasoning & training leakage** | — |

Recommended order: 01 → 08. If you read only two: **07** (what the numbers are)
and **08** (why to be careful).

---

## The tools at a glance

| Tool | Type | Score → pathogenic | Learns from clinical labels? | Paper |
|---|---|---|---|---|
| **AlphaMissense** | missense | ≥ 0.564 | No (sequence/structure) | Cheng 2023, PMID 37733863 |
| **EVE** | missense | ≥ 0.50 | No (MSA) | Frazer 2021, PMID 34707284 |
| **ESM1b** | missense | LLR ≤ −7.5 | No (protein LM) | Brandes 2023, PMID 37563329 |
| **PrimateAI** | missense | ≥ 0.803 | Proxy (primate/common) | Sundaram 2018, PMID 30038395 |
| **REVEL** | missense | ≥ 0.75 (graded) | **Yes (HGMD+ESP)** ⚠ | Ioannidis 2016, PMID 27666373 |
| **SpliceAI** | splice | DS_max ≥ 0.5 | No (GTEx junctions) | Jaganathan 2019, PMID 30661751 |
| **Pangolin** | splice | ≥ 0.5 | No | Zeng 2021, PMID 34753498 |
| **CADD** | general | PHRED ≥ 15 | Proxy | Rentzsch 2021, PMID 33618777 |

⚠ REVEL is the one to distrust when benchmarking against ClinVar (it may have
trained on the same labels). Note ESM1b runs the *opposite* direction — more
negative = more damaging. All thresholds are single-cut simplifications; the ACMG
calibration (Pejaver 2022, PMID 36413997) uses *graded* thresholds — see notebook
04/08.

---

## Setup

```bash
pip install -r requirements.txt
# then, from this folder:
jupyter lab           # or: jupyter notebook
```

Open `notebooks/00_overview_and_setup.ipynb` first.

### The REAL loaders need the download cache

The `load_*` REAL loaders read cached extracts from the project's `_tmp_fetch/`
directory (gnomAD, AlphaMissense-for-CFTR, ClinVar). Those large files are **not**
committed. To regenerate them, each loader's docstring in `toolkit.py` gives the
exact source and filter (gnomAD GraphQL API; tabix AlphaMissense to CFTR + filter
UniProt P13569; filter ClinVar `variant_summary.txt.gz` to CFTR). The DEMO loaders
and the live CADD API work with no cache.

### Run all notebooks headless

```bash
cd notebooks
for nb in *.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
```

---

## Files

```
cftr_variant_toolkit/
├── README.md              ← you are here
├── requirements.txt
├── toolkit.py             ← the library: loaders (REAL+DEMO), thresholds,
│                            tool registry, scoring helpers — all documented
├── _nbutil.py             ← tiny helper used to build the notebooks
├── notebooks/             ← 00–08 (see table above)
└── outputs/               ← results written by notebook 07
```

## Known limitations (by design / honesty)

- Four of five missense predictors and both splice predictors are **DEMO** — swap
  in real data (per-notebook instructions) before treating any of it as findings.
- The 9 curated splice variants have hand-entered genomic coordinates; **only ~4
  of 9 validate against the GRCh38 reference** (notebook 06 shows the live check).
  One VUS (`c.2657+120C>T`) is an explicitly *synthetic* teaching example.
- The A1 discordance list did **not** apply a training-cutoff temporal hold-out or
  validate against functional (CFTR2) truth — notebook 08 explains how to.

## References

AlphaMissense (Cheng 2023, *Science*) · EVE (Frazer 2021, *Nature*) · ESM1b
(Brandes 2023, *Nat Genet*) · REVEL (Ioannidis 2016, *AJHG*) · PrimateAI
(Sundaram 2018, *Nat Genet*) · SpliceAI (Jaganathan 2019, *Cell*) · Pangolin
(Zeng 2021, *Genome Biol*) · CADD-Splice (Rentzsch 2021, *Genome Med*) · REVEL
ACMG calibration (Pejaver 2022, *AJHG*) · gnomAD v4 · ClinVar · CFTR2 (cftr2.org).
