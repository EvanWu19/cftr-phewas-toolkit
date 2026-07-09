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
| **CFTR2** (30 Jan 2026) | **REAL** | ~2,097 variants / 780 missense keys | 05 |
| CADD | **REAL** (live API) | per-variant | 06 |
| EVE, ESM1b | **DEMO** | ~13 curated variants | 03 |
| REVEL, PrimateAI | **DEMO** | ~13 curated variants | 04 |
| SpliceAI, Pangolin | **DEMO** | 9 curated splice variants | 06 |

**The only genome-wide real *predictor* here is AlphaMissense** — CFTR2 and
ClinVar are real *truth sets* (databases), not predictors. Every DataFrame a
loader returns carries a `source` column (`REAL` / `DEMO`) so the two are never
confused. Never quote a DEMO value as a finding.

> **CFTR2 is now REAL.** The full public CFTR2 variant list (30 January 2026
> release, ~2,097 variants) ships at `data/cftr2_2026-01-30.csv`, built from the
> cftr2.org download by `build_cftr2.py`. `tk.load_cftr2()` returns it with a
> 1-letter `protein_variant` key for the ~780 simple-missense variants (779 of
> which join to real AlphaMissense). Redistributed under cftr2.org's public
> data-use terms — please cite CFTR2 if you use it.

---

## 📋 The one-page summary (A1 / A2 / B1)

This is the combined one-page summary these notebooks document — reproduced here
so a reader lands on the headline results and their source. Generated 2026-07-01
as part of a CFTR variant-interpretation collaboration.

> **Read the numbers with the REAL/DEMO table above in mind.** The five-predictor
> framing below is the *presentation*; the section right after this one
> (“What the headline numbers actually mean”) and notebook 07 show that only
> AlphaMissense is genome-wide-real and the other predictors are DEMO on ~13
> variants. (CFTR2, one of the two truth sets, *is* now real — see above.)

### Dashboard

| Block | Metric | Value |
|---|---|---|
| **A1 · Missense** | CFTR missense variants scored | **2,496** |
| **A1 · Discordant** | Predictor↔database disagreements | **413** |
| **A1 · Priority 1** | VUS, ≥3/5 tools pathogenic | **4** |
| **A2 · Splice** | Splice variants scored | **1,094** |
| **A2 · High impact** | HIGH SpliceAI (+1 MODERATE) | **7** |
| **A2 · VUS worklist** | VUS with high splice risk | **2** |
| **B1 · Pipeline** | Nextflow processes succeeded | **6/6** |
| **B1 · Modules** | Containerised, CI-green | **5** |
| **Total worklist** | Variants for expert curation | **415** |

### A1 — Missense VUS triage / predictor discordance

Every CFTR missense variant scored by five orthogonal predictors (AlphaMissense,
EVE, ESM1b, REVEL, PrimateAI), then cross-checked against its CFTR2 class and
ClinVar assertion → a **413-variant discordance worklist** where computational
evidence disagrees with the curated classification (403 upgrade + 10 downgrade
candidates; 0 reverse discordance).

*Pathogenic cutoffs:* AlphaMissense ≥ 0.564 · EVE ≥ 0.50 · ESM1b ≤ −7.5 · REVEL ≥ 0.75 · PrimateAI ≥ 0.803.

**Priority 1 — VUS but ≥3/5 tools predict pathogenic** (primary upgrade candidates):

| Variant | HGVS c. | CFTR2 | ClinVar | AM | EVE | ESM1b | REVEL | PAI | Votes |
|---|---|---|---|---|---|---|---|---|---|
| **Tyr161Cys** | c.482A>G | VUS | Uncertain | 0.891 | 0.832 | −7.20 | 0.872 | 0.841 | **4/5** |
| **Gly970Asp** | c.2909G>A | VUS | Uncertain | 0.831 | 0.773 | −6.40 | 0.812 | 0.782 | **3/5** |
| **Ser912Leu** | c.2735C>T | VUS | Uncertain | 0.805 | 0.742 | −6.20 | 0.782 | 0.751 | **3/5** |
| **Val520Phe** | c.1558G>T | VUS | Uncertain | 0.778 | 0.718 | −6.00 | 0.755 | 0.731 | **3/5** |

> With the **real CFTR2** loader now available, notebook 05 also computes a fully
> real upgrade set: **256** variants that CFTR2 calls *"No interpretation available"*
> or *"Varying clinical consequence"* while AlphaMissense scores ≥ 0.564.

Source: `outputs/A1_upgrade_worklist_REAL.csv` (real AlphaMissense-vs-ClinVar upgrades) and the notebook 07 reconstruction.

### A2 — Splice-variant discordance

Deep-intronic, synonymous, and splice-site CFTR variants scored with SpliceAI +
Pangolin delta scores and CADD-Splice PHRED — **invisible to the A1 missense
tools** — then cross-checked against CFTR2/ClinVar. *Thresholds:* SpliceAI/Pangolin
DS_max ≥ 0.5 = HIGH, ≥ 0.2 = MODERATE; CADD-PHRED ≥ 15 = top 3%.

**Splice-risk VUS (primary worklist):**

| Variant | Type | SpliceAI | Pangolin | CADD | Tier |
|---|---|---|---|---|---|
| c.2657+120C>T | deep intronic | 0.540 | 0.510 | 17.9 | HIGH |
| IVS8 5T (c.1210-34TG(12)T(5)) | deep intronic | 0.310 | 0.220 | 0.0 | MODERATE |

**Known CF splice variants (positive controls)** — any real SpliceAI/Pangolin run
should recover HIGH here: 2988+1G>A, 2789+5G>A, 2657+3A>G, 3849+10kb C>T,
3272-26A>G, 1811+1.6kb A>G.

Source: `outputs/A2_splice_DEMO.csv` (all rows `source=DEMO`).

### B1 — cftr-varqc reproducible pipeline

An nf-core-style Nextflow workflow packaging the A1 + A2 logic into portable,
containerised, CI-tested infrastructure: takes a CFTR VCF, normalises to HGVS
(MANE NM_000492.4), joins precomputed effect scores, cross-checks CFTR2 + ClinVar,
and emits a MultiQC-style HTML report plus the discordance worklists. Last run:
**6/6 processes succeeded, 0 failed**; 5 local modules; CI stub-run green on every push.

```
INPUT VCF
   ├─ HGVS_NORMALIZE       bcftools norm + SnpEff HGVS (MANE NM_000492.4)
   ├─ SCORE_JOIN_MISSENSE  join AlphaMissense/EVE/ESM1b/REVEL/PrimateAI
   │      └─ CLINVAR_CHECK  → missense discordance worklist   (A1)
   ├─ SPLICE_SCORE         join SpliceAI/Pangolin + CADD-Splice PHRED
   │      └─ CLINVAR_CHECK  → splice discordance worklist      (A2)
   └─ MULTIQC_REPORT       aggregate summaries → HTML + TSV
```

*The `cftr-varqc` pipeline itself is not included in this repository — this
toolkit is the notebook companion. B1 is summarised here for context only.*

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
| 05 | `notebooks/05_clinvar_cftr2.ipynb` | ClinVar + CFTR2 — clinical & functional truth sets | **both REAL** |
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
UniProt P13569; filter ClinVar `variant_summary.txt.gz` to CFTR). The DEMO loaders,
the **real CFTR2 loader** (`data/cftr2_2026-01-30.csv` is committed), and the live
CADD API all work with **no cache**.

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
├── build_cftr2.py         ← rebuilds data/cftr2_*.csv from a CFTR2 release xlsx
├── _nbutil.py             ← tiny helper used to build the notebooks
├── data/                  ← REAL CFTR2 variant list (cftr2_2026-01-30.csv)
├── notebooks/             ← 00–08 (see table above)
└── outputs/               ← results written by notebook 07
```

## Known limitations (by design / honesty)

- Four of five missense predictors and both splice predictors are **DEMO** — swap
  in real data (per-notebook instructions) before treating any of it as findings.
  (gnomAD, AlphaMissense, ClinVar, **CFTR2**, and CADD are already REAL.)
- The 9 curated splice variants have hand-entered genomic coordinates; **only ~4
  of 9 validate against the GRCh38 reference** (notebook 06 shows the live check).
  One VUS (`c.2657+120C>T`) is an explicitly *synthetic* teaching example.
- The A1 discordance list did **not** apply a training-cutoff temporal hold-out;
  notebook 08 explains how to. Functional (CFTR2) truth is now available as a REAL
  loader — notebook 05 cross-checks it against ClinVar over 654 shared missense variants.

## References

AlphaMissense (Cheng 2023, *Science*) · EVE (Frazer 2021, *Nature*) · ESM1b
(Brandes 2023, *Nat Genet*) · REVEL (Ioannidis 2016, *AJHG*) · PrimateAI
(Sundaram 2018, *Nat Genet*) · SpliceAI (Jaganathan 2019, *Cell*) · Pangolin
(Zeng 2021, *Genome Biol*) · CADD-Splice (Rentzsch 2021, *Genome Med*) · REVEL
ACMG calibration (Pejaver 2022, *AJHG*) · gnomAD v4 · ClinVar · CFTR2 (cftr2.org).
