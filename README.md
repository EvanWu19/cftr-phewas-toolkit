# CFTR Variant Toolkit

A beginner-friendly, **provenance-honest** walkthrough of the computational tools
used to interpret CFTR variants — the missense pathogenicity predictors, the
splice predictors, the population-frequency reference, and the clinical/functional
truth sets. Built as a teaching companion to the project's A1 (missense triage)
and A2 (splice discordance) analyses.

Each tool gets its own Jupyter notebook explaining **what it is → what the score
means → the threshold and why → how to get the real data**. A final pair of
notebooks reproduces the combined worklists (12) and covers the methodology traps
— circular reasoning and training-data leakage (13).

---

> ## 📦 Data is NOT included in this repo
> This repository ships **code + notebooks + build scripts + a data manifest only**.
> Every dataset (raw sources *and* the derived per-CFTR score extracts) is
> license-restricted / non-commercial — SpliceAI CC BY-NC, REVEL non-commercial,
> **PrimateAI/dbNSFP CC BY-NC-ND**, AlphaMissense CC BY-NC-SA, EVE per-publication,
> CFTR2 data-use terms — and is **not redistributed here**. Regenerate each extract
> locally with the `build_*.py` scripts from your own downloads; `data_manifest.json`
> lists the exact source, version, and checksum. Committed notebooks keep their
> outputs so you can *read* the results; to *re-run*, fetch the data first.

## ⚠️ Read this first: REAL vs DEMO — and what a fresh clone ships

Two things matter most about this toolkit. First, which numbers come from *real*
data and which come from small *hand-curated demo tables*. Second — and just as
important — **this repo ships no data at all.** `data/`, `outputs/`, and
`_tmp_fetch/` are gitignored, so "REAL" below means *what a loader returns once you
have built or cached its extract locally*, **not** what you get on `git clone`.

| Source | REAL once you… | Coverage (once built) | Fresh clone gives you | Notebook |
|---|---|---|---|---|
| gnomAD v4 (allele freq) | cache `_tmp_fetch/` | ~2,466 missense / ~4,717 non-coding | ❌ `FileNotFoundError` | 01 |
| **AlphaMissense** | cache `_tmp_fetch/` | genome-wide CFTR missense | ❌ `FileNotFoundError` | 02 |
| ClinVar | cache `_tmp_fetch/` | genome-wide | ❌ `FileNotFoundError` | 07 |
| **CFTR2** (30 Jan 2026) | build `data/…csv` | ~2,097 variants / ~780 missense keys | ⚠️ DEMO fallback | 08 |
| **EVE** | build `data/…csv` | ~26,809 CFTR variants | ⚠️ DEMO fallback | 03 |
| **ESM1b** | build `data/…csv` | ~28,120 CFTR (saturation) | ⚠️ DEMO fallback | 04 |
| **REVEL** | build `data/…csv` | ~10,826 CFTR (coord-keyed; non-commercial) | ⚠️ DEMO fallback | 05 |
| **PrimateAI** | build `data/…csv` | ~1,976 CFTR (dbNSFP subset; non-commercial) | ⚠️ DEMO fallback | 06 |
| **SpliceAI** | build `data/…csv` | ~566k CFTR SNVs (Illumina v1.3, CC BY-NC) | ⚠️ DEMO fallback | 09 |
| CADD | — (live API) | per-variant | ✅ live (cache it) | 11 |
| Pangolin | — | 9 curated splice variants | 🟡 DEMO always | 10 |

There are **three** states, not two:

- **Cache-backed REAL** (gnomAD, AlphaMissense, ClinVar) — the loader **raises
  `FileNotFoundError`** until you populate `_tmp_fetch/` (each docstring says how).
  No demo fallback.
- **Build-backed REAL** (CFTR2, EVE, ESM1b, REVEL, PrimateAI, SpliceAI) — the loader
  returns REAL once you've built `data/<tool>.csv` with the matching `build_*.py`,
  otherwise it **falls back to a small DEMO table** (`source='DEMO'`) *with a
  warning*. Pass `strict=True` to raise instead of degrading silently.
- **Live / always-DEMO** — CADD is a live API (cache the responses for
  reproducibility); Pangolin ships as 9 curated DEMO variants only.

CFTR2 and ClinVar are real *truth sets* (databases), not predictors. Every DataFrame
a loader returns carries a `source` column (`REAL` / `DEMO`) so the two are never
confused. **Never quote a DEMO value as a finding.**

> ### 📍 Current status in *this* checkout
> This repo ships **code + notebooks + build scripts + manifest only** — no data.
> On a clean clone: the three cache-backed loaders raise `FileNotFoundError`; the six
> build-backed loaders return `source='DEMO'` (with a warning) until you run their
> `build_*.py`. See **[`data/README.md`](data/README.md)** for exactly how to fetch
> and build each extract, and **[`data_manifest.json`](data_manifest.json)** for the
> source URL, version, checksum, and license of every dataset.

---

## 🔰 Beginner primer (read before the numbers)

New to CFTR or variant prediction? These four notes unlock the rest of the README.

**Variant vocabulary.** A gene is read in three-letter codons that spell out a
protein. A **missense** variant changes one amino acid; a **synonymous** variant
changes the DNA but *not* the amino acid (yet can still disrupt splicing); a
**splice** variant hits the signals that cut introns out and join exons together; a
**deep intronic** variant sits far inside an intron but can create a *cryptic* splice
site; **non-coding** is the umbrella for everything that isn't a protein-coding
change (intronic, UTR, splice-region, synonymous). The A1 analysis is about
missense; A2 is about splice/non-coding — the variants the missense tools can't see.

**Two ways variants are named, and why joins break.** CFTR variants travel under
several keys, and most "biological disagreements" are really key mismatches:

| Key type | Example | Where it's used | Gotcha |
|---|---|---|---|
| Protein (1-letter) | `G551D` | AlphaMissense, EVE, ESM1b, gnomAD `protein_variant` | only exists for missense |
| Protein (3-letter) | `Gly551Asp` | some curated tables, ClinVar `Name` | must convert to 1-letter to join (`three_to_one()` in `toolkit.py`) |
| HGVS coding | `c.1652G>A` | CFTR2, clinical reports | needs the MANE transcript (`NM_000492.4`) |
| Legacy CFTR name | `2789+5G>A` | CFTR2 history, older literature | no formula — kept as a lookup column |
| Genomic coordinate | `7-117587799-G-A` (`chrom,pos,ref,alt`) | REVEL, SpliceAI, CADD, gnomAD `variant_id` | **CFTR is on the minus strand**, so the genomic `ref/alt` is the *complement* of the coding change (coding C>T shows as G>A) |

REVEL and SpliceAI are keyed by **coordinate** (no protein position), so they join by
`chrom,pos,ref,alt` — mind the strand. AlphaMissense/EVE/ESM1b join by
**protein_variant**. `hgvsp_to_short()` and `three_to_one()` in `toolkit.py` are the
helpers that normalise between them.

**Why "the predictor disagrees with ClinVar" isn't always evidence (circularity).**
Testing a predictor against labels it *trained on* is like grading a student with the
exact questions they studied — a high score proves memorisation, not understanding.
REVEL was trained on curated pathogenic/benign labels that share lineage with ClinVar,
so "REVEL vs ClinVar" is partly **circular**. AlphaMissense/EVE/ESM1b never saw
clinical labels, so comparing *them* to ClinVar is closer to independent evidence.
Notebook 13 quantifies this.

> ### ⚠️ A predictor score is not a clinical diagnosis
> Every threshold in this README (AlphaMissense ≥ 0.564, SpliceAI ≥ 0.5, …) is a
> deliberately simple single cut-point used to build **teaching worklists** — lists
> of variants worth a human's attention. They are **not** ACMG classifications and
> **not** diagnoses. Real clinical use applies *graded* thresholds and multiple lines
> of evidence (Pejaver 2022; notebooks 05/13). `score ≥ cutoff` ≠ "pathogenic".

---

## 📋 The one-page summary (A1 / A2 / B1)

This is the combined one-page summary these notebooks document — reproduced here
so a reader lands on the headline results and their source. Generated 2026-07-01
as part of a CFTR variant-interpretation collaboration.

> ### 🕰️ These are HISTORICAL / demo-reproduced numbers
> Everything in this "one-page summary" block (the Dashboard, the A1 discordance
> figure, and the Priority-1 table) reproduces the **original project webpage**,
> which was computed on a **~13-variant demo** footing before the real extracts
> existed. They are kept for provenance — *do not cite them as current results.*
> The **current real-data rerun** is the next section,
> **[“What the headline numbers actually mean”](#what-the-headline-numbers-in-the-summary-report-actually-mean)**,
> and the **[count glossary](#count-glossary--every-number-in-one-place)** maps each
> historical number to its corrected real value.

### Dashboard *(historical / demo — see banner above)*

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

> **Update (real EVE):** with EVE now real, **S912L scores benign (0.085)** and drops
> out, so the honest ≥3/5 count is **3, not 4** (archived integration notebook). The `4` above is the
> original webpage's demo-based figure, kept here as the reproduced summary.

> With the **real CFTR2** loader (notebook 08) now available, you can compute a fully
> real upgrade set: **256** variants that CFTR2 calls *"No interpretation available"*
> or *"Varying clinical consequence"* while AlphaMissense scores ≥ 0.564.

Source: the archived integration notebook wrote `outputs/A1_upgrade_worklist_REAL.csv`
(real AlphaMissense-vs-ClinVar upgrades). `outputs/` is gitignored and the notebook is
in `archive/` — this worklist will be regenerated by the future `predict/` pipeline.

### A2 — Splice-variant discordance

Deep-intronic, synonymous, and splice-site CFTR variants scored with SpliceAI +
Pangolin delta scores and CADD-Splice PHRED — **invisible to the A1 missense
tools** — then cross-checked against CFTR2/ClinVar. *Thresholds:* SpliceAI/Pangolin
DS_max ≥ 0.5 = HIGH, ≥ 0.2 = MODERATE; CADD-PHRED ≥ 15 = top 3%.

**Splice-risk VUS (primary worklist):**

| Variant | Type | SpliceAI | Pangolin | CADD | Tier |
|---|---|---|---|---|---|
| c.2657+120C>T † | deep intronic | 0.540 | 0.510 | 17.9 | HIGH |
| IVS8 5T (c.1210-34TG(12)T(5)) | deep intronic | 0.310 | 0.220 | 0.0 | MODERATE |

> **† `c.2657+120C>T` is a synthetic teaching example**, not a confirmed real
> observation, and all scores in this table are `source=DEMO` (hand-authored, not a
> real SpliceAI/Pangolin run). Do not treat this row as a real patient or database
> variant.

**Known CF splice variants (positive controls)** — any real SpliceAI/Pangolin run
should recover HIGH here: 2988+1G>A, 2789+5G>A, 2657+3A>G, 3849+10kb C>T,
3272-26A>G, 1811+1.6kb A>G.

Source: the archived integration notebook wrote `outputs/A2_splice_DEMO.csv` (all rows
`source=DEMO`). `outputs/` is gitignored and the notebook is in `archive/`.

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

> ### ✅ These are the CURRENT real-data rerun numbers
> Computed by notebook **12** (missense) and notebook **09** (splice) *with the real
> extracts built*. Where they differ from the historical block above, **these are the
> ones to cite.** (Reproduce them yourself only after building the extracts — see
> [`data/README.md`](data/README.md).)

The earlier one-page summary reported the historical `2496 / 413 / 403 / 10 / 4 / 1094`.
Here is what each one actually is, and its corrected real value:

- **2,496 → 2,466.** The real gnomAD v4 CFTR missense backbone is **2,466** variants;
  the historical **2,496** added ~30 hand-curated famous alleles (G551D, …). Of the
  2,466, **2,430** have an AlphaMissense score and **2,437** have ≥1 real predictor.
- **413 = 403 + 10 → 402 = 392 + 10.** This is a **two-source comparison —
  AlphaMissense vs ClinVar** (AM pathogenic while ClinVar is uncertain = *upgrade*, or
  AM benign while ClinVar is pathogenic = *downgrade*), **not** a five-tool vote. The
  live rerun over real data gives **392 upgrade + 10 downgrade = 402** (archived integration notebook);
  the historical `403/10/413` was the webpage's figure and drifts by a few with the
  ClinVar release used.
- **4 → 473.** "VUS but ≥3/5 tools pathogenic" was a demo-only figure (≥3/5 over ~13
  demo variants; it becomes **3** once real EVE drops S912L). With the real missense
  extracts built, the consensus runs over the **2,466 observed** variants and the real
  A1 Priority-1 worklist is **473** (archived integration notebook). (PrimateAI covers only ~53% of
  sites, so some variants are voted by 4 tools not 5; REVEL/PrimateAI carry
  circularity — notebook 12.)
- **1,094 → 4,260 / 164 HIGH.** "1,094 splice variants scored" originally meant **9
  demo variants** (the other ~1,085 unscored). With the real SpliceAI extract built
  (notebook 09), **4,260 of the 4,717** observed gnomAD non-coding CFTR variants get a
  real SpliceAI score, of which **164 are HIGH-impact** (≥0.5) and 82 MODERATE — the
  real A2 worklist. (Pangolin's `build_pangolin.py` now scores it for real too — see
  notebook 10; the archived integration notebook's A2 section still used the 9
  curated DEMO rows, so build the real join via notebook 09.)

### Count glossary — every number in one place

Every headline number that appears in this README, what it actually counts, and
whether it is current or historical. "Source" points at where it is computed.

| Number | What it counts | Status | Source |
|---|---|---|---|
| **2,466** | gnomAD v4 CFTR **missense** variants (no PASS/AC filter) — the real backbone | ✅ current | nb 01; `gnomad_missense.rows` |
| 2,133 | subset of those that are PASS + AC>0 in the gnomAD browser (stricter view) | ✅ current (alt filter) | manifest note |
| 2,430 / 2,437 | of the 2,466: have an AlphaMissense score / have ≥1 real predictor | ✅ current | archived integration nb |
| **2,496** | 2,466 + ~30 hand-curated famous alleles injected by the original script | 🕰️ historical | archived integration nb |
| **413** = 403 + 10 | AlphaMissense-vs-ClinVar discordance on the **original webpage** | 🕰️ historical | webpage |
| **402** = 392 + 10 | same comparison on the **current real rerun** (upgrade + downgrade) | ✅ current | archived integration nb |
| **473** | observed VUS with ≥3/5 tools pathogenic — the real A1 Priority-1 worklist | ✅ current | archived integration nb |
| 4 / 3 | historical Priority-1 (≥3/5 over ~13 demo variants; 3 after real EVE drops S912L) | 🕰️ historical/demo | webpage / archived integration nb |
| 256 | CFTR2 "no interpretation" or "varying consequence" **and** AM ≥ 0.564 (fully-real upgrade set) | ✅ current | nb 08/12 |
| **1,085** | older stated gnomAD non-coding count | 🕰️ stale | old table |
| **4,717** | gnomAD v4 CFTR **non-coding** variants (intron + synonymous + UTR + splice-region) | ✅ current | nb 01; `gnomad_noncoding.rows` |
| **1,094** | historical "splice variants scored" (really 9 DEMO scored + ~1,085 unscored) | 🕰️ historical/demo | webpage |
| **4,260 / 4,717** | non-coding variants that get a **real SpliceAI** score | ✅ current | nb 09 |
| **164 / 82** | of those 4,260: real SpliceAI HIGH (≥0.5) / MODERATE (0.2–0.5) | ✅ current | nb 09 |
| **566k** (566,106) | all precomputed SpliceAI CFTR SNVs in the built extract | ✅ current | nb 09; `spliceai.rows` |
| 9 | hand-curated DEMO splice variants (the A2 teaching table) | 🟡 DEMO | `toolkit.py` |

Coverage counts for the built extracts (saturation unless noted): EVE ~26,809 ·
ESM1b ~28,120 · REVEL ~10,826 raw (~10,127 after per-site dedup) · PrimateAI ~1,976
(dbNSFP ClinVar subset, ~53%, **not** saturation) · CFTR2 ~2,097.

See notebook **13** for why "predictor disagrees with ClinVar" is only evidence
when the predictor never trained on ClinVar-lineage labels (REVEL did; the
unsupervised tools did not).

---

## Notebooks

| # | File | Covers | Data on a fresh clone |
|---|---|---|---|
| 00 | `notebooks/00_overview_and_setup.ipynb` | setup + the provenance map | — |
| 01 | `notebooks/01_gnomad.ipynb` | gnomAD — population allele frequency as orthogonal evidence | REAL if cached, else error |
| 02 | `notebooks/02_alphamissense.ipynb` | AlphaMissense — genome-wide missense predictor | REAL if cached, else error |
| 03 | `notebooks/03_eve.ipynb` | EVE — unsupervised evolutionary model | REAL if built, else DEMO |
| 04 | `notebooks/04_esm1b.ipynb` | ESM1b — protein language model (backwards scale) | REAL if built, else DEMO |
| 05 | `notebooks/05_revel.ipynb` | REVEL — supervised ensemble + **circularity** | REAL if built, else DEMO |
| 06 | `notebooks/06_primateai.ipynb` | PrimateAI — semi-supervised (subset) | REAL if built, else DEMO |
| 07 | `notebooks/07_clinvar.ipynb` | ClinVar — crowd-sourced clinical truth set | REAL if cached, else error |
| 08 | `notebooks/08_cftr2.ipynb` | CFTR2 — disease-specific functional truth set | REAL if built, else DEMO |
| 09 | `notebooks/09_spliceai.ipynb` | SpliceAI — splice deltas (all CFTR SNVs) | REAL if built, else DEMO |
| 10 | `notebooks/10_pangolin.ipynb` | Pangolin — independent splice model (`build_pangolin.py` runs it for real) | DEMO scope |
| 11 | `notebooks/11_cadd.ipynb` | CADD — live deleteriousness score | REAL (live API) |
| 12 | `notebooks/12_decircularization_benchmark.ipynb` | **circularity & temporal-leakage reference** (which tool can you trust vs which truth set) | — |

One tool per notebook (03–11); notebook **12** is the circularity / temporal-leakage
reference. Recommended order: 01 → 12. If you read only one, read **12** (why a
predictor "disagreeing with ClinVar" is only sometimes evidence).

> **Roadmap / archived work:** the old A1/A2 *integration* notebook (which reproduced
> the discordance worklists) has moved to `archive/` (kept locally, **not** shipped).
> Real genome-wide prediction is being rebuilt in a dedicated **`predict/`** pipeline
> (next step). The current-real A1/A2 numbers below (402 / 473 / 256) were computed by
> that archived integration notebook and will be regenerated by `predict/`.

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
| **Pangolin** | splice | ≥ 0.5 | No | Zeng & Li 2022, PMID 35449021 |
| **CADD** | general | PHRED ≥ 15 | Proxy | Rentzsch 2021, PMID 33618777 |

⚠ REVEL is the one to distrust when benchmarking against ClinVar (it may have
trained on the same labels). Note ESM1b runs the *opposite* direction — more
negative = more damaging. All thresholds are single-cut simplifications; the ACMG
calibration (Pejaver 2022, PMID 36413997) uses *graded* thresholds — see notebooks
05/13.

---

## Setup

```bash
pip install -r requirements.txt
# then, from this folder:
jupyter lab           # or: jupyter notebook
```

Open `notebooks/00_overview_and_setup.ipynb` first.

### The REAL loaders need data you build/cache yourself

**Nothing under `data/`, `outputs/`, or `_tmp_fetch/` ships in the repo** — see
**[`data/README.md`](data/README.md)** for how to fetch and build every extract.
Concretely:

- **Cache-backed** (gnomAD, AlphaMissense, ClinVar) — read from `_tmp_fetch/`; each
  loader's docstring in `toolkit.py` gives the exact source and filter (gnomAD GraphQL
  API; tabix AlphaMissense to CFTR + filter UniProt P13569; filter ClinVar
  `variant_summary.txt.gz` to CFTR). Missing → `FileNotFoundError`.
- **Build-backed** (CFTR2, EVE, ESM1b, REVEL, PrimateAI, SpliceAI) — read from
  `data/<tool>.csv` built by the matching `build_*.py`. Missing → **DEMO fallback**
  (warning), or `strict=True` to raise.
- **No data needed** — the DEMO tables and the live CADD API work on a bare clone.

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
├── build_*.py             ← rebuild each data/ extract from its raw download
├── verify_data.py         ← checks locally-built extracts vs data_manifest.json
├── data_manifest.json     ← source/version/checksum/license for every dataset
├── build_pangolin.py      ← RUNS the Pangolin model locally (no genome download)
├── dev/_nbutil.py         ← author-side helper used to build the notebooks
├── notebooks/             ← 00–12 (see table above) — committed WITH outputs
├── data/                  ← gitignored; only data/README.md ships (build guide)
├── outputs/               ← gitignored; worklists (regenerated by predict/)
└── archive/               ← gitignored; superseded integration notebook
```

> Only the plain files and `notebooks/` are committed. `data/` (except its
> `README.md`), `outputs/`, and `_tmp_fetch/` are gitignored — a clone must rebuild
> them (see [`data/README.md`](data/README.md)).

## Known limitations (by design / honesty)

- **On a fresh clone, every predictor is DEMO or errors** — because no data ships
  (see the REAL/DEMO table). The six build-backed loaders (CFTR2, EVE, ESM1b, REVEL,
  PrimateAI, SpliceAI) become REAL once you run their `build_*.py`; gnomAD,
  AlphaMissense, and ClinVar need `_tmp_fetch/`. Pangolin has no precomputed release,
  but `build_pangolin.py` runs the real model locally (a curated run stays DEMO by
  scope). Build the extracts before treating any output as a finding, and check the
  `source` column.
- The 9 curated splice variants have hand-entered genomic coordinates; **only ~4
  of 9 validate against the GRCh38 reference** (notebook 11 shows the live check).
  One VUS (`c.2657+120C>T`) is an explicitly *synthetic* teaching example.
- **Reproducibility caveats:** ClinVar's release is **unpinned** (it updates ~weekly
  and drives the A1 VUS/discordance counts — record the exact release you use), and
  CADD is a **live API** (cache responses, or a rerun can change/fail on network
  behaviour rather than biology). Both are noted in `data_manifest.json`.
- The A1 discordance list did **not** apply a training-cutoff temporal hold-out;
  **notebook 12** (the circularity / temporal-leakage reference) explains how, and flags
  which tool×variant pairs risk leakage (e.g. F508del, reported 1989, vs any
  label-supervised tool trained after). Note CFTR2 is **not** an independent truth set —
  it shares clinical evidence with ClinVar and they cross-cite, so it is only *partially*
  orthogonal (its functional-assay axis), not circularity-free.

## Related work — we *reproduce*, not discover

Aggregating predictors + cross-checking ClinVar/CFTR2 is well established
(**OpenCRAVAT / OakVar**, **dbNSFP**, Ensembl VEP + plugins). The A1/A2 analyses
here **reproduce** published CFTR results rather than discovering them:
McDonald et al. 2023 (*PLOS ONE*, AlphaMissense's high false-positive rate vs
CFTR2), Tordai et al. 2024 (*Sci Data*), Bergougnoux et al. 2022 (*J Cyst Fibros*,
splice VUS), and the ACMG CFTR standard (Deignan et al. 2021, *Genet Med*). The
toolkit's contribution is the **honest REAL/DEMO provenance + CFTR teaching**.

## References

AlphaMissense (Cheng 2023, *Science*) · EVE (Frazer 2021, *Nature*) · ESM1b
(Brandes 2023, *Nat Genet*) · REVEL (Ioannidis 2016, *AJHG*) · PrimateAI
(Sundaram 2018, *Nat Genet*) · SpliceAI (Jaganathan 2019, *Cell*) · Pangolin
(Zeng & Li 2022, *Genome Biol*) · CADD-Splice (Rentzsch 2021, *Genome Med*) · REVEL
ACMG calibration (Pejaver 2022, *AJHG*) · gnomAD v4 · ClinVar · CFTR2 (cftr2.org).
