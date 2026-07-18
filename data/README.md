# `data/` — how to fetch and build every extract

> **Nothing in this folder is committed.** `data/`, `outputs/`, and `_tmp_fetch/`
> are gitignored (see the repo `.gitignore`), so a fresh clone contains **only this
> README** — no CSVs, no parquet, no raw archives. This file is the recipe for
> rebuilding each dataset yourself. Every dataset below is license-restricted or
> non-commercial and is **not redistributed** here.

Once you have rebuilt an extract, `python verify_data.py` checks its `sha256`/row
count against [`../data_manifest.json`](../data_manifest.json) (the machine-readable
version of everything below).

There are two kinds of local data:

- **`data/` extracts** — small per-CFTR CSVs the `build_*.py` scripts distil from a
  large raw download. The `load_*()` loaders in `toolkit.py` read these; if the CSV
  is absent they fall back to a small DEMO table (or raise under `strict=True`).
- **`_tmp_fetch/` caches** — TSVs pulled live from an API/FTP by the `fetch_*`
  helpers. The gnomAD / AlphaMissense / ClinVar loaders **raise `FileNotFoundError`**
  until these exist (no DEMO fallback).

---

## `data/` extracts (build these with `build_*.py`)

| Tool | Save-as in `data/` | Rows | Build script | Raw download (fetch yourself) | Source | Retrieved | License |
|---|---|---|---|---|---|---|---|
| **CFTR2** | `cftr2_2026-01-30.csv` | 2,097 | `build_cftr2.py` | `CFTR2_30January2026.xlsx` | cftr2.org variant-list history tab | 2026-07-01 | CFTR2 public data-use terms — **cite CFTR2** |
| **EVE** | `eve_cftr_2021-08.csv` | 26,809 | `build_eve.py` | `EVE_all_data.zip` → `variant_files/CFTR_HUMAN.csv` | https://evemodel.org (release 2021-08, UniProt P13569) | 2026-07-12 | EVE / evemodel.org terms — **CONFIRM before publishing** |
| **ESM1b** | `esm1b_cftr.csv` | 28,120 | `build_esm1b.py` | `ALL_hum_isoforms_ESM1b_LLR.zip` → `…/P13569_LLR.csv` | HuggingFace Space `ntranoslab/esm_variants` | 2026-07-13 | MIT (code); scores per publication |
| **REVEL** | `revel_cftr_v1.3.csv` | 10,826 | `build_revel.py` | `revel-v1.3_all_chromosomes.zip` → `revel_with_transcript_ids` (6.5 GB) | https://sites.google.com/site/revelgenomics | 2026-07-13 | **Non-commercial** (contact authors otherwise) |
| **PrimateAI** | `primateai_cftr.csv` | 1,976 | `build_primateai.py` | `dbNSFP5.0a_variant.…parquet` | dbNSFP v5.0a (Zenodo 15131632), `PrimateAI_score` column | 2026-07-13 | dbNSFP **CC BY-NC-ND**; PrimateAI non-commercial |
| **SpliceAI** | `spliceai_cftr_2021_v1.3.csv` | 566,106 | `build_spliceai.py` | `spliceai_scores.masked.snv.hg38.vcf.gz` (~28.6 GB) | Illumina BaseSpace share https://basespace.illumina.com/s/otSPW8hnhaZR (`genome_scores_v1.3`) | 2026-07-12 | **CC BY-NC 4.0** — attribute SpliceAI + Illumina |

Notes:
- **EVE** carries an unconfirmed publishing license — treat the extract as local-only
  until you confirm evemodel.org's terms.
- **PrimateAI** comes from dbNSFP's ClinVar-re-annotated subset, so it is ~1,976
  *observed* CFTR variants, **not** a saturation set like EVE/ESM1b.
- **SpliceAI** ships as ~566k precomputed SNV scores (masked, Gencode v24, distance
  50 bp); the raw VCF stays external — never commit it.

## `_tmp_fetch/` caches (pull these with the `fetch_*` helpers)

| Dataset | Save-as in `_tmp_fetch/` | Rows | Fetcher | Source | Notes |
|---|---|---|---|---|---|
| **gnomAD missense** | `gnomad_cftr_missense.tsv` | 2,466 | `fetch_scores.py::fetch_gnomad_cftr` | gnomAD v4.1.1 GraphQL API (`ENSG00000001626`, `gnomad_r4`) | **No PASS/AC filter** → 2,466 incl. AC0/VQSR-filtered; the browser's PASS+AC>0 view is 2,133 |
| **gnomAD non-coding** | `gnomad_cftr_noncoding.tsv` | 4,717 | `fetch_splice_scores.py` | gnomAD v4.1.1 GraphQL API | intron + synonymous + UTR + splice-region classes |
| **AlphaMissense** | `alphamissense_cftr.tsv` | 9,721 | `fetch_scores.py::fetch_alphamissense` | `gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz` (bulk file), filtered to UniProt P13569 | **CONFIRM AlphaMissense data-use license** |
| **ClinVar** | `clinvar_cftr.tsv` | 6,101 | `fetch_scores.py::fetch_clinvar_cftr` | `ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`, filtered to CFTR | **Release date UNPINNED** — ClinVar updates ~weekly; record the exact release you used, since it drives the A1 VUS set |

## Live API (no local file)

- **CADD v1.7** — queried per-variant at
  `https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/` by `toolkit.fetch_cadd()`.
  Not reproducible unless you cache the responses; results can change or fail on
  network/API behaviour rather than biology.

---

For the exact query strings, checksums, and build provenance of every item above,
see [`../data_manifest.json`](../data_manifest.json).
