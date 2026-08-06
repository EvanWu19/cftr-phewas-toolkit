# `data/` — how to fetch and build every extract

> **Nothing in this folder is committed.** `data/` and `outputs/` are gitignored
> (see the repo `.gitignore`), so a fresh clone contains **only this README** — no
> CSVs, no parquet, no raw archives. Every dataset below is license-restricted or
> non-commercial and is **not redistributed** here.

**There is one data folder now, not two.** Every dataset — whether it's pulled
live from an API, filtered from a bulk download, or produced by running a model
locally — lands in `data/`. There is no separate cache directory.

**The fetch/build code lives in the notebook that owns the tool, not in a
separate script.** Open that notebook, find the "Fetching/Building the REAL
data" cell near the top, and run it — it either queries a live source directly,
or tells you exactly what to manually download and where to put it before
re-running. This table is the summary; the notebook cell is the actual recipe.

Once you have rebuilt an extract, `python verify_data.py` checks its `sha256`/row
count against [`../data_manifest.json`](../data_manifest.json) (the machine-readable
version of everything below).

---

## Live-fetched (the notebook queries a public API/FTP directly — no manual download)

| Tool | Save-as in `data/` | Rows | Notebook | Source | Notes |
|---|---|---|---|---|---|
| **gnomAD missense** | `gnomad_cftr_missense.tsv` | 2,466 | `tools/01_gnomad.ipynb` | gnomAD v4.1.1 GraphQL API (`ENSG00000001626`, `gnomad_r4`) | **No PASS/AC filter** → 2,466 incl. AC0/VQSR-filtered; the browser's PASS+AC>0 view is 2,133 |
| **gnomAD non-coding** | `gnomad_cftr_noncoding.tsv` | 4,717 | `tools/01_gnomad.ipynb` (same cell) | gnomAD v4.1.1 GraphQL API | intron + synonymous + UTR + splice-region classes |
| **AlphaMissense** | `alphamissense_cftr.tsv` | 9,721 (→ 8,597 unique) | `tools/02_alphamissense.ipynb` | `AlphaMissense_hg38.tsv.gz` (genome-coordinate keyed), streamed from Google's public `dm_alphamissense` GCS bucket and filtered to UniProt P13569 | **CC BY 4.0** (DeepMind relicensed from CC BY-NC-SA on 2024-03-13) |
| **ClinVar** | `clinvar_cftr.tsv` + `clinvar_cftr.release.json` | ~6,100+ (drifts by design) | `benchmark/00_clinvar.ipynb` | `variant_summary.txt.gz` (default) or a pinned `archive/variant_summary_YYYY-MM.txt.gz`, streamed from `ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/` and filtered to `GeneSymbol=='CFTR', Assembly=='GRCh38'` | **Pinnable** — set `CLINVAR_RELEASE` in the fetch cell to `'latest'` (default; dated via HTTP `Last-Modified`, recorded in the `.release.json` sidecar) or `'YYYY-MM'` to reproduce a specific past month |
| **CADD** | *(no file — live per-call)* | — | `tools/09_cadd.ipynb` | `https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/` | Not reproducible unless you cache responses; a CADD version bump changes scores |

## Manual download + build (the notebook has a cell that reads the file once you provide it)

| Tool | Save-as in `data/` | Rows | Notebook | Manual download | Source | License |
|---|---|---|---|---|---|---|
| **CFTR2** | `cftr2_cftr.csv` + `cftr2_cftr.release.json` | 2,097 | `benchmark/01_cftr2.ipynb` | `CFTR2_30January2026.xlsx` (filename configurable in the build cell) | cftr2.org variant-list history tab | CFTR2 public data-use terms — **cite CFTR2** |
| **EVE** | `eve_cftr_2021-08.csv` | 26,809 | `tools/03_eve.ipynb` | `EVE_all_data.zip` → reads only `variant_files/CFTR_HUMAN.csv` from inside it | https://evemodel.org (release 2021-08, UniProt P13569) | EVE / evemodel.org terms — **CONFIRM before publishing** |
| **ESM1b** | `esm1b_cftr.csv` | 28,120 | `tools/04_esm1b.ipynb` | `ALL_hum_isoforms_ESM1b_LLR.zip` → reads only `…/P13569_LLR.csv` from inside it | HuggingFace Space `ntranoslab/esm_variants` | MIT (code); scores per publication |
| **REVEL** | `revel_cftr_v1.3.csv` | 10,826 (→ 10,127 unique coordinates) | `tools/05_revel.ipynb` | `revel-v1.3_all_chromosomes.zip` → streams `revel_with_transcript_ids` (6.5 GB member, stops after chr7) | https://sites.google.com/site/revelgenomics | **Non-commercial** (contact authors otherwise) |
| **PrimateAI** | `primateai_cftr.csv` | 1,976 | `tools/06_primateai.ipynb` | `dbNSFP5.0a_variant.…parquet` (reads only chromosome 7's columns) | dbNSFP v5.0a (Zenodo 15131632), `PrimateAI_score` column | dbNSFP **CC BY-NC-ND**; PrimateAI non-commercial |
| **SpliceAI** | `spliceai_cftr_2021_v1.3.csv` | ~2.08M (566,106 SNVs + 1,509,624 indels) | `tools/07_spliceai.ipynb` | `spliceai_scores.masked.snv.hg38.vcf.gz` + `.tbi` (~28.6 GB) and `spliceai_scores.raw.indel.hg38.vcf.gz` + `.tbi` (~69.3 GB) — the notebook seeks directly to the CFTR region via the `.tbi` index, never reading the full files | Illumina BaseSpace share https://basespace.illumina.com/s/otSPW8hnhaZR (`genome_scores_v1.3`) | **CC BY-NC 4.0** — attribute SpliceAI + Illumina |
| **Pangolin** | `pangolin_cftr.csv` | ~1,892 scored / 2,097 targets | `tools/08_pangolin.ipynb` | No data file — `pip install` the model package; the notebook auto-fetches+caches the ~215 kb CFTR reference region from Ensembl on first run (no whole-genome download) and needs `data/cftr2_cftr.csv` built first | github.com/tkzeng/Pangolin (Zeng & Li 2022, PMID 35449021) | non-commercial |

Notes:
- **EVE** carries an unconfirmed publishing license — treat the extract as local-only
  until you confirm evemodel.org's terms.
- **PrimateAI** comes from dbNSFP's ClinVar-re-annotated subset, so it is ~1,976
  *observed* CFTR variants, **not** a saturation set like EVE/ESM1b.
- **AlphaMissense**'s raw extract has 9,721 rows because the genome-coordinate file
  lists a protein change once per genomic codon path that produces it; the loader
  collapses that to 8,597 unique `protein_variant` values (highest score kept).
- **SpliceAI** is usually built **mixed masked/raw**: Illumina's `masked.indel`
  release is commonly a 0-byte failed download, so indels fall back to `raw.indel`
  while SNVs come from `masked.snv`. Every row carries `score_type`.
- **Pangolin**'s default scope (`SCOPE = "cftr2"` in the build cell) scores every
  CFTR2 variant with GRCh38 coordinates and labels the result `source='REAL'`;
  `SCOPE = "curated"` scores just 5 classic alleles and stays `source='DEMO'` —
  the label follows coverage, never the model.
- **CFTR2 has no historical archive** (checked directly against cftr2.org —
  there is no dated-release listing, unlike ClinVar). The build cell reads
  whatever release date is in the workbook's own header and records it in
  `cftr2_cftr.release.json`; reproducing a past run means manually sourcing
  that older workbook from cftr2.org yourself.

---

For the exact query strings, checksums, and build provenance of every item above,
see [`../data_manifest.json`](../data_manifest.json).
