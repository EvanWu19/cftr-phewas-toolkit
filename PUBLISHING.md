# PUBLISHING.md — how to publish as much as possible, safely

*Companion to [`DATA_LICENSING_CONCERNS.md`](DATA_LICENSING_CONCERNS.md) and
[`data_manifest.json`](data_manifest.json). Reviewed 2026-07-31. Not legal advice.*

This file answers one question: **given the goal of publishing as much CFTR
in-silico data as possible, what can go where?** It separates two very different
distribution channels, because the same dataset can be fine in one and forbidden in
the other:

- **Channel A — the public GitHub repo.** Redistribution *to the entire world*, with
  no ability to enforce who uses it or how (GitHub is a commercial platform). Only
  **permissive** licenses (CC BY, CC0, ODbL, public domain) belong here.
- **Channel B — a strictly non-commercial hosted search tool.** You control the terms
  of service, can require non-commercial use, and can display attribution/NC notices.
  This channel may *additionally* serve **NonCommercial** data — but **never**
  NoDerivatives (ND) data, and never data whose redistribution right is merely absent.

---

## The decision table

| Dataset | License (verified) | Channel A: commit to repo? | Channel B: host in NC tool? | Why |
|---|---|---|---|---|
| **AlphaMissense** | CC BY 4.0 | ✅ **Yes** | ✅ Yes | NC lifted Mar 2024; attribution only |
| **gnomAD** (freq cols) | ODbL + MIT | ✅ **Yes** | ✅ Yes | Attribute, keep-open, ODbL share-alike; no reidentification |
| **ClinVar** | CC0 / public domain | ✅ **Yes** | ✅ Yes | No restrictions; attribution requested |
| **CFTR2** | Public data-use (cite); also CC0 via ClinVar | ✅ **Yes** (cite) | ✅ Yes | Prefer the ClinVar-sourced calls (CC0) if you want zero ambiguity |
| **ESM1b** | MIT code; scores "per publication" | ⚠️ **Confirm first** | ⚠️ Confirm | Likely fine; no explicit redistribution license — verify with ntranoslab repo |
| **Pangolin** | Non-commercial (your own model output) | ❌ No | ✅ Yes (cite) | Your output, but model is NC → NC channel only |
| **SpliceAI** | CC BY-NC 4.0 | ❌ No | ✅ Yes (+ NC notice) **or** live-API | NC forbids the open repo; a genuine NC tool may host with attribution |
| **REVEL** | Non-commercial; redistribution unaddressed | ❌ No | ⚠️ Not without permission | Grants *use*, not *republish* — email authors or live-query |
| **PrimateAI / dbNSFP** | CC BY-NC-**ND** | ❌ No | ❌ **No** | ND blocks redistributing a derived extract in *any* channel |
| **EVE** | Unconfirmed (site is JS-only) | ❌ No | ❌ Not yet | Genuine blank — confirm evemodel.org terms before either |
| **CADD** | Live API; separate terms | ❌ No (don't cache-republish) | ✅ Live-query only | Keep querying live; don't host a cached copy |

**Your own derived columns** — `category`, `mechanism`, `n_missense_tools`, calls,
coverage flags — are **always publishable in both channels**, because they are your
work product, not a licensor's scores. This is the key to publishing the benchmark:
ship the variant key + your derived calls, omit the third-party score columns, and
anyone can rebuild the scores with `build_*.py`.

---

## What to actually publish, in order

### 1. Commit a `data/publishable/` set to the repo (permissive only)
The green rows above are all small (AlphaMissense ~9.7k, gnomAD ~2.5k+4.7k, ClinVar
~6.1k, CFTR2 ~2.1k rows). They are safe to commit. Recommended structure:

```
data/publishable/
  alphamissense_cftr.csv     # CC BY 4.0  — cite Cheng et al. 2023
  gnomad_cftr_missense.csv   # ODbL       — cite gnomAD v4.1.1 (freq cols only)
  gnomad_cftr_noncoding.csv  # ODbL
  clinvar_cftr.csv           # CC0        — pin the release date in the filename
  cftr2_2026-01-30.csv       # cite CFTR2
  LICENSES.md                # per-file license + attribution text (below)
```

Then relax `.gitignore` for exactly this folder (keep everything else ignored):

```gitignore
# publishable, permissively-licensed extracts (see PUBLISHING.md)
!data/publishable/
!data/publishable/**
```

Each file needs its attribution to travel *with* it. Put this in
`data/publishable/LICENSES.md`:

- **AlphaMissense** — © Google DeepMind, CC BY 4.0. Cheng et al., *Science* 2023
  (PMID 37733863). Predictions relicensed to CC BY 4.0 on 2024-03-13.
- **gnomAD** — Genome Aggregation Database, ODbL. gnomAD v4.1.1, Broad Institute.
  Frequency columns only; no attempt to reidentify individuals.
- **ClinVar** — NCBI ClinVar, public domain (CC0). Release date: `<pin it>`.
- **CFTR2** — Data from CFTR2 (cftr2.org); cite CFTR2. Release 2026-01-30.

### 2. Mint a citable release (Zenodo → DOI)
Committing CSVs bloats git history. Cleaner: attach the publishable set as a **GitHub
Release** asset, and enable the **Zenodo–GitHub integration** so each release gets a
**DOI**. That gives other bioinformaticians a stable, citable download and keeps the
repo lean. Put the DOI badge in the README.

### 3. Publish aggregates freely (already safe)
`outputs/predict_tool_performance.csv`, `predict_category_summary.csv`, coverage
figures, and the notebook narratives are summary statistics *about* the data, not the
data — publishable regardless of source license.

### 4. Publish the benchmark as calls-only
Ship `predict_cftr2_benchmark_ALL.csv` **without** the SpliceAI/REVEL/PrimateAI score
columns — keep the variant key + your `category`/`mechanism`/`n_missense_tools`. Fully
publishable, fully reproducible via `build_*.py`.

### 5. Close the two remaining blanks
- **ESM1b** — confirm the score-redistribution terms on the ntranoslab/esm-variants
  repo. If permissive, promote it to Channel A.
- **EVE** — confirm evemodel.org's terms (or email the Marks lab). Until then, hold.

---

## Permission emails (copy, fill the brackets, send)

Getting a one-line "yes" replaces all guesswork and would let you host REVEL, EVE, or
even a merged table. Several groups grant this readily for teaching resources.

### To the REVEL authors (redistribution of a CFTR subset)
> Subject: Permission to redistribute a CFTR-only REVEL subset (non-commercial, educational)
>
> Dear REVEL authors,
>
> I maintain an open, non-commercial teaching resource on CFTR variant
> interpretation (github.com/EvanWu19/cftr-phewas-toolkit). I would like to publish a
> **CFTR-gene-restricted subset** of REVEL v1.3 scores (~10,826 rows) so learners can
> reproduce the analysis without downloading the 6.5 GB genome-wide file.
>
> May I (1) redistribute this CFTR-only subset in a public repository, and (2) serve
> it through a strictly non-commercial variant-lookup web tool, in both cases with
> attribution and a link back to REVEL? Happy to add any notice you require.
>
> Thank you — [name, affiliation]

### To the EVE / Marks lab (confirm terms + redistribution)
> Subject: License terms for redistributing a CFTR EVE extract
>
> Dear EVE team,
>
> I'd like to include EVE scores for CFTR (UniProt P13569, ~26,809 variants, 2021-08
> release) in an open, non-commercial CFTR teaching toolkit and a non-commercial
> lookup tool. I couldn't find a machine-readable license on evemodel.org. Under what
> terms may I redistribute this gene-restricted extract, with attribution? — [name]

### To Illumina / SpliceAI (optional — clarify NC hosting)
> Subject: Hosting a CFTR SpliceAI subset in a non-commercial tool (CC BY-NC)
>
> Dear SpliceAI team,
>
> I understand the precomputed SpliceAI v1.3 scores are CC BY-NC 4.0. For a strictly
> non-commercial, educational CFTR lookup tool, may I serve a CFTR-region subset with
> attribution to SpliceAI + Illumina and the CC BY-NC notice displayed? — [name]

---

## The questions to put to each provider (if you write)
1. May I publish a **gene-restricted subset** of your scores in a **public,
   non-commercial, educational** repository?
2. Does a **column-renamed / re-keyed / de-duplicated** extract count as a
   **derivative** under your terms? *(critical for dbNSFP's ND)*
3. Is hosting on **GitHub** (commercial platform, free tier) compatible with your
   **non-commercial** clause?
4. May I publish a **merged multi-tool table**, or must each tool's scores be
   distributed separately under its own terms?
