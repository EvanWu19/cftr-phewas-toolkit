# Data licensing & sharing — concerns before publishing data or results

**Status: partly resolved (July 2026). See the update box below.**
**I am not a lawyer and this is not legal advice.** This is an engineering-level
review of the data-use terms recorded in [`data_manifest.json`](data_manifest.json),
written so the open questions are visible instead of implicit. Every item marked
⚠️ should be confirmed against the provider's own current terms before anything is
published — the terms below are what the manifest records, not what I verified with
each licensor.

---

## ✅ Update — 2026-07-31: three blanks resolved, all in your favour

Three of the datasets this document treated as restricted or unconfirmed were
re-checked against the providers' current public terms. **All three moved toward
*more* permissive**, which materially widens what you can publish. See
[`PUBLISHING.md`](PUBLISHING.md) for the resulting publish/hold decision per file.

| Dataset | This doc originally said | Verified current term (2026-07-31) | Effect |
|---|---|---|---|
| **AlphaMissense** | CC BY-NC-SA (unconfirmed) | **CC BY 4.0** — DeepMind lifted the NonCommercial restriction on **13 Mar 2024** | ✅ **Publishable** with attribution. Also **dissolves License Conflict #3** — there is no ShareAlike term left to clash with PrimateAI's ND. |
| **gnomAD** | *(no license recorded)* | **ODbL** (Open Database License) + MIT terms | ✅ **Redistributable** with attribution, keep-open, and ShareAlike-under-ODbL. gnomAD explicitly states there are "no restrictions or embargoes" on publishing derived results. |
| **ClinVar** | *(no license recorded)* | **CC0 1.0 / U.S. public domain**; NCBI places "no restrictions on use or distribution" | ✅ **Redistributable**; attribution merely *requested*. CFTR2's expert-panel calls reach ClinVar under the same CC0. |

**Still restricted (unchanged, correctly recorded):** SpliceAI (CC BY-NC 4.0),
REVEL (non-commercial; redistribution unaddressed), PrimateAI via dbNSFP
(**CC BY-NC-ND** — the one hard blocker). **Still a genuine blank:** EVE — the
evemodel.org site is JavaScript-only and states no machine-readable license; its
terms must be confirmed on the site or by emailing the Marks lab before any EVE
extract is redistributed.

> One caveat on the gnomAD/ClinVar wins: gnomAD ships **SpliceAI annotations under
> CC BY-NC**, so if you ever republish a gnomAD *slice that includes its SpliceAI
> columns*, that column inherits the NC restriction even though the allele
> frequencies do not. Publish gnomAD frequency columns, not gnomAD's SpliceAI column.

---

## The short version

**The code, the notebooks, and the build scripts are yours to publish freely.**
The concern is only about **the data extracts** (`data/*.csv`) and **the derived
result tables** (`outputs/*.csv`) — because they contain, or are computed from,
third-party scores whose licenses restrict redistribution.

Right now the repo is **safe by construction**: `.gitignore` excludes `data/`
(except `data/README.md`), `outputs/`, `_tmp_fetch/`, and `archive/`, so none of it
has ever been pushed. **The concerns below are about changing that.**

**The single riskiest file** is
`outputs/predict_cftr2_benchmark_ALL.csv` — 2,097 variants × per-variant scores from
**seven** different tools (AlphaMissense, EVE, ESM1b, REVEL, PrimateAI, SpliceAI,
Pangolin) merged into one table. It is simultaneously subject to *every* restriction
listed below, and at least two of those restrictions appear to be mutually
incompatible (see [License conflict](#3-two-licenses-that-may-not-be-combinable)).

---

## Why this is not one question but three

Not all "data" carries the same risk. It helps to separate three layers:

| Layer | Example | Risk | Reasoning |
|---|---|---|---|
| **1. Raw third-party downloads** | `spliceai_scores.masked.snv.hg38.vcf.gz` (28.6 GB), `revel-v1.3_all_chromosomes.zip`, dbNSFP parquet | **Highest — do not publish** | Straight redistribution of the licensor's product. Never in the repo; keep it that way. |
| **2. Derived per-CFTR extracts** | `data/spliceai_cftr_2021_v1.3.csv` (~2.08 M rows, ~116 MB), `data/revel_cftr_v1.3.csv`, `data/primateai_cftr.csv` | **High — currently my main concern** | These are *filtered subsets*, not transformations. A CFTR-region slice of SpliceAI is still SpliceAI's scores; the substance is the licensor's. |
| **3. Derived aggregate results** | `outputs/predict_tool_performance.csv` (17 lines: n / TP / FP / sensitivity), `outputs/predict_category_summary.csv` (counts) | **Lowest — likely publishable** | Summary statistics *about* the data, not the data. Widely treated as ordinary research reporting. |

**The per-variant benchmark table sits in layer 2, not layer 3**, even though it
lives in `outputs/`. A table with one row per variant and a column per tool
reconstitutes a usable copy of each tool's scores for CFTR — that is redistribution,
regardless of the folder it is in.

---

## Per-dataset concerns

License text is as recorded in `data_manifest.json`. **"Confidence" is my confidence
in the *restriction*, not permission** — low confidence means *check before relying
on it*, never *assume it's fine*.

| Dataset | Recorded terms | My concern | Confidence |
|---|---|---|---|
| **PrimateAI** (via dbNSFP v5.0a) | dbNSFP **CC BY-NC-ND**; PrimateAI non-commercial | ⚠️ **The sharpest problem.** **ND = NoDerivatives.** A CC BY-NC-ND license permits sharing the work *verbatim* but prohibits distributing **adapted** material. Our `primateai_cftr.csv` is a filtered, re-keyed, column-renamed, de-duplicated extract — that is plausibly an adaptation. The merged benchmark table certainly is. | High that ND applies; unsure whether a filtered subset counts as a "derivative" in this context |
| **SpliceAI** (Illumina precomputed v1.3) | **CC BY-NC 4.0** — attribute SpliceAI + Illumina | ⚠️ **NonCommercial** is genuinely ambiguous for a public GitHub repo. Publishing to a free public repo is not obviously "commercial," but GitHub is a commercial platform and downstream users may be commercial. CC BY-NC also **requires** the license notice + attribution to travel with the data. Our extract is ~2.08 million rows / ~116 MB — substantial, not a token sample. | High that NC + attribution apply; unsure whether public-repo hosting violates NC |
| **REVEL** | "free for **non-commercial** use (contact authors otherwise)" | ⚠️ This is **not a standard license** — it is a usage permission with no stated *redistribution* right. Permission to *use* data is not permission to *republish* it. The instruction to "contact the authors" suggests case-by-case licensing. | High that redistribution is unaddressed |
| **EVE** | "EVE / evemodel.org terms — **CONFIRM before publishing this extract**" | ⚠️ **The terms were never confirmed.** The manifest has flagged this as unresolved since the extract was built. I cannot tell you whether publishing 26,809 EVE scores is permitted, because nobody has checked. | Unknown — this is a genuine blank |
| **AlphaMissense** | ~~CC BY-NC-SA (unconfirmed)~~ → **CC BY 4.0** (relicensed 13 Mar 2024) | ✅ **Resolved — publishable.** Attribution only; no NonCommercial, no ShareAlike. Removes the Conflict #3 blocker. Cite Cheng et al. 2023 + note the CC BY 4.0 release. | High — verified against DeepMind's release |
| **CFTR2** | "CFTR2 public data-use terms (cite CFTR2)" | Mildest of the restricted set — CFTR2 publishes this variant list openly and asks for citation. Still has *its own terms*; "public" is not the same as "unrestricted." Worth reading the actual data-use statement rather than assuming. | Medium-low concern, but read the terms |
| **gnomAD** | ~~*(no license recorded)*~~ → **ODbL** + MIT | ✅ **Resolved — redistributable.** Attribute gnomAD, keep derived DB open, share adaptations under ODbL. Do not reidentify participants. Publish frequency columns, **not** gnomAD's bundled CC BY-NC SpliceAI column. | High — verified against gnomAD policies page |
| **ClinVar** | ~~*(no license recorded)*~~ → **CC0 1.0 / public domain** | ✅ **Resolved — redistributable.** NCBI places no use/distribution restrictions; attribution requested. Don't imply endorsement. | High — verified against NCBI policy |
| **CADD** | *(no license recorded)*; live API | We store no bulk CADD data — scores are fetched live per variant. If you ever cache CADD responses into a published file, this needs checking (CADD has separate academic/commercial terms). | Low now; would rise if cached |
| **Pangolin** | non-commercial; cite Zeng & Li 2022 | `pangolin_cftr.csv` is **model output we computed ourselves**, which is a better position than redistributing someone's published table — but the model's license is still non-commercial, and output-ownership terms are worth confirming. | Medium — better footing, still NC |

---

## The three specific problems that worry me most

### 1. "Non-commercial" is undefined for a public GitHub repo

Five of these datasets are NC-restricted. Nobody is paying you, so intuitively this
feels non-commercial — but a public repo is a *distribution channel* to anyone,
including commercial users, hosted on a commercial platform. Different licensors read
this differently, and CC's own NC definition turns on whether use is "primarily
intended for or directed toward commercial advantage," which is about *your* purpose
but says little about downstream recipients. **I cannot tell you where this lands,
and I don't think you should guess.**

### 2. A "filtered subset" may still be redistribution

It is tempting to think that slicing SpliceAI down to CFTR (~2 M of billions of rows)
makes it a new thing. I don't think that reasoning is safe: the *value* in the file is
entirely the licensor's predictions, and the selection is mechanical (a coordinate
range). The strongest version of this concern is PrimateAI's **ND** term, which
targets exactly this.

### 3. Two licenses that may not be combinable

**➡️ 2026-07-31: this conflict is largely dissolved.** AlphaMissense is now
**CC BY 4.0** (no ShareAlike), so the direct SA-vs-ND contradiction is gone. What
remains is simpler but still real: `predict_cftr2_benchmark_ALL.csv` still merges
**PrimateAI (CC BY-NC-ND)**, **REVEL (NC)**, and **SpliceAI (CC BY-NC)** into one
table. ND means you may not distribute a derivative of the PrimateAI column *at all*,
so the merged table remains **hold / do-not-publish** as-is.

Original framing, for the record:

- ~~**AlphaMissense: CC BY-NC-SA** → derivatives must be shared alike~~ (no longer applies)
- **dbNSFP/PrimateAI: CC BY-NC-ND** → adapted material *must not* be distributed.

**Fix, unchanged:** split results per-source, drop the ND/NC columns, or publish only
aggregates. A benchmark table carrying **only your own derived calls** (`category`,
`mechanism`, `n_missense_tools`) plus the variant key is fully publishable and lets
anyone rebuild the scores with `build_*.py`.

---

## ~~Also: the repo itself has no LICENSE file~~ → ✅ Added 2026-07-31

A root [`LICENSE`](LICENSE) (MIT) now covers the code, notebooks, and docs, with an
explicit **scope notice** that it does **not** re-license any third-party data. A
[`CITATION.cff`](CITATION.cff) makes the toolkit citable. The original reasoning
below still explains *why* this matters:

- **For your code:** with no license, default copyright applies — others technically
  have no right to use, modify, or redistribute your notebooks, even though the repo
  is public. If you want this to be a teaching resource people can actually use, it
  needs an explicit license (MIT/Apache-2.0/CC-BY are common for this kind of work).
- **For the data question:** a repo-level license cannot grant rights you don't hold.
  If restricted data were ever committed, an MIT LICENSE at the root would be making
  a promise about third-party data you're not in a position to make.

---

## What I'd suggest (in order)

1. **Keep the current posture as the default.** Code + notebooks + build scripts +
   manifest + `data/README.md` ship; data and per-variant outputs stay gitignored.
   This is already implemented and it is the reason there's no problem today.
2. **Publish aggregates freely.** Tool performance stats, category counts, coverage
   figures, and the notebook *narrative* carry the scientific message without
   redistributing anyone's scores. `predict_tool_performance.csv` and
   `predict_category_summary.csv` are in this category.
3. **Add a LICENSE for your own work**, with an explicit note that it covers the code
   and documentation **only**, not any third-party data the scripts download.
4. **Resolve EVE and AlphaMissense**, which the manifest itself flags as unconfirmed.
   These are known blanks, not judgment calls — they can be closed by reading the
   current terms.
5. **If you want to publish per-variant tables, ask the licensors.** A short email to
   the REVEL authors, dbNSFP, and Illumina/SpliceAI describing exactly what you want
   to publish (a CFTR-only extract, for teaching, on a public repo) would replace all
   of the guessing above with a real answer. Several groups grant this readily.
6. **If you need per-variant results published before that resolves**, consider
   publishing **only the variant identifiers plus your own derived calls/categories**
   (e.g. `category`, `mechanism`, `n_missense_tools`) and *omitting the third-party
   score columns*. That preserves reproducibility — anyone can rebuild the scores with
   the `build_*.py` scripts — without you redistributing the scores.

---

## Questions worth putting to each provider

If you do write to them, these are the specific things left unanswered here:

- May I publish a **gene-restricted subset** of your scores in a **public,
  non-commercial, educational** repository?
- Does a **column-renamed / re-keyed / de-duplicated** extract count as a
  **derivative** under your terms? *(critical for dbNSFP's ND)*
- Is hosting on **GitHub** (a commercial platform, free tier) compatible with your
  **non-commercial** clause?
- May I publish a **merged multi-tool table**, or must each tool's scores be
  distributed separately under its own terms?

---

*Written by Claude (Anthropic) during a repo review, from the terms recorded in
`data_manifest.json`. Not legal advice — a checklist of what I could not resolve on
my own, so it can be resolved deliberately rather than by omission.*
