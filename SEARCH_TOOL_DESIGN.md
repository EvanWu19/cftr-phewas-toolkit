# SEARCH_TOOL_DESIGN.md — a CFTR variant in-silico lookup tool

*Design spec for a **strictly non-commercial** web tool that lets a bioinformatician
type a CFTR variant and see every in-silico prediction for it. Written 2026-07-31 to
be built against [`toolkit.py`](toolkit.py), [`data_manifest.json`](data_manifest.json),
and the publish/hold decisions in [`PUBLISHING.md`](PUBLISHING.md).*

---

## The core design principle

**Serving a score through a web UI is still redistribution.** So the licensing split
from `PUBLISHING.md` maps directly onto *how each score reaches the user*. Four
delivery modes, chosen per tool by its license — never one blanket approach:

| Mode | For which tools | Why it's clean |
|---|---|---|
| **A. Bundled** — stored in the tool's own DB | AlphaMissense (CC BY), gnomAD (ODbL), ClinVar (CC0), CFTR2, your derived calls | Permissive; ship with attribution |
| **B. Computed on demand** — run the model yourself | SpliceAI + Pangolin | *Your* model output (like `build_pangolin.py` already is), not Illumina's precomputed table — best footing |
| **C. Live-API proxy** — query provider at request time, cache briefly, never re-publish | CADD; optionally SpliceAI-lookup | Provider serves its own data; you never host it |
| **D. Link-out** — deep-link to the provider's browser | REVEL, PrimateAI/dbNSFP (ND), EVE (until confirmed) | You redistribute nothing; user gets the score from the source |

This is exactly how established browsers stay compliant: host what's permissive,
compute what you can, proxy or link the rest.

---

## Architecture

```
                    ┌──────────────────────────────────────────┐
   user types  ───► │  Frontend (variant box + results table)  │
   "G551D" or       │  shows a source + license badge per cell │
   "c.1652G>A" or   └───────────────────┬──────────────────────┘
   "7-117587799-G-A"                     │  /lookup?variant=...
                                         ▼
                    ┌──────────────────────────────────────────┐
                    │  API (FastAPI) — reuses toolkit.py         │
                    │  1. normalize key (three_to_one,           │
                    │     hgvsp_to_short, strand-aware coords)   │
                    │  2. fan out to the four modes ▼            │
                    └───┬──────────┬───────────┬────────────┬───┘
             (A) bundled│   (B) run│   (C) live│   (D) link  │
                        ▼          ▼           ▼             ▼
              ┌───────────────┐ ┌────────┐ ┌────────┐ ┌────────────────┐
              │ DuckDB/SQLite │ │SpliceAI│ │ CADD   │ │ deep-link URLs │
              │ AlphaMissense │ │+Pangoli│ │  API   │ │ REVEL / dbNSFP │
              │ gnomAD ClinVar│ │n model │ │        │ │ EVE browsers   │
              │ CFTR2 + calls │ │(docker)│ │(cache) │ │                │
              └───────────────┘ └────────┘ └────────┘ └────────────────┘
```

### Backend
- **FastAPI** (Python) so you reuse `toolkit.py` directly — the loaders, the key
  normalizers (`three_to_one()`, `hgvsp_to_short()`), and the strand handling are
  already written and tested. Don't reimplement them in JS.
- **Storage: DuckDB** (or SQLite) over the bundled permissive extracts. Two indexes:
  a genomic key `(chrom, pos, ref, alt)` for REVEL/SpliceAI/CADD/gnomAD-style joins
  and a `protein_variant` key for AlphaMissense/EVE/ESM1b. DuckDB reads Parquet
  directly and is fast enough to serve per-variant lookups with no server tuning.
- Every stored row keeps its `source` (`REAL`/`DEMO`) and a `license` string, so the
  API can return them and the UI can render provenance — never silently mix DEMO into
  a result.

### Splice scores as your own output (mode B)
The Broad's **SpliceAI-lookup** service (`github.com/broadinstitute/SpliceAI-lookup`)
runs *both* SpliceAI and Pangolin and ships **public Docker images**. Self-host it (or
run the models directly, as `build_pangolin.py` already does) and your splice numbers
are model outputs you generated — which sidesteps redistributing Illumina's CC BY-NC
precomputed file entirely. This is the single most useful move for splice coverage.
- Alternative for an MVP: proxy the public API at
  `https://spliceailookup-api.broadinstitute.org/spliceai/` — but it's rate-limited to
  a few requests/user/minute (interactive use only), so self-hosting is the real
  answer for a public tool.

### Live-API (mode C)
CADD stays live via `toolkit.fetch_cadd()`. Cache responses in-process to be polite,
but **do not** persist-and-republish a CADD cache (separate terms). Show a "fetched
live from CADD" note.

### Link-out (mode D)
For REVEL, PrimateAI/dbNSFP, and EVE (until its terms are confirmed), render a button
that deep-links to the provider's own variant page. The user gets the score from the
licensor; you host nothing. If you later get email permission (templates in
`PUBLISHING.md`), promote those tools from mode D to mode A.

### Frontend
- MVP: **Streamlit** or a single FastAPI + HTML page — fastest path to a working demo.
- Polished: a small **React/Next** front against the FastAPI JSON. A results table
  with one row per tool, columns for score / threshold-call / source / license badge,
  and a footer carrying all attributions + the non-commercial Terms of Use.

### Hosting (must stay non-commercial)
Good non-commercial fits: **Hugging Face Spaces** (free, clearly academic — great for
this), Google Cloud Run free tier, Render, or Fly.io. Whatever you pick, publish a
**Terms of Use** page stating the tool is for non-commercial research/education, which
is what lets mode-B splice output and any future NC-licensed data be served here but
not in the open repo.

---

## Licensing guardrails to bake into the UI (not optional)
1. **Per-cell provenance** — every score shows its source + license (e.g. a small
   "CC BY 4.0 · AlphaMissense" badge). This *is* the attribution the licenses require.
2. **Footer attributions** — AlphaMissense (Cheng 2023), gnomAD, ClinVar, CFTR2,
   SpliceAI+Illumina, Pangolin (Zeng & Li 2022), CADD — with links.
3. **Non-commercial Terms of Use** — the ToS is what distinguishes this channel from
   the public repo and makes NC hosting defensible.
4. **Never render DEMO as a finding** — carry the `source` column through to the UI
   and grey out / label DEMO rows, matching the toolkit's existing discipline.
5. **No PrimateAI/dbNSFP or unconfirmed-EVE data in the DB** — link-out only.

---

## Build phases

**Phase 0 — MVP (bundled + link-out).** FastAPI + `toolkit.py`, DuckDB over the four
permissive extracts, a variant box, a results table, mode-D link-outs for the
restricted tools. Ships a genuinely useful tool with zero licensing risk.

**Phase 1 — splice on demand.** Add the self-hosted SpliceAI-lookup container (SpliceAI
+ Pangolin) behind mode B. Now every variant, coding or splice, gets a splice call.

**Phase 2 — live CADD + caching + batch.** Wire mode-C CADD, add a short-TTL cache,
and a "paste a list of variants" batch endpoint.

**Phase 3 — promote on permission.** As REVEL/EVE emails come back "yes", move those
from link-out to bundled, with their required notices.

---

## Reuse checklist (already in this repo)
- Key normalization + strand handling → `toolkit.py` (`three_to_one`, `hgvsp_to_short`)
- Per-tool loaders + `source`/DEMO discipline → `toolkit.py` `load_*()`
- Provenance, versions, checksums, licenses → `data_manifest.json`
- Data integrity check → `verify_data.py`
- Splice model run → `build_pangolin.py` (pattern for self-hosting the models)
