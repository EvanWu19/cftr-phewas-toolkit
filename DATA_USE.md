# Data sources, licenses & attribution

This toolkit's **code, notebooks, and documentation** are released under the
[MIT License](LICENSE). That license covers our own work only — it does **not**
re-license any third-party dataset or model output that the `build_*.py` / `fetch_*`
scripts download or compute. Each source below carries its own license and its own
attribution requirement. No third-party data is redistributed in this repository;
see [`data_manifest.json`](data_manifest.json) for the exact source, version, and
checksum of every dataset, and [`data/README.md`](data/README.md) for how to rebuild
each extract yourself.

*License terms were re-verified on 2026-07-31 against each provider's public terms.
This is a summary for reuse, not legal advice — confirm the current terms with the
provider before redistributing any extract.*

## How to read this

- 🟢 **Permissive** — redistributable with attribution (safe to share/host).
- 🟡 **Confirm** — terms unresolved; treat as do-not-redistribute until checked.
- 🔴 **Restricted** — do not redistribute the scores (non-commercial and/or
  no-derivatives, or redistribution not granted). Rebuild locally from the source.

| Dataset | License | Reuse | Cite |
|---|---|---|---|
| **AlphaMissense** | CC BY 4.0 (relicensed from CC BY-NC-SA on 2024-03-13) | 🟢 Permissive — attribution required | Cheng et al., *Science* 2023 (PMID 37733863); © Google DeepMind |
| **gnomAD** | ODbL + MIT terms | 🟢 Permissive — attribute, keep-open, ShareAlike-under-ODbL; no reidentification. Publish frequency columns only (gnomAD's bundled SpliceAI column is CC BY-NC) | gnomAD v4.1.1, Broad Institute |
| **ClinVar** | CC0 1.0 / U.S. public domain | 🟢 Permissive — attribution requested; don't imply endorsement | NCBI ClinVar (pin the release date) |
| **CFTR2** | Public data-use terms (cite CFTR2); expert-panel calls also in ClinVar under CC0 | 🟢 Permissive — cite CFTR2 | CFTR2 (cftr2.org), release 2026-01-30 |
| **ESM1b** | MIT (code); scores "per publication" | 🟡 Confirm redistribution of scores | Brandes et al. 2023 (PMID 37563329); ntranoslab/esm-variants |
| **EVE** | evemodel.org terms (no machine-readable license found) | 🟡 Confirm before redistributing | Frazer et al., *Nature* 2021 (PMID 34707284) |
| **Pangolin** | Non-commercial (our own model output) | 🔴 Non-commercial — do not commit to a public repo | Zeng & Li, *Genome Biol* 2022 (PMID 35449021) |
| **SpliceAI** | CC BY-NC 4.0 | 🔴 Non-commercial — attribute SpliceAI + Illumina; do not redistribute the precomputed extract | Jaganathan et al., *Cell* 2019 (PMID 30661751) |
| **REVEL** | Free for non-commercial use (contact authors otherwise) | 🔴 Restricted — grants use, not redistribution | Ioannidis et al. 2016 (PMID 27666373) |
| **PrimateAI** (via dbNSFP) | dbNSFP CC BY-NC-**ND**; PrimateAI non-commercial | 🔴 Restricted — NoDerivatives blocks redistributing a filtered extract | Sundaram et al. 2018 (PMID 30038395); dbNSFP v5.0a |
| **CADD** | Queried live per-variant; separate academic/commercial terms | 🔴 Do not cache-and-republish | Rentzsch et al. 2021 (PMID 33618777) |

## Attribution, if you reuse the 🟢 data

- **AlphaMissense** — © Google DeepMind, CC BY 4.0. Cheng et al., *Science* 2023.
- **gnomAD** — Genome Aggregation Database (Broad Institute), ODbL; v4.1.1.
- **ClinVar** — NCBI ClinVar, public domain (CC0); release date: *pin it*.
- **CFTR2** — Data from CFTR2 (cftr2.org); please cite CFTR2.

If you use this toolkit itself, please also cite it — see [`CITATION.cff`](CITATION.cff).
