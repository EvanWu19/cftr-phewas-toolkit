"""
cftr_variant_toolkit / toolkit.py
=================================
The small set of constants and helpers genuinely SHARED across notebooks:
published thresholds, the tool/circularity registry, score->call logic, the
curated DEMO panel, and amino-acid key helpers — plus a thin load_<tool>()
reader per dataset (read data/<file>, tag source, DEMO fallback).

This file deliberately does NOT contain fetch/build logic anymore. Each
dataset's fetch-or-build code (how data/<file> gets created in the first
place — a live API call, a filtered download, or a locally-run model) lives
in the ONE notebook that owns that tool (tools/00-10, benchmark/00-01), right
next to the loader call and the license/provenance note for that source. That
used to be split across this file's docstrings, seven build_*.py scripts, and
two fetch_*.py scripts that were referenced everywhere but never committed —
opaque provenance masquerading as documentation. Read a tool's own notebook
to see exactly how its data was obtained; there is nowhere else to look.

Design goals
------------
1.  **Honest about provenance.** Every DataFrame this module returns carries a
    ``source`` column whose value is either ``"REAL"`` (downloaded from the
    primary source / queried live) or ``"DEMO"`` (a small, hand-curated table of
    illustrative values baked into this file). Never mix them silently.

2.  **One data/ folder, no hidden cache.** Every dataset — whether pulled live
    from an API, filtered from a bulk download, or produced by running a model
    locally — lands in ``data/`` once its owning notebook's fetch cell has run.
    There is no separate ``_tmp_fetch/`` cache; a loader either finds its file
    in ``data/`` or falls back to DEMO (or raises, for the three live-queried
    sources that have no DEMO fallback).

3.  **Beginner-readable.** Each function documents *what the tool is*, *what the
    score means*, *the threshold and why*, and *the primary reference*.

Real vs demo — and what a fresh clone actually ships
----------------------------------------------------
    IMPORTANT: this repo ships CODE + NOTEBOOKS + manifest ONLY. ``data/`` and
    ``outputs/`` are gitignored, so a fresh clone contains NONE of the datasets
    below. "REAL" describes what a loader returns *once its owning notebook's
    fetch cell has been run* — not what is in the box.

    On a fresh clone, each loader behaves as:
      REAL only after you run the fetch cell in its notebook (else FileNotFoundError):
        gnomAD v4 missense + non-coding (tools/01), AlphaMissense (tools/02),
        ClinVar (benchmark/00) — each queries/downloads live, no DEMO fallback
      REAL only after you run the build cell in its notebook, else DEMO fallback:
        CFTR2 (~2,097, benchmark/01), EVE (~26,809, tools/03),
        ESM1b (~28,120 saturation, tools/04), REVEL (~10,127 unique coordinates,
        non-commercial, tools/05), PrimateAI (~1,976, dbNSFP ClinVar subset,
        non-commercial, tools/06), SpliceAI (~566k SNVs; CC BY-NC 4.0, tools/07)
      REAL, queried live per-call (no local data needed):
        CADD v1.7 REST API (tools/09)
      DEMO always (hand-curated illustrative values — NOT real predictions):
        Pangolin (9 curated splice variants only; tools/08 also has a REAL
        model-run path over the full CFTR2 list)

    => The six build-locally loaders fall back to a tiny DEMO table when their
       extract is missing. Pass ``strict=True`` to raise instead of silently
       degrading; the default emits a warning. See ``data/README.md`` for the
       provenance summary and ``data_manifest.json`` for exact hashes/licenses.

References
----------
    AlphaMissense : Cheng et al. 2023   Science    PMID 37733863
    EVE           : Frazer et al. 2021  Nature     PMID 34707284
    ESM1b         : Brandes et al. 2023 Nat Genet  PMID 37563329
    REVEL         : Ioannidis et al. 2016 AJHG     PMID 27666373
    PrimateAI     : Sundaram et al. 2018 Nat Genet PMID 30038395
    SpliceAI      : Jaganathan et al. 2019 Cell    PMID 30661751
    Pangolin      : Zeng & Li 2022  Genome Biol PMID 35449021 (github.com/tkzeng/Pangolin)
    CADD-Splice   : Rentzsch et al. 2021 Genome Med  PMID 33618777
    REVEL thresholds : Pejaver et al. 2022 AJHG (ACMG calibration) PMID 36413997
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Paths — ONE data/ folder for every dataset (no _tmp_fetch/ cache anymore)
# ─────────────────────────────────────────────────────────────────────────────
PKG_DIR  = Path(__file__).resolve().parent
DATA_DIR = PKG_DIR / "data"           # gitignored; built by each notebook's fetch/build cell
OUT_DIR  = PKG_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

CFTR2_CSV     = DATA_DIR / "cftr2_2026-01-30.csv"           # benchmark/01
EVE_CSV       = DATA_DIR / "eve_cftr_2021-08.csv"            # tools/03
SPLICEAI_CSV  = DATA_DIR / "spliceai_cftr_2021_v1.3.csv"     # tools/07
ESM1B_CSV     = DATA_DIR / "esm1b_cftr.csv"                  # tools/04
REVEL_CSV     = DATA_DIR / "revel_cftr_v1.3.csv"              # tools/05
PRIMATEAI_CSV = DATA_DIR / "primateai_cftr.csv"               # tools/06
PANGOLIN_CSV  = DATA_DIR / "pangolin_cftr.csv"                # tools/08

# CFTR locus, GRCh38. CFTR is on the PLUS (forward) strand of chromosome 7 (7q31.2),
# so cDNA/coding alleles read the same as the genomic ref/alt — no complementing.
CFTR_CHR, CFTR_START, CFTR_END = "7", 117_470_098, 117_667_108
CFTR_UNIPROT = "P13569"            # AlphaMissense protein key
CFTR_MANE_TX = "NM_000492.4"       # MANE Select transcript


# ─────────────────────────────────────────────────────────────────────────────
# Published decision thresholds (one binary cut per tool)
# ─────────────────────────────────────────────────────────────────────────────
# NB: these are deliberately *simple* single cut-points taken from each tool's
# calibration paper. Real ACMG use (esp. REVEL) applies GRADED thresholds for
# different evidence strengths — see tools/05 and Pejaver 2022.
THRESHOLDS = {
    # higher score = more pathogenic
    "am":         {"path": 0.564, "benign": 0.340},   # AlphaMissense class cuts
    "eve":        {"path": 0.500, "benign": 0.500},   # EVE posterior midpoint
    "revel":      {"path": 0.750, "benign": 0.290},   # REVEL (ACMG moderate / BP4)
    "primate_ai": {"path": 0.803, "benign": 0.483},   # PrimateAI
    # LOWER (more negative) = more damaging
    "esm1b":      {"path": -7.5,  "benign": -7.5},    # ESM1b LLR cut
    # splice deltas: higher = more impact
    "spliceai":   {"high": 0.5,   "moderate": 0.2},   # SpliceAI DS_max
    "pangolin":   {"high": 0.5,   "moderate": 0.2},   # Pangolin
    "cadd":       {"path": 15.0},                     # CADD PHRED (top ~3%)
}

# Which direction is "more damaging", per tool. Used by call_from_score().
HIGHER_IS_WORSE = {"am": True, "eve": True, "revel": True, "primate_ai": True,
                   "esm1b": False}

# ─────────────────────────────────────────────────────────────────────────────
# Tool registry — metadata for the de-circularization / benchmark notebook (13)
# ─────────────────────────────────────────────────────────────────────────────
# `learning` is the key field for circularity reasoning:
#   "unsupervised"  — trained only on sequences/MSAs, NEVER on clinical labels
#                     → safe(ish) to benchmark against ClinVar
#   "semi"          — trained on a proxy for benignity (e.g. common variants)
#   "supervised"    — trained on curated pathogenic/benign labels (HGMD/ClinVar
#                     lineage) → benchmarking against ClinVar is PARTLY CIRCULAR
TOOL_REGISTRY = {
    "AlphaMissense": dict(kind="missense", learning="unsupervised",
        signal="protein language model + structure (AlphaFold), weakly calibrated on population frequency",
        circularity="low", pmid="37733863"),
    "EVE": dict(kind="missense", learning="unsupervised",
        signal="deep generative model over the multiple-sequence alignment (evolutionary constraint)",
        circularity="low", pmid="34707284"),
    "ESM1b": dict(kind="missense", learning="unsupervised",
        signal="protein language model log-likelihood ratio (variant vs wild-type)",
        circularity="low", pmid="37563329"),
    "PrimateAI": dict(kind="missense", learning="semi-supervised",
        signal="deep net trained on common human/primate missense as a benign proxy",
        circularity="medium", pmid="30038395"),
    "REVEL": dict(kind="missense", learning="supervised",
        signal="random-forest ENSEMBLE of 13 scores, trained on curated pathogenic/benign labels",
        circularity="HIGH (label lineage overlaps ClinVar/HGMD)", pmid="27666373"),
    "SpliceAI": dict(kind="splice", learning="supervised (on GTEx splice junctions, not clinical labels)",
        signal="deep net predicting acceptor/donor gain/loss probability from sequence",
        circularity="low-for-clinvar", pmid="30661751"),
    "Pangolin": dict(kind="splice", learning="supervised (on splice usage across species/tissues)",
        signal="deep net predicting splice-site usage change",
        circularity="low-for-clinvar", pmid="35449021"),
    "CADD": dict(kind="general", learning="semi-supervised",
        signal="SVM/logistic on many annotations; trained on simulated-vs-observed variants",
        circularity="medium", pmid="33618777"),
}

# Publication / training-freeze year per tool — the anchor for the temporal-leakage
# reference (tools/10). A variant's clinical label can only leak into a tool if the
# tool's training data postdates the variant's first pathogenic report AND the tool
# learned from clinical labels. `label_supervised` marks the tools where a post-report
# training year is a *direct* leakage risk (REVEL); unsupervised/proxy tools carry only
# INDIRECT risk (benchmarks/frequency calibration), so a date flag there is weaker.
# NOTE: CADD is NOT trained on clinical labels — it contrasts observed vs simulated
# variants (proxy). CFTR2 is NOT independent of ClinVar (they cross-cite); see nb12.
TOOL_YEAR = {
    "AlphaMissense": 2023, "EVE": 2021, "ESM1b": 2023, "PrimateAI": 2018,
    "REVEL": 2016, "SpliceAI": 2019, "Pangolin": 2022, "CADD": 2021,  # CADD v1.7
}
LABEL_SUPERVISED = {  # trained directly on curated clinical pathogenic/benign labels?
    "AlphaMissense": False, "EVE": False, "ESM1b": False, "PrimateAI": False,
    "REVEL": True, "SpliceAI": False, "Pangolin": False, "CADD": False,
}


# =============================================================================
# 1. gnomAD — population allele frequencies (REAL, cached)
# =============================================================================
def load_gnomad_missense() -> pd.DataFrame:
    """CFTR missense variants + allele frequency from gnomAD v4 (REAL).

    gnomAD (Genome Aggregation Database) is a reference set of human variation
    from ~800k individuals. Allele frequency (AF) is the single most useful
    *orthogonal* filter for pathogenicity: a truly CF-causing recessive allele
    should be rare (AF typically < 1e-3), so a "pathogenic" prediction on a
    common variant is a red flag.

    Returns columns: variant_id, hgvs_c, hgvs_p, protein_variant, consequence,
    gnomad_af, source. One row per missense variant (~2,466).

    Thin reader only: the gnomAD GraphQL fetch that builds
    data/gnomad_cftr_missense.tsv lives in tools/01_gnomad.ipynb (no precomputed
    per-gene download exists, so the notebook queries the live API directly).
    """
    fp = DATA_DIR / "gnomad_cftr_missense.tsv"
    if not fp.exists():
        raise FileNotFoundError(
            f"{fp} missing. Run the fetch cell in tools/01_gnomad.ipynb "
            "(queries the gnomAD v4 GraphQL API live; no download needed).")
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    df = df.rename(columns={"hgvsc": "hgvs_c", "hgvsp": "hgvs_p"})
    df["protein_variant"] = df["hgvs_p"].apply(hgvsp_to_short)
    df["source"] = "REAL"
    return df[["variant_id", "hgvs_c", "hgvs_p", "protein_variant",
               "consequence", "gnomad_af", "source"]]


def load_gnomad_noncoding() -> pd.DataFrame:
    """CFTR non-coding / synonymous / splice-region variants from gnomAD v4 (REAL).

    Same source as load_gnomad_missense(), but the consequence classes that the
    missense predictors CANNOT see: intronic, synonymous, UTR, splice-region.
    These are the substrate for the A2 splice analysis.

    Thin reader only: fetched by the same cell as load_gnomad_missense(), in
    tools/01_gnomad.ipynb -> data/gnomad_cftr_noncoding.tsv (~4,717 variants).
    """
    fp = DATA_DIR / "gnomad_cftr_noncoding.tsv"
    if not fp.exists():
        raise FileNotFoundError(
            f"{fp} missing. Run the fetch cell in tools/01_gnomad.ipynb.")
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    df["source"] = "REAL"
    return df


# =============================================================================
# 2. AlphaMissense — the one REAL genome-wide missense predictor here
# =============================================================================
def load_alphamissense() -> pd.DataFrame:
    """AlphaMissense pathogenicity for every CFTR missense change (REAL).

    AlphaMissense (Cheng 2023, DeepMind) adapts AlphaFold into a variant-effect
    predictor. It is *unsupervised* w.r.t. clinical labels — trained on protein
    sequences/structures plus weak population-frequency calibration, NOT on
    ClinVar pathogenic/benign labels. That is why it is a good tool to compare
    *against* ClinVar without circular reasoning (see tools/10).

    Score `am_pathogenicity` in [0,1]; AlphaMissense's own 3-class cut-points:
        >= 0.564  -> "likely_pathogenic"
        <= 0.340  -> "likely_benign"
        else      -> "ambiguous"

    Returns: protein_variant, am_score, am_class, source (~8,597 rows).
    Thin reader only: the fetch (streams DeepMind's genome-wide
    AlphaMissense_hg38.tsv.gz and filters to UniProt P13569) lives in
    tools/02_alphamissense.ipynb -> data/alphamissense_cftr.tsv.
    """
    fp = DATA_DIR / "alphamissense_cftr.tsv"
    if not fp.exists():
        raise FileNotFoundError(
            f"{fp} missing. Run the fetch cell in tools/02_alphamissense.ipynb.")
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    # collapse multiple codon changes giving the same protein change → 1 row
    df = (df.sort_values("am_pathogenicity", ascending=False)
            .drop_duplicates("protein_variant", keep="first"))
    df = df.rename(columns={"am_pathogenicity": "am_score"})
    df["source"] = "REAL"
    return df[["protein_variant", "am_score", "am_class", "source"]]


# =============================================================================
# 3. ClinVar — clinical assertions (REAL, cached)
# =============================================================================
def load_clinvar() -> pd.DataFrame:
    """ClinVar clinical significance for CFTR variants (REAL).

    ClinVar aggregates clinical assertions (Pathogenic / Benign / Uncertain /
    Conflicting) submitted by labs. It is the de-facto clinical "truth" set —
    but treat it with care: assertions vary in review status (star rating) and,
    crucially, some predictors were TRAINED on ClinVar-lineage labels, so
    comparing those predictors to ClinVar is partly circular (see tools/10).

    Returns: protein_variant, clinvar_sig, review_status, clinvar_call, source.
    `clinvar_call` collapses the free-text significance to pathogenic/benign/
    uncertain via cv_class().
    Thin reader only: the fetch (streams NCBI's variant_summary.txt.gz and
    filters to GeneSymbol=='CFTR', Assembly=='GRCh38') lives in
    benchmark/00_clinvar.ipynb -> data/clinvar_cftr.tsv. ClinVar updates
    ~weekly and is otherwise unpinned — record the exact retrieval date.
    """
    fp = DATA_DIR / "clinvar_cftr.tsv"
    if not fp.exists():
        raise FileNotFoundError(f"{fp} missing. Run the fetch cell in benchmark/00_clinvar.ipynb.")
    df = pd.read_csv(fp, sep="\t", low_memory=False)
    df["protein_variant"] = df["Name"].apply(extract_hgvsp_from_name)
    df = df.dropna(subset=["protein_variant"]).drop_duplicates("protein_variant")
    out = pd.DataFrame({
        "protein_variant": df["protein_variant"],
        "clinvar_sig":     df["ClinicalSignificance"],
        "review_status":   df.get("ReviewStatus"),
    })
    out["clinvar_call"] = out["clinvar_sig"].apply(cv_class)
    out["source"] = "REAL"
    return out.reset_index(drop=True)


# =============================================================================
# 4. Missense predictors — EVE / ESM1b / REVEL / PrimateAI
#    REAL when the data/*.csv extract exists; DEMO fallback otherwise.
# =============================================================================
# Each loader returns REAL scores once you have run its owning notebook's build
# cell to produce `data/<tool>.csv` locally (see data/README.md). None of those CSVs ship in
# the repo (data/ is gitignored), so on a fresh clone every loader falls back to
# the small curated DEMO table below. The returned `source` column ("REAL"/"DEMO")
# always tells you which you got; pass strict=True to raise instead of falling back.
#
# demo columns: protein_variant, eve_score, esm1b_score, revel_score, primate_ai_score
_DEMO_MISSENSE = [
    # protein_variant, am, eve, esm1b, revel, primate_ai, cftr2_class, clinvar_sig, gnomad_af
    ("G551D", 0.999, 0.98, -12.1, 0.94, 0.93, "CF-causing",  "Pathogenic",             1.3e-4),
    ("F508del", None, None, None,  None, None, "CF-causing",  "Pathogenic",             None ),  # not a missense
    ("R117H", 0.34, 0.42, -4.1,  0.52, 0.44, "CF-causing (mild)", "Pathogenic",         2.5e-3),
    ("R334W", 0.92, 0.90, -9.2,  0.88, 0.85, "CF-causing",  "Pathogenic",               1.1e-4),
    ("G85E",  0.95, 0.93, -10.1, 0.90, 0.88, "CF-causing",  "Pathogenic",               8.0e-5),
    ("D1152H",0.695,0.61, -5.4,  0.672,0.60, "VUS",         "Conflicting interpretations", 4.2e-5),
    ("R668C", 0.20, 0.18, -3.1,  0.22, 0.28, "Not CF-causing","Benign",                 3.0e-3),
    ("Tyr161Cys", 0.891, 0.832, -7.2, 0.872, 0.841, "VUS", "Uncertain significance",    6.4e-6),
    ("Gly970Asp", 0.831, 0.773, -6.4, 0.812, 0.782, "VUS", "Uncertain significance",    4.1e-6),
    ("Ser912Leu", 0.805, 0.742, -6.2, 0.782, 0.751, "VUS", "Uncertain significance",    2.3e-6),
    ("Val520Phe", 0.778, 0.718, -6.0, 0.755, 0.731, "VUS", "Uncertain significance",    5.4e-6),
    ("His949Tyr", 0.762, 0.741, -5.6, 0.741, 0.700, "VUS", "Uncertain significance",    1.1e-5),
    ("Pro205Ser", 0.741, 0.700, -5.1, 0.719, 0.660, "VUS", "Uncertain significance",    8.7e-6),
    ("M470V", 0.09, 0.12, -1.8, 0.10, 0.30, "Not CF-causing", "Benign",                 0.45  ),
]
_DEMO_COLS = ["protein_variant", "am_score", "eve_score", "esm1b_score",
              "revel_score", "primate_ai_score", "cftr2_class", "clinvar_sig", "gnomad_af"]


def _demo_frame() -> pd.DataFrame:
    df = pd.DataFrame(_DEMO_MISSENSE, columns=_DEMO_COLS)
    df["source"] = "DEMO"
    return df


def _missing_extract(name: str, path: Path, strict: bool) -> None:
    """Signal that a REAL extract CSV is absent, before falling back to DEMO.

    Under ``strict=True`` this raises FileNotFoundError; otherwise it warns so the
    silent demo fallback becomes visible. Keeps the teaching convenience (the repo
    ships no data, so a fresh clone still runs) without letting a broken/missing
    extract masquerade as a successful REAL load.
    """
    msg = (f"{name}: real extract not found at {path} — returning the DEMO "
           f"fallback table (source='DEMO'). Run the build cell in that tool's "
           f"own notebook (see data/README.md), or pass strict=True to raise.")
    if strict:
        raise FileNotFoundError(msg)
    warnings.warn(msg, stacklevel=3)


_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
_THREE_RE = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def three_to_one(protein_variant: str) -> str:
    """Normalise a missense key to the 1-letter form used by the REAL tables.

    'Tyr161Cys' -> 'Y161C'. Already-1-letter keys (e.g. 'G551D') and anything that
    isn't a simple 3-letter missense are returned unchanged. Use this to join the
    curated DEMO variants (some keyed 3-letter) onto REAL EVE/AlphaMissense/ClinVar,
    which use 1-letter keys — the "join hygiene" issue notebooks 12/13 discuss.
    """
    if not isinstance(protein_variant, str):
        return protein_variant
    m = _THREE_RE.match(protein_variant.strip())
    if not m:
        return protein_variant
    wt, pos, mt = m.group(1), m.group(2), m.group(3)
    if wt in _AA3TO1 and mt in _AA3TO1:
        return f"{_AA3TO1[wt]}{pos}{_AA3TO1[mt]}"
    return protein_variant


def load_eve(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """EVE score per CFTR missense variant.

    EVE (Evolutionary model of Variant Effect, Frazer 2021) is an UNSUPERVISED
    deep generative model trained on the multiple-sequence alignment of the
    protein family — it scores how well a variant fits evolutionary constraint.
    Score in [0,1]; >= 0.5 ~ pathogenic. No clinical labels used → low
    circularity vs ClinVar.

    REAL if the extract exists (`data/eve_cftr_2021-08.csv`, ~26,809 scored
    variants), built by the manual-download build cell in tools/03_eve.ipynb
    from the EVE release (evemodel.org, CFTR = UniProt P13569); keyed by the
    1-letter protein_variant. Columns:
    protein_variant, eve_score, eve_class, source. The extract is gitignored, so
    on a fresh clone this falls back to the tiny curated DEMO table (source='DEMO')
    with a warning — pass strict=True to raise instead, or demo=True to request the
    DEMO table silently.
    """
    if not demo and EVE_CSV.exists():
        df = pd.read_csv(EVE_CSV, dtype={"protein_variant": "string"})
        df["source"] = "REAL"
        return df[["protein_variant", "eve_score", "eve_class", "source"]]
    if not demo:
        _missing_extract("EVE", EVE_CSV, strict)
    d = _demo_frame()
    return d[["protein_variant", "eve_score", "source"]].dropna(subset=["eve_score"])


def load_esm1b(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """ESM1b LLR per CFTR missense variant — REAL if the extract exists.

    ESM1b (Brandes 2023) is a protein LANGUAGE model. It scores a variant by the
    log-likelihood ratio (LLR) of the mutant vs wild-type amino acid — a more
    NEGATIVE LLR means the model finds the mutation more surprising/damaging.
    Cut: LLR <= -7.5 ~ pathogenic. Unsupervised → low circularity.

    REAL if the extract exists: full CFTR **saturation** LLR (~28,120 variants, all
    1,480 residues) from `data/esm1b_cftr.csv`, built by the manual-download build
    cell in tools/04_esm1b.ipynb from the ntranoslab esm_variants release
    (canonical UniProt **P13569**). protein_variant
    keyed. The extract is gitignored → fresh clone falls back to the DEMO table
    (source='DEMO') with a warning; strict=True raises, demo=True is silent.
    """
    if not demo and ESM1B_CSV.exists():
        df = pd.read_csv(ESM1B_CSV, dtype={"protein_variant": "string"})
        df["source"] = "REAL"
        return df[["protein_variant", "esm1b_score", "source"]]
    if not demo:
        _missing_extract("ESM1b", ESM1B_CSV, strict)
    d = _demo_frame()
    return d[["protein_variant", "esm1b_score", "source"]].dropna(subset=["esm1b_score"])


def load_revel(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """REVEL score per CFTR missense variant — REAL if the extract exists.

    REVEL (Ioannidis 2016) is a SUPERVISED random-forest ENSEMBLE of 13 other
    predictors, trained on curated pathogenic/benign variants. Score in [0,1];
    the ACMG calibration (Pejaver 2022) gives GRADED cut-points, of which 0.75
    is a common single 'likely pathogenic' point.

    ⚠ CIRCULARITY: REVEL's training labels share lineage with ClinVar/HGMD, so
    'REVEL disagrees with ClinVar' can partly reflect label leakage, not
    independent evidence. Handled in tools/10.

    REAL if the extract exists: genome-wide REVEL v1.3 for CFTR (~10,127 variants)
    from `data/revel_cftr_v1.3.csv`, built by the manual-download build cell in
    tools/05_revel.ipynb. **Keyed by genomic
    coordinate** (chrom,pos,ref,alt) — the REVEL table has no protein position, so
    join it onto observed variants by coordinate, not protein_variant. (CFTR is
    plus-strand, so the coding alleles match the genomic ref/alt directly.) Non-commercial license. The extract is gitignored → fresh
    clone falls back to the DEMO table (source='DEMO') with a warning; strict=True
    raises, demo=True is silent.
    """
    if not demo and REVEL_CSV.exists():
        df = pd.read_csv(REVEL_CSV)
        df["chrom"] = df["chrom"].astype(str)
        # REVEL can list >1 transcript row per genomic site; keep the max per site
        df = (df.sort_values("revel_score", ascending=False)
                .drop_duplicates(["chrom", "pos", "ref", "alt"]).reset_index(drop=True))
        df["source"] = "REAL"
        return df
    if not demo:
        _missing_extract("REVEL", REVEL_CSV, strict)
    d = _demo_frame()
    return d[["protein_variant", "revel_score", "source"]].dropna(subset=["revel_score"])


def load_primateai(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """PrimateAI score per CFTR missense variant — REAL if the extract exists.

    PrimateAI (Sundaram 2018) is a deep net trained SEMI-supervised on common
    human & non-human primate missense variants as a proxy for benignity.
    Score in [0,1]; >= 0.803 ~ pathogenic. Medium circularity.

    REAL if the extract exists: PrimateAI for CFTR from `data/primateai_cftr.csv`,
    built by the manual-download build cell in tools/06_primateai.ipynb from the
    **dbNSFP v5.0a** parquet. ⚠ COVERAGE:
    dbNSFP's ClinVar-re-annotated subset, so ~1,976 observed CFTR variants (NOT
    saturation). protein_variant + coordinate keyed. Non-commercial. The extract is
    gitignored → fresh clone falls back to the DEMO table (source='DEMO') with a
    warning; strict=True raises, demo=True is silent.
    """
    if not demo and PRIMATEAI_CSV.exists():
        df = pd.read_csv(PRIMATEAI_CSV, dtype={"protein_variant": "string"})
        df["source"] = "REAL"
        return df
    if not demo:
        _missing_extract("PrimateAI", PRIMATEAI_CSV, strict)
    d = _demo_frame()
    return d[["protein_variant", "primate_ai_score", "source"]].dropna(subset=["primate_ai_score"])


def load_cftr2(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """CFTR2 clinical-functional class per CFTR variant — REAL if the extract exists.

    CFTR2 (cftr2.org) is the clinical-functional reference for CF: it labels
    variants as 'CF-causing', 'Varying clinical consequence', 'Non CF-causing',
    or 'No interpretation available', based on patient data + in-vitro CFTR
    function. Its *in-vitro functional-assay* axis is evidence ClinVar largely
    lacks, so CFTR2 is **partially** more orthogonal than ClinVar for benchmarking.
    ⚠ It is NOT an independent gold standard: CFTR2 and ClinVar share clinical/
    patient evidence, ClinVar entries cite CFTR2, and CFTR2 informs the ACMG CFTR
    guidance that ClinVar submitters follow — so benchmarking against it is not
    circularity-free (tools/10).

    REAL if the extract exists: the full public CFTR2 variant list (30 January 2026
    release, ~2,097 variants) built from a manually-downloaded cftr2.org workbook
    by the build cell in benchmark/01_cftr2.ipynb into
    ``data/cftr2_2026-01-30.csv``. Returns a 1-letter
    ``protein_variant`` key (e.g. 'G551D') for the ~780 simple-missense variants
    so it joins onto the AlphaMissense/gnomAD tables; non-missense rows carry an
    empty key but keep their legacy/cDNA names and genomic coordinates.

    The extract is gitignored (CFTR2's data-use terms allow local use; rebuild it
    yourself from cftr2.org — see data/README.md), so a fresh clone falls back to
    the tiny embedded curated set (source='DEMO') with a warning; strict=True
    raises, demo=True is silent. Please cite CFTR2 (cftr2.org) if you use it.
    """
    if not demo and CFTR2_CSV.exists():
        df = pd.read_csv(CFTR2_CSV, dtype={"protein_variant": "string"})
        df["protein_variant"] = df["protein_variant"].fillna("")
        df["source"] = "REAL"
        return df
    if not demo:
        _missing_extract("CFTR2", CFTR2_CSV, strict)
    d = _demo_frame()
    return d[["protein_variant", "cftr2_class", "source"]]


def load_cftr2_demo() -> pd.DataFrame:
    """Tiny embedded CFTR2 class table (DEMO). Prefer load_cftr2() for the full
    REAL list; kept for the early teaching cells that predate the real loader."""
    d = _demo_frame()
    return d[["protein_variant", "cftr2_class", "source"]]


# =============================================================================
# 5. SpliceAI / Pangolin — DEMO curated splice variants
# =============================================================================
# 9 hand-curated CFTR splice variants. The SpliceAI/Pangolin numbers below are
# ILLUSTRATIVE (source='DEMO'), NOT the output of a real SpliceAI/Pangolin run.
# The 6 CF-causing ones are genuine, well-known splice alleles; the 3 VUS
# (incl. the deep-intronic c.2657+120C>T) are teaching examples — do not treat
# c.2657+120C>T as a confirmed real observation.
#
# tuple: (variant_id, hgvs_c, legacy, type, cftr2_class, clinvar_sig,
#         chrom, pos, ref, alt, DS_AG, DS_AL, DS_DG, DS_DL, pangolin, cadd, note)
KNOWN_CF_SPLICE_VARIANTS = [
    ("7-117645994-C-T","c.3849+10246C>T","3849+10kb C>T","deep_intronic","CF-causing","Pathogenic",
     "7",117_645_994,"C","T", 0.81,0.03,0.02,0.01, 0.79, 24.2,
     "Creates cryptic exon 6b; 10.2 kb downstream of exon 22. Common in UK/N-European."),
    ("7-117587799-G-A","c.2657+5G>A","2789+5G>A","splice_site","CF-causing","Pathogenic",
     "7",117_587_799,"G","A", 0.02,0.01,0.91,0.05, 0.88, 27.8,
     "IVS14b+5G>A; weakens intron-14b donor. 5th most common CF variant in some populations."),
    ("7-117616814-A-G","c.3140-26A>G","3272-26A>G","deep_intronic","CF-causing","Pathogenic",
     "7",117_616_814,"A","G", 0.76,0.02,0.03,0.01, 0.73, 22.4,
     "New acceptor 26 bp upstream of exon 20. Common in Hispanic/Latino."),
    ("7-117587800-T-C","c.2657+3A>G","2657+3A>G","splice_site","CF-causing","Pathogenic",
     "7",117_587_800,"T","C", 0.02,0.01,0.88,0.04, 0.85, 26.1,
     "Donor splice-site variant intron 14b. Less common than 2657+5."),
    ("7-117548628-G-A","c.1210-34TG(12)T(5)","IVS8_5T","deep_intronic","VUS","Uncertain significance",
     "7",117_548_628,"G","A", 0.28,0.12,0.31,0.08, 0.22, 14.1,
     "Poly-T tract (5T/7T/9T) affects exon 9 skipping. 5T alone insufficient; 5T+TG12 raises risk."),
    ("7-117592260-C-T","c.2988+1G>A","2988+1G>A","splice_site","CF-causing","Pathogenic",
     "7",117_592_260,"C","T", 0.02,0.01,0.97,0.03, 0.96, 32.1,
     "Canonical +1 donor site, intron 16. Abolishes splicing."),
    ("7-117559590-G-A","c.1680-886A>G","1811+1.6kbA>G","deep_intronic","CF-causing","Pathogenic",
     "7",117_559_590,"G","A", 0.68,0.04,0.02,0.01, 0.65, 20.8,
     "New acceptor in intron 11, 1.6 kb from canonical donor."),
    ("7-117548735-A-G","c.1210A>G","syn context","synonymous","VUS","Uncertain significance",
     "7",117_548_735,"A","G", 0.18,0.05,0.14,0.06, 0.12, 11.3,
     "Coding synonymous near exon-intron boundary; possible ESR disruption."),
    ("7-117588032-C-T","c.2657+120C>T","c.2657+120C>T","deep_intronic","VUS","Uncertain significance",
     "7",117_588_032,"C","T", 0.54,0.03,0.02,0.01, 0.51, 17.9,
     "TEACHING EXAMPLE (synthetic): illustrative deep-intronic VUS. Not a confirmed real observation."),
]
_SPLICE_COLS = ["variant_id","hgvs_c","legacy_name","variant_type","cftr2_class","clinvar_sig",
                "chrom","pos","ref","alt","DS_AG","DS_AL","DS_DG","DS_DL",
                "pangolin_score","cadd_phred","note"]


def load_splice_demo() -> pd.DataFrame:
    """Curated CFTR splice variants with DEMO SpliceAI/Pangolin/CADD scores.

    SpliceAI (Jaganathan 2019) predicts, per position, the probability that a
    variant creates/destroys an acceptor or donor site: four deltas DS_AG
    (acceptor gain), DS_AL (acceptor loss), DS_DG (donor gain), DS_DL
    (donor loss), each 0–1. DS_max = max of the four. DS_max >= 0.5 = high
    impact, >= 0.2 = moderate. Pangolin (Zeng & Li 2022) is a similar model giving
    one 0–1 score.

    ⚠ The DS_/pangolin numbers here are DEMO (hand-authored). For real scores:
    SpliceAI — precomputed VCF from Illumina BaseSpace (login) or the Broad
    SpliceAI-lookup app; Pangolin — run the model locally (GPU). See README.
    """
    df = pd.DataFrame(KNOWN_CF_SPLICE_VARIANTS, columns=_SPLICE_COLS)
    df["spliceai_ds_max"] = df[["DS_AG", "DS_AL", "DS_DG", "DS_DL"]].max(axis=1)
    df["source"] = "DEMO"
    return df


def load_spliceai(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """SpliceAI delta scores for CFTR — REAL if the extract exists.

    REAL if the extract exists: the precomputed Illumina **SpliceAI v1.3** scores for
    the whole CFTR region (`data/spliceai_cftr_2021_v1.3.csv`, ~2.08 M records), built
    by the manual-download build cell in tools/07_spliceai.ipynb (a hand-rolled .tbi
    seek, since pysam doesn't build on Windows). Keyed by genomic coordinate (chrom,pos,ref,alt); columns
    DS_AG/DS_AL/DS_DG/DS_DL and spliceai_ds_max (>= 0.5 high, >= 0.2 moderate). Join
    onto observed variants (e.g. gnomAD non-coding) by coordinate to build the real A2
    splice worklist.

    **SNVs *and* indels** (~566 k + ~1.51 M). The indels matter: no missense predictor
    can score one, so without them every CFTR2 deletion/insertion was invisible to every
    tool — including F508del, which SpliceAI scores DS_max = 0.01 (correct: it is a
    folding defect, not a splice one). See predict/13.

    ⚠ The extract is usually **MIXED masked/raw**: Illumina's `masked.indel` release is
    very often a 0-byte failed download, so the builder falls back to `raw.indel` while
    SNVs come from `masked.snv`. Each row carries `score_type` ('masked'/'raw') and
    `variant_class` ('snv'/'indel'). Raw does not zero out biologically implausible
    directions, so raw >= masked for the same variant: across the CFTR SNVs the two
    agree exactly 95.8% of the time (mean |diff| 0.0019), but 113 cross the 0.5 cut.
    Don't compare a masked score against a raw one and call the difference biology.

    The extract is gitignored (CC BY-NC 4.0) → fresh clone falls back to the 9
    curated variants (load_splice_demo) with a warning; strict=True raises,
    demo=True is silent. Note those hand-entered coordinates mostly do NOT
    reproduce against real precomputed SpliceAI (coordinate errors + limited
    deep-intronic coverage in the precomputed release).

    LICENSE: SpliceAI scores are CC BY-NC 4.0 (Jaganathan et al. 2019, PMID
    30661751). The 28.6 GB source VCF stays external; cite SpliceAI + Illumina.
    """
    if not demo and SPLICEAI_CSV.exists():
        df = pd.read_csv(SPLICEAI_CSV)
        df["source"] = "REAL"
        return df
    if not demo:
        _missing_extract("SpliceAI", SPLICEAI_CSV, strict)
    return load_splice_demo()


def load_pangolin(demo: bool = False, strict: bool = False) -> pd.DataFrame:
    """Pangolin splice scores for CFTR — REAL if you have run the model.

    Pangolin (Zeng & Li 2022, Genome Biol 23:103, PMID 35449021,
    github.com/tkzeng/Pangolin) has no precomputed per-gene release and is not in
    dbNSFP, so real scores require RUNNING the model: the build cell in
    tools/08_pangolin.ipynb does that locally (weights bundled with the pip
    package; only the ~215 kb CFTR reference region is needed, no whole-genome
    download). Score 0-1, >= 0.5 high, >= 0.2 moderate; keyed by genomic coordinate.

    REAL if the extract exists (``data/pangolin_cftr.csv``). The build cell's
    default ``SCOPE = "cftr2"`` scores every CFTR2 variant carrying GRCh38
    coordinates — ~1,892 of 2,097, **SNVs and indels** — and labels the result
    ``source='REAL'``; ``SCOPE = "curated"`` scores just the 5 classic splice alleles
    and stays ``source='DEMO'``, because that coverage is a teaching subset rather
    than a worklist. The label follows the *scope*, never the model. Rows that could
    not be scored are kept with an empty score and a ``skip_reason``.

    Because Pangolin is run locally it reaches variant classes the precomputed
    masked-SNV SpliceAI release cannot (notably indels) — but see predict/13: a
    Pangolin score on a frameshift is a *splice* verdict, and "no splice impact" is
    usually correct and rarely the reason the variant is pathogenic.

    If the file is absent this falls back to the hand-authored splice-demo pangolin
    values (with a warning; strict=True raises).
    """
    if not demo and PANGOLIN_CSV.exists():
        df = pd.read_csv(PANGOLIN_CSV)
        return df  # source column set by the build cell (DEMO for curated scope)
    if not demo:
        _missing_extract("Pangolin", PANGOLIN_CSV, strict)
    d = load_splice_demo()
    return d[["variant_id", "hgvs_c", "legacy_name", "pangolin_score", "source"]]


# =============================================================================
# 6. Shared curated "worked-example" panels (A1 missense / A2 splice)
# =============================================================================
# A fixed panel of famous CFTR variants scored by EVERY tool, so a reader can
# follow the SAME variants through all the per-tool notebooks. Each notebook
# scores this panel itself, inline, with its own loader -- see e.g.
# tools/02_alphamissense.ipynb's "shared missense worked-example panel" section.

# The shared A1 panel membership (famous CFTR missense variants, 1-letter keys) —
# single source of truth for WHICH variants; each notebook shows its own SCORING.
A1_PANEL_VARIANTS = [three_to_one(v[0]) for v in _DEMO_MISSENSE]

# Known CF splice alleles for the A2 panel — looked up by CFTR2 cDNA name so each
# notebook can use the AUTHORITATIVE GRCh38 coordinates (not the demo's hand-entered ones).
A2_KNOWN_CDNA = ["c.2988+1G>A", "c.2657+5G>A", "c.3718-2477C>T", "c.3140-26A>G", "c.1680-886A>G"]


# =============================================================================
# 7. Scoring helpers — turn a raw score into a 3-class call
# =============================================================================
def call_from_score(score, tool: str) -> str:
    """Map a raw score to {'pathogenic','uncertain','benign'} using THRESHOLDS.

    Respects each tool's direction (ESM1b: lower = worse). Returns 'na' if the
    score is missing. This is the single place the binary cut-points live, so a
    notebook can change a threshold and re-score consistently.
    """
    if score is None or (isinstance(score, float) and np.isnan(score)):
        return "na"
    t = THRESHOLDS[tool]
    if HIGHER_IS_WORSE.get(tool, True):
        if score >= t["path"]:
            return "pathogenic"
        if score < t["benign"]:
            return "benign"
        return "uncertain"
    else:                                   # lower = worse (ESM1b)
        if score <= t["path"]:
            return "pathogenic"
        return "benign" if score > t["benign"] else "uncertain"


def cv_class(sig) -> str:
    """Collapse ClinVar free-text significance to pathogenic/benign/uncertain."""
    if sig is None or (isinstance(sig, float) and np.isnan(sig)):
        return "unknown"
    s = str(sig).lower()
    if "conflicting" in s:
        return "uncertain"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    return "uncertain"


# ─────────────────────────────────────────────────────────────────────────────
# Amino-acid helpers (HGVS 3-letter → 1-letter short key, e.g. p.Gly551Asp→G551D)
# ─────────────────────────────────────────────────────────────────────────────
_AA3_TO_1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Glu":"E","Gln":"Q",
             "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
             "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}


def hgvsp_to_short(hgvsp) -> str | None:
    """p.Gly551Asp -> G551D  (the join key AlphaMissense/EVE/ESM1b use)."""
    if not hgvsp or (isinstance(hgvsp, float) and np.isnan(hgvsp)):
        return None
    s = str(hgvsp).replace("p.", "").replace("(", "").replace(")", "")
    for three, one in _AA3_TO_1.items():
        s = s.replace(three, one)
    return s or None


def extract_hgvsp_from_name(name) -> str | None:
    """Pull the protein change out of a ClinVar 'Name' field → short key."""
    import re
    if name is None:
        return None
    m = re.search(r"p\.([A-Z][a-z]{2}\d+[A-Z][a-z]{2}|[A-Z][a-z]{2}\d+Ter)", str(name))
    return hgvsp_to_short("p." + m.group(1)) if m else None


if __name__ == "__main__":
    # smoke test — prints each loader's row count and its ACTUAL source column
    # (REAL when the extract/cache exists, DEMO fallback otherwise), so the label
    # can never disagree with what was loaded. The three live-fetched loaders
    # (gnomAD, AlphaMissense, ClinVar) raise if data/ is missing; catch that so
    # the smoke test still completes.
    def _src(df):
        return df["source"].iloc[0] if len(df) and "source" in df else "n/a"

    def _try(label, fn):
        try:
            df = fn()
            print(f"{label:16}: {len(df):>7} rows  source={_src(df)}")
        except FileNotFoundError as exc:
            print(f"{label:16}: MISSING — {str(exc).splitlines()[0]}")

    _try("gnomAD missense", load_gnomad_missense)
    _try("gnomAD noncoding", load_gnomad_noncoding)
    _try("AlphaMissense", load_alphamissense)
    _try("ClinVar", load_clinvar)
    _try("EVE", load_eve)          # REAL if data/eve_cftr_2021-08.csv exists, else DEMO
    _try("ESM1b", load_esm1b)
    _try("REVEL", load_revel)
    _try("PrimateAI", load_primateai)
    _try("CFTR2", load_cftr2)
    _try("SpliceAI", load_spliceai)
    _try("Pangolin", load_pangolin)
    _try("splice demo", load_splice_demo)
    # CADD has no thin reader here -- fetch_cadd() now lives in tools/09_cadd.ipynb
    # (it needs no local file at all, so there was nothing left to smoke-test here).
