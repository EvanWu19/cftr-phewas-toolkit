"""
cftr_variant_toolkit / toolkit.py
=================================
A small, heavily-documented library of loaders for the variant-effect
prediction tools used in the CFTR PheWAS variant-interpretation work
(the A1 missense-triage and A2 splice-discordance analyses).

Design goals
------------
1.  **Honest about provenance.** Every DataFrame this module returns carries a
    ``source`` column whose value is either ``"REAL"`` (downloaded from the
    primary source / queried live) or ``"DEMO"`` (a small, hand-curated table of
    illustrative values baked into this file). Never mix them silently.

2.  **Reuse the cache.** The heavy real datasets (AlphaMissense genome-wide,
    ClinVar variant_summary, gnomAD queries) were already downloaded into the
    project ``_tmp_fetch/`` cache. The ``load_*`` functions read those cached
    extracts so the notebooks produce genuine results without re-downloading
    ~1 GB. If the cache is missing, each function's docstring tells you exactly
    how to regenerate it.

3.  **Beginner-readable.** Each function documents *what the tool is*, *what the
    score means*, *the threshold and why*, and *the primary reference*.

Real vs demo, as shipped
------------------------
    REAL (cached extracts of primary data):
        gnomAD v4 missense + non-coding, AlphaMissense (CFTR), ClinVar (CFTR)
    REAL (shipped with the toolkit — small, public, static):
        CFTR2 variant list (30 January 2026 release, ~2,097 variants)
    REAL (queried live per-call):
        CADD v1.7 REST API
    DEMO (hand-curated illustrative values — NOT real predictions):
        EVE, ESM1b, REVEL, PrimateAI, SpliceAI, Pangolin

    => The only genome-wide REAL predictor here is AlphaMissense. EVE/ESM1b/
       REVEL/PrimateAI cover only the ~30 curated demo variants, and SpliceAI/
       Pangolin only the 9 curated splice variants. See README.md.

References
----------
    AlphaMissense : Cheng et al. 2023   Science    PMID 37733863
    EVE           : Frazer et al. 2021  Nature     PMID 34707284
    ESM1b         : Brandes et al. 2023 Nat Genet  PMID 37563329
    REVEL         : Ioannidis et al. 2016 AJHG     PMID 27666373
    PrimateAI     : Sundaram et al. 2018 Nat Genet PMID 30038395
    SpliceAI      : Jaganathan et al. 2019 Cell    PMID 30661751
    Pangolin      : Zeng & Bromberg 2021 Genome Biol PMID 34753498
    CADD-Splice   : Rentzsch et al. 2021 Genome Med  PMID 33618777
    REVEL thresholds : Pejaver et al. 2022 AJHG (ACMG calibration) PMID 36413997
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
# toolkit.py lives in  <project>/05_analysis/cftr_variant_toolkit/
# the shared download cache lives in  <project>/_tmp_fetch/
PKG_DIR    = Path(__file__).resolve().parent
PROJECT    = PKG_DIR.parents[1]
CACHE_DIR  = PROJECT / "_tmp_fetch"
OUT_DIR    = PKG_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)
# REAL CFTR2 release, shipped with the toolkit (small, public, static). Built
# from the official cftr2.org variant list by build_cftr2.py.
CFTR2_CSV  = PKG_DIR / "data" / "cftr2_2026-01-30.csv"

# CFTR locus, GRCh38. CFTR is on the minus strand of chromosome 7.
CFTR_CHR, CFTR_START, CFTR_END = "7", 117_470_098, 117_667_108
CFTR_UNIPROT = "P13569"            # AlphaMissense protein key
CFTR_MANE_TX = "NM_000492.4"       # MANE Select transcript


# ─────────────────────────────────────────────────────────────────────────────
# Published decision thresholds (one binary cut per tool)
# ─────────────────────────────────────────────────────────────────────────────
# NB: these are deliberately *simple* single cut-points taken from each tool's
# calibration paper. Real ACMG use (esp. REVEL) applies GRADED thresholds for
# different evidence strengths — see notebook 13 and Pejaver 2022.
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
        circularity="low-for-clinvar", pmid="34753498"),
    "CADD": dict(kind="general", learning="semi-supervised",
        signal="SVM/logistic on many annotations; trained on simulated-vs-observed variants",
        circularity="medium", pmid="33618777"),
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

    Cache: _tmp_fetch/gnomad_cftr_missense.tsv
    Regenerate: gnomAD GraphQL API, gene ENSG00000001626, dataset gnomad_r4,
      keep consequence == 'missense_variant'. (See fetch_scores.fetch_gnomad_cftr)
    """
    fp = CACHE_DIR / "gnomad_cftr_missense.tsv"
    if not fp.exists():
        raise FileNotFoundError(
            f"{fp} missing. Query the gnomAD v4 GraphQL API for CFTR "
            "(ENSG00000001626) missense variants — see docstring.")
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

    Cache: _tmp_fetch/gnomad_cftr_noncoding.tsv  (~1,085 variants)
    """
    fp = CACHE_DIR / "gnomad_cftr_noncoding.tsv"
    if not fp.exists():
        raise FileNotFoundError(f"{fp} missing — see docstring / fetch_splice_scores.py")
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
    *against* ClinVar without circular reasoning (see notebook 13).

    Score `am_pathogenicity` in [0,1]; AlphaMissense's own 3-class cut-points:
        >= 0.564  -> "likely_pathogenic"
        <= 0.340  -> "likely_benign"
        else      -> "ambiguous"

    Returns: protein_variant, am_score, am_class, source (~2,460 rows).
    Cache: _tmp_fetch/alphamissense_cftr.tsv  (extracted from the 0.6 GB
      AlphaMissense_hg38.tsv.gz, filtered to UniProt P13569).
    Regenerate: tabix the genome-wide AlphaMissense_hg38.tsv.gz to the CFTR
      region and keep uniprot_id == 'P13569'.
    """
    fp = CACHE_DIR / "alphamissense_cftr.tsv"
    if not fp.exists():
        raise FileNotFoundError(
            f"{fp} missing. Download AlphaMissense_hg38.tsv.gz and filter to "
            "UniProt P13569 — see docstring.")
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
    comparing those predictors to ClinVar is partly circular (notebook 13).

    Returns: protein_variant, clinvar_sig, review_status, clinvar_call, source.
    `clinvar_call` collapses the free-text significance to pathogenic/benign/
    uncertain via cv_class().
    Cache: _tmp_fetch/clinvar_cftr.tsv  (filtered from ClinVar variant_summary).
    """
    fp = CACHE_DIR / "clinvar_cftr.tsv"
    if not fp.exists():
        raise FileNotFoundError(f"{fp} missing — filter ClinVar variant_summary.txt.gz to CFTR.")
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
# 4. CADD — live REST API (REAL, per-variant)
# =============================================================================
def fetch_cadd(chrom: str, pos: int, ref: str, alt: str,
               delay_sec: float = 0.3) -> dict:
    """Score ONE variant with the CADD v1.7 REST API (REAL, live).

    CADD (Combined Annotation Dependent Depletion) integrates dozens of
    annotations into one deleteriousness score. The PHRED-scaled value is the
    handy one: PHRED >= 15 means "top ~3% most deleterious of all possible
    variants"; >= 20 means top 1%. CADD v1.7 folds in some splice features.

    NOTE ON STRAND: CFTR is on the genomic minus strand, so a coding-strand
    change (e.g. C>T) appears on the plus strand as its complement (G>A). This
    helper tries both orientations.

    API: https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chr}:{pos}-{pos}
    Returns dict(cadd_raw, cadd_phred) or None values on miss.
    """
    url = f"https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chrom}:{pos}-{pos}"
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    try:
        data = requests.get(url, timeout=15).json()
    except Exception as exc:                                    # network / JSON
        return {"cadd_raw": None, "cadd_phred": None, "error": str(exc)}
    header = data[0] if data else []
    for rec in data[1:]:
        if len(rec) < 6:
            continue
        r, a = rec[2], rec[3]
        if (r == ref and a == alt) or (r == comp.get(ref) and a == comp.get(alt)):
            time.sleep(delay_sec)
            return {"cadd_raw": float(rec[4]), "cadd_phred": float(rec[5])}
    time.sleep(delay_sec)
    return {"cadd_raw": None, "cadd_phred": None}


# =============================================================================
# 5. DEMO predictors — EVE / ESM1b / REVEL / PrimateAI  (NOT real predictions)
# =============================================================================
# These four tools are gated (login) or heavy (GPU / large files) to obtain, so
# the toolkit ships a tiny curated table of ~30 illustrative variants. The
# `source` column is "DEMO" so these can never be mistaken for real output.
#
# To use REAL scores instead, download the per-tool CFTR files and join on
# protein_variant. Each function's docstring gives the URL / method.
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


def load_eve(demo: bool = True) -> pd.DataFrame:
    """EVE score per CFTR missense variant.

    EVE (Evolutionary model of Variant Effect, Frazer 2021) is an UNSUPERVISED
    deep generative model trained on the multiple-sequence alignment of the
    protein family — it scores how well a variant fits evolutionary constraint.
    Score in [0,1]; >= 0.5 ~ pathogenic. No clinical labels used → low
    circularity vs ClinVar.

    demo=True (default) returns the curated DEMO table (source='DEMO').
    REAL: download CFTR (P13569) scores from https://evemodel.org (per-protein
    CSV) and join on protein_variant; set demo=False and point EVE_CSV at it.
    """
    d = _demo_frame()
    return d[["protein_variant", "eve_score", "source"]].dropna(subset=["eve_score"])


def load_esm1b(demo: bool = True) -> pd.DataFrame:
    """ESM1b LLR per CFTR missense variant.

    ESM1b (Brandes 2023) is a protein LANGUAGE model. It scores a variant by the
    log-likelihood ratio (LLR) of the mutant vs wild-type amino acid — a more
    NEGATIVE LLR means the model finds the mutation more surprising/damaging.
    Cut: LLR <= -7.5 ~ pathogenic. Unsupervised → low circularity.

    REAL: bulk ESM1b variant-effect files are published per-protein (Brandes
    2023 supplement / HuggingFace); join on protein_variant.
    """
    d = _demo_frame()
    return d[["protein_variant", "esm1b_score", "source"]].dropna(subset=["esm1b_score"])


def load_revel(demo: bool = True) -> pd.DataFrame:
    """REVEL score per CFTR missense variant.

    REVEL (Ioannidis 2016) is a SUPERVISED random-forest ENSEMBLE of 13 other
    predictors, trained on curated pathogenic/benign variants. Score in [0,1];
    the ACMG calibration (Pejaver 2022) gives GRADED cut-points, of which 0.75
    is a common single 'likely pathogenic' point.

    ⚠ CIRCULARITY: REVEL's training labels share lineage with ClinVar/HGMD, so
    'REVEL disagrees with ClinVar' can partly reflect label leakage, not
    independent evidence. Handle in notebook 13.

    REAL: download the genome-wide REVEL table (sites.google.com/site/revelgenomics)
    and join on genomic coordinate (chr,pos,ref,alt).
    """
    d = _demo_frame()
    return d[["protein_variant", "revel_score", "source"]].dropna(subset=["revel_score"])


def load_primateai(demo: bool = True) -> pd.DataFrame:
    """PrimateAI score per CFTR missense variant.

    PrimateAI (Sundaram 2018) is a deep net trained SEMI-supervised on common
    human & non-human primate missense variants as a proxy for benignity.
    Score in [0,1]; >= 0.803 ~ pathogenic. Medium circularity.

    REAL: PrimateAI-3D / PrimateAI scores are distributed by Illumina; join on
    genomic coordinate.
    """
    d = _demo_frame()
    return d[["protein_variant", "primate_ai_score", "source"]].dropna(subset=["primate_ai_score"])


def load_cftr2(demo: bool = False) -> pd.DataFrame:
    """CFTR2 clinical-functional class per CFTR variant — REAL by default.

    CFTR2 (cftr2.org) is the clinical-functional reference for CF: it labels
    variants as 'CF-causing', 'Varying clinical consequence', 'Non CF-causing',
    or 'No interpretation available', based on patient data + in-vitro CFTR
    function. Because CFTR2's call includes FUNCTIONAL assay evidence, it is a
    useful ORTHOGONAL truth set for benchmarking sequence predictors (notebook
    08) — less circular than ClinVar for supervised tools.

    REAL (default): the full public CFTR2 variant list (30 January 2026 release,
    ~2,097 variants) shipped at ``data/cftr2_2026-01-30.csv`` and built from the
    official cftr2.org download by ``build_cftr2.py``. Returns a 1-letter
    ``protein_variant`` key (e.g. 'G551D') for the ~780 simple-missense variants
    so it joins onto the AlphaMissense/gnomAD tables; non-missense rows carry an
    empty key but keep their legacy/cDNA names and genomic coordinates.

    demo=True returns the tiny embedded curated set (source='DEMO') used by the
    early teaching cells.

    NOTE: CFTR2 data is redistributed here under CFTR2's public data-use terms —
    please cite CFTR2 (cftr2.org) if you use it.
    """
    if not demo and CFTR2_CSV.exists():
        df = pd.read_csv(CFTR2_CSV, dtype={"protein_variant": "string"})
        df["protein_variant"] = df["protein_variant"].fillna("")
        df["source"] = "REAL"
        return df
    d = _demo_frame()
    return d[["protein_variant", "cftr2_class", "source"]]


def load_cftr2_demo() -> pd.DataFrame:
    """Tiny embedded CFTR2 class table (DEMO). Prefer load_cftr2() for the full
    REAL list; kept for the early teaching cells that predate the real loader."""
    d = _demo_frame()
    return d[["protein_variant", "cftr2_class", "source"]]


# =============================================================================
# 6. SpliceAI / Pangolin — DEMO curated splice variants
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
    impact, >= 0.2 = moderate. Pangolin (Zeng 2021) is a similar model giving
    one 0–1 score.

    ⚠ The DS_/pangolin numbers here are DEMO (hand-authored). For real scores:
    SpliceAI — precomputed VCF from Illumina BaseSpace (login) or the Broad
    SpliceAI-lookup app; Pangolin — run the model locally (GPU). See README.
    """
    df = pd.DataFrame(KNOWN_CF_SPLICE_VARIANTS, columns=_SPLICE_COLS)
    df["spliceai_ds_max"] = df[["DS_AG", "DS_AL", "DS_DG", "DS_DL"]].max(axis=1)
    df["source"] = "DEMO"
    return df


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
    # smoke test — prints how many REAL vs DEMO rows each loader yields
    print("gnomAD missense :", len(load_gnomad_missense()), "(REAL)")
    print("gnomAD noncoding:", len(load_gnomad_noncoding()), "(REAL)")
    print("AlphaMissense   :", len(load_alphamissense()), "(REAL)")
    print("ClinVar         :", len(load_clinvar()), "(REAL)")
    print("EVE demo        :", len(load_eve()), "(DEMO)")
    print("splice demo     :", len(load_splice_demo()), "(DEMO)")
    print("CADD live 2988+1:", fetch_cadd("7", 117_592_260, "C", "T"))
