# Lay-reader review: unclear explanations and logic breakpoints

Perspective: I read this repository as a new learner in bioinformatics, CFTR, and variant-effect prediction. This is not a scientific validation; it is a clarity and internal-consistency review focused on places where a beginner could get confused or where the repository's stated logic appears to break.

## High-impact clarity issues

1. **The repository repeatedly says real data files are committed, but they are absent from the checkout.**
   - `README.md` says CFTR2 ships at `data/cftr2_2026-01-30.csv` and lists a `data/` directory.
   - `data_manifest.json` lists multiple `derived_extract` CSVs under `data/`.
   - In this checkout, there is no `data/` directory at repository depth 2, so the default loaders for CFTR2, EVE, ESM1b, REVEL, PrimateAI, and SpliceAI will silently fall back to demo data or fail depending on the function.
   - Beginner impact: the reader cannot tell whether they are reading a data-included toolkit, a code-only toolkit, or a notebook-output-only artifact.

2. **The REAL/DEMO story contradicts itself.**
   - The top table says every predictor except Pangolin is REAL.
   - The A1 summary immediately says the five-predictor framing uses AlphaMissense real data and other predictors are demo on about 13 variants.
   - Later, the explanation says all five predictors are now REAL and the real priority worklist is 473.
   - The Known limitations section still says four of five missense predictors and both splice predictors are DEMO.
   - Beginner impact: I cannot tell which numerical claims are current, which are historical, and which are intentionally preserved for provenance.

3. **Several headline counts are hard to reconcile.**
   - gnomAD missense is described as about 2,466 in the REAL/DEMO table, but A1 says 2,496 missense variants scored after adding about 30 curated variants.
   - Splice variants are described as about 1,085 non-coding in the REAL/DEMO table, 1,094 in the dashboard, 4,717 in `data_manifest.json`, and about 4,260 observed gnomAD non-coding variants with SpliceAI scores later in the README.
   - Beginner impact: without a short count glossary, it is unclear whether these are raw variants, filtered variants, scored variants, observed variants, or historical webpage numbers.

4. **The README mixes reproduced historical results with current results without a visual boundary.**
   - Some tables are “kept here as reproduced summary,” while later paragraphs update them.
   - The reader needs a clear label such as “Historical/demo reproduced numbers” versus “Current real-data rerun numbers.”
   - Beginner impact: it is easy to cite the wrong number, especially the Priority 1 count of 4 versus 3 versus 473.

5. **`toolkit.py` still contains stale comments and smoke-test labels.**
   - Section 5 is titled “DEMO predictors — EVE / ESM1b / REVEL / PrimateAI (NOT real predictions),” but those loaders are now documented as REAL by default when files exist.
   - The `__main__` smoke test prints `EVE demo` while calling `load_eve()` with default `demo=False`.
   - Beginner impact: the code comments disagree with the loader behavior, so a learner cannot confidently infer what will run.

## Bioinformatics concepts that need more lay explanation

6. **“Missense,” “splice,” “synonymous,” “deep intronic,” and “non-coding” need a one-paragraph glossary near the top.**
   - The notebooks may explain them individually, but the README uses them before defining them.
   - Suggested beginner framing: missense changes one amino acid; splice variants affect exon/intron cutting and joining; synonymous changes do not change the amino acid but can still affect splicing; deep intronic variants sit far inside introns but can create cryptic splice sites.

7. **The CFTR strand issue is easy to miss.**
   - `fetch_cadd()` mentions that CFTR is on the minus strand and tries complements.
   - A reader needs a central explanation that genomic `ref/alt` may look like the complement of the cDNA/HGVS coding change.
   - This matters because several tools are keyed by protein variant, while others are keyed by genomic coordinate.

8. **Join keys are a recurring source of hidden complexity.**
   - Some tables join by `protein_variant` like `G551D`; others join by `chrom,pos,ref,alt`; CFTR2 has cDNA names, legacy names, and sometimes protein keys.
   - The repo would benefit from a “join key map” table: source, key type, example, gotcha.
   - Beginner impact: many biological disagreements may actually be data-joining mismatches.

9. **“Orthogonal evidence,” “circularity,” and “label leakage” are important but abstract.**
   - The README mentions them often, but a lay reader may not understand why REVEL-vs-ClinVar can be circular while AlphaMissense-vs-ClinVar is less circular.
   - Suggested short analogy: testing a student with questions they studied from is not independent evidence.

10. **Thresholds are presented as if binary, then caveated later.**
    - The table gives single pathogenic cutoffs, while text later says ACMG uses graded thresholds.
    - Beginner impact: a reader may treat `score >= cutoff` as a diagnosis.
    - Suggest explicitly saying “these thresholds create teaching worklists, not clinical classifications.”

## Logic breakpoints or implementation surprises

11. **Missing data can cause silent DEMO fallback for some loaders.**
    - `load_eve()`, `load_esm1b()`, `load_revel()`, `load_primateai()`, `load_cftr2()`, and `load_spliceai()` return demo data if their real CSV is missing.
    - That is convenient for teaching, but it can hide a broken data setup.
    - Suggested improvement: add a `strict=True` option or print/warn when default real loading falls back to demo.

12. **The README says every returned DataFrame carries `source`, but real/demo fallback can make coverage ambiguous.**
    - `source` marks rows as REAL or DEMO, which is good.
    - However, a beginner may not inspect the column and may only see a successful notebook execution.
    - Suggested improvement: notebooks should display a `source` count after every load.

13. **ClinVar release is unpinned in the manifest.**
    - `data_manifest.json` explicitly says the ClinVar release date is unpinned.
    - Because ClinVar changes over time, A1 VUS and discordance counts can drift.
    - Beginner impact: results may not be reproducible even if the code is.

14. **CADD is live API by default, so it is not fully reproducible.**
    - The manifest notes that CADD responses should be cached.
    - Beginner impact: a notebook can change or fail due to network/API behavior rather than biology.

15. **Pangolin remains demo, but it is placed beside real SpliceAI in A2 language.**
    - The README says A2 uses SpliceAI + Pangolin + CADD, but Pangolin is still demo.
    - Suggested wording: “real SpliceAI, demo Pangolin, live CADD unless cached.”

16. **The synthetic splice example is too close to the real worklist language.**
    - `c.2657+120C>T` is described as a splice-risk VUS in the dashboard and later as a synthetic teaching example.
    - Beginner impact: it can be mistaken for a real patient or database variant.
    - Suggest marking synthetic variants directly in every table where they appear.

17. **The stated repository layout includes `outputs/` and `data/`, but this checkout does not include them.**
    - The README says notebook 12 writes outputs and references `outputs/A1_upgrade_worklist_REAL.csv` and `outputs/A2_splice_DEMO.csv`.
    - In this checkout, those directories/files are absent.
    - Beginner impact: source citations in the README point to files that are not available.

18. **Notebook outputs may contain the only readable results, but source data is absent.**
    - The README says notebooks keep outputs so readers can read results.
    - A learner who wants to verify a table cannot rerun without external downloads and cannot inspect the referenced output CSVs.
    - Suggested improvement: commit small non-restricted summary CSVs or a frozen report with only aggregate counts and non-licensed examples.

## Suggested beginner-facing additions

1. Add a “current status in this checkout” section that says exactly which files are present, which must be downloaded, and which loaders will fall back to demo.
2. Add a count glossary explaining every major number: 2,466, 2,496, 413, 473, 1,085, 1,094, 4,260, 4,717, 566k.
3. Add a keying/joining glossary for `protein_variant`, HGVS c., HGVS p., legacy CFTR names, and genomic coordinates.
4. Add one warning box: “Predictor score does not equal clinical diagnosis.”
5. Add a strict data-loading mode or warnings for demo fallback.
6. Clearly split historical/demo reproduced results from current real-data rerun results.
