# Coffee Consumption and Atrial Fibrillation — Analysis Code

Reproducible analysis code for the integrated Mendelian randomization (MR) + FAERS
pharmacovigilance study of coffee intake and atrial fibrillation (AF).

All analyses use **publicly available data only**. No patient-level data were accessed.

## Environment

- R 4.6.0 (Windows x86_64). Required packages:
  `TwoSampleMR`, `MRPRESSO`, `MVMR`, `RadialMR`, `ieugwasr`, `dplyr`, `ggplot2`.
  - `TwoSampleMR`, `MRPRESSO`, `ieugwasr`: CRAN.
  - `MVMR`: `remotes::install_github("WSpiller/MVMR")`
  - `RadialMR`: `remotes::install_github("WSpiller/RadialMR")`
- Python 3 (managed) with `pandas`, `requests` for FAERS scripts.
- OpenGWAS API token required (free): register at https://api.opengwas.io/profile/
  and replace the placeholder `jwt <- "YOUR_OPENGWAS_TOKEN"` at the top of each
  R script that queries OpenGWAS. Tokens expire ~13 days after issue.

## Data sources

| Role | Source | Access |
|---|---|---|
| Exposure (primary) | Coffee intake, cups/day — UK Biobank, OpenGWAS `ukb-b-5237` (n = 428,860) | OpenGWAS API |
| Exposure (sensitivity) | Coffee intake — PAGE / Wojcik 2019, OpenGWAS `ebi-a-GCST008028` (n = 15,837, multi-ancestry) | OpenGWAS API |
| Outcome (primary) | AF — Nielsen et al. 2018, OpenGWAS `ebi-a-GCST006414` (n = 1,030,836) | OpenGWAS API |
| Outcome (replication) | FinnGen R11 `I9_AF` (55,853 cases / 231,952 controls) | https://r11.finngen.fi/ (free registration) |
| Outcome (cross-ancestry replication) | AF and atrial flutter — BioBank Japan (Sakaue & Kanai 2021; 4,150 cases / 155,540 controls, East Asian; GRCh37) | pheweb.jp, file `hum0197.v3.BBJ.AF.v1` (redirects to humandbs.dbcls.jp) |
| Mediators (MVMR) | SBP, BMI, smoking initiation, total cholesterol (OpenGWAS IDs inside `mr_mvmr*.R`) | OpenGWAS API |
| Pharmacovigilance | FDA Adverse Event Reporting System via openFDA API | https://api.fda.gov/drug/event.json |
| Outcome (spectrum / sensitivity) | FinnGen R12 `I9_AF` and arrhythmia-spectrum endpoints (`CARDIAC_ARRHYTM`, `I9_PAROXTAC`, `I9_OTHARR`, `I9_AVBLOCK`, `I9_CONDUCTIO`) | https://r12.finngen.fi/ REST API |
| Outcome (sensitivity) | FinnGen R13 `I9_AF` summary statistics | https://www.finngen.fi/en/access_results |
| Mediators (metabolites) | Chen et al. 2023 plasma metabolome GWAS — 9 caffeine-family metabolites (GCST90199644–GCST90200436) | EBI GWAS Catalog FTP (harmonised files) |
| Intermediate phenotypes | QT interval (GCST90165290), resting heart rate (GCST007609), PR interval (GCST010321) | EBI GWAS Catalog FTP (harmonised files) |

Place the FinnGen file as `finngen_R11_I9_AF.gz` in the working directory
(see download URL inside `mr_r11.R`).

## Pipeline order

**Univariable MR (primary)**
1. `mr_main.R` — instrument selection (P < 5e-8, clumped), harmonisation, IVW /
   weighted-median / MR-Egger / MR-PRESSO vs Nielsen outcome.
2. `mr_wave2.R`, `mr_wave2b.R` — instrument/robustness extensions (leave-one-out,
   heterogeneity).
3. `mr_wave3.R`, `mr_wave3b.R` — additional robustness and plotting.
4. `mr_r11.R` — replication vs FinnGen R11.
5. `mr_step1.R`, `mr_step1b.R`, `mr_p0_fixes.R` — exploratory/first-pass runs.

**Multivariable MR**
6. `mr_mvmr_exp.R`, `mr_mvmr_outcome.R`, `mr_mvmr.R` — MVMR vs Nielsen outcome
   (SBP, BMI, smoking initiation, TC).
7. `mr_mvmr_r11.R` — MVMR vs FinnGen R11.
8. `condf.R` — conditional F-statistics for MVMR instruments.

**P2 robustness**
9. `radial.R`, `radial2.R` — Radial MR (modified second-order weights) vs both
   outcomes; outlier detection.
10. `wojcik.R`, `wojcik_r11.R` — independent-instrument sensitivity using the
    PAGE/Wojcik exposure (P < 1e-6 set) vs both outcomes.
11. Equivalence testing (TOST, bounds OR 0.80–1.25) and power analysis were
    computed from the harmonised estimates above (script snippets embedded in
    the supplementary material, Table S8).

**Pharmacovigilance**
12. `faers_rebuild.py` — openFDA query of caffeine-containing products and AF
    reports; PRR / ROR / IC computation.
13. `faers_sens_a.py` — FAERS sensitivity analyses (brand-name phrasing,
    seriousness restriction).

**Colocalization and PheWAS** (`colocalization/`)
14. `fetch_api.py` — batch OpenGWAS association lookup for window variants
    (coffee intake) and per-SNP phenome-wide scans (PheWAS, 40 instruments).
15. `fetch_ahr.py` — OpenGWAS lookup for the post hoc AHR region (rs4410790
    ±250 kb; rsid list in `data/ahr_rsids.txt`).
16. `coloc_run.R` — coloc.abf at the outlier-locus windows (coffee intake vs
    AF, FinnGen R11 and Nielsen), default priors.
17. `coloc_prior_sensitivity.R` — prior-sensitivity re-fits
    (p12 = 5e-6 / 1e-5 / 5e-5).
18. `coloc_ahr.R` — post hoc AHR-region colocalization across priors.

**Cross-ancestry replication** (`bbj_replication/`)
19. Download the BBJ AF summary statistics (`GWASsummary_Atrial_Flutter_Japanese_SakaueKanai2020.auto.txt.gz`
    inside `hum0197.v3.BBJ.AF.v1.zip`); harmonise the 40 instruments by rsID
    (34/40 matched: 33 by rsID + rs442355 by position; 6 unavailable,
    see `why_missing.py` output and manuscript Supplementary Table S20);
    align BBJ BETA (Allele2) to the coffee-intake effect allele; compute
    the Wald ratio theta = by/bx with weights w = (bx/sy)^2.
20. `fix_weights_bbj.py` — recomputes IVW (multiplicative random-effects),
    Cochran's Q / I², weighted MR-Egger intercept, and the three-dataset
    pooled estimate from `data/bbj_harmonized.csv`.
21. `verify_bbj.py` — independent recomputation plus single-variant lookup
    of rs762551 / rs4410790 in the raw BBJ file.
22. `why_missing.py` — coordinate-level scan documenting why 6 variants
    are absent from the BBJ release (point-level imputation-QC filtering,
    not regional coverage gaps).

**Figures**
23. `mr_figures_final.R`, `mr_figures_r11.R` — forest/scatter/funnel/leave-one-out
    figures (300 dpi TIFF).

**2026-09 revision: mechanism, spectrum and synthesis analyses**
(`ep_20260901/`, `spectrum_20260901/`, `r13/`, `review_20260901/`)
24. `ep_20260901/cluster_mr.py` — mechanism-cluster assignment of the 40
    instruments (M = caffeine pharmacokinetics, A = adiposity/intake
    propensity, O = other/unassigned) and cluster-stratified IVW estimates
    (`cluster_assignments.csv`, `cluster_mr_results.csv`).
25. `ep_20260901/ep_mr.py` — MR on electrophysiological intermediate
    phenotypes (PR interval, QT interval, resting heart rate)
    (`ep_mr_results.csv`).
26. `ep_20260901/query_qt.py` — extraction of the 40 instruments from the
    QT-interval GWAS via the shared pure-Python remote tabix client
    (`ep_20260901/remote_tabix.py`; no pysam/htslib required).
27. `ep_20260901/metabolite_panel.py`, `ep_20260901/analyze_metabolites.py` —
    extraction of the instruments from 9 caffeine-family metabolite GWAS and
    overall/cluster-stratified IVW with BH FDR (`metab_results.csv`).
28. `spectrum_20260901/fetch_spectrum.py` — FinnGen R12 REST extraction of the
    40 instruments across the arrhythmia endpoints;
    `spectrum_20260901/fig_spectrum.py` — spectrum forest figure. Harmonised
    per-SNP dataset: `spectrum_20260901/spectrum_harmonized_snps.csv`; pooled
    estimates: `spectrum_20260901/spectrum_mr_results.csv`.
29. `r13/run_r13.py` — FinnGen R13 I9_AF sensitivity: rsID harmonisation of
    the 40 instruments, IVW (multiplicative random-effects) and weighted
    MR-Egger.
30. `review_20260901/reproduce_recurrence_synthesis.py` — fixed-effects
    pooling of the DECAF hazard ratio and the post-PVI adjusted odds ratio
    (output `recurrence_synthesis.json`; inputs in manuscript Supplementary
    Table S26). `params.json` stores the PubMed query used for the 2026
    literature update; raw API dumps are not redistributed.

**Manuscript audit** (`audit/`)
31. `audit/final_audit.py` — automated consistency battery (citation order,
    key-number recurrence, word counts, section numbering) run against the
    manuscript file.
32. `audit/final_audit_fixes.py`, `audit/conclusions_certainty.py` — scripted
    text-level revisions applied during finalisation (kept for provenance).

## Notes

- R scripts in this environment occasionally print a harmless segfault after the
  final `DONE` line (Windows/R 4.6.0 shutdown issue); all outputs are written
  before that point.
- `MRlap` (sample-overlap correction) was attempted but is not installable in
  this environment (requires Rtools + `GenomicSEM` compilation). Overlap bias is
  addressed analytically in the manuscript; the primary exposure GWAS
  (`ukb-b-5237`) and the FinnGen R11 outcome share no samples.
- FAERS queries use the openFDA public endpoint; query strings are embedded in
  the Python scripts and in Supplementary Table S5.

## Key outputs (`key_outputs/`)

Pre-computed results allowing table/figure-level verification without rerunning
the pipeline: `instruments.csv` (40 exposure instruments), per-dataset MR
(`results_*.csv`), MVMR (`mvmr_results_*.csv`), radial-MR outliers
(`radial_*.csv`), conditional F (`condf_*.csv`), PAGE instruments
(`instruments_wojcik_page.csv`), equivalence/power (`tost_power.json`), pooled
meta-analysis (`pooled_meta_nielsen_r11.json`), colocalization
(`coloc_results.csv`, `coloc_ahr_results.csv`, `coloc_prior_sensitivity.csv`),
PheWAS (`phewas_40snps.csv`, `phewas_per_snp_summary.csv`), and the BBJ
harmonised dataset (`bbj_harmonized.csv`, with corrected Wald-ratio weights).

Large raw inputs (FinnGen R11 file, BBJ zip, Nielsen full sumstats) are not
redistributed; download links are given above.
