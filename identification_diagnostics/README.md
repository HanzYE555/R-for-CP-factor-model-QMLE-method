# Identification diagnostics for matrix-variate CP-factor models

Complete numerical results accompanying **Identification problem and
quasi-maximum likelihood estimation for matrix-variate CP-factor models**.

The paper's supplement presents means for **p=q=40** and **T=100,200,400**.
This directory retains the complete grid: all nine ordered pairs
p,q in {20,40,80}, T in {100,200,400,800}, and the interior, boundary
and non-identifiable ("exterior") designs. There are **10,796 available
estimates across 108 cells**. Only interior p=q=80,T=100,seeds1–4 are
missing (96 results); all other cells contain 100 seeds.

## Start here

| Location | Contents |
|---|---|
| `paper/p40_q40_T100_200_400_means.csv` | Nine mean-only rows selected for the supplement |
| `paper/p40_q40_T100_200_400_table.tex` | Copyable mean table, matching the supplement |
| `data/paper_40x40_T100_200_400_means.csv` | Extended mean extract, additionally including the absolute J margin |
| `data/geometry/tables/key_metrics_mean_var.csv` | All108 cells: means, sample variances, counts, contact/detection rates |
| `data/geometry/primary_by_seed.csv.gz` | All10,796 primary new-domain geometric results |
| `data/geometry/all_rho_by_seed/` | All50,577 rho candidates and their5°/10° cap searches, losslessly split/gzipped |
| `data/geometry/minima_pair_angles_by_seed/` | All552,507 independent candidate pairs, losslessly split/gzipped |
| `data/spectra/all_J_singular_values_mean_variance.csv.gz` | All3,456 dimension/design/T/ell/singular-index rows, includingorig/nor summaries |
| `data/fit_reports/J_svd_by_seed/` | All individual J singular values, seed-normalized ratios and Gram eigenvalues |
| `data/fit_reports/all_matrix_spectra_by_seed/` | Complete per-seed Ar,Br,Z,C and J spectra |
| `data/geometry/iteration_sensitivity/` | Saved-fit continuation comparisons re-evaluated on the new domain |
| `data/geometry/validation/` | Export/matrix checks, source checkpoint hashes and portable provenance |
| `report/supplementary_data.tex` | Optional standalone full table/figure report; not the main supplement |
| `code/` | Original simulation, geometric searches, analyses and validators |

CSV files larger than2MB are gzip-compressed; files larger than20MB are
split at row boundaries with repeated CSV headers. Parts remain in original
source order. `DATA_MANIFEST.json` records logical filenames, storage parts,
row counts, uncompressed SHA256 hashes and original source hashes. No stored
file exceeds10MB. The preparation script reconstructs byte-identical canonical
CSVs and checks every uncompressed SHA256.

## Definitions and experiment conventions

* New geometric searches use max_j|lambda_j|/||lambda||<=0.90, so lambda is at
  least25.8419° from every coefficient axis. The combined matrix has Frobenius
  norm1. Lambda itself is not also constrained to have Euclidean norm1.
* The physical cosine to each rank-one loading is an auxiliary value, not an
  additional exclusion constraint. It may exceed0.90.
* `rho_ind` uses each distinct2x2 minor once and equals half the ordered-minor
  ratio in the paper. E1 and orig/nor J singular values are unsquared; E1 is
  not multiplied by100. Each J's nor ratios are computed within each seed.
* `Delta` uses the nearest independently found E1 candidate; this need not be
  the independently lowest-E1 candidate. The5° joint score Q is infinite if
  no eligible independent pair is found. Counts accompany conditional means.
* Noise is generated with covariance I; fitting estimates diagonal noise
  variances in[0.05,20]. QMLE uses one exact and five near-truth starts,
  with per-start caps300/400/500 for interior/boundary/exterior. This is a
  truth-assisted numerical experiment. Later T fits also use an earlier fit.
* Separate multistart SLSQP searches use30 random starts, feasible truth-based
  and warm starts,400 evaluations/run and two cross-refinement rounds.
  The5°/10° cap searches use up to5 starts and300 evaluations/run, followed
  by400-evaluation refinement of the10° solution without the angular cap.
* Source fitting records and J spectra are unchanged. All primary geometry
  under `data/geometry/` was recomputed on the coefficient-angle domain.
  Files under `data/archive_matrix_domain/`, and geometric columns in the
  original `data/fit_reports/metrics_by_seed.csv.gz`, refer to the archived
  physical/loading-angle domain and must not be relabelled as new results.
* All valid saved fits are retained, irrespective of strict QMLE stopping flags.
  Sample variances use n-1. Shared seeds across dimensions are not900 independent
  replications. AUC files are retained as ancillary outputs, not core findings.

## Reproduce summaries, tables and figures without refitting

From this directory, install Python requirements and prepare an ignored workspace:

```sh
python -m pip install -r requirements.txt
python code/prepare_reproduction.py
cd .reproduce
python work/analyze_lambda_geometry.py --bootstrap-reps 1000
python work/build_lambda_manuscript_tables.py
python work/build_lambda_denominator_table.py
python work/plot_lambda_geometry.py
python work/analyze_lambda_iteration_sensitivity.py
python work/validate_lambda_geometry.py --aggregate
```

These commands analyze saved CSV exports and create output in `.reproduce/`;
they do not rerun QMLE. The1000-replicate bootstrap is ancillary. The original
analysis scripts are preserved byte-for-byte in `code/work/`. Their original
working paths are recreated by `prepare_reproduction.py` rather than changing
their mathematical implementations. The optional report can be compiled from
`report/` with `supplementary_data.tex` as the main file; figures are included.

## Recompute model fits or geometry

The complete original simulation is
`code/CP_QMLE_grid_20_40_80_100seeds.R`. It includes its required functions and
can be run from a chosen working directory with R4.6 and package `nloptr`:

```sh
Rscript code/CP_QMLE_grid_20_40_80_100seeds.R --cores=4 --out=PATH/CP_QMLE_grid_results
```

This is computationally expensive and is not required to read or reanalyze the
released results. The self-contained original simulation includes the archived
geometry; apply the new-domain runner to the fitted checkpoints afterward.
After preparing `.reproduce/`, run from that directory (R packages `nloptr`
and `jsonlite` are needed):

```sh
Rscript work/run_lambda_geometry_portable.R --source PATH/CP_QMLE_grid_results --output outputs/recomputed_lambda --pilot
```

Remove `--pilot` for all10,796 available keys; use `--shard i 12`
for i=1,...,12 to run separate processes. Do not run two
processes with the same shard. The source QMLE checkpoint tree and new per-case
RDS files are not tracked in Git; the portable source paths and MD5 hashes are
retained. The new-domain engine and optimization signature are unchanged.

## Provenance and coverage

Full export validation covers10,796 fits,50,577 rho centers and552,507 pairs.
The direct108-cell matrix audit covers4,169 endpoints. The continuation audit
covers28 common configurations plus three separately saved adaptive records.
All source checkpoint MD5 hashes were independently rechecked. Original local
absolute source paths were replaced with their portable checkpoint paths in the
released provenance CSV; the hashes of the fitted files were not changed.

`DATA_MANIFEST.json` is the logical-data manifest. `FILES_SHA256SUMS.txt` verifies
every payload file; `BUNDLE_AUDIT.json` summarizes completeness and sizes.
Large source checkpoint objects and generated time-series matrices are excluded;
all reported per-seed metrics, singular values and candidate-pair exports are
included. This repository publishes numerical search results, not a certificate
of global optimization or a proof of a convergence rate.
