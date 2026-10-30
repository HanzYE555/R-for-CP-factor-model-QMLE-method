# Complete results for all nine dimension pairs

**The repository contains all nine ordered pairs, not only 40 × 40.**
The `paper/` folder is the small selection used in the article. Use the links
below to browse the full mean results directly on GitHub.

[**All dimension pairs in one CSV: 108 rows**](all_dimensions_means.csv)

| p \ q | 20 | 40 | 80 |
|---|---|---|---|
| 20 | [20 × 20](p20_q20_means.csv) | [20 × 40](p20_q40_means.csv) | [20 × 80](p20_q80_means.csv) |
| 40 | [40 × 20](p40_q20_means.csv) | [40 × 40](p40_q40_means.csv) | [40 × 80](p40_q80_means.csv) |
| 80 | [80 × 20](p80_q20_means.csv) | [80 × 40](p80_q40_means.csv) | [80 × 80](p80_q80_means.csv) |

Each per-dimension CSV has **12 rows**: the interior, identifiable boundary and
non-identifiable (`exterior`) designs, each at **T = 100, 200, 400, 800**.
The CSVs contain the same ten primary mean metrics as the article's selected
mean table, with `p`, `q`, `regime`, `T` and `n` identifying each cell.
The values are copied from the complete summary CSV without recalculation or
additional rounding.

## Coverage and counts

There are **10,796 available fits across 108 cells**. Only interior
`p=80,q=80,T=100` has `n=96`, because seeds 1–4 are missing; every other cell
has 100 fits. `n` is the number of available fits. The mean joint score
`pair_joint_max_mean` uses finite Q values: exterior `p=80,q=20,T=200` has
99 finite scores, and all other cells have 100 (or 96 in the missing-seed cell).
The full summaries below retain the per-metric finite counts.

## Meaning of the primary mean columns

| Column | Quantity averaged across seeds |
|---|---|
| `rho_min_mean` | Lowest found independent-minor polynomial ratio |
| `E1_min_mean` | Independently lowest found relative rank-one error |
| `E1_at_rho_mean` | E1 at the selected lowest-rho center |
| `rho_at_E1_mean` | Polynomial ratio at the independently lowest-E1 candidate |
| `E1_cap10_mean` | Lowest found E1 within 10 degrees of the rho center |
| `pair_joint_max_mean` | Finite joint score Q for independent pairs within 5 degrees |
| `J_relative_min_mean` | Per-seed minimum normalized smallest J singular value |
| `center_to_E1_nearest_deg_090_mean` | Coefficient angle, in degrees, to the nearest independent E1 candidate |
| `E1_max_loading_cosine_mean` | Physical loading cosine at the independently lowest-E1 candidate |
| `cap10_max_loading_cosine_mean` | Physical loading cosine at the 10-degree cap solution |

E1 and J singular values are unsquared; E1 is not multiplied by 100. J ratios
are normalized within each seed before averaging. The polynomial ratio is
`rho_ind`, one half of the complete ordered-minor ratio in the paper. The
coefficient-axis cosine is constrained to at most 0.90; the physical loading
cosine is auxiliary and may exceed 0.90.

## Complete summaries, individual spectra and saved search results

* [All 108 cells, including means, variances, finite counts and detection/contact rates](../data/geometry/tables/key_metrics_mean_var.csv)
* [All 10,796 primary per-seed geometric results](../data/geometry/primary_by_seed.csv.gz)
* [All rho candidates and their cap searches](../data/geometry/all_rho_by_seed/)
* [All independent candidate-pair angles](../data/geometry/minima_pair_angles_by_seed/)
* [All J singular-value summaries, every ell and singular index](../data/spectra/all_J_singular_values_mean_variance.csv.gz)
* [All individual per-seed J spectra](../data/fit_reports/J_svd_by_seed/)
* [Standalone complete numerical report source](../report/supplementary_data.tex)
* [Reproduction instructions and code overview](../README.md)
