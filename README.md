# R-for-CP-factor-model-QMLE-method
This is the R.code for simulation (int/bdry) and real data analysis via QMLE of CP factor model time series

## SIMULATION PROCEDURE
### CLI Guide (Small Models for Checking Convergence-Rate Trends)

This repository provides a command-line interface (CLI) to run the simulation script:

- Script path: `scripts/run_sim.R`



### Arguments
--case: scenario
int: interior case
bdry: boundary case
--p: dimension p (integer, required)
--q: dimension q (integer, required)
--T: sample size T (integer, required)
--d: factor dimension (default: 4; note that bdry requires d = 5)
--seeds: comma-separated random seeds (e.g., 1,2,3,4,5)
--maxit: maximum number of iterations
--out: output CSV path (directories will be created automatically if needed)
The output is a CSV file containing results for each seed plus one aggregated (average) row (seed = NA).

### Recommended “Small-Model” Settings (for Quick Trend Checks)
To quickly verify that error metrics decrease as T increases, we recommend using small dimensions: 
--p,q∈{5,10,20}
--T∈{50,100,200}. 
Suggested --d 5 (for bdry, make sure d = 5)\\
Suggested iteration budgets: 
int case recommended --maxit 80;

bdry case recommended --maxit 150 to 200

### Example Commands (Copy & Run)
```bash
Rscript scripts/run_sim.R --case int --p 5 --q 5 --T 50  --d 4 --seeds 1,2,3,4,5 --maxit 80 --out out/int_p5_q5_T50.csv

Rscript scripts/run_sim.R --case bdry --p 5 --q 5 --T 50 --d 5 --seeds 1,2,3,4,5 --maxit 180 --out out/bdry_p5_q5_T50.csv
```
## REAL DATA ANALYSIS
### Data (Fama–French, monthly)

This project uses two monthly datasets from the Kenneth R. French Data Library:

1) **Fama/French 3 Factors (monthly)**  
   File (ZIP): `F-F_Research_Data_Factors_CSV.zip`

2) **100 Portfolios Formed on Size and Book-to-Market (10×10, monthly)**  
   File (ZIP): `100_Portfolios_10x10_CSV.zip`

#### How to obtain the data
Download the ZIP files manually from the official data library (to avoid redistributing the raw data in this repository), then place them in:

- `data/raw/F-F_Research_Data_Factors_CSV.zip`
- `data/raw/100_Portfolios_10x10_CSV.zip`

#### Build the processed dataset
From the repository root, run:

```bash  
Rscript scripts/prepare_data.R
```

This creates:
- `data/processed/ten_by_ten_prepared.rds`

The RDS contains:

- date: monthly dates
- R_mat: 10×10 portfolio returns matrix (T×100), in decimal returns (e.g., 0.01 = 1%)
- MktRF: market excess return (decimal)
- RF: risk-free rate (decimal)
- E_capm: CAPM residuals of portfolio excess returns on MktRF (T×100, decimal)
- colnames_10x10: portfolio column names






## Reproducibility: Parallel Analysis + Rolling (CP.Unified vs QMLE)

This repo contains two scripts for the experiments in the same paper section:

1. **Parallel analysis** (real vs null, with 95% threshold).
2. **Rolling comparison** (CP.Unified vs QMLE), including rolling-window control and QMLE/VAR settings.

Both workflows require the data-prep step first.

---

## 1) Parallel analysis (real vs null, 95%)

### Step 0: prepare data (required)
```bash 
Rscript prepare_data.R
```
### Step 1: run parallel analysis
```bash
Rscript parallel_analysis_entry.R
```

#### Output 
This produces a figure similar to:
- `parallel analysis real vs null 95%` (saved under `outputs/parallel`)

## 2) Rolling: CP.Unified vs QMLE (with rolling eigengap-ready outputs)

### Step 0: prepare data (required)
```bash
Rscript prepare_data.R
```
### Step 1: run rolling entry
```bash
Rscript rolling_entry.R
```

#### Parameters (command line)
`rolling_entry.R` reads the following arguments (defaults shown are the script defaults):
```r
d       <- as_int(get_arg("d", "2"))      # latent dimension, suggested latent dimension: d\in [2,6]
train_T <- as_int(get_arg("train_T", "456"))
horizon <- c(1, 2)

# rolling window control
S_run <- get_arg("S_run", NULL)           # NULL => run ALL windows; otherwise run first S_run windows

# VAR + stability
pmax_var      <- as_int(get_arg("pmax_var", "24"))
var_ic        <- get_arg("var_ic", "AIC") # "AIC" or "BIC"
stable_thresh <- as_num(get_arg("stable_thresh", "0.999"))

# QMLE
qmle_maxit <- as_int(get_arg("qmle_maxit", "120"))

# output + verbosity
verbose_every <- as_int(get_arg("verbose_every", "20"))
out_dir <- get_arg("out_dir", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(d %in% 1:20)
stopifnot(var_ic %in% c("AIC", "BIC"))
```
you can use the following R code example to get results: 
```r
# assume you already ran prepare_data.R and loaded/created R_npq in the session

# all windows
res_rec_2 <- run_one(R_npq, d = 2, S_run = NULL)

# first 40 windows
res_rec_5_40 <- run_one(R_npq, d = 5, S_run = 40)

# eigengap plot for d=5, first 40 windows (QMLE A/B)
out_dir <- "outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "eigengap_logratio_AB_d5_S40.png"), width=1200, height=900, res=150)
par(mfrow=c(2,1), mar=c(4,4,3,1))
plot_eigengap_logratio(res_rec_5_40$store$sv_A_q,
                       main = "QMLE A: log(s_k/s_{k+1}) (d=5, S_run=40)")
plot_eigengap_logratio(res_rec_5_40$store$sv_B_q,
                       main = "QMLE B: log(s_k/s_{k+1}) (d=5, S_run=40)")
dev.off()
par(mfrow=c(1,1))
```

We also provide a choice for bash: 
```bash
Rscript rolling_entry.R --d=5 --S_run=40 --out_dir=outputs
```


















