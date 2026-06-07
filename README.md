# R-for-CP-factor-model-QMLE-method
This is the R.code for simulation (int/bdry) and real data analysis via QMLE of CP factor model time series

## SIMULATION PROCEDURE
### CLI Guide (Small Models for Checking Convergence-Rate Trends)

This repository provides a command-line interface (CLI) to run the simulation script:

- Script path: `scripts/run_sim.R`
- How to run: open a terminal at the repository root (same level as `scripts/`) and run `Rscript ...`


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
--T∈{50,100,200}
Suggested --d 5 (for bdry, make sure d = 5)
Suggested iteration budgets
int case: recommended --maxit 80
bdry case: recommended --maxit 150 to 200

### Example Commands (Copy & Run)
Rscript scripts/run_sim.R --case int --p 5 --q 5 --T 50  --d 4 --seeds 1,2,3,4,5 --maxit 80 --out out/int_p5_q5_T50.csv

Rscript scripts/run_sim.R --case bdry --p 5 --q 5 --T 50 --d 5 --seeds 1,2,3,4,5 --maxit 180 --out out/bdry_p5_q5_T50.csv

## REAL DATA ANALYSIS


