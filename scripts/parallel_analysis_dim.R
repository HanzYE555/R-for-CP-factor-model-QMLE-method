# rolling_compare_cp_vs_qmle.R
# ============================================================
# Rolling-window forecast comparison:
#   CP.Unified vs QMLE (CP factor model)
# With option roll_n: only run first K rolling windows
#
# Inputs expected in your repo:
#   - E_capm: [T x N] numeric matrix (CAPM residuals)
#   - R_npq:  [T x p x q] numeric array (your constructed matrix time series)
#
# Requires user-provided model functions (source them before running):
#   - CP_MTS(Y_array, d, ...) -> list(A, B, Fhat, ...) (or adapt in wrapper below)
#   - QMLE_marm(Y_array, d, init = NULL, maxit = ...) -> list(A, B, ...)
#   - odotprod(A, B, p, q, d): factor extraction mapping
#   - vec(M), invvec(v, p, q): vectorize helpers
#
# ============================================================
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
})

# -------------------------
# Config
# -------------------------
in_rds  <- file.path("data", "processed", "ten_by_ten_prepared.rds")
out_dir <- file.path("outputs", "parallel")
out_rds <- file.path(out_dir, "parallel_analysis_results.rds")
out_png <- file.path(out_dir, "parallel_analysis_real_vs_null95.png")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load
# -------------------------
dat <- readRDS(in_rds)
E_capm <- dat$E_capm

stopifnot(is.matrix(E_capm), is.numeric(E_capm))
Tn <- nrow(E_capm)
Nn <- ncol(E_capm)

cat(sprintf("Loaded E_capm: T=%d, N=%d\n", Tn, Nn))

# If you really want to enforce 10x10:
if (Nn != 100) stop("Expected N=100 portfolios but got N=", Nn, call. = FALSE)

# -------------------------
# Helper
# -------------------------
eig_of <- function(X, scale_to_percent = TRUE) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- cov(Xc)  # uses 1/(T-1)
  if (scale_to_percent) S <- S * 10000
  eigen(S, symmetric = TRUE, only.values = TRUE)$values
}

# 1) Real eigenvalues
eig_real <- eig_of(E_capm, scale_to_percent = TRUE)

# 2) Null eigenvalues via column-wise time permutation
B <- 500
eig_null <- matrix(NA_real_, nrow = B, ncol = Nn)

set.seed(1)
for (b in seq_len(B)) {
  Xb <- E_capm
  for (j in seq_len(Nn)) Xb[, j] <- sample(Xb[, j], size = Tn, replace = FALSE)
  eig_null[b, ] <- eig_of(Xb, scale_to_percent = TRUE)
  if (b %% 50 == 0) cat(sprintf("perm %d / %d\n", b, B))
}

# 3) Thresholds and k-hat
thr90 <- apply(eig_null, 2, quantile, probs = 0.90, na.rm = TRUE)
thr95 <- apply(eig_null, 2, quantile, probs = 0.95, na.rm = TRUE)
thr99 <- apply(eig_null, 2, quantile, probs = 0.99, na.rm = TRUE)

k_out <- c(
  k90 = sum(eig_real > thr90),
  k95 = sum(eig_real > thr95),
  k99 = sum(eig_real > thr99)
)

print(k_out)

cross_k <- which(eig_real <= thr95)[1]
cat("First k where eig_real <= thr95:", cross_k, "\n")

# 4) Save results
saveRDS(list(
  eig_real = eig_real,
  eig_null = eig_null,
  thr90 = thr90,
  thr95 = thr95,
  thr99 = thr99,
  k_out = k_out,
  cross_k_95 = cross_k,
  B = B,
  in_rds = in_rds
), out_rds)
cat("Saved: ", out_rds, "\n", sep = "")

# -------------------------
# Plot (works under Rscript by saving to PNG)
# -------------------------
png(out_png, width = 1200, height = 800, res = 150)
op <- par(mar = c(4,4,2,1))
plot(eig_real, type = "b", pch = 16, cex = 0.7,
     xlab = "k (eigenvalue rank)", ylab = "eigenvalue (cov, % units)",
     main = "Parallel analysis: real vs null 95% threshold")
lines(thr95, type = "b", pch = 1, col = "red", cex = 0.7)
legend("topright", legend = c("real", "null 95%"),
       col = c("black", "red"), lty = 1, pch = c(16,1), bty = "n")
par(op)
dev.off()

cat("Saved plot: ", out_png, "\n", sep = "")

