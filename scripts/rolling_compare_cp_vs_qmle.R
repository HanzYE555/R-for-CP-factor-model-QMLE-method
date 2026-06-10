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
R_t <- vector("list", length = 696)
for (t in 1:696) {
  R_t[[t]] <- 100 * matrix(E_capm[t,],nrow = 10)
}
R_arr_pqn <- simplify2array(R_t)          # p x q x Tn
R_npq <- aperm(R_arr_pqn, c(3, 1, 2))           # Tn x p x q  (n x p x q)

eig_real <- eigen(t(E_capm) %*% E_capm *(10000/696))$values

## ===== Parallel analysis for latent dimension (CAPM residuals) =====
## Assumes: E_capm is a numeric matrix with dim 696 x 100 (T x N)
stopifnot(is.matrix(E_capm), is.numeric(E_capm))
stopifnot(nrow(E_capm) == 696, ncol(E_capm) == 100)

Tn <- nrow(E_capm)
Nn <- ncol(E_capm)

# Helper: eigenvalues of covariance, optionally scaled to "% returns" units
eig_of <- function(X, scale_to_percent = TRUE) {
  # center columns (residuals should already be ~0 mean, but do it anyway)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- cov(Xc)  # N x N, uses 1/(T-1)
  if (scale_to_percent) S <- S * 10000  # (x100)^2
  eigen(S, symmetric = TRUE, only.values = TRUE)$values
}

# 1) Real eigenvalues
eig_real <- eig_of(E_capm, scale_to_percent = TRUE)

# 2) Null eigenvalues via column-wise time permutation
B <- 500  # increase to 1000+ if you want smoother thresholds
eig_null <- matrix(NA_real_, nrow = B, ncol = Nn)

set.seed(1)
for (b in 1:B) {
  Xb <- E_capm
  # independently permute time within each asset/portfolio column
  for (j in 1:Nn) {
    Xb[, j] <- sample(Xb[, j], size = Tn, replace = FALSE)
  }
  eig_null[b, ] <- eig_of(Xb, scale_to_percent = TRUE)
}

# 3) Thresholds and k-hat at different confidence levels
thr90 <- apply(eig_null, 2, quantile, probs = 0.90, na.rm = TRUE)
thr95 <- apply(eig_null, 2, quantile, probs = 0.95, na.rm = TRUE)
thr99 <- apply(eig_null, 2, quantile, probs = 0.99, na.rm = TRUE)

k_out <- c(
  k90 = sum(eig_real > thr90),
  k95 = sum(eig_real > thr95),
  k99 = sum(eig_real > thr99)
)
print(k_out)

# 4) Plot: real vs null threshold (95%)
op <- par(mar = c(4,4,2,1))
plot(eig_real, type = "b", pch = 16, cex = 0.7,
     xlab = "k (eigenvalue rank)", ylab = "eigenvalue (cov, % units)",
     main = "Parallel analysis: real vs null 95% threshold")
lines(thr95, type = "b", pch = 1, col = "red", cex = 0.7)
legend("topright", legend = c("real", "null 95%"),
       col = c("black", "red"), lty = 1, pch = c(16,1), bty = "n")
par(op)

# Optional: where do they cross?
cross_k <- which(eig_real <= thr95)[1]
cat("First k where eig_real <= thr95:", cross_k, "\n")


