#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDTSA)
  library(MASS)
  library(MTS)
})

# ============================================================
# 0) Minimal command-line args
# ============================================================
get_arg <- function(key, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}
as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}
as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

# ============================================================
# 1) User-configurable (reviewer can override by CLI)
# ============================================================
data_path <- get_arg("data", file.path("data", "processed", "ten_by_ten_prepared.rds"))

d         <- as_int(get_arg("d", "4"))              # latent dimension
train_T   <- as_int(get_arg("train_T", "456"))
horizon   <- c(1, 2)

# rolling window control
s_start   <- as_int(get_arg("s_start", "1"))
n_windows <- get_arg("n_windows", NULL)
n_windows <- if (is.null(n_windows) || n_windows == "") NULL else as_int(n_windows)

# VAR + stability
pmax_var      <- as_int(get_arg("pmax_var", "24"))
var_ic        <- get_arg("var_ic", "AIC")           # "AIC" or "BIC"
stable_thresh <- as_num(get_arg("stable_thresh", "0.999"))

# QMLE
qmle_maxit <- as_int(get_arg("qmle_maxit", "120"))

# output + verbosity
verbose_every <- as_int(get_arg("verbose_every", "20"))
out_dir <- get_arg("out_dir", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(d %in% 1:20)
stopifnot(var_ic %in% c("AIC", "BIC"))

cat("============================================================\n")
cat("Rolling CP.Unified vs QMLE\n")
cat(sprintf("data: %s\n", data_path))
cat(sprintf("d=%d, train_T=%d, horizon={%s}\n", d, train_T, paste(horizon, collapse = ",")))
cat(sprintf("s_start=%d, n_windows=%s\n", s_start, if (is.null(n_windows)) "ALL" else n_windows))
cat(sprintf("VAR: pmax=%d, IC=%s, stable_thresh=%.4f\n", pmax_var, var_ic, stable_thresh))
cat(sprintf("QMLE: maxit=%d\n", qmle_maxit))
cat("============================================================\n")

# ============================================================
safe_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-12) return(1.0)
  s
}


svd_report <- function(M, name = "M", tol = NULL, k = NULL) {
  s <- svd(M, nu = 0, nv = 0)$d
  if (is.null(tol)) tol <- max(dim(M)) * .Machine$double.eps * max(s)
  r_num <- sum(s > tol)
  
  cat("\n--- SVD report:", name, "---\n")
  cat("dim =", paste(dim(M), collapse = " x "), "\n")
  cat("tol =", format(tol, digits = 6), "\n")
  cat("numerical rank =", r_num, "\n")
  cat("singular values (descending):\n")
  if (!is.null(k)) s <- head(s, k)
  print(format(s, digits = 6), quote = FALSE)
  
  invisible(list(singular_values = s, tol = tol, numerical_rank = r_num))
}

col_normalize_to1 <- function(M, eps = 1e-12) {
  cn <- sqrt(colSums(M^2))
  cn[cn < eps] <- 1
  sweep(M, 2, cn, "/")
}

reconstruct_Y_from_ABf <- function(A, B, f_vec) {
  A %*% diag(as.numeric(f_vec), nrow = length(f_vec)) %*% t(B)
}

## -------------------------
## Khatri-Rao loading
## -------------------------
odotprod <- function(A, B, p, q, d) {
  K <- matrix(0, nrow = p * q, ncol = d)
  for (c in 1:d) K[, c] <- kronecker(B[, c], A[, c])
  K
}
loading_KR <- function(A, B) {
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  odotprod(A, B, p, q, d)
}
loading_diagnostics <- function(A, B, name = "", print = TRUE) {
  L <- loading_KR(A, B)
  sv <- svd(L, nu = 0, nv = 0)$d
  cond <- if (min(sv) > 0) max(sv) / min(sv) else Inf
  rnk <- qr(L)$rank
  out <- list(singular_values = sv, condition_number = cond, rank = rnk, dim = dim(L))
  
  if (print) {
    cat("\n--- Loading diagnostics", if (nzchar(name)) paste0(" [", name, "]") else "", " ---\n", sep = "")
    cat("dim(L) =", nrow(L), "x", ncol(L), "\n")
    cat("rank(L) =", rnk, "\n")
    cat("singular values =", paste(format(sv, digits = 6), collapse = ", "), "\n")
    cat("condition number =", format(cond, digits = 6), "\n")
  }
  out
}

## -------------------------
## CPMTS wrapper
## -------------------------
fit_cpmts_unified_safe <- function(R_npq, d, lag.k = 20, lag.ktilde = 10) {
  Tn <- dim(R_npq)[1]; p <- dim(R_npq)[2]; q <- dim(R_npq)[3]
  delta1 <- 2 * sqrt(log(p * q) / Tn)
  
  res <- CP_MTS(
    R_npq,
    xi = NULL,
    Rank = list(d = d, d1 = d, d2 = d),
    lag.k = lag.k,
    lag.ktilde = lag.ktilde,
    method = "CP.Unified",
    thresh1 = TRUE, thresh2 = TRUE, thresh3 = TRUE,
    delta1 = delta1, delta2 = delta1, delta3 = delta1
  )
  
  
  
  A <- col_normalize_to1(res$A)
  B <- col_normalize_to1(res$B)
  
  f <- NULL
  if (!is.null(res$f)) f <- res$f
  if (is.null(f) && !is.null(res$F)) f <- res$F
  if (is.null(f) && !is.null(res$fhat)) f <- res$fhat
  if (is.null(f) && !is.null(res$factor)) f <- res$factor
  if (is.null(f)) stop("CP_MTS result does not contain latent factors field (f/F/fhat/factor not found).")
  
  if (is.list(f)) fmat <- do.call(rbind, lapply(f, function(x) as.numeric(x))) else fmat <- as.matrix(f)
  if (ncol(fmat) != d) {
    if (nrow(fmat) == d && ncol(fmat) == Tn) fmat <- t(fmat)
  }
  stopifnot(nrow(fmat) == Tn, ncol(fmat) == d)
  
  list(res = res, A = A, B = B, f = fmat)
}

## -------------------------
## Data helpers
## -------------------------
as_Ylist_from_Rnpq <- function(R_npq) {
  Tn <- dim(R_npq)[1]; p <- dim(R_npq)[2]; q <- dim(R_npq)[3]
  Ylist <- vector("list", Tn)
  for (t in 1:Tn) Ylist[[t]] <- matrix(R_npq[t,,], nrow = p, ncol = q)
  Ylist
}

empirical_second_moment <- function(Ylist, p, q, Tn) {
  M <- matrix(0, nrow = p * q, ncol = p * q)
  for (t in 1:Tn) {
    v <- as.vector(Ylist[[t]])
    M <- M + tcrossprod(v)
  }
  M / Tn
}

## -------------------------
## QMLE with SPD Mx via Cholesky parametrization
## -------------------------
pack_lower_raw <- function(Lraw) {
  d <- nrow(Lraw)
  v <- numeric(d * (d + 1) / 2)
  idx <- 1
  for (i in 1:d) {
    for (j in 1:i) {
      v[idx] <- Lraw[i, j]
      idx <- idx + 1
    }
  }
  v
}

unpack_lower_raw <- function(v, d) {
  Lraw <- matrix(0, d, d)
  idx <- 1
  for (i in 1:d) {
    for (j in 1:i) {
      Lraw[i, j] <- v[idx]
      idx <- idx + 1
    }
  }
  Lraw
}

build_L_from_Lraw <- function(Lraw) {
  L <- Lraw
  diag(L) <- exp(diag(Lraw))
  L
}

pack_params <- function(A, B, Lraw_lower_vec, log_diag_Se) {
  c(as.vector(A), as.vector(B), as.vector(Lraw_lower_vec), as.vector(log_diag_Se))
}

unpack_params <- function(theta, p, q, d) {
  nA <- p * d
  nB <- q * d
  nL <- d * (d + 1) / 2
  nSe <- p * q
  
  A <- matrix(theta[1:nA], nrow = p, ncol = d)
  B <- matrix(theta[(nA + 1):(nA + nB)], nrow = q, ncol = d)
  
  off1 <- nA + nB
  Lraw_vec <- theta[(off1 + 1):(off1 + nL)]
  Lraw <- unpack_lower_raw(Lraw_vec, d)
  
  off2 <- off1 + nL
  log_diag_Se <- theta[(off2 + 1):(off2 + nSe)]
  
  list(A = A, B = B, Lraw = Lraw, log_diag_Se = log_diag_Se)
}

q_function <- function(M_y, A, B, Mx, log_diag_Se, ridge = 1e-8) {
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  N <- p * q
  
  K <- odotprod(A, B, p, q, d)
  
  Se_diag <- exp(log_diag_Se)
  Se_eff <- Se_diag + ridge
  Dinv <- 1 / Se_eff
  logdet_D <- sum(log(Se_eff))
  
  # G = K' D^{-1} K
  W <- K * Dinv
  G <- crossprod(K, W)
  
  # Mx SPD by construction; ridge for numerical stability
  Mx <- Mx + ridge * diag(d)
  chol_Mx <- chol(Mx)
  Mx_inv <- chol2inv(chol_Mx)
  
  S <- Mx_inv + G
  chol_S <- try(chol(S), silent = TRUE)
  if (inherits(chol_S, "try-error")) {
    S <- S + 1e-8 * diag(d)
    chol_S <- chol(S)
  }
  
  logdet_S <- 2 * sum(log(diag(chol_S)))
  logdet_Mx <- 2 * sum(log(diag(chol_Mx)))
  logdet_SigmaY <- logdet_D + (logdet_S + logdet_Mx)
  
  # trace term using Woodbury identity:
  # tr(SigmaY^{-1} My) = tr(D^{-1} My) - tr(S^{-1} K' D^{-1} My D^{-1} K)
  tr_Dinv_My <- sum(Dinv * diag(M_y))
  
  # More stable: compute T = K' D^{-1} My D^{-1} K without forming D^{-1}MyD^{-1} explicitly
  # Let Z = (My %*% (K * Dinv)) ; then T = (K * Dinv)' %*% Z
  Z <- M_y %*% (K * Dinv)              # (pq) x d
  Tmat <- crossprod(K * Dinv, Z)       # d x d
  
  SinvT <- backsolve(chol_S, forwardsolve(t(chol_S), Tmat))
  tr_term <- tr_Dinv_My - sum(diag(SinvT))
  
  (logdet_SigmaY + tr_term) / N
}

eig_report <- function(M, name = "M", tol = 1e-8) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  ev <- sort(as.numeric(ev), decreasing = TRUE)
  
  out <- list(
    name = name,
    eigenvalues = ev,
    min = min(ev),
    max = max(ev),
    cond = max(ev) / min(ev),
    n_small = sum(ev < tol),
    tol = tol
  )
  
  cat("\n--- Eigen report:", name, "---\n")
  cat("max eig =", format(out$max, digits = 8), "\n")
  cat("min eig =", format(out$min, digits = 8), "\n")
  cat("condition number (max/min) =", format(out$cond, digits = 8), "\n")
  cat("count < tol (", tol, ") =", out$n_small, "\n", sep = "")
  cat("eigenvalues (desc):", paste(format(ev, digits = 8), collapse = ", "), "\n")
  
  invisible(out)
}

objective_fn <- function(theta, M_y, p, q, d) {
  pars <- unpack_params(theta, p, q, d)
  A <- col_normalize_to1(pars$A)
  B <- col_normalize_to1(pars$B)
  L <- build_L_from_Lraw(pars$Lraw)
  Mx <- L %*% t(L)
  q_function(M_y, A, B, Mx = Mx, log_diag_Se = pars$log_diag_Se, ridge = 1e-8)
}

fit_qmle <- function(M_y, p, q, d, init, maxit = 120) {
  Ainit <- col_normalize_to1(init$A)
  Binit <- col_normalize_to1(init$B)
  
  Mx0 <- init$M_x
  Mx0 <- Mx0 + 1e-6 * diag(d)
  L0 <- t(chol(Mx0))
  Lraw0 <- L0
  diag(Lraw0) <- log(pmax(diag(L0), 1e-12))
  Lraw_vec0 <- pack_lower_raw(Lraw0)
  
  logSe0 <- init$log_diag_Se
  
  theta0 <- pack_params(Ainit, Binit, Lraw_vec0, logSe0)
  
  width <- 5
  lower <- theta0 - width
  upper <- theta0 + width
  
  opt <- optim(
    par = theta0,
    fn = function(th) objective_fn(th, M_y, p, q, d),
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = maxit, factr = 1e7)
  )
  
  pars <- unpack_params(opt$par, p, q, d)
  Ahat <- col_normalize_to1(pars$A)
  Bhat <- col_normalize_to1(pars$B)
  Lhat <- build_L_from_Lraw(pars$Lraw)
  Mxhat <- Lhat %*% t(Lhat)
  
  list(
    Ahat = Ahat,
    Bhat = Bhat,
    Mxhat = Mxhat,
    logSehat = pars$log_diag_Se,
    opt = opt,
    conv_code = opt$convergence,
    conv_msg  = opt$message
  )
}

## -------------------------
## QMLE factor extraction: GLS/BLUP under working covariance
## f_hat(t) = (Mx^{-1} + K' D^{-1} K)^{-1} K' D^{-1} y_t
## -------------------------
extract_f_blup <- function(Ylist, A, B, Mx, log_diag_Se, ridge = 1e-8) {
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  K <- loading_KR(A, B)               # (pq) x d
  
  Se_diag <- exp(log_diag_Se)
  Se_eff <- Se_diag + ridge
  Dinv <- 1 / Se_eff
  
  # Precompute S = Mx^{-1} + K' D^{-1} K
  Mx <- Mx + ridge * diag(d)
  chol_Mx <- chol(Mx)
  Mx_inv <- chol2inv(chol_Mx)
  
  KD <- K * Dinv                      # (pq) x d == D^{-1}K rowwise
  G <- crossprod(K, KD)               # d x d
  S <- Mx_inv + G
  chol_S <- try(chol(S), silent = TRUE)
  if (inherits(chol_S, "try-error")) {
    S <- S + 1e-8 * diag(d)
    chol_S <- chol(S)
  }
  
  Fhat <- matrix(NA_real_, nrow = length(Ylist), ncol = d)
  for (t in 1:length(Ylist)) {
    y <- as.vector(Ylist[[t]])
    rhs <- crossprod(K, Dinv * y)     # K' D^{-1} y
    fhat <- backsolve(chol_S, forwardsolve(t(chol_S), rhs))
    Fhat[t, ] <- as.numeric(fhat)
  }
  Fhat
}

## -------------------------
## VAR forecasting (iterated h-step, no peeking) with stability constraint
## -------------------------
var_pred_step <- function(state_mat, Bhat, include_const = TRUE) {
  p <- nrow(state_mat)
  x <- NULL
  for (lag in 1:p) x <- c(x, state_mat[p - lag + 1, ])
  if (include_const) x <- c(1, x)
  as.numeric(x %*% Bhat)
}

companion_spectral_radius <- function(Bhat, d, p, include_const = TRUE) {
  if (include_const) Astack <- t(Bhat[-1, , drop = FALSE]) else Astack <- t(Bhat)
  # Astack: d x (d*p)
  if (p == 1) {
    Fmat <- Astack
  } else {
    top <- Astack
    bottom <- cbind(diag(d * (p - 1)), matrix(0, nrow = d * (p - 1), ncol = d))
    Fmat <- rbind(top, bottom)
  }
  max(Mod(eigen(Fmat, only.values = TRUE)$values))
}

fit_var_ols_ic_stable <- function(F, pmax = 12, include_const = TRUE,
                                  ic = c("BIC","AIC"),
                                  stable_thresh = 0.999) {
  ic <- match.arg(ic)
  n <- nrow(F); d <- ncol(F)
  best <- list(ic = Inf, p = NA_integer_, Bhat = NULL, rho = NA_real_)
  
  for (p in 1:pmax) {
    if (n <= p + 5) next
    
    Y <- F[(p + 1):n, , drop = FALSE]
    X <- NULL
    for (lag in 1:p) X <- cbind(X, F[(p + 1 - lag):(n - lag), , drop = FALSE])
    if (include_const) X <- cbind(1, X)
    
    XtX <- crossprod(X)
    if (qr(XtX)$rank < ncol(XtX)) next
    
    Bhat <- solve(XtX, crossprod(X, Y))
    E <- Y - X %*% Bhat
    Sigma <- crossprod(E) / nrow(E)
    
    ld <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
    if (!is.finite(ld)) next
    
    n_eff <- nrow(E)
    k_params <- d * (d * p + if (include_const) 1 else 0)
    ic_val <- if (ic == "AIC") (ld + 2 * k_params / n_eff) else (ld + log(n_eff) * k_params / n_eff)
    
    rho <- companion_spectral_radius(Bhat, d = d, p = p, include_const = include_const)
    if (!is.finite(rho) || rho >= stable_thresh) next
    
    if (ic_val < best$ic) best <- list(ic = ic_val, p = p, Bhat = Bhat, rho = rho)
  }
  
  if (!is.finite(best$ic)) stop("Stable VAR selection failed: no stable model found under pmax.")
  best
}

generate_fhat_h_seq_stable <- function(f_train, Tn, train_T, h,
                                       pmax = 12, include_const = TRUE,
                                       ic = "BIC", stable_thresh = 0.999) {
  # Standardize factors for VAR fitting (reduces near-unit-root and scaling issues)
  mu <- colMeans(f_train)
  sdv <- apply(f_train, 2, sd)
  sdv[sdv < 1e-12] <- 1
  Fz <- scale(f_train, center = mu, scale = sdv)
  
  fit <- fit_var_ols_ic_stable(Fz, pmax = pmax, include_const = include_const,
                               ic = ic, stable_thresh = stable_thresh)
  p <- fit$p; Bhat <- fit$Bhat
  
  state <- Fz[(train_T - p + 1):train_T, , drop = FALSE]
  n_out <- (Tn - h) - train_T + 1
  d <- ncol(Fz)
  out_z <- matrix(NA_real_, nrow = n_out, ncol = d)
  
  state_roll <- state
  for (k in 1:n_out) {
    tmp <- state_roll
    f_new <- NULL
    for (step in 1:h) {
      f_new <- var_pred_step(tmp, Bhat, include_const)
      tmp <- rbind(tmp[-1, , drop = FALSE], matrix(f_new, nrow = 1))
    }
    out_z[k, ] <- f_new
    
    f1 <- var_pred_step(state_roll, Bhat, include_const)
    state_roll <- rbind(state_roll[-1, , drop = FALSE], matrix(f1, nrow = 1))
  }
  
  # Inverse-transform forecasts
  out <- sweep(out_z, 2, sdv, "*")
  out <- sweep(out, 2, mu, "+")
  
  list(fhat_h = out, p_selected = p, rho = fit$rho, ic = fit$ic)
}

## -------------------------
## Evaluation: pointwise rRMSE / rMAE in Y-space
## -------------------------
eval_pointwise_rmetrics <- function(R_npq, A, B, fhat_h, train_T, h) {
  Tn <- dim(R_npq)[1]; p <- dim(R_npq)[2]; q <- dim(R_npq)[3]
  tau_vec <- (train_T + h):Tn
  stopifnot(nrow(fhat_h) == length(tau_vec))
  
  rRMSE <- numeric(length(tau_vec))
  rMAE  <- numeric(length(tau_vec))
  
  for (k in seq_along(tau_vec)) {
    tau <- tau_vec[k]
    Y_true <- matrix(R_npq[tau, , ], nrow = p, ncol = q)
    Y_hat  <- A %*% diag(as.numeric(fhat_h[k, ]), nrow = ncol(A)) %*% t(B)
    
    E <- Y_hat - Y_true
    rRMSE[k] <- sqrt(mean(E^2))
    rMAE[k]  <- mean(abs(E))
  }
  
  list(
    tau = tau_vec,
    rRMSE = rRMSE,
    rMAE = rMAE,
    mean_rRMSE = mean(rRMSE),
    sd_rRMSE   = sd(rRMSE),
    mean_rMAE  = mean(rMAE),
    sd_rMAE    = sd(rMAE)
  )
}

## -------------------------
## Main: compare in Y-space RMSE
## -------------------------
compare_cpmts_vs_qmle_rolling_rrmse <- function(
    R_npq,
    d = 4,
    train_T = 456,              # rolling window length
    horizon = c(1, 2),
    pmax_var = 24,
    include_const = TRUE,
    var_ic = "BIC",
    stable_thresh = 0.999,
    qmle_maxit = 80,
    verbose_every = 20,
    use_warm_start = TRUE,
    store_cp = FALSE,           # whether also store CP A,B and singular values
    eps_denom = 1e-12,          # avoid divide-by-zero in denominators
    store_abs_errors = TRUE     # store absolute RMSE/MAE and some scale diagnostics
) {
  Tn <- dim(R_npq)[1]; p <- dim(R_npq)[2]; q <- dim(R_npq)[3]
  stopifnot(train_T >= 2, train_T < Tn)
  
  # Rolling: estimation window is [s, t_end] with fixed length train_T
  S <- Tn - train_T
  
  out <- list(
    d = d, train_T = train_T, T = Tn, S = S,
    horizon = horizon,
    rolling = TRUE,
    r_definition = "SD-normalized: rRMSE = sqrt(mean((E/sd_ij)^2)); rMAE = mean(|E/sd_ij|), sd_ij from in-window elementwise sd",
    per_h = setNames(vector("list", length(horizon)), paste0("h", horizon))
  )
  
  for (h in horizon) {
    n_eval <- S - h + 1
    out$per_h[[paste0("h", h)]] <- list(
      rRMSE_cp = rep(NA_real_, n_eval),
      rMAE_cp  = rep(NA_real_, n_eval),
      rRMSE_q  = rep(NA_real_, n_eval),
      rMAE_q   = rep(NA_real_, n_eval),
      var_p_cp = rep(NA_integer_, n_eval),
      var_p_q  = rep(NA_integer_, n_eval)
    )
    if (isTRUE(store_abs_errors)) {
      out$per_h[[paste0("h", h)]]$RMSE_cp <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$MAE_cp  <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$RMSE_q  <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$MAE_q   <- rep(NA_real_, n_eval)
      
      # diagnostics about sd scale in the rolling window
      out$per_h[[paste0("h", h)]]$sd_mean <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$sd_med  <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$sd_min  <- rep(NA_real_, n_eval)
      out$per_h[[paste0("h", h)]]$sd_p10  <- rep(NA_real_, n_eval)
    }
  }
  
  # storage for per-origin estimates
  out$store <- list(
    s = rep(NA_integer_, S),
    t_end = rep(NA_integer_, S),
    
    A_q = array(NA_real_, dim = c(p, d, S)),
    B_q = array(NA_real_, dim = c(q, d, S)),
    sv_A_q = matrix(NA_real_, nrow = S, ncol = min(p, d)),
    sv_B_q = matrix(NA_real_, nrow = S, ncol = min(q, d))
  )
  
  if (isTRUE(store_cp)) {
    out$store$A_cp    <- array(NA_real_, dim = c(p, d, S))
    out$store$B_cp    <- array(NA_real_, dim = c(q, d, S))
    out$store$sv_A_cp <- matrix(NA_real_, nrow = S, ncol = min(p, d))
    out$store$sv_B_cp <- matrix(NA_real_, nrow = S, ncol = min(q, d))
  }
  
  # warm-start holder for QMLE across rolling windows
  last_init <- NULL
  
  # ---- helpers ----
  .RMSE_abs <- function(E) sqrt(mean(E^2))
  .MAE_abs  <- function(E) mean(abs(E))
  
  # elementwise sd from rolling window (train window)
  .sd_ij_from_window <- function(R_window, eps = 1e-12) {
    # R_window: [train_T, p, q]
    sd_ij <- apply(R_window, c(2, 3), sd)
    pmax(sd_ij, eps)
  }
  
  # sd-normalized relative metrics
  .rRMSE_sd <- function(E, sd_ij) sqrt(mean((E / sd_ij)^2))
  .rMAE_sd  <- function(E, sd_ij) mean(abs(E / sd_ij))
  # -----------------
  
  for (s in 1:S) {
    t_end <- s + train_T - 1
    if (t_end + 1 > Tn) break
    
    R_window <- R_npq[s:t_end, , , drop = FALSE]  # [train_T, p, q]
    
    # precompute sd_ij once per rolling window (same for all horizons h)
    sd_ij <- .sd_ij_from_window(R_window, eps = eps_denom)
    
    # --- Fit CP_MTS on rolling window
    cp_tr <- fit_cpmts_unified_safe(R_window, d = d)
    A_cp <- cp_tr$A; B_cp <- cp_tr$B; f_cp <- cp_tr$f
    
    # --- Fit QMLE on rolling window
    Ylist <- as_Ylist_from_Rnpq(R_window)
    M_y <- empirical_second_moment(Ylist, p, q, train_T)
    
    if (use_warm_start && !is.null(last_init)) {
      init_q <- last_init
    } else {
      init_q <- list(
        A = A_cp,
        B = B_cp,
        M_x = diag(d),
        log_diag_Se = rep(log(1), p * q)
      )
    }
    
    qmle <- fit_qmle(
      M_y = M_y, p = p, q = q, d = d,
      init = init_q, maxit = qmle_maxit
    )
    
    A_q <- qmle$Ahat
    B_q <- qmle$Bhat
    
    f_q <- extract_f_blup(
      Ylist, A_q, B_q,
      Mx = qmle$Mxhat,
      log_diag_Se = qmle$logSehat
    )
    
    if (use_warm_start) {
      last_init <- list(
        A = A_q,
        B = B_q,
        M_x = qmle$Mxhat,
        log_diag_Se = qmle$logSehat
      )
    }
    
    # store per-origin paths
    out$store$s[s] <- s
    out$store$t_end[s] <- t_end
    
    out$store$A_q[, , s] <- A_q
    out$store$B_q[, , s] <- B_q
    out$store$sv_A_q[s, ] <- svd(A_q, nu = 0, nv = 0)$d
    out$store$sv_B_q[s, ] <- svd(B_q, nu = 0, nv = 0)$d
    
    if (isTRUE(store_cp)) {
      out$store$A_cp[, , s] <- A_cp
      out$store$B_cp[, , s] <- B_cp
      out$store$sv_A_cp[s, ] <- svd(A_cp, nu = 0, nv = 0)$d
      out$store$sv_B_cp[s, ] <- svd(B_cp, nu = 0, nv = 0)$d
    }
    
    # --- Forecast for each horizon
    for (h in horizon) {
      if (t_end + h > Tn) next
      idx_eval <- s
      
      # CP: VAR on factors
      cp_path <- generate_fhat_h_seq_stable(
        f_train = f_cp, Tn = train_T + h, train_T = train_T, h = h,
        pmax = pmax_var, include_const = include_const,
        ic = var_ic, stable_thresh = stable_thresh
      )
      fhat_cp_h <- cp_path$fhat_h[nrow(cp_path$fhat_h), , drop = FALSE]
      Yhat_cp <- reconstruct_Y_from_ABf(A_cp, B_cp, fhat_cp_h[1, ])
      
      # QMLE: VAR on factors
      q_path <- generate_fhat_h_seq_stable(
        f_train = f_q, Tn = train_T + h, train_T = train_T, h = h,
        pmax = pmax_var, include_const = include_const,
        ic = var_ic, stable_thresh = stable_thresh
      )
      fhat_q_h <- q_path$fhat_h[nrow(q_path$fhat_h), , drop = FALSE]
      Yhat_q <- reconstruct_Y_from_ABf(A_q, B_q, fhat_q_h[1, ])
      
      # Truth at global time t_end + h
      Ytrue <- matrix(R_npq[t_end + h, , ], nrow = p, ncol = q)
      
      Ecp <- Yhat_cp - Ytrue
      Eq  <- Yhat_q  - Ytrue
      
      # --- SD-normalized relative metrics (elementwise) ---
      out$per_h[[paste0("h", h)]]$rRMSE_cp[idx_eval] <- .rRMSE_sd(Ecp, sd_ij)
      out$per_h[[paste0("h", h)]]$rMAE_cp[idx_eval]  <- .rMAE_sd(Ecp, sd_ij)
      out$per_h[[paste0("h", h)]]$rRMSE_q[idx_eval]  <- .rRMSE_sd(Eq,  sd_ij)
      out$per_h[[paste0("h", h)]]$rMAE_q[idx_eval]   <- .rMAE_sd(Eq,  sd_ij)
      
      # VAR orders
      out$per_h[[paste0("h", h)]]$var_p_cp[idx_eval] <- cp_path$p_selected
      out$per_h[[paste0("h", h)]]$var_p_q[idx_eval]  <- q_path$p_selected
      
      # optional: store absolute errors and sd diagnostics
      if (isTRUE(store_abs_errors)) {
        out$per_h[[paste0("h", h)]]$RMSE_cp[idx_eval] <- .RMSE_abs(Ecp)
        out$per_h[[paste0("h", h)]]$MAE_cp[idx_eval]  <- .MAE_abs(Ecp)
        out$per_h[[paste0("h", h)]]$RMSE_q[idx_eval]  <- .RMSE_abs(Eq)
        out$per_h[[paste0("h", h)]]$MAE_q[idx_eval]   <- .MAE_abs(Eq)
        
        sds <- as.numeric(sd_ij)
        out$per_h[[paste0("h", h)]]$sd_mean[idx_eval] <- mean(sds)
        out$per_h[[paste0("h", h)]]$sd_med[idx_eval]  <- median(sds)
        out$per_h[[paste0("h", h)]]$sd_min[idx_eval]  <- min(sds)
        out$per_h[[paste0("h", h)]]$sd_p10[idx_eval]  <- as.numeric(stats::quantile(sds, probs = 0.10, names = FALSE))
      }
    }
    
    if (!is.null(verbose_every) && verbose_every > 0 && (s %% verbose_every == 0)) {
      cat(sprintf("Rolling done: s=%d / %d (window=[%d,%d])\n", s, S, s, t_end))
      for (h in horizon) {
        v <- out$per_h[[paste0("h", h)]]
        cat(sprintf(
          "  h=%d: CP mean rRMSE=%.4f, QMLE mean rRMSE=%.4f (so far)\n",
          h,
          mean(v$rRMSE_cp[1:s], na.rm = TRUE),
          mean(v$rRMSE_q[1:s],  na.rm = TRUE)
        ))
      }
    }
  }
  
  out$summary <- lapply(horizon, function(h) {
    v <- out$per_h[[paste0("h", h)]]
    res <- list(
      h = h,
      CP_mean_rRMSE = mean(v$rRMSE_cp, na.rm = TRUE),
      CP_sd_rRMSE   = sd(v$rRMSE_cp,   na.rm = TRUE),
      CP_mean_rMAE  = mean(v$rMAE_cp,  na.rm = TRUE),
      CP_sd_rMAE    = sd(v$rMAE_cp,    na.rm = TRUE),
      
      Q_mean_rRMSE  = mean(v$rRMSE_q,  na.rm = TRUE),
      Q_sd_rRMSE    = sd(v$rRMSE_q,    na.rm = TRUE),
      Q_mean_rMAE   = mean(v$rMAE_q,   na.rm = TRUE),
      Q_sd_rMAE     = sd(v$rMAE_q,     na.rm = TRUE)
    )
    if (isTRUE(store_abs_errors)) {
      res$CP_mean_RMSE <- mean(v$RMSE_cp, na.rm = TRUE)
      res$Q_mean_RMSE  <- mean(v$RMSE_q,  na.rm = TRUE)
      res$mean_sd_mean <- mean(v$sd_mean, na.rm = TRUE)
      res$mean_sd_med  <- mean(v$sd_med,  na.rm = TRUE)
      res$mean_sd_p10  <- mean(v$sd_p10,  na.rm = TRUE)
    }
    res
  })
  names(out$summary) <- paste0("h", horizon)
  
  out
}
