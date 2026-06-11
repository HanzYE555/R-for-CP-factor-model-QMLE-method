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
# 2) Helpers
# ============================================================
trim_to_NA <- function(x) x[is.finite(x)]

safe_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-12) return(1.0)
  s
}

# CP wrapper (robust)
fit_cpmts_unified_safe <- function(X_npq, d, max_tries = 3) {
  # X_npq: n x p x q (rolling train window)
  # returns list(A, B, Fhat)
  ok <- FALSE
  last_err <- NULL

  for (k in 1:max_tries) {
    out <- tryCatch({
      # CP_MTS expects time x p x q
      fit <- HDTSA::CP_MTS(X_npq, k = d, method = "CP.Unified")

      # common outputs in HDTSA: A, B, F or Fhat (names vary by version)
      A <- fit$A
      B <- fit$B
      F <- NULL
      if (!is.null(fit$Fhat)) F <- fit$Fhat
      if (is.null(F) && !is.null(fit$F)) F <- fit$F

      if (is.null(A) || is.null(B)) stop("CP_MTS did not return A/B.")

      # Ensure dimensions: A pxd, B qxd
      list(A = A, B = B, Fhat = F)
    }, error = function(e) {
      last_err <<- e
      NULL
    })

    if (!is.null(out)) { ok <- TRUE; return(out) }
  }

  stop(sprintf("CP_MTS failed after %d tries. Last error: %s",
               max_tries, if (is.null(last_err)) "NA" else last_err$message))
}

# extract BLUP factors given A,B and observation R_t (p x q)
extract_f_blup <- function(R_t, A, B) {
  # vec(R) approx (B \kron A) vec(F)
  # so vec(F) = solve( (B'B \kron A'A), (B' \kron A') vec(R) )
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  stopifnot(all(dim(R_t) == c(p, q)), ncol(B) == d)

  AtA <- crossprod(A)   # dxd
  BtB <- crossprod(B)   # dxd
  K   <- kronecker(BtB, AtA)               # d^2 x d^2
  rhs <- as.vector(t(A) %*% R_t %*% B)     # d^2 (since A' R B is dxd, then vec)
  fvec <- solve(K, rhs)
  matrix(fvec, nrow = d, ncol = d)
}

# stability check for VAR coefficients (companion matrix spectral radius)
is_var_stable <- function(coef_mat, k, p, thresh = 0.999) {
  # coef_mat: k x (k*p + const?) from MTS::VAR order
  # We will build companion using first k*p columns
  if (p <= 0) return(TRUE)
  A <- coef_mat[, 1:(k*p), drop = FALSE]
  # companion: kp x kp
  C <- matrix(0, nrow = k*p, ncol = k*p)
  C[1:k, ] <- A
  if (p > 1) C[(k+1):(k*p), 1:(k*(p-1))] <- diag(k*(p-1))
  rho <- max(Mod(eigen(C, only.values = TRUE)$values))
  is.finite(rho) && (rho < thresh)
}

# Fit VAR with IC selection + stability filter
fit_var_select <- function(Y, pmax = 24, ic = "AIC", include_const = TRUE, stable_thresh = 0.999) {
  # Y: T x k
  stopifnot(is.matrix(Y))
  k <- ncol(Y)
  Tn <- nrow(Y)
  if (Tn < 10) stop("Not enough data to fit VAR.")

  best <- NULL
  best_ic <- Inf
  best_p <- NA_integer_

  for (p in 1:pmax) {
    fit <- tryCatch({
      MTS::VAR(Y, p = p, include.mean = include_const)
    }, error = function(e) NULL)
    if (is.null(fit)) next

    # IC from fit: fit$aic, fit$bic (MTS typically provides)
    val <- if (ic == "AIC") fit$aic else fit$bic
    if (!is.finite(val)) next

    # stability check
    coef_mat <- fit$coef
    if (!is.null(coef_mat)) {
      if (!is_var_stable(coef_mat, k = k, p = p, thresh = stable_thresh)) next
    }

    if (val < best_ic) {
      best_ic <- val
      best <- fit
      best_p <- p
    }
  }

  if (is.null(best)) stop("No stable VAR fit found up to pmax.")
  list(fit = best, p = best_p, ic_value = best_ic)
}

# multi-step forecast from VAR (recursive)
var_forecast_h <- function(fit_obj, Y_hist, h) {
  # fit_obj from MTS::VAR
  # We use MTS::VARpred for 1-step? To keep simple: use VARpred which returns multi-step.
  pr <- MTS::VARpred(fit_obj, h = h)
  # VARpred returns $pred : h x k
  pr$pred[h, ]
}

# ============================================================
# 3) QMLE (your idea: SPD Mx via Cholesky parametrization)
#     NOTE: This is a placeholder wrapper around your existing QMLE.
#     If you already have a working fit_qmle(...) in your repo, paste it here.
# ============================================================

# ---- paste your existing fit_qmle here ----
# Must return list(A_q, B_q, Mx, converged, ...)
# For now, we provide a minimal "fallback" that uses CP init only (so code runs),
# but you SHOULD replace this with your actual QMLE implementation.
fit_qmle <- function(X_npq, d, init_A, init_B, maxit = 120) {
  # Replace with real QMLE.
  # Fallback: just return init (so reviewer can run end-to-end)
  list(A_q = init_A, B_q = init_B, Mx = diag(d), converged = TRUE)
}

# ============================================================
# 4) Rolling evaluation: CP.Unified vs QMLE
# ============================================================
compare_cpmts_vs_qmle_rolling_rrmse <- function(
    R_npq,
    d = 4,
    train_T = 456,
    horizon = c(1, 2),
    pmax_var = 24,
    include_const = TRUE,
    var_ic = "AIC",
    stable_thresh = 0.999,
    qmle_maxit = 120,
    verbose_every = 20,
    use_warm_start = TRUE,
    eps_denom = 1e-12,
    store_abs_errors = TRUE,
    s_start = 1,
    n_windows = NULL
) {
  stopifnot(length(dim(R_npq)) == 3)
  Tn <- dim(R_npq)[1]; p <- dim(R_npq)[2]; q <- dim(R_npq)[3]
  S_total <- Tn - train_T
  if (S_total < 1) stop("train_T too large for T.")

  s_start <- as.integer(s_start)
  if (s_start < 1 || s_start > S_total) stop("Invalid s_start.")

  s_end <- S_total
  if (!is.null(n_windows)) {
    n_windows <- as.integer(n_windows)
    if (n_windows < 1) stop("n_windows must be >= 1.")
    s_end <- min(S_total, s_start + n_windows - 1)
  }

  s_grid <- s_start:s_end
  S_run <- length(s_grid)

  # Precompute sd_ij over the whole sample (or you can use training-only; keep your current definition)
  sd_ij <- matrix(NA_real_, p, q)
  for (i in 1:p) for (j in 1:q) sd_ij[i, j] <- safe_sd(R_npq[, i, j])

  out <- list(
    d = d, train_T = train_T, T = Tn,
    p = p, q = q,
    horizon = horizon,
    S_total = S_total,
    S_run = S_run,
    s_start = s_start,
    s_end = s_end,
    per_h = list(),
    summary = list()
  )

  # store abs errors if requested
  if (store_abs_errors) {
    out$abs_errors <- list()
    for (h in horizon) {
      out$abs_errors[[paste0("h", h)]] <- list(
        cp = array(NA_real_, dim = c(S_run, p, q)),
        q  = array(NA_real_, dim = c(S_run, p, q))
      )
    }
  }

  for (h in horizon) {
    out$per_h[[paste0("h", h)]] <- list(
      rRMSE_cp = rep(NA_real_, S_run),
      rMAE_cp  = rep(NA_real_, S_run),
      rRMSE_q  = rep(NA_real_, S_run),
      rMAE_q   = rep(NA_real_, S_run),
      var_p_cp = rep(NA_integer_, S_run),
      var_p_q  = rep(NA_integer_, S_run)
    )
  }

  # warm start holders
  last_Aq <- NULL; last_Bq <- NULL

  for (kk in seq_along(s_grid)) {
    s <- s_grid[kk]
    t_end <- s + train_T - 1
    X_train <- R_npq[s:t_end, , , drop = FALSE]  # train_T x p x q

    # ---- CP.Unified on train window
    cp <- fit_cpmts_unified_safe(X_train, d = d)
    A_cp <- cp$A
    B_cp <- cp$B

    # ---- QMLE on train window (init by CP or warm start)
    init_A <- if (use_warm_start && !is.null(last_Aq)) last_Aq else A_cp
    init_B <- if (use_warm_start && !is.null(last_Bq)) last_Bq else B_cp

    qmle <- fit_qmle(X_train, d = d, init_A = init_A, init_B = init_B, maxit = qmle_maxit)
    A_q <- qmle$A_q
    B_q <- qmle$B_q

    last_Aq <- A_q; last_Bq <- B_q

    # ---- build factor time series over training (BLUP)
    Fcp <- array(NA_real_, dim = c(train_T, d, d))
    Fq  <- array(NA_real_, dim = c(train_T, d, d))

    for (tt in 1:train_T) {
      Rt <- X_train[tt, , ]
      Fcp[tt, , ] <- extract_f_blup(Rt, A_cp, B_cp)
      Fq[tt, , ]  <- extract_f_blup(Rt, A_q,  B_q)
    }

    # ---- vectorize factors for VAR: vec(F_t) length d^2
    Ycp <- matrix(NA_real_, nrow = train_T, ncol = d*d)
    Yq  <- matrix(NA_real_, nrow = train_T, ncol = d*d)
    for (tt in 1:train_T) {
      Ycp[tt, ] <- as.vector(Fcp[tt, , ])
      Yq[tt, ]  <- as.vector(Fq[tt, , ])
    }

    # ---- VAR fit + forecast for each h
    var_cp <- fit_var_select(Ycp, pmax = pmax_var, ic = var_ic, include_const = include_const, stable_thresh = stable_thresh)
    var_q  <- fit_var_select(Yq,  pmax = pmax_var, ic = var_ic, include_const = include_const, stable_thresh = stable_thresh)

    for (h in horizon) {
      t_pred <- t_end + h
      if (t_pred > Tn) next

      yhat_cp <- var_forecast_h(var_cp$fit, Ycp, h)
      yhat_q  <- var_forecast_h(var_q$fit,  Yq,  h)

      Fhat_cp <- matrix(yhat_cp, nrow = d, ncol = d)
      Fhat_q  <- matrix(yhat_q,  nrow = d, ncol = d)

      Rhat_cp <- A_cp %*% Fhat_cp %*% t(B_cp)
      Rhat_q  <- A_q  %*% Fhat_q  %*% t(B_q)

      R_true <- R_npq[t_pred, , ]

      err_cp <- Rhat_cp - R_true
      err_q  <- Rhat_q  - R_true

      # sd-normalized rRMSE / rMAE
      z_cp <- err_cp / (sd_ij + eps_denom)
      z_q  <- err_q  / (sd_ij + eps_denom)

      rrmse_cp <- sqrt(mean(z_cp^2, na.rm = TRUE))
      rrmse_q  <- sqrt(mean(z_q^2,  na.rm = TRUE))
      rmae_cp  <- mean(abs(z_cp), na.rm = TRUE)
      rmae_q   <- mean(abs(z_q),  na.rm = TRUE)

      out$per_h[[paste0("h", h)]]$rRMSE_cp[kk] <- rrmse_cp
      out$per_h[[paste0("h", h)]]$rRMSE_q[kk]  <- rrmse_q
      out$per_h[[paste0("h", h)]]$rMAE_cp[kk]  <- rmae_cp
      out$per_h[[paste0("h", h)]]$rMAE_q[kk]   <- rmae_q
      out$per_h[[paste0("h", h)]]$var_p_cp[kk] <- var_cp$p
      out$per_h[[paste0("h", h)]]$var_p_q[kk]  <- var_q$p

      if (store_abs_errors) {
        out$abs_errors[[paste0("h", h)]]$cp[kk, , ] <- abs(err_cp)
        out$abs_errors[[paste0("h", h)]]$q[kk, , ]  <- abs(err_q)
      }
    }

    if (!is.null(verbose_every) && verbose_every > 0 && (kk %% verbose_every == 0)) {
      cat(sprintf("Done kk=%d/%d (s=%d, t_end=%d)\n", kk, S_run, s, t_end))
      for (h in horizon) {
        a <- out$per_h[[paste0("h", h)]]$rRMSE_cp[1:kk]
        b <- out$per_h[[paste0("h", h)]]$rRMSE_q[1:kk]
        cat(sprintf("  h=%d: mean rRMSE CP=%.4f, Q=%.4f\n",
                    h, mean(trim_to_NA(a)), mean(trim_to_NA(b))))
      }
    }
  }

  # summary
  for (h in horizon) {
    a <- out$per_h[[paste0("h", h)]]$rRMSE_cp
    b <- out$per_h[[paste0("h", h)]]$rRMSE_q
    ma <- out$per_h[[paste0("h", h)]]$rMAE_cp
    mb <- out$per_h[[paste0("h", h)]]$rMAE_q

    out$summary[[paste0("h", h)]] <- list(
      h = h,
      CP_mean_rRMSE = mean(trim_to_NA(a)),
      CP_sd_rRMSE   = sd(trim_to_NA(a)),
      Q_mean_rRMSE  = mean(trim_to_NA(b)),
      Q_sd_rRMSE    = sd(trim_to_NA(b)),
      CP_mean_rMAE  = mean(trim_to_NA(ma)),
      Q_mean_rMAE   = mean(trim_to_NA(mb)),
      n_eval        = sum(is.finite(a) & is.finite(b))
    )
  }

  out
}

# ============================================================
# 5) Load prepared data (R_npq or E_capm)
# ============================================================
dat <- readRDS(data_path)

if (!is.null(dat$R_npq)) {
  R_npq <- dat$R_npq
} else if (!is.null(dat$E_capm)) {
  E_capm <- dat$E_capm
  stopifnot(is.matrix(E_capm), ncol(E_capm) == 100)
  Tn <- nrow(E_capm)
  R_t <- vector("list", length = Tn)
  for (t in seq_len(Tn)) R_t[[t]] <- 100 * matrix(E_capm[t, ], nrow = 10)
  R_arr_pqn <- simplify2array(R_t)
  R_npq <- aperm(R_arr_pqn, c(3, 1, 2))
} else {
  stop("RDS must contain R_npq or E_capm.")
}

stopifnot(length(dim(R_npq)) == 3)
cat(sprintf("Loaded R_npq dim: %s\n", paste(dim(R_npq), collapse = " x ")))

# ============================================================
# 6) Run + save
# ============================================================
res <- compare_cpmts_vs_qmle_rolling_rrmse(
  R_npq = R_npq,
  d = d,
  train_T = train_T,
  horizon = horizon,
  pmax_var = pmax_var,
  include_const = TRUE,
  var_ic = var_ic,
  stable_thresh = stable_thresh,
  qmle_maxit = qmle_maxit,
  verbose_every = verbose_every,
  use_warm_start = TRUE,
  store_abs_errors = TRUE,
  s_start = s_start,
  n_windows = n_windows
)

tag <- sprintf("d%d_trainT%d_s%03d_n%s_%s",
               d, train_T, s_start,
               if (is.null(n_windows)) "ALL" else sprintf("%03d", n_windows),
               var_ic)

rds_out <- file.path(out_dir, paste0("rolling_", tag, ".rds"))
saveRDS(res, rds_out)

# summary CSV
sum_df <- do.call(rbind, lapply(names(res$summary), function(nm) {
  x <- res$summary[[nm]]
  data.frame(
    horizon = x$h,
    CP_mean_rRMSE = x$CP_mean_rRMSE,
    CP_sd_rRMSE   = x$CP_sd_rRMSE,
    Q_mean_rRMSE  = x$Q_mean_rRMSE,
    Q_sd_rRMSE    = x$Q_sd_rRMSE,
    CP_mean_rMAE  = x$CP_mean_rMAE,
    Q_mean_rMAE   = x$Q_mean_rMAE,
    n_eval = x$n_eval
  )
}))
csv_out <- file.path(out_dir, paste0("summary_", tag, ".csv"))
write.csv(sum_df, csv_out, row.names = FALSE)

cat("============================================================\n")
cat("Saved outputs:\n")
cat("  ", rds_out, "\n", sep = "")
cat("  ", csv_out, "\n", sep = "")
cat("Summary:\n")
print(sum_df)
cat("============================================================\n")
