#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MASS)
  library(MTS)
})

# =========================
# Utilities (shared)
# =========================

# Column-normalize matrix to target L2 norm
col_normalize <- function(M, target_norm) {
  cn <- sqrt(colSums(M^2))
  cn[cn == 0] <- 1
  M <- sweep(M, 2, cn, "/")
  M * target_norm
}

# Project matrix to rank r via SVD, keep top-r components
proj_rank <- function(M, r) {
  s <- svd(M, nu = nrow(M), nv = ncol(M))
  if (r > 0) {
    U_r <- s$u[, 1:r, drop = FALSE]
    D_r <- diag(s$d[1:r], nrow = r, ncol = r)
    V_r <- s$v[, 1:r, drop = FALSE]
    U_r %*% D_r %*% t(V_r)
  } else {
    matrix(0, nrow = nrow(M), ncol = ncol(M))
  }
}

# Force (numerical) rank = r by truncation, then re-pad columns to desired norms
force_rank_and_normalize <- function(M, r, target_norm) {
  Mr <- proj_rank(M, r)
  col_normalize(Mr, target_norm)
}

# Khatri-Rao columns: vec(a_c b_c^T) = b_c ⊗ a_c
odotprod <- function(A, B) {
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  K <- matrix(0, nrow = p * q, ncol = d)
  for (c in 1:d) K[, c] <- kronecker(B[, c], A[, c])
  K
}

# Upper-triangular packing for symmetric Mx
n_upper <- function(d) d * (d + 1) / 2

pack_upper <- function(M) {
  d <- nrow(M)
  v <- numeric(n_upper(d))
  idx <- 1
  for (i in 1:d) {
    for (j in i:d) {
      v[idx] <- M[i, j]
      idx <- idx + 1
    }
  }
  v
}

unpack_upper <- function(v, d) {
  M <- matrix(0, d, d)
  idx <- 1
  for (i in 1:d) {
    for (j in i:d) {
      M[i, j] <- v[idx]
      idx <- idx + 1
    }
  }
  M <- M + t(M) - diag(diag(M))
  M
}

pack_params <- function(A, B, Mx_upper, log_diag_Se) {
  c(as.vector(A), as.vector(B), as.vector(Mx_upper), as.vector(log_diag_Se))
}

unpack_params <- function(theta, p, q, d) {
  nA <- p * d
  nB <- q * d
  nMxU <- n_upper(d)
  nSe <- p * q

  A <- matrix(theta[1:nA], nrow = p, ncol = d)
  B <- matrix(theta[(nA + 1):(nA + nB)], nrow = q, ncol = d)

  off1 <- nA + nB
  Mx_upper <- theta[(off1 + 1):(off1 + nMxU)]

  off2 <- off1 + nMxU
  log_diag_Se <- theta[(off2 + 1):(off2 + nSe)]

  list(A = A, B = B, Mx_upper = Mx_upper, log_diag_Se = log_diag_Se)
}

empirical_second_moment <- function(Ylist, p, q, Tn) {
  M <- matrix(0, nrow = p * q, ncol = p * q)
  for (t in 1:Tn) {
    v <- as.vector(Ylist[[t]])
    M <- M + tcrossprod(v)
  }
  M / Tn
}

# White noise: entrywise independent, diag covariance for vec(E_t)
gnrt_matrix_normal_ts <- function(p, q, Tn, diag_Se_entry, s) {
  set.seed(s + 12345)
  pq <- p * q
  sd_vec <- sqrt(diag_Se_entry)
  Elist <- vector("list", Tn)
  for (t in 1:Tn) {
    z <- rnorm(pq) * sd_vec
    Elist[[t]] <- matrix(z, nrow = p, ncol = q)
  }
  Elist
}

# VARMA latent factor generator: returns list length Tn of d-vectors
gnrt_varmalatent <- function(d, s, Tn) {
  set.seed(s)
  k <- d
  n_obs <- Tn + 100

  generate_stable_matrix <- function(k, scale = 0.3) {
    M <- matrix(rnorm(k^2, sd = scale), k, k)
    eig <- max(Mod(eigen(M)$values))
    if (eig >= 1) M <- M / (eig + 0.2)
    M
  }

  phi <- generate_stable_matrix(k)
  theta <- generate_stable_matrix(k)

  diagvalue_sigma_latent <- sample(c(1, 2), size = k, replace = TRUE)
  sigma_latent <- diag(diagvalue_sigma_latent)

  sim_data <- VARMAsim(
    nobs = n_obs, arlags = 1, malags = 1,
    phi = phi, theta = theta, sigma = sigma_latent
  )
  F_all <- sim_data$series

  Flist <- vector("list", length = Tn)
  for (t in 1:Tn) Flist[[t]] <- 3 * F_all[t + 100, ]
  Flist
}

# QMLE core (shared)
q_function <- function(M_y, A, B, M_x, log_diag_Se, ridge = 1e-8) {
  p <- nrow(A); q <- nrow(B); d <- ncol(A)
  N <- p * q

  K <- odotprod(A, B)

  Se_diag <- exp(log_diag_Se)
  Se_eff <- Se_diag + ridge

  Dinv <- 1 / Se_eff
  logdet_D <- sum(log(Se_eff))

  # G = K' D^{-1} K
  W <- K * Dinv
  G <- crossprod(K, W)

  # Cholesky for Mx (ensure PD)
  Mx <- M_x
  chol_Mx <- try(chol(Mx), silent = TRUE)
  if (inherits(chol_Mx, "try-error")) {
    Mx <- Mx + 1e-8 * diag(d)
    chol_Mx <- chol(Mx)
  }
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

  tr_Dinv_My <- sum(Dinv * diag(M_y))

  U <- W
  V <- M_y %*% U
  V <- V * Dinv
  Tmat <- crossprod(K, V)

  SinvT <- backsolve(chol_S, forwardsolve(t(chol_S), Tmat))
  tr_term <- tr_Dinv_My - sum(diag(SinvT))

  (logdet_SigmaY + tr_term) / N
}

objective_fn <- function(theta, M_y, p, q, d) {
  pars <- unpack_params(theta, p, q, d)
  A <- col_normalize(pars$A, sqrt(p))
  B <- col_normalize(pars$B, sqrt(q))
  M_x <- unpack_upper(pars$Mx_upper, d) + 1e-8 * diag(d)
  q_function(M_y, A, B, M_x = M_x, log_diag_Se = pars$log_diag_Se, ridge = 1e-8)
}

# =========================
# Boundary-case specific A,B construction
# =========================

make_Ar_Br <- function() {
  MA <- rbind(
    c( 1,  5,  2, -1, 0),
    c( 2,  1, -2,  5, 1),
    c(-1, -5, -2,  1, 0),
    c( 0,  2,  1, -1, 0),
    c( 3,  4, -1,  5, 1)
  )
  MB <- rbind(
    c( 2,  1,  0,  0, -1),
    c( 0, -1,  0, -2,  5),
    c( 1,  0, -1, -2,  3),
    c(-1,  0,  0,  1, -2),
    c( 3,  0,  3,  0,  3)
  )

  DA <- diag(c(1/sqrt(15), 1/sqrt(71), 1/sqrt(14), 1/sqrt(53), 1/sqrt(2)), 5, 5)
  DB <- diag(c(1/sqrt(15), 1/sqrt(2), 1/sqrt(10), 1/3, 1/sqrt(48)), 5, 5)

  Ar <- MA %*% DA
  Br <- MB %*% DB
  list(Ar = Ar, Br = Br)
}

make_AB_from_QR <- function(p, q, d, seed) {
  if (d > 5) stop("Boundary construction requires d <= 5 (A_r, B_r are 5x5).")
  set.seed(seed + 202403)

  Atilde <- matrix(rnorm(p * d, sd = sqrt(3)), p, d)
  Btilde <- matrix(rnorm(q * d, sd = sqrt(3)), q, d)

  QA <- qr.Q(qr(Atilde))
  QB <- qr.Q(qr(Btilde))

  if (ncol(QA) < d) QA <- cbind(QA, matrix(0, nrow = p, ncol = d - ncol(QA)))
  if (ncol(QB) < d) QB <- cbind(QB, matrix(0, nrow = q, ncol = d - ncol(QB)))

  QA <- QA[, 1:d, drop = FALSE]
  QB <- QB[, 1:d, drop = FALSE]

  base <- make_Ar_Br()
  Ar <- base$Ar[, 1:d, drop = FALSE]
  Br <- base$Br[, 1:d, drop = FALSE]

  A <- sqrt(p) * (QA %*% Ar)
  B <- sqrt(q) * (QB %*% Br)

  list(A = A, B = B)
}

# =========================
# Data generation / fit per case
# =========================

simulate_data_bdry <- function(p, q, d, Tn, s) {
  set.seed(s)
  AB <- make_AB_from_QR(p = p, q = q, d = d, seed = s)
  A0 <- AB$A
  B0 <- AB$B

  F_t <- gnrt_varmalatent(d, s, Tn)

  M_x0 <- matrix(0, d, d)
  for (t in 1:Tn) {
    ft <- matrix(F_t[[t]], nrow = d, ncol = 1)
    M_x0 <- M_x0 + (ft %*% t(ft))
  }
  M_x0 <- M_x0 / Tn

  pq <- p * q
  diag_Se_entry0 <- runif(pq, 1, 2)

  Elist <- gnrt_matrix_normal_ts(p, q, Tn, diag_Se_entry0, s)

  Ylist <- vector("list", Tn)
  for (t in 1:Tn) {
    ft <- F_t[[t]]
    mean_mat <- A0 %*% diag(ft, nrow = d, ncol = d) %*% t(B0)
    Ylist[[t]] <- mean_mat + Elist[[t]]
  }

  M_y <- empirical_second_moment(Ylist, p, q, Tn)
  list(A0 = A0, B0 = B0, M_x0 = M_x0, diag_Se_entry0 = diag_Se_entry0, M_y = M_y, Ylist = Ylist)
}

simulate_data_int <- function(p, q, d, Tn, s) {
  set.seed(s)

  A0 <- matrix(rnorm(p * d, sd = sqrt(3)), nrow = p, ncol = d)
  B0 <- matrix(rnorm(q * d, sd = sqrt(3)), nrow = q, ncol = d)
  A0 <- col_normalize(A0, sqrt(p))
  B0 <- col_normalize(B0, sqrt(q))

  # interior case uses rank (d-1)
  A0 <- force_rank_and_normalize(A0, d - 1, sqrt(p))
  B0 <- force_rank_and_normalize(B0, d - 1, sqrt(q))

  F_t <- gnrt_varmalatent(d, s, Tn)

  M_x0 <- matrix(0, d, d)
  for (t in 1:Tn) {
    ft <- matrix(F_t[[t]], nrow = d, ncol = 1)
    M_x0 <- M_x0 + (ft %*% t(ft))
  }
  M_x0 <- M_x0 / Tn

  pq <- p * q
  diag_Se_entry0 <- runif(pq, 1, 2)

  Elist <- gnrt_matrix_normal_ts(p, q, Tn, diag_Se_entry0, s)

  Ylist <- vector("list", Tn)
  for (t in 1:Tn) {
    ft <- F_t[[t]]
    mean_mat <- A0 %*% diag(ft, nrow = d, ncol = d) %*% t(B0)
    Ylist[[t]] <- mean_mat + Elist[[t]]
  }

  M_y <- empirical_second_moment(Ylist, p, q, Tn)
  list(A0 = A0, B0 = B0, M_x0 = M_x0, diag_Se_entry0 = diag_Se_entry0, M_y = M_y, Ylist = Ylist)
}

fit_qmle_bdry <- function(M_y, p, q, d, init, maxit = 50) {
  Ainit <- init$A
  Binit <- init$B
  M_x_init <- init$M_x
  log_Se_init <- init$log_diag_Se

  theta0 <- pack_params(Ainit, Binit, pack_upper(M_x_init), log_Se_init)

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

  # boundary case: enforce rank d and normalize
  A1 <- force_rank_and_normalize(pars$A, d, sqrt(p))
  B1 <- force_rank_and_normalize(pars$B, d, sqrt(q))

  list(
    Ahat = A1, Bhat = B1,
    Mxhat = unpack_upper(pars$Mx_upper, d),
    logSehat = pars$log_diag_Se,
    opt = opt,
    conv_code = opt$convergence,
    conv_msg = opt$message
  )
}

fit_qmle_int <- function(M_y, p, q, d, init, maxit = 50) {
  Ainit <- init$A
  Binit <- init$B
  M_x_init <- init$M_x
  log_Se_init <- init$log_diag_Se

  theta0 <- pack_params(Ainit, Binit, pack_upper(M_x_init), log_Se_init)

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

  # interior case: enforce rank (d-1) and normalize
  A1 <- force_rank_and_normalize(pars$A, d - 1, sqrt(p))
  B1 <- force_rank_and_normalize(pars$B, d - 1, sqrt(q))

  list(
    Ahat = A1, Bhat = B1,
    Mxhat = unpack_upper(pars$Mx_upper, d),
    logSehat = pars$log_diag_Se,
    opt = opt,
    conv_code = opt$convergence,
    conv_msg = opt$message
  )
}

# =========================
# Run once / many (shared)
# =========================

extract_errors <- function(one_res) {
  c(
    A_relF_err = one_res$A_relF_err,
    B_relF_err = one_res$B_relF_err,
    Se_relRMSE = one_res$Se_relRMSE,
    Mx_relF_err = one_res$Mx_relF_err
  )
}

run_once_case <- function(case = c("bdry", "int"), p = 20, q = 20, d = 4, Tn = 100, s = 1, maxit = 50) {
  case <- match.arg(case)

  sim <- switch(case,
    bdry = simulate_data_bdry(p, q, d, Tn, s),
    int  = simulate_data_int(p, q, d, Tn, s)
  )

  init <- list(
    A = sim$A0,
    B = sim$B0,
    M_x = sim$M_x0,
    log_diag_Se = log(sim$diag_Se_entry0)
  )

  fit <- switch(case,
    bdry = fit_qmle_bdry(sim$M_y, p, q, d, init, maxit),
    int  = fit_qmle_int(sim$M_y, p, q, d, init, maxit)
  )

  Ahat <- col_normalize(fit$Ahat, sqrt(p))
  Bhat <- col_normalize(fit$Bhat, sqrt(q))

  # Errors
  A_err <- norm(sim$A0 - Ahat, "F") / max(1, norm(sim$A0, "F"))
  B_err <- norm(sim$B0 - Bhat, "F") / max(1, norm(sim$B0, "F"))
  Se_true <- sim$diag_Se_entry0
  Se_hat <- exp(fit$logSehat)
  Se_relRMSE <- sqrt(mean((Se_hat - Se_true)^2)) / max(1, mean(Se_true))
  Mx_err <- norm(sim$M_x0 - fit$Mxhat, "F") / max(1, norm(sim$M_x0, "F"))

  list(
    case = case,
    A0 = sim$A0, B0 = sim$B0,
    Ahat = Ahat, Bhat = Bhat,
    Mx0 = sim$M_x0, Mxhat = fit$Mxhat,
    Se_true = Se_true, Se_hat = Se_hat,
    M_y = sim$M_y,
    opt = fit$opt,
    converged = isTRUE(fit$conv_code == 0),
    conv_code = fit$conv_code,
    conv_msg = fit$conv_msg,
    A_relF_err = A_err,
    B_relF_err = B_err,
    Se_relRMSE = Se_relRMSE,
    Mx_relF_err = Mx_err
  )
}

run_many <- function(case, p, q, d, Tn, seeds, maxit = 50, verbose = TRUE) {
  out_mat <- matrix(NA_real_, nrow = length(seeds), ncol = 5)
  colnames(out_mat) <- c("conv_code", "A_relF_err", "B_relF_err", "Se_relRMSE", "Mx_relF_err")

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    if (verbose) cat(sprintf("[%s] Running seed %d ...\n", case, s))
    one <- run_once_case(case = case, p = p, q = q, d = d, Tn = Tn, s = s, maxit = maxit)

    out_mat[i, "conv_code"] <- one$conv_code
    out_mat[i, c("A_relF_err", "B_relF_err", "Se_relRMSE", "Mx_relF_err")] <- extract_errors(one)
  }

  out_avg <- colMeans(out_mat, na.rm = TRUE)
  list(seeds = seeds, out_each = out_mat, out_avg = out_avg)
}

# =========================
# CLI
# =========================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}

if (length(args) == 0) {
  cat("Usage:\n")
  cat("  Rscript scripts/run_sim.R --case bdry --p 60 --q 60 --T 100 --d 4 --seeds 1,2,3,4,5 --maxit 30 --out out/bdry_p60_q60_T100.csv\n")
  cat("  Rscript scripts/run_sim.R --case int  --p 60 --q 60 --T 100 --d 4 --seeds 1,2,3,4,5 --maxit 30 --out out/int_p60_q60_T100.csv\n")
  quit(status = 0)
}

case  <- get_arg("--case", "bdry")
p     <- as.integer(get_arg("--p"))
q     <- as.integer(get_arg("--q"))
Tn    <- as.integer(get_arg("--T"))
d     <- as.integer(get_arg("--d", "4"))
seeds <- get_arg("--seeds", "1,2,3,4,5")
maxit <- as.integer(get_arg("--maxit", "50"))
out   <- get_arg("--out", sprintf("out/%s_p%d_q%d_T%d.csv", case, p, q, Tn))

if (is.na(p) || is.na(q) || is.na(Tn)) {
  stop("必须提供 --p --q --T，例如 --case bdry --p 60 --q 60 --T 100")
}

case <- tolower(case)
if (!case %in% c("bdry", "int")) stop("--case must be one of: bdry, int")

if (case == "bdry" && d > 5) stop("Boundary case requires d <= 5 (A_r, B_r are 5x5).")

seed_vec <- as.integer(strsplit(seeds, ",")[[1]])
agg <- run_many(case = case, p = p, q = q, d = d, Tn = Tn, seeds = seed_vec, maxit = maxit, verbose = TRUE)

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

df_each <- data.frame(
  case = case,
  p = p, q = q, T = Tn,
  seed = seed_vec,
  conv_code = agg$out_each[, "conv_code"],
  A_relF_err = agg$out_each[, "A_relF_err"],
  B_relF_err = agg$out_each[, "B_relF_err"],
  Se_relRMSE = agg$out_each[, "Se_relRMSE"],
  Mx_relF_err = agg$out_each[, "Mx_relF_err"]
)

df_avg <- data.frame(
  case = case,
  p = p, q = q, T = Tn,
  seed = NA_integer_,
  conv_code = agg$out_avg["conv_code"],
  A_relF_err = agg$out_avg["A_relF_err"],
  B_relF_err = agg$out_avg["B_relF_err"],
  Se_relRMSE = agg$out_avg["Se_relRMSE"],
  Mx_relF_err = agg$out_avg["Mx_relF_err"]
)

res <- rbind(df_each, df_avg)
write.csv(res, out, row.names = FALSE)
cat("Wrote:", out, "\n")
