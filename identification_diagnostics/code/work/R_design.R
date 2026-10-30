# Data-generating designs for the requested p=q=40 CP-QMLE experiment.
# Base R only. No packages, external data, or hidden fitted/truth parameters.
# Rows of Y index time; columns use column-major vec of a p by q matrix.

simulation_seeds <- seq_len(100L)
simulation_T <- c(100L, 200L, 400L, 800L)
simulation_regimes <- c("interior", "boundary", "exterior")

normalize_columns <- function(A, target_norm = 1) {
  norms <- sqrt(colSums(A * A))
  if (any(!is.finite(norms)) || any(norms <= 0)) stop("Zero or invalid loading column")
  sweep(A, 2L, norms / target_norm, "/")
}

khatri_rao <- function(A, B) {
  if (ncol(A) != ncol(B)) stop("A and B must have the same number of columns")
  vapply(seq_len(ncol(A)), function(j) as.vector(outer(A[, j], B[, j])),
         numeric(nrow(A) * nrow(B)))
}

with_local_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  eval.parent(substitute(expr))
}

seed_key <- function(seed, tag, dimension = 0L) {
  as.integer(seed + 100003L * tag + 1009L * dimension)
}

haar_thin_basis <- function(n, r, seed) {
  if (r > n || r < 1L) stop("The embedding rank must lie between 1 and n")
  with_local_seed(seed, {
    G <- matrix(rnorm(n * r), n, r)
    decomp <- qr(G)
    Q <- qr.Q(decomp)[, seq_len(r), drop = FALSE]
    signs <- sign(diag(qr.R(decomp)))
    signs[signs == 0] <- 1
    sweep(Q, 2L, signs, "*")
  })
}

core_loading_design <- function(regime, n = 0.5) {
  regime <- match.arg(regime, simulation_regimes)
  if (regime == "interior") {
    A <- rbind(c(-1, 1, 0, 1, 1, -1, 0, 0),
               c(1, 0, -1, -1, -1, 0, 1, 0),
               c(1, 0, 0, -1, 1, 1, -1, -1),
               c(-1, 1, 1, 0, 1, 0, 0, 1))
    B <- rbind(c(0, 1, 0, 1, 0, 1, 1, -1),
               c(-1, 1, 1, -1, 1, 0, -1, 1),
               c(0, 0, 1, 0, 0, 1, 0, -1),
               c(1, 1, 1, -1, 1, 0, 0, -1))
  } else {
    A <- rbind(c(1, 1, 0, 1, 0), c(0, 2, 1, -1, 0), c(0, 1, -1, 2, 1))
    B <- rbind(c(2, 1, -1, -1, 0), c(0, -1, 5, 3, 0), c(0, -1, 3, 1, 2))
    if (regime == "exterior") {
      if (!is.finite(n) || n <= 0) stop("n must be positive")
      B[, 1L] <- B[, 1L] - c(0, 2, 0) / n
    }
  }
  extra_lambda <- NULL
  if (regime == "exterior") {
    extra_lambda <- c(1, rep(1 / n, ncol(A) - 1L)) *
      sqrt(colSums(A * A)) * sqrt(colSums(B * B))
    extra_lambda <- extra_lambda / sqrt(sum(extra_lambda^2))
  }
  list(raw_A = A, raw_B = B, Ar_unit = normalize_columns(A),
       Br_unit = normalize_columns(B), d = ncol(A), rank_A = nrow(A),
       rank_B = nrow(B), regime = regime, extra_lambda = extra_lambda,
       n = if (regime == "exterior") n else NA_real_)
}

make_loading_design <- function(regime, seed, p = 40L, q = 40L, n = 0.5) {
  core <- core_loading_design(regime, n)
  # Boundary and exterior share exactly the same embedding for the same seed.
  QA <- haar_thin_basis(p, core$rank_A, seed_key(seed, 1L))
  QB <- haar_thin_basis(q, core$rank_B, seed_key(seed, 2L))
  A_unit <- QA %*% core$Ar_unit
  B_unit <- QB %*% core$Br_unit
  A <- sqrt(p) * A_unit
  B <- sqrt(q) * B_unit
  out <- c(core, list(A = A, B = B, A_unit = A_unit, B_unit = B_unit,
                     QA = QA, QB = QB, p = as.integer(p), q = as.integer(q),
                     seed = as.integer(seed)))
  if (!is.null(core$extra_lambda)) {
    zm <- drop(khatri_rao(A_unit, B_unit) %*% core$extra_lambda)
    singular <- svd(matrix(zm, p, q), nu = 0L, nv = 0L)$d
    stopifnot(singular[2L] < singular[1L] * 1e-10)
    out$extra_max_loading_cosine <- max(abs(drop(crossprod(khatri_rao(A_unit, B_unit), zm)))) /
      sqrt(sum(zm * zm))
  }
  out
}

stable_random_matrix <- function(d, radius_range) {
  A <- matrix(rnorm(d * d), d, d) / sqrt(d)
  radius <- max(Mod(eigen(A, only.values = TRUE)$values))
  A * (runif(1L, radius_range[1L], radius_range[2L]) / radius)
}

stationary_varma_covariance <- function(Phi, Theta, innovation_variance) {
  d <- nrow(Phi)
  zero <- matrix(0, d, d)
  F <- rbind(cbind(Phi, Theta), cbind(zero, zero))
  G <- rbind(diag(d), diag(d))
  W <- G %*% diag(innovation_variance, d) %*% t(G)
  m <- 2L * d
  V <- matrix(solve(diag(m * m) - kronecker(F, F), as.vector(W)), m, m)
  V <- (V + t(V)) / 2
  list(F = F, G = G, V = V, Mx = V[seq_len(d), seq_len(d), drop = FALSE])
}

simulate_varma <- function(d, seed, Tmax = 800L, burn = 100L,
                           factor_multiplier = 1) {
  if (Tmax < 1L || burn < 0L || factor_multiplier <= 0) stop("Invalid simulation settings")
  with_local_seed(seed_key(seed, 3L, d), {
    Phi <- stable_random_matrix(d, c(0.3, 0.7))
    Theta <- stable_random_matrix(d, c(0.1, 0.5))
    innovation_variance <- sample(c(1, 2), d, replace = TRUE)
    stationary <- stationary_varma_covariance(Phi, Theta, innovation_variance)
    ev <- eigen(stationary$V, symmetric = TRUE)
    state <- drop(ev$vectors %*% (sqrt(pmax(ev$values, 0)) * rnorm(2L * d)))
    innovations <- sweep(matrix(rnorm((burn + Tmax) * d), burn + Tmax, d, byrow = TRUE),
                         2L, sqrt(innovation_variance), "*")
    Xall <- matrix(NA_real_, burn + Tmax, d)
    for (tt in seq_len(burn + Tmax)) {
      state <- drop(stationary$F %*% state + stationary$G %*% innovations[tt, ])
      Xall[tt, ] <- factor_multiplier * state[seq_len(d)]
    }
    list(X = Xall[burn + seq_len(Tmax), , drop = FALSE], Phi = Phi, Theta = Theta,
         Mx = stationary$Mx * factor_multiplier^2,
         innovation_variance = innovation_variance)
  })
}

generate_data <- function(regime, seed, Tmax = 800L, p = 40L, q = 40L,
                          burn = 100L, factor_multiplier = 1, n = 0.5) {
  loadings <- make_loading_design(regime, seed, p, q, n)
  factors <- simulate_varma(loadings$d, seed, Tmax, burn, factor_multiplier)
  # A separate noise stream makes all regimes share the same standard-normal
  # observation noise for a given seed; boundary/exterior also share X and Q.
  noise <- with_local_seed(seed_key(seed, 4L), matrix(rnorm(Tmax * p * q), Tmax, p * q, byrow = TRUE))
  Z <- khatri_rao(loadings$A, loadings$B)
  Y <- factors$X %*% t(Z) + noise
  metadata <- list(regime = regime, seed = as.integer(seed), p = as.integer(p),
                   q = as.integer(q), d = loadings$d, rank_A = loadings$rank_A,
                   rank_B = loadings$rank_B, Tmax = as.integer(Tmax), burn = as.integer(burn),
                   factor_multiplier = factor_multiplier, noise_covariance = "identity",
                   noise_variance = 1, embedding = "Gaussian QR Haar thin basis, fixed within seed",
                   AR_radius = max(Mod(eigen(factors$Phi, only.values = TRUE)$values)),
                   MA_radius = max(Mod(eigen(factors$Theta, only.values = TRUE)$values)),
                   AR_target_range = c(0.3, 0.7), MA_target_range = c(0.1, 0.5),
                   factor_process = "x_t = Phi x_(t-1) + eps_t + Theta eps_(t-1)",
                   nested_T = TRUE, stationary_initialization = TRUE,
                   n = loadings$n,
                   signal_noise_trace_ratio = sum((Z %*% factors$Mx) * Z) / (p * q))
  c(loadings, factors, list(Y = Y, Z = Z, noise = noise,
                           noise_variance = rep(1, p * q), metadata = metadata))
}

centered_prefix <- function(design, T, make_S = FALSE) {
  if (T > nrow(design$Y) || T < 1L) stop("T is outside the generated path")
  Y <- design$Y[seq_len(T), , drop = FALSE]
  Y <- sweep(Y, 2L, colMeans(Y), "-")
  out <- list(Y = Y, T = as.integer(T))
  if (make_S) out$S <- crossprod(Y) / T
  out
}

# Keep all thin-QR rows of estimated loadings. Do not truncate an estimate to
# the population rank 3 or 4: doing so changes its fitted matrix subspace.
normalized_qr_cores <- function(A, B) {
  qa <- qr(A)
  qb <- qr(B)
  # Undo any rank-revealing column pivot so component labels are preserved.
  Ra <- qr.R(qa)[, order(qa$pivot), drop = FALSE]
  Rb <- qr.R(qb)[, order(qb$pivot), drop = FALSE]
  list(Ar = normalize_columns(Ra), Br = normalize_columns(Rb),
       QA = qr.Q(qa), QB = qr.Q(qb),
       A_numerical_rank = qa$rank, B_numerical_rank = qb$rank)
}

validate_R_design <- function() {
  bd <- generate_data("boundary", 1L, Tmax = 800L)
  ex <- generate_data("exterior", 1L, Tmax = 800L)
  it <- generate_data("interior", 1L, Tmax = 800L)
  bd_short <- generate_data("boundary", 1L, Tmax = 100L)
  stopifnot(identical(bd$X, ex$X), identical(bd$noise, ex$noise),
            identical(bd$QA, ex$QA), identical(bd$QB, ex$QB),
            nrow(bd$Y) == 800L, ncol(bd$Y) == 1600L,
            max(abs(colSums(it$A^2) - 40)) < 1e-10,
            max(abs(colSums(it$B^2) - 40)) < 1e-10,
            qr(it$A)$rank == 4L, qr(bd$A)$rank == 3L,
            qr(it$Z)$rank == 8L, qr(bd$Z)$rank == 5L,
            ex$extra_max_loading_cosine < 0.9,
            max(abs(bd$Y[seq_len(100L), ] - bd_short$Y)) < 1e-12)
  cores <- normalized_qr_cores(it$A, it$B)
  stopifnot(max(abs(crossprod(cores$Ar) - crossprod(it$A_unit))) < 1e-10)
  invisible(TRUE)
}
