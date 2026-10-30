# ================================================================
# CP-QMLE grid simulation, 2026-10-27
# p,q in {20,40,80}; T=100,200,400,800; 100 seeds; three regimes.
# Iteration caps per start: interior=300, boundary=400, exterior=500.
# Uses simulation truth for optimizer starts; this is an oracle-start experiment.
# DGP noise covariance=I. QMLE still ESTIMATES diagonal noise variances.
# Geometry uses jointly matched, unitary-removed, unit-column Ar and Br.
# J orig/nor are UNSQUARED singular values. Gram eigenvalues are separate.
# Every local candidate and its position is saved; no fake-zero classification.
# ================================================================
if(.Platform$OS.type=="windows" && !isTRUE(l10n_info()[["UTF-8"]]))
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE)))
.cp_script_dir<-local({
  paths<-lapply(sys.frames(),function(frame)frame$ofile)
  paths<-Filter(function(x)is.character(x)&&length(x)==1L&&nzchar(x)&&file.exists(x),paths)
  if(length(paths))dirname(normalizePath(paths[[length(paths)]],winslash="/",mustWork=FALSE))
  else {
    arg<-grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
    if(length(arg))dirname(normalizePath(sub("^--file=","",arg[1L]),winslash="/",mustWork=FALSE))
    else normalizePath(getwd(),winslash="/",mustWork=FALSE)
  }
})
.cp_r_minor<-strsplit(R.version$minor,".",fixed=TRUE)[[1L]][1L]
.cp_default_dependency_lib<-file.path(.cp_script_dir,"CP_QMLE_R_library",
  paste(R.version$major,.cp_r_minor,sep="."),R.version$platform)
.cp_grid_cpu<-suppressWarnings(parallel::detectCores(logical=FALSE))
if(length(.cp_grid_cpu)!=1L||!is.finite(.cp_grid_cpu)).cp_grid_cpu<-2L
.cp_grid_workers<-max(1L,min(4L,as.integer(.cp_grid_cpu)-1L))

# -------- USER CONFIGURATION --------
GRID_CONFIG<-list(
  p_values=c(20L,40L,80L),
  q_values=c(20L,40L,80L),
  seeds=1:100,
  sample_sizes=c(100L,200L,400L,800L),
  regimes=c("interior","boundary","exterior"),
  maxit_by_regime=c(interior=300L,boundary=400L,exterior=500L),
  truth_starts=6L,              # 1 exact truth + 5 near-truth starts
  truth_angle_deg=3,            # Max angle for each A/B loading perturbation
  truth_chol_sd=.05,
  truth_noise_log_sd=.05,
  factr=1000,
  pgtol=1e-8,
  lmm=15L,
  estimate_noise=TRUE,
  noise_bounds=c(.05,20),
  diagnostic_random_starts=30L, # Truth-direction starts are added separately
  diagnostic_maxeval=400L,
  diagnostic_cosine_cap=.90,
  diagnostic_cross_start_limit=5L,
  diagnostic_cross_rounds=2L,   # Recheck both objectives from each other's candidates
  alignment_angle_deg=5,
  alignment_cosine_caps=c(.95,.90),
  alignment_starts=5L,
  alignment_maxeval=300L,
  rank_tolerance=1e-10,
  n_exterior=.5,
  factor_multiplier=1,
  burn=100L,
  ncores=.cp_grid_workers,      # Set 2, 4, 6, ... to suit CPU and memory
  batch_jobs=2L*.cp_grid_workers,
  blas_threads_per_worker=1L,
  resume=TRUE,
  save_search_details=FALSE,    # All distinct local candidates are still saved
  save_all_start_fits=FALSE,    # All QMLE-start scalar logs are always saved
  report_detailed_summaries=TRUE,
  report_chunk_rows=20000L,
  validate_gradients=TRUE,
  verbose=TRUE,
  auto_install_dependencies=TRUE,
  dependency_repos="https://cloud.r-project.org",
  dependency_lib=.cp_default_dependency_lib,
  output_dir=file.path(.cp_script_dir,"CP_QMLE_grid_results")
)

# Dependency discovery and optional project-local installation.
# Defining this function performs no installation or other side effect.

cp_ensure_dependencies <- function(cfg) {
  # Keep a working package already visible in the current R session.
  if (requireNamespace("nloptr", quietly = TRUE)) return(invisible(TRUE))

  script_dir <- get0(".cp_script_dir", envir = environment(cp_ensure_dependencies),
                     inherits = TRUE, ifnotfound = NULL)
  if (!is.character(script_dir) || length(script_dir) != 1L ||
      is.na(script_dir) || !nzchar(script_dir)) {
    script_dir <- getwd()
  }
  script_dir <- normalizePath(path.expand(script_dir), winslash = "/", mustWork = FALSE)

  project_lib <- cfg$dependency_lib
  if (!is.character(project_lib) || length(project_lib) != 1L ||
      is.na(project_lib) || !nzchar(project_lib)) {
    version_minor <- strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L]
    version_tag <- paste(R.version$major, version_minor, sep = ".")
    project_lib <- file.path(script_dir, "CP_QMLE_R_library", version_tag,
                             R.version$platform)
  }
  project_lib <- normalizePath(path.expand(project_lib), winslash = "/", mustWork = FALSE)
  repos <- cfg$dependency_repos
  if (!is.character(repos) || !length(repos) || anyNA(repos) || any(!nzchar(repos))) {
    repos <- "https://cloud.r-project.org"
  }

  extra_lib <- Sys.getenv("CP_QMLE_R_LIB", unset = "")
  candidates <- c(project_lib, extra_lib,
                  file.path(script_dir, "R_lib"),
                  file.path(script_dir, "..", "work", "R_lib"),
                  file.path(getwd(), "R_lib"),
                  file.path(getwd(), "work", "R_lib"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- unique(vapply(candidates, function(path) {
    normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  }, character(1L)))
  existing <- candidates[dir.exists(candidates)]
  original_paths <- .libPaths()
  if (length(existing)) {
    # Try each candidate first. This also avoids an incompatible package in
    # an earlier library masking a usable build in a later project library.
    for (library_path in existing) {
      .libPaths(unique(c(library_path, existing, original_paths)))
      if (requireNamespace("nloptr", quietly = TRUE)) {
        message("Using nloptr from an existing R library: ",
                dirname(find.package("nloptr")))
        return(invisible(TRUE))
      }
    }
  }

  lib_literal <- encodeString(project_lib, quote = '"')
  repos_literal <- paste(deparse(repos), collapse = "")
  install_command <- paste0('install.packages("nloptr", lib = ', lib_literal,
                             ", repos = ", repos_literal, ")")
  library_command <- paste0(".libPaths(c(", lib_literal, ", .libPaths()))")
  instructions <- paste0("\nProject library: ", project_lib,
                          "\nRun in the R console:\n  ", install_command,
                          "\n  ", library_command)

  auto_install <- cfg$auto_install_dependencies
  if (is.null(auto_install)) auto_install <- TRUE
  if (!isTRUE(auto_install)) {
    stop("Package 'nloptr' is not available in this R session.", instructions,
         "\nThen rerun the script, or enable CONFIG$auto_install_dependencies.",
         call. = FALSE)
  }

  # .libPaths() silently ignores nonexistent folders, so create it first.
  if (!dir.exists(project_lib)) {
    dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(project_lib)) {
    stop("Could not create the project R library: ", project_lib,
         "\nSet CONFIG$dependency_lib to a writable project folder.",
         instructions, call. = FALSE)
  }
  .libPaths(unique(c(project_lib, existing, original_paths)))
  message("Installing the missing package 'nloptr' into the project library: ", project_lib)
  install_error <- NULL
  tryCatch(
    utils::install.packages("nloptr", lib = project_lib, repos = repos),
    error = function(e) install_error <<- conditionMessage(e)
  )
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    detail <- if (is.null(install_error)) "" else paste0("\nInstallation error: ", install_error)
    stop("The installation did not make package 'nloptr' available.", detail,
         "\nCheck the CRAN connection and whether this library is writable.",
         instructions, "\nAfter installation succeeds, rerun the script.",
         call. = FALSE)
  }
  invisible(TRUE)
}


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


# CP covariance QMLE backend. Only base R is required.
# Y has T rows and p*q columns; each row uses vec_F of a p by q matrix.
# The objective uses Woodbury algebra and never forms/inverts a (p*q)^2 matrix.

cp_normalize_columns <- function(A, norm = 1) {
  sweep(A, 2L, pmax(sqrt(colSums(A * A)), 1e-12) / norm, "/")
}

cp_khatri_rao <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  vapply(seq_len(ncol(A)), function(j) as.vector(outer(A[, j], B[, j])),
         numeric(nrow(A) * nrow(B)))
}

cp_qmle_objective <- function(Y, p, q, d, estimate_noise = TRUE,
                              fixed_noise = 1, noise_bounds = c(0.05, 20),
                              chol_diag_bounds = c(1e-3, 20),
                              chol_offdiag_bound = 20) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  N <- p * q
  stopifnot(ncol(Y) == N, nrow(Y) >= 2L, all(is.finite(Y)), d < N)
  if (length(fixed_noise) == 1L) fixed_noise <- rep(fixed_noise, N)
  stopifnot(length(fixed_noise) == N, all(fixed_noise > 0))
  Y <- sweep(Y, 2L, colMeans(Y), "-")
  TT <- nrow(Y)
  y2 <- colMeans(Y * Y)
  nA <- p * d
  nB <- q * d
  lower_mask <- lower.tri(matrix(0, d, d), diag = TRUE)
  lower_rc <- which(lower_mask, arr.ind = TRUE)
  diagonal_positions <- which(lower_rc[, 1L] == lower_rc[, 2L])
  nL <- nrow(lower_rc)
  loading_end <- nA + nB
  chol_end <- loading_end + nL
  lower <- c(rep(-Inf, loading_end), rep(-chol_offdiag_bound, nL))
  upper <- c(rep( Inf, loading_end), rep( chol_offdiag_bound, nL))
  lower[loading_end + diagonal_positions] <- log(chol_diag_bounds[1L])
  upper[loading_end + diagonal_positions] <- log(chol_diag_bounds[2L])
  if (estimate_noise) {
    lower <- c(lower, rep(log(noise_bounds[1L]), N))
    upper <- c(upper, rep(log(noise_bounds[2L]), N))
  }

  unpack <- function(par) {
    raw_A <- matrix(par[seq_len(nA)], p, d)
    raw_B <- matrix(par[nA + seq_len(nB)], q, d)
    norms_A <- pmax(sqrt(colSums(raw_A * raw_A)), 1e-12)
    norms_B <- pmax(sqrt(colSums(raw_B * raw_B)), 1e-12)
    unit_A <- sweep(raw_A, 2L, norms_A, "/")
    unit_B <- sweep(raw_B, 2L, norms_B, "/")
    L <- matrix(0, d, d)
    lv <- par[loading_end + seq_len(nL)]
    L[lower_mask] <- lv
    diag(L) <- exp(lv[diagonal_positions])
    psi <- if (estimate_noise) exp(par[chol_end + seq_len(N)]) else fixed_noise
    list(A = sqrt(p) * unit_A, B = sqrt(q) * unit_B,
         unit_A = unit_A, unit_B = unit_B,
         norms_A = norms_A, norms_B = norms_B, L = L, psi = psi)
  }

  pack <- function(A, B, M, psi = fixed_noise) {
    ee <- eigen((M + t(M)) / 2, symmetric = TRUE)
    M <- tcrossprod(sweep(ee$vectors, 2L, sqrt(pmax(ee$values, 1e-6)), "*"))
    L <- t(chol(M))
    lv <- L[lower_mask]
    lv[diagonal_positions] <- log(diag(L))
    par <- c(as.vector(cp_normalize_columns(A)),
             as.vector(cp_normalize_columns(B)), lv)
    if (estimate_noise) par <- c(par, log(pmax(noise_bounds[1L], pmin(noise_bounds[2L], psi))))
    pmax(lower, pmin(upper, par))
  }

  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$result <- NULL
  evaluate <- function(par) {
    if (identical(par, cache$par)) return(cache$result)
    parts <- unpack(par)
    A <- parts$A
    B <- parts$B
    L <- parts$L
    psi <- parts$psi
    Z <- cp_khatri_rao(A, B)
    F <- Z %*% L
    DiF <- sweep(F, 1L, psi, "/")
    K <- diag(d) + crossprod(F, DiF)
    ck <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(ck)) {
      result <- list(value = 1e30, gradient = rep(0, length(par)))
    } else {
      Ki <- chol2inv(ck)
      Q <- Y %*% DiF
      YV <- Q %*% Ki
      V <- DiF %*% Ki                         # Sigma^{-1} F
      trace_SW <- sum(y2 / psi) - sum(Q * YV) / TT
      value <- (sum(log(psi)) + 2 * sum(log(diag(ck))) + trace_SW) / (2 * N)
      X <- crossprod(Y, YV) / TT
      grad_F <- (V - sweep(X, 1L, psi, "/") +
                   V %*% (crossprod(Q, YV) / TT)) / N
      grad_Z <- grad_F %*% t(L)
      grad_A <- matrix(0, p, d)
      grad_B <- matrix(0, q, d)
      for (j in seq_len(d)) {
        gj <- matrix(grad_Z[, j], p, q)
        grad_A[, j] <- gj %*% B[, j]
        grad_B[, j] <- t(gj) %*% A[, j]
      }
      grad_raw_A <- sweep(grad_A - sweep(parts$unit_A, 2L,
                                  colSums(parts$unit_A * grad_A), "*"),
                          2L, parts$norms_A / sqrt(p), "/")
      grad_raw_B <- sweep(grad_B - sweep(parts$unit_B, 2L,
                                  colSums(parts$unit_B * grad_B), "*"),
                          2L, parts$norms_B / sqrt(q), "/")
      grad_L <- crossprod(Z, grad_F)
      glv <- grad_L[lower_mask]
      glv[diagonal_positions] <- glv[diagonal_positions] * diag(L)
      gradient <- c(as.vector(grad_raw_A), as.vector(grad_raw_B), glv)
      if (estimate_noise) {
        diag_W <- 1 / psi - rowSums(V * DiF)
        YW <- sweep(Y, 2L, psi, "/") - YV %*% t(DiF)
        diag_WSW <- colMeans(YW * YW)
        gradient <- c(gradient, (diag_W - diag_WSW) * psi / (2 * N))
      }
      result <- list(value = as.numeric(value), gradient = as.vector(gradient))
    }
    cache$par <- par
    cache$result <- result
    result
  }
  projected_gradient <- function(par, gradient) {
    gradient[(par <= lower + 1e-8) & gradient > 0] <- 0
    gradient[(par >= upper - 1e-8) & gradient < 0] <- 0
    gradient
  }
  list(fn = function(par) evaluate(par)$value,
       gr = function(par) evaluate(par)$gradient,
       evaluate = evaluate, unpack = unpack, pack = pack,
       projected_gradient = projected_gradient, lower = lower, upper = upper,
       Y = Y, p = p, q = q, d = d, N = N, y2 = y2,
       estimate_noise = estimate_noise, fixed_noise = fixed_noise)
}

cp_principal_subspace <- function(Y, d, method = "auto", oversample = 8L,
                                  power_iterations = 2L) {
  # The deterministic seed is supplied by cp_qmle before calling this helper.
  if (method %in% c("auto", "irlba") && requireNamespace("irlba", quietly = TRUE)) {
    fit <- tryCatch(irlba::irlba(Y, nv = d, nu = 0, tol = 1e-7), error = function(e) NULL)
    if (!is.null(fit)) return(list(U = fit$v, values = fit$d, method = "irlba"))
  }
  if (method %in% c("auto", "RSpectra") && requireNamespace("RSpectra", quietly = TRUE)) {
    fit <- tryCatch(RSpectra::svds(Y, k = d, nu = 0, nv = d,
                                  opts = list(tol = 1e-7)), error = function(e) NULL)
    if (!is.null(fit)) return(list(U = fit$v, values = fit$d, method = "RSpectra"))
  }
  if (method == "base_svd") {
    fit <- svd(Y, nu = 0, nv = d)
    return(list(U = fit$v[, seq_len(d), drop = FALSE], values = fit$d[seq_len(d)], method = "base_svd"))
  }
  ell <- min(d + oversample, nrow(Y) - 1L, ncol(Y))
  if (ell < d) stop("Too few observations for the requested factor dimension")
  Q <- qr.Q(qr(crossprod(Y, matrix(rnorm(nrow(Y) * ell), nrow(Y), ell))))
  for (iteration in seq_len(power_iterations)) Q <- qr.Q(qr(crossprod(Y, Y %*% Q)))
  small <- svd(Y %*% Q, nu = 0, nv = d)
  list(U = Q %*% small$v[, seq_len(d), drop = FALSE],
       values = small$d[seq_len(d)], method = "randomized_svd_2_power_iterations")
}

cp_rankone_candidate_pool <- function(U, p, q, count = 4L * ncol(U),
                                      alternating_steps = 20L) {
  d <- ncol(U)
  blocks <- lapply(seq_len(d), function(j) matrix(U[, j], p, q))
  candidates <- list()
  for (k in seq_len(count)) {
    co <- if (k <= d) diag(d)[, k] else rnorm(d)
    initial <- matrix(U %*% co, p, q)
    ss <- svd(initial, nu = 1, nv = 1)
    a <- ss$u[, 1L]
    b <- ss$v[, 1L]
    for (iteration in seq_len(alternating_steps)) {
      GG <- vapply(blocks, function(M) as.vector(M %*% b), numeric(p))
      a <- svd(GG, nu = 1, nv = 0)$u[, 1L]
      HH <- vapply(blocks, function(M) as.vector(crossprod(M, a)), numeric(q))
      b <- svd(HH, nu = 1, nv = 0)$u[, 1L]
    }
    z <- as.vector(outer(a, b))
    duplicate <- length(candidates) > 0L && any(vapply(candidates,
                         function(old) abs(sum(z * old$z)) > 0.9999, logical(1)))
    if (!duplicate) candidates[[length(candidates) + 1L]] <-
      list(a = a, b = b, z = z, residual = max(0, 1 - sum(crossprod(U, z)^2)))
  }
  candidates
}

cp_qmle <- function(Y, p, q, d, seed = 1L, n_starts = 3L, maxit = 300L,
                    estimate_noise = TRUE, fixed_noise = 1,
                    previous_fit = NULL, noise_bounds = c(0.05, 20),
                    pca_method = "auto", factr = 1e7, pgtol = 1e-6) {
  # Each start has exactly the supplied maxit cap; no hidden extra polishing run.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  n_starts <- as.integer(n_starts)
  stopifnot(n_starts >= 1L, maxit >= 1L)
  began <- proc.time()[["elapsed"]]
  obj <- cp_qmle_objective(Y, p, q, d, estimate_noise, fixed_noise, noise_bounds)
  pca <- cp_principal_subspace(obj$Y, d, method = pca_method)
  pool <- cp_rankone_candidate_pool(pca$U, p, q, count = max(4L*d, 24L))
  initial_noise <- max(0.1, min(5, (sum(obj$y2) - sum(pca$values^2) / nrow(obj$Y)) / (p*q-d)))
  starts <- list()
  sources <- character()
  if (!is.null(previous_fit)) {
    starts[[1L]] <- obj$pack(previous_fit$A, previous_fit$B, previous_fit$M, previous_fit$psi)
    sources <- "previous_T"
  }
  while (length(starts) < n_starts) {
    s <- length(starts) + 1L
    if (length(pool) >= d) {
      Zpool <- vapply(pool, function(x) x$z, numeric(p*q))
      residual <- vapply(pool, function(x) x$residual, numeric(1))
      weight <- exp(-3 * residual)
      if (s > 1L) weight <- weight * exp(rnorm(length(weight), 0, 0.35))
      chosen <- qr(sweep(Zpool, 2L, weight, "*"), LAPACK = TRUE)$pivot[seq_len(d)]
      A <- sqrt(p) * vapply(pool[chosen], function(x) x$a, numeric(p))
      B <- sqrt(q) * vapply(pool[chosen], function(x) x$b, numeric(q))
    } else {
      A <- cp_normalize_columns(matrix(rnorm(p*d), p, d), sqrt(p))
      B <- cp_normalize_columns(matrix(rnorm(q*d), q, d), sqrt(q))
    }
    if (s > 1L) {
      A <- cp_normalize_columns(A + matrix(rnorm(p*d, 0, 0.05), p, d), sqrt(p))
      B <- cp_normalize_columns(B + matrix(rnorm(q*d, 0, 0.05), q, d), sqrt(q))
    }
    Z <- cp_khatri_rao(A, B)
    gram <- crossprod(Z)
    Zi <- solve(gram + diag(max(diag(gram)) * 1e-8, d), t(Z))
    psi <- if (estimate_noise) rep(initial_noise, p*q) else obj$fixed_noise
    scores <- obj$Y %*% t(Zi)
    M <- crossprod(scores) / nrow(scores) - tcrossprod(sweep(Zi, 2L, sqrt(psi), "*"))
    em <- eigen((M + t(M))/2, symmetric = TRUE)
    M <- tcrossprod(sweep(em$vectors, 2L, sqrt(pmax(0.01, pmin(20, em$values))), "*"))
    starts[[s]] <- obj$pack(A, B, M, psi)
    sources[s] <- "sample_PCA_rankone"
  }
  fits <- vector("list", n_starts)
  logs <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    tick <- proc.time()[["elapsed"]]
    fit <- optim(starts[[s]], obj$fn, obj$gr, method = "L-BFGS-B",
                 lower = obj$lower, upper = obj$upper,
                 control = list(maxit = maxit, factr = factr, pgtol = pgtol, lmm = 15L))
    ev <- obj$evaluate(fit$par)
    pg <- obj$projected_gradient(fit$par, ev$gradient)
    logs[[s]] <- data.frame(start = s, source = sources[s], nll = ev$value,
      convergence = fit$convergence, converged = fit$convergence == 0L,
      function_evaluations = unname(fit$counts["function"]),
      gradient_evaluations = unname(fit$counts["gradient"]),
      gradient_norm = max(abs(ev$gradient)), projected_gradient_norm = max(abs(pg)),
      elapsed_seconds = proc.time()[["elapsed"]] - tick,
      message = if (is.null(fit$message)) "" else fit$message,
      stringsAsFactors = FALSE)
    fits[[s]] <- fit
  }
  logs <- do.call(rbind, logs)
  best <- which.min(logs$nll)
  selected <- fits[[best]]
  parts <- obj$unpack(selected$par)
  list(A = parts$A, B = parts$B, M = tcrossprod(parts$L), psi = parts$psi,
       nll = logs$nll[best], converged = logs$converged[best],
       convergence = logs$convergence[best], message = logs$message[best],
       gradient_norm = logs$gradient_norm[best],
       projected_gradient_norm = logs$projected_gradient_norm[best],
       best_start = best, all_starts = logs, maxit_per_start = maxit,
       n_starts = n_starts, pca_method = pca$method,
       noise_mode = if (estimate_noise) "estimated_diagonal" else "known_diagonal",
       elapsed_seconds = proc.time()[["elapsed"]] - began,
       noise_at_lower_bound = if (estimate_noise) sum(parts$psi < noise_bounds[1L]*(1+1e-4)) else 0L,
       noise_at_upper_bound = if (estimate_noise) sum(parts$psi > noise_bounds[2L]*(1-1e-4)) else 0L)
}

cp_check_qmle_gradient <- function(seed = 101L, estimate_noise = TRUE) {
  set.seed(seed)
  p <- 3L; q <- 4L; d <- 5L
  Y <- matrix(rnorm(120L * p*q), 120L, p*q)
  obj <- cp_qmle_objective(Y, p, q, d, estimate_noise, fixed_noise = 1)
  A <- matrix(rnorm(p*d), p, d)
  B <- matrix(rnorm(q*d), q, d)
  LL <- diag(d) + matrix(rnorm(d*d, 0, 0.05), d, d) * lower.tri(matrix(0,d,d))
  par <- obj$pack(A, B, tcrossprod(LL), rep(1.2,p*q))
  analytical <- obj$gr(par)
  numerical <- numeric(length(par))
  for (j in seq_along(par)) {
    h <- 1e-5 * max(1, abs(par[j]))
    plus <- minus <- par
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    numerical[j] <- (obj$fn(plus)-obj$fn(minus))/(2*h)
  }
  absolute <- max(abs(analytical - numerical))
  relative <- sqrt(sum((analytical - numerical)^2)) / max(1e-12, sqrt(sum(analytical^2)))
  stopifnot(absolute < 2e-7, relative < 2e-6)
  # Independently compare Woodbury loss with direct dense covariance evaluation.
  pieces <- obj$unpack(par)
  Z <- cp_khatri_rao(pieces$A, pieces$B)
  Sigma <- Z %*% tcrossprod(pieces$L) %*% t(Z) + diag(pieces$psi)
  S <- crossprod(obj$Y)/nrow(obj$Y)
  direct <- (as.numeric(determinant(Sigma, logarithm=TRUE)$modulus) +
               sum(S * solve(Sigma)))/(2*p*q)
  stopifnot(abs(direct - obj$fn(par)) < 1e-10)
  list(estimate_noise = estimate_noise, max_absolute_error = absolute,
       relative_error = relative, dense_loss_difference = direct-obj$fn(par))
}


# Core geometry and Proposition 12 diagnostics for CP QMLE simulations.
#
# Each fitted loading is first jointly matched to the population columns
# (when population loadings are supplied), then reduced by economy QR.
# The default keeps min(nrow(A), d) and min(nrow(B), d) rows: estimated
# loading ranks are NOT forced to equal the population loading ranks.
# Columns of the resulting cores are normalized to unit length.
#
# rho and J use each independent 2 x 2 minor once.  The complete ordered
# minor vector has twice the Euclidean norm. J_orig and J_nor contain
# singular values directly, sorted increasingly; neither is squared.
# For each ell, nor = orig / max(orig), before averaging across seeds.
#
# Required package for constrained diagnostic searches: nloptr.
# Optional package for joint assignment: clue; a base-R exact dynamic
# programming fallback is included for d <= 16.
# No installation or filesystem mutation is performed by this module.

cp_assert_loadings <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  if (!is.numeric(A) || !is.numeric(B) || any(!is.finite(A)) ||
      any(!is.finite(B)) || ncol(A) != ncol(B) || ncol(A) < 2L) {
    stop("A and B must be finite numeric matrices with the same d >= 2 columns.")
  }
  if (any(colSums(A * A) <= 0) || any(colSums(B * B) <= 0)) {
    stop("Every loading column must be nonzero.")
  }
  list(A = A, B = B)
}

cp_unit_columns <- function(A) {
  norms <- sqrt(colSums(A * A))
  if (any(!is.finite(norms) | norms <= 0)) stop("Cannot normalize zero columns.")
  sweep(A, 2L, norms, "/")
}

cp_max_assignment <- function(score) {
  d <- nrow(score)
  if (d != ncol(score) || any(!is.finite(score))) stop("Invalid assignment score.")
  if (requireNamespace("clue", quietly = TRUE)) {
    return(as.integer(clue::solve_LSAP(score, maximum = TRUE)))
  }
  if (d > 16L) stop("Install package 'clue' for assignment with d > 16.")
  # dp[mask] is the best assignment of the first popcount(mask) rows to
  # columns selected by mask.  Returns the population column for each
  # estimated column, exactly as solve_LSAP would.
  nstates <- bitwShiftL(1L, d)
  dp <- rep(-Inf, nstates); dp[1L] <- 0
  previous <- integer(nstates); chosen <- integer(nstates)
  popcount <- integer(nstates)
  if (nstates > 1L) for (mask in seq_len(nstates - 1L)) {
    popcount[mask + 1L] <- popcount[bitwAnd(mask, mask - 1L) + 1L] + 1L
  }
  for (mask in 0L:(nstates - 2L)) {
    row <- popcount[mask + 1L] + 1L
    for (column in seq_len(d)) {
      bit <- bitwShiftL(1L, column - 1L)
      if (bitwAnd(mask, bit) != 0L) next
      next_mask <- bitwOr(mask, bit)
      value <- dp[mask + 1L] + score[row, column]
      if (value > dp[next_mask + 1L]) {
        dp[next_mask + 1L] <- value
        previous[next_mask + 1L] <- mask
        chosen[next_mask + 1L] <- column
      }
    }
  }
  assignment <- integer(d); mask <- nstates - 1L
  for (row in d:1L) {
    assignment[row] <- chosen[mask + 1L]
    mask <- previous[mask + 1L]
  }
  assignment
}

cp_joint_alignment <- function(A, B, A_true, B_true) {
  inputs <- cp_assert_loadings(A, B); A <- inputs$A; B <- inputs$B
  truth <- cp_assert_loadings(A_true, B_true)
  if (!identical(dim(A), dim(truth$A)) || !identical(dim(B), dim(truth$B))) {
    stop("Estimated and population loading dimensions must agree for matching.")
  }
  ca <- crossprod(cp_unit_columns(A), cp_unit_columns(truth$A))
  cb <- crossprod(cp_unit_columns(B), cp_unit_columns(truth$B))
  # This is absolute cosine agreement of rank-one loading matrices.
  # A and B always receive the SAME component permutation.
  score <- abs(ca * cb)
  estimated_to_truth <- cp_max_assignment(score)
  estimated_for_truth <- order(estimated_to_truth)
  A <- A[, estimated_for_truth, drop = FALSE]
  B <- B[, estimated_for_truth, drop = FALSE]
  signs_A <- sign(colSums(A * truth$A)); signs_A[signs_A == 0] <- 1
  signs_B <- sign(colSums(B * truth$B)); signs_B[signs_B == 0] <- 1
  list(A = sweep(A, 2L, signs_A, "*"), B = sweep(B, 2L, signs_B, "*"),
       estimated_for_truth = estimated_for_truth,
       estimated_to_truth = estimated_to_truth,
       signs_A = signs_A, signs_B = signs_B,
       matched_rankone_cosines = score[cbind(estimated_for_truth, seq_len(ncol(A)))],
       matching_rule = "maximize sum(abs(cos(Ahat_i,A_j)*cos(Bhat_i,B_j)))")
}

cp_qr_core <- function(A, retained_rows = NULL, rank_tolerance = 1e-10) {
  A <- as.matrix(A); d <- ncol(A)
  factor <- base::qr(A, tol = rank_tolerance, LAPACK = FALSE)
  Q <- base::qr.Q(factor, complete = FALSE)
  R <- base::qr.R(factor, complete = FALSE)
  # qr() may pivot when columns are numerically dependent.  Undo column
  # pivoting so labels continue to agree with the jointly matched columns.
  R <- R[, order(factor$pivot), drop = FALSE]
  available <- min(nrow(A), d)
  if (is.null(retained_rows)) retained_rows <- available
  retained_rows <- as.integer(retained_rows)
  if (length(retained_rows) != 1L || retained_rows < 1L || retained_rows > available) {
    stop("Configured retained_rows is outside the economy-QR dimensions.")
  }
  singular_values <- base::svd(A, nu = 0L, nv = 0L)$d
  numerical_rank <- sum(singular_values > max(singular_values) * rank_tolerance)
  core <- R[seq_len(retained_rows), , drop = FALSE]
  list(core = core, Q = Q[, seq_len(retained_rows), drop = FALSE],
       numerical_rank = numerical_rank, retained_rows = retained_rows,
       available_rows = available, configured_truncation = retained_rows < available,
       reconstruction_relative_error = sqrt(sum((A - Q %*% R)^2) / sum(A^2)))
}

cp_prepare_core <- function(A, B, A_true = NULL, B_true = NULL,
                            core_ranks = NULL, rank_tolerance = 1e-10) {
  inputs <- cp_assert_loadings(A, B); A <- inputs$A; B <- inputs$B
  matching <- NULL
  if (xor(is.null(A_true), is.null(B_true))) stop("Supply both population loadings or neither.")
  if (!is.null(A_true)) {
    matching <- cp_joint_alignment(A, B, A_true, B_true)
    A <- matching$A; B <- matching$B
    matching$A <- NULL; matching$B <- NULL
  }
  if (!is.null(core_ranks) && length(core_ranks) != 2L) {
    stop("core_ranks must be NULL or an explicit two-element truncation configuration.")
  }
  qa <- cp_qr_core(A, if (is.null(core_ranks)) NULL else core_ranks[1L], rank_tolerance)
  qb <- cp_qr_core(B, if (is.null(core_ranks)) NULL else core_ranks[2L], rank_tolerance)
  list(A = cp_unit_columns(qa$core), B = cp_unit_columns(qb$core),
       original_p = nrow(A), original_q = nrow(B), d = ncol(A),
       core_p = nrow(qa$core), core_q = nrow(qb$core),
       numerical_rank_A = qa$numerical_rank, numerical_rank_B = qb$numerical_rank,
       truncated_A = qa$configured_truncation, truncated_B = qb$configured_truncation,
       qr_reconstruction_error_A = qa$reconstruction_relative_error,
       qr_reconstruction_error_B = qb$reconstruction_relative_error,
       matching = matching)
}

cp_minor_matrix <- function(A, B) {
  d <- ncol(A)
  pairs <- utils::combn(seq_len(d), 2L)
  ap <- utils::combn(seq_len(nrow(A)), 2L)
  bp <- utils::combn(seq_len(nrow(B)), 2L)
  CA <- A[ap[1L, ], pairs[1L, ], drop = FALSE] * A[ap[2L, ], pairs[2L, ], drop = FALSE] -
        A[ap[1L, ], pairs[2L, ], drop = FALSE] * A[ap[2L, ], pairs[1L, ], drop = FALSE]
  CB <- B[bp[1L, ], pairs[1L, ], drop = FALSE] * B[bp[2L, ], pairs[2L, ], drop = FALSE] -
        B[bp[1L, ], pairs[2L, ], drop = FALSE] * B[bp[2L, ], pairs[1L, ], drop = FALSE]
  C <- vapply(seq_len(ncol(pairs)), function(j) {
    as.vector(kronecker(CA[, j], CB[, j]))
  }, numeric(nrow(CA) * nrow(CB)))
  list(C = C, pairs = pairs)
}

cp_jacobian_spectra <- function(A, B, C = NULL, pairs = NULL,
                                numerical_clip_tolerance = 1e-12) {
  A <- cp_unit_columns(A); B <- cp_unit_columns(B)
  d <- ncol(A); GA <- crossprod(A); GB <- crossprod(B)
  values <- matrix(NA_real_, d, d - 1L)
  relative <- values; grams <- vector("list", d)
  negative_roundoff <- numeric(d)
  if (is.null(C) || is.null(pairs)) {
    minors <- cp_minor_matrix(A, B)
    C <- minors$C; pairs <- minors$pairs
  }
  for (ell in seq_len(d)) {
    others <- setdiff(seq_len(d), ell)
    gram_A <- GA[ell, ell] * GA[others, others, drop = FALSE] -
              outer(GA[ell, others], GA[ell, others])
    gram_B <- GB[ell, ell] * GB[others, others, drop = FALSE] -
              outer(GB[ell, others], GB[ell, others])
    gram <- gram_A * gram_B
    gram <- (gram + t(gram)) / 2
    gram_ev <- sort(base::eigen(gram, symmetric = TRUE, only.values = TRUE)$values)
    negative_roundoff[ell] <- min(0, min(gram_ev))
    if (min(gram_ev) < -numerical_clip_tolerance * max(1, max(abs(gram_ev)))) {
      warning("A Jacobian Gram eigenvalue is materially negative; check numerical conditioning.")
    }
    # Report the SVD values directly, with no squaring. The Gram above is
    # retained only as an independently constructed algebraic check.
    incident <- which(colSums(pairs == ell) > 0)
    singular_values <- base::svd(C[, incident, drop = FALSE], nu = 0L, nv = 0L)$d
    sv <- sort(c(singular_values, rep(0, d - 1L - length(singular_values))))
    values[ell, ] <- sv
    relative[ell, ] <- if (max(sv) > 0) sv / max(sv) else rep(NA_real_, d - 1L)
    grams[[ell]] <- gram
  }
  dimnames(values) <- list(paste0("component_", seq_len(d)), paste0("svd_", seq_len(d - 1L)))
  dimnames(relative) <- dimnames(values)
  list(orig = values, nor = relative,
       gram_matrices = grams,
       J_min = min(values),
       J_relative_min = min(relative, na.rm = TRUE),
       negative_roundoff = negative_roundoff,
       svd_convention = "ascending svd(J_ell); orig=s; nor=s/max(s); not squared")
}

cp_diagnostic_problem <- function(core, cosine_cap = .90, rank_tolerance = 1e-10) {
  A <- core$A; B <- core$B; d <- ncol(A)
  minors <- cp_minor_matrix(A, B)
  Z <- vapply(seq_len(d), function(j) as.vector(outer(A[, j], B[, j])),
              numeric(nrow(A) * nrow(B)))
  singular_values <- base::svd(Z, nu = 0L, nv = 0L)$d
  if (min(singular_values) <= max(singular_values) * rank_tolerance) {
    stop("B odot A is numerically column-rank deficient; diagnostic coordinates are not invertible.")
  }
  factor <- base::qr(Z, LAPACK = FALSE)
  Q <- base::qr.Q(factor, complete = FALSE)
  R <- base::qr.R(factor, complete = FALSE)[, order(factor$pivot), drop = FALSE]
  Rinv <- solve(R); G <- crossprod(Z, Q)
  C <- minors$C; ii <- minors$pairs[1L, ]; jj <- minors$pairs[2L, ]
  rho_squared <- function(x) {
    lam <- drop(Rinv %*% x); k <- lam[ii] * lam[jj]
    den <- sum(k * k)
    if (!is.finite(den) || den < 1e-28) return(list(value = 1e10, gradient = rep(0, d)))
    h <- drop(C %*% k); num <- sum(h * h)
    gk <- 2 * (drop(crossprod(C, h)) * den - k * num) / den^2
    gl <- numeric(d)
    for (j in seq_along(ii)) {
      gl[ii[j]] <- gl[ii[j]] + gk[j] * lam[jj[j]]
      gl[jj[j]] <- gl[jj[j]] + gk[j] * lam[ii[j]]
    }
    list(value = num / den, gradient = drop(crossprod(Rinv, gl)))
  }
  E_squared <- function(x) {
    matrix_value <- matrix(Q %*% x, nrow(A), nrow(B))
    sv <- base::svd(matrix_value, nu = 1L, nv = 1L)
    norm2 <- sum(x * x); top2 <- sv$d[1L]^2
    rankone_vec <- as.vector(outer(sv$u[, 1L], sv$v[, 1L]))
    gradtop <- 2 * sv$d[1L] * drop(crossprod(Q, rankone_vec))
    gradient <- -(gradtop * norm2 - top2 * 2 * x) / norm2^2
    list(value = sum(sv$d[-1L]^2) / norm2, gradient = gradient)
  }
  g_ineq <- function(x) drop(G %*% x)^2 - cosine_cap^2 * sum(x * x)
  jac_ineq <- function(x) 2 * sweep(G, 1L, drop(G %*% x), "*") -
                         matrix(2 * cosine_cap^2 * x, d, d, byrow = TRUE)
  summarize <- function(x) {
    x <- x / sqrt(sum(x * x))
    list(x = x, lambda = drop(Rinv %*% x),
         rho = sqrt(max(0, rho_squared(x)$value)),
         E1 = sqrt(max(0, E_squared(x)$value)),
         max_cosine = max(abs(G %*% x)), matrix_norm = sqrt(sum(x * x)))
  }
  list(A = A, B = B, C = C, pairs = minors$pairs, Z = Z, Q = Q, R = R,
       Rinv = Rinv, G = G, d = d, rho_squared = rho_squared,
       E_squared = E_squared, g_ineq = g_ineq, jac_ineq = jac_ineq,
       summarize = summarize, cosine_cap = cosine_cap,
       Z_singular_values = singular_values)
}

cp_cluster_minima <- function(candidates, objective,
                              cluster_cosine_tolerance = 1e-5) {
  if (length(candidates) == 0L) return(list())
  values <- vapply(candidates, function(x) x$objective_value, numeric(1L))
  candidates <- candidates[order(values)]
  clusters <- list()
  for (candidate in candidates) {
    similarities <- if (length(clusters)) vapply(clusters, function(x) {
      abs(sum(candidate$x * x$x))
    }, numeric(1L)) else numeric()
    hit <- which(similarities >= 1 - cluster_cosine_tolerance)
    if (!length(hit)) {
      candidate$hit_count <- 1L
      candidate$converged_hit_count <- as.integer(candidate$success)
      clusters[[length(clusters) + 1L]] <- candidate
    } else {
      j <- hit[which.max(similarities[hit])]
      clusters[[j]]$hit_count <- clusters[[j]]$hit_count + 1L
      clusters[[j]]$converged_hit_count <- clusters[[j]]$converged_hit_count + as.integer(candidate$success)
    }
  }
  for (j in seq_along(clusters)) {
    clusters[[j]]$minimum_id <- j
    clusters[[j]]$objective <- objective
  }
  clusters
}

cp_diagnostics <- function(A, B, A_true = NULL, B_true = NULL,
                           seed = 123L, nstarts = 30L, cosine_cap = .90,
                           core_ranks = NULL, rank_tolerance = 1e-10,
                           maxeval = 400L, xtol_rel = 1e-9,
                           ftol_abs = 1e-12, feasibility_tolerance = 2e-6,
                           cluster_cosine_tolerance = 1e-5,
                           cross_start_limit = 5L, extra_lambdas = NULL) {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("Package 'nloptr' is required: install.packages('nloptr') before running diagnostic searches.")
  }
  if (nstarts < 1L || cosine_cap <= 0 || cosine_cap >= 1) stop("Invalid search configuration.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(as.integer(seed))
  core <- cp_prepare_core(A, B, A_true, B_true, core_ranks, rank_tolerance)
  problem <- cp_diagnostic_problem(core, cosine_cap, rank_tolerance)
  spectra <- cp_jacobian_spectra(core$A, core$B, problem$C, problem$pairs)
  d <- problem$d
  starts <- list()
  if (!is.null(extra_lambdas)) {
    if (is.numeric(extra_lambdas)) extra_lambdas <- list(extra_lambdas)
    for (lam in extra_lambdas) {
      if (length(lam) != d || any(!is.finite(lam))) stop("Invalid extra lambda start.")
      x <- drop(problem$R %*% lam); x <- x / sqrt(sum(x * x))
      if (max(problem$g_ineq(x)) <= 0) starts[[length(starts) + 1L]] <- x
    }
  }
  attempts <- 0L
  while (length(starts) < nstarts && attempts < max(10000L, nstarts * 1000L)) {
    x <- stats::rnorm(d); x <- x / sqrt(sum(x * x))
    if (max(problem$g_ineq(x)) <= -1e-10) starts[[length(starts) + 1L]] <- x
    attempts <- attempts + 1L
  }
  if (length(starts) < nstarts) stop("Could not generate enough feasible starts.")
  optimize_start <- function(x0, objective, source) {
    fun <- if (objective == "rho") problem$rho_squared else problem$E_squared
    fit <- tryCatch(nloptr::nloptr(
      x0 = x0,
      eval_f = function(x) fun(x)$value,
      eval_grad_f = function(x) fun(x)$gradient,
      eval_g_ineq = problem$g_ineq, eval_jac_g_ineq = problem$jac_ineq,
      eval_g_eq = function(x) sum(x * x) - 1,
      eval_jac_g_eq = function(x) matrix(2 * x, nrow = 1L),
      opts = list(algorithm = "NLOPT_LD_SLSQP", maxeval = as.integer(maxeval),
                  xtol_rel = xtol_rel, ftol_abs = ftol_abs, print_level = 0L)
    ), error = function(e) list(solution = x0, status = -999L,
                                message = conditionMessage(e), iterations = 0L))
    raw_norm <- sqrt(sum(fit$solution * fit$solution))
    result <- problem$summarize(fit$solution)
    result$feasible <- is.finite(raw_norm) && abs(raw_norm - 1) <= feasibility_tolerance &&
                       result$max_cosine <= cosine_cap + feasibility_tolerance
    result$objective_value <- if (objective == "rho") result$rho else result$E1
    result$optimized_objective_squared <- result$objective_value^2
    result$success <- fit$status %in% 1:4
    result$status <- fit$status; result$message <- fit$message
    result$iterations <- if (is.null(fit$iterations)) NA_integer_ else fit$iterations
    result$source <- source
    result
  }
  rho_candidates <- lapply(seq_along(starts), function(i) {
    optimize_start(starts[[i]], "rho", paste0("random_or_extra_", i))
  })
  rho_valid <- Filter(function(x) x$feasible, rho_candidates)
  rho_clusters <- cp_cluster_minima(rho_valid, "rho", cluster_cosine_tolerance)
  E_starts <- starts
  if (length(rho_clusters) && cross_start_limit > 0L) {
    E_starts <- c(E_starts, lapply(utils::head(rho_clusters, cross_start_limit), function(x) x$x))
  }
  E_candidates <- lapply(seq_along(E_starts), function(i) {
    optimize_start(E_starts[[i]], "E1", paste0(if (i <= length(starts)) "random_or_extra_" else "rho_warm_", i))
  })
  E_valid <- Filter(function(x) x$feasible, E_candidates)
  E_clusters <- cp_cluster_minima(E_valid, "E1", cluster_cosine_tolerance)
  if (length(E_clusters) && cross_start_limit > 0L) {
    extra_rho <- lapply(utils::head(E_clusters, cross_start_limit), function(x) {
      optimize_start(x$x, "rho", "E1_warm")
    })
    rho_candidates <- c(rho_candidates, extra_rho)
    rho_valid <- Filter(function(x) x$feasible, rho_candidates)
    rho_clusters <- cp_cluster_minima(rho_valid, "rho", cluster_cosine_tolerance)
  }
  if (!length(rho_clusters) || !length(E_clusters)) stop("No feasible diagnostic minimum was returned.")
  best_rho <- rho_clusters[[1L]]; best_E <- E_clusters[[1L]]
  # Save lambda, not x, in delivered local-minimum records.  Lambda is in
  # matched unit-column core coordinates; reconstruction uses these cores.
  strip_x <- function(x) { x$x <- NULL; x }
  list(
    rho_min = best_rho$rho, E_at_rho = best_rho$E1,
    E_min = best_E$E1, rho_at_E = best_E$rho,
    lambda_rho = best_rho$lambda, lambda_E = best_E$lambda,
    cosine_at_rho = best_rho$max_cosine, cosine_at_E = best_E$max_cosine,
    rho_best_success = best_rho$success, E_best_success = best_E$success,
    rho_local_minima = lapply(rho_clusters, strip_x),
    E_local_minima = lapply(E_clusters, strip_x),
    n_rho_local_minima = length(rho_clusters), n_E_local_minima = length(E_clusters),
    rho_search_runs = length(rho_candidates), E_search_runs = length(E_candidates),
    rho_search_feasible = length(rho_valid), E_search_feasible = length(E_valid),
    rho_search_converged = sum(vapply(rho_candidates, function(x) x$success, logical(1L))),
    E_search_converged = sum(vapply(E_candidates, function(x) x$success, logical(1L))),
    J_min = spectra$J_min, J_relative_min = spectra$J_relative_min,
    J_orig = spectra$orig,
    J_nor = spectra$nor,
    J_negative_roundoff = spectra$negative_roundoff,
    Z_singular_values = problem$Z_singular_values,
    C_singular_values = base::svd(problem$C, nu = 0L, nv = 0L)$d,
    core = core,
    configuration = list(seed = seed, random_or_extra_starts = length(starts),
                         requested_nstarts = nstarts, cosine_cap = cosine_cap,
                         maxeval = maxeval, xtol_rel = xtol_rel, ftol_abs = ftol_abs,
                         feasibility_tolerance = feasibility_tolerance,
                         cluster_cosine_tolerance = cluster_cosine_tolerance,
                         cross_start_limit = cross_start_limit,
                         rank_tolerance = rank_tolerance,
                         minor_convention = "independent 2x2 minors counted once",
                         svd_convention = spectra$svd_convention,
                         minimum_note = "Distinct feasible numerical local minima found; not a certificate that every local minimum was found.")
  )
}


# Local E1 searches in projective lambda caps about numerical rho minima.
# Requires the functions in R_diagnostics.R and package nloptr.
# Angles below are Euclidean lambda angles, modulo a common sign.

cp_alignment_angle_deg <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  den <- sqrt(sum(a*a)*sum(b*b))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  value <- min(1, max(0, abs(sum(a*b))/den))
  acos(value)*180/pi
}

cp_alignment_change_cap <- function(problem, cosine_cap) {
  G <- problem$G; d <- problem$d; cap <- cosine_cap
  problem$g_ineq <- function(x) drop(G %*% x)^2 - cap^2*sum(x*x)
  problem$jac_ineq <- function(x) 2*sweep(G,1L,drop(G %*% x),"*") -
    matrix(2*cap^2*x,d,d,byrow=TRUE)
  problem$cosine_cap <- cap
  problem
}

cp_alignment_x_from_lambda <- function(problem, lambda) {
  x <- drop(problem$R %*% lambda)
  x/sqrt(sum(x*x))
}

cp_alignment_cap_constraints <- function(problem, lambda0, angle_deg) {
  unit0 <- lambda0/sqrt(sum(lambda0^2))
  cosine <- cos(angle_deg*pi/180)
  L <- problem$Rinv
  linear <- drop(crossprod(L,unit0))
  list(
    value=function(x) {
      lambda <- drop(L %*% x)
      c(problem$g_ineq(x), cosine*sqrt(sum(lambda^2))-sum(unit0*lambda))
    },
    jacobian=function(x) {
      lambda <- drop(L %*% x)
      norm <- sqrt(sum(lambda^2))
      if (!is.finite(norm) || norm<=0) stop("Degenerate lambda in angular constraint.")
      angular <- cosine*drop(crossprod(L,lambda))/norm-linear
      rbind(problem$jac_ineq(x),angular)
    },
    unit0=unit0,
    cosine=cosine
  )
}

cp_alignment_cap_search <- function(problem, lambda0, angle_deg=5,
                                    nstarts=5L,maxeval=300L,
                                    extra_lambdas=list(),
                                    feasibility_tolerance=2e-6,
                                    angle_tolerance_deg=1e-4,
                                    edge_tolerance_deg=.05,
                                    xtol_rel=1e-9,ftol_abs=1e-13) {
  stopifnot(angle_deg>0,angle_deg<90,nstarts>=1L,maxeval>=1L)
  constraints <- cp_alignment_cap_constraints(problem,lambda0,angle_deg)
  unit0 <- constraints$unit0
  center_x <- cp_alignment_x_from_lambda(problem,lambda0)
  summarize <- function(solution,status,message,source,iterations=NA_integer_) {
    raw_norm <- sqrt(sum(solution^2))
    point <- problem$summarize(solution)
    point$angle_deg <- cp_alignment_angle_deg(point$lambda,lambda0)
    point$matrix_angle_deg <- cp_alignment_angle_deg(point$x,center_x)
    point$eq_residual <- abs(raw_norm^2-1)
    point$angle_violation_deg <- max(0,point$angle_deg-angle_deg)
    point$cosine_violation <- max(0,point$max_cosine-problem$cosine_cap)
    point$oriented_cap_constraint <- utils::tail(constraints$value(point$x),1L)
    point$feasible <- is.finite(raw_norm) && abs(raw_norm-1)<=feasibility_tolerance &&
      point$max_cosine<=problem$cosine_cap+feasibility_tolerance &&
      point$angle_deg<=angle_deg+angle_tolerance_deg &&
      sum(point$lambda*unit0)>=-feasibility_tolerance
    point$success <- status %in% 1:4
    point$status <- as.integer(status); point$message <- as.character(message)
    point$source <- source; point$iterations <- iterations
    point$cap_edge <- point$angle_deg >= angle_deg-edge_tolerance_deg
    point$loading_constraint_edge <- point$max_cosine >= problem$cosine_cap-1e-5
    point
  }
  # The center is an explicitly labelled fallback, never optimizer convergence.
  baseline <- summarize(center_x,0L,"Feasible center baseline; not an optimization run.","center_baseline",0L)
  starts <- list(center_x)
  add_start <- function(lambda) {
    if (length(lambda)!=problem$d || any(!is.finite(lambda))) return(invisible(NULL))
    if (sum(lambda*unit0)<0) lambda <- -lambda
    x <- cp_alignment_x_from_lambda(problem,lambda)
    angle <- cp_alignment_angle_deg(lambda,unit0)
    if (angle>angle_deg+1e-8 || max(problem$g_ineq(x))>feasibility_tolerance) return(invisible(NULL))
    duplicate <- any(vapply(starts,function(y) abs(sum(x*y))>=1-1e-12,logical(1L)))
    if (!duplicate) starts[[length(starts)+1L]] <<- x
    invisible(NULL)
  }
  if (length(extra_lambdas)) for (lambda in extra_lambdas) {
    if (length(starts)>=nstarts) break
    add_start(lambda)
  }
  attempts <- 0L
  while (length(starts)<nstarts && attempts<max(1000L,200L*nstarts)) {
    v <- stats::rnorm(problem$d); v <- v-sum(v*unit0)*unit0
    norm <- sqrt(sum(v*v)); attempts <- attempts+1L
    if (norm<1e-12) next
    v <- v/norm
    theta <- angle_deg*pi/180*stats::runif(1L)^.5
    add_start(cos(theta)*unit0+sin(theta)*v)
  }
  runs <- lapply(seq_along(starts),function(i) {
    x0 <- starts[[i]]
    fit <- tryCatch(nloptr::nloptr(
      x0=x0,
      eval_f=function(x) problem$E_squared(x)$value,
      eval_grad_f=function(x) problem$E_squared(x)$gradient,
      eval_g_ineq=constraints$value,eval_jac_g_ineq=constraints$jacobian,
      eval_g_eq=function(x) sum(x*x)-1,
      eval_jac_g_eq=function(x) matrix(2*x,nrow=1L),
      opts=list(algorithm="NLOPT_LD_SLSQP",maxeval=as.integer(maxeval),
                xtol_rel=xtol_rel,ftol_abs=ftol_abs,print_level=0L,
                tol_constraints_ineq=rep(1e-10,problem$d+1L),
                tol_constraints_eq=1e-10)
    ),error=function(e) list(solution=x0,status=-999L,message=conditionMessage(e),iterations=0L))
    summarize(fit$solution,fit$status,fit$message,paste0("start_",i),fit$iterations)
  })
  feasible <- Filter(function(x) isTRUE(x$feasible),c(runs,list(baseline)))
  # On near ties retain an actually converged optimizer instead of the baseline.
  best <- if (length(feasible)) {
    values <- vapply(feasible,function(x)x$E1,numeric(1L))
    near <- which(values^2<=min(values^2)+ftol_abs & values<=min(values)+1e-10)
    successful <- near[vapply(feasible[near],function(x)isTRUE(x$success),logical(1L))]
    pool <- if(length(successful))successful else seq_along(feasible)
    feasible[[pool[which.min(values[pool])]]]
  } else NULL
  list(best=best,runs=runs,baseline=baseline,
       n_candidates_feasible=length(feasible),
       n_candidates_successful=sum(vapply(runs,function(x)isTRUE(x$feasible)&&isTRUE(x$success),logical(1L))),
       n_optimizer_runs=length(runs),configuration=list(angle_deg=angle_deg,
       cosine_cap=problem$cosine_cap,nstarts=nstarts,maxeval=maxeval,
       feasibility_tolerance=feasibility_tolerance,angle_tolerance_deg=angle_tolerance_deg))
}

cp_local_alignment_case <- function(checkpoint,angle_deg=5,
                                    cap_cosines=c(.95,.90),nstarts=5L,
                                    maxeval=300L,seed=123L,
                                    feasibility_tolerance=2e-6,
                                    angle_tolerance_deg=1e-4,
                                    edge_tolerance_deg=.05) {
  if (!requireNamespace("nloptr",quietly=TRUE)) stop("Package 'nloptr' is required.")
  if (any(!is.finite(cap_cosines)) || any(cap_cosines<=0 | cap_cosines>=1)) stop("Invalid cosine caps.")
  dg <- checkpoint$diagnostics
  if (is.null(dg) || !length(dg$rho_local_minima) || !length(dg$E_local_minima)) {
    stop("Checkpoint must contain saved rho and E1 numerical local minima.")
  }
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  previous_seed <- if(had_seed)get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit(if(had_seed)assign(".Random.seed",previous_seed,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(as.integer(seed))
  core <- dg$core
  if(is.null(core))core <- list(A=checkpoint$Ar,B=checkpoint$Br)
  base_problem <- cp_diagnostic_problem(core,max(cap_cosines))
  caps <- sort(unique(as.numeric(cap_cosines)))
  problems <- lapply(caps,function(cap)cp_alignment_change_cap(base_problem,cap))
  E_pool <- lapply(seq_along(dg$E_local_minima),function(j) {
    old <- dg$E_local_minima[[j]]
    point <- base_problem$summarize(cp_alignment_x_from_lambda(base_problem,old$lambda))
    point$id <- if(is.null(old$minimum_id))j else old$minimum_id
    point$success <- isTRUE(old$success); point$status <- old$status
    point
  })
  E_pool <- Filter(function(x)x$max_cosine<=.90+feasibility_tolerance,E_pool)
  if(!length(E_pool))stop("No saved E1 candidate is feasible at cosine <=0.90.")
  E_best <- E_pool[[which.min(vapply(E_pool,function(x)x$E1,numeric(1L)))]]
  rows <- list(); details <- list()
  for (rho_index in seq_along(dg$rho_local_minima)) {
    rho_old <- dg$rho_local_minima[[rho_index]]
    center <- base_problem$summarize(cp_alignment_x_from_lambda(base_problem,rho_old$lambda))
    rho_order <- if(is.null(rho_old$minimum_id))rho_index else rho_old$minimum_id
    E_center_angles <- vapply(E_pool,function(x)cp_alignment_angle_deg(center$lambda,x$lambda),numeric(1L))
    nearest_center <- E_pool[[which.min(E_center_angles)]]
    extra <- lapply(E_pool[order(vapply(E_pool,function(x)x$E1,numeric(1L)))],function(x)x$lambda)
    previous_best <- NULL
    for (cap_index in seq_along(caps)) {
      cap <- caps[cap_index]
      warm <- if(is.null(previous_best))extra else c(list(previous_best$lambda),extra)
      search <- cp_alignment_cap_search(problems[[cap_index]],center$lambda,angle_deg,nstarts,maxeval,
        extra_lambdas=warm,feasibility_tolerance=feasibility_tolerance,
        angle_tolerance_deg=angle_tolerance_deg,edge_tolerance_deg=edge_tolerance_deg)
      best <- search$best
      if(is.null(best))stop("No feasible cap point; inspect checkpoint rho minima and tolerances.")
      # A feasible result for a tighter cap remains a valid candidate for a looser cap.
      if(!is.null(previous_best) && previous_best$E1 < best$E1-1e-10) {
        best <- previous_best
        best$source <- "tighter_cap_warm_candidate"
        best$success <- FALSE; best$status <- 0L
        best$message <- "Feasible tighter-cap candidate retained; no convergence claim for this cap."
        best$loading_constraint_edge <- best$max_cosine>=cap-1e-5
        best$cosine_violation <- max(0,best$max_cosine-cap)
        search$best <- best
      }
      previous_best <- best
      E_cap_angles <- vapply(E_pool,function(x)cp_alignment_angle_deg(best$lambda,x$lambda),numeric(1L))
      nearest_cap <- E_pool[[which.min(E_cap_angles)]]
      row <- data.frame(
        regime=checkpoint$regime,seed=checkpoint$seed,T=checkpoint$T,
        qmle_converged=isTRUE(checkpoint$fit$converged),rho_order=rho_order,
        cap_cosine=cap,angle_deg=angle_deg,
        E1_center=center$E1,rho_center=center$rho,E1_cap=best$E1,rho_cap=best$rho,
        center_to_cap_deg=best$angle_deg,center_to_cap_matrix_deg=best$matrix_angle_deg,
        center_to_E1_best_deg=cp_alignment_angle_deg(center$lambda,E_best$lambda),
        cap_to_E1_best_deg=cp_alignment_angle_deg(best$lambda,E_best$lambda),
        center_to_E1_best_matrix_deg=cp_alignment_angle_deg(center$x,E_best$x),
        cap_to_E1_best_matrix_deg=cp_alignment_angle_deg(best$x,E_best$x),
        center_to_E1_nearest_deg=min(E_center_angles),cap_to_E1_nearest_deg=min(E_cap_angles),
        center_to_E1_nearest_matrix_deg=cp_alignment_angle_deg(center$x,nearest_center$x),
        cap_to_E1_nearest_matrix_deg=cp_alignment_angle_deg(best$x,nearest_cap$x),
        E1_best=E_best$E1,E1_nearest_to_center=nearest_center$E1,E1_nearest_to_cap=nearest_cap$E1,
        E1_gap_best=best$E1-E_best$E1,
        E1_gap_nearest_to_center=best$E1-nearest_center$E1,
        E1_gap_nearest_to_cap=best$E1-nearest_cap$E1,
        E1_improvement_over_center=center$E1-best$E1,
        E1_best_id=E_best$id,E1_nearest_center_id=nearest_center$id,E1_nearest_cap_id=nearest_cap$id,
        E1_candidate_in_5deg=any(E_center_angles<=angle_deg+angle_tolerance_deg),
        cap_within_5deg_of_E1_candidate=min(E_cap_angles)<=angle_deg+angle_tolerance_deg,
        cap_max_cosine=best$max_cosine,cap_feasible_at_090=best$max_cosine<=.90+feasibility_tolerance,
        cap_edge=best$cap_edge,cap_loading_edge=best$loading_constraint_edge,
        optimizer_success=isTRUE(best$success),optimizer_status=best$status,
        optimizer_message=best$message,optimizer_source=best$source,
        n_candidates_feasible=search$n_candidates_feasible,
        n_candidates_successful=search$n_candidates_successful,n_optimizer_runs=search$n_optimizer_runs,
        cap_eq_residual=best$eq_residual,cap_angle_violation_deg=best$angle_violation_deg,
        cap_cosine_violation=max(0,best$max_cosine-cap),
        rho_center_success=isTRUE(rho_old$success),E1_best_success=isTRUE(E_best$success),
        E1_nearest_center_success=isTRUE(nearest_center$success),E1_nearest_cap_success=isTRUE(nearest_cap$success),
        n_E1_candidates=length(E_pool),stringsAsFactors=FALSE)
      # Current experiments have d=5 or d=8. Pad to eight for ordinary
      # rbind across regimes; preserve every coordinate for larger d.
      for(j in seq_len(max(8L,base_problem$d))) {
        row[[sprintf("lambda_%02d",j)]] <- if(j<=base_problem$d)best$lambda[j] else NA_real_
        row[[sprintf("lambda_unit_%02d",j)]] <- if(j<=base_problem$d)best$lambda[j]/sqrt(sum(best$lambda^2)) else NA_real_
        row[[sprintf("center_lambda_%02d",j)]] <- if(j<=base_problem$d)center$lambda[j] else NA_real_
      }
      rows[[length(rows)+1L]] <- row
      details[[length(details)+1L]] <- list(rho_order=rho_order,cap_cosine=cap,
        center=center,search=search,nearest_E1_to_center=nearest_center,
        nearest_E1_to_cap=nearest_cap,E1_best=E_best)
    }
  }
  result <- do.call(rbind,rows)
  result <- result[order(result$rho_order,match(result$cap_cosine,cap_cosines)),,drop=FALSE]
  rownames(result) <- NULL
  list(rows=result,details=details,configuration=list(angle_deg=angle_deg,
    cap_cosines=cap_cosines,nstarts=nstarts,maxeval=maxeval,seed=seed,
    angle_coordinates="unit lambda in matched unit-column Ar/Br coordinates; sign ignored",
    comparison="saved E1 numerical local candidates at cosine <=0.90; not certified exhaustive"))
}


# Truth-near multiple starts for the p,q grid simulation.
# Requires R_qmle_backend.R. Only base R is used here.
# These starts use simulation truth and are NOT a blind-data estimator.
# The objective still estimates every diagonal noise variance by default.

cp_grid_option <- function(cfg, name, default) {
  value <- cfg[[name]]
  if (is.null(value)) default else value
}

cp_grid_perturb_directions <- function(A, max_angle_deg=3) {
  A <- cp_normalize_columns(as.matrix(A))
  if (!is.finite(max_angle_deg) || max_angle_deg<0 || max_angle_deg>=90)
    stop("truth_angle_deg must be in [0,90).")
  if (max_angle_deg==0) return(A)
  for(j in seq_len(ncol(A))) {
    unit <- A[,j]
    direction <- stats::rnorm(nrow(A))
    direction <- direction-sum(direction*unit)*unit
    norm <- sqrt(sum(direction^2))
    if (!is.finite(norm) || norm<1e-14) {
      direction <- numeric(nrow(A));direction[which.min(abs(unit))] <- 1
      direction <- direction-sum(direction*unit)*unit
      norm <- sqrt(sum(direction^2))
    }
    theta <- stats::runif(1L,0,max_angle_deg)*pi/180
    A[,j] <- cos(theta)*unit+sin(theta)*direction/norm
  }
  A
}

cp_grid_pack_exact <- function(obj,A,B,M,psi=rep(1,obj$N),repair_floor=1e-10) {
  A <- as.matrix(A);B <- as.matrix(B);M <- as.matrix(M)
  if (!identical(dim(A),c(obj$p,obj$d)) ||
      !identical(dim(B),c(obj$q,obj$d)) ||
      !identical(dim(M),c(obj$d,obj$d))) stop("Invalid dimensions in exact QMLE packing.")
  if(any(!is.finite(c(A,B,M,psi))) || any(colSums(A*A)<=0) ||
     any(colSums(B*B)<=0)) stop("Invalid or zero loading columns in QMLE start.")
  if(length(psi)==1L)psi <- rep(psi,obj$N)
  if(length(psi)!=obj$N || any(psi<=0))stop("Invalid positive noise variances in QMLE start.")
  M <- (M+t(M))/2
  initial_eigen <- eigen(M,symmetric=TRUE,only.values=TRUE)$values
  L <- tryCatch(t(chol(M)),error=function(e)NULL)
  repaired <- is.null(L)
  used_floor <- 0
  if(repaired) {
    ev <- eigen(M,symmetric=TRUE)
    used_floor <- max(repair_floor, .Machine$double.eps*max(1,max(abs(ev$values)))*10)
    M <- tcrossprod(sweep(ev$vectors,2L,sqrt(pmax(ev$values,used_floor)),"*"))
    L <- t(chol((M+t(M))/2))
  }
  mask <- lower.tri(L,diag=TRUE)
  positions <- which(mask,arr.ind=TRUE)
  lv <- L[mask]
  diagonal <- which(positions[,1L]==positions[,2L])
  lv[diagonal] <- log(diag(L))
  unbounded <- c(as.vector(cp_normalize_columns(A)),as.vector(cp_normalize_columns(B)),lv)
  if(obj$estimate_noise)unbounded <- c(unbounded,log(psi))
  par <- pmax(obj$lower,pmin(obj$upper,unbounded))
  list(par=par,covariance_repaired=repaired,covariance_repair_floor=used_floor,
       covariance_min_eigen_before=min(initial_eigen),
       clipped_parameters=sum(par!=unbounded),
       exact_before_bounds=!repaired && !any(par!=unbounded))
}

cp_grid_initial_angles <- function(candidate,truth) {
  cosine <- colSums(cp_normalize_columns(candidate)*cp_normalize_columns(truth))
  acos(pmax(-1,pmin(1,cosine)))*180/pi
}

cp_grid_qmle <- function(Y,design,cfg,seed,previous_fit=NULL) {
  p <- as.integer(design$p);q <- as.integer(design$q);d <- as.integer(design$d)
  regime <- as.character(design$regime)
  maxit_by_regime <- cp_grid_option(cfg,"maxit_by_regime",
    c(interior=300L,boundary=400L,exterior=500L))
  maxit <- as.integer(maxit_by_regime[[regime]])
  truth_starts <- as.integer(cp_grid_option(cfg,"truth_starts",6L))
  angle_deg <- cp_grid_option(cfg,"truth_angle_deg",3)
  chol_sd <- cp_grid_option(cfg,"truth_chol_sd",.05)
  noise_sd <- cp_grid_option(cfg,"truth_noise_log_sd",.05)
  estimate_noise <- isTRUE(cp_grid_option(cfg,"estimate_noise",TRUE))
  noise_bounds <- cp_grid_option(cfg,"noise_bounds",c(.05,20))
  factr <- cp_grid_option(cfg,"factr",1000)
  pgtol <- cp_grid_option(cfg,"pgtol",1e-8)
  lmm <- as.integer(cp_grid_option(cfg,"lmm",15L))
  save_all <- isTRUE(cp_grid_option(cfg,"save_all_start_fits",FALSE))
  if(length(maxit)!=1L || !is.finite(maxit) || maxit<1L)stop("Missing positive regime-specific maxit.")
  if(length(truth_starts)!=1L || truth_starts<1L || !is.finite(truth_starts))stop("truth_starts must be positive.")
  if(any(!is.finite(c(chol_sd,noise_sd))) || min(chol_sd,noise_sd)<0)stop("Invalid truth perturbation SD.")
  if(!is.finite(factr) || factr<0 || !is.finite(pgtol) || pgtol<0 || lmm<1L)
    stop("Invalid optimization tolerances.")
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_seed <- if(had_seed)get(".Random.seed",envir=.GlobalEnv)else NULL
  on.exit(if(had_seed)assign(".Random.seed",old_seed,envir=.GlobalEnv)else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(as.integer(seed))
  began <- proc.time()[["elapsed"]]
  obj <- cp_qmle_objective(Y,p,q,d,estimate_noise=estimate_noise,fixed_noise=1,
    noise_bounds=noise_bounds)
  starts <- list();source_names <- character();initial_info <- list()
  add_start <- function(A,B,M,psi,source) {
    packed <- cp_grid_pack_exact(obj,A,B,M,psi)
    i <- length(starts)+1L
    starts[[i]] <<- packed$par
    source_names[i] <<- source
    angle_A <- cp_grid_initial_angles(A,design$A)
    angle_B <- cp_grid_initial_angles(B,design$B)
    initial_info[[i]] <<- list(
      initial_A_angle_mean=mean(angle_A),initial_A_angle_max=max(angle_A),
      initial_B_angle_mean=mean(angle_B),initial_B_angle_max=max(angle_B),
      initial_M_relative_error=sqrt(sum((M-design$Mx)^2)/sum(design$Mx^2)),
      initial_noise_log_rms=sqrt(mean(log(psi)^2)),
      covariance_repaired=packed$covariance_repaired,
      covariance_repair_floor=packed$covariance_repair_floor,
      covariance_min_eigen_before=packed$covariance_min_eigen_before,
      clipped_parameters=packed$clipped_parameters,
      exact_before_bounds=packed$exact_before_bounds)
    invisible(NULL)
  }
  add_start(design$A,design$B,design$Mx,rep(1,p*q),"truth_exact")
  L0 <- t(chol((design$Mx+t(design$Mx))/2))
  if(truth_starts>1L)for(i in 2:truth_starts) {
    A <- sqrt(p)*cp_grid_perturb_directions(design$A,angle_deg)
    B <- sqrt(q)*cp_grid_perturb_directions(design$B,angle_deg)
    L <- L0
    diag(L) <- diag(L0)*exp(stats::rnorm(d,0,chol_sd))
    if(d>1L)for(row in 2:d)for(column in seq_len(row-1L))
      L[row,column] <- L0[row,column]+stats::rnorm(1L,0,chol_sd)*sqrt(design$Mx[row,row])
    psi <- exp(stats::rnorm(p*q,0,noise_sd))
    add_start(A,B,tcrossprod(L),psi,paste0("truth_near_",i-1L))
  }
  warm_index <- NULL
  if(!is.null(previous_fit)) {
    add_start(previous_fit$A,previous_fit$B,previous_fit$M,previous_fit$psi,"previous_T")
    warm_index <- length(starts)
  }
  n_starts <- length(starts)
  logs <- vector("list",n_starts)
  all_fits <- if(save_all)vector("list",n_starts)else NULL
  best <- NULL;best_value <- Inf;best_start <- NA_integer_
  for(i in seq_along(starts)) {
    tick <- proc.time()[["elapsed"]]
    initial_loss <- obj$fn(starts[[i]])
    warm_loss_delta <- if(!is.null(warm_index) && i==warm_index && length(previous_fit$nll)==1L)
      initial_loss-as.numeric(previous_fit$nll)else NA_real_
    # The warm loss changes with T. The parameters are packed exactly before
    # bounds; old and current-T objective values must not be equated.
    optimized <- tryCatch(optim(starts[[i]],obj$fn,obj$gr,method="L-BFGS-B",
      lower=obj$lower,upper=obj$upper,
      control=list(maxit=maxit,factr=factr,pgtol=pgtol,lmm=lmm)),error=function(e)
        list(par=starts[[i]],convergence=9999L,message=paste("Optimizer error; start retained:",conditionMessage(e)),
             counts=c("function"=0L,gradient=0L)))
    evaluated <- obj$evaluate(optimized$par)
    pg <- obj$projected_gradient(optimized$par,evaluated$gradient)
    counts <- optimized$counts
    row <- data.frame(start=i,source=source_names[i],nll=evaluated$value,
      initial_nll=initial_loss,convergence=optimized$convergence,
      converged=optimized$convergence==0L,
      function_evaluations=if(is.null(counts["function"]))NA_real_ else unname(counts["function"]),
      gradient_evaluations=if(is.null(counts["gradient"]))NA_real_ else unname(counts["gradient"]),
      gradient_norm=max(abs(evaluated$gradient)),projected_gradient_norm=max(abs(pg)),
      elapsed_seconds=proc.time()[["elapsed"]]-tick,
      message=if(is.null(optimized$message))""else optimized$message,
      maxit_per_start=maxit,factr=factr,pgtol=pgtol,
      warm_nll_delta_across_T=warm_loss_delta,
      stringsAsFactors=FALSE)
    logs[[i]] <- cbind(row,as.data.frame(initial_info[[i]],stringsAsFactors=FALSE))
    if(save_all) {
      parts_i <- obj$unpack(optimized$par)
      all_fits[[i]] <- list(A=parts_i$A,B=parts_i$B,M=tcrossprod(parts_i$L),L=parts_i$L,
                           psi=parts_i$psi,nll=evaluated$value,convergence=optimized$convergence)
    }
    # A caught optimizer exception must not turn an unoptimized truth start
    # into a reported estimator, even if that fallback has a finite loss.
    usable<-optimized$convergence!=9999L && is.finite(evaluated$value) && evaluated$value<1e29
    if(usable && evaluated$value<best_value) {
      best_value <- evaluated$value;best <- optimized;best_start <- i
    }
  }
  if(is.null(best))stop("Every truth-near QMLE run failed or returned an invalid objective; no estimator selected.")
  logs <- do.call(rbind,logs)
  parts <- obj$unpack(best$par)
  M <- tcrossprod(parts$L)
  Z <- cp_khatri_rao(parts$A,parts$B)
  signal_singular <- svd(Z%*%parts$L,nu=0L,nv=0L)$d
  list(A=parts$A,B=parts$B,M=M,L=parts$L,psi=parts$psi,
    nll=best_value,converged=logs$converged[best_start],
    convergence=logs$convergence[best_start],message=logs$message[best_start],
    gradient_norm=logs$gradient_norm[best_start],
    projected_gradient_norm=logs$projected_gradient_norm[best_start],
    best_start=best_start,all_starts=logs,maxit_per_start=maxit,n_starts=n_starts,
    n_truth_starts=truth_starts,n_previous_T_starts=as.integer(!is.null(previous_fit)),
    pca_method="not_used_truth_near_starts",noise_mode=if(estimate_noise)"estimated_diagonal"else"known_diagonal",
    noise_at_lower_bound=if(estimate_noise)sum(parts$psi<=noise_bounds[1L]*(1+1e-4))else 0L,
    noise_at_upper_bound=if(estimate_noise)sum(parts$psi>=noise_bounds[2L]*(1-1e-4))else 0L,
    M_min_eigen=min(eigen(M,symmetric=TRUE,only.values=TRUE)$values),
    signal_singular_values=signal_singular,
    elapsed_seconds=proc.time()[["elapsed"]]-began,
    truth_near_initialization=TRUE,
    initialization=list(truth_starts=truth_starts,truth_angle_deg=angle_deg,
      truth_chol_sd=chol_sd,truth_noise_log_sd=noise_sd,
      statement="Simulation-truth initialization; not a blind-data estimator. Noise remains fitted."),
    all_start_estimates=all_fits)
}


# Geometry for the p/q grid experiment. Requires R_design.R,
# R_diagnostics.R and R_local_alignment_core.R. No QMLE fitting here.

cp_grid_get <- function(cfg, name, default) {
  if (is.null(cfg[[name]])) default else cfg[[name]]
}

cp_grid_bind_frames <- function(frames) {
  frames <- Filter(function(x)!is.null(x) && nrow(x)>0L,frames)
  if(!length(frames))return(data.frame())
  fields<-unique(unlist(lapply(frames,names),use.names=FALSE))
  frames<-lapply(frames,function(x){
    for(field in setdiff(fields,names(x)))x[[field]]<-NA
    x[,fields,drop=FALSE]
  })
  do.call(rbind,frames)
}

cp_grid_true_core_directions <- function(regime,n=.5) {
  truth<-core_loading_design(regime,n)
  A<-truth$Ar_unit;B<-truth$Br_unit
  combination_norm<-rep(1,ncol(A))
  if(regime=="exterior") {
    M<-A%*%diag(truth$extra_lambda)%*%t(B)
    sv<-svd(M,nu=1L,nv=1L)
    combination_norm<-c(combination_norm,sqrt(sum(M*M)))
    A<-cbind(A,sv$u[,1L]);B<-cbind(B,sv$v[,1L])
  }
  list(A=A,B=B,Z=khatri_rao(A,B),n_original=truth$d,
    n_known=ncol(A),combination_norm=combination_norm,truth=truth)
}

cp_grid_population_bases <- function(regime,n=.5,rank_tolerance=1e-10) {
  # The lazy cache is local to each R worker, so it need not be exported.
  cache_name<-".cp_grid_population_reference_cache"
  if(!exists(cache_name,envir=.GlobalEnv,inherits=FALSE))
    assign(cache_name,new.env(parent=emptyenv()),envir=.GlobalEnv)
  cache<-get(cache_name,envir=.GlobalEnv,inherits=FALSE)
  key<-paste(regime,format(n,digits=17),format(rank_tolerance,digits=17),sep=";")
  if(exists(key,envir=cache,inherits=FALSE))return(get(key,envir=cache,inherits=FALSE))
  truth<-cp_grid_true_core_directions(regime,n)
  selections<-if(regime=="exterior")utils::combn(seq_len(truth$n_known),truth$n_original,simplify=FALSE) else
    list(seq_len(truth$n_original))
  rows<-list()
  for(ids in selections) {
    A<-truth$A[,ids,drop=FALSE];B<-truth$B[,ids,drop=FALSE]
    minors<-cp_minor_matrix(A,B)
    js<-cp_jacobian_spectra(A,B,minors$C,minors$pairs)
    Z<-khatri_rao(A,B)
    zsv<-svd(Z,nu=0L,nv=0L)$d
    csv<-svd(minors$C,nu=0L,nv=0L)$d
    c_rank<-sum(csv>max(csv)*rank_tolerance)
    basis_key<-paste(sort(ids),collapse="_")
    for(ell in seq_along(ids))for(index in seq_len(length(ids)-1L)) {
      values<-js$orig[ell,]
      rows[[length(rows)+1L]]<-data.frame(regime=regime,
        population_basis_key=basis_key,ell=ell,assigned_true_direction_id=ids[ell],
        svd_index=index,orig=js$orig[ell,index],nor=js$nor[ell,index],
        gram_eigen_orig=js$orig[ell,index]^2,gram_eigen_nor=js$nor[ell,index]^2,
        J_rank=sum(values>max(values)*rank_tolerance),basis_Z_min_svd=min(zsv),
        basis_C_rank=c_rank,basis_C_ncol=ncol(minors$C),stringsAsFactors=FALSE)
    }
  }
  result<-do.call(rbind,rows)
  assign(key,result,envir=cache)
  result
}

cp_grid_known_matching <- function(core,matched_A,matched_B,design) {
  truth<-cp_grid_true_core_directions(design$regime,
    if(is.null(design$n) || !is.finite(design$n)).5 else design$n)
  true_A<-design$QA%*%truth$A;true_B<-design$QB%*%truth$B
  true_Z<-khatri_rao(true_A,true_B)
  fitted_Z<-khatri_rao(cp_unit_columns(matched_A),cp_unit_columns(matched_B))
  score<-abs(crossprod(fitted_Z,true_Z))
  d<-nrow(score);m<-ncol(score)
  square<-rbind(score,matrix(0,m-d,m))
  assignment<-cp_max_assignment(square)[seq_len(d)]
  values<-score[cbind(seq_len(d),assignment)]
  missing<-setdiff(seq_len(m),assignment)
  basis_key<-paste(sort(assignment),collapse="_")
  original_values<-diag(score[,seq_len(d),drop=FALSE])
  matching<-data.frame(ell=seq_len(d),assigned_true_direction_id=assignment,
    matched_true_direction_cosine=values,original_true_direction_cosine=original_values,
    population_basis_key=basis_key,
    missing_true_direction_id=if(length(missing)==1L)missing else NA_integer_,
    extra_selected=design$regime=="exterior" && m%in%assignment,
    matching_to_original_mean=mean(original_values),
    matching_to_known_set_mean=mean(values),matching_to_known_set_min=min(values),
    is_original_basis_set=identical(sort(assignment),seq_len(d)),
    near_original_basis=mean(original_values)>=.98,n_known_true_directions=m,
    stringsAsFactors=FALSE)
  list(rows=matching,truth=truth,true_A=true_A,true_B=true_B,true_Z=true_Z,
    fitted_Z=fitted_Z,assignment=assignment,score=score)
}

cp_grid_project_known <- function(problem,qa,qb,matching,cosine_cap=.90) {
  rows<-list();points<-list();extra_lambdas<-list();extra_labels<-character()
  for(j in seq_len(ncol(matching$true_A))) {
    ca<-drop(crossprod(qa$Q,matching$true_A[,j]))
    cb<-drop(crossprod(qb$Q,matching$true_B[,j]))
    x_projection<-drop(crossprod(problem$Q,as.vector(outer(ca,cb))))
    projection_norm<-sqrt(sum(x_projection*x_projection))
    if(projection_norm>1+1e-8)stop("Projection norm exceeds one; core/physical mapping is inconsistent.")
    point<-if(projection_norm>1e-14)problem$summarize(x_projection/projection_norm)else NULL
    similarities<-abs(drop(crossprod(matching$fitted_Z,matching$true_Z[,j])))
    assigned<-which(matching$assignment==j)
    true_sv<-svd(outer(matching$truth$A[,j],matching$truth$B[,j]),nu=0L,nv=0L)$d
    rho_defined<-!is.null(point) && sqrt(sum((point$lambda[problem$pairs[1L,]]*
      point$lambda[problem$pairs[2L,]])^2))>1e-14
    feasible<-!is.null(point) && point$max_cosine<=cosine_cap+1e-10
    row<-data.frame(true_direction_id=j,
      direction_type=if(j>matching$truth$n_original)"extra_rankone"else"original_loading",
      is_extra=j>matching$truth$n_original,population_projection_norm=1,
      true_combination_norm=matching$truth$combination_norm[j],
      population_rankone_error=sqrt(sum(true_sv[-1L]^2)/sum(true_sv^2)),
      projection_norm=projection_norm,projection_residual=sqrt(max(0,1-projection_norm^2)),
      E1_projected=if(is.null(point))NA_real_ else point$E1,
      rho_projected=if(rho_defined)point$rho else NA_real_,rho_defined=rho_defined,
      max_loading_cosine=if(is.null(point))NA_real_ else point$max_cosine,
      feasible_primary=feasible,nearest_estimated_ell=which.max(similarities),
      nearest_estimated_cosine=max(similarities),
      assigned_estimated_ell=if(length(assigned))assigned[1L]else NA_integer_,
      stringsAsFactors=FALSE)
    for(k in seq_len(max(8L,problem$d))) {
      row[[sprintf("lambda_%02d",k)]]<-if(!is.null(point)&&k<=problem$d)point$lambda[k]else NA_real_
      row[[sprintf("lambda_unit_%02d",k)]]<-if(!is.null(point)&&k<=problem$d)
        point$lambda[k]/sqrt(sum(point$lambda^2))else NA_real_
    }
    rows[[length(rows)+1L]]<-row
    points[[j]]<-list(id=j,x_projection=x_projection,point=point)
    if(feasible) {
      duplicate<-length(extra_lambdas)>0L && any(vapply(extra_lambdas,function(lam)
        cp_alignment_angle_deg(lam,point$lambda)<1e-5,logical(1L)))
      if(!duplicate) {
        extra_lambdas[[length(extra_lambdas)+1L]]<-point$lambda
        extra_labels<-c(extra_labels,paste0("truth_direction_",j,"_projection"))
      }
    }
  }
  list(rows=do.call(rbind,rows),points=points,extra_lambdas=extra_lambdas,extra_labels=extra_labels)
}

cp_grid_search_geometry <- function(core,problem,cfg,seed,extra_lambdas=list(),extra_labels=character()) {
  had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  previous_seed<-if(had_seed)get(".Random.seed",envir=.GlobalEnv)else NULL
  on.exit(if(had_seed)assign(".Random.seed",previous_seed,envir=.GlobalEnv)else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(as.integer(seed))
  nrandom<-as.integer(cp_grid_get(cfg,"diagnostic_random_starts",30L))
  maxeval<-as.integer(cp_grid_get(cfg,"diagnostic_maxeval",400L))
  cross_limit<-as.integer(cp_grid_get(cfg,"diagnostic_cross_start_limit",5L))
  cross_rounds<-as.integer(cp_grid_get(cfg,"diagnostic_cross_rounds",2L))
  if(length(cross_rounds)!=1L || is.na(cross_rounds) || cross_rounds<0L)
    stop("diagnostic_cross_rounds must be one nonnegative integer.")
  ftol<-cp_grid_get(cfg,"diagnostic_ftol_abs",1e-13)
  cluster_tol<-cp_grid_get(cfg,"diagnostic_cluster_cosine_tolerance",1e-5)
  feasibility_tol<-cp_grid_get(cfg,"diagnostic_feasibility_tolerance",2e-6)
  d<-problem$d
  starts<-lapply(extra_lambdas,function(lam)cp_alignment_x_from_lambda(problem,lam))
  labels<-extra_labels
  if(length(labels)!=length(starts))labels<-paste0("truth_projection_",seq_along(starts))
  count_random<-0L;attempts<-0L
  while(count_random<nrandom && attempts<max(10000L,nrandom*1000L)) {
    x<-stats::rnorm(d);x<-x/sqrt(sum(x*x));attempts<-attempts+1L
    if(max(problem$g_ineq(x))<=-1e-10) {
      starts[[length(starts)+1L]]<-x;count_random<-count_random+1L
      labels<-c(labels,sprintf("random_%03d",count_random))
    }
  }
  if(count_random<nrandom)stop("Could not generate enough random geometry starts.")
  optimize_one<-function(x0,objective,source) {
    fun<-if(objective=="rho")problem$rho_squared else problem$E_squared
    fit<-tryCatch(nloptr::nloptr(x0=x0,
      eval_f=function(x)fun(x)$value,eval_grad_f=function(x)fun(x)$gradient,
      eval_g_ineq=problem$g_ineq,eval_jac_g_ineq=problem$jac_ineq,
      eval_g_eq=function(x)sum(x*x)-1,eval_jac_g_eq=function(x)matrix(2*x,nrow=1L),
      opts=list(algorithm="NLOPT_LD_SLSQP",maxeval=maxeval,xtol_rel=1e-9,
        ftol_abs=ftol,print_level=0L,tol_constraints_ineq=rep(1e-10,d),tol_constraints_eq=1e-10)),
      error=function(e)list(solution=x0,status=-999L,message=conditionMessage(e),iterations=0L))
    norm<-sqrt(sum(fit$solution^2))
    result<-problem$summarize(fit$solution)
    result$feasible<-is.finite(norm)&&abs(norm-1)<=feasibility_tol &&
      result$max_cosine<=problem$cosine_cap+feasibility_tol
    result$eq_residual<-abs(norm^2-1)
    result$objective_value<-if(objective=="rho")result$rho else result$E1
    result$optimized_objective_squared<-result$objective_value^2
    result$success<-fit$status%in%1:4;result$status<-fit$status
    result$message<-fit$message;result$source<-source
    result$iterations<-if(is.null(fit$iterations))NA_integer_ else fit$iterations
    result$loading_constraint_edge<-result$max_cosine>=problem$cosine_cap-1e-5
    result
  }
  rho_runs<-lapply(seq_along(starts),function(i)optimize_one(starts[[i]],"rho",labels[i]))
  cluster<-function(runs,objective) {
    valid<-Filter(function(x)isTRUE(x$feasible),runs)
    out<-cp_cluster_minima(valid,objective,cluster_tol)
    for(j in seq_along(out)) {
      hits<-valid[vapply(valid,function(x)abs(sum(x$x*out[[j]]$x))>=1-cluster_tol,logical(1L))]
      out[[j]]$all_hit_sources<-paste(unique(vapply(hits,function(x)x$source,character(1L))),collapse=";")
      out[[j]]$truth_guided_hit<-any(grepl("^truth",vapply(hits,function(x)x$source,character(1L))))
    }
    out
  }
  rho_local<-cluster(rho_runs,"rho")
  E_starts<-starts;E_labels<-labels
  if(length(rho_local)&&cross_limit>0L)for(j in seq_len(min(length(rho_local),cross_limit))) {
    E_starts[[length(E_starts)+1L]]<-rho_local[[j]]$x
    E_labels<-c(E_labels,paste0("rho_candidate_",j))
  }
  E_runs<-lapply(seq_along(E_starts),function(i)optimize_one(E_starts[[i]],"E1",E_labels[i]))
  E_local<-cluster(E_runs,"E1")
  # A fixed number of cross-objective rounds, with the SAME maxeval per run.
  # Candidate witnesses are optimized before they enter the other local pool.
  # Keep every original run and all discovered local candidates/provenance.
  if(cross_rounds>0L)for(round in seq_len(cross_rounds)) {
    if(length(E_local)&&cross_limit>0L) {
      rho_extra<-lapply(seq_len(min(length(E_local),cross_limit)),function(j)
        optimize_one(E_local[[j]]$x,"rho",paste0("refinement_round_",round,
          "_E1_candidate_",j,"_to_rho")))
      rho_runs<-c(rho_runs,rho_extra);rho_local<-cluster(rho_runs,"rho")
    }
    if(length(rho_local)) {
      E_extra<-lapply(seq_along(rho_local),function(j)
        optimize_one(rho_local[[j]]$x,"E1",paste0("refinement_round_",round,
          "_rho_candidate_",j,"_to_E1")))
      E_runs<-c(E_runs,E_extra);E_local<-cluster(E_runs,"E1")
    }
  }
  if(!length(rho_local)||!length(E_local))stop("No feasible rho or E1 geometry candidates found.")
  spectra<-cp_jacobian_spectra(core$A,core$B,problem$C,problem$pairs)
  strip<-function(x){x$x<-NULL;x}
  best_rho<-rho_local[[1L]];best_E<-E_local[[1L]]
  # Every feasible point supplies a valid upper bound on a global infimum,
  # even when it was optimized for the other objective. Do not relabel such
  # a witness as a stationary point/local minimum of the other objective.
  rho_valid<-Filter(function(x)isTRUE(x$feasible),rho_runs)
  E_valid<-Filter(function(x)isTRUE(x$feasible),E_runs)
  witness_pool<-c(lapply(rho_valid,function(x){x$optimized_for<-"rho";x}),
                  lapply(E_valid,function(x){x$optimized_for<-"E1";x}))
  E_witness<-witness_pool[[which.min(vapply(witness_pool,function(x)x$E1,numeric(1L)))]]
  rho_witness<-witness_pool[[which.min(vapply(witness_pool,function(x)x$rho,numeric(1L)))]]
  dg<-list(rho_min=best_rho$rho,E_at_rho=best_rho$E1,E_min=best_E$E1,rho_at_E=best_E$rho,
    E1_feasible_upper_bound=E_witness$E1,rho_feasible_upper_bound=rho_witness$rho,
    E1_search_consistency_gap=max(0,best_E$E1-E_witness$E1),
    rho_search_consistency_gap=max(0,best_rho$rho-rho_witness$rho),
    E1_upper_bound_witness_source=E_witness$source,
    E1_upper_bound_witness_optimized_for=E_witness$optimized_for,
    rho_upper_bound_witness_source=rho_witness$source,
    rho_upper_bound_witness_optimized_for=rho_witness$optimized_for,
    feasible_witnesses=list(E1=strip(E_witness),rho=strip(rho_witness)),
    lambda_rho=best_rho$lambda,lambda_E=best_E$lambda,
    cosine_at_rho=best_rho$max_cosine,cosine_at_E=best_E$max_cosine,
    rho_best_success=best_rho$success,E_best_success=best_E$success,
    rho_local_minima=lapply(rho_local,strip),E_local_minima=lapply(E_local,strip),
    n_rho_local_minima=length(rho_local),n_E_local_minima=length(E_local),
    rho_search_runs=length(rho_runs),E_search_runs=length(E_runs),
    rho_search_feasible=sum(vapply(rho_runs,function(x)isTRUE(x$feasible),logical(1L))),
    E_search_feasible=sum(vapply(E_runs,function(x)isTRUE(x$feasible),logical(1L))),
    rho_search_converged=sum(vapply(rho_runs,function(x)isTRUE(x$success),logical(1L))),
    E_search_converged=sum(vapply(E_runs,function(x)isTRUE(x$success),logical(1L))),
    J_min=spectra$J_min,J_relative_min=spectra$J_relative_min,
    J_orig=spectra$orig,J_nor=spectra$nor,J_negative_roundoff=spectra$negative_roundoff,
    Z_singular_values=problem$Z_singular_values,
    C_singular_values=svd(problem$C,nu=0L,nv=0L)$d,core=core,
    configuration=list(seed=seed,random_starts=nrandom,truth_projection_starts=length(extra_lambdas),
      cosine_cap=problem$cosine_cap,maxeval=maxeval,ftol_abs=ftol,
      cluster_cosine_tolerance=cluster_tol,feasibility_tolerance=feasibility_tol,
      cross_start_limit=cross_limit,cross_refinement_rounds=cross_rounds,
      per_run_maxeval=maxeval,svd_convention="orig=s; nor=s/max(s); neither squared",
      feasible_bound_note="Bounds use all feasible rho/E1 runs; cross-objective witnesses are not relabelled as local minima.",
      minimum_note="Numerical candidates found; no automatic true/false-zero classification."))
  run_rows<-list()
  for(objective in c("rho","E1")) {
    runs<-if(objective=="rho")rho_runs else E_runs
    for(i in seq_along(runs)) {
      x<-runs[[i]]
      run_rows[[length(run_rows)+1L]]<-data.frame(objective=objective,run=i,source=x$source,
        feasible=x$feasible,converged=x$success,status=x$status,message=x$message,
        rho=x$rho,E1=x$E1,max_loading_cosine=x$max_cosine,matrix_norm=x$matrix_norm,
        eq_residual=x$eq_residual,iterations=x$iterations,stringsAsFactors=FALSE)
    }
  }
  list(diag=dg,all_search_runs=do.call(rbind,run_rows))
}

cp_grid_spectrum_rows <- function(M,name,rank_tolerance=1e-10,ell=NA_integer_,
                                  assigned_id=NA_integer_,basis_key="",reference=NULL) {
  singular<-sort(c(svd(M,nu=0L,nv=0L)$d,rep(0,max(0,ncol(M)-min(dim(M))))))
  relative<-if(max(singular)>0)singular/max(singular)else rep(NA_real_,length(singular))
  rank<-sum(singular>max(singular)*rank_tolerance)
  out<-data.frame(matrix_name=name,ell=ell,svd_index=seq_along(singular),orig=singular,nor=relative,
    gram_eigen_orig=singular^2,gram_eigen_nor=relative^2,numerical_rank=rank,
    n_svd=length(singular),n_rows=nrow(M),n_cols=ncol(M),
    assigned_true_direction_id=assigned_id,population_basis_key=basis_key,
    reference_orig=NA_real_,reference_nor=NA_real_,reference_gram_eigen_orig=NA_real_,
    reference_gram_eigen_nor=NA_real_,reference_numerical_rank=NA_integer_,stringsAsFactors=FALSE)
  if(!is.null(reference)&&nrow(reference)) {
    ref<-reference[match(out$svd_index,reference$svd_index),,drop=FALSE]
    out$reference_orig<-ref$orig;out$reference_nor<-ref$nor
    out$reference_gram_eigen_orig<-ref$gram_eigen_orig
    out$reference_gram_eigen_nor<-ref$gram_eigen_nor
    out$reference_numerical_rank<-ref$J_rank
  }
  out
}

cp_grid_candidate_pairs <- function(dg,problem) {
  R<-dg$rho_local_minima;E<-dg$E_local_minima
  lambda_angles<-matrix(NA_real_,length(R),length(E))
  matrix_angles<-lambda_angles
  for(i in seq_along(R))for(j in seq_along(E)) {
    lambda_angles[i,j]<-cp_alignment_angle_deg(R[[i]]$lambda,E[[j]]$lambda)
    matrix_angles[i,j]<-cp_alignment_angle_deg(cp_alignment_x_from_lambda(problem,R[[i]]$lambda),
      cp_alignment_x_from_lambda(problem,E[[j]]$lambda))
  }
  nearE<-apply(lambda_angles,1L,which.min);nearR<-apply(lambda_angles,2L,which.min)
  rows<-list()
  for(i in seq_along(R))for(j in seq_along(E))rows[[length(rows)+1L]]<-data.frame(
    rho_id=R[[i]]$minimum_id,E1_id=E[[j]]$minimum_id,
    lambda_angle_deg=lambda_angles[i,j],matrix_angle_deg=matrix_angles[i,j],
    rho_at_rho=R[[i]]$rho,E1_at_rho=R[[i]]$E1,rho_at_E1=E[[j]]$rho,E1_at_E1=E[[j]]$E1,
    E1_gap=R[[i]]$E1-E[[j]]$E1,is_nearest_E1_to_rho=nearE[i]==j,
    is_nearest_rho_to_E1=nearR[j]==i,is_rho_best=i==1L,is_E1_best=j==1L,
    rho_success=R[[i]]$success,E1_success=E[[j]]$success,stringsAsFactors=FALSE)
  do.call(rbind,rows)
}

cp_grid_annotate_candidates <- function(candidates,known,problem) {
  for(i in seq_along(candidates)) {
    point<-candidates[[i]]
    x<-cp_alignment_x_from_lambda(problem,point$lambda)
    physical_cosines<-vapply(known$points,function(p)abs(sum(x*p$x_projection)),numeric(1L))
    lambda_angles<-vapply(known$points,function(p)if(is.null(p$point))NA_real_ else
      cp_alignment_angle_deg(point$lambda,p$point$lambda),numeric(1L))
    nearest<-which.max(physical_cosines)
    candidates[[i]]$nearest_known_direction_id<-nearest
    candidates[[i]]$nearest_known_matrix_angle_deg<-acos(min(1,max(0,physical_cosines[nearest])))*180/pi
    nearest_lambda<-which.min(lambda_angles)
    candidates[[i]]$nearest_projected_known_lambda_id<-nearest_lambda
    candidates[[i]]$nearest_projected_known_lambda_angle_deg<-lambda_angles[nearest_lambda]
    candidates[[i]]$rankone_error_percent<-100*point$E1
  }
  candidates
}

cp_grid_geometry <- function(fit,design,cfg,seed) {
  if(!requireNamespace("nloptr",quietly=TRUE))stop("Package nloptr is required for grid geometry.")
  rank_tolerance<-cp_grid_get(cfg,"rank_tolerance",1e-10)
  cosine_cap<-cp_grid_get(cfg,"diagnostic_cosine_cap",.90)
  core<-cp_prepare_core(fit$A,fit$B,design$A,design$B,rank_tolerance=rank_tolerance)
  aligned<-cp_joint_alignment(fit$A,fit$B,design$A,design$B)
  qa<-cp_qr_core(aligned$A,rank_tolerance=rank_tolerance)
  qb<-cp_qr_core(aligned$B,rank_tolerance=rank_tolerance)
  stopifnot(max(abs(cp_unit_columns(qa$core)-core$A))<1e-10,
            max(abs(cp_unit_columns(qb$core)-core$B))<1e-10)
  problem<-cp_diagnostic_problem(core,cosine_cap,rank_tolerance)
  matching<-cp_grid_known_matching(core,aligned$A,aligned$B,design)
  known<-cp_grid_project_known(problem,qa,qb,matching,cosine_cap)
  search<-cp_grid_search_geometry(core,problem,cfg,seed,known$extra_lambdas,known$extra_labels)
  dg<-search$diag
  dg$rho_local_minima<-cp_grid_annotate_candidates(dg$rho_local_minima,known,problem)
  dg$E_local_minima<-cp_grid_annotate_candidates(dg$E_local_minima,known,problem)
  n<-if(is.null(design$n)||!is.finite(design$n)).5 else design$n
  all_reference<-cp_grid_population_bases(design$regime,n,rank_tolerance)
  basis_key<-matching$rows$population_basis_key[1L]
  reference<-all_reference[all_reference$population_basis_key==basis_key,,drop=FALSE]
  reference$ell<-match(reference$assigned_true_direction_id,matching$assignment)
  reference<-reference[order(reference$ell,reference$svd_index),setdiff(names(reference),"regime"),drop=FALSE]
  spectra<-list(cp_grid_spectrum_rows(core$A,"Ar",rank_tolerance,basis_key=basis_key),
    cp_grid_spectrum_rows(core$B,"Br",rank_tolerance,basis_key=basis_key),
    cp_grid_spectrum_rows(problem$Z,"Z",rank_tolerance,basis_key=basis_key),
    cp_grid_spectrum_rows(problem$C,"C",rank_tolerance,basis_key=basis_key))
  for(ell in seq_len(problem$d)) {
    incident<-which(colSums(problem$pairs==ell)>0L)
    spectra[[length(spectra)+1L]]<-cp_grid_spectrum_rows(problem$C[,incident,drop=FALSE],"J",rank_tolerance,
      ell=ell,assigned_id=matching$assignment[ell],basis_key=basis_key,
      reference=reference[reference$ell==ell,,drop=FALSE])
  }
  fake_checkpoint<-list(regime=design$regime,seed=design$seed,
    T=if(is.null(fit$T))NA_integer_ else fit$T,fit=fit,diagnostics=dg,Ar=core$A,Br=core$B)
  alignment<-cp_local_alignment_case(fake_checkpoint,
    angle_deg=cp_grid_get(cfg,"alignment_angle_deg",5),
    cap_cosines=cp_grid_get(cfg,"alignment_cosine_caps",c(.95,.90)),
    nstarts=cp_grid_get(cfg,"alignment_starts",5L),
    maxeval=cp_grid_get(cfg,"alignment_maxeval",300L),seed=seed+1907L)
  alignment_rows<-alignment$rows[,setdiff(names(alignment$rows),c("regime","seed","T")),drop=FALSE]
  details<-if(isTRUE(cp_grid_get(cfg,"save_search_details",FALSE)))alignment$details else NULL
  list(diag=dg,alignment_rows=alignment_rows,alignment_details=details,
    pair_distances=cp_grid_candidate_pairs(dg,problem),known_directions=known$rows,
    direction_matching=matching$rows,spectra_extra=do.call(rbind,spectra),
    population_reference=reference,all_search_runs=search$all_search_runs,
    configuration=list(angle_coordinates="lambda vector projective angle; sign ignored",
      primary_cosine_cap=cosine_cap,uses_truth_projection_starts=TRUE,
      random_starts=cp_grid_get(cfg,"diagnostic_random_starts",30L),
      truth_projection_starts=length(known$extra_lambdas),
      gram_eigenvalue_note="gram_eigen_orig=orig^2; gram_eigen_nor=nor^2; orig/nor remain unsquared SVD",
      classification_note="No automatic fake-zero label is assigned from angular disagreement."))
}


# Streaming reports for the p/q grid experiment. Only one checkpoint is retained.
cp_grid_bind_missing <- function(frames) {
  frames <- Filter(function(x)is.data.frame(x) && nrow(x)>0L,frames)
  if(!length(frames))return(data.frame())
  columns <- unique(unlist(lapply(frames,names),use.names=FALSE))
  frames <- lapply(frames,function(x) {
    for(column in setdiff(columns,names(x)))x[[column]]<-NA
    x[,columns,drop=FALSE]
  })
  do.call(rbind,frames)
}

cp_grid_safe_mean <- function(x) {
  x<-as.numeric(x);x<-x[is.finite(x)]
  if(length(x))mean(x)else NA_real_
}

cp_grid_safe_var <- function(x) {
  x<-as.numeric(x);x<-x[is.finite(x)]
  if(length(x)>1L)stats::var(x)else NA_real_
}

cp_grid_report_scalar <- function(x, fallback=NA_real_) {
  if(length(x)==1L)x else fallback
}

cp_grid_report_metadata <- function(frame, record) {
  if(is.null(frame) || !is.data.frame(frame) || !nrow(frame))return(data.frame())
  # Strip embedded newlines so each CSV record occupies one physical line.
  for(nm in names(frame))if(is.character(frame[[nm]]))
    frame[[nm]]<-gsub("[\r\n]+"," ",frame[[nm]])
  for(nm in c("p","q","regime","seed","T"))frame[[nm]]<-record[[nm]]
  frame$qmle_converged<-isTRUE(record$fit$converged)
  frame$checkpoint_status<-if(is.null(record$status))"unknown"else record$status
  first<-c("p","q","regime","seed","T","qmle_converged","checkpoint_status")
  frame[,c(first,setdiff(names(frame),first)),drop=FALSE]
}

cp_grid_report_append <- function(state, name, frame) {
  if(is.null(frame) || !is.data.frame(frame) || !nrow(frame))return(invisible(NULL))
  for(nm in names(frame))if(is.character(frame[[nm]]))
    frame[[nm]]<-gsub("[\r\n]+"," ",frame[[nm]])
  path<-file.path(state$directory,name)
  known<-state$schemas[[name]]
  numeric_now<-names(frame)[vapply(frame,function(x)is.numeric(x)||is.logical(x),logical(1L))]
  state$numeric[[name]]<-union(state$numeric[[name]],numeric_now)
  if(is.null(known)) {
    state$schemas[[name]]<-names(frame)
    utils::write.table(frame,path,sep=",",row.names=FALSE,col.names=TRUE,na="",quote=TRUE)
  } else {
    added<-setdiff(names(frame),known)
    if(length(added)) {
      # Rare schema extension (for example d=5 followed by d=8): append
      # empty cells to existing physical records without loading the CSV.
      columns<-c(known,added)
      temporary<-paste0(path,".schema.tmp")
      empty<-as.data.frame(setNames(rep(list(character()),length(columns)),columns),check.names=FALSE)
      utils::write.table(empty,temporary,sep=",",row.names=FALSE,col.names=TRUE,na="",quote=TRUE)
      input<-file(path,open="rt",encoding="UTF-8")
      output<-file(temporary,open="at",encoding="UTF-8")
      readLines(input,n=1L,warn=FALSE)
      suffix<-paste(rep(",",length(added)),collapse="")
      repeat {
        lines<-readLines(input,n=20000L,warn=FALSE)
        if(!length(lines))break
        writeLines(paste0(lines,suffix),output,useBytes=TRUE)
      }
      close(input);close(output)
      if(!file.copy(temporary,path,overwrite=TRUE))stop("Cannot extend CSV schema: ",path)
      unlink(temporary)
      known<-columns;state$schemas[[name]]<-columns
    }
    for(nm in setdiff(known,names(frame)))frame[[nm]]<-NA
    utils::write.table(frame[,known,drop=FALSE],path,sep=",",row.names=FALSE,
      col.names=FALSE,append=TRUE,na="",quote=TRUE)
  }
  state$counts[[name]]<-if(is.null(state$counts[[name]]))nrow(frame)else state$counts[[name]]+nrow(frame)
  invisible(NULL)
}

cp_grid_report_metrics <- function(r,cfg) {
  f<-r$fit;g<-r$geometry;dg<-g$diag
  if(is.null(dg))dg<-r$diagnostics
  core<-dg$core
  scalar<-cp_grid_report_scalar
  maxit<-if(!is.null(f$maxit_per_start))f$maxit_per_start else cfg$maxit_by_regime[[r$regime]]
  row<-data.frame(p=r$p,q=r$q,regime=r$regime,seed=r$seed,T=r$T,d=r$d,
    status=if(is.null(r$status))"unknown"else r$status,
    noise_estimated=if(is.null(cfg$estimate_noise))TRUE else isTRUE(cfg$estimate_noise),
    qmle_converged=isTRUE(f$converged),nll=scalar(f$nll),
    projected_gradient=scalar(f$projected_gradient_norm),gradient_norm=scalar(f$gradient_norm),
    maxit_per_start=scalar(maxit),qmle_starts=scalar(f$n_starts),
    rho_min=scalar(dg$rho_min),E1_at_rho=scalar(dg$E_at_rho),
    E1_min=scalar(dg$E_min),rho_at_E1=scalar(dg$rho_at_E),
    E1_feasible_upper_bound=scalar(dg$E1_feasible_upper_bound),
    rho_feasible_upper_bound=scalar(dg$rho_feasible_upper_bound),
    E1_search_consistency_gap=scalar(dg$E1_search_consistency_gap),
    rho_search_consistency_gap=scalar(dg$rho_search_consistency_gap),
    E1_at_rho_percent=100*scalar(dg$E_at_rho),E1_min_percent=100*scalar(dg$E_min),
    J_min=scalar(dg$J_min),J_relative_min=scalar(dg$J_relative_min),
    rho_best_success=isTRUE(dg$rho_best_success),E1_best_success=isTRUE(dg$E_best_success),
    rho_local_count=scalar(dg$n_rho_local_minima),E1_local_count=scalar(dg$n_E_local_minima),
    subspace_error=scalar(r$quality$subspace_error),
    covariance_relative_error=scalar(r$quality$covariance_relative_error),
    core_p=scalar(core$core_p),core_q=scalar(core$core_q),
    fitted_A_numerical_rank=scalar(core$numerical_rank_A),
    fitted_B_numerical_rank=scalar(core$numerical_rank_B),
    matching_cosine_mean=if(is.null(core$matching))NA_real_ else cp_grid_safe_mean(core$matching$matched_rankone_cosines),
    noise_at_lower_bound=scalar(f$noise_at_lower_bound),noise_at_upper_bound=scalar(f$noise_at_upper_bound),
    seconds=scalar(r$seconds),fit_seconds=scalar(f$elapsed_seconds),
    error_message=if(is.null(r$error_message))""else r$error_message,stringsAsFactors=FALSE)
  matching<-g$direction_matching
  if(is.data.frame(matching) && nrow(matching)) {
    for(nm in intersect(c("matching_to_original_mean","matching_to_known_set_mean",
      "matching_to_known_set_min","missing_true_direction_id","extra_selected",
      "is_original_basis_set","near_original_basis","n_known_true_directions"),names(matching)))
      row[[nm]]<-matching[[nm]][1L]
  }
  if(!is.null(f$M)) {
    ev<-eigen((f$M+t(f$M))/2,symmetric=TRUE,only.values=TRUE)$values
    row$Mx_min_eigenvalue<-min(ev);row$Mx_max_eigenvalue<-max(ev)
    row$Mx_condition_number<-if(min(ev)>0)max(ev)/min(ev)else Inf
  }
  row
}

cp_grid_report_local_minima <- function(r) {
  dg<-r$geometry$diag
  if(is.null(dg))return(data.frame())
  rows<-list()
  d<-if(is.null(r$d))ncol(dg$core$A)else r$d
  for(objective in c("rho","E1")) {
    points<-if(objective=="rho")dg$rho_local_minima else dg$E_local_minima
    for(j in seq_along(points)) {
      point<-points[[j]]
      row<-data.frame(objective=objective,local_minimum_order=if(is.null(point$minimum_id))j else point$minimum_id,
        objective_value=point$objective_value,rho=point$rho,E1=point$E1,
        rankone_error_percent=100*point$E1,max_loading_cosine=point$max_cosine,
        matrix_norm=point$matrix_norm,hit_count=point$hit_count,
        converged_hits=point$converged_hit_count,solver_converged=isTRUE(point$success),
        solver_status=point$status,solver_message=point$message,
        source=if(is.null(point$source))""else point$source,stringsAsFactors=FALSE)
      for(nm in setdiff(names(point),names(row))) {
        value<-point[[nm]]
        if(length(value)==1L && (is.numeric(value)||is.logical(value)||is.character(value)))row[[nm]]<-value
      }
      if(!is.null(point$all_hit_sources))row$all_hit_sources<-paste(point$all_hit_sources,collapse=";")
      unit<-point$lambda/sqrt(sum(point$lambda^2))
      for(k in seq_len(max(8L,d))) {
        row[[sprintf("lambda_%02d",k)]]<-if(k<=d)point$lambda[k]else NA_real_
        row[[sprintf("lambda_unit_%02d",k)]]<-if(k<=d)unit[k]else NA_real_
      }
      rows[[length(rows)+1L]]<-row
    }
  }
  cp_grid_report_metadata(cp_grid_bind_missing(rows),r)
}

cp_grid_report_J <- function(r) {
  dg<-r$geometry$diag
  if(is.null(dg))return(data.frame())
  orig<-dg$J_orig;nor<-dg$J_nor
  if(is.null(orig))return(data.frame())
  if(is.null(nor))nor<-sweep(orig,1L,apply(orig,1L,max),"/")
  rows<-list()
  extra<-r$geometry$spectra_extra
  for(ell in seq_len(nrow(orig)))for(j in seq_len(ncol(orig))) {
    row<-data.frame(ell=ell,svd_index=j,orig=orig[ell,j],nor=nor[ell,j],
      gram_eigen_orig=orig[ell,j]^2,gram_eigen_nor=nor[ell,j]^2)
    if(is.data.frame(extra) && nrow(extra) && "matrix_name" %in% names(extra)) {
      ref<-extra[extra$matrix_name=="J" & !is.na(extra$ell) & extra$ell==ell & extra$svd_index==j,,drop=FALSE]
      if(nrow(ref))for(nm in intersect(c("assigned_true_direction_id","population_basis_key",
        "reference_orig","reference_nor","reference_gram_eigen_orig","reference_gram_eigen_nor",
        "reference_numerical_rank","numerical_rank"),names(ref)))row[[nm]]<-ref[[nm]][1L]
    }
    rows[[length(rows)+1L]]<-row
  }
  cp_grid_report_metadata(cp_grid_bind_missing(rows),r)
}

cp_grid_report_fit_spectra <- function(r) {
  fit<-r$fit
  if(is.null(fit$M))return(data.frame())
  eigenvalues<-sort(eigen((fit$M+t(fit$M))/2,symmetric=TRUE,only.values=TRUE)$values)
  signal<-fit$signal_singular_values
  if(is.null(signal) && !is.null(fit$A) && !is.null(fit$B)) {
    Z<-vapply(seq_len(ncol(fit$A)),function(j)as.vector(outer(fit$A[,j],fit$B[,j])),
      numeric(nrow(fit$A)*nrow(fit$B)))
    L<-tryCatch(t(chol(fit$M)),error=function(e)NULL)
    if(!is.null(L))signal<-svd(Z%*%L,nu=0L,nv=0L)$d
  }
  signal<-sort(signal)
  n<-max(length(eigenvalues),length(signal))
  rows<-data.frame(spectrum_index=seq_len(n),
    Mx_eigenvalue=c(eigenvalues,rep(NA_real_,n-length(eigenvalues))),
    signal_singular_value=c(signal,rep(NA_real_,n-length(signal))))
  rows$signal_singular_value_nor<-rows$signal_singular_value/sqrt(r$p*r$q)
  cp_grid_report_metadata(rows,r)
}

cp_grid_report_numeric <- function(x) {
  x[x=="TRUE"]<-"1";x[x=="FALSE"]<-"0"
  suppressWarnings(as.numeric(x))
}

cp_grid_report_summarize <- function(path,output_path,keys,fields,cfg,
                                     wide=FALSE,filter_rho_best=FALSE) {
  if(!file.exists(path))return(invisible(0L))
  header<-names(utils::read.csv(path,nrows=0L,check.names=FALSE))
  keys<-intersect(keys,header);fields<-intersect(fields,header)
  if(!length(fields))return(invisible(0L))
  connection<-file(path,open="rt",encoding="UTF-8")
  on.exit(close(connection),add=TRUE)
  readLines(connection,n=1L,warn=FALSE)
  accumulator<-new.env(hash=TRUE,parent=emptyenv())
  chunk_size<-if(is.null(cfg$report_chunk_rows))20000L else as.integer(cfg$report_chunk_rows)
  repeat {
    lines<-readLines(connection,n=chunk_size,warn=FALSE)
    if(!length(lines))break
    frame<-utils::read.csv(text=paste(lines,collapse="\n"),header=FALSE,col.names=header,
      colClasses="character",check.names=FALSE,na.strings=c("","NA"),stringsAsFactors=FALSE)
    if(filter_rho_best && "rho_order" %in% names(frame))frame<-frame[frame$rho_order=="1",,drop=FALSE]
    if(!nrow(frame))next
    for(subset in c("all_valid","qmle_converged")) {
      x<-if(subset=="qmle_converged")frame[!is.na(frame$qmle_converged) & frame$qmle_converged=="TRUE",,drop=FALSE]else frame
      if(!nrow(x))next
      key_values<-lapply(x[,keys,drop=FALSE],function(v)ifelse(is.na(v),"<NA>",v))
      group_keys<-do.call(paste,c(key_values,list(sep="\034")))
      indices<-split(seq_len(nrow(x)),group_keys)
      for(group_key in names(indices)) {
        ii<-indices[[group_key]];key<-paste(subset,group_key,sep="\034")
        state<-accumulator[[key]]
        if(is.null(state))state<-list(keys=as.list(x[ii[1L],keys,drop=FALSE]),subset=subset,
          n_saved=0L,n_qmle_converged=0L,n_diagnostic_ok=0L,values=list())
        state$n_saved<-state$n_saved+length(ii)
        state$n_qmle_converged<-state$n_qmle_converged+sum(x$qmle_converged[ii]=="TRUE",na.rm=TRUE)
        if("checkpoint_status" %in% names(x))
          state$n_diagnostic_ok<-state$n_diagnostic_ok+sum(x$checkpoint_status[ii]=="ok",na.rm=TRUE)
        else if("status" %in% names(x))
          state$n_diagnostic_ok<-state$n_diagnostic_ok+sum(x$status[ii]=="ok",na.rm=TRUE)
        for(field in fields) {
          value<-cp_grid_report_numeric(x[[field]][ii]);value<-value[is.finite(value)]
          n<-length(value);if(!n)next
          block_mean<-mean(value);block_M2<-sum((value-block_mean)^2)
          old<-state$values[[field]]
          if(is.null(old))state$values[[field]]<-list(n=n,mean=block_mean,M2=block_M2)
          else {
            delta<-block_mean-old$mean;total<-old$n+n
            state$values[[field]]<-list(n=total,mean=old$mean+delta*n/total,
              M2=old$M2+block_M2+delta^2*old$n*n/total)
          }
        }
        accumulator[[key]]<-state
      }
    }
  }
  rows<-list()
  for(key in ls(accumulator,all.names=TRUE)) {
    state<-accumulator[[key]]
    base<-as.data.frame(state$keys,stringsAsFactors=FALSE,check.names=FALSE)
    base$subset<-state$subset;base$n_saved<-state$n_saved;base$n_seeds_present<-state$n_saved
    base$n_requested<-length(cfg$seeds);base$n_qmle_converged<-state$n_qmle_converged
    base$n_diagnostic_ok<-state$n_diagnostic_ok
    if(wide) {
      row<-base
      for(field in fields) {
        value<-state$values[[field]]
        row[[paste0(field,"_mean")]]<-if(is.null(value))NA_real_ else value$mean
        row[[paste0(field,"_variance")]]<-if(is.null(value)||value$n<2L)NA_real_ else value$M2/(value$n-1L)
        row[[paste0(field,"_n_finite")]]<-if(is.null(value))0L else value$n
      }
      rows[[length(rows)+1L]]<-row
    } else for(field in fields) {
      row<-base;value<-state$values[[field]];row$metric<-field
      row$mean<-if(is.null(value))NA_real_ else value$mean
      row$variance<-if(is.null(value)||value$n<2L)NA_real_ else value$M2/(value$n-1L)
      row$n_finite<-if(is.null(value))0L else value$n
      rows[[length(rows)+1L]]<-row
    }
  }
  result<-cp_grid_bind_missing(rows)
  if(nrow(result))utils::write.csv(result,output_path,row.names=FALSE,na="")
  invisible(nrow(result))
}

cp_grid_write_reports <- function(cfg) {
  if(is.null(cfg$signature))cfg$signature<-cp_grid_config_signature(cfg)
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  state<-new.env(parent=emptyenv());state$directory<-cfg$output_dir
  state$schemas<-list();state$numeric<-list();state$counts<-list()
  filenames<-c("metrics_by_seed.csv","J_svd_by_seed.csv","J_gram_eigenvalues_by_seed.csv",
    "local_minima_by_seed.csv","minima_pair_angles_by_seed.csv","rho_neighborhood_alignment_by_seed.csv",
    "known_rankone_directions_by_seed.csv","direction_matching_by_seed.csv","all_matrix_spectra_by_seed.csv",
    "population_basis_reference_by_seed.csv","diagnostic_search_runs_by_seed.csv","qmle_all_starts.csv",
    "nonconvergence_and_errors.csv","fit_covariance_spectra_by_seed.csv",
    "checkpoint_read_errors.csv","checkpoint_signature_mismatches.csv")
  # Only replace generated report files; checkpoints are never modified here.
  for(name in filenames) {
    path<-file.path(cfg$output_dir,name)
    if(file.exists(path))unlink(path)
  }
  p_values<-cfg$p_values;q_values<-cfg$q_values
  n_read<-0L;n_signature_mismatch<-0L;n_read_error<-0L
  for(p in p_values)for(q in q_values)for(regime in cfg$regimes)for(T in cfg$sample_sizes)for(seed in cfg$seeds) {
    path<-cp_grid_checkpoint_file(cfg,p,q,regime,seed,T)
    if(!file.exists(path))next
    r<-tryCatch(readRDS(path),error=function(e)e)
    if(inherits(r,"error")) {
      n_read_error<-n_read_error+1L
      cp_grid_report_append(state,"checkpoint_read_errors.csv",data.frame(p=p,q=q,regime=regime,seed=seed,T=T,
        path=path,error_message=conditionMessage(r),stringsAsFactors=FALSE))
      next
    }
    if(!identical(r$config_signature,cfg$signature)) {
      n_signature_mismatch<-n_signature_mismatch+1L
      cp_grid_report_append(state,"checkpoint_signature_mismatches.csv",data.frame(p=p,q=q,regime=regime,seed=seed,T=T,
        path=path,stored_signature=if(is.null(r$config_signature))""else r$config_signature,
        expected_signature=cfg$signature,stringsAsFactors=FALSE))
      next
    }
    n_read<-n_read+1L
    for(nm in c("p","q","regime","seed","T"))if(is.null(r[[nm]]))r[[nm]]<-get(nm)
    if(is.null(r$geometry))r$geometry<-list(diag=r$diagnostics)
    metrics<-cp_grid_report_metrics(r,cfg)
    cp_grid_report_append(state,"metrics_by_seed.csv",metrics)
    if(!identical(r$status,"ok") || !isTRUE(r$fit$converged))
      cp_grid_report_append(state,"nonconvergence_and_errors.csv",metrics)
    cp_grid_report_append(state,"qmle_all_starts.csv",cp_grid_report_metadata(r$fit$all_starts,r))
    cp_grid_report_append(state,"fit_covariance_spectra_by_seed.csv",cp_grid_report_fit_spectra(r))
    J<-cp_grid_report_J(r)
    cp_grid_report_append(state,"J_svd_by_seed.csv",J)
    if(nrow(J)) {
      gram<-J[,setdiff(names(J),c("orig","nor")),drop=FALSE]
      cp_grid_report_append(state,"J_gram_eigenvalues_by_seed.csv",gram)
    }
    cp_grid_report_append(state,"local_minima_by_seed.csv",cp_grid_report_local_minima(r))
    tables<-c(alignment_rows="rho_neighborhood_alignment_by_seed.csv",
      pair_distances="minima_pair_angles_by_seed.csv",known_directions="known_rankone_directions_by_seed.csv",
      direction_matching="direction_matching_by_seed.csv",spectra_extra="all_matrix_spectra_by_seed.csv",
      population_reference="population_basis_reference_by_seed.csv",all_search_runs="diagnostic_search_runs_by_seed.csv")
    for(field in names(tables))cp_grid_report_append(state,tables[[field]],cp_grid_report_metadata(r$geometry[[field]],r))
    rm(r)
  }
  if(!n_read)stop("No compatible grid checkpoint was available for reporting.")
  base_keys<-c("p","q","regime","T")
  summarize<-function(input,output,keys,fields=NULL,wide=FALSE,best=FALSE) {
    if(is.null(fields)) {
      fields<-setdiff(state$numeric[[input]],c("p","q","seed","T","d","ell","svd_index",
        "rho_id","E1_id","rho_order","angle_deg","cap_cosine","true_direction_id",
        "assigned_true_direction_id","missing_true_direction_id","nearest_estimated_ell",
        "assigned_estimated_ell","E1_best_id","E1_nearest_center_id","E1_nearest_cap_id",
        "n_known_true_directions",keys))
      fields<-fields[!grepl("^(lambda_|lambda_unit_|center_lambda_|true_lambda_|coef_)",fields)]
    }
    cp_grid_report_summarize(file.path(cfg$output_dir,input),file.path(cfg$output_dir,output),
      keys,fields,cfg,wide=wide,filter_rho_best=best)
  }
  summarize("metrics_by_seed.csv","metrics_mean_variance.csv",base_keys)
  summarize("fit_covariance_spectra_by_seed.csv","fit_covariance_spectra_mean_variance.csv",
    c(base_keys,"spectrum_index"),c("Mx_eigenvalue","signal_singular_value","signal_singular_value_nor"),wide=TRUE)
  summarize("J_svd_by_seed.csv","J_svd_mean_variance.csv",c(base_keys,"ell","svd_index"),
    c("orig","nor","reference_orig","reference_nor"),wide=TRUE)
  summarize("J_svd_by_seed.csv","J_svd_by_basis_mean_variance.csv",
    c(base_keys,"population_basis_key","ell","assigned_true_direction_id","svd_index"),
    c("orig","nor","reference_orig","reference_nor"),wide=TRUE)
  summarize("J_gram_eigenvalues_by_seed.csv","J_gram_eigenvalues_mean_variance.csv",
    c(base_keys,"ell","svd_index"),c("gram_eigen_orig","gram_eigen_nor",
      "reference_gram_eigen_orig","reference_gram_eigen_nor"),wide=TRUE)
  summarize("local_minima_by_seed.csv","local_minima_mean_variance.csv",
    c(base_keys,"objective","local_minimum_order"),
    c("objective_value","rho","E1","rankone_error_percent","max_loading_cosine",
      "hit_count","converged_hits","solver_converged","truth_guided_hit",
      "nearest_known_matrix_angle_deg","nearest_projected_known_lambda_angle_deg"),wide=TRUE)
  summarize("rho_neighborhood_alignment_by_seed.csv","rho_neighborhood_alignment_mean_variance.csv",
    c(base_keys,"rho_order","cap_cosine"))
  summarize("rho_neighborhood_alignment_by_seed.csv","rho_neighborhood_best_mean_variance.csv",
    c(base_keys,"cap_cosine"),best=TRUE)
  summarize("known_rankone_directions_by_seed.csv","known_rankone_directions_mean_variance.csv",
    c(base_keys,"true_direction_id","direction_type"))
  summarize("direction_matching_by_seed.csv","direction_matching_mean_variance.csv",
    c(base_keys,"ell","assigned_true_direction_id"))
  detailed<-is.null(cfg$report_detailed_summaries)||isTRUE(cfg$report_detailed_summaries)
  if(detailed) {
    summarize("minima_pair_angles_by_seed.csv","minima_pair_angles_mean_variance.csv",
      c(base_keys,"rho_id","E1_id"))
    summarize("all_matrix_spectra_by_seed.csv","all_matrix_spectra_mean_variance.csv",
      c(base_keys,"matrix_name","ell","svd_index"),
      c("orig","nor","gram_eigen_orig","gram_eigen_nor","reference_orig","reference_nor",
        "reference_gram_eigen_orig","reference_gram_eigen_nor","numerical_rank",
        "reference_numerical_rank"),wide=TRUE)
    summarize("population_basis_reference_by_seed.csv","population_basis_reference_mean_variance.csv",
      c(base_keys,"population_basis_key","ell","assigned_true_direction_id","svd_index"),
      c("orig","nor","gram_eigen_orig","gram_eigen_nor","J_rank","basis_Z_min_svd"),wide=TRUE)
  }
  manifest<-list(created_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),
    n_checkpoints_read=n_read,n_signature_mismatch=n_signature_mismatch,n_read_errors=n_read_error,
    n_requested=length(p_values)*length(q_values)*length(cfg$regimes)*length(cfg$sample_sizes)*length(cfg$seeds),
    raw_row_counts=state$counts,raw_columns=state$schemas,
    variance_convention="sample variance across seeds, denominator n-1",
    local_order_note="k-th smallest candidate within each seed, not a matched cross-seed branch",
    classification_note="Angles and minima are recorded; no candidate is automatically labelled a fake zero.")
  saveRDS(manifest,file.path(cfg$output_dir,"report_manifest.rds"))
  writeLines(capture.output(str(manifest,max.level=2L)),file.path(cfg$output_dir,"report_manifest.txt"))
  message("Grid reports completed: ",cfg$output_dir," (",n_read," compatible checkpoints).")
  invisible(manifest)
}


# Grid orchestration. All source data are generated; each result has its own file.
cp_grid_checkpoint_file <- function(cfg,p,q,regime,seed,T) {
  file.path(cfg$output_dir,"checkpoints",sprintf("p%03d_q%03d",p,q),regime,
            sprintf("seed_%03d",seed),paste0("T",T,".rds"))
}

cp_grid_atomic_rds <- function(value,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- tempfile(pattern=".writing-",tmpdir=dirname(path),fileext=".rds")
  on.exit(if(file.exists(temporary))unlink(temporary),add=TRUE)
  saveRDS(value,temporary,compress=TRUE)
  if(file.exists(path) && !file.remove(path))stop("Cannot replace result file: ",path)
  if(!file.rename(temporary,path)) {
    if(!file.copy(temporary,path,overwrite=TRUE))stop("Cannot save result: ",path)
    unlink(temporary)
  }
  invisible(path)
}

cp_grid_config_signature <- function(cfg) {
  fields <- c("maxit_by_regime","truth_starts","truth_angle_deg","truth_chol_sd",
    "truth_noise_log_sd","estimate_noise","noise_bounds","factr","pgtol","lmm",
    "diagnostic_random_starts","diagnostic_maxeval","diagnostic_cosine_cap",
    "diagnostic_cross_start_limit","diagnostic_cross_rounds","alignment_angle_deg","alignment_cosine_caps",
    "alignment_starts","alignment_maxeval","rank_tolerance","n_exterior",
    "factor_multiplier","burn","save_search_details","save_all_start_fits")
  payload <- list(method_version="CP_grid_truth_near_2026-10-27_v1",
    parameters=cfg[fields],
    designs=lapply(c("interior","boundary","exterior"),core_loading_design,n=cfg$n_exterior),
    generating_functions=lapply(list(with_local_seed,seed_key,haar_thin_basis,make_loading_design,
      stable_random_matrix,stationary_varma_covariance,simulate_varma,generate_data),
      function(f)paste(deparse(body(f)),collapse="\n")))
  path<-tempfile(fileext=".rds")
  on.exit(unlink(path),add=TRUE)
  saveRDS(payload,path,compress=FALSE)
  unname(tools::md5sum(path))
}

cp_grid_quality <- function(fit,design) {
  Z<-cp_khatri_rao(fit$A,fit$B);Z0<-design$Z;d<-design$d
  orth<-function(X){f<-qr(X,tol=1e-10);qr.Q(f)[,seq_len(f$rank),drop=FALSE]}
  Q<-orth(Z);Q0<-orth(Z0)
  space<-sqrt(max(0,(ncol(Q)+ncol(Q0)-2*sum(crossprod(Q,Q0)^2))/(2*d)))
  loading_chol<-if(!is.null(fit$L))fit$L else t(chol(fit$M))
  F<-Z%*%loading_chol;F0<-Z0%*%t(chol(design$Mx))
  diagonal<-rowSums(F^2);diagonal0<-rowSums(F0^2)
  delta_noise<-fit$psi-1
  numerator<-sum(crossprod(F)^2)+sum(crossprod(F0)^2)-2*sum(crossprod(F,F0)^2)+
    sum(delta_noise^2)+2*sum(delta_noise*(diagonal-diagonal0))
  denominator<-sum(crossprod(F0)^2)+nrow(F0)+2*sum(diagonal0)
  list(subspace_error=space,covariance_relative_error=sqrt(max(0,numerator)/denominator),
    estimated_Z_rank=ncol(Q),noise_mean=mean(fit$psi),noise_variance_mse=mean((fit$psi-1)^2))
}

cp_grid_read_checkpoint <- function(path) {
  if(!file.exists(path))return(NULL)
  tryCatch(readRDS(path),error=function(e)NULL)
}

cp_grid_run_job <- function(job,cfg) {
  p<-as.integer(job$p);q<-as.integer(job$q);seed<-as.integer(job$seed)
  regime<-as.character(job$regime)
  id<-sprintf("p%03d_q%03d_%s_seed%03d",p,q,regime,seed)
  progress_file<-file.path(cfg$output_dir,"job_progress",paste0(id,".rds"))
  counts<-list(p=p,q=q,regime=regime,seed=seed,completed=0L,errors=0L,reused=0L)
  saved<-lapply(cfg$sample_sizes,function(T)cp_grid_read_checkpoint(
    cp_grid_checkpoint_file(cfg,p,q,regime,seed,T)))
  for(old in saved)if(!is.null(old) && !identical(old$config_signature,cfg$signature))
    stop("A checkpoint uses a different design or budget. Choose a new output_dir: ",id)
  complete<-vapply(saved,function(old)!is.null(old)&&identical(old$status,"ok"),logical(1L))
  if(cfg$resume && all(complete)) {
    counts$completed<-length(cfg$sample_sizes);counts$reused<-counts$completed
    cp_grid_atomic_rds(counts,progress_file)
    return(as.data.frame(counts,stringsAsFactors=FALSE))
  }
  if(cfg$verbose)message("[job] ",id)
  design<-generate_data(regime,seed,Tmax=max(cfg$sample_sizes),p=p,q=q,
    burn=cfg$burn,factor_multiplier=cfg$factor_multiplier,n=cfg$n_exterior)
  # Only Y is needed for fitting. Do not retain another full noise matrix per worker.
  design$noise<-NULL;design$X<-NULL
  metadata_path<-file.path(dirname(cp_grid_checkpoint_file(cfg,p,q,regime,seed,cfg$sample_sizes[1L])),
                           "design_metadata.rds")
  cp_grid_atomic_rds(list(metadata=design$metadata,Mx=design$Mx,
    raw_A=design$raw_A,raw_B=design$raw_B,extra_lambda=design$extra_lambda),metadata_path)
  previous_fit<-NULL
  for(i in seq_along(cfg$sample_sizes)) {
    T<-cfg$sample_sizes[i];old<-saved[[i]]
    path<-cp_grid_checkpoint_file(cfg,p,q,regime,seed,T)
    if(cfg$resume && !is.null(old) && identical(old$status,"ok")) {
      previous_fit<-old$fit;counts$completed<-counts$completed+1L
      counts$reused<-counts$reused+1L
      cp_grid_atomic_rds(counts,progress_file)
      next
    }
    tick<-proc.time()[[3L]]
    result<-list(p=p,q=q,regime=regime,seed=seed,T=T,d=design$d,
      config_signature=cfg$signature,status="fit_error",error_message="",
      truth_metadata=design$metadata,truth_Mx=design$Mx,reused_fit=FALSE)
    if(cfg$verbose)message(sprintf("[%s p=%d q=%d] seed=%d T=%d",regime,p,q,seed,T))
    if(cfg$resume && !is.null(old) && !is.null(old$fit)) {
      fit<-old$fit;result$reused_fit<-TRUE
    } else {
      Y<-design$Y[seq_len(T),,drop=FALSE]
      fit_seed<-as.integer(seed+17L*T+100003L*match(regime,simulation_regimes)+1009L*p+1013L*q)
      fit<-tryCatch(cp_grid_qmle(Y,design,cfg,seed=fit_seed,previous_fit=previous_fit),error=function(e)e)
      rm(Y)
    }
    if(inherits(fit,"error")) {
      result$error_message<-conditionMessage(fit)
      counts$errors<-counts$errors+1L
    } else {
      result$fit<-fit;previous_fit<-fit
      result$quality<-cp_grid_quality(fit,design)
      result$status<-"diagnostic_error"
      diagnostic_seed<-as.integer(seed+31L*T+500003L*match(regime,simulation_regimes)+10007L*p+10009L*q)
      geometry<-tryCatch(cp_grid_geometry(fit,design,cfg,seed=diagnostic_seed),error=function(e)e)
      if(inherits(geometry,"error")) {
        result$error_message<-conditionMessage(geometry)
        counts$errors<-counts$errors+1L
      } else {
        result$geometry<-geometry;result$status<-"ok"
        result$Ar<-geometry$diag$core$A;result$Br<-geometry$diag$core$B
        stopifnot(max(abs(colSums(result$Ar^2)-1))<1e-8,
                  max(abs(colSums(result$Br^2)-1))<1e-8)
        counts$completed<-counts$completed+1L
      }
    }
    result$seconds<-proc.time()[[3L]]-tick
    result$completed_at<-format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")
    cp_grid_atomic_rds(result,path)
    cp_grid_atomic_rds(counts,progress_file)
    if(!identical(result$status,"ok"))message("  ",result$status,": ",result$error_message)
    rm(result,fit)
    gc(verbose=FALSE)
  }
  rm(design)
  gc(verbose=FALSE)
  as.data.frame(counts,stringsAsFactors=FALSE)
}

cp_grid_progress <- function(output_dir) {
  files<-list.files(file.path(output_dir,"job_progress"),pattern="\\.rds$",full.names=TRUE)
  values<-lapply(files,cp_grid_read_checkpoint)
  values<-Filter(Negate(is.null),values)
  if(!length(values))return(data.frame(jobs=0L,completed=0L,errors=0L,reused=0L))
  data.frame(jobs=length(values),completed=sum(vapply(values,function(x)x$completed,numeric(1L))),
    errors=sum(vapply(values,function(x)x$errors,numeric(1L))),
    reused=sum(vapply(values,function(x)x$reused,numeric(1L))))
}

cp_grid_write_population <- function(cfg) {
  output<-list();rows<-list();spectra<-list();minima<-list()
  for(regime in cfg$regimes) {
    truth<-core_loading_design(regime,cfg$n_exterior)
    dg<-cp_diagnostics(truth$Ar_unit,truth$Br_unit,
      seed=1819L,
      nstarts=max(8L,cfg$diagnostic_random_starts),
      cosine_cap=cfg$diagnostic_cosine_cap,maxeval=cfg$diagnostic_maxeval,
      extra_lambdas=truth$extra_lambda)
    output[[regime]]<-dg
    sv<-dg$C_singular_values
    rows[[regime]]<-data.frame(regime=regime,d=truth$d,
      rho_min=dg$rho_min,E1_at_rho=dg$E_at_rho,E1_min=dg$E_min,rho_at_E1=dg$rho_at_E,
      J_min=dg$J_min,J_relative_min=dg$J_relative_min,
      C_rank=sum(sv>max(sv)*cfg$rank_tolerance),C_columns=choose(truth$d,2L),
      rho_local_count=dg$n_rho_local_minima,E1_local_count=dg$n_E_local_minima)
    spectra[[regime]]<-cp_grid_population_bases(regime,n=cfg$n_exterior,
                                              rank_tolerance=cfg$rank_tolerance)
    for(objective in c("rho","E1")) {
      points<-if(objective=="rho")dg$rho_local_minima else dg$E_local_minima
      for(point in points) {
        value<-data.frame(regime=regime,objective=objective,
          local_minimum_order=point$minimum_id,rho=point$rho,E1=point$E1,
          max_loading_cosine=point$max_cosine,solver_converged=point$success,
          stringsAsFactors=FALSE)
        for(j in seq_len(8L))value[[sprintf("lambda_%02d",j)]]<-
          if(j<=length(point$lambda))point$lambda[j]else NA_real_
        minima[[length(minima)+1L]]<-value
      }
    }
  }
  cp_grid_atomic_rds(output,file.path(cfg$output_dir,"population_benchmarks.rds"))
  utils::write.csv(do.call(rbind,rows),file.path(cfg$output_dir,"population_benchmarks.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,spectra),file.path(cfg$output_dir,"population_J_all_bases.csv"),row.names=FALSE)
  utils::write.csv(do.call(rbind,minima),file.path(cfg$output_dir,"population_local_minima.csv"),row.names=FALSE)
  invisible(NULL)
}

cp_grid_make_cluster <- function(ncores,log,threads) {
  variables<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","LC_ALL")
  previous<-Sys.getenv(variables,unset=NA_character_)
  on.exit({
    absent<-is.na(previous)
    if(any(absent))Sys.unsetenv(variables[absent])
    if(any(!absent))do.call(Sys.setenv,as.list(setNames(previous[!absent],variables[!absent])))
  },add=TRUE)
  values<-c(rep(as.character(threads),3L),Sys.getlocale("LC_CTYPE"))
  do.call(Sys.setenv,as.list(setNames(values,variables)))
  parallel::makePSOCKcluster(ncores,outfile=log)
}

cp_run_grid <- function(cfg) {
  cp_ensure_dependencies(cfg)
  cfg$p_values<-sort(unique(as.integer(cfg$p_values)))
  cfg$q_values<-sort(unique(as.integer(cfg$q_values)))
  cfg$seeds<-sort(unique(as.integer(cfg$seeds)))
  cfg$sample_sizes<-sort(unique(as.integer(cfg$sample_sizes)))
  stopifnot(length(cfg$p_values)>0L,length(cfg$q_values)>0L,length(cfg$sample_sizes)>0L,
    length(cfg$regimes)>0L,all(cfg$p_values>=8L),all(cfg$q_values>=8L),all(cfg$sample_sizes>=2L),
    all(cfg$regimes %in% simulation_regimes),length(cfg$seeds)>0L,
    length(cfg$ncores)==1L,is.finite(cfg$ncores),cfg$ncores>=1L,
    all(cfg$regimes %in% names(cfg$maxit_by_regime)),
    all(cfg$maxit_by_regime[cfg$regimes]>=1L))
  cfg$output_dir<-normalizePath(cfg$output_dir,winslash="/",mustWork=FALSE)
  cfg$signature<-cp_grid_config_signature(cfg)
  config_path<-file.path(cfg$output_dir,"configuration.rds")
  existing<-cp_grid_read_checkpoint(config_path)
  if(!is.null(existing) && !identical(existing$signature,cfg$signature))
    stop("This output directory has another design/budget. Choose a new output_dir before running.")
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  cp_grid_atomic_rds(cfg,config_path)
  writeLines(capture.output(str(cfg)),file.path(cfg$output_dir,"configuration.txt"))
  writeLines(capture.output(sessionInfo()),file.path(cfg$output_dir,"sessionInfo.txt"))
  cp_grid_write_population(cfg)
  if(isTRUE(cfg$validate_gradients)) {
    message("Checking analytic QMLE gradients...")
    invisible(cp_check_qmle_gradient(estimate_noise=cfg$estimate_noise))
  }
  jobs<-expand.grid(seed=cfg$seeds,regime=cfg$regimes,p=cfg$p_values,q=cfg$q_values,
    KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  jobs<-jobs[order(-jobs$p*jobs$q,match(jobs$regime,cfg$regimes),jobs$p,jobs$q,jobs$seed),,drop=FALSE]
  job_list<-lapply(seq_len(nrow(jobs)),function(i)as.list(jobs[i,,drop=FALSE]))
  total<-nrow(jobs)*length(cfg$sample_sizes)
  message("Configurations: ",total,"; jobs: ",nrow(jobs),"; local workers: ",cfg$ncores)
  message("Checkpoints: ",file.path(cfg$output_dir,"checkpoints"))
  cluster<-NULL
  if(cfg$ncores>1L) {
    log<-file.path(cfg$output_dir,"worker_log.txt")
    message("Worker progress log: ",log)
    cluster<-cp_grid_make_cluster(cfg$ncores,log,cfg$blas_threads_per_worker)
    on.exit(if(!is.null(cluster))parallel::stopCluster(cluster),add=TRUE)
    parallel::clusterCall(cluster,function(paths,threads) {
      if(.Platform$OS.type=="windows" && !isTRUE(l10n_info()[["UTF-8"]]))
        invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE)))
      .libPaths(paths)
      Sys.setenv(OMP_NUM_THREADS=as.character(threads),OPENBLAS_NUM_THREADS=as.character(threads),
        MKL_NUM_THREADS=as.character(threads))
      stopifnot(requireNamespace("nloptr",quietly=TRUE));NULL
    },.libPaths(),cfg$blas_threads_per_worker)
    parallel::clusterExport(cluster,GRID_EXPORTS,envir=environment(cp_run_grid))
  }
  batch_size<-max(cfg$ncores,as.integer(cfg$batch_jobs))
  batches<-split(seq_along(job_list),ceiling(seq_along(job_list)/batch_size))
  all_jobs<-list()
  for(batch_index in seq_along(batches)) {
    indices<-batches[[batch_index]]
    answer<-if(is.null(cluster))lapply(job_list[indices],cp_grid_run_job,cfg=cfg)else
      parallel::parLapplyLB(cluster,job_list[indices],cp_grid_run_job,cfg=cfg)
    all_jobs<-c(all_jobs,answer)
    status<-do.call(rbind,all_jobs)
    utils::write.csv(status,file.path(cfg$output_dir,"job_status.csv"),row.names=FALSE)
    message(sprintf("Progress: %d/%d jobs; %d/%d configurations saved; errors=%d",
      nrow(status),nrow(jobs),sum(status$completed),total,sum(status$errors)))
  }
  if(!is.null(cluster)){parallel::stopCluster(cluster);cluster<-NULL}
  result<-cp_grid_write_reports(cfg)
  metric_path<-file.path(cfg$output_dir,"metrics_by_seed.csv")
  if(file.exists(metric_path)) {
    metrics<-utils::read.csv(metric_path,stringsAsFactors=FALSE)
    message("Selected QMLE convergence flags: ",sum(metrics$qmle_converged,na.rm=TRUE),
      "/",nrow(metrics),". See nonconvergence_and_errors.csv for the remaining cases.")
  }
  message("Completed grid. Results: ",cfg$output_dir)
  invisible(result)
}

cp_grid_main <- function() {
  cfg<-GRID_CONFIG
  args<-commandArgs(trailingOnly=TRUE)
  if("--smoke" %in% args) {
    cfg$p_values<-c(20L,40L);cfg$q_values<-c(20L,40L)
    cfg$seeds<-1L;cfg$sample_sizes<-c(100L,200L)
    cfg$maxit_by_regime<-c(interior=5L,boundary=5L,exterior=5L)
    cfg$truth_starts<-2L;cfg$diagnostic_random_starts<-4L
    cfg$diagnostic_maxeval<-40L;cfg$alignment_starts<-2L;cfg$alignment_maxeval<-40L
    cfg$ncores<-min(2L,cfg$ncores)
    cfg$output_dir<-paste0(cfg$output_dir,"_SMOKE")
    message("Smoke test only: 24 configurations with deliberately small optimization budgets.")
  }
  out<-grep("^--out=",args,value=TRUE)
  cores<-grep("^--cores=",args,value=TRUE)
  if(length(out))cfg$output_dir<-sub("^--out=","",out[1L])
  if(length(cores))cfg$ncores<-as.integer(sub("^--cores=","",cores[1L]))
  cfg$batch_jobs<-max(cfg$ncores,cfg$batch_jobs)
  if("--summary-only" %in% args) {
    saved<-readRDS(file.path(cfg$output_dir,"configuration.rds"))
    return(invisible(cp_grid_write_reports(saved)))
  }
  cp_run_grid(cfg)
}


GRID_EXPORTS<-c(
  "cp_ensure_dependencies",
  "normalize_columns",
  "khatri_rao",
  "with_local_seed",
  "seed_key",
  "haar_thin_basis",
  "core_loading_design",
  "make_loading_design",
  "stable_random_matrix",
  "stationary_varma_covariance",
  "simulate_varma",
  "generate_data",
  "centered_prefix",
  "normalized_qr_cores",
  "validate_R_design",
  "cp_normalize_columns",
  "cp_khatri_rao",
  "cp_qmle_objective",
  "cp_principal_subspace",
  "cp_rankone_candidate_pool",
  "cp_qmle",
  "cp_check_qmle_gradient",
  "cp_assert_loadings",
  "cp_unit_columns",
  "cp_max_assignment",
  "cp_joint_alignment",
  "cp_qr_core",
  "cp_prepare_core",
  "cp_minor_matrix",
  "cp_jacobian_spectra",
  "cp_diagnostic_problem",
  "cp_cluster_minima",
  "cp_diagnostics",
  "cp_alignment_angle_deg",
  "cp_alignment_change_cap",
  "cp_alignment_x_from_lambda",
  "cp_alignment_cap_constraints",
  "cp_alignment_cap_search",
  "cp_local_alignment_case",
  "cp_grid_option",
  "cp_grid_perturb_directions",
  "cp_grid_pack_exact",
  "cp_grid_initial_angles",
  "cp_grid_qmle",
  "cp_grid_get",
  "cp_grid_bind_frames",
  "cp_grid_true_core_directions",
  "cp_grid_population_bases",
  "cp_grid_known_matching",
  "cp_grid_project_known",
  "cp_grid_search_geometry",
  "cp_grid_spectrum_rows",
  "cp_grid_candidate_pairs",
  "cp_grid_annotate_candidates",
  "cp_grid_geometry",
  "cp_grid_bind_missing",
  "cp_grid_safe_mean",
  "cp_grid_safe_var",
  "cp_grid_report_scalar",
  "cp_grid_report_metadata",
  "cp_grid_report_append",
  "cp_grid_report_metrics",
  "cp_grid_report_local_minima",
  "cp_grid_report_J",
  "cp_grid_report_fit_spectra",
  "cp_grid_report_numeric",
  "cp_grid_report_summarize",
  "cp_grid_write_reports",
  "cp_grid_checkpoint_file",
  "cp_grid_atomic_rds",
  "cp_grid_config_signature",
  "cp_grid_quality",
  "cp_grid_read_checkpoint",
  "cp_grid_run_job",
  "cp_grid_progress",
  "cp_grid_write_population",
  "cp_grid_make_cluster",
  "cp_run_grid",
  "cp_grid_main",
  "simulation_seeds",
  "simulation_T",
  "simulation_regimes"
)

# Load functions without starting the grid:
# options(cp_grid.no_autorun=TRUE); source("CP_QMLE_grid_20_40_80_100seeds.R")
# Then edit GRID_CONFIG and call cp_run_grid(GRID_CONFIG).
.cp_grid_no_autorun<-isTRUE(getOption("cp_grid.no_autorun",FALSE)) ||
  identical(Sys.getenv("CP_GRID_NO_AUTORUN"),"1")
if(!.cp_grid_no_autorun)cp_grid_main()
