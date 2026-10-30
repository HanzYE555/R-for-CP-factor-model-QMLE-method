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
