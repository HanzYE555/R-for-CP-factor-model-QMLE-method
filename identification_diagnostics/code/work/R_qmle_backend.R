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
