# Experiment driver and CSV/RDS reporting. Modules are combined into one file.

cp_config_signature <- function(cfg) {
  # Seeds and T are excluded so a later extension can reuse existing fits.
  fields <- c("p", "q", "maxit", "qmle_starts", "diagnostic_starts",
              "diagnostic_maxeval", "cosine_cap", "estimate_noise",
              "noise_bounds", "factor_multiplier", "burn", "n_exterior")
  paste(vapply(fields, function(nm) paste(nm, paste(cfg[[nm]], collapse=","), sep="="),
               character(1L)), collapse=";")
}

cp_quality_metrics <- function(fit, design) {
  Z <- cp_khatri_rao(fit$A, fit$B); Z0 <- design$Z
  orth <- function(X) {
    z <- qr(X, tol=1e-10); qr.Q(z)[, seq_len(z$rank), drop=FALSE]
  }
  Q <- orth(Z); Q0 <- orth(Z0); d <- ncol(Z0)
  subspace_error <- sqrt(max(0, (ncol(Q)+ncol(Q0)-2*sum(crossprod(Q,Q0)^2))/(2*d)))
  F <- Z %*% t(chol(fit$M)); F0 <- Z0 %*% t(chol(design$Mx))
  diag_signal <- rowSums(F^2); diag_signal0 <- rowSums(F0^2)
  delta_noise <- fit$psi - 1
  numerator2 <- sum(crossprod(F)^2)+sum(crossprod(F0)^2)-2*sum(crossprod(F,F0)^2)+
    sum(delta_noise^2)+2*sum(delta_noise*(diag_signal-diag_signal0))
  denominator2 <- sum(crossprod(F0)^2)+nrow(F0)+2*sum(diag_signal0)
  list(subspace_error=subspace_error,
       covariance_relative_error=sqrt(max(0,numerator2)/denominator2),
       estimated_Z_rank=ncol(Q))
}

cp_checkpoint_file <- function(cfg, regime, seed, T) {
  file.path(cfg$output_dir, "checkpoints", regime, sprintf("seed_%03d",seed), paste0("T",T,".rds"))
}

cp_run_seed <- function(seed, cfg) {
  completed <- 0L; errors <- 0L
  for (regime in cfg$regimes) {
    design <- generate_data(regime, seed, Tmax=max(cfg$sample_sizes), p=cfg$p, q=cfg$q,
                            burn=cfg$burn, factor_multiplier=cfg$factor_multiplier,
                            n=cfg$n_exterior)
    seed_dir <- dirname(cp_checkpoint_file(cfg,regime,seed,cfg$sample_sizes[1L]))
    dir.create(seed_dir, recursive=TRUE, showWarnings=FALSE)
    saveRDS(design$metadata, file.path(seed_dir,"design_metadata.rds"))
    if (cfg$save_data) saveRDS(design, file.path(seed_dir,"generated_data.rds"))
    previous_fit <- NULL
    for (T in sort(cfg$sample_sizes)) {
      path <- cp_checkpoint_file(cfg,regime,seed,T)
      if (cfg$resume && file.exists(path)) {
        old <- readRDS(path)
        if (!identical(old$config_signature,cfg$signature)) {
          stop("Checkpoint configuration differs. Choose a new output_dir or set resume=FALSE: ",path)
        }
        if (identical(old$status,"ok")) {
          previous_fit <- old$fit; completed <- completed+1L; next
        }
      }
      if (cfg$verbose) message(sprintf("[%s] seed=%d, T=%d",regime,seed,T))
      tick <- proc.time()[[3L]]
      Y <- centered_prefix(design,T)$Y
      result <- list(regime=regime,seed=seed,T=T,p=cfg$p,q=cfg$q,d=design$d,
                     config_signature=cfg$signature, status="fit_error",
                     error_message="",truth_metadata=design$metadata)
      fit <- tryCatch(cp_qmle(Y,cfg$p,cfg$q,design$d,
                             seed=seed+17L*T+100000L*match(regime,simulation_regimes),
                             n_starts=cfg$qmle_starts,maxit=cfg$maxit,
                             estimate_noise=cfg$estimate_noise,fixed_noise=1,
                             previous_fit=previous_fit,noise_bounds=cfg$noise_bounds),
                      error=function(e)e)
      if (inherits(fit,"error")) {
        result$error_message <- conditionMessage(fit); errors <- errors+1L
      } else {
        result$fit <- fit; previous_fit <- fit
        result$quality <- cp_quality_metrics(fit,design)
        result$status <- "diagnostic_error"
        diag <- tryCatch(cp_diagnostics(fit$A,fit$B,design$A,design$B,
                         seed=seed+31L*T+100000L*match(regime,simulation_regimes),
                         nstarts=cfg$diagnostic_starts,cosine_cap=cfg$cosine_cap,
                         maxeval=cfg$diagnostic_maxeval,core_ranks=NULL),error=function(e)e)
        if (inherits(diag,"error")) {
          result$error_message <- conditionMessage(diag); errors <- errors+1L
          result$core_fallback <- tryCatch(cp_prepare_core(fit$A,fit$B,design$A,design$B),
                                           error=function(e)NULL)
        } else {
          result$diagnostics <- diag; result$status <- "ok"; completed <- completed+1L
          result$Ar <- diag$core$A; result$Br <- diag$core$B
          # Ar/Br are the ONLY coordinates used in all reported geometry.
          stopifnot(max(abs(colSums(diag$core$A^2)-1))<1e-8,
                    max(abs(colSums(diag$core$B^2)-1))<1e-8)
        }
      }
      result$seconds <- proc.time()[[3L]]-tick
      result$completed_at <- format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")
      saveRDS(result,path,compress=TRUE)
      if (!identical(result$status,"ok")) message("  ",result$status,": ",result$error_message)
      rm(Y,result); gc(verbose=FALSE)
    }
    rm(design); gc(verbose=FALSE)
  }
  data.frame(seed=seed,completed=completed,errors=errors)
}

cp_bind_rows <- function(rows) if (length(rows)) do.call(rbind,rows) else data.frame()
cp_safe_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
cp_safe_var <- function(x) if (sum(is.finite(x))>1L) stats::var(x[is.finite(x)]) else NA_real_

cp_reported_J_svd <- function(diagnostics) {
  if (!is.null(diagnostics$J_orig)) {
    orig <- diagnostics$J_orig
  } else if (!is.null(diagnostics$J_eigenvalues)) {
    # Legacy checkpoint: convert each seed before taking means/variances.
    orig <- sqrt(pmax(diagnostics$J_eigenvalues,0))
  } else {
    stop("Checkpoint contains neither J_orig nor legacy J_eigenvalues.")
  }
  maxima <- apply(orig,1L,max)
  nor <- sweep(orig,1L,maxima,"/")
  nor[maxima<=0,] <- NA_real_
  list(orig=orig,nor=nor)
}

cp_collect_results <- function(cfg) {
  metrics <- list(); spectra <- list(); minima <- list(); starts <- list()
  for (regime in cfg$regimes) for (T in cfg$sample_sizes) for (seed in cfg$seeds) {
    path <- cp_checkpoint_file(cfg,regime,seed,T)
    if (!file.exists(path)) next
    r <- readRDS(path)
    if (!identical(r$config_signature,cfg$signature)) next
    fit <- r$fit; dg <- r$diagnostics
    scalar <- function(x) if (length(x)==1L) x else NA_real_
    row <- data.frame(regime=regime,seed=seed,T=T,p=r$p,q=r$q,d=r$d,
       status=r$status,noise_estimated=cfg$estimate_noise,
       qmle_converged=if (is.null(fit)) FALSE else fit$converged,
       nll=scalar(fit$nll),projected_gradient=scalar(fit$projected_gradient_norm),
       maxit_per_start=cfg$maxit,qmle_starts=cfg$qmle_starts,
       rho_min=scalar(dg$rho_min),E1_at_rho=scalar(dg$E_at_rho),
       E1_min=scalar(dg$E_min),rho_at_E1=scalar(dg$rho_at_E),
       E1_at_rho_percent=100*scalar(dg$E_at_rho),E1_min_percent=100*scalar(dg$E_min),
       J_min=scalar(dg$J_min),J_relative_min=scalar(dg$J_relative_min),
       rho_best_success=if (is.null(dg)) FALSE else dg$rho_best_success,
       E1_best_success=if (is.null(dg)) FALSE else dg$E_best_success,
       rho_local_count=scalar(dg$n_rho_local_minima),E1_local_count=scalar(dg$n_E_local_minima),
       subspace_error=scalar(r$quality$subspace_error),
       covariance_relative_error=scalar(r$quality$covariance_relative_error),
       core_p=scalar(dg$core$core_p),core_q=scalar(dg$core$core_q),
       fitted_A_numerical_rank=scalar(dg$core$numerical_rank_A),
       fitted_B_numerical_rank=scalar(dg$core$numerical_rank_B),
       matching_cosine_mean=if (is.null(dg)) NA_real_ else mean(dg$core$matching$matched_rankone_cosines),
       noise_at_lower_bound=if (is.null(fit)) NA_integer_ else sum(fit$psi<=cfg$noise_bounds[1L]*(1+1e-4)),
       noise_at_upper_bound=if (is.null(fit)) NA_integer_ else sum(fit$psi>=cfg$noise_bounds[2L]*(1-1e-4)),
       seconds=scalar(r$seconds),error_message=r$error_message,stringsAsFactors=FALSE)
    metrics[[length(metrics)+1L]] <- row
    if (!is.null(fit)) {
      st <- fit$all_starts
      st$regime <- regime; st$seed <- seed; st$T <- T
      starts[[length(starts)+1L]] <- st
    }
    if (is.null(dg)) next
    reported_svd <- cp_reported_J_svd(dg)
    for (ell in seq_len(r$d)) for (j in seq_len(r$d-1L)) {
      spectra[[length(spectra)+1L]] <- data.frame(regime=regime,seed=seed,T=T,
        ell=ell,svd_index=j,orig=reported_svd$orig[ell,j],
        nor=reported_svd$nor[ell,j],qmle_converged=fit$converged)
    }
    for (objective in c("rho","E1")) {
      local <- if (objective=="rho") dg$rho_local_minima else dg$E_local_minima
      for (point in local) {
        mr <- data.frame(regime=regime,seed=seed,T=T,objective=objective,
          local_minimum_order=point$minimum_id,objective_value=point$objective_value,
          rho=point$rho,E1=point$E1,rankone_error_percent=100*point$E1,
          max_loading_cosine=point$max_cosine,matrix_norm=point$matrix_norm,
          hit_count=point$hit_count,converged_hits=point$converged_hit_count,
          solver_converged=point$success,solver_status=point$status,
          solver_message=point$message,qmle_converged=fit$converged,stringsAsFactors=FALSE)
        unit <- point$lambda/sqrt(sum(point$lambda^2))
        for (j in seq_len(8L)) {
          mr[[sprintf("lambda_%02d",j)]] <- if(j<=r$d)point$lambda[j] else NA_real_
          mr[[sprintf("lambda_unit_%02d",j)]] <- if(j<=r$d)unit[j] else NA_real_
        }
        minima[[length(minima)+1L]] <- mr
      }
    }
  }
  list(metrics=cp_bind_rows(metrics),spectra=cp_bind_rows(spectra),
       minima=cp_bind_rows(minima),starts=cp_bind_rows(starts))
}

cp_summarize_metrics <- function(df,cfg) {
  rows <- list()
  fields <- c("rho_min","E1_at_rho","E1_min","rho_at_E1",
              "E1_at_rho_percent","E1_min_percent","J_min","J_relative_min",
              "subspace_error","covariance_relative_error","matching_cosine_mean",
              "rho_local_count","E1_local_count","projected_gradient","seconds")
  for (regime in cfg$regimes) for (T in cfg$sample_sizes) {
    group <- df[df$regime==regime & df$T==T,,drop=FALSE]
    for (subset in c("all_valid","qmle_converged")) {
      x <- if(subset=="qmle_converged")group[group$qmle_converged,,drop=FALSE] else group
      for (field in fields) {
        value <- x[[field]]
        rows[[length(rows)+1L]] <- data.frame(regime=regime,T=T,subset=subset,metric=field,
          mean=cp_safe_mean(value),variance=cp_safe_var(value),sd=sqrt(cp_safe_var(value)),
          n_finite=sum(is.finite(value)),n_saved=nrow(group),n_requested=length(cfg$seeds),
          n_qmle_converged=sum(group$qmle_converged),n_diagnostic_ok=sum(group$status=="ok"))
      }
    }
  }
  cp_bind_rows(rows)
}

cp_summarize_spectra <- function(df,cfg) {
  if (!nrow(df)) return(data.frame())
  groups <- split(seq_len(nrow(df)),interaction(df$regime,df$T,df$ell,df$svd_index,drop=TRUE))
  rows <- lapply(groups,function(ii) {
    x <- df[ii,,drop=FALSE]
    data.frame(regime=x$regime[1L],T=x$T[1L],ell=x$ell[1L],svd_index=x$svd_index[1L],
      orig_mean=cp_safe_mean(x$orig),orig_variance=cp_safe_var(x$orig),
      nor_mean=cp_safe_mean(x$nor),nor_variance=cp_safe_var(x$nor),
      n_seeds=nrow(x),n_requested=length(cfg$seeds),n_qmle_converged=sum(x$qmle_converged))
  })
  ans <- cp_bind_rows(rows)
  ans[order(match(ans$regime,cfg$regimes),ans$T,ans$ell,ans$svd_index),,drop=FALSE]
}

cp_summarize_minima <- function(df,cfg) {
  if (!nrow(df)) return(data.frame())
  groups <- split(seq_len(nrow(df)),interaction(df$regime,df$T,df$objective,df$local_minimum_order,drop=TRUE))
  rows <- lapply(groups,function(ii) {
    x <- df[ii,,drop=FALSE]
    data.frame(regime=x$regime[1L],T=x$T[1L],objective=x$objective[1L],
      local_minimum_order=x$local_minimum_order[1L],
      objective_mean=cp_safe_mean(x$objective_value),objective_variance=cp_safe_var(x$objective_value),
      rho_mean=cp_safe_mean(x$rho),rho_variance=cp_safe_var(x$rho),
      E1_mean=cp_safe_mean(x$E1),E1_variance=cp_safe_var(x$E1),
      rankone_error_percent_mean=cp_safe_mean(x$rankone_error_percent),
      rankone_error_percent_variance=cp_safe_var(x$rankone_error_percent),
      n_seeds_present=nrow(x),n_requested=length(cfg$seeds),
      n_solver_converged=sum(x$solver_converged),n_qmle_converged=sum(x$qmle_converged),
      aggregation_note="k-th smallest numerical local candidate per seed; not a matched geometric branch")
  })
  ans <- cp_bind_rows(rows)
  ans[order(match(ans$regime,cfg$regimes),ans$T,ans$objective,ans$local_minimum_order),,drop=FALSE]
}

cp_write_reports <- function(cfg) {
  raw <- cp_collect_results(cfg)
  if (!nrow(raw$metrics)) stop("No checkpoint was available for reporting.")
  summarize_subsets <- function(frame,fun) {
    all <- fun(frame,cfg)
    conv <- fun(frame[frame$qmle_converged,,drop=FALSE],cfg)
    all$subset <- rep("all_valid",nrow(all))
    conv$subset <- rep("qmle_converged",nrow(conv))
    if (!nrow(all)) return(conv)
    if (!nrow(conv)) return(all)
    rbind(all,conv)
  }
  summaries <- list(metrics=cp_summarize_metrics(raw$metrics,cfg),
                    spectra=summarize_subsets(raw$spectra,cp_summarize_spectra),
                    minima=summarize_subsets(raw$minima,cp_summarize_minima))
  write <- function(x,name) utils::write.csv(x,file.path(cfg$output_dir,name),row.names=FALSE,na="")
  write(raw$metrics,"metrics_by_seed.csv")
  write(summaries$metrics,"metrics_mean_variance.csv")
  write(raw$spectra,"J_svd_by_seed.csv")
  write(summaries$spectra,"J_svd_mean_variance.csv")
  write(raw$minima,"local_minima_by_seed.csv")
  write(summaries$minima,"local_minima_mean_variance.csv")
  write(raw$starts,"qmle_all_starts.csv")
  write(raw$metrics[raw$metrics$status!="ok" | !raw$metrics$qmle_converged,,drop=FALSE],
        "nonconvergence_and_errors.csv")
  saveRDS(list(config=cfg,raw=raw,summaries=summaries),file.path(cfg$output_dir,"all_results.rds"))
  invisible(list(raw=raw,summaries=summaries))
}

cp_population_benchmarks <- function(cfg) {
  file <- file.path(cfg$output_dir,"population_benchmarks.rds")
  result <- list()
  for (regime in cfg$regimes) {
    truth <- core_loading_design(regime,cfg$n_exterior)
    result[[regime]] <- cp_diagnostics(truth$Ar_unit,truth$Br_unit,seed=1819L,
       nstarts=cfg$diagnostic_starts,cosine_cap=cfg$cosine_cap,
       maxeval=cfg$diagnostic_maxeval,extra_lambdas=truth$extra_lambda)
  }
  saveRDS(result,file)
  table <- do.call(rbind,lapply(names(result),function(regime) {
    r <- result[[regime]]
    data.frame(regime=regime,rho_min=r$rho_min,E1_at_rho=r$E_at_rho,E1_min=r$E_min,
               rho_at_E1=r$rho_at_E,J_min=r$J_min,rho_local_count=r$n_rho_local_minima,
               E1_local_count=r$n_E_local_minima)
  }))
  utils::write.csv(table,file.path(cfg$output_dir,"population_benchmarks.csv"),row.names=FALSE)
  invisible(result)
}

cp_run_experiment <- function(cfg) {
  cp_ensure_dependencies(cfg)
  stopifnot(cfg$p>=8L,cfg$q>=8L,cfg$maxit>=1L,length(cfg$seeds)>0L,
            all(cfg$regimes %in% simulation_regimes),all(cfg$sample_sizes>=2L))
  cfg$output_dir <- normalizePath(cfg$output_dir,winslash="/",mustWork=FALSE)
  cfg$signature <- cp_config_signature(cfg)
  dir.create(cfg$output_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(cfg,file.path(cfg$output_dir,"configuration.rds"))
  writeLines(capture.output(str(cfg)),file.path(cfg$output_dir,"configuration.txt"))
  writeLines(capture.output(sessionInfo()),file.path(cfg$output_dir,"sessionInfo.txt"))
  cp_population_benchmarks(cfg)
  if (cfg$ncores>1L) {
    cl <- parallel::makePSOCKcluster(cfg$ncores,outfile="")
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(paths) {
      .libPaths(paths)
      if (!requireNamespace("nloptr",quietly=TRUE)) stop("Worker cannot load package 'nloptr'.")
      NULL
    },.libPaths())
    parallel::clusterExport(cl,CP_EXPORTS,envir=environment(cp_run_experiment))
    jobs <- parallel::parLapplyLB(cl,cfg$seeds,function(seed,cfg) cp_run_seed(seed,cfg),cfg=cfg)
  } else {
    jobs <- lapply(cfg$seeds,function(seed) cp_run_seed(seed,cfg))
  }
  utils::write.csv(cp_bind_rows(jobs),file.path(cfg$output_dir,"job_status.csv"),row.names=FALSE)
  result <- cp_write_reports(cfg)
  message("Completed. Results: ",cfg$output_dir)
  message("Mean/variance tables include sample variance var(x), denominator n-1.")
  invisible(result)
}
