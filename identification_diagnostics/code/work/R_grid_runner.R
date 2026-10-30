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
