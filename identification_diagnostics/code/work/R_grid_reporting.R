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
