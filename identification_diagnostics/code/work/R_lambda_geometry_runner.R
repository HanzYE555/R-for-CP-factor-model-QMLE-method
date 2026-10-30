# New lambda-axis geometry; saved QMLE files are read-only.
invisible(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE))
Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
source("work/R_diagnostics.R",encoding="UTF-8")
source("work/R_grid_geometry.R",encoding="UTF-8")
source("work/R_local_alignment_core.R",encoding="UTF-8")
source("work/R_lambda_geometry_engine.R",encoding="UTF-8")
stopifnot(requireNamespace("nloptr",quietly=TRUE))
LAMBDA_SOURCE <- "C:/Users/Lenovo/AppData/Local/Temp/7a0a08b7-1844-4bc8-bfe8-652f5efa9228_CP_QMLE_grid_20_40_80_100seeds.zip.228/CP_QMLE_grid_results"
LAMBDA_OUT <- normalizePath("outputs/lambda_axis_geometry_20261030",winslash="/",mustWork=FALSE)
LAMBDA_CFG <- list(diagnostic_random_starts=30L,diagnostic_maxeval=400L,
  diagnostic_cross_rounds=2L,diagnostic_cross_start_limit=5L,
  cap_starts=5L,cap_maxeval=300L,refine_maxeval=400L,cosine_cap=.90)
LAMBDA_SIGNATURE <- paste("lambda-axis-v1;cos=.90;random30;eval400;cross2;caps5x300;refine400",
  unname(tools::md5sum("work/R_lambda_geometry_engine.R")),sep=";")
dir.create(file.path(LAMBDA_OUT,"cases"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(LAMBDA_OUT,"parts"),recursive=TRUE,showWarnings=FALSE)

lg_keyframe <- function(job,n=1L) {
  data.frame(p=rep(job$p,n),q=rep(job$q,n),regime=rep(job$regime,n),
             seed=rep(job$seed,n),T=rep(job$T,n),stringsAsFactors=FALSE)
}
lg_append <- function(x,path) {
  if(is.null(x)||!nrow(x))return(invisible(NULL))
  write.table(x,path,sep=",",quote=TRUE,row.names=FALSE,col.names=!file.exists(path),append=file.exists(path),na="NA")
}
lg_target <- function(job)file.path(LAMBDA_OUT,"cases",sprintf("p%03d_q%03d_%s_s%03d_T%d.rds",job$p,job$q,job$regime,job$seed,job$T))
lg_make_primary <- function(job,result,record) {
  dg<-result$diag; pr<-lg_problem(dg$core,.90)
  cap<-result$cap_rows
  cap<-cap[cap$rho_order==1L,,drop=FALSE];stopifnot(nrow(cap)==1L)
  pairs<-cp_grid_candidate_pairs(dg,pr)
  near<-pairs[pairs$lambda_angle_deg<=5+1e-8,,drop=FALSE]
  pairs$joint_max<-pmax(pairs$rho_at_rho,pairs$E1_at_rho,pairs$rho_at_E1,pairs$E1_at_E1)
  if(nrow(near)){
    near$joint_max<-pmax(near$rho_at_rho,near$E1_at_rho,near$rho_at_E1,near$E1_at_E1)
    best<-near[which.min(near$joint_max),,drop=FALSE]
  }else best<-NULL
  r<-pr$summarize(cp_alignment_x_from_lambda(pr,dg$lambda_rho))
  e<-pr$summarize(cp_alignment_x_from_lambda(pr,dg$lambda_E))
  angles<-vapply(dg$E_local_minima,function(z)cp_alignment_angle_deg(dg$lambda_rho,z$lambda),numeric(1L))
  row<-cbind(lg_keyframe(job),data.frame(d=pr$d,case_index=job$case_index,
    qmle_converged=isTRUE(record$fit$converged),rho_min=dg$rho_min,E1_at_rho=dg$E_at_rho,
    E1_min=dg$E_min,rho_at_E1=dg$rho_at_E,
    rho_max_lambda_cosine=r$max_cosine,E1_max_lambda_cosine=e$max_cosine,
    rho_max_loading_cosine=r$physical_max_cosine,E1_max_loading_cosine=e$physical_max_cosine,
    center_to_E1_nearest_deg_090=min(angles),center_to_independent_E1_le5_090=min(angles)<=5,
    n_rho_local=length(dg$rho_local_minima),n_E_local=length(dg$E_local_minima),
    n_pairs=nrow(pairs),n_pairs_angle_le5=nrow(near),pair_exists=nrow(near)>0L,
    pair_joint_max=if(is.null(best))Inf else best$joint_max,
    pair_joint_angle=if(is.null(best))NA_real_ else best$lambda_angle_deg,
    pair_joint_rho_id=if(is.null(best))NA_integer_ else best$rho_id,
    pair_joint_E1_id=if(is.null(best))NA_integer_ else best$E1_id,
    pair_joint_solvers_success=if(is.null(best))FALSE else best$rho_success&&best$E1_success,
    E1_cap_090=cap$E1_cap5_saved,rho_cap_090=cap$rho_cap5_saved,
    center_to_cap_deg_090=cap$center_to_cap5_deg,
    E1_cap10=cap$E1_cap10,rho_cap10=cap$rho_cap10,
    center_to_cap10_deg=cap$center_to_cap10_deg,
    cap10_angle_edge=cap$cap10_angle_edge,cap10_axis_edge=cap$cap10_axis_edge,
    cap10_max_lambda_cosine=cap$cap10_max_lambda_cosine,
    cap10_max_loading_cosine=cap$cap10_max_loading_cosine,
    local_candidate_exists10=cap$local_candidate_exists10,
    E1_local_within10_found=cap$E1_local_within10_found,
    local_best_angle10=cap$local_best_angle10,
    rho_best_success=isTRUE(dg$rho_best_success),E_best_success=isTRUE(dg$E_best_success),
    cap10_selected_converged=cap$cap10_selected_converged,
    cap10_selected_status=cap$cap10_selected_status,
    domain="lambda_axis_cosine_le_0.90",stringsAsFactors=FALSE))
  for(j in seq_len(8L)){
    row[[sprintf("rho_lambda_%02d",j)]]<-if(j<=pr$d)dg$lambda_rho[j]else NA_real_
    row[[sprintf("E_lambda_%02d",j)]]<-if(j<=pr$d)dg$lambda_E[j]else NA_real_
  }
  list(primary=row,pairs=cbind(lg_keyframe(job,nrow(pairs)),pairs))
}
lg_case <- function(job) {
  target<-lg_target(job)
  if(file.exists(target)){
    x<-readRDS(target)
    if(!identical(x$signature,LAMBDA_SIGNATURE))stop("Incompatible new-domain checkpoint: ",target)
    if(identical(x$status,"ok"))return(x)
  }
  path<-file.path(LAMBDA_SOURCE,"checkpoints",sprintf("p%03d_q%03d",job$p,job$q),
                 job$regime,sprintf("seed_%03d",job$seed),paste0("T",job$T,".rds"))
  started<-proc.time()[[3L]]
  x<-tryCatch({
    old<-readRDS(path);stopifnot(identical(old$status,"ok"))
    result<-lg_search(old,as.integer(20261101+job$case_index*101L),LAMBDA_CFG)
    if(!is.null(result$independent)&&is.null(result$diag))result$diag<-result$independent$diag
    if(!is.null(result$caps)&&is.null(result$cap_rows))result$cap_rows<-result$caps$rows
    summary<-lg_make_primary(job,result,old)
    result$primary<-summary$primary;result$pairs<-summary$pairs
    result$cap_rows<-cbind(lg_keyframe(job,nrow(result$cap_rows)),result$cap_rows)
    result$problem<-NULL
    c(list(status="ok",error=""),result)
  },error=function(e)list(status="error",error=conditionMessage(e)))
  x$job<-job;x$signature<-LAMBDA_SIGNATURE;x$source_path<-path
  x$source_md5<-unname(tools::md5sum(path));x$source_bytes<-file.info(path)$size
  x$seconds<-proc.time()[[3L]]-started
  tmp<-paste0(target,".tmp");saveRDS(x,tmp,compress=TRUE)
  if(file.exists(target)){
    previous<-paste0(target,".previous");if(!file.exists(previous))file.rename(target,previous)
  }
  stopifnot(file.rename(tmp,target));x
}
lg_main <- function() {
  base<-read.csv("work/grid_final_geometry_by_seed.csv",stringsAsFactors=FALSE)
  jobs<-base[,c("p","q","regime","seed","T")]
  jobs<-jobs[order(jobs$p,jobs$q,match(jobs$regime,c("interior","boundary","exterior")),jobs$seed,jobs$T),]
  jobs$case_index<-seq_len(nrow(jobs));rownames(jobs)<-NULL
  stopifnot(nrow(jobs)==10796L)
  args<-commandArgs(trailingOnly=TRUE)
  pilot<-"--pilot"%in%args
  if(pilot)jobs<-jobs[jobs$p==40&jobs$q==40&jobs$seed==1&jobs$T%in%c(100,800),]
  at<-match("--shard",args);shard<-if(is.na(at))1L else as.integer(args[at+1L])
  ns<-if(is.na(at))1L else as.integer(args[at+2L])
  if(shard==1L)write.csv(jobs,file.path(LAMBDA_OUT,if(pilot)"pilot_manifest.csv"else"manifest.csv"),row.names=FALSE)
  jobs<-jobs[(jobs$case_index-1L)%%ns==shard-1L,,drop=FALSE]
  prefix<-if(pilot)"pilot"else sprintf("shard%02d",shard)
  # Part files are regenerated on restart from validated checkpoints.
  output<-list(primary=NULL,caps=NULL,pairs=NULL)
  initialized<-c(primary=FALSE,caps=FALSE,pairs=FALSE)
  schemas<-list(primary=NULL,caps=NULL,pairs=NULL)
  statuses<-list()
  for(i in seq_len(nrow(jobs))){
    result<-lg_case(as.list(jobs[i,]))
    statuses[[i]]<-data.frame(case_index=jobs$case_index[i],status=result$status,
      n_rho=if(result$status=="ok")nrow(result$cap_rows)else 0L,
      seconds=result$seconds,error=result$error,stringsAsFactors=FALSE)
    if(result$status=="ok"){
      output$primary[[length(output$primary)+1L]]<-result$primary
      output$caps[[length(output$caps)+1L]]<-result$cap_rows
      output$pairs[[length(output$pairs)+1L]]<-result$pairs
    }
    if(i%%8L==0L||i==nrow(jobs)){
      write.csv(do.call(rbind,statuses),file.path(LAMBDA_OUT,paste0(prefix,"_progress.csv")),row.names=FALSE)
      for(what in names(output))if(length(output[[what]])){
        tab<-cp_grid_bind_frames(output[[what]])
        if(is.null(schemas[[what]]))schemas[[what]]<-names(tab)
        stopifnot(identical(names(tab),schemas[[what]]))
        write.table(tab,file.path(LAMBDA_OUT,"parts",paste0(prefix,"_",what,".csv")),
          sep=",",quote=TRUE,row.names=FALSE,col.names=!initialized[what],append=initialized[what],na="NA")
        initialized[what]<-TRUE;output[[what]]<-NULL
      }
      cat(sprintf("LAMBDA %s: %d/%d cases; errors=%d; last_seconds=%.2f\n",prefix,i,nrow(jobs),
        sum(vapply(statuses,function(z)z$status!="ok",logical(1L))),result$seconds));flush.console()
    }
  }
  if(any(vapply(statuses,function(z)z$status!="ok",logical(1L))))stop("Inspect new-domain errors.")
}
if(!isTRUE(getOption("lambda.runner.no_autorun",FALSE)))lg_main()
