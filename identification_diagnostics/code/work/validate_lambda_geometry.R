# Independent direct-matrix validation of the lambda-axis geometry.
# This file never writes source fits, QMLE checkpoints, or archived results.
if(.Platform$OS.type=="windows")invisible(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE))
Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")

lv_angle <- function(u,v) {
  a<-abs(sum(u*v))/sqrt(sum(u*u)*sum(v*v))
  acos(min(1,max(0,a)))*180/pi
}

lv_direct <- function(A,B,lambda) {
  A<-as.matrix(A);B<-as.matrix(B);lambda<-as.numeric(lambda)
  stopifnot(ncol(A)==ncol(B),length(lambda)==ncol(A),all(is.finite(lambda)))
  M<-A%*%diag(lambda,ncol(A))%*%t(B)
  norm2<-sum(M*M);stopifnot(is.finite(norm2),norm2>0)
  ss<-svd(M,nu=0L,nv=0L)$d
  minor2<-0
  for(rr in combn(nrow(M),2L,simplify=FALSE))
    for(cc in combn(ncol(M),2L,simplify=FALSE)) {
      mm<-M[rr[1L],cc[1L]]*M[rr[2L],cc[2L]]-M[rr[1L],cc[2L]]*M[rr[2L],cc[1L]]
      minor2<-minor2+mm*mm
    }
  kk<-combn(lambda,2L,FUN=prod);den<-sum(kk*kk)
  cnorm<-sqrt(colSums(A*A)*colSums(B*B))
  loading_cosine<-vapply(seq_len(ncol(A)),function(j)
    abs(sum(M*outer(A[,j],B[,j])))/sqrt(norm2)/cnorm[j],numeric(1L))
  axis_cosine<-max(abs(lambda))/sqrt(sum(lambda*lambda))
  list(E1=sqrt(sum(ss[-1L]^2)/norm2),rho=if(den>0)sqrt(minor2/den)else NA_real_,
       H=sqrt(minor2),matrix_norm=sqrt(norm2),lambda_norm=sqrt(sum(lambda*lambda)),
       axis_cosine=axis_cosine,axis_angle_deg=acos(min(1,axis_cosine))*180/pi,
       max_loading_cosine=max(loading_cosine),loading_cosines=loading_cosine,
       top_squared_gap=if(length(ss)>1L)ss[1L]^2-ss[2L]^2 else ss[1L]^2)
}

lv_population <- function(out="work/lambda_geometry_validation") {
  source("work/R_design.R",encoding="UTF-8")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  rows<-list()
  for(regime in c("interior","boundary","exterior")) {
    core<-core_loading_design(regime,n=.5)
    A<-core$Ar_unit;B<-core$Br_unit
    Z<-vapply(seq_len(ncol(A)),function(j)as.vector(outer(A[,j],B[,j])),numeric(nrow(A)*nrow(B)))
    stopifnot(qr(Z)$rank==ncol(A))
    # A coordinate axis is excluded by the coefficient-domain rule even
    # though it is an exact rank-one direction.
    ax<-lv_direct(A,B,c(1,rep(0,ncol(A)-1L)))
    stopifnot(ax$axis_cosine>.9,ax$E1<1e-12)
    rows[[length(rows)+1L]]<-data.frame(regime=regime,kind="given_axis",E1=ax$E1,
      rho=ax$rho,axis_cosine=ax$axis_cosine,axis_angle_deg=ax$axis_angle_deg,
      max_loading_cosine=ax$max_loading_cosine,lambda_domain_feasible=FALSE)
    if(regime=="exterior") {
      lam<-core$extra_lambda
      p<-lv_direct(A,B,lam)
      stopifnot(p$E1<1e-12,p$rho<1e-12,p$axis_cosine<.9,
                abs(p$axis_angle_deg-49.0558760194184)<1e-10,
                abs(p$max_loading_cosine-.5270462766947305)<1e-12)
      rows[[length(rows)+1L]]<-data.frame(regime=regime,kind="extra_rank_one",E1=p$E1,
        rho=p$rho,axis_cosine=p$axis_cosine,axis_angle_deg=p$axis_angle_deg,
        max_loading_cosine=p$max_loading_cosine,lambda_domain_feasible=TRUE)
    }
  }
  # A physical cosine exceeding the old threshold is allowed in the new
  # domain: orthogonal coefficient axes and matrix angles are distinct.
  A<-rbind(c(1,1),c(0,0));B<-rbind(c(1,.99),c(0,sqrt(1-.99^2)))
  p<-lv_direct(A,B,c(1,1))
  stopifnot(p$axis_cosine<.9,p$max_loading_cosine>.9,p$E1<1e-12)
  rows[[length(rows)+1L]]<-data.frame(regime="constructed_nonidentifiable",kind="physical_cosine_over_old_cap",
    E1=p$E1,rho=p$rho,axis_cosine=p$axis_cosine,axis_angle_deg=p$axis_angle_deg,
    max_loading_cosine=p$max_loading_cosine,lambda_domain_feasible=TRUE)
  tab<-do.call(rbind,rows)
  write.csv(tab,file.path(out,"population_direct_validation.csv"),row.names=FALSE)
  jsonlite::write_json(list(status="passed",population_checks=nrow(tab),
    new_domain="max(abs(lambda))/norm(lambda)<=0.90",axis_exclusion_deg=acos(.9)*180/pi,
    old_physical_cap_not_applied=TRUE,QMLE_refitted=FALSE),
    file.path(out,"population_direct_validation.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
  tab
}

lv_fd <- function(f,x,h=1e-6) {
  y<-as.numeric(f(x));ans<-matrix(NA_real_,length(y),length(x))
  for(j in seq_along(x)) {
    delta<-rep(0,length(x));delta[j]<-h*max(1,abs(x[j]))
    ans[,j]<-(as.numeric(f(x+delta))-as.numeric(f(x-delta)))/(2*delta[j])
  }
  if(length(y)==1L)drop(ans)else ans
}

# Objective adapter arguments are explicit so that this audit does not
# depend on the engine's naming conventions or reuse its value formulas.
lv_gradient_audit <- function(problem,rho_value_gradient,E_value_gradient,
                              axis_value,axis_jacobian,n=20L,seed=4901L,
                              L=NULL,A=NULL,B=NULL,cap=.90) {
  if(is.null(L))L<-if(!is.null(problem$L))problem$L else problem$Rinv
  if(is.null(A))A<-problem$A
  if(is.null(B))B<-problem$B
  stopifnot(is.matrix(L),nrow(L)==ncol(A))
  d<-ncol(A);set.seed(seed);rows<-list();accepted<-0L;attempts<-0L
  while(accepted<n && attempts<10000L) {
    attempts<-attempts+1L
    lam<-rnorm(d);lam<-lam/sqrt(sum(lam*lam))
    if(max(abs(lam))>=cap-.04)next
    x<-drop(solve(L,lam));x<-x/sqrt(sum(x*x))
    p<-lv_direct(A,B,drop(L%*%x))
    if(p$top_squared_gap<1e-4)next
    accepted<-accepted+1L
    direct_rho<-function(xx)lv_direct(A,B,drop(L%*%xx))$rho^2
    direct_E<-function(xx)lv_direct(A,B,drop(L%*%xx))$E1^2
    direct_g<-function(xx){ll<-drop(L%*%xx);ll^2-cap^2*sum(ll*ll)}
    fd_rho<-lv_fd(direct_rho,x);fd_E<-lv_fd(direct_E,x);fd_g<-lv_fd(direct_g,x)
    rho<-rho_value_gradient(x);ee<-E_value_gradient(x)
    ja<-axis_jacobian(x)
    rows[[accepted]]<-data.frame(point=accepted,
      rho_value_error=abs(rho$value-direct_rho(x)),
      E1_squared_value_error=abs(ee$value-direct_E(x)),
      axis_constraint_value_error=max(abs(axis_value(x)-direct_g(x))),
      rho_gradient_error=max(abs(rho$gradient-fd_rho)),
      E1_squared_gradient_error=max(abs(ee$gradient-fd_E)),
      axis_constraint_gradient_error=max(abs(ja-fd_g)),
      rho_gradient_scaled_error=max(abs(rho$gradient-fd_rho))/(1+max(abs(fd_rho))),
      E1_gradient_scaled_error=max(abs(ee$gradient-fd_E))/(1+max(abs(fd_E))),
      axis_gradient_scaled_error=max(abs(ja-fd_g))/(1+max(abs(fd_g))),
      axis_cosine=p$axis_cosine,physical_max_cosine=p$max_loading_cosine)
  }
  stopifnot(accepted==n)
  tab<-do.call(rbind,rows)
  stopifnot(max(tab$rho_value_error)<1e-9,max(tab$E1_squared_value_error)<1e-10,
            max(tab$axis_constraint_value_error)<1e-10,
            max(tab$rho_gradient_scaled_error)<2e-6,
            max(tab$E1_gradient_scaled_error)<2e-6,
            max(tab$axis_gradient_scaled_error)<2e-6)
  tab
}

lv_point_check <- function(A,B,point,center=NULL,angle_limit=NULL,
                          cap=.9,numeric_tol=1e-7,angle_tol=1e-4) {
  stopifnot(is.list(point),!is.null(point$lambda))
  p<-lv_direct(A,B,point$lambda)
  result<-data.frame(E1_error=if(is.null(point$E1))NA_real_ else abs(point$E1-p$E1),
    rho_error=if(is.null(point$rho))NA_real_ else abs(point$rho-p$rho),
    reported_axis_error=if(is.null(point$axis_max_cosine))NA_real_ else abs(point$axis_max_cosine-p$axis_cosine),
    reported_physical_error=if(is.null(point$physical_max_cosine))NA_real_ else abs(point$physical_max_cosine-p$max_loading_cosine),
    reported_axis_angle_error=if(is.null(point$axis_min_angle_deg))NA_real_ else abs(point$axis_min_angle_deg-p$axis_angle_deg),
    matrix_norm_error=abs(p$matrix_norm-1),axis_violation=max(0,p$axis_cosine-cap),
    axis_cosine=p$axis_cosine,physical_max_cosine=p$max_loading_cosine,
    lambda_angle_deg=if(is.null(center))NA_real_ else lv_angle(point$lambda,center$lambda))
  if(!is.null(angle_limit))result$angle_violation<-max(0,result$lambda_angle_deg-angle_limit)
  stopifnot(all(result$E1_error[is.finite(result$E1_error)]<numeric_tol),
            all(result$rho_error[is.finite(result$rho_error)]<numeric_tol),
            all(result$reported_axis_error[is.finite(result$reported_axis_error)]<numeric_tol),
            all(result$reported_physical_error[is.finite(result$reported_physical_error)]<numeric_tol),
            all(result$reported_axis_angle_error[is.finite(result$reported_axis_angle_error)]<1e-5),
            result$matrix_norm_error<2e-6,result$axis_violation<2e-6)
  if(!is.null(angle_limit))stopifnot(result$angle_violation<angle_tol)
  # Intentionally no check that physical_max_cosine <= .90.
  result
}

lv_pair_pool_check <- function(A,B,rho_pool,E_pool,alpha=5,cap=.9) {
  if(!length(rho_pool)||!length(E_pool))return(list(Q=Inf,min_angle=Inf,n_eligible=0L))
  scored<-list();all_angles<-numeric()
  for(i in seq_along(rho_pool))for(j in seq_along(E_pool)) {
    r<-rho_pool[[i]];e<-E_pool[[j]]
    pr<-lv_direct(A,B,r$lambda);pe<-lv_direct(A,B,e$lambda)
    stopifnot(pr$axis_cosine<=cap+2e-6,pe$axis_cosine<=cap+2e-6)
    aa<-lv_angle(r$lambda,e$lambda);all_angles<-c(all_angles,aa)
    if(aa<=alpha+1e-8)scored[[length(scored)+1L]]<-c(i=i,j=j,angle=aa,
      Q=max(pr$rho,pr$E1,pe$rho,pe$E1))
  }
  if(!length(scored))return(list(Q=Inf,min_angle=min(all_angles),n_eligible=0L))
  mat<-do.call(rbind,scored);best<-mat[which.min(mat[,"Q"]),]
  list(Q=unname(best["Q"]),min_angle=min(all_angles),n_eligible=nrow(mat),best=best)
}

lv_engine_audit <- function(out="work/lambda_geometry_validation",n=20L) {
  source("work/R_diagnostics.R",encoding="UTF-8")
  source("work/R_grid_geometry.R",encoding="UTF-8")
  source("work/R_local_alignment_core.R",encoding="UTF-8")
  source("work/R_lambda_geometry_engine.R",encoding="UTF-8")
  source("work/R_design.R",encoding="UTF-8")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  rows<-list()
  for(regime in c("interior","boundary","exterior")) {
    truth<-core_loading_design(regime,n=.5)
    pr<-lg_problem(list(A=truth$Ar_unit,B=truth$Br_unit),.9)
    tab<-lv_gradient_audit(pr,pr$rho_squared,pr$E_squared,pr$g_ineq,pr$jac_ineq,
                           n=n,seed=4901L+match(regime,c("interior","boundary","exterior")))
    tab$regime<-regime;rows[[length(rows)+1L]]<-tab
  }
  ans<-do.call(rbind,rows)
  write.csv(ans,file.path(out,"gradient_direct_validation.csv"),row.names=FALSE)
  jsonlite::write_json(list(status="passed",points=nrow(ans),
    maximum_rho_value_error=max(ans$rho_value_error),
    maximum_E1_squared_value_error=max(ans$E1_squared_value_error),
    maximum_axis_constraint_value_error=max(ans$axis_constraint_value_error),
    maximum_rho_gradient_error=max(ans$rho_gradient_error),
    maximum_E1_gradient_error=max(ans$E1_squared_gradient_error),
    maximum_axis_gradient_error=max(ans$axis_constraint_gradient_error),
    values_reconstructed_from="literal A diag(lambda) B^T, matrix SVD, independent determinants",
    physical_cosine_not_a_constraint=TRUE),file.path(out,"gradient_direct_validation.json"),
    auto_unbox=TRUE,pretty=TRUE,digits=17)
  ans
}

lv_case_audit <- function(path) {
  x<-readRDS(path);stopifnot(identical(x$status,"ok"))
  dg<-x$diag;A<-dg$core$A;B<-dg$core$B
  stopifnot(identical(dg$configuration$domain,"lambda-axis"),
            isFALSE(dg$configuration$physical_cosine_is_constraint))
  rows<-list()
  add<-function(point,kind,id,center=NULL,angle_limit=NULL) {
    row<-lv_point_check(A,B,point,center,angle_limit)
    row$kind<-kind;row$id<-id
    if(!"angle_violation"%in%names(row))row$angle_violation<-NA_real_
    rows[[length(rows)+1L]]<<-row
  }
  for(i in seq_along(dg$rho_local_minima))add(dg$rho_local_minima[[i]],"rho_local",i)
  for(i in seq_along(dg$E_local_minima))add(dg$E_local_minima[[i]],"E_local",i)
  for(i in seq_along(x$cap_details)) {
    detail<-x$cap_details[[i]]
    add(detail$center,"cap_center",i)
    add(detail$selected5,"cap5",i,detail$center,5)
    add(detail$selected10,"cap10",i,detail$center,10)
    if(isTRUE(detail$refined_without_angle$feasible))add(detail$refined_without_angle,"refined",i)
    stopifnot(detail$selected10$E1<=detail$selected5$E1+1e-8,
              detail$selected5$E1<=detail$center$E1+1e-8)
    if(!is.null(detail$local_best)) {
      add(detail$local_best,"local_within10",i,detail$center,10)
      stopifnot(isTRUE(detail$local_best$success))
    }
  }
  result<-do.call(rbind,rows)
  for(name in c("p","q","regime","seed","T"))result[[name]]<-x$job[[name]]
  pairs<-lv_pair_pool_check(A,B,dg$rho_local_minima,dg$E_local_minima)
  if(!is.null(x$primary)) {
    p<-x$primary
    if(is.finite(pairs$Q))stopifnot(abs(p$pair_joint_max-pairs$Q)<1e-7)else stopifnot(is.infinite(p$pair_joint_max))
    stopifnot(p$n_pairs_angle_le5==pairs$n_eligible)
    closest<-min(vapply(dg$E_local_minima,function(e)lv_angle(dg$lambda_rho,e$lambda),numeric(1L)))
    stopifnot(abs(p$center_to_E1_nearest_deg_090-closest)<1e-7)
  }
  if(!is.null(x$source_path)) {
    stopifnot(unname(tools::md5sum(x$source_path))==x$source_md5)
    old<-readRDS(x$source_path)
    stopifnot(identical(A,old$geometry$diag$core$A),identical(B,old$geometry$diag$core$B))
    stopifnot(abs(dg$J_min-old$geometry$diag$J_min)<1e-12,
              abs(dg$J_relative_min-old$geometry$diag$J_relative_min)<1e-12)
  }
  result
}

lv_cases_audit <- function(paths,out="work/lambda_geometry_validation") {
  stopifnot(length(paths)>0L)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  tab<-do.call(rbind,lapply(paths,lv_case_audit))
  write.csv(tab,file.path(out,"cases_direct_validation.csv"),row.names=FALSE)
  jsonlite::write_json(list(status="passed",cases=length(paths),points=nrow(tab),
    maximum_E1_error=max(tab$E1_error,na.rm=TRUE),maximum_rho_error=max(tab$rho_error,na.rm=TRUE),
    maximum_norm_error=max(tab$matrix_norm_error),maximum_axis_violation=max(tab$axis_violation),
    maximum_reported_axis_cosine_error=max(tab$reported_axis_error,na.rm=TRUE),
    maximum_reported_physical_cosine_error=max(tab$reported_physical_error,na.rm=TRUE),
    maximum_reported_axis_angle_error=max(tab$reported_axis_angle_error,na.rm=TRUE),
    maximum_angle_violation=max(tab$angle_violation,na.rm=TRUE),
    physical_cosine_above_old_cap_count=sum(tab$physical_max_cosine>.9+2e-6),
    Q_reconstructed_from_all_independent_pairs=TRUE,unchanged_QMLE_core_and_spectra=TRUE,
    source_checkpoint_hashes_unchanged=TRUE),file.path(out,"cases_direct_validation.json"),
    auto_unbox=TRUE,pretty=TRUE,digits=17)
  tab
}

lv_stratified_audit <- function(case_dir="outputs/lambda_axis_geometry_20261030/cases",
                                out="outputs/lambda_axis_geometry_20261030/validation",
                                require_complete=TRUE) {
  paths<-list.files(case_dir,pattern="[.]rds$",full.names=TRUE)
  fields<-strcapture("p([0-9]+)_q([0-9]+)_(interior|boundary|exterior)_s([0-9]+)_T([0-9]+)[.]rds$",
    basename(paths),proto=data.frame(p=integer(),q=integer(),regime=character(),seed=integer(),T=integer()))
  stopifnot(!anyNA(fields));fields$path<-paths
  fields<-fields[order(fields$p,fields$q,fields$regime,fields$T,fields$seed),]
  selected<-fields[!duplicated(fields[c("p","q","regime","T")]),]
  if(require_complete)stopifnot(nrow(selected)==108L)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(selected,file.path(out,"stratified_direct_manifest.csv"),row.names=FALSE)
  lv_cases_audit(selected$path,out)
}

lv_truth_reference <- function(out="outputs/lambda_axis_geometry_20261030/validation") {
  source("work/R_design.R",encoding="UTF-8")
  source("work/R_diagnostics.R",encoding="UTF-8")
  source("work/R_grid_geometry.R",encoding="UTF-8")
  source("work/R_local_alignment_core.R",encoding="UTF-8")
  source("work/R_lambda_geometry_engine.R",encoding="UTF-8")
  cfg<-list(axis_cosine_cap=.90,diagnostic_random_starts=30L,
    diagnostic_maxeval=400L,diagnostic_cross_rounds=2L,
    diagnostic_cross_start_limit=5L,diagnostic_old_start_limit=0L,
    alignment_starts=5L,alignment_maxeval=300L,refine_maxeval=400L)
  rows<-list();details<-list()
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(regime in c("interior","boundary","exterior")) {
    truth<-core_loading_design(regime,n=.5)
    core<-list(A=truth$Ar_unit,B=truth$Br_unit)
    record<-list(diag=list(core=core))
    if(regime=="exterior") {
      known<-as.data.frame(as.list(setNames(truth$extra_lambda,sprintf("lambda_%02d",seq_along(truth$extra_lambda)))))
      record$geometry<-list(known_directions=known)
    }
    result<-lg_search(record,593071L+match(regime,c("interior","boundary","exterior")),cfg)
    dg<-result$diag;br<-dg$rho_local_minima[[1L]];be<-dg$E_local_minima[[1L]]
    # The numerical reference is checked by literal matrix reconstruction,
    # without using C^T C to evaluate the saved objectives.
    lv_point_check(core$A,core$B,br);lv_point_check(core$A,core$B,be)
    row<-data.frame(regime=regime,d=ncol(core$A),rho_min=br$rho,E1_at_rho=br$E1,
      E1_min=be$E1,rho_at_E1=be$rho,
      rho_axis_cosine=br$axis_max_cosine,E1_axis_cosine=be$axis_max_cosine,
      rho_physical_cosine=br$physical_max_cosine,E1_physical_cosine=be$physical_max_cosine,
      rho_axis_angle_deg=br$axis_min_angle_deg,E1_axis_angle_deg=be$axis_min_angle_deg,
      rho_to_E1_angle_deg=lv_angle(br$lambda,be$lambda),
      E1_cap5=result$cap_rows$E1_cap5_saved[1L],E1_cap10=result$cap_rows$E1_cap10[1L],
      rho_success=br$success,E1_success=be$success,
      rho_local_count=length(dg$rho_local_minima),E1_local_count=length(dg$E_local_minima),
      random_starts=30L,cross_rounds=2L,maxeval=400L)
    for(j in seq_len(8L)) {
      row[[sprintf("rho_lambda_%02d",j)]]<-if(j<=ncol(core$A))br$lambda[j]else NA_real_
      row[[sprintf("E1_lambda_%02d",j)]]<-if(j<=ncol(core$A))be$lambda[j]else NA_real_
    }
    rows[[regime]]<-row;details[[regime]]<-result
  }
  tab<-do.call(rbind,rows);rownames(tab)<-NULL
  write.csv(tab,file.path(out,"population_lambda_search_reference.csv"),row.names=FALSE)
  saveRDS(details,file.path(out,"population_lambda_search_reference.rds"),version=2L)
  jsonlite::write_json(list(domain="lambda-axis cosine <= .90",axis_exclusion_deg=acos(.9)*180/pi,
    data="exact unit-column truth cores from R_design.R, n=.5",
    search_configuration=cfg,numerical_reference_only=TRUE,global_optimality_proved=FALSE,
    identified_cases_strict_positive="Compact fixed lambda-axis domain, injective Z, and only coordinate rank-one directions imply strictly positive E1/rho minima.",
    QMLE_optimizations=0L),file.path(out,"population_lambda_search_reference.json"),
    auto_unbox=TRUE,pretty=TRUE,digits=17)
  tab
}

if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  out<-if(length(args))args[1L]else"work/lambda_geometry_validation"
  tab<-lv_population(out)
  cat("Population/direct-domain checks passed:",nrow(tab),"\n")
}
