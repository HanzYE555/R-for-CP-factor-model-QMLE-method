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
