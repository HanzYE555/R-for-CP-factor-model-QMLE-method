# Lambda-axis geometry on saved QMLE cores. No QMLE fitting or file I/O.
# Source R_diagnostics.R, R_grid_geometry.R and R_local_alignment_core.R first.
# Coordinates x satisfy ||x||=||M(lambda)||_F=1 and lambda=Rinv*x.
# max_cosine is deliberately the lambda-axis cosine for old cap-helper
# compatibility; physical_max_cosine is a separate reported quantity.

lg_get <- function(cfg, name, default) if (is.null(cfg[[name]])) default else cfg[[name]]

lg_problem <- function(core, c=.90, rank_tolerance=1e-10) {
  pr <- cp_diagnostic_problem(core, c, rank_tolerance)
  d <- pr$d
  if (!is.finite(c) || c<=1/sqrt(d) || c>=1) stop("Invalid nonempty lambda-axis cosine cutoff.")
  L <- pr$Rinv; C <- pr$C; ii <- pr$pairs[1L,]; jj <- pr$pairs[2L,]
  CC <- crossprod(C); CC <- (CC+t(CC))/2
  ec <- new.env(parent=emptyenv()); rc <- new.env(parent=emptyenv())
  old_E <- pr$E_squared; physical_G <- pr$G
  E_squared <- function(x) {
    if (exists("x",ec,inherits=FALSE) && identical(x,ec$x)) return(ec$value)
    value <- old_E(x); ec$x <- x; ec$value <- value; value
  }
  rho_squared <- function(x) {
    if (exists("x",rc,inherits=FALSE) && identical(x,rc$x)) return(rc$value)
    lambda <- drop(L %*% x); k <- lambda[ii]*lambda[jj]
    den <- sum(k*k)
    if (!is.finite(den) || den<1e-28) return(list(value=1e10,gradient=rep(0,d)))
    u <- drop(CC %*% k); num <- sum(k*u)
    # Direct minors protect near-kernel values from Gram cancellation.
    if (!is.finite(num) || num<1e-12*max(1,den)) {
      h <- drop(C %*% k); num <- sum(h*h); u <- drop(crossprod(C,h))
    }
    gk <- 2*(u*den-k*num)/den^2
    gl <- numeric(d)
    for (j in seq_along(ii)) {
      gl[ii[j]] <- gl[ii[j]]+gk[j]*lambda[jj[j]]
      gl[jj[j]] <- gl[jj[j]]+gk[j]*lambda[ii[j]]
    }
    value <- list(value=max(0,num/den),gradient=drop(crossprod(L,gl)))
    rc$x <- x; rc$value <- value; value
  }
  g_ineq <- function(x) {
    lambda <- drop(L %*% x)
    lambda^2-c^2*sum(lambda^2)
  }
  jac_ineq <- function(x) {
    lambda <- drop(L %*% x)
    2*sweep(L,1L,lambda,"*")-
      matrix(2*c^2*drop(crossprod(L,lambda)),d,d,byrow=TRUE)
  }
  summarize <- function(x) {
    nr <- sqrt(sum(x*x))
    if (!is.finite(nr) || nr<=0) stop("Cannot summarize a zero/nonfinite point.")
    x <- as.numeric(x/nr); lambda <- drop(L %*% x)
    axis_cos <- min(1,max(abs(lambda))/sqrt(sum(lambda^2)))
    list(x=x,lambda=lambda,rho=sqrt(max(0,rho_squared(x)$value)),
         E1=sqrt(max(0,E_squared(x)$value)),max_cosine=axis_cos,
         axis_max_cosine=axis_cos,physical_max_cosine=max(abs(physical_G%*%x)),
         axis_min_angle_deg=acos(axis_cos)*180/pi,matrix_norm=1)
  }
  pr$rho_squared <- rho_squared; pr$E_squared <- E_squared
  pr$g_ineq <- g_ineq; pr$jac_ineq <- jac_ineq; pr$summarize <- summarize
  pr$axis_cosine_cap <- c; pr$cosine_cap <- c
  pr$axis_angle_deg <- acos(c)*180/pi; pr$domain <- "lambda-axis"
  pr$C_gram <- CC
  pr
}

lg_record_diag <- function(record) {
  if (!is.null(record$geometry$diag)) return(record$geometry$diag)
  if (!is.null(record$diagnostics)) return(record$diagnostics)
  if (!is.null(record$diag)) return(record$diag)
  stop("Record has no saved diagnostic core.")
}

lg_seed <- function(seed) set.seed(as.integer((abs(as.double(seed)) %% 2147483646)+1))

lg_usable <- function(point) isTRUE(point$feasible) &&
  is.finite(point$rho) && is.finite(point$E1) &&
  (point$status>0L || point$status==-4L)

lg_optimize <- function(pr,x0,objective,source,maxeval=400L,cfg=list()) {
  fun <- if (objective=="rho") pr$rho_squared else pr$E_squared
  fit <- tryCatch(nloptr::nloptr(x0=x0,
    eval_f=function(x)fun(x)$value,eval_grad_f=function(x)fun(x)$gradient,
    eval_g_ineq=pr$g_ineq,eval_jac_g_ineq=pr$jac_ineq,
    eval_g_eq=function(x)sum(x*x)-1,
    eval_jac_g_eq=function(x)matrix(2*x,nrow=1L),
    opts=list(algorithm="NLOPT_LD_SLSQP",maxeval=as.integer(maxeval),
      xtol_rel=lg_get(cfg,"diagnostic_xtol_rel",1e-9),
      ftol_abs=lg_get(cfg,"diagnostic_ftol_abs",1e-13),print_level=0L,
      tol_constraints_ineq=rep(1e-10,pr$d),tol_constraints_eq=1e-10)),
    error=function(e)list(solution=x0,status=-999L,message=conditionMessage(e),iterations=0L))
  raw_norm <- sqrt(sum(fit$solution^2)); tol <- lg_get(cfg,"diagnostic_feasibility_tolerance",2e-6)
  point <- tryCatch(pr$summarize(fit$solution),error=function(e)NULL)
  if (is.null(point)) point <- list(x=rep(NA_real_,pr$d),lambda=rep(NA_real_,pr$d),
    rho=NA_real_,E1=NA_real_,max_cosine=Inf,axis_max_cosine=Inf,
    physical_max_cosine=NA_real_,axis_min_angle_deg=NA_real_,matrix_norm=NA_real_)
  point$feasible <- is.finite(raw_norm) && abs(raw_norm-1)<=tol &&
    is.finite(point$max_cosine) && point$max_cosine<=pr$cosine_cap+tol
  point$eq_residual <- abs(raw_norm^2-1)
  point$objective_value <- if(objective=="rho")point$rho else point$E1
  point$optimized_objective_squared <- point$objective_value^2
  point$success <- isTRUE(fit$status%in%1:4)
  point$status <- as.integer(fit$status);point$message <- as.character(fit$message)
  point$source <- source;point$objective <- objective
  point$iterations <- if(is.null(fit$iterations))NA_integer_ else as.integer(fit$iterations)
  point$axis_constraint_edge <- point$max_cosine>=pr$cosine_cap-1e-5
  point
}

lg_cluster <- function(runs,objective,tolerance=1e-5) {
  valid <- Filter(lg_usable,runs)
  if (!length(valid)) return(list())
  valid <- valid[order(vapply(valid,function(x)x$objective_value,numeric(1L)))]
  out <- list()
  for (point in valid) {
    unit <- point$lambda/sqrt(sum(point$lambda^2))
    similarities <- if(length(out))vapply(out,function(q)
      abs(sum(unit*q$lambda))/sqrt(sum(q$lambda^2)),numeric(1L))else numeric()
    hits <- which(similarities>=1-tolerance)
    if (!length(hits)) {
      point$hit_count <- 1L;point$converged_hit_count <- as.integer(point$success)
      point$all_hit_sources <- point$source
      point$truth_guided_hit <- grepl("^truth",point$source)
      out[[length(out)+1L]] <- point
    } else {
      j <- hits[which.max(similarities[hits])]
      out[[j]]$hit_count <- out[[j]]$hit_count+1L
      out[[j]]$converged_hit_count <- out[[j]]$converged_hit_count+as.integer(point$success)
      out[[j]]$all_hit_sources <- paste(out[[j]]$all_hit_sources,point$source,sep=";")
      out[[j]]$truth_guided_hit <- out[[j]]$truth_guided_hit||grepl("^truth",point$source)
    }
  }
  for (j in seq_along(out)) {out[[j]]$minimum_id<-j;out[[j]]$objective<-objective}
  out
}

lg_run_rows <- function(runs,objective=NULL,extra=list()) {
  if (!length(runs)) return(data.frame())
  rows <- lapply(seq_along(runs),function(i) {
    x <- runs[[i]]
    z <- data.frame(objective=if(is.null(objective))x$objective else objective,
      run=i,source=x$source,feasible=isTRUE(x$feasible),converged=isTRUE(x$success),
      status=x$status,message=if(is.null(x$message))""else x$message,
      rho=x$rho,E1=x$E1,max_axis_cosine=x$max_cosine,
      max_lambda_cosine=x$max_cosine,max_loading_cosine=x$physical_max_cosine,
      physical_max_cosine=x$physical_max_cosine,axis_min_angle_deg=x$axis_min_angle_deg,
      matrix_norm=x$matrix_norm,eq_residual=if(is.null(x$eq_residual))NA_real_ else x$eq_residual,
      iterations=if(is.null(x$iterations))NA_integer_ else x$iterations,
      axis_constraint_edge=isTRUE(x$axis_constraint_edge),stringsAsFactors=FALSE)
    for (j in seq_len(8L))z[[sprintf("lambda_%02d",j)]] <- if(j<=length(x$lambda))x$lambda[j]else NA_real_
    for (name in names(extra))z[[name]] <- extra[[name]]
    z
  })
  do.call(rbind,rows)
}

lg_independent <- function(record,pr,seed,cfg=list()) {
  stopifnot(requireNamespace("nloptr",quietly=TRUE))
  had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  old_seed <- if(had)get(".Random.seed",.GlobalEnv)else NULL
  on.exit(if(had)assign(".Random.seed",old_seed,.GlobalEnv)else
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  lg_seed(seed)
  nrandom <- as.integer(lg_get(cfg,"diagnostic_random_starts",30L))
  maxeval <- as.integer(lg_get(cfg,"diagnostic_maxeval",400L))
  rounds <- as.integer(lg_get(cfg,"diagnostic_cross_rounds",2L))
  limit <- as.integer(lg_get(cfg,"diagnostic_cross_start_limit",5L))
  old_limit <- as.integer(lg_get(cfg,"diagnostic_old_start_limit",5L))
  cluster_tol <- lg_get(cfg,"diagnostic_cluster_cosine_tolerance",1e-5)
  starts <- list(); labels <- character(); kinds <- character()
  add <- function(lambda,label,kind) {
    if(length(lambda)!=pr$d || any(!is.finite(lambda)) || sum(lambda^2)<=0)return(FALSE)
    x <- cp_alignment_x_from_lambda(pr,lambda);p <- pr$summarize(x)
    if(p$max_cosine>pr$cosine_cap+1e-10)return(FALSE)
    u <- p$lambda/sqrt(sum(p$lambda^2))
    if(length(starts)&&any(vapply(starts,function(y){v<-drop(pr$Rinv%*%y);
      abs(sum(u*v))/sqrt(sum(v*v))>=1-1e-12},logical(1L))))return(FALSE)
    starts[[length(starts)+1L]] <<- x;labels <<- c(labels,label);kinds <<- c(kinds,kind)
    TRUE
  }
  known <- record$geometry$known_directions
  if(!is.null(known)&&nrow(known))for(i in seq_len(nrow(known))) {
    fields <- sprintf("lambda_%02d",seq_len(pr$d))
    if(all(fields%in%names(known)))add(as.numeric(unlist(known[i,fields,drop=FALSE],use.names=FALSE)),paste0("truth_projection_",i),"truth")
  }
  truth_n <- sum(kinds=="truth")
  old <- lg_record_diag(record);nold <- 0L
  for(i in seq_len(max(length(old$rho_local_minima),length(old$E_local_minima),0L))) {
    for(key in c("rho_local_minima","E_local_minima"))if(nold<old_limit&&i<=length(old[[key]])) {
      if(add(old[[key]][[i]]$lambda,paste0("old_",key,"_",i,"_reoptimized"),"old"))nold<-nold+1L
    }
    if(nold>=old_limit)break
  }
  # Draw in the same unit-M QR coordinates as the original search. Extras
  # are bounded warm starts and do not replace any of the random starts.
  attempts <- 0L;nr <- 0L
  while(nr<nrandom&&attempts<max(10000L,1000L*nrandom)) {
    x <- rnorm(pr$d);x<-x/sqrt(sum(x*x));attempts<-attempts+1L
    lambda <- drop(pr$Rinv%*%x)
    if(max(abs(lambda))/sqrt(sum(lambda^2))<pr$cosine_cap-1e-9) {
      starts[[length(starts)+1L]]<-x;nr<-nr+1L
      labels<-c(labels,sprintf("random_%03d",nr));kinds<-c(kinds,"random")
    }
  }
  if(nr<nrandom)stop("Could not generate requested feasible lambda-axis starts.")
  rho_runs <- lapply(seq_along(starts),function(i)lg_optimize(pr,starts[[i]],"rho",labels[i],maxeval,cfg))
  rho_local <- lg_cluster(rho_runs,"rho",cluster_tol)
  E_starts <- starts;E_labels<-labels
  if(length(rho_local)&&limit>0L)for(j in seq_len(min(length(rho_local),limit))) {
    E_starts[[length(E_starts)+1L]]<-rho_local[[j]]$x
    E_labels<-c(E_labels,paste0("initial_rho_candidate_",j,"_to_E1"))
  }
  E_runs <- lapply(seq_along(E_starts),function(i)lg_optimize(pr,E_starts[[i]],"E1",E_labels[i],maxeval,cfg))
  E_local <- lg_cluster(E_runs,"E1",cluster_tol)
  if(rounds>0L)for(round in seq_len(rounds)) {
    if(length(E_local)&&limit>0L) {
      extra <- lapply(seq_len(min(length(E_local),limit)),function(j)
        lg_optimize(pr,E_local[[j]]$x,"rho",paste0("round_",round,"_E1_",j,"_to_rho"),maxeval,cfg))
      rho_runs<-c(rho_runs,extra);rho_local<-lg_cluster(rho_runs,"rho",cluster_tol)
    }
    if(length(rho_local)&&limit>0L) {
      extra<-lapply(seq_len(min(length(rho_local),limit)),function(j)
        lg_optimize(pr,rho_local[[j]]$x,"E1",paste0("round_",round,"_rho_",j,"_to_E1"),maxeval,cfg))
      E_runs<-c(E_runs,extra);E_local<-lg_cluster(E_runs,"E1",cluster_tol)
    }
  }
  if(!length(rho_local)||!length(E_local))stop("No feasible optimized lambda-axis candidates.")
  br<-rho_local[[1L]];be<-E_local[[1L]]
  witnesses<-Filter(function(x)isTRUE(x$feasible)&&is.finite(x$rho)&&is.finite(x$E1),c(rho_runs,E_runs))
  ew<-witnesses[[which.min(vapply(witnesses,function(x)x$E1,numeric(1L)))]]
  rw<-witnesses[[which.min(vapply(witnesses,function(x)x$rho,numeric(1L)))]]
  spectra<-if(!is.null(old$J_orig))list(J_min=old$J_min,J_relative_min=old$J_relative_min,
    orig=old$J_orig,nor=old$J_nor,negative_roundoff=old$J_negative_roundoff)else
    cp_jacobian_spectra(pr$A,pr$B,pr$C,pr$pairs)
  dg<-list(rho_min=br$rho,E_at_rho=br$E1,E_min=be$E1,rho_at_E=be$rho,
    lambda_rho=br$lambda,lambda_E=be$lambda,
    cosine_at_rho=br$max_cosine,cosine_at_E=be$max_cosine,
    axis_cosine_at_rho=br$max_cosine,axis_cosine_at_E=be$max_cosine,
    physical_cosine_at_rho=br$physical_max_cosine,physical_cosine_at_E=be$physical_max_cosine,
    rho_best_success=br$success,E_best_success=be$success,
    E1_feasible_upper_bound=ew$E1,rho_feasible_upper_bound=rw$rho,
    E1_search_consistency_gap=max(0,be$E1-ew$E1),rho_search_consistency_gap=max(0,br$rho-rw$rho),
    feasible_witnesses=list(E1=ew,rho=rw),
    rho_local_minima=rho_local,E_local_minima=E_local,
    n_rho_local_minima=length(rho_local),n_E_local_minima=length(E_local),
    rho_search_runs=length(rho_runs),E_search_runs=length(E_runs),
    rho_search_feasible=sum(vapply(rho_runs,function(x)isTRUE(x$feasible),logical(1L))),
    E_search_feasible=sum(vapply(E_runs,function(x)isTRUE(x$feasible),logical(1L))),
    rho_search_converged=sum(vapply(rho_runs,function(x)isTRUE(x$success),logical(1L))),
    E_search_converged=sum(vapply(E_runs,function(x)isTRUE(x$success),logical(1L))),
    J_min=spectra$J_min,J_relative_min=spectra$J_relative_min,J_orig=spectra$orig,J_nor=spectra$nor,
    J_negative_roundoff=spectra$negative_roundoff,Z_singular_values=pr$Z_singular_values,
    C_singular_values=old$C_singular_values,core=old$core,
    configuration=list(seed=seed,domain="lambda-axis",axis_cosine_cap=pr$cosine_cap,
      axis_angle_deg=pr$axis_angle_deg,random_starts=nrandom,truth_projection_starts=truth_n,
      old_feasible_warm_starts=nold,maxeval=maxeval,cross_start_limit=limit,cross_refinement_rounds=rounds,
      cluster_metric="projective lambda angle",cluster_cosine_tolerance=cluster_tol,
      physical_cosine_is_constraint=FALSE,minimum_note="Feasible numerical optimized candidates; no exhaustive/local-optimality certificate."))
  list(diag=dg,all_search_runs=rbind(lg_run_rows(rho_runs,"rho"),lg_run_rows(E_runs,"E1")))
}

lg_cap_point <- function(pr,lambda,center,angle,source,status=0L,success=FALSE) {
  p<-pr$summarize(cp_alignment_x_from_lambda(pr,lambda))
  p$angle_deg<-cp_alignment_angle_deg(p$lambda,center$lambda)
  p$matrix_angle_deg<-cp_alignment_angle_deg(p$x,center$x)
  p$feasible<-p$max_cosine<=pr$cosine_cap+2e-6&&p$angle_deg<=angle+1e-4
  p$eq_residual<-0;p$source<-source;p$status<-as.integer(status);p$success<-isTRUE(success)
  p$message<-"Retained feasible witness; not an additional optimizer run."
  p$objective<-"E1";p$iterations<-0L
  p$cap_edge<-p$angle_deg>=angle-.05
  p$axis_constraint_edge<-p$max_cosine>=pr$cosine_cap-1e-5
  p
}

lg_choose_cap <- function(points) {
  points<-Filter(function(x)isTRUE(x$feasible)&&is.finite(x$E1),points)
  if(!length(points))stop("No feasible cap point.")
  values<-vapply(points,function(x)x$E1,numeric(1L))
  near<-which(values<=min(values)+1e-10)
  successful<-near[vapply(points[near],function(x)isTRUE(x$success),logical(1L))]
  pool<-if(length(successful))successful else seq_along(points)
  points[[pool[which.min(values[pool])]]]
}

lg_one_cap <- function(pr,center,Epool,angle,cfg,warm=NULL) {
  in_cap<-which(vapply(Epool,function(e)cp_alignment_angle_deg(center$lambda,e$lambda)<=angle+1e-4,logical(1L)))
  witnesses<-lapply(in_cap,function(j)lg_cap_point(pr,Epool[[j]]$lambda,center,angle,
    paste0("independent_E1_witness_",j),Epool[[j]]$status,Epool[[j]]$success))
  extra<-lapply(Epool[in_cap],function(e)e$lambda)
  if(length(extra))extra<-extra[order(vapply(Epool[in_cap],function(e)e$E1,numeric(1L)))]
  warm_point<-NULL
  if(!is.null(warm)) {
    warm_point<-lg_cap_point(pr,warm$lambda,center,angle,"cap5_feasible_witness")
    extra<-c(list(warm$lambda),extra)
  }
  search<-cp_alignment_cap_search(pr,center$lambda,angle_deg=angle,
    nstarts=as.integer(lg_get(cfg,"alignment_starts",5L)),
    maxeval=as.integer(lg_get(cfg,"alignment_maxeval",300L)),extra_lambdas=extra)
  for(i in seq_along(search$runs)) {
    search$runs[[i]]$source<-paste0("cap",angle,"_",search$runs[[i]]$source)
    search$runs[[i]]$objective<-"E1"
    search$runs[[i]]$axis_constraint_edge<-search$runs[[i]]$loading_constraint_edge
    search$runs[[i]]$loading_constraint_edge<-NULL
  }
  baseline<-lg_cap_point(pr,center$lambda,center,angle,"center_feasible_witness")
  points<-c(search$runs,list(baseline),witnesses,if(!is.null(warm_point))list(warm_point)else list())
  best<-lg_choose_cap(points)
  list(best=best,runs=search$runs,witnesses=witnesses,baseline=baseline,
       previous5=warm_point,search_configuration=search$configuration)
}

lg_caps <- function(pr,dg,seed,cfg=list()) {
  had<-exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  old_seed<-if(had)get(".Random.seed",.GlobalEnv)else NULL
  on.exit(if(had)assign(".Random.seed",old_seed,.GlobalEnv)else
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  rows<-list();details<-list();run_rows<-list();Epool<-dg$E_local_minima
  for(i in seq_along(dg$rho_local_minima)) {
    lg_seed(as.double(seed)+i*1009)
    old<-dg$rho_local_minima[[i]]
    center<-pr$summarize(cp_alignment_x_from_lambda(pr,old$lambda))
    five<-lg_one_cap(pr,center,Epool,5,cfg)
    ten<-lg_one_cap(pr,center,Epool,10,cfg,five$best)
    b5<-five$best;b10<-ten$best
    if(b10$E1>b5$E1+1e-9||b5$E1>center$E1+1e-9)stop("Cap witness nesting failed.")
    refined<-lg_optimize(pr,b10$x,"E1","uncapped_E1_refinement_from_cap10",
      as.integer(lg_get(cfg,"refine_maxeval",400L)),cfg)
    refined$angle_deg<-cp_alignment_angle_deg(refined$lambda,center$lambda)
    refined$matrix_angle_deg<-cp_alignment_angle_deg(refined$x,center$x)
    refined$domain_feasible<-refined$feasible
    refined$inside10<-isTRUE(refined$domain_feasible)&&is.finite(refined$angle_deg)&&refined$angle_deg<=10+1e-4
    local<-Filter(function(p)isTRUE(p$success)&&isTRUE(p$feasible),ten$witnesses)
    if(isTRUE(refined$success)&&isTRUE(refined$inside10))local<-c(local,list(refined))
    lb<-if(length(local))local[[which.min(vapply(local,function(e)e$E1,numeric(1L)))]]else NULL
    e_angles<-vapply(Epool,function(e)cp_alignment_angle_deg(center$lambda,e$lambda),numeric(1L))
    r<-data.frame(rho_order=i,rho_center=center$rho,E1_center=center$E1,
      E1_cap5_saved=b5$E1,rho_cap5_saved=b5$rho,center_to_cap5_deg=b5$angle_deg,
      E1_cap10=b10$E1,rho_cap10=b10$rho,center_to_cap10_deg=b10$angle_deg,
      center_to_cap10_matrix_deg=b10$matrix_angle_deg,
      center_max_axis_cosine=center$max_cosine,center_physical_max_cosine=center$physical_max_cosine,
      cap5_max_axis_cosine=b5$max_cosine,cap5_physical_max_cosine=b5$physical_max_cosine,
      cap10_max_axis_cosine=b10$max_cosine,cap10_physical_max_cosine=b10$physical_max_cosine,
      cap5_max_lambda_cosine=b5$max_cosine,cap5_max_loading_cosine=b5$physical_max_cosine,
      cap10_max_lambda_cosine=b10$max_cosine,cap10_max_loading_cosine=b10$physical_max_cosine,
      cap5_axis_min_angle_deg=b5$axis_min_angle_deg,cap10_axis_min_angle_deg=b10$axis_min_angle_deg,
      cap5_angle_edge=isTRUE(b5$cap_edge),cap5_axis_edge=isTRUE(b5$axis_constraint_edge),
      cap10_angle_edge=isTRUE(b10$cap_edge),cap10_axis_edge=isTRUE(b10$axis_constraint_edge),
      cap10_selected_converged=isTRUE(b10$success),cap10_selected_status=b10$status,
      cap10_selected_source=b10$source,cap10_optimizer_runs=length(ten$runs),
      cap10_optimizer_success_count=sum(vapply(ten$runs,function(e)isTRUE(e$success)&&isTRUE(e$feasible),logical(1L))),
      cap10_fallback_selected=b10$source%in%c("center_feasible_witness","cap5_feasible_witness"),
      E1_min_saved=dg$E_min,rho_at_E1_saved=dg$rho_at_E,
      center_to_nearest_saved_E1_deg=min(e_angles),saved_E1_within10_exists=length(ten$witnesses)>0L,
      E1_refined_from10=refined$E1,rho_refined_from10=refined$rho,
      refined_center_angle_deg=refined$angle_deg,refined_success=isTRUE(refined$success),
      refined_status=refined$status,refined_domain_feasible=isTRUE(refined$domain_feasible),
      refined_inside10=isTRUE(refined$inside10),local_candidate_exists10=!is.null(lb),
      E1_local_within10_found=if(is.null(lb))NA_real_ else lb$E1,
      local_best_angle10=if(is.null(lb))NA_real_ else lb$angle_deg,
      selected10_to_local_deg=if(is.null(lb))NA_real_ else cp_alignment_angle_deg(b10$lambda,lb$lambda),
      E1_gain_center_to10=center$E1-b10$E1,E1_gain5_to10=b5$E1-b10$E1,
      search_gap_vs_saved_E1=max(0,dg$E_min-b10$E1),stringsAsFactors=FALSE)
    for(j in seq_len(8L))for(key in c("center","cap5","cap10","refined")) {
      p<-switch(key,center=center,cap5=b5,cap10=b10,refined=refined)
      r[[sprintf("%s_lambda_%02d",key,j)]]<-if(j<=pr$d)p$lambda[j]else NA_real_
    }
    rows[[i]]<-r;details[[i]]<-list(rho_order=i,center=center,cap5=five,cap10=ten,
      refined_without_angle=refined,selected5=b5,selected10=b10,local_best=lb)
    run_rows[[length(run_rows)+1L]]<-lg_run_rows(five$runs,"E1",list(rho_order=i,phase="cap5"))
    run_rows[[length(run_rows)+1L]]<-lg_run_rows(ten$runs,"E1",list(rho_order=i,phase="cap10"))
    run_rows[[length(run_rows)+1L]]<-lg_run_rows(list(refined),"E1",list(rho_order=i,phase="uncapped_refinement"))
  }
  list(rows=do.call(rbind,rows),details=details,all_cap_runs=do.call(rbind,run_rows))
}

lg_search <- function(record,jobseed,cfg=list()) {
  old<-lg_record_diag(record);core<-old$core
  if(is.null(core)&&!is.null(record$Ar)&&!is.null(record$Br))core<-list(A=record$Ar,B=record$Br)
  if(is.null(core))stop("Saved aligned, unit-column diagnostic core is required.")
  cap<-lg_get(cfg,"axis_cosine_cap",lg_get(cfg,"diagnostic_cosine_cap",.90))
  pr<-lg_problem(core,cap,lg_get(cfg,"rank_tolerance",1e-10))
  independent<-lg_independent(record,pr,jobseed,cfg)
  caps<-lg_caps(pr,independent$diag,as.double(jobseed)+1907,cfg)
  pairs<-cp_grid_candidate_pairs(independent$diag,pr)
  list(diag=independent$diag,cap_rows=caps$rows,cap_details=caps$details,
       pair_distances=pairs,all_search_runs=independent$all_search_runs,
       all_cap_runs=caps$all_cap_runs,
       configuration=list(domain="lambda-axis",axis_cosine_cap=cap,
         axis_angle_deg=acos(cap)*180/pi,physical_cosine_is_constraint=FALSE,
         n_QMLE_optimizations=0L,source="saved aligned unit-column cores"))
}
