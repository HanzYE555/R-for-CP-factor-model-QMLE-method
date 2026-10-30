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
