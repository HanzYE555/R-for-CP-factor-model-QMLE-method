# Truth-near multiple starts for the p,q grid simulation.
# Requires R_qmle_backend.R. Only base R is used here.
# These starts use simulation truth and are NOT a blind-data estimator.
# The objective still estimates every diagonal noise variance by default.

cp_grid_option <- function(cfg, name, default) {
  value <- cfg[[name]]
  if (is.null(value)) default else value
}

cp_grid_perturb_directions <- function(A, max_angle_deg=3) {
  A <- cp_normalize_columns(as.matrix(A))
  if (!is.finite(max_angle_deg) || max_angle_deg<0 || max_angle_deg>=90)
    stop("truth_angle_deg must be in [0,90).")
  if (max_angle_deg==0) return(A)
  for(j in seq_len(ncol(A))) {
    unit <- A[,j]
    direction <- stats::rnorm(nrow(A))
    direction <- direction-sum(direction*unit)*unit
    norm <- sqrt(sum(direction^2))
    if (!is.finite(norm) || norm<1e-14) {
      direction <- numeric(nrow(A));direction[which.min(abs(unit))] <- 1
      direction <- direction-sum(direction*unit)*unit
      norm <- sqrt(sum(direction^2))
    }
    theta <- stats::runif(1L,0,max_angle_deg)*pi/180
    A[,j] <- cos(theta)*unit+sin(theta)*direction/norm
  }
  A
}

cp_grid_pack_exact <- function(obj,A,B,M,psi=rep(1,obj$N),repair_floor=1e-10) {
  A <- as.matrix(A);B <- as.matrix(B);M <- as.matrix(M)
  if (!identical(dim(A),c(obj$p,obj$d)) ||
      !identical(dim(B),c(obj$q,obj$d)) ||
      !identical(dim(M),c(obj$d,obj$d))) stop("Invalid dimensions in exact QMLE packing.")
  if(any(!is.finite(c(A,B,M,psi))) || any(colSums(A*A)<=0) ||
     any(colSums(B*B)<=0)) stop("Invalid or zero loading columns in QMLE start.")
  if(length(psi)==1L)psi <- rep(psi,obj$N)
  if(length(psi)!=obj$N || any(psi<=0))stop("Invalid positive noise variances in QMLE start.")
  M <- (M+t(M))/2
  initial_eigen <- eigen(M,symmetric=TRUE,only.values=TRUE)$values
  L <- tryCatch(t(chol(M)),error=function(e)NULL)
  repaired <- is.null(L)
  used_floor <- 0
  if(repaired) {
    ev <- eigen(M,symmetric=TRUE)
    used_floor <- max(repair_floor, .Machine$double.eps*max(1,max(abs(ev$values)))*10)
    M <- tcrossprod(sweep(ev$vectors,2L,sqrt(pmax(ev$values,used_floor)),"*"))
    L <- t(chol((M+t(M))/2))
  }
  mask <- lower.tri(L,diag=TRUE)
  positions <- which(mask,arr.ind=TRUE)
  lv <- L[mask]
  diagonal <- which(positions[,1L]==positions[,2L])
  lv[diagonal] <- log(diag(L))
  unbounded <- c(as.vector(cp_normalize_columns(A)),as.vector(cp_normalize_columns(B)),lv)
  if(obj$estimate_noise)unbounded <- c(unbounded,log(psi))
  par <- pmax(obj$lower,pmin(obj$upper,unbounded))
  list(par=par,covariance_repaired=repaired,covariance_repair_floor=used_floor,
       covariance_min_eigen_before=min(initial_eigen),
       clipped_parameters=sum(par!=unbounded),
       exact_before_bounds=!repaired && !any(par!=unbounded))
}

cp_grid_initial_angles <- function(candidate,truth) {
  cosine <- colSums(cp_normalize_columns(candidate)*cp_normalize_columns(truth))
  acos(pmax(-1,pmin(1,cosine)))*180/pi
}

cp_grid_qmle <- function(Y,design,cfg,seed,previous_fit=NULL) {
  p <- as.integer(design$p);q <- as.integer(design$q);d <- as.integer(design$d)
  regime <- as.character(design$regime)
  maxit_by_regime <- cp_grid_option(cfg,"maxit_by_regime",
    c(interior=300L,boundary=400L,exterior=500L))
  maxit <- as.integer(maxit_by_regime[[regime]])
  truth_starts <- as.integer(cp_grid_option(cfg,"truth_starts",6L))
  angle_deg <- cp_grid_option(cfg,"truth_angle_deg",3)
  chol_sd <- cp_grid_option(cfg,"truth_chol_sd",.05)
  noise_sd <- cp_grid_option(cfg,"truth_noise_log_sd",.05)
  estimate_noise <- isTRUE(cp_grid_option(cfg,"estimate_noise",TRUE))
  noise_bounds <- cp_grid_option(cfg,"noise_bounds",c(.05,20))
  factr <- cp_grid_option(cfg,"factr",1000)
  pgtol <- cp_grid_option(cfg,"pgtol",1e-8)
  lmm <- as.integer(cp_grid_option(cfg,"lmm",15L))
  save_all <- isTRUE(cp_grid_option(cfg,"save_all_start_fits",FALSE))
  if(length(maxit)!=1L || !is.finite(maxit) || maxit<1L)stop("Missing positive regime-specific maxit.")
  if(length(truth_starts)!=1L || truth_starts<1L || !is.finite(truth_starts))stop("truth_starts must be positive.")
  if(any(!is.finite(c(chol_sd,noise_sd))) || min(chol_sd,noise_sd)<0)stop("Invalid truth perturbation SD.")
  if(!is.finite(factr) || factr<0 || !is.finite(pgtol) || pgtol<0 || lmm<1L)
    stop("Invalid optimization tolerances.")
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_seed <- if(had_seed)get(".Random.seed",envir=.GlobalEnv)else NULL
  on.exit(if(had_seed)assign(".Random.seed",old_seed,envir=.GlobalEnv)else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(as.integer(seed))
  began <- proc.time()[["elapsed"]]
  obj <- cp_qmle_objective(Y,p,q,d,estimate_noise=estimate_noise,fixed_noise=1,
    noise_bounds=noise_bounds)
  starts <- list();source_names <- character();initial_info <- list()
  add_start <- function(A,B,M,psi,source) {
    packed <- cp_grid_pack_exact(obj,A,B,M,psi)
    i <- length(starts)+1L
    starts[[i]] <<- packed$par
    source_names[i] <<- source
    angle_A <- cp_grid_initial_angles(A,design$A)
    angle_B <- cp_grid_initial_angles(B,design$B)
    initial_info[[i]] <<- list(
      initial_A_angle_mean=mean(angle_A),initial_A_angle_max=max(angle_A),
      initial_B_angle_mean=mean(angle_B),initial_B_angle_max=max(angle_B),
      initial_M_relative_error=sqrt(sum((M-design$Mx)^2)/sum(design$Mx^2)),
      initial_noise_log_rms=sqrt(mean(log(psi)^2)),
      covariance_repaired=packed$covariance_repaired,
      covariance_repair_floor=packed$covariance_repair_floor,
      covariance_min_eigen_before=packed$covariance_min_eigen_before,
      clipped_parameters=packed$clipped_parameters,
      exact_before_bounds=packed$exact_before_bounds)
    invisible(NULL)
  }
  add_start(design$A,design$B,design$Mx,rep(1,p*q),"truth_exact")
  L0 <- t(chol((design$Mx+t(design$Mx))/2))
  if(truth_starts>1L)for(i in 2:truth_starts) {
    A <- sqrt(p)*cp_grid_perturb_directions(design$A,angle_deg)
    B <- sqrt(q)*cp_grid_perturb_directions(design$B,angle_deg)
    L <- L0
    diag(L) <- diag(L0)*exp(stats::rnorm(d,0,chol_sd))
    if(d>1L)for(row in 2:d)for(column in seq_len(row-1L))
      L[row,column] <- L0[row,column]+stats::rnorm(1L,0,chol_sd)*sqrt(design$Mx[row,row])
    psi <- exp(stats::rnorm(p*q,0,noise_sd))
    add_start(A,B,tcrossprod(L),psi,paste0("truth_near_",i-1L))
  }
  warm_index <- NULL
  if(!is.null(previous_fit)) {
    add_start(previous_fit$A,previous_fit$B,previous_fit$M,previous_fit$psi,"previous_T")
    warm_index <- length(starts)
  }
  n_starts <- length(starts)
  logs <- vector("list",n_starts)
  all_fits <- if(save_all)vector("list",n_starts)else NULL
  best <- NULL;best_value <- Inf;best_start <- NA_integer_
  for(i in seq_along(starts)) {
    tick <- proc.time()[["elapsed"]]
    initial_loss <- obj$fn(starts[[i]])
    warm_loss_delta <- if(!is.null(warm_index) && i==warm_index && length(previous_fit$nll)==1L)
      initial_loss-as.numeric(previous_fit$nll)else NA_real_
    # The warm loss changes with T. The parameters are packed exactly before
    # bounds; old and current-T objective values must not be equated.
    optimized <- tryCatch(optim(starts[[i]],obj$fn,obj$gr,method="L-BFGS-B",
      lower=obj$lower,upper=obj$upper,
      control=list(maxit=maxit,factr=factr,pgtol=pgtol,lmm=lmm)),error=function(e)
        list(par=starts[[i]],convergence=9999L,message=paste("Optimizer error; start retained:",conditionMessage(e)),
             counts=c("function"=0L,gradient=0L)))
    evaluated <- obj$evaluate(optimized$par)
    pg <- obj$projected_gradient(optimized$par,evaluated$gradient)
    counts <- optimized$counts
    row <- data.frame(start=i,source=source_names[i],nll=evaluated$value,
      initial_nll=initial_loss,convergence=optimized$convergence,
      converged=optimized$convergence==0L,
      function_evaluations=if(is.null(counts["function"]))NA_real_ else unname(counts["function"]),
      gradient_evaluations=if(is.null(counts["gradient"]))NA_real_ else unname(counts["gradient"]),
      gradient_norm=max(abs(evaluated$gradient)),projected_gradient_norm=max(abs(pg)),
      elapsed_seconds=proc.time()[["elapsed"]]-tick,
      message=if(is.null(optimized$message))""else optimized$message,
      maxit_per_start=maxit,factr=factr,pgtol=pgtol,
      warm_nll_delta_across_T=warm_loss_delta,
      stringsAsFactors=FALSE)
    logs[[i]] <- cbind(row,as.data.frame(initial_info[[i]],stringsAsFactors=FALSE))
    if(save_all) {
      parts_i <- obj$unpack(optimized$par)
      all_fits[[i]] <- list(A=parts_i$A,B=parts_i$B,M=tcrossprod(parts_i$L),L=parts_i$L,
                           psi=parts_i$psi,nll=evaluated$value,convergence=optimized$convergence)
    }
    # A caught optimizer exception must not turn an unoptimized truth start
    # into a reported estimator, even if that fallback has a finite loss.
    usable<-optimized$convergence!=9999L && is.finite(evaluated$value) && evaluated$value<1e29
    if(usable && evaluated$value<best_value) {
      best_value <- evaluated$value;best <- optimized;best_start <- i
    }
  }
  if(is.null(best))stop("Every truth-near QMLE run failed or returned an invalid objective; no estimator selected.")
  logs <- do.call(rbind,logs)
  parts <- obj$unpack(best$par)
  M <- tcrossprod(parts$L)
  Z <- cp_khatri_rao(parts$A,parts$B)
  signal_singular <- svd(Z%*%parts$L,nu=0L,nv=0L)$d
  list(A=parts$A,B=parts$B,M=M,L=parts$L,psi=parts$psi,
    nll=best_value,converged=logs$converged[best_start],
    convergence=logs$convergence[best_start],message=logs$message[best_start],
    gradient_norm=logs$gradient_norm[best_start],
    projected_gradient_norm=logs$projected_gradient_norm[best_start],
    best_start=best_start,all_starts=logs,maxit_per_start=maxit,n_starts=n_starts,
    n_truth_starts=truth_starts,n_previous_T_starts=as.integer(!is.null(previous_fit)),
    pca_method="not_used_truth_near_starts",noise_mode=if(estimate_noise)"estimated_diagonal"else"known_diagonal",
    noise_at_lower_bound=if(estimate_noise)sum(parts$psi<=noise_bounds[1L]*(1+1e-4))else 0L,
    noise_at_upper_bound=if(estimate_noise)sum(parts$psi>=noise_bounds[2L]*(1-1e-4))else 0L,
    M_min_eigen=min(eigen(M,symmetric=TRUE,only.values=TRUE)$values),
    signal_singular_values=signal_singular,
    elapsed_seconds=proc.time()[["elapsed"]]-began,
    truth_near_initialization=TRUE,
    initialization=list(truth_starts=truth_starts,truth_angle_deg=angle_deg,
      truth_chol_sd=chol_sd,truth_noise_log_sd=noise_sd,
      statement="Simulation-truth initialization; not a blind-data estimator. Noise remains fitted."),
    all_start_estimates=all_fits)
}
