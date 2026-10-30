options(lambda.runner.no_autorun=TRUE)
source("work/R_lambda_geometry_runner.R",encoding="UTF-8")
v_out <- file.path(LAMBDA_OUT,"iteration_sensitivity")
v_files <- sort(list.files(file.path(v_out,"cases"),pattern="\\.rds$",full.names=TRUE))
stopifnot(length(v_files)==31L)
v_rows <- list();v_pairs <- list()
for(v_file in v_files) {
  v_case <- readRDS(v_file);v_dg<-v_case$diag;v_pr<-lg_problem(v_dg$core,.90)
  v_check <- function(lambda,E1,rho,kind,order) {
    M<-v_pr$A%*%diag(lambda,nrow=v_pr$d)%*%t(v_pr$B)
    normM<-sqrt(sum(M*M));sv<-svd(M,nu=0L,nv=0L)$d
    direct_E<-sqrt(sum(sv[-1L]^2)/sum(sv^2))
    kv<-lambda[v_pr$pairs[1L,]]*lambda[v_pr$pairs[2L,]]
    direct_rho<-sqrt(sum((v_pr$C%*%kv)^2)/sum(kv^2))
    axis<-max(abs(lambda))/sqrt(sum(lambda^2))
    physical<-max(abs(drop(crossprod(v_pr$Z,as.vector(M)))))/normM
    stopifnot(abs(normM-1)<1e-5,axis<=.90+2e-6,abs(direct_E-E1)<2e-6,abs(direct_rho-rho)<2e-6)
    data.frame(case=basename(v_file),kind=kind,order=order,
      unit_Frobenius_error=abs(normM-1),E1_direct_error=abs(direct_E-E1),
      rho_direct_error=abs(direct_rho-rho),max_lambda_cosine=axis,
      max_physical_cosine=physical,stringsAsFactors=FALSE)
  }
  for(v_kind in c("rho","E")) {
    v_candidates<-if(v_kind=="rho")v_dg$rho_local_minima else v_dg$E_local_minima
    for(v_i in seq_along(v_candidates)) {
      v_z<-v_candidates[[v_i]]
      v_rows[[length(v_rows)+1L]]<-v_check(v_z$lambda,v_z$E1,v_z$rho,v_kind,v_i)
    }
  }
  for(v_i in seq_len(nrow(v_case$cap_rows))) {
    v_z<-v_case$cap_rows[v_i,]
    for(v_kind in c("center","cap5","cap10")) {
      v_lam<-unlist(v_z[paste0(v_kind,sprintf("_lambda_%02d",seq_len(v_pr$d)))],use.names=FALSE)
      v_E<-v_z[[switch(v_kind,center="E1_center",cap5="E1_cap5_saved",cap10="E1_cap10")]]
      v_rho<-v_z[[switch(v_kind,center="rho_center",cap5="rho_cap5_saved",cap10="rho_cap10")]]
      v_rows[[length(v_rows)+1L]]<-v_check(v_lam,v_E,v_rho,v_kind,v_i)
    }
    stopifnot(v_z$E1_cap10<=v_z$E1_cap5_saved+2e-6,v_z$E1_cap5_saved<=v_z$E1_center+2e-6,
              v_z$center_to_cap5_deg<=5+1e-4,v_z$center_to_cap10_deg<=10+1e-4)
  }
  v_angles<-vapply(v_dg$E_local_minima,function(z)acos(min(1,abs(sum(v_dg$lambda_rho*z$lambda))/sqrt(sum(v_dg$lambda_rho^2)*sum(z$lambda^2))))*180/pi,numeric(1L))
  stopifnot(abs(min(v_angles)-v_case$primary$center_to_E1_nearest_deg_090)<1e-4)
  v_all<-list()
  for(v_r in v_dg$rho_local_minima)for(v_e in v_dg$E_local_minima) {
    v_angle<-acos(min(1,abs(sum(v_r$lambda*v_e$lambda))/sqrt(sum(v_r$lambda^2)*sum(v_e$lambda^2))))*180/pi
    if(v_angle<=5+1e-8)v_all[[length(v_all)+1L]]<-max(v_r$rho,v_r$E1,v_e$rho,v_e$E1)
  }
  v_Q<-if(length(v_all))min(unlist(v_all))else Inf
  v_Qerror<-if(is.finite(v_Q)&&is.finite(v_case$primary$pair_joint_max))abs(v_Q-v_case$primary$pair_joint_max)else if(identical(is.infinite(v_Q),is.infinite(v_case$primary$pair_joint_max)))0 else Inf
  stopifnot(v_Qerror<1e-10)
  v_pairs[[length(v_pairs)+1L]]<-data.frame(case=basename(v_file),
    Q_direct_error=v_Qerror,Delta_direct_error=abs(min(v_angles)-v_case$primary$center_to_E1_nearest_deg_090),
    source_unchanged=identical(v_case$source_md5,unname(tools::md5sum(v_case$source_path))))
}
write.csv(do.call(rbind,v_rows),file.path(v_out,"lambda_iteration_direct_point_validation.csv"),row.names=FALSE)
write.csv(do.call(rbind,v_pairs),file.path(v_out,"lambda_iteration_direct_case_validation.csv"),row.names=FALSE)
cat(sprintf("Validated %d direct candidate/cap points in %d continued-fit cases.\n",length(v_rows),length(v_files)))
