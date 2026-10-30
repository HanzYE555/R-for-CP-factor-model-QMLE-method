# Re-search saved continued QMLE fits on the new coefficient-angle domain.
# No likelihood optimization is performed; source checkpoints are read-only.
options(lambda.runner.no_autorun=TRUE)
source("work/R_lambda_geometry_runner.R",encoding="UTF-8")

iteration_source <- "outputs/final_summary_20261030/iteration_cap10_audit/source_records"
iteration_out <- file.path(LAMBDA_OUT,"iteration_sensitivity")
dir.create(file.path(iteration_out,"cases"),recursive=TRUE,showWarnings=FALSE)
iteration_files <- sort(list.files(iteration_source,pattern="\\.rds$",recursive=TRUE,full.names=TRUE))
stopifnot(length(iteration_files)==31L)
iteration_manifest <- read.csv(file.path(LAMBDA_OUT,"manifest.csv"),stringsAsFactors=FALSE)
iteration_metrics <- c("rho_min","E1_at_rho","E1_min","rho_at_E1","E1_cap_090","E1_cap10",
  "rho_cap10","center_to_E1_nearest_deg_090","pair_joint_max","center_to_cap10_deg",
  "cap10_max_loading_cosine","cap10_max_lambda_cosine","E1_local_within10_found",
  "local_best_angle10","J_min","J_relative_min")
iteration_flags <- c("qmle_converged","pair_exists","cap10_angle_edge","cap10_axis_edge",
                     "local_candidate_exists10","rho_best_success","E_best_success")
iteration_rows <- list();iteration_primary <- list();iteration_old_primary <- list()
for(iteration_i in seq_along(iteration_files)) {
  iteration_path <- iteration_files[iteration_i]
  iteration_record <- readRDS(iteration_path)
  stopifnot(identical(iteration_record$status,"ok"))
  iteration_stage <- if(grepl("/adaptive_final/",gsub("\\\\","/",iteration_path)))"adaptive_final"else"common_stage2"
  iteration_job <- iteration_manifest[iteration_manifest$p==iteration_record$p &
    iteration_manifest$q==iteration_record$q & iteration_manifest$regime==iteration_record$regime &
    iteration_manifest$seed==iteration_record$seed & iteration_manifest$T==iteration_record$T,,drop=FALSE]
  stopifnot(nrow(iteration_job)==1L)
  iteration_job <- as.list(iteration_job[1L,])
  iteration_target <- file.path(iteration_out,"cases",sprintf("%s_p%03d_q%03d_%s_s%03d_T%d.rds",
    iteration_stage,iteration_job$p,iteration_job$q,iteration_job$regime,iteration_job$seed,iteration_job$T))
  iteration_baseline <- readRDS(lg_target(iteration_job))
  stopifnot(identical(iteration_baseline$status,"ok"),identical(iteration_baseline$signature,LAMBDA_SIGNATURE))
  iteration_original <- readRDS(iteration_baseline$source_path)
  iteration_md5 <- unname(tools::md5sum(iteration_path))
  iteration_core <- lg_record_diag(iteration_record)$core
  iteration_matching <- iteration_core$matching
  iteration_A <- cp_unit_columns(sweep(iteration_record$fit$A[,iteration_matching$estimated_for_truth,drop=FALSE],2L,iteration_matching$signs_A,"*"))
  iteration_B <- cp_unit_columns(sweep(iteration_record$fit$B[,iteration_matching$estimated_for_truth,drop=FALSE],2L,iteration_matching$signs_B,"*"))
  iteration_gram_error <- max(abs(crossprod(iteration_A)-crossprod(iteration_core$A)),
                             abs(crossprod(iteration_B)-crossprod(iteration_core$B)))
  stopifnot(iteration_gram_error<1e-10)
  iteration_started <- proc.time()[[3L]]
  if(file.exists(iteration_target)) {
    iteration_result <- readRDS(iteration_target)
    stopifnot(identical(iteration_result$signature,LAMBDA_SIGNATURE),identical(iteration_result$source_md5,iteration_md5))
  } else {
    iteration_result <- lg_search(iteration_record,as.integer(20261101+iteration_job$case_index*101L),LAMBDA_CFG)
    iteration_summary <- lg_make_primary(iteration_job,iteration_result,iteration_record)
    iteration_result$primary <- iteration_summary$primary
    iteration_result$pairs <- iteration_summary$pairs
    iteration_result$signature <- LAMBDA_SIGNATURE
    iteration_result$source_md5 <- iteration_md5
    iteration_result$source_path <- iteration_path
    iteration_result$baseline_path <- lg_target(iteration_job)
    iteration_result$seconds <- proc.time()[[3L]]-iteration_started
    iteration_result$status <- "ok"
    saveRDS(iteration_result,iteration_target,compress=TRUE)
  }
  stopifnot(identical(iteration_md5,unname(tools::md5sum(iteration_path))))
  iteration_old <- iteration_baseline$primary
  iteration_new <- iteration_result$primary
  iteration_old$J_min <- iteration_baseline$diag$J_min
  iteration_old$J_relative_min <- iteration_baseline$diag$J_relative_min
  iteration_new$J_min <- iteration_result$diag$J_min
  iteration_new$J_relative_min <- iteration_result$diag$J_relative_min
  iteration_row <- cbind(lg_keyframe(iteration_job),data.frame(
    comparison_group=iteration_stage,case_index=iteration_job$case_index,
    core_fit_gram_max_error=iteration_gram_error,source_md5=iteration_md5,
    old_nll=iteration_original$fit$nll,new_nll=iteration_record$fit$nll,
    old_projected_gradient=iteration_original$fit$projected_gradient_norm,
    new_projected_gradient=iteration_record$fit$projected_gradient_norm,
    seconds=iteration_result$seconds,stringsAsFactors=FALSE))
  for(iteration_metric in iteration_metrics) {
    iteration_a <- iteration_old[[iteration_metric]]
    iteration_b <- iteration_new[[iteration_metric]]
    stopifnot(length(iteration_a)==1L,length(iteration_b)==1L)
    iteration_row[[paste0("old_",iteration_metric)]] <- iteration_a
    iteration_row[[paste0("new_",iteration_metric)]] <- iteration_b
    iteration_row[[paste0("change_",iteration_metric)]] <- if(is.finite(iteration_a)&&is.finite(iteration_b))iteration_b-iteration_a else NA_real_
  }
  for(iteration_flag in iteration_flags) {
    iteration_row[[paste0("old_",iteration_flag)]] <- iteration_old[[iteration_flag]]
    iteration_row[[paste0("new_",iteration_flag)]] <- iteration_new[[iteration_flag]]
  }
  iteration_old$comparison_group <- iteration_stage
  iteration_new$comparison_group <- iteration_stage
  iteration_rows[[iteration_i]] <- iteration_row
  iteration_old_primary[[iteration_i]] <- iteration_old
  iteration_primary[[iteration_i]] <- iteration_new
  cat(sprintf("LAMBDA continuation %d/%d: %s p%dq%d %s seed%d T%d (%.2fs)\n",
    iteration_i,length(iteration_files),iteration_stage,iteration_job$p,iteration_job$q,
    iteration_job$regime,iteration_job$seed,iteration_job$T,iteration_result$seconds));flush.console()
}
write.csv(cp_grid_bind_frames(iteration_rows),file.path(iteration_out,"paired_lambda_iteration_by_seed.csv"),row.names=FALSE)
write.csv(cp_grid_bind_frames(iteration_old_primary),file.path(iteration_out,"original_lambda_primary.csv"),row.names=FALSE)
write.csv(cp_grid_bind_frames(iteration_primary),file.path(iteration_out,"continued_lambda_primary.csv"),row.names=FALSE)
cat("Saved new-domain paired continuation geometry; no QMLE performed.\n")
