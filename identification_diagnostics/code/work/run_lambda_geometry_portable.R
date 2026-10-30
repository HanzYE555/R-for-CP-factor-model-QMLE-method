# Portable entry point around the unchanged, audited geometry runner.
# Run from the extracted package's root directory.
options(lambda.runner.no_autorun=TRUE)
source("work/R_lambda_geometry_runner.R",encoding="UTF-8")
portable_args<-commandArgs(trailingOnly=TRUE)
source_at<-match("--source",portable_args)
if(is.na(source_at)||source_at==length(portable_args))
  stop("Provide --source followed by the CP_QMLE_grid_results directory containing checkpoints/.")
LAMBDA_SOURCE<-normalizePath(portable_args[source_at+1L],winslash="/",mustWork=TRUE)
stopifnot(dir.exists(file.path(LAMBDA_SOURCE,"checkpoints")))
output_at<-match("--output",portable_args)
if(!is.na(output_at)) {
  if(output_at==length(portable_args))stop("Missing --output directory.")
  LAMBDA_OUT<-normalizePath(portable_args[output_at+1L],winslash="/",mustWork=FALSE)
}
dir.create(file.path(LAMBDA_OUT,"cases"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(LAMBDA_OUT,"parts"),recursive=TRUE,showWarnings=FALSE)
lg_main()
