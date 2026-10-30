# Read-only provenance for the compact lambda-geometry data package.
out<-"outputs/lambda_axis_geometry_20261030"
paths<-sort(list.files(file.path(out,"cases"),pattern="[.]rds$",full.names=TRUE))
stopifnot(length(paths)==10796L)
rows<-lapply(paths,function(path) {
  x<-readRDS(path);stopifnot(identical(x$status,"ok"))
  data.frame(p=x$job$p,q=x$job$q,regime=x$job$regime,seed=x$job$seed,T=x$job$T,
    source_path=x$source_path,source_md5=x$source_md5,source_bytes=x$source_bytes,
    geometry_signature=x$signature,stringsAsFactors=FALSE)
})
tab<-do.call(rbind,rows)
stopifnot(!anyDuplicated(tab[c("p","q","regime","seed","T")]))
tab$source_md5_verified<-unname(tools::md5sum(tab$source_path))==tab$source_md5
stopifnot(all(tab$source_md5_verified))
tab$relative_checkpoint_path<-sub("^.*?/checkpoints/","checkpoints/",gsub("\\\\","/",tab$source_path))
dir.create(file.path(out,"validation"),showWarnings=FALSE)
write.csv(tab,file.path(out,"validation/source_checkpoint_provenance.csv"),row.names=FALSE)
jsonlite::write_json(list(status="passed",source_fits=nrow(tab),
  total_source_bytes=sum(tab$source_bytes),source_hashes_independently_verified=all(tab$source_md5_verified),
  source_files_copied_or_modified=FALSE,new_geometry_RDS_files_in_compact_package=FALSE),
  file.path(out,"validation/source_checkpoint_provenance.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Verified source checkpoint hashes:",nrow(tab),"source bytes",sum(tab$source_bytes),"\n")
