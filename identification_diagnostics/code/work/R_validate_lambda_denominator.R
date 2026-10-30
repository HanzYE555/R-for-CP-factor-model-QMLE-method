# Algebraically independent check of all denominator summaries and displayed rows.
out<-"outputs/lambda_axis_geometry_20261030"
d<-read.csv(file.path(out,"primary_by_seed.csv"),stringsAsFactors=FALSE)
L<-as.matrix(d[sprintf("rho_lambda_%02d",1:8)]);L[is.na(L)]<-0
sum2<-rowSums(L^2)
d$rho_lambda_norm<-sqrt(sum2)
# Use the closed formula here; the Python reporter sums pairwise products.
d$rho_k_norm<-sqrt(pmax(0,(sum2^2-rowSums(L^4))/2))
d$rho_H<-d$rho_min*d$rho_k_norm
d$rho_H_over_E1<-d$rho_H/d$E1_at_rho
metrics<-c("rho_min","E1_at_rho","rho_lambda_norm","rho_k_norm","rho_H","rho_H_over_E1","rho_max_loading_cosine")
by<-d[c("p","q","regime","T")]
mu<-aggregate(d[metrics],by,mean)
va<-aggregate(d[metrics],by,var)
names(mu)[match(metrics,names(mu))]<-paste0(metrics,"_mean")
names(va)[match(metrics,names(va))]<-paste0(metrics,"_variance")
direct<-merge(mu,va,by=c("p","q","regime","T"))
reported<-read.csv(file.path(out,"tables/rho_denominator_mean_variance.csv"),stringsAsFactors=FALSE)
merged<-merge(direct,reported,by=c("p","q","regime","T"),suffixes=c("_direct","_reported"))
errors<-vapply(c(paste0(metrics,"_mean"),paste0(metrics,"_variance")),function(name)
  max(abs(merged[[paste0(name,"_direct")]]-merged[[paste0(name,"_reported")]])),numeric(1L))
stopifnot(nrow(merged)==108L,max(errors)<1e-10)
small<-reported[reported$p==40&reported$q==40,]
tex<-paste(readLines("outputs/Id_Prob_QMLE_CP_Model_JASA_S23/additions/lambda/denominator_table.tex",warn=FALSE),collapse="\n")
expected<-0L
for(i in seq_len(nrow(small)))for(name in c("rho_min","E1_at_rho","rho_k_norm","rho_H_over_E1")) {
  mean_token<-sprintf(if(name=="rho_k_norm")"%.4f"else"%.6f",small[[paste0(name,"_mean")]][i])
  var_token<-sprintf("%.2e",small[[paste0(name,"_variance")]][i])
  stopifnot(grepl(mean_token,tex,fixed=TRUE),grepl(var_token,tex,fixed=TRUE))
  expected<-expected+1L
}
write.csv(data.frame(metric=names(errors),maximum_absolute_error=errors),
  file.path(out,"validation/denominator_independent_R_validation.csv"),row.names=FALSE)
jsonlite::write_json(list(status="passed",cells=nrow(merged),displayed_mean_variance_pairs=expected,
  maximum_summary_error=max(errors),method="R closed formula (sum(lambda^2)^2-sum(lambda^4))/2 vs Python pairwise products",
  paper_denominator_table_values_checked=TRUE),
  file.path(out,"validation/denominator_independent_R_validation.json"),auto_unbox=TRUE,pretty=TRUE)
print(small[c("regime","T","rho_min_mean","E1_at_rho_mean","rho_k_norm_mean","rho_H_over_E1_mean")])
cat("All108cell independent closed-form checks passed; max absolute error",max(errors),"\n")
