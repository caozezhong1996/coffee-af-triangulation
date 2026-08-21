jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"wave3_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

exp <- extract_instruments("ukb-b-17814", p1=5e-8, clump=TRUE, r2=0.001, kb=10000)
exp$exposure <- "Caffeine drink within last hour (UKB)"
cat("Instruments:", nrow(exp), "\n")
for (oid in c("ebi-a-GCST006414","finn-b-I9_AF")) {
  oname <- ifelse(oid=="ebi-a-GCST006414","AF Nielsen 2018","AF FinnGen R7")
  cat("\n===== ", oname, "=====\n")
  out <- extract_outcome_data(snps=exp$SNP, outcomes=oid)
  out$outcome <- oname
  dat <- harmonise_data(exp, out)
  res <- mr(dat)
  res_or <- generate_odds_ratios(res)
  print(res_or[, c("method","nsnp","or","or_lci95","or_uci95","pval")])
  write.csv(res_or, file.path(outdir, paste0("wave3_", gsub("-","_",oid), ".csv")), row.names=FALSE)
  cat("-- het --\n"); print(mr_heterogeneity(dat))
  cat("-- pleio --\n"); print(mr_pleiotropy_test(dat))
}
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
