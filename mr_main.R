jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
dir.create(outdir, showWarnings=FALSE)
log <- file(file.path(outdir, "run_log.txt"), "w")
sink(log, type="output"); sink(log, type="message")

# 1. 暴露工具变量：Coffee intake (UKB, n=428,860)
exp_dat <- extract_instruments(outcomes="ukb-b-5237", p1=5e-8, clump=TRUE, r2=0.001, kb=10000)
exp_dat$exposure <- "Coffee intake (cups/day, UKB)"
cat("Instruments:", nrow(exp_dat), "\n")
print(exp_dat[, c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure")])
write.csv(exp_dat, file.path(outdir,"instruments.csv"), row.names=FALSE)

run_mr <- function(outcome_id, outcome_name){
  cat("\n===== OUTCOME:", outcome_name, "=====\n")
  out_dat <- extract_outcome_data(snps=exp_dat$SNP, outcomes=outcome_id)
  if(is.null(out_dat)){ cat("NO OUTCOME DATA\n"); return(NULL) }
  out_dat$outcome <- outcome_name
  cat("SNPs in outcome:", nrow(out_dat), "\n")
  dat <- harmonise_data(exp_dat, out_dat)
  write.csv(dat, file.path(outdir, paste0("harmonised_", outcome_id, ".csv")), row.names=FALSE)
  res <- mr(dat)
  res_or <- generate_odds_ratios(res)
  print(res_or[, c("outcome","method","nsnp","or","or_lci95","or_uci95","pval")])
  write.csv(res_or, file.path(outdir, paste0("results_", outcome_id, ".csv")), row.names=FALSE)
  cat("\n-- Heterogeneity --\n"); print(mr_heterogeneity(dat))
  cat("\n-- Pleiotropy (Egger intercept) --\n"); print(tryCatch(mr_pleiotropy_test(dat), error=function(e) conditionMessage(e)))
  loo <- mr_leaveoneout(dat); write.csv(loo, file.path(outdir, paste0("loo_", outcome_id, ".csv")), row.names=FALSE)
  singlesnp <- mr_singlesnp(dat); write.csv(singlesnp, file.path(outdir, paste0("singlesnp_", outcome_id, ".csv")), row.names=FALSE)
  return(list(dat=dat, res=res))
}

r1 <- run_mr("ebi-a-GCST006414", "AF (Nielsen 2018, AFGen+)")
r2 <- run_mr("finn-b-I9_AF", "AF & flutter (FinnGen R7)")
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
