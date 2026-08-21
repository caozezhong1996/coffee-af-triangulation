jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
args <- commandArgs(trailingOnly=TRUE)
outcome_id <- args[1]; outcome_name <- args[2]
rds <- file.path(outdir, "mv_exp.rds")
if (file.exists(rds)) { mv_exp <- readRDS(rds) } else {
  mv_exp <- mv_extract_exposures(c("ukb-b-5237","ukb-b-19953","ieu-b-4877","ukb-b-5779"),
                                 pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000)
  saveRDS(mv_exp, rds)
}
log <- file(file.path(outdir, paste0("mvmr_", outcome_id, "_log.txt")),"w"); sink(log,type="output"); sink(log,type="message")
cat("===== MVMR OUTCOME:", outcome_name, "=====\n")
out_dat <- extract_outcome_data(snps = mv_exp$SNP, outcomes = outcome_id)
out_dat$outcome <- outcome_name
mvdat <- mv_harmonise_data(mv_exp, out_dat)
res <- mv_multiple(mvdat)
res_or <- generate_odds_ratios(res$result)
print(res_or[, c("exposure","nsnp","b","se","pval","or","or_lci95","or_uci95")])
write.csv(res_or, file.path(outdir, paste0("mvmr_results_", outcome_id, ".csv")), row.names=FALSE)
cat("\n-- Instrument strength (conditional F) --\n")
print(tryCatch(strength_mvmr(mvdat), error=function(e) conditionMessage(e)))
cat("\n-- Pleiotropy (Q) --\n")
print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
