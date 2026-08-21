jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"mvmr_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

ids <- c("ukb-b-5237",   # Coffee intake
         "ukb-b-19953",  # BMI
         "ieu-b-4877",   # Smoking initiation (GSCAN)
         "ukb-b-5779")   # Alcohol intake frequency / drinks per week
cat("Extracting MVMR exposures...\n")
mv_exp <- mv_extract_exposures(ids, pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000)
cat("Combined SNPs:", nrow(mv_exp), "\n")
print(table(mv_exp$exposure))

run_mvmr <- function(outcome_id, outcome_name){
  cat("\n===== MVMR OUTCOME:", outcome_name, "=====\n")
  out_dat <- extract_outcome_data(snps = mv_exp$SNP, outcomes = outcome_id)
  out_dat$outcome <- outcome_name
  mvdat <- mv_harmonise_data(mv_exp, out_dat)
  res <- mv_multiple(mvdat)
  res_or <- generate_odds_ratios(res$result)
  print(res_or[, c("exposure","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  write.csv(res_or, file.path(outdir, paste0("mvmr_results_", outcome_id, ".csv")), row.names=FALSE)
  cat("\n-- Pleiotropy (MVMR Egger) --\n")
  print(tryCatch(mv_residual(mvdat), error=function(e) conditionMessage(e)))
  cat("\n-- Instrument strength --\n")
  print(tryCatch({s <- strength_mvmr(mvdat); s}, error=function(e) conditionMessage(e)))
  cat("\n-- Pleiotropy test --\n")
  print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
}

run_mvmr("ebi-a-GCST006414", "AF (Nielsen 2018)")
run_mvmr("finn-b-I9_AF", "AF & flutter (FinnGen R7)")
cat("\nMVMR DONE\n")
sink(type="message"); sink(); close(log)
