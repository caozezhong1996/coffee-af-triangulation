jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"wave2_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

# ---- 1. 咖啡类型 MR：instant / ground / decaf -> AF ----
types <- c("ebi-a-GCST90096914"="Instant coffee", "ebi-a-GCST90096913"="Ground coffee", "ebi-a-GCST90096908"="Decaffeinated coffee")
outcomes <- c("ebi-a-GCST006414"="AF Nielsen 2018", "finn-b-I9_AF"="AF FinnGen R7")
for (tid in names(types)) {
  cat("\n##### EXPOSURE:", types[tid], "#####\n")
  exp <- tryCatch(extract_instruments(tid, p1=5e-8, clump=TRUE, r2=0.001, kb=10000),
                  error=function(e){cat("ERR:",conditionMessage(e),"\n"); NULL})
  if (is.null(exp) || nrow(exp)==0) { cat("No instruments\n"); next }
  exp$exposure <- paste0(types[tid], " consumption")
  cat("Instruments:", nrow(exp), "\n")
  for (oid in names(outcomes)) {
    out <- tryCatch(extract_outcome_data(snps=exp$SNP, outcomes=oid), error=function(e) NULL)
    if (is.null(out)) { cat("  no outcome data for", oid, "\n"); next }
    out$outcome <- outcomes[oid]
    dat <- harmonise_data(exp, out)
    res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
    res_or <- generate_odds_ratios(res)
    print(res_or[, c("exposure","outcome","method","nsnp","or","or_lci95","or_uci95","pval")])
    write.csv(res_or, file.path(outdir, paste0("wave2_", gsub("-","_",tid), "_", gsub("-","_",oid), ".csv")), row.names=FALSE)
    print(mr_pleiotropy_test(dat))
  }
}

# ---- 2. 反向 MR：AF -> 咖啡摄入 ----
cat("\n##### REVERSE MR: AF -> Coffee intake #####\n")
exp_af <- tryCatch(extract_instruments("ebi-a-GCST006414", p1=5e-8, clump=TRUE, r2=0.001, kb=10000),
                   error=function(e){cat("ERR:",conditionMessage(e),"\n"); NULL})
if (!is.null(exp_af)) {
  exp_af$exposure <- "Atrial fibrillation (genetic liability)"
  cat("AF instruments:", nrow(exp_af), "\n")
  out_c <- extract_outcome_data(snps=exp_af$SNP, outcomes="ukb-b-5237")
  out_c$outcome <- "Coffee intake (cups/day, UKB)"
  dat2 <- harmonise_data(exp_af, out_c)
  res2 <- mr(dat2, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  print(res2[, c("method","nsnp","b","se","pval")])
  write.csv(res2, file.path(outdir,"wave2_reverse_AF_to_coffee.csv"), row.names=FALSE)
  print(mr_heterogeneity(dat2))
  print(mr_pleiotropy_test(dat2))
}
cat("\nWAVE2 DONE\n")
sink(type="message"); sink(); close(log)
