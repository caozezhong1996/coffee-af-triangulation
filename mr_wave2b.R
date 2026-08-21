jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"wave2b_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")
cat("##### REVERSE MR: AF -> Coffee intake #####\n")
exp_af <- extract_instruments("ebi-a-GCST006414", p1=5e-8, clump=TRUE, r2=0.001, kb=10000)
exp_af$exposure <- "Atrial fibrillation (genetic liability)"
cat("AF instruments:", nrow(exp_af), "\n")
out_c <- extract_outcome_data(snps=exp_af$SNP, outcomes="ukb-b-5237")
out_c$outcome <- "Coffee intake (cups/day, UKB)"
dat2 <- harmonise_data(exp_af, out_c)
res2 <- mr(dat2, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
print(res2[, c("method","nsnp","b","se","pval")])
write.csv(res2, file.path(outdir,"wave2_reverse_AF_to_coffee.csv"), row.names=FALSE)
cat("\n-- heterogeneity --\n"); print(mr_heterogeneity(dat2))
cat("\n-- pleiotropy --\n"); print(mr_pleiotropy_test(dat2))
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
