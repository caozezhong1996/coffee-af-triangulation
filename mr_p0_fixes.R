jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(MRPRESSO)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"p0_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

cat("##### 1. DECAF x FinnGen R7 #####\n")
exp_decaf <- extract_instruments("ebi-a-GCST90096908", p1=5e-8, clump=TRUE, r2=0.001, kb=10000)
exp_decaf$exposure <- "Decaffeinated coffee consumption"
cat("decaf instruments:", nrow(exp_decaf), "\n")
out_fg <- extract_outcome_data(snps=exp_decaf$SNP, outcomes="finn-b-I9_AF")
out_fg$outcome <- "AF & flutter (FinnGen R7)"
dat_d <- harmonise_data(exp_decaf, out_fg)
res_d <- mr(dat_d, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
print(res_d[, c("method","nsnp","b","se","pval")])
write.csv(res_d, file.path(outdir,"p0_decaf_finn.csv"), row.names=FALSE)

cat("\n##### 2. MR-PRESSO outlier SNP identities (main analysis) #####\n")
for (ds in c("ebi-a-GCST006414","finn-b-I9_AF")) {
  f <- file.path(outdir, paste0("harmonised_", gsub("-","-",ds), ".csv"))
  cat("\n--- PRESSO outlier test:", ds, "---\n")
  d <- read.csv(file.path(outdir, paste0("harmonised_", ds, ".csv")))
  d <- d[d$mr_keep==TRUE, ]
  cat("SNPs kept:", nrow(d), "\n")
  pr <- mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
                  SdOutcome="se.outcome", SdExposure="se.exposure",
                  OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=d,
                  NbDistribution=2000, SignifThreshold=0.05)
  print(pr$`Main MR results`)
  ot <- pr$`MR-PRESSO results`$`Outlier Test`
  if (!is.null(ot)) {
    ot <- ot[order(ot$Pvalue), ]
    print(head(ot[ot$Pvalue<0.05, ], 20))
    write.csv(ot, file.path(outdir, paste0("p0_presso_outliers_", ds, ".csv")), row.names=FALSE)
  }
}
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
