suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"r11_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

inst <- read.csv("mr_results/instruments.csv")
cat("instruments.csv cols:", paste(names(inst), collapse=","), "\n")
r11 <- read.delim("r11_af_hits.tsv", comment.char="", check.names=FALSE)
names(r11)[1] <- "chrom"
cat("r11 rows:", nrow(r11), "\n")

# resolve which rsid per row matches instruments
pick <- sapply(strsplit(r11$rsids, ","), function(x){ m <- x[x %in% inst$SNP]; if(length(m)) m[1] else NA })
r11$SNP <- pick
r11 <- r11[!is.na(r11$SNP), ]
cat("r11 rows with resolved rsid:", nrow(r11), "\n")

exp <- format_data(inst, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                   effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                   eaf_col="eaf.exposure", pval_col="pval.exposure")
exp$exposure <- "Coffee intake (cups/day, UKB)"
out <- format_data(data.frame(SNP=r11$SNP, beta=r11$beta, se=r11$sebeta, eaf=r11$af_alt,
                              effect_allele=toupper(r11$alt), other_allele=toupper(r11$ref), pval=r11$pval),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", pval_col="pval")
out$outcome <- "AF & flutter (FinnGen R11)"
dat <- harmonise_data(exp, out, action=2)
cat("harmonised:", nrow(dat), " kept:", sum(dat$mr_keep), "\n")
write.csv(dat, file.path(outdir,"harmonised_finngenR11.csv"), row.names=FALSE)

d <- dat[dat$mr_keep==TRUE, ]
res <- mr(d, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
print(res[, c("method","nsnp","b","se","pval")])
write.csv(res, file.path(outdir,"results_finngenR11.csv"), row.names=FALSE)
cat("\n-- heterogeneity --\n"); print(mr_heterogeneity(d))
cat("\n-- pleiotropy --\n"); print(mr_pleiotropy_test(d))
cat("\n-- PRESSO --\n")
pr <- mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
                SdOutcome="se.outcome", SdExposure="se.exposure",
                OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=d,
                NbDistribution=2000, SignifThreshold=0.05)
print(pr$`Main MR results`)
ot <- pr$`MR-PRESSO results`$`Outlier Test`
if (!is.null(ot)) { ot$SNP <- d$SNP; ot <- ot[order(ot$Pvalue),]; print(head(ot[ot$Pvalue<0.05,],10))
  write.csv(ot, file.path(outdir,"p0_presso_outliers_finngenR11.csv"), row.names=FALSE) }
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
