jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"wojcik_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")
exp_w <- read.csv(file.path(outdir,"instruments_wojcik_page.csv"))
cat("Wojcik instruments:", nrow(exp_w), " mean F:", round(mean((exp_w$beta.exposure/exp_w$se.exposure)^2),1), "\n")
exp_w$exposure <- "Coffee intake (cups/day, PAGE/Wojcik 2019)"
# Nielsen outcome via API
out_n <- extract_outcome_data(snps=exp_w$SNP, outcomes="ebi-a-GCST006414")
out_n$outcome <- "AF (Nielsen 2018)"
d1 <- harmonise_data(exp_w, out_n); d1 <- d1[d1$mr_keep==TRUE, ]
cat("Nielsen kept:", nrow(d1), "\n")
if (nrow(d1) >= 3) {
  r1 <- mr(d1, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
  print(r1[, c("method","nsnp","b","se","pval")])
  write.csv(r1, file.path(outdir,"wojcik_nielsen.csv"), row.names=FALSE)
}
# R11 outcome locally
x <- read.delim("finngen_R11_I9_AF_w.tsv", comment.char="", check.names=FALSE) if (file.exists("finngen_R11_I9_AF_w.tsv")) else NULL
cat("R11 hits file present:", !is.null(x), "\n")
sink(type="message"); sink(); close(log)
