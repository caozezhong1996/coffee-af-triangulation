suppressMessages(library(TwoSampleMR))
outdir <- "./mr_results"
exp_w <- read.csv(file.path(outdir,"instruments_wojcik_page.csv"))
exp_w <- format_data(exp_w, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                     effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                     eaf_col="eaf.exposure", pval_col="pval.exposure")
exp_w$exposure <- "Coffee intake (cups/day, PAGE/Wojcik 2019)"
x <- read.delim("r11_wojcik_hits.tsv", comment.char="", check.names=FALSE); names(x)[1]<-"chrom"
pick <- sapply(strsplit(x$rsids, ","), function(v){ m <- v[v %in% exp_w$SNP]; if(length(m)) m[1] else NA })
x$SNP <- pick; x <- x[!is.na(x$SNP), ]
o <- format_data(data.frame(SNP=x$SNP, beta=x$beta, se=x$sebeta, eaf=x$af_alt,
                 effect_allele=toupper(x$alt), other_allele=toupper(x$ref), pval=x$pval),
                 type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                 effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="eaf", pval_col="pval")
o$outcome <- "AF & flutter (FinnGen R11)"
d <- harmonise_data(exp_w, o, action=2); d <- d[d$mr_keep==TRUE, ]
cat("kept:", nrow(d), "\n")
r <- mr(d, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
print(r[, c("method","nsnp","b","se","pval")])
write.csv(r, file.path(outdir,"wojcik_finnR11.csv"), row.names=FALSE)
cat("DONE\n")
