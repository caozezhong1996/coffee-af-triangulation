jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(MVMR)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"condf_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")
mv_exp <- readRDS(file.path(outdir,"mv_exp.rds"))

mk_r11_out <- function(){
  x <- read.delim("r11_mvmr_hits.tsv", comment.char="", check.names=FALSE); names(x)[1]<-"chrom"
  pick <- sapply(strsplit(x$rsids, ","), function(v){ m <- v[v %in% unique(mv_exp$SNP)]; if(length(m)) m[1] else NA })
  x$SNP <- pick; x <- x[!is.na(x$SNP), ]
  o <- format_data(data.frame(SNP=x$SNP, beta=x$beta, se=x$sebeta, eaf=x$af_alt,
                   effect_allele=toupper(x$alt), other_allele=toupper(x$ref), pval=x$pval),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="eaf", pval_col="pval")
  o$outcome <- "AF & flutter (FinnGen R11)"; o
}

do_cf <- function(mvdat, tag){
  snps <- rownames(mvdat$exposure_beta)
  bx <- as.data.frame(mvdat$exposure_beta); se_x <- as.data.frame(mvdat$exposure_se)
  by <- as.numeric(mvdat$outcome_beta); se_y <- as.numeric(mvdat$outcome_se)
  f <- format_mvmr(BXGs=bx, BYG=by, seBXGs=se_x, seBYG=se_y, RSID=snps)
  cat("\n== Conditional F:", tag, "==\n")
  s <- strength_mvmr(f); print(s)
  cat("-- MVMR Egger pleiotropy --\n"); print(tryCatch(pleiotropy_mvmr(f), error=function(e) conditionMessage(e)))
  row.names(s) <- NULL
  write.csv(data.frame(outcome=tag, s), file.path(outdir, paste0("condf_", tag, ".csv")), row.names=FALSE)
}

mvdat11 <- mv_harmonise_data(mv_exp, mk_r11_out()); do_cf(mvdat11, "finngenR11")
for (oid in c("ebi-a-GCST006414","finn-b-I9_AF")) {
  od <- extract_outcome_data(snps=unique(mv_exp$SNP), outcomes=oid); od$outcome <- oid
  do_cf(mv_harmonise_data(mv_exp, od), oid)
}
cat("\nCF DONE\n")
sink(type="message"); sink(); close(log)
