jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(MVMR)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"mvmr_r11_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

mv_exp <- readRDS(file.path(outdir,"mv_exp.rds"))

read_r11 <- function(f){
  x <- read.delim(f, comment.char="", check.names=FALSE)
  names(x)[1] <- "chrom"
  x
}
mk_outcome <- function(x, ids, name){
  pick <- sapply(strsplit(x$rsids, ","), function(v){ m <- v[v %in% ids]; if(length(m)) m[1] else NA })
  x$SNP <- pick; x <- x[!is.na(x$SNP), ]
  o <- format_data(data.frame(SNP=x$SNP, beta=x$beta, se=x$sebeta, eaf=x$af_alt,
                              effect_allele=toupper(x$alt), other_allele=toupper(x$ref), pval=x$pval),
                   type="outcome", snp_col="SNP", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele", eaf_col="eaf", pval_col="pval")
  o$outcome <- name
  o
}

run_mvmr_local <- function(out_dat, tag){
  cat("\n===== MVMR OUTCOME:", tag, "=====\n")
  mvdat <- mv_harmonise_data(mv_exp, out_dat)
  res <- mv_multiple(mvdat)
  res_or <- generate_odds_ratios(res$result)
  print(res_or[, c("exposure","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  write.csv(res_or, file.path(outdir, paste0("mvmr_results_", tag, ".csv")), row.names=FALSE)
  cat("\n-- Conditional F (strength_mvmr) --\n")
  print(tryCatch(strength_mvmr(mvdat), error=function(e) conditionMessage(e)))
  cat("\n-- MVMR Egger pleiotropy --\n")
  print(tryCatch(pleiotropy_mvmr(mvdat), error=function(e) conditionMessage(e)))
}

# R11 outcome (local)
r11m <- read_r11("r11_mvmr_hits.tsv")
out_r11 <- mk_outcome(r11m, unique(mv_exp$SNP), "AF & flutter (FinnGen R11)")
cat("R11 outcome SNPs:", nrow(out_r11), "\n")
run_mvmr_local(out_r11, "finngenR11")

# Nielsen & R7 outcomes (API) for conditional F
for (oid in c("ebi-a-GCST006414","finn-b-I9_AF")) {
  cat("\n== re-extract outcome for conditional F:", oid, "==\n")
  od <- extract_outcome_data(snps=unique(mv_exp$SNP), outcomes=oid)
  od$outcome <- oid
  run_mvmr_local(od, oid)
}

# Subtype MR vs R11
subs <- list("ebi-a-GCST90096908"="Decaffeinated", "ebi-a-GCST90096913"="Ground", "ebi-a-GCST90096914"="Instant")
r11s <- read_r11("r11_subtype_hits.tsv")
for (id in names(subs)) {
  cat("\n===== SUBTYPE:", subs[[id]], "x FinnGen R11 =====\n")
  e <- read.csv(paste0(outdir,"/subtype_inst_", id, ".csv"))
  e <- format_data(e, type="exposure", snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure",
                   effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure",
                   eaf_col="eaf.exposure", pval_col="pval.exposure")
  e$exposure <- paste0(subs[[id]], " coffee")
  o <- mk_outcome(r11s, e$SNP, "AF & flutter (FinnGen R11)")
  d <- harmonise_data(e, o, action=2)
  d <- d[d$mr_keep==TRUE, ]
  cat("kept:", nrow(d), "\n")
  if (nrow(d) >= 2) {
    res <- mr(d, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))
    print(res[, c("method","nsnp","b","se","pval")])
    write.csv(res, file.path(outdir, paste0("r11_subtype_", id, ".csv")), row.names=FALSE)
  }
}
cat("\nALL DONE\n")
sink(type="message"); sink(); close(log)
