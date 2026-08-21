jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr)})
outdir <- "./mr_results"
mv_exp <- mv_extract_exposures(c("ukb-b-5237","ukb-b-19953","ieu-b-4877","ukb-b-5779"),
                               pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000)
saveRDS(mv_exp, file.path(outdir, "mv_exp.rds"))
cat("SAVED", nrow(mv_exp), "SNPs\n")
