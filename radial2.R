jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages({library(TwoSampleMR); library(ieugwasr); library(RadialMR)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"radial2_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")

d <- read.csv(file.path(outdir,"harmonised_finngenR11.csv")); d <- d[d$mr_keep==TRUE, ]
rw <- format_radial(BXG=d$beta.exposure, BYG=d$beta.outcome, seBXG=d$se.exposure, seBYG=d$se.outcome, RSID=d$SNP)
ivw <- ivw_radial(rw, alpha=0.05, weights=3, tol=0.0001, summary=TRUE)
cat("===== RADIAL FinnGen R11 =====\n")
print(ivw$coef); cat("Q:", ivw$qstatistic, "\n")
cat("-- outliers --\n"); print(ivw$outliers)
write.csv(ivw$outliers, file.path(outdir,"radial_outliers_FinnGenR11.csv"), row.names=FALSE)
egg <- egger_radial(rw, alpha=0.05, weights=3, summary=TRUE)
cat("-- radial egger --\n"); print(egg$coef)

cat("\n===== Wojcik 2019 (PAGE, non-UKB) instruments =====\n")
w <- extract_instruments("ebi-a-GCST008028", p1=5e-8, clump=TRUE, r2=0.001, kb=10000)
cat("instruments at 5e-8:", nrow(w), "\n")
if (nrow(w) < 3) {
  w <- extract_instruments("ebi-a-GCST008028", p1=1e-6, clump=TRUE, r2=0.001, kb=10000)
  cat("instruments at 1e-6 (fallback):", nrow(w), "\n")
}
write.csv(w, file.path(outdir,"instruments_wojcik_page.csv"), row.names=FALSE)
cat("\nDONE\n")
sink(type="message"); sink(); close(log)
