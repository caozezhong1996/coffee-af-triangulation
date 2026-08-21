suppressMessages({library(TwoSampleMR); library(RadialMR)})
outdir <- "./mr_results"
log <- file(file.path(outdir,"radial_log.txt"),"w"); sink(log,type="output"); sink(log,type="message")
run_radial <- function(csvfile, tag){
  d <- read.csv(csvfile); d <- d[d$mr_keep==TRUE, ]
  rw <- format_radial(BXG=d$beta.exposure, BYG=d$beta.outcome, seBXG=d$se.exposure, seBYG=d$se.outcome, RSID=d$SNP)
  cat("\n===== RADIAL:", tag, "=====\n")
  ivw <- ivw_radial(rw, alpha=0.05, weights=3, tol=0.0001, summary=TRUE)
  print(ivw)
  cat("-- Radial outliers --\n")
  if (!is.null(ivw$outliers) && nrow(ivw$outliers)>0) print(ivw$outliers) else cat("none flagged\n")
  egg <- tryCatch(egger_radial(rw, alpha=0.05, weights=3, summary=TRUE), error=function(e) conditionMessage(e))
  print(egg)
  res <- data.frame(dataset=tag,
    ivw_coef=ivw$coef["Estimate",1], ivw_se=ivw$coef["Estimate",2], ivw_p=ivw$coef["Estimate",4],
    ivw_qstat=ivw$qstatistic[1], ivw_qp=ivw$qstatistic[2])
  write.csv(res, file.path(outdir, paste0("radial_", tag, ".csv")), row.names=FALSE)
  pdf(file.path(outdir, paste0("radial_plot_", tag, ".pdf")), width=7, height=7)
  print(plot_radial(rw, radial_est=FALSE, show_outliers=TRUE, scale=FALSE))
  dev.off()
}
run_radial(file.path(outdir,"harmonised_ebi-a-GCST006414.csv"), "Nielsen2018")
run_radial(file.path(outdir,"harmonised_finngenR11.csv"), "FinnGenR11")
cat("\nRADIAL DONE\n")
sink(type="message"); sink(); close(log)
