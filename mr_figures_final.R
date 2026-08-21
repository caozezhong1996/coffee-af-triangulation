suppressMessages({library(TwoSampleMR); library(ggplot2); library(patchwork)})
outdir <- "./mr_results"
ids  <- c("ebi-a-GCST006414", "finn-b-I9_AF")
labs <- c("A", "B")

load_dat <- function(oid){
  dat <- read.csv(file.path(outdir, paste0("harmonised_", oid, ".csv")), stringsAsFactors=FALSE)
  dat <- dat[dat$mr_keep==TRUE, ]
  exp <- data.frame(SNP=dat$SNP, beta.exposure=dat$beta.exposure, se.exposure=dat$se.exposure,
                    effect_allele.exposure=dat$effect_allele.exposure, other_allele.exposure=dat$other_allele.exposure,
                    eaf.exposure=dat$eaf.exposure, exposure=dat$exposure[1], id.exposure=dat$id.exposure[1],
                    pval.exposure=dat$pval.exposure, mr_keep.exposure=TRUE, stringsAsFactors=FALSE)
  out <- data.frame(SNP=dat$SNP, beta.outcome=dat$beta.outcome, se.outcome=dat$se.outcome,
                    effect_allele.outcome=dat$effect_allele.outcome, other_allele.outcome=dat$other_allele.outcome,
                    eaf.outcome=dat$eaf.outcome, outcome=dat$outcome[1], id.outcome=dat$id.outcome[1],
                    pval.outcome=dat$pval.outcome, mr_keep.outcome=TRUE, stringsAsFactors=FALSE)
  harmonise_data(exp, out)
}

clean_theme <- theme_bw(base_size=9) + theme(legend.position="bottom", legend.text=element_text(size=7),
  legend.title=element_text(size=7), plot.tag=element_text(face="bold", size=11))

mk <- function(fn, oid, ...){
  d <- load_dat(oid)
  switch(fn,
    scatter = mr_scatter_plot(mr(d), d)[[1]],
    forest  = mr_forest_plot(mr_singlesnp(d))[[1]],
    funnel  = mr_funnel_plot(mr_singlesnp(d))[[1]],
    loo     = mr_leaveoneout_plot(mr_leaveoneout(d))[[1]])
}

save_panel <- function(fn, outfile, w, h){
  p1 <- mk(fn, ids[1]) + clean_theme + labs(tag="A", title="Nielsen 2018 (AFGen+)")
  p2 <- mk(fn, ids[2]) + clean_theme + labs(tag="B", title="FinnGen R7")
  g <- p1 + p2 + plot_layout(ncol=2)
  ggsave(file.path(outdir, outfile), g, width=w, height=h, dpi=300, device="tiff", compression="lzw")
  ggsave(file.path(outdir, sub(".tiff", ".png", outfile)), g, width=w, height=h, dpi=300)
  cat("saved", outfile, "\n")
}

save_panel("scatter", "Figure1_scatter_300dpi.tiff", 10, 5)
save_panel("forest",  "Figure2_forest_300dpi.tiff", 12, 10)
save_panel("funnel",  "Figure3_funnel_300dpi.tiff", 10, 5)
save_panel("loo",     "Figure4_loo_300dpi.tiff", 12, 10)
cat("ALL FIGURES DONE\n")
