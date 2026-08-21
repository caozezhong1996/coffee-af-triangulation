suppressMessages({library(data.table); library(coloc)})
base <- "data"
fg  <- fread(file.path(base, "ahr_fg_window.tsv"))
api <- fread(file.path(base, "ahr_api_associations.csv"))
nei <- fread(file.path(base, "ahr_nielsen_window.tsv"))

comp <- function(a) chartr("ACGT", "TGCA", a)
is_pal <- function(a1, a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")

harm <- function(eff_a, e1, e2, fr, ref, alt) {
  n <- length(eff_a); beta <- rep(NA_real_, n); eaf <- rep(NA_real_, n); ok <- rep(FALSE, n)
  for (i in seq_len(n)) {
    a1 <- e1[i]; a2 <- e2[i]; f <- fr[i]; rf <- ref[i]; al <- alt[i]
    if (is.na(a1) || is.na(a2) || is.na(rf) || is.na(al)) next
    if (a1 == al && a2 == rf)      { beta[i] <-  eff_a[i]; eaf[i] <- f;       ok[i] <- TRUE }
    else if (a1 == rf && a2 == al) { beta[i] <- -eff_a[i]; eaf[i] <- if (is.na(f)) NA else 1 - f; ok[i] <- TRUE }
    else if (comp(a1) == al && comp(a2) == rf) { beta[i] <-  eff_a[i]; eaf[i] <- f; ok[i] <- TRUE }
    else if (comp(a1) == rf && comp(a2) == al) { beta[i] <- -eff_a[i]; eaf[i] <- if (is.na(f)) NA else 1 - f; ok[i] <- TRUE }
  }
  list(beta = beta, eaf = eaf, ok = ok)
}
maf_from <- function(e) ifelse(is.na(e), NA_real_, pmin(e, 1 - e))

# coffee sdY from the full original api set (consistent with main analysis)
api_all <- fread(file.path(base, "api_associations.csv"))
cf0 <- api_all[id == "ukb-b-5237"]; cf0[, m := maf_from(eaf)]
sdy_est <- median(cf0$se * sqrt(2 * cf0$n * cf0$m * (1 - cf0$m)), na.rm = TRUE)
cat("coffee sdY estimate:", round(sdy_est, 3), "\n")

sub <- fg  # whole AHR window
sub <- sub[!duplicated(rsids)]
fg_maf <- maf_from(sub$af_alt)
csub <- api[id == "ukb-b-5237"][match(sub$rsids, rsid)]
h <- harm(csub$beta, csub$ea, csub$nea, csub$eaf, sub$ref, sub$alt)
keep <- h$ok & !is.na(csub$beta) & !is.na(sub$beta)
pal <- is_pal(sub$ref, sub$alt)
caf_maf <- maf_from(h$eaf); caf_maf[is.na(caf_maf)] <- fg_maf[is.na(caf_maf)]
ambig <- pal & (caf_maf > 0.42)
keep <- keep & !ambig & !is.na(caf_maf)

run_coloc <- function(af_beta, af_se, af_maf, N_af, s_af, tag, p12 = 5e-6) {
  k <- keep & !is.na(af_beta) & !is.na(af_se) & af_se > 0 & !is.na(af_maf) & af_maf > 0 & af_maf < 0.5
  idx <- which(k)
  if (length(idx) < 50) return(data.table(outcome = tag, nsnps = length(idx), note = "too few SNPs"))
  d1 <- list(snp = sub$rsids[idx], beta = csub$beta[idx], varbeta = (csub$se[idx])^2,
             type = "quant", N = 428860, sdY = sdy_est, MAF = caf_maf[idx])
  d2 <- list(snp = sub$rsids[idx], beta = af_beta[idx], varbeta = (af_se[idx])^2,
             type = "cc", N = N_af, s = s_af, MAF = af_maf[idx])
  res <- tryCatch(coloc.abf(d1, d2, p12 = p12), error = function(e) e)
  if (inherits(res, "error")) return(data.table(outcome = tag, nsnps = length(idx), note = conditionMessage(res)))
  sm <- res$summary; tr <- res$results
  ppcol <- grep("SNP.PP", names(tr), value = TRUE)
  top <- tr[which.max(tr[[ppcol[1]]]), ]
  data.table(outcome = tag, nsnps = length(idx),
             PP.H0 = sm["PP.H0.abf"], PP.H1 = sm["PP.H1.abf"], PP.H2 = sm["PP.H2.abf"],
             PP.H3 = sm["PP.H3.abf"], PP.H4 = sm["PP.H4.abf"],
             top_snp = as.character(top$snp), top_SNP.PP.H4 = as.numeric(top[[ppcol[1]]]))
}

res <- list()
res[["fg_default"]]  <- run_coloc(sub$beta, sub$sebeta, fg_maf, 287805, 55853/287805, "FinnGenR11")
nsub <- nei[match(sub$rsids, rsid)]
hn <- harm(nsub$Effect_A2, nsub$A2, nsub$A1, suppressWarnings(as.numeric(nsub$Freq_A2)), sub$ref, sub$alt)
n_maf <- maf_from(hn$eaf); n_maf[is.na(n_maf)] <- caf_maf[is.na(n_maf)]
res[["nei_default"]] <- run_coloc(hn$beta, nsub$StdErr, n_maf, 1030836, 60620/1030836, "Nielsen2018")
# prior sensitivity
for (p12 in c(1e-5, 5e-5)) {
  res[[paste0("fg_", p12)]]  <- run_coloc(sub$beta, sub$sebeta, fg_maf, 287805, 55853/287805, paste0("FinnGenR11_p12=", p12), p12)
  res[[paste0("nei_", p12)]] <- run_coloc(hn$beta, nsub$StdErr, n_maf, 1030836, 60620/1030836, paste0("Nielsen2018_p12=", p12), p12)
}
out <- rbindlist(res, fill = TRUE)
out <- cbind(window = "AHR_rs4410790", out)
fwrite(out, file.path(base, "coloc_ahr_results.csv"))
print(out[, .(outcome, nsnps, PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3), top_snp, top_SNP.PP.H4 = round(top_SNP.PP.H4, 3))])
