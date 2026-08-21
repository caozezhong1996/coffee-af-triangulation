suppressMessages({library(data.table); library(coloc)})

base <- "data"
fg  <- fread(file.path(base, "windows_thinned.tsv"))          # FinnGen: #chrom,pos,ref,alt,rsids,beta,sebeta,af_alt,pval,window
api <- fread(file.path(base, "api_associations.csv"))          # OpenGWAS: id,rsid,ea,nea,eaf,beta,se,p,n
nei <- fread(file.path(base, "nielsen_windows.tsv"))           # Nielsen: rsid,A1,A2,Freq_A2,Effect_A2,StdErr,Pvalue,CHR,POS_GRCh37

comp <- function(a) chartr("ACGT", "TGCA", a)
is_pal <- function(a1, a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")

# harmonize one source to FinnGen alt allele. returns beta aligned to alt, eaf aligned to alt (or NA)
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

# coffee sdY estimate from se: sdY = se * sqrt(2*N*maf*(1-maf))
cf <- api[id == "ukb-b-5237"]
cf[, m := maf_from(eaf)]
sdy_est <- median(cf$se * sqrt(2 * cf$n * cf$m * (1 - cf$m)), na.rm = TRUE)
cat("coffee sdY estimate:", round(sdy_est, 3), "\n")

windows <- unique(fg$window)
results <- list()

for (wn in windows) {
  sub <- fg[window == wn]
  # FinnGen AF dataset (alt-aligned already)
  fg_maf <- maf_from(sub$af_alt)
  # --- coffee ---
  csub <- api[id == "ukb-b-5237"][match(sub$rsids, rsid)]
  h <- harm(csub$beta, csub$ea, csub$nea, csub$eaf, sub$ref, sub$alt)
  keep <- h$ok & !is.na(csub$beta) & !is.na(sub$beta)
  # drop palindromes with ambiguous freq
  pal <- is_pal(sub$ref, sub$alt)
  caf_maf <- maf_from(h$eaf); caf_maf[is.na(caf_maf)] <- fg_maf[is.na(caf_maf)]
  ambig <- pal & (caf_maf > 0.42)
  keep <- keep & !ambig & !is.na(caf_maf)

  run_coloc <- function(af_beta, af_se, af_maf, N_af, s_af, tag) {
    k <- keep & !is.na(af_beta) & !is.na(af_se) & af_se > 0 & !is.na(af_maf) & af_maf > 0 & af_maf < 0.5
    idx <- which(k)
    if (length(idx) < 50) return(data.table(window = wn, outcome = tag, nsnps = length(idx), note = "too few SNPs"))
    d1 <- list(snp = sub$rsids[idx], beta = csub$beta[idx], varbeta = (csub$se[idx])^2,
               type = "quant", N = 428860, sdY = sdy_est, MAF = caf_maf[idx])
    d2 <- list(snp = sub$rsids[idx], beta = af_beta[idx], varbeta = (af_se[idx])^2,
               type = "cc", N = N_af, s = s_af, MAF = af_maf[idx])
    res <- tryCatch(coloc.abf(d1, d2), error = function(e) e)
    if (inherits(res, "error")) return(data.table(window = wn, outcome = tag, nsnps = length(idx), note = conditionMessage(res)))
    sm <- res$summary
    tr <- res$results
    ppcol <- grep("SNP.PP", names(tr), value = TRUE)
    if (length(ppcol) > 0) {
      top <- tr[which.max(tr[[ppcol[1]]]), ]
      top_snp <- as.character(top$snp); top_pp <- as.numeric(top[[ppcol[1]]])
    } else { top_snp <- NA_character_; top_pp <- NA_real_ }
    data.table(window = wn, outcome = tag, nsnps = length(idx),
               PP.H0 = sm["PP.H0.abf"], PP.H1 = sm["PP.H1.abf"], PP.H2 = sm["PP.H2.abf"],
               PP.H3 = sm["PP.H3.abf"], PP.H4 = sm["PP.H4.abf"],
               top_snp = top_snp, top_SNP.PP.H4 = top_pp)
  }

  # FinnGen AF (alt-aligned): beta= sub$beta, se=sub$sebeta
  results[[paste0(wn, "_fg")]] <- run_coloc(sub$beta, sub$sebeta, fg_maf, 287805, 55853/287805, "FinnGenR11")

  # Nielsen AF (local): align Effect_A2 (allele A2) to alt
  nsub <- nei[match(sub$rsids, rsid)]
  hn <- harm(nsub$Effect_A2, nsub$A2, nsub$A1, suppressWarnings(as.numeric(nsub$Freq_A2)), sub$ref, sub$alt)
  n_maf <- maf_from(hn$eaf); n_maf[is.na(n_maf)] <- caf_maf[is.na(n_maf)]
  results[[paste0(wn, "_nei")]] <- run_coloc(hn$beta, nsub$StdErr, n_maf, 1030836, 60620/1030836, "Nielsen2018")
}

out <- rbindlist(results, fill = TRUE)
fwrite(out, file.path(base, "coloc_results.csv"))
print(out[, .(window, outcome, nsnps, PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3), top_snp, top_SNP.PP.H4 = round(top_SNP.PP.H4, 3))])
