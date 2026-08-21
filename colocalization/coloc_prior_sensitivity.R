suppressMessages({library(data.table); library(coloc)})

base <- "data"
fg  <- fread(file.path(base, "windows_thinned.tsv"))
api <- fread(file.path(base, "api_associations.csv"))
nei <- fread(file.path(base, "nielsen_windows.tsv"))

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

cf <- api[id == "ukb-b-5237"]
cf[, m := maf_from(eaf)]
sdy_est <- median(cf$se * sqrt(2 * cf$n * cf$m * (1 - cf$m)), na.rm = TRUE)

# ---- build harmonized dataset lists once per window x outcome ----
build_pair <- function(wn) {
  sub <- fg[window == wn]
  fg_maf <- maf_from(sub$af_alt)
  csub <- api[id == "ukb-b-5237"][match(sub$rsids, rsid)]
  h <- harm(csub$beta, csub$ea, csub$nea, csub$eaf, sub$ref, sub$alt)
  keep <- h$ok & !is.na(csub$beta) & !is.na(sub$beta)
  pal <- is_pal(sub$ref, sub$alt)
  caf_maf <- maf_from(h$eaf); caf_maf[is.na(caf_maf)] <- fg_maf[is.na(caf_maf)]
  ambig <- pal & (caf_maf > 0.42)
  keep <- keep & !ambig & !is.na(caf_maf)

  mk <- function(af_beta, af_se, af_maf, N_af, s_af) {
    k <- keep & !is.na(af_beta) & !is.na(af_se) & af_se > 0 & !is.na(af_maf) & af_maf > 0 & af_maf < 0.5
    idx <- which(k)
    if (length(idx) < 50) return(NULL)
    list(d1 = list(snp = sub$rsids[idx], beta = csub$beta[idx], varbeta = (csub$se[idx])^2,
                   type = "quant", N = 428860, sdY = sdy_est, MAF = caf_maf[idx]),
         d2 = list(snp = sub$rsids[idx], beta = af_beta[idx], varbeta = (af_se[idx])^2,
                   type = "cc", N = N_af, s = s_af, MAF = af_maf[idx]),
         nsnps = length(idx))
  }

  nsub <- nei[match(sub$rsids, rsid)]
  hn <- harm(nsub$Effect_A2, nsub$A2, nsub$A1, suppressWarnings(as.numeric(nsub$Freq_A2)), sub$ref, sub$alt)
  n_maf <- maf_from(hn$eaf); n_maf[is.na(n_maf)] <- caf_maf[is.na(n_maf)]

  list(fg  = mk(sub$beta, sub$sebeta, fg_maf, 287805, 55853/287805),
       nei = mk(hn$beta, nsub$StdErr, n_maf, 1030836, 60620/1030836))
}

P12 <- c(5e-6, 1e-5, 5e-5)
results <- list()
for (wn in unique(fg$window)) {
  pair <- build_pair(wn)
  for (tag in c("fg", "nei")) {
    obj <- pair[[tag]]
    if (is.null(obj)) next
    outcome <- ifelse(tag == "fg", "FinnGenR11", "Nielsen2018")
    for (p12 in P12) {
      res <- tryCatch(coloc.abf(obj$d1, obj$d2, p12 = p12), error = function(e) e)
      if (inherits(res, "error")) {
        results[[length(results)+1]] <- data.table(window = wn, outcome = outcome, p12 = p12,
                                                   nsnps = obj$nsnps, note = conditionMessage(res))
        next
      }
      sm <- res$summary
      results[[length(results)+1]] <- data.table(
        window = wn, outcome = outcome, p12 = p12, nsnps = obj$nsnps,
        PP.H0 = sm["PP.H0.abf"], PP.H1 = sm["PP.H1.abf"], PP.H2 = sm["PP.H2.abf"],
        PP.H3 = sm["PP.H3.abf"], PP.H4 = sm["PP.H4.abf"])
    }
  }
}

out <- rbindlist(results, fill = TRUE)
out[, class := fifelse(PP.H4 >= 0.7, "shared (H4>=0.7)",
                fifelse(PP.H3 >= 0.7, "distinct (H3>=0.7)",
                 fifelse(PP.H4 < 0.05, "no shared (H4<0.05)", "indeterminate")))]
fwrite(out, file.path(base, "coloc_prior_sensitivity.csv"))
show <- out[, .(window, outcome, nsnps, p12,
                PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3), class)]
print(show)
