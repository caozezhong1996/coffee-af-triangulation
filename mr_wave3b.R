jwt <- "YOUR_OPENGWAS_TOKEN"
Sys.setenv(OPENGWAS_JWT = jwt)
suppressMessages(library(ieugwasr))
th <- tryCatch(tophits("ukb-b-17814", pval=5e-8, clump=1), error=function(e) paste("ERR:", conditionMessage(e)))
if (is.character(th)) cat(th, "\n") else { cat("hits:", nrow(th), "\n"); print(head(th[, c("rsid","chr","position","p")], 15)) }
th2 <- tryCatch(tophits("ukb-b-17814", pval=5e-6, clump=1), error=function(e) paste("ERR:", conditionMessage(e)))
if (is.character(th2)) cat(th2, "\n") else cat("at 5e-6:", nrow(th2), "hits\n")
