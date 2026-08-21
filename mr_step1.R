Sys.setenv(OPENGWAS_JWT = "YOUR_OPENGWAS_TOKEN")
suppressMessages(library(ieugwasr))
log <- file(".//mr_search_log.txt", "w")
sink(log, type="output"); sink(log, type="message")
g <- tryCatch(gwasinfo(), error=function(e){cat("ERR:",conditionMessage(e),"\n"); NULL})
if(is.data.frame(g)){
  cat("datasets:", nrow(g), "\n\n=== COFFEE ===\n")
  print(g[grepl("coffee", g$trait, ignore.case=TRUE), c("id","trait","year","sample_size","author")])
  cat("\n=== ATRIAL FIBRILLATION ===\n")
  print(g[grepl("atrial fibrillation", g$trait, ignore.case=TRUE), c("id","trait","year","sample_size","ncase")])
}
sink(type="message"); sink(); close(log)
