setwd("c:/Users/jpaige/git/jittering")
load("savedOutput/global/BYM2_models_test_results.RData")
for (nm in c("fitMM","fitMD","fitMDM")) {
  res <- out[[nm]]
  nll <- res$opt$value
  if (is.null(nll)) nll <- res$opt$objective
  cat("NLL", nm, nll, "\n")
}
