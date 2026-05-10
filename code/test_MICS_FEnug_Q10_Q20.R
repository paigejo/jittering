#!/usr/bin/env Rscript
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_FEnug.R")

covs <- c("urban","access","elev","distRiversLakes","normPop")
KMICS <- 100

cat("Preparing inputs with KMICS=", KMICS, "\n")
inputsMDM <- makeInputsMDM(datDHS=ed, datMICS=edMICS, KMICS=KMICS)
cat("Available MICS covariates:", paste(colnames(inputsMDM$intPtsMICS$XUrb), collapse=", "), "\n")

run_fit <- function(Qgh) {
  cat("\nRunning fit Qgh=", Qgh, "\n")
  fit <- fitMFEM(datDHS=ed, datMICS=edMICS, inputsMDM=inputsMDM, covariates=covs, KMICS=KMICS, Qgh=Qgh, fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
  if(is.null(fit$TMBsd)) stop(paste0("sdreport failed for Q=", Qgh))
  pe <- summary(fit$TMBsd, "fixed")
  alpha <- pe["alpha","Estimate"]; alpha_se <- pe["alpha","Std. Error"]
  beta_idx <- grep("^beta", rownames(pe))
  beta_est <- pe[beta_idx,"Estimate"]; beta_se <- pe[beta_idx,"Std. Error"]
  beta_names <- fit$covNames
  sigmaEps <- NA
  if("log_tauEps" %in% rownames(pe)) {
    sigmaEps <- exp(-0.5 * pe["log_tauEps","Estimate"])
  }
  NLL <- fit$opt$objective
  list(alpha=alpha, alpha_se=alpha_se, beta_est=beta_est, beta_se=beta_se, beta_names=beta_names, sigmaEps=sigmaEps, NLL=NLL)
}

res10 <- run_fit(10)
res20 <- run_fit(20)

params <- c("alpha", paste0("beta[", covs, "]"))
est10 <- c(res10$alpha, res10$beta_est)
se10  <- c(res10$alpha_se, res10$beta_se)
est20 <- c(res20$alpha, res20$beta_est)
se20  <- c(res20$alpha_se, res20$beta_se)

df <- data.frame(Parameter=params,
                 Q10_Est=round(est10,4), Q10_SE=round(se10,4),
                 Q20_Est=round(est20,4), Q20_SE=round(se20,4),
                 stringsAsFactors=FALSE)

cat("\nParameter estimates (Q=10 vs Q=20):\n")
print(df, row.names=FALSE)

cat("\nSummary stats:\n")
sumdf <- data.frame(Parameter=c("sigmaEps","NLL"),
                    Q10=round(c(res10$sigmaEps, res10$NLL),4),
                    Q20=round(c(res20$sigmaEps, res20$NLL),4))
print(sumdf, row.names=FALSE)

cat("\nDone.\n")
