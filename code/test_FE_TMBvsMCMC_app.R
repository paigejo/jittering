# Comparison test: FE+nugget+GH for MICS-only and DHS-only,
# fitted with both TMB optimization and TMBstan MCMC,
# using the application data (Nigeria DHS+MICS women's secondary education).
#
# 4 fits:
#   (1) MICS-only FE+nugget GH,  TMB     (fitFEM,  doMCMC=FALSE)
#   (2) DHS-only  FE+nugget GH,  TMB     (fitFED,  doMCMC=FALSE)
#   (3) MICS-only FE+nugget GH,  TMBstan (fitFEM,  doMCMC=TRUE)
#   (4) DHS-only  FE+nugget GH,  TMBstan (fitFED,  doMCMC=TRUE)
#
# Each fit is cached individually so the script is resumable.

source("code/setup.R")

# ---- config ----
Qgh    <- 10
KMICS  <- 100
KDHSu  <- 16
KDHSr  <- 21
mcIter <- 2000
mcWarm <- 1000
mcChn  <- 2

cacheDir <- "savedOutput/test"
if(!dir.exists(cacheDir)) dir.create(cacheDir, recursive=TRUE)

cacheTMBM  <- "savedOutput/global/fitFEM_final.RData"   # reuse if present
cacheTMBD  <- "savedOutput/global/fitFED_final.RData"   # reuse if present
cacheMCMCM <- file.path(cacheDir, "fitFEM_MCMC_app.RData")
cacheMCMCD <- file.path(cacheDir, "fitFED_MCMC_app.RData")

# ---- load datasets ----
load("savedOutput/global/ed.RData")
load("savedOutput/global/edMICS.RData")
edMICS <- sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

cat("DHS  n =", nrow(ed), "    MICS n =", nrow(edMICS), "\n")
cat("KMICS=", KMICS, "  KDHSu=", KDHSu, "  KDHSr=", KDHSr, "  Qgh=", Qgh, "\n",
    "MCMC: chains=", mcChn, " iter=", mcIter, " warmup=", mcWarm, "\n\n", sep="")

fit_or_load <- function(file, fitExpr, label) {
    if(file.exists(file)) {
        cat("[", label, "] loading cache:", file, "\n")
        e <- new.env(); load(file, envir=e)
        return(e$res)
    }
    cat("[", label, "] fitting ...\n")
    res <- eval(fitExpr)
    save(res, file=file)
    cat("[", label, "] saved to:", file, "\n")
    res
}

# ---- (1) MICS-only TMB ----
res_FEM_TMB <- fit_or_load(
    cacheTMBM,
    quote(fitFEM(datDHS=ed, datMICS=edMICS,
                 KMICS=KMICS, Qgh=Qgh,
                 getSDs=TRUE, verbose=FALSE, doMCMC=FALSE)),
    "FEM TMB")

# ---- (2) DHS-only TMB ----
res_FED_TMB <- fit_or_load(
    cacheTMBD,
    quote(fitFED(datDHS=ed,
                 KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                 getSDs=TRUE, verbose=FALSE, doMCMC=FALSE)),
    "FED TMB")

# ---- (3) MICS-only MCMC ----
res_FEM_MCMC <- fit_or_load(
    cacheMCMCM,
    quote(fitFEM(datDHS=ed, datMICS=edMICS,
                 KMICS=KMICS, Qgh=Qgh,
                 getSDs=FALSE, verbose=FALSE,
                 doMCMC=TRUE,
                 mcmcIter=mcIter, mcmcWarmup=mcWarm, mcmcChains=mcChn)),
    "FEM MCMC")

# ---- (4) DHS-only MCMC ----
res_FED_MCMC <- fit_or_load(
    cacheMCMCD,
    quote(fitFED(datDHS=ed,
                 KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                 getSDs=FALSE, verbose=FALSE,
                 doMCMC=TRUE,
                 mcmcIter=mcIter, mcmcWarmup=mcWarm, mcmcChains=mcChn)),
    "FED MCMC")

# ---- timings ----
cat("\nTimings (minutes):\n")
tt <- data.frame(
    model = c("FEM_TMB", "FED_TMB", "FEM_MCMC", "FED_MCMC"),
    fit_min   = c(res_FEM_TMB$totalTime, res_FED_TMB$totalTime,
                  res_FEM_MCMC$totalTime, res_FED_MCMC$totalTime) / 60,
    sdrep_min = c(res_FEM_TMB$sdTime, res_FED_TMB$sdTime,
                  res_FEM_MCMC$sdTime, res_FED_MCMC$sdTime) / 60
)
print(tt, digits=4, row.names=FALSE)

# ---- helpers ----
tmb_fixed <- function(res) {
    SD <- res$TMBsd
    if(!inherits(SD, "sdreport")) return(NULL)
    ss <- summary(SD, "fixed")
    data.frame(par=rownames(ss), est=ss[,"Estimate"], sd=ss[,"Std. Error"],
               row.names=NULL, stringsAsFactors=FALSE)
}

mcmc_summary <- function(res) {
    sf <- res$TMBsd
    if(!inherits(sf, "stanfit")) return(NULL)
    s <- rstan::summary(sf)$summary
    keep <- !grepl("^lp__", rownames(s))
    s <- s[keep, , drop=FALSE]
    data.frame(par=rownames(s),
               mean = s[,"mean"], sd = s[,"sd"],
               q025 = s[,"2.5%"], q975 = s[,"97.5%"],
               n_eff = s[,"n_eff"], Rhat = s[,"Rhat"],
               row.names=NULL, stringsAsFactors=FALSE)
}

label_betas <- function(df, covNames) {
    isB <- grepl("^beta", df$par)
    if(any(isB) && !is.null(covNames) && length(covNames) == sum(isB)) {
        df$par[isB] <- paste0("beta[", covNames, "]")
    }
    df
}

# ---- per-model comparison: TMB vs MCMC ----
cmp <- function(tmb, mcmc, covNames, label) {
    cat("\n========== ", label, "  (TMB vs MCMC)  ==========\n", sep="")
    if(is.null(tmb) || is.null(mcmc)) { cat("  missing piece, skipping\n"); return(NULL) }
    tmb  <- label_betas(tmb,  covNames)
    mcmc <- label_betas(mcmc, covNames)
    keep <- intersect(tmb$par, mcmc$par)
    tmb  <- tmb[match(keep, tmb$par), , drop=FALSE]
    mcmc <- mcmc[match(keep, mcmc$par), , drop=FALSE]
    out <- data.frame(par   = keep,
                      tmb_est  = tmb$est,
                      tmb_sd   = tmb$sd,
                      mcmc_mean= mcmc$mean,
                      mcmc_sd  = mcmc$sd,
                      mcmc_q025= mcmc$q025,
                      mcmc_q975= mcmc$q975,
                      Rhat     = mcmc$Rhat,
                      n_eff    = mcmc$n_eff,
                      diff_est = mcmc$mean - tmb$est,
                      ratio_sd = mcmc$sd / tmb$sd)
    print(out, digits=4, row.names=FALSE)
    invisible(out)
}

cmpM <- cmp(tmb_fixed(res_FEM_TMB), mcmc_summary(res_FEM_MCMC),
            res_FEM_TMB$covNames, "MICS-only FE+nug+GH")
cmpD <- cmp(tmb_fixed(res_FED_TMB), mcmc_summary(res_FED_MCMC),
            res_FED_TMB$covNames, "DHS-only FE+nug+GH")

# ---- save combined comparison ----
saveRDS(list(timings=tt, MICS=cmpM, DHS=cmpD,
             config=list(Qgh=Qgh, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                         mcmcIter=mcIter, mcmcWarmup=mcWarm, mcmcChains=mcChn)),
        file=file.path(cacheDir, "FE_TMBvsMCMC_app_summary.rds"))

cat("\nDone. Summary saved to ",
    file.path(cacheDir, "FE_TMBvsMCMC_app_summary.rds"), "\n", sep="")
