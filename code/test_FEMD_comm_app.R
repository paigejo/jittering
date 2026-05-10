# Fit the commensurate-priors fusion model on the application data and compare
# beta_D vs beta_M, plus the EB-selected sigma_comm.
#
# Cached fit: savedOutput/test/fitFEMD_comm_app.RData
# Compares against the cached MICS-only and DHS-only FE+nug+GH fits when available.

source("code/setup.R")
source("code/modFEMD_comm.R")

Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21

cacheDir  <- "savedOutput/test"
if(!dir.exists(cacheDir)) dir.create(cacheDir, recursive=TRUE)
cacheFile <- file.path(cacheDir, "fitFEMD_comm_app.RData")

# ---- load datasets ----
load("savedOutput/global/ed.RData")
load("savedOutput/global/edMICS.RData")
edMICS <- sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

cat("DHS  n =", nrow(ed), "    MICS n =", nrow(edMICS), "\n")
cat("KMICS=", KMICS, "  KDHSu=", KDHSu, "  KDHSr=", KDHSr, "  Qgh=", Qgh, "\n",
    sep="")

# ---- fit (or load) ----
if(file.exists(cacheFile)) {
    cat("[fitFEMD_comm] loading cache:", cacheFile, "\n")
    load(cacheFile)
    res <- res_comm
} else {
    cat("[fitFEMD_comm] fitting ...\n")
    startTime <- proc.time()[3]
    res <- fitFEMD_comm(datDHS=ed, datMICS=edMICS,
                        KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                        getSDs=TRUE, verbose=FALSE)
    cat("Total wall time:", round((proc.time()[3] - startTime)/60, 2), "min\n")
    res_comm <- res
    save(res_comm, file=cacheFile)
    cat("[fitFEMD_comm] saved to:", cacheFile, "\n")
}

# ---- summary ----
cat("\n========== fitFEMD_comm =====\n")
cat("method =", res$method,
    "  fitTime(min) =", round(res$totalTime/60, 2),
    "  sdTime(min) =",  round(res$sdTime/60, 2), "\n")
opt <- res$opt
if(!is.null(opt)) {
    cat("convergence =", opt$convergence,
        "   nll =", round(opt$objective, 4), "\n")
    cat("Outer optimum:\n"); print(round(opt$par, 4))
    if("log_sigma_comm" %in% names(opt$par)) {
        ls <- opt$par["log_sigma_comm"]
        cat("\nlog_sigma_comm =", round(ls, 4),
            "  =>  sigma_comm =", round(exp(ls), 4), "\n")
    }
}

if(inherits(res$TMBsd, "sdreport")) {
    SD <- res$TMBsd
    cat("pdHess =", SD$pdHess, "\n")

    cv <- res$covNames
    nB <- length(cv)

    # outer hyperparameters
    fp <- summary(SD, "fixed")
    cat("\n----- Outer hyperparameters (with SE) -----\n")
    print(round(fp, 4))

    # random-effect MAPs (alpha, alpha_M, beta, beta_M)
    re <- summary(SD, "random")
    rn <- rownames(re)

    pickRows <- function(name, n) {
        idx <- which(rn == name)
        if(length(idx) == 0) return(NULL)
        re[idx, , drop=FALSE]
    }

    alphaD <- pickRows("alpha", 1)
    alphaM <- pickRows("alpha_M", 1)
    betaD  <- pickRows("beta",   nB)
    betaM  <- pickRows("beta_M", nB)

    if(!is.null(betaD) && nrow(betaD) == nB) {
        tab <- data.frame(
            par             = c("alpha", paste0("beta[", cv, "]")),
            DHS_est         = c(alphaD[1,1], betaD[,1]),
            DHS_sd          = c(alphaD[1,2], betaD[,2]),
            MICS_est        = c(alphaM[1,1], betaM[,1]),
            MICS_sd         = c(alphaM[1,2], betaM[,2]),
            diff_M_minus_D  = c(alphaM[1,1] - alphaD[1,1],
                                betaM[,1]  - betaD[,1])
        )
        cat("\n----- DHS-side vs MICS-side MAP estimates -----\n")
        print(tab, row.names=FALSE, digits=4)
    }
}

# ---- comparison vs separate-survey fits ----
compareCacheM <- "savedOutput/global/fitFEM_final.RData"
compareCacheD <- "savedOutput/global/fitFED_final.RData"
if(file.exists(compareCacheM) && file.exists(compareCacheD) &&
   inherits(res$TMBsd, "sdreport")) {
    eM <- new.env(); load(compareCacheM, envir=eM); resM <- eM$res
    eD <- new.env(); load(compareCacheD, envir=eD); resD <- eD$res
    ssM <- summary(resM$TMBsd, "fixed")
    ssD <- summary(resD$TMBsd, "fixed")

    cv  <- res$covNames
    nB  <- length(cv)
    re  <- summary(res$TMBsd, "random")
    rn  <- rownames(re)
    betaD_comm <- re[rn == "beta",   , drop=FALSE]
    betaM_comm <- re[rn == "beta_M", , drop=FALSE]
    alphaD_comm <- re[rn == "alpha",   , drop=FALSE]
    alphaM_comm <- re[rn == "alpha_M", , drop=FALSE]

    # ssM/ssD have rownames "alpha","beta","beta",... pull alpha first then 5 beta rows
    rnM <- rownames(ssM); rnD <- rownames(ssD)
    aM  <- ssM[rnM == "alpha", , drop=FALSE]
    aD  <- ssD[rnD == "alpha", , drop=FALSE]
    bM  <- ssM[rnM == "beta",  , drop=FALSE]
    bD  <- ssD[rnD == "beta",  , drop=FALSE]

    cmp <- data.frame(
        par           = c("alpha", paste0("beta[", cv, "]")),
        MICSonly_est  = c(aM[1,1], bM[,1]),
        MICSonly_sd   = c(aM[1,2], bM[,2]),
        DHSonly_est   = c(aD[1,1], bD[,1]),
        DHSonly_sd    = c(aD[1,2], bD[,2]),
        comm_DHS_est  = c(alphaD_comm[1,1], betaD_comm[,1]),
        comm_DHS_sd   = c(alphaD_comm[1,2], betaD_comm[,2]),
        comm_MICS_est = c(alphaM_comm[1,1], betaM_comm[,1]),
        comm_MICS_sd  = c(alphaM_comm[1,2], betaM_comm[,2])
    )
    cat("\n----- All four sets of estimates side-by-side -----\n")
    print(cmp, row.names=FALSE, digits=4)
}

cat("\nDone.\n")
