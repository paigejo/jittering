# Fit the commensurate-priors fusion model on the application data with BYM2 enabled.
# This is the "real" run — MICS now contributes to the shared spatial random field u
# while the commensurate prior protects beta from MICS misspecification.
#
# Cached fit: savedOutput/test/fitFEMD_comm_BYM2_app.RData

source("code/setup.R")
source("code/modFEMD_comm.R")

Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21

cacheDir  <- "savedOutput/test"
if(!dir.exists(cacheDir)) dir.create(cacheDir, recursive=TRUE)
cacheFile <- file.path(cacheDir, "fitFEMD_comm_BYM2_app.RData")

load("savedOutput/global/ed.RData")
load("savedOutput/global/edMICS.RData")
edMICS <- sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

cat("DHS  n =", nrow(ed), "    MICS n =", nrow(edMICS), "\n")
cat("KMICS=", KMICS, "  KDHSu=", KDHSu, "  KDHSr=", KDHSr, "  Qgh=", Qgh, "\n")

if(file.exists(cacheFile)) {
    cat("Loading cached fit:", cacheFile, "\n")
    load(cacheFile)
    res <- res_comm_bym2
} else {
    cat("Fitting BYM2-enabled commensurate-priors fusion (fixedEffectsOnly=FALSE)...\n")
    startTime <- proc.time()[3]
    res <- fitFEMD_comm(datDHS=ed, datMICS=edMICS,
                        KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                        fixedEffectsOnly=FALSE,
                        getSDs=TRUE, verbose=TRUE)
    cat("Total wall time:", round((proc.time()[3] - startTime)/60, 2), "min\n")
    res_comm_bym2 <- res
    save(res_comm_bym2, file=cacheFile)
    cat("Saved to:", cacheFile, "\n")
}

# ---- summary ----
cat("\n========== fitFEMD_comm (BYM2 enabled) =====\n")
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
    if("log_tau" %in% names(opt$par)) {
        lt <- opt$par["log_tau"]
        cat("log_tau =", round(lt, 4),
            "  =>  sigmaBYM2 =", round(exp(-0.5*lt), 4), "\n")
    }
    if("logit_phi" %in% names(opt$par)) {
        lp <- opt$par["logit_phi"]
        cat("logit_phi =", round(lp, 4),
            "  =>  phi =", round(1/(1+exp(-lp)), 4), "\n")
    }
}

if(inherits(res$TMBsd, "sdreport")) {
    SD <- res$TMBsd
    cat("pdHess =", SD$pdHess, "\n")

    fp <- summary(SD, "fixed")
    cat("\n----- Outer hyperparameters (with SE) -----\n")
    print(round(fp, 4))

    cv <- res$covNames
    nB <- length(cv)
    re <- summary(SD, "random")
    rn <- rownames(re)

    pickRows <- function(name) re[rn == name, , drop=FALSE]
    aD <- pickRows("alpha");   aM <- pickRows("alpha_M")
    bD <- pickRows("beta");    bM <- pickRows("beta_M")

    if(!is.null(bD) && nrow(bD) == nB) {
        tab <- data.frame(
            par             = c("alpha", paste0("beta[", cv, "]")),
            DHS_est         = c(aD[1,1], bD[,1]),
            DHS_sd          = c(aD[1,2], bD[,2]),
            MICS_est        = c(aM[1,1], bM[,1]),
            MICS_sd         = c(aM[1,2], bM[,2]),
            diff_M_minus_D  = c(aM[1,1] - aD[1,1], bM[,1] - bD[,1])
        )
        cat("\n----- DHS-side vs MICS-side MAP estimates (BYM2 enabled) -----\n")
        print(tab, row.names=FALSE, digits=4)
    }
}

# ---- comparison vs FE-only commensurate ----
feOnlyCache <- "savedOutput/test/fitFEMD_comm_app.RData"
if(file.exists(feOnlyCache) && inherits(res$TMBsd, "sdreport")) {
    eF <- new.env(); load(feOnlyCache, envir=eF); resFE <- eF$res_comm
    if(inherits(resFE$TMBsd, "sdreport")) {
        rF  <- summary(resFE$TMBsd, "random"); rnF <- rownames(rF)
        rB  <- summary(res$TMBsd,   "random"); rnB <- rownames(rB)

        pickF <- function(nm) rF[rnF == nm, , drop=FALSE]
        pickB <- function(nm) rB[rnB == nm, , drop=FALSE]

        cv <- res$covNames
        cmp <- data.frame(
            par           = c("alpha", paste0("beta[", cv, "]")),
            FEonly_DHS    = c(pickF("alpha")[1,1],   pickF("beta")[,1]),
            FEonly_MICS   = c(pickF("alpha_M")[1,1], pickF("beta_M")[,1]),
            BYM2_DHS      = c(pickB("alpha")[1,1],   pickB("beta")[,1]),
            BYM2_MICS     = c(pickB("alpha_M")[1,1], pickB("beta_M")[,1])
        )
        cat("\n----- FE-only vs BYM2-enabled commensurate fusion -----\n")
        print(cmp, row.names=FALSE, digits=4)

        cat("\nsigma_comm (FE-only) =",
            round(exp(resFE$opt$par["log_sigma_comm"]), 4), "\n")
        cat("sigma_comm (BYM2)    =",
            round(exp(res$opt$par["log_sigma_comm"]), 4), "\n")
    }
}

cat("\nDone.\n")
