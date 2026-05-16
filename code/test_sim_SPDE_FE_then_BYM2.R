# Two-stage sanity test on SPDE-simulated data:
#   Stage 1: FE+nug+GH for DHS, MICS, DHS+MICS (fitFED, fitFEM, fitFEMD).
#   Gate:    if any FE model badly misses a covariate effect, STOP and report.
#   Stage 2: BYM2+nug+GH for DHS, MICS, DHS+MICS (fitMD, fitMM, fitMDM).
#
# Truth (simStudy1 DGP per CLAUDE.md):
#   alpha=-1.25, beta=(urban=1.0, access=0, elev=0, distRiversLakes=0, normPop=0.5)
#   sigmaEpsilon = sqrt(1.5)  ;  for BYM2 spatial hyperparams: not part of the FE pass criterion
#
# Run from code/:  Rscript test_sim_SPDE_FE_then_BYM2.R
# Wall time: ~10-30 min for FE stage; +30-60 min if BYM2 stage runs.

source("setup.R")
options(error=recover)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/modM_DSep.R")     # fitMD  (BYM2 DHS-only)
source("code/modM_MSep.R")     # fitMM  (BYM2 MICS-only)
source("code/modM_DMSep.R")    # fitMDM (BYM2 fusion)
source("code/makeInputsTMB.R")

simIndex <- 1
Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
covs  <- c("urban","access","elev","distRiversLakes","normPop")

# Truth for the simStudy1 DGP.
true_alpha    <- -1.25
true_beta     <- c(urban=1.0, access=0, elev=0, distRiversLakes=0, normPop=0.5)
true_sigmaEps <- sqrt(1.5)
ZTHRESH       <- 3   # |z| > ZTHRESH on any beta => FE failure

# ---- load SPDE-simulated surveys (the unsuffixed file) ----
cat("Loading SPDE-simulated surveys (savedOutput/simStudy1/simPopsSurveys.RData)\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")
datDHS  <- surveysDHS[[simIndex]]
datMICS <- surveysMICS[[simIndex]]

nameTab <- rbind(c("N","n"), c("N","ns"), c("Z","y"), c("Z","ys"))
for(i in 1:nrow(nameTab)) {
    fromN <- nameTab[i,1]; toN <- nameTab[i,2]
    if(!(toN %in% names(datDHS)))  datDHS[[toN]]  <- datDHS[[fromN]]
    if(!(toN %in% names(datMICS))) datMICS[[toN]] <- datMICS[[fromN]]
}
if(!("Stratum" %in% names(datMICS))) datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

# ---- load integration points directly and pass them into makeInputsMDM ----
# Otherwise makeInputsMDM would (a) load DHS int pts only if datDHS is identical
# to the application `ed` (it isn't here — it's a sim survey), so it would fall
# through to makeAllIntegrationPointsDHS and trigger GDAL/terra. Passing the
# pre-built integration points explicitly avoids GDAL entirely.
# Prefer the K-suffixed DHS int pts file (built at the requested K=KDHSu/KDHSr).
# Unsuffixed files in savedOutput/simStudy1/ are the legacy K=11/16 build, which
# CLAUDE.md flags as biased; only the *_K<u>_<r>.RData variants are at K=16/21.
dhsIntPtsFile_K <- sprintf(
    "savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K%d_%d.RData",
    simIndex, KDHSu, KDHSr)
dhsIntPtsFile_default <- sprintf(
    "savedOutput/simStudy1/intPtsDHS_simStudy1_%d.RData", simIndex)
dhsIntPtsFile <- if(file.exists(dhsIntPtsFile_K)) dhsIntPtsFile_K else dhsIntPtsFile_default
micsIntPtsFile <- sprintf("savedOutput/global/intPtsMICS_%d.RData", KMICS)

cat("Loading DHS int pts: ", dhsIntPtsFile, "\n", sep="")
load(dhsIntPtsFile)    # expected to populate object `intPtsDHS`
cat("Loading MICS int pts:", micsIntPtsFile, "\n", sep="")
load(micsIntPtsFile)   # expected to populate object `intPtsMICS`
stopifnot(exists("intPtsDHS"), exists("intPtsMICS"))

# Hard-check that the loaded resolution matches the requested K, so we don't
# silently fit at the deprecated K=11/16.
gotKU <- ncol(intPtsDHS$wUrban); gotKR <- ncol(intPtsDHS$wRural)
if(gotKU != KDHSu || gotKR != KDHSr) {
    stop(sprintf(
        "DHS int pts at WRONG resolution: loaded K=%d/%d but requested KDHSu=%d, KDHSr=%d.\n  File: %s\n  Regenerate this file at the correct K or use a different simIndex.",
        gotKU, gotKR, KDHSu, KDHSr, dhsIntPtsFile))
}

# ---- build / load inputsMDM once ----
inputsFile <- sprintf(
    "savedOutput/simStudy1/inputs_SPDE_sim%d_KMICS%03d_KDHS%02d_%02d.RData",
    simIndex, KMICS, KDHSu, KDHSr)
if(file.exists(inputsFile)) {
    cat("Loading precomputed inputs:", inputsFile, "\n")
    load(inputsFile)
} else {
    cat("Building inputsMDM (one-shot, integration points passed in directly)\n")
    inputsMDM <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS,
                               intPtsDHS=intPtsDHS, intPtsMICS=intPtsMICS,
                               KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                               saveNewIntPts=FALSE)
    save(inputsMDM, file=inputsFile)
}
fixDhsNames <- function(mat) {
    if(is.null(mat)) return(mat)
    cn <- colnames(mat)
    if(!is.null(cn)) {
        cn <- gsub("^urb$", "urban", cn); cn <- gsub("^pop$", "normPop", cn)
        colnames(mat) <- cn
    }
    mat
}
inputsMDM$intPtsDHS$covsUrb <- fixDhsNames(inputsMDM$intPtsDHS$covsUrb)
inputsMDM$intPtsDHS$covsRur <- fixDhsNames(inputsMDM$intPtsDHS$covsRur)

# ---- run-one helper ----
runOne <- function(label, fitExpr) {
    cat("\n========== ", label, " ==========\n", sep="")
    t0 <- proc.time()[3]
    res <- try(eval(fitExpr), silent=FALSE)
    elapsed <- (proc.time()[3] - t0)/60
    if(inherits(res, "try-error")) {
        cat(sprintf("  FAILED to fit (%.2f min): %s\n", elapsed, attr(res,"condition")$message))
        return(NULL)
    }
    pdh <- if(inherits(res$TMBsd, "sdreport")) as.character(res$TMBsd$pdHess) else "n/a"
    conv <- if(!is.null(res$opt)) as.character(res$opt$convergence) else "n/a"
    cat(sprintf("  fit time: %.2f min   convergence: %s   pdHess: %s\n", elapsed, conv, pdh))
    res
}

# ---- pull (alpha, beta, sigmaEps) for FE models ----
pullFEests <- function(res) {
    if(is.null(res) || !inherits(res$TMBsd, "sdreport")) return(NULL)
    fe <- summary(res$TMBsd, "fixed")
    a  <- fe["alpha", , drop=FALSE]
    b  <- fe[grep("^beta", rownames(fe)), , drop=FALSE]
    sigEps <- if("log_tauEps" %in% rownames(fe)) exp(-0.5*fe["log_tauEps","Estimate"]) else NA
    list(alpha=a, beta=b, sigmaEps=sigEps,
         pdHess=res$TMBsd$pdHess, conv=res$opt$convergence)
}

# ============== STAGE 1 : FE+nug+GH ==============
cat("\n##############################################")
cat("\n##  STAGE 1: FE+nug+GH on SPDE-simulated data ")
cat("\n##############################################\n")

resD <- runOne("M_D    fitFED",
               quote(fitFED(datDHS=datDHS, inputsMDM=inputsMDM,
                            covariates=covs, Qgh=Qgh,
                            fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))
resM <- runOne("M_M    fitFEM",
               quote(fitFEM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                            covariates=covs, KMICS=KMICS, Qgh=Qgh,
                            fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))
resDM <- runOne("M_DM   fitFEMD",
                quote(fitFEMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                              covariates=covs, KMICS=KMICS, Qgh=Qgh,
                              fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))

# Build a side-by-side table of parameter recovery + max |z| per model.
parNames <- c("alpha", paste0("beta[", covs, "]"))
truthVec <- c(true_alpha, true_beta)

estTab <- data.frame(par=parNames, truth=round(truthVec, 3))
maxZ   <- list()
for(nm in c("M_D", "M_M", "M_DM")) {
    res <- list(M_D=resD, M_M=resM, M_DM=resDM)[[nm]]
    e   <- pullFEests(res)
    if(is.null(e)) {
        estTab[[paste0(nm,"_est")]] <- NA
        estTab[[paste0(nm,"_se")]]  <- NA
        estTab[[paste0(nm,"_z")]]   <- NA
        maxZ[[nm]] <- NA
        next
    }
    est <- c(e$alpha[,"Estimate"], e$beta[,"Estimate"])
    se  <- c(e$alpha[,"Std. Error"], e$beta[,"Std. Error"])
    z   <- (est - truthVec) / se
    estTab[[paste0(nm,"_est")]] <- round(est, 3)
    estTab[[paste0(nm,"_se")]]  <- round(se,  3)
    estTab[[paste0(nm,"_z")]]   <- round(z,   2)
    maxZ[[nm]] <- max(abs(z[-1]), na.rm=TRUE)   # exclude alpha from the gate; the user said covariate effects
}

cat("\n----- FE parameter recovery (truth vs estimate +/- SE; z = (est-truth)/SE) -----\n")
print(estTab, row.names=FALSE)

cat("\n----- nugget SD (truth =", round(true_sigmaEps,3), ") -----\n")
for(nm in c("M_D","M_M","M_DM")) {
    res <- list(M_D=resD, M_M=resM, M_DM=resDM)[[nm]]
    e   <- pullFEests(res)
    if(!is.null(e)) cat(sprintf("  %-5s  sigmaEps = %.3f\n", nm, e$sigmaEps))
}

cat(sprintf("\n----- max |z| on covariate betas per FE model (gate = %g) -----\n", ZTHRESH))
for(nm in names(maxZ)) cat(sprintf("  %-5s  max|z| = %s\n", nm, format(maxZ[[nm]], digits=3)))

passed <- !any(is.na(unlist(maxZ))) && all(unlist(maxZ) <= ZTHRESH)
if(!passed) {
    cat("\n>>> FE-stage FAILED: at least one model has max |z| >", ZTHRESH,
        "(or did not converge).\n>>> Stopping before BYM2 stage as instructed.\n\n")
    quit(save="no", status=0)
}

cat("\n>>> FE-stage PASSED. Proceeding to BYM2 stage.\n")

# ============== STAGE 2 : BYM2+nug+GH ==============
cat("\n#####################################################")
cat("\n##  STAGE 2: BYM2+nug+GH on SPDE-simulated data     ")
cat("\n#####################################################\n")

# fitMD / fitMM / fitMDM are the BYM2 GH variants. They expect the same
# inputsMDM the FE wrappers use; pass it through. They DO NOT take a
# fixedEffectsOnly arg — BYM2 is always on.
resD2 <- runOne("M_D BYM2  fitMD",
                quote(fitMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                            KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                            Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))
resM2 <- runOne("M_M BYM2  fitMM",
                quote(fitMM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                            KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))
resDM2 <- runOne("M_DM BYM2  fitMDM",
                 quote(fitMDM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                              KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                              Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))

# Pull (alpha, beta, sigmaEps, sigmaBYM2, phi) for BYM2 models.
pullBYM2ests <- function(res) {
    if(is.null(res) || !inherits(res$TMBsd, "sdreport")) return(NULL)
    fe <- summary(res$TMBsd, "fixed")
    rn <- rownames(fe)
    aRow <- if("alpha" %in% rn) fe["alpha", , drop=FALSE] else NULL
    bRow <- fe[grep("^beta", rn), , drop=FALSE]
    sigEps   <- if("log_tauEps" %in% rn)    exp(-0.5*fe["log_tauEps","Estimate"]) else NA
    sigBYM2  <- if("log_tau" %in% rn)       exp(-0.5*fe["log_tau","Estimate"])    else NA
    phiHat   <- if("logit_phi" %in% rn)     1/(1+exp(-fe["logit_phi","Estimate"])) else NA
    list(alpha=aRow, beta=bRow,
         sigmaEps=sigEps, sigmaBYM2=sigBYM2, phi=phiHat,
         pdHess=res$TMBsd$pdHess)
}

estTab2 <- data.frame(par=parNames, truth=round(truthVec, 3))
for(nm in c("M_D","M_M","M_DM")) {
    res <- list(M_D=resD2, M_M=resM2, M_DM=resDM2)[[nm]]
    e   <- pullBYM2ests(res)
    if(is.null(e) || is.null(e$alpha)) {
        estTab2[[paste0(nm,"_est")]] <- NA
        estTab2[[paste0(nm,"_se")]]  <- NA
        next
    }
    est <- c(e$alpha[,"Estimate"], e$beta[,"Estimate"])
    se  <- c(e$alpha[,"Std. Error"], e$beta[,"Std. Error"])
    estTab2[[paste0(nm,"_est")]] <- round(est, 3)
    estTab2[[paste0(nm,"_se")]]  <- round(se,  3)
}

cat("\n----- BYM2 parameter recovery -----\n")
print(estTab2, row.names=FALSE)

cat("\n----- BYM2 spatial / nugget hyperparameters -----\n")
for(nm in c("M_D","M_M","M_DM")) {
    res <- list(M_D=resD2, M_M=resM2, M_DM=resDM2)[[nm]]
    e   <- pullBYM2ests(res)
    if(is.null(e)) next
    cat(sprintf("  %-5s  sigmaEps=%.3f  sigmaBYM2=%.3f  phi=%.3f  pdHess=%s\n",
                nm, e$sigmaEps, e$sigmaBYM2, e$phi, as.character(e$pdHess)))
}

cat("\nDone.\n")
