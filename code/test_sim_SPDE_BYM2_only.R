# Stage-2-only follow-up to test_sim_SPDE_FE_then_BYM2.R: run BYM2+nug+GH for
# DHS, MICS, and DHS+MICS fusion on the SPDE-simulated data and report covariate
# recovery. Reuses the cached inputsMDM from the FE-stage run.

source("setup.R")
options(error=traceback)
source("code/modM_DSep.R")
source("code/modM_MSep.R")
source("code/modM_DMSep.R")
source("code/makeInputsTMB.R")

simIndex <- 1
Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
covs  <- c("urban","access","elev","distRiversLakes","normPop")

true_alpha    <- -1.25
true_beta     <- c(urban=1.0, access=0, elev=0, distRiversLakes=0, normPop=0.5)
true_sigmaEps <- sqrt(1.5)

# ---- load SPDE-simulated surveys (must be the same dataset as the FE stage) ----
cat("Loading SPDE-simulated surveys ...\n")
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

inputsFile <- sprintf(
    "savedOutput/simStudy1/inputs_SPDE_sim%d_KMICS%03d_KDHS%02d_%02d.RData",
    simIndex, KMICS, KDHSu, KDHSr)
if(file.exists(inputsFile)) {
    cat("Loading cached inputsMDM:", inputsFile, "\n")
    load(inputsFile)
} else {
    cat("No cached inputsMDM — building from scratch.\n")
    inputsMDM <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS, KMICS=KMICS,
                               KDHSurb=KDHSu, KDHSrur=KDHSr, saveNewIntPts=TRUE)
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

runOne <- function(label, fitExpr) {
    cat("\n========== ", label, " ==========\n", sep="")
    t0 <- proc.time()[3]
    res <- try(eval(fitExpr), silent=FALSE)
    elapsed <- (proc.time()[3] - t0)/60
    if(inherits(res, "try-error")) {
        cat(sprintf("  FAILED to fit (%.2f min): %s\n", elapsed, attr(res,"condition")$message))
        return(NULL)
    }
    pdh  <- if(inherits(res$TMBsd, "sdreport")) as.character(res$TMBsd$pdHess) else "n/a"
    conv <- if(!is.null(res$opt)) as.character(res$opt$convergence) else "n/a"
    cat(sprintf("  fit time: %.2f min   convergence: %s   pdHess: %s\n", elapsed, conv, pdh))
    res
}

cat("\n#####################################################")
cat("\n##  STAGE 2: BYM2+nug+GH on SPDE-simulated data     ")
cat("\n#####################################################\n")

resD <- runOne("M_D BYM2  fitMD",
               quote(fitMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                           KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                           Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))
resM <- runOne("M_M BYM2  fitMM",
               quote(fitMM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                           KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))
resDM <- runOne("M_DM BYM2  fitMDM",
                quote(fitMDM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                             KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                             Qgh=Qgh, getSDs=TRUE, verbose=FALSE)))

# ---- pull (alpha, beta, sigmaEps, sigmaBYM2, phi) for BYM2 models ----
pullBYM2 <- function(res) {
    if(is.null(res) || !inherits(res$TMBsd, "sdreport")) return(NULL)
    fe <- summary(res$TMBsd, "fixed")
    rn <- rownames(fe)
    aRow <- if("alpha" %in% rn) fe["alpha", , drop=FALSE] else NULL
    bRow <- fe[grep("^beta", rn), , drop=FALSE]
    sigEps  <- if("log_tauEps" %in% rn) exp(-0.5*fe["log_tauEps","Estimate"]) else NA
    sigBYM2 <- if("log_tau"    %in% rn) exp(-0.5*fe["log_tau","Estimate"])    else NA
    phiHat  <- if("logit_phi"  %in% rn) 1/(1+exp(-fe["logit_phi","Estimate"])) else NA
    list(alpha=aRow, beta=bRow, sigmaEps=sigEps, sigmaBYM2=sigBYM2, phi=phiHat,
         pdHess=res$TMBsd$pdHess)
}

parNames <- c("alpha", paste0("beta[", covs, "]"))
truthVec <- c(true_alpha, true_beta)

estTab <- data.frame(par=parNames, truth=round(truthVec, 3))
maxZ   <- list()
for(nm in c("M_D","M_M","M_DM")) {
    res <- list(M_D=resD, M_M=resM, M_DM=resDM)[[nm]]
    e   <- pullBYM2(res)
    if(is.null(e) || is.null(e$alpha)) {
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
    maxZ[[nm]] <- max(abs(z[-1]), na.rm=TRUE)
}

cat("\n----- BYM2 parameter recovery (truth vs estimate +/- SE; z = (est-truth)/SE) -----\n")
print(estTab, row.names=FALSE)

cat("\n----- BYM2 spatial / nugget hyperparameters -----\n")
for(nm in c("M_D","M_M","M_DM")) {
    res <- list(M_D=resD, M_M=resM, M_DM=resDM)[[nm]]
    e   <- pullBYM2(res)
    if(is.null(e)) next
    cat(sprintf("  %-5s  sigmaEps=%.3f  sigmaBYM2=%.3f  phi=%.3f  pdHess=%s\n",
                nm, e$sigmaEps, e$sigmaBYM2, e$phi, as.character(e$pdHess)))
}

cat("\n----- max |z| on covariate betas per BYM2 model (FE stage had 3.22 / 3.29 / 4.49) -----\n")
for(nm in names(maxZ)) cat(sprintf("  %-5s  max|z| = %s\n", nm, format(maxZ[[nm]], digits=3)))

cat("\nDone.\n")
