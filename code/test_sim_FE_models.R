# Cluster-friendly sanity test: run all four FE+nug+GH models on simulation
# replicate 1 and report estimates vs known truth.
#
# Models:  M_D (fitFED), M_M (fitFEM), M_DM (fitFEMD), M_DM_comm (fitFEMD_comm)
# Truth (simStudy1, per CLAUDE.md):
#   alpha = -1.25,  beta = (urban=1.0, access=0, elev=0, distRiversLakes=0, normPop=0.5)
#   sigmaEpsilon = sqrt(1.5)
#
# Run from code/:  Rscript test_sim_FE_models.R
# Wall time: ~10-20 min.

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/modFEMD_comm.R")
source("code/makeInputsTMB.R")

simIndex <- 1
Qgh   <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
covs  <- c("urban","access","elev","distRiversLakes","normPop")

# truth from simStudy1 DGP
true_alpha    <- -1.25
true_beta     <- c(urban=1.0, access=0, elev=0, distRiversLakes=0, normPop=0.5)
true_sigmaEps <- sqrt(1.5)

# ---- load simulated surveys ----
cat("Loading simulated surveys...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
datMICS <- surveysMICS[[simIndex]]
datDHS  <- surveysDHS[[simIndex]]

# normalize column names
nameTab <- rbind(c("N", "n"), c("N", "ns"), c("Z", "y"), c("Z", "ys"))
for(i in 1:nrow(nameTab)) {
    fromN <- nameTab[i,1]; toN <- nameTab[i,2]
    if(!(toN %in% names(datDHS)))  datDHS[[toN]]  <- datDHS[[fromN]]
    if(!(toN %in% names(datMICS))) datMICS[[toN]] <- datMICS[[fromN]]
}
if(!("Stratum" %in% names(datMICS))) datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

# ---- build / load inputs (one-shot) ----
inputsFile <- sprintf("savedOutput/simStudy1/inputs_sim%d_KMICS%03d_KDHS%02d_%02d.RData",
                      simIndex, KMICS, KDHSu, KDHSr)
if(file.exists(inputsFile)) {
    cat("Loading precomputed inputs:", inputsFile, "\n")
    load(inputsFile)
} else {
    cat("Building inputs (one-shot, not timed) ...\n")
    inputsMDM <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS, KMICS=KMICS,
                               KDHSurb=KDHSu, KDHSrur=KDHSr, saveNewIntPts=TRUE)
    save(inputsMDM, file=inputsFile)
}

# DHS column-name canonicalization (urb -> urban, pop -> normPop)
fixDhsNames <- function(mat) {
    if(is.null(mat)) return(mat)
    cn <- colnames(mat)
    if(!is.null(cn)) {
        cn <- gsub("^urb$", "urban", cn)
        cn <- gsub("^pop$", "normPop", cn)
        colnames(mat) <- cn
    }
    mat
}
inputsMDM$intPtsDHS$covsUrb <- fixDhsNames(inputsMDM$intPtsDHS$covsUrb)
inputsMDM$intPtsDHS$covsRur <- fixDhsNames(inputsMDM$intPtsDHS$covsRur)

# ---- helper: run one model and pull fixed-effect estimates ----
runOne <- function(label, fitExpr) {
    cat("\n========== ", label, " ==========\n", sep="")
    t0 <- proc.time()[3]
    res <- eval(fitExpr)
    elapsed <- proc.time()[3] - t0
    cat(sprintf("  fit time: %.2f min   convergence: %s   pdHess: %s\n",
                elapsed/60,
                ifelse(is.null(res$opt), "(MCMC)", as.character(res$opt$convergence)),
                if(inherits(res$TMBsd, "sdreport")) as.character(res$TMBsd$pdHess) else "n/a"))
    res
}

res_D <- runOne("M_D    (fitFED)",
                quote(fitFED(datDHS=datDHS, inputsMDM=inputsMDM,
                             covariates=covs, Qgh=Qgh,
                             fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))

res_M <- runOne("M_M    (fitFEM)",
                quote(fitFEM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                             covariates=covs, KMICS=KMICS, Qgh=Qgh,
                             fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))

res_DM <- runOne("M_DM   (fitFEMD)",
                 quote(fitFEMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                               covariates=covs, KMICS=KMICS, Qgh=Qgh,
                               fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))

res_C <- runOne("M_DM_comm (fitFEMD_comm)",
                quote(fitFEMD_comm(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                                   covariates=covs, KMICS=KMICS, Qgh=Qgh,
                                   fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)))

# ---- collect estimates ----
pickFixed <- function(res, alphaName="alpha", betaName="beta") {
    if(!inherits(res$TMBsd, "sdreport")) return(NULL)
    re <- summary(res$TMBsd, "random")
    if(nrow(re) > 0 && (alphaName %in% rownames(re) || betaName %in% rownames(re))) {
        # commensurate model: alpha/beta are random effects
        rn <- rownames(re)
        a  <- re[rn == alphaName, , drop=FALSE]
        b  <- re[rn == betaName,  , drop=FALSE]
        return(list(alpha=as.numeric(a[1,]), beta=b))
    }
    fe <- summary(res$TMBsd, "fixed")
    a  <- fe["alpha", , drop=FALSE]
    b  <- fe[grep("^beta", rownames(fe)), , drop=FALSE]
    list(alpha=as.numeric(a[1,]), beta=b)
}

estD  <- pickFixed(res_D)
estM  <- pickFixed(res_M)
estDM <- pickFixed(res_DM)
estC  <- pickFixed(res_C, alphaName="alpha", betaName="beta")  # DHS-side params

# log_tauEps -> sigmaEps
getSigmaEps <- function(res) {
    if(!inherits(res$TMBsd, "sdreport")) return(NA)
    fe <- summary(res$TMBsd, "fixed")
    if("log_tauEps" %in% rownames(fe)) exp(-0.5 * fe["log_tauEps", "Estimate"]) else NA
}
sig_D  <- getSigmaEps(res_D)
sig_M  <- getSigmaEps(res_M)
sig_DM <- getSigmaEps(res_DM)
sig_C  <- getSigmaEps(res_C)

# sigma_comm for the commensurate model (outer hyperparameter)
sigComm <- NA
if(inherits(res_C$TMBsd, "sdreport")) {
    fe <- summary(res_C$TMBsd, "fixed")
    if("log_sigma_comm" %in% rownames(fe)) sigComm <- exp(fe["log_sigma_comm","Estimate"])
}

# ---- print summary table ----
parNames <- c("alpha", paste0("beta[", covs, "]"))
truth    <- c(true_alpha, true_beta)

asRow <- function(est) {
    if(is.null(est)) return(rep(NA_real_, length(parNames)))
    c(est$alpha[1], est$beta[,"Estimate"])
}
seRow <- function(est) {
    if(is.null(est)) return(rep(NA_real_, length(parNames)))
    se_a <- if(length(est$alpha) >= 2) est$alpha[2] else est$alpha["Std. Error"]
    c(se_a, est$beta[,"Std. Error"])
}

tab <- data.frame(
    par   = parNames,
    truth = round(truth, 3),
    M_D_est  = round(asRow(estD),  3),  M_D_se  = round(seRow(estD),  3),
    M_M_est  = round(asRow(estM),  3),  M_M_se  = round(seRow(estM),  3),
    M_DM_est = round(asRow(estDM), 3),  M_DM_se = round(seRow(estDM), 3),
    Mcomm_est= round(asRow(estC),  3),  Mcomm_se= round(seRow(estC),  3)
)

cat("\n========== Parameter recovery on sim ", simIndex, " ==========\n", sep="")
print(tab, row.names=FALSE)

cat("\nNugget SD (truth = ", round(true_sigmaEps, 3), "):\n", sep="")
cat(sprintf("  M_D       sigmaEps = %.3f\n", sig_D))
cat(sprintf("  M_M       sigmaEps = %.3f\n", sig_M))
cat(sprintf("  M_DM      sigmaEps = %.3f\n", sig_DM))
cat(sprintf("  M_DM_comm sigmaEps = %.3f\n", sig_C))
cat(sprintf("  M_DM_comm sigma_comm = %.3f  (EB-selected)\n", sigComm))

cat("\nDone.\n")
