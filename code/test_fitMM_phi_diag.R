# Local diagnostic for fitMM phi -> 1 boundary pathology on BYM2-simulated data.
# Three configurations compared, plus NLL profile in phi for each.
#
#   r1: fitMM, default PC prior on phi
#   r2: fitMM, uniform prior on phi
#   r3: fitMM, MICS weights collapsed to K=1 (all mass on int pt 1)
#
# After each fit prints: (logit_phi, phi, sigma_BYM2, sigma_eps, pdHess, NLL).
# Then profiles NLL at phi = 0.5..0.999 holding other params at that fit's MLE.

source("setup.R")
options(error = recover)
source("code/modM_MSep.R")
source("code/makeInputsTMB.R")

# --- config (match the runSimStudyScores.R defaults) ---
simIdx <- 1
KMICS  <- 100
KDHSu  <- 16
KDHSr  <- 21
Qgh    <- 10

# --- load data + build inputsMDM for sim simIdx ---
load("savedOutput/simStudy1/simPopsSurveys.RData")     # BYM2-simulated (10 sims)
load("savedOutput/global/intPtsMICS_100.RData")        # MICS int pts
load(sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K16_21.RData", simIdx))

datDHS  <- surveysDHS [[simIdx]]
datMICS <- surveysMICS[[simIdx]]
for(p in list(c("N","n"), c("N","ns"), c("Z","y"), c("Z","ys"))) {
    if(!(p[2] %in% names(datDHS)))  datDHS[[p[2]]]  <- datDHS[[p[1]]]
    if(!(p[2] %in% names(datMICS))) datMICS[[p[2]]] <- datMICS[[p[1]]]
}
if(!("Stratum" %in% names(datMICS))) datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

inp <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS,
                     intPtsDHS=intPtsDHS, intPtsMICS=intPtsMICS,
                     KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                     saveNewIntPts=FALSE)
# Canonicalize DHS covariate names (urb -> urban, pop -> normPop)
canon <- function(m) {
    cn <- colnames(m)
    if(!is.null(cn)) { cn <- gsub("^urb$","urban",cn); cn <- gsub("^pop$","normPop",cn); colnames(m) <- cn }
    m
}
inp$intPtsDHS$covsUrb <- canon(inp$intPtsDHS$covsUrb)
inp$intPtsDHS$covsRur <- canon(inp$intPtsDHS$covsRur)

# --- helper: pretty-print one fit's summary ---
summFit <- function(label, res) {
    p  <- res$opt$par
    cat(sprintf("\n%-44s   logit_phi=%6.3f  phi=%.4f  sigmaBYM2=%.3f  sigmaEps=%.3f  pdHess=%s  NLL=%.4f\n",
        label, p["logit_phi"], plogis(p["logit_phi"]),
        exp(-0.5*p["log_tau"]), exp(-0.5*p["log_tauEps"]),
        as.character(res$TMBsd$pdHess), res$opt$objective))
}

# --- helper: NLL profile in phi at fit MLE ---
profilePhi <- function(label, res) {
    obj <- res$TMBobj; mle <- res$opt$par
    phiVals <- c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 0.999)
    nlls <- sapply(phiVals, function(p) { par2 <- mle; par2["logit_phi"] <- qlogis(p); obj$fn(par2) })
    cat("---- NLL profile:", label, "----\n")
    prof <- data.frame(phi=phiVals, logit_phi=round(qlogis(phiVals),2),
                       NLL=round(nlls,4), diffNLL=round(nlls - min(nlls),4))
    print(prof, row.names=FALSE)
}

# --- run 1: default PC prior ---
cat("\n========== Run 1: PC prior on phi ==========\n")
t0 <- proc.time()[3]
r1 <- fitMM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inp,
            KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
cat(sprintf("fit time: %.1f min\n", (proc.time()[3]-t0)/60))
summFit("PC prior", r1); profilePhi("PC prior", r1)

# --- run 2: uniform prior on phi ---
cat("\n========== Run 2: uniform prior on phi ==========\n")
t0 <- proc.time()[3]
r2 <- fitMM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inp,
            KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE,
            uniform_phi_prior=TRUE)
cat(sprintf("fit time: %.1f min\n", (proc.time()[3]-t0)/60))
summFit("uniform prior", r2); profilePhi("uniform prior", r2)

# --- run 3: MICS weights collapsed to K=1 (degenerate MICS integration) ---
cat("\n========== Run 3: degenerate MICS weights (K=1, mass on int pt 1) ==========\n")
inpMd_M <- inp
inpMd_M$intPtsMICS$wUrban[, 1] <- 1; inpMd_M$intPtsMICS$wUrban[, -1] <- 0
inpMd_M$intPtsMICS$wRural[, 1] <- 1; inpMd_M$intPtsMICS$wRural[, -1] <- 0
t0 <- proc.time()[3]
r3 <- fitMM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inpMd_M,
            KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
cat(sprintf("fit time: %.1f min\n", (proc.time()[3]-t0)/60))
summFit("MICS K=1", r3); profilePhi("MICS K=1", r3)

cat("\nDone.\n")
