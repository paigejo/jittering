#!/usr/bin/env Rscript
# Compare FE+nug models: Combined (MICS+DHS), MICS-only, DHS-only on simulated data
# - Uses simulated surveys in savedOutput/simStudy1/simPopsSurveys_BYM2.RData
# - Uses KMICS=100, KDHSurb=16, KDHSrur=21, Qgh=10
# - Precomputes and saves inputs, compiles DLLs, then times only optimization + SE computation

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/makeInputsTMB.R")
source("code/modFEM.R")
source("code/modFED.R")
source("code/modFEMD.R")

Qgh <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
covs <- c("urban","access","elev","distRiversLakes","normPop")
simIndex <- 1

# Load simulated surveys
cat("Loading simulated surveys...\n")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
datMICS <- surveysMICS[[simIndex]]
datDHS <- surveysDHS[[simIndex]]

# Ensure expected column names (N->n/ns, Z->y/ys) and Stratum for MICS
nameTab <- rbind(c("N", "n"), c("N", "ns"), c("Z", "y"), c("Z", "ys"))
for(i in 1:nrow(nameTab)) {
  fromN <- nameTab[i,1]; toN <- nameTab[i,2]
  if(!(toN %in% names(datDHS)))     datDHS[[toN]]     <- datDHS[[fromN]]
  if(!(toN %in% names(datMICS))) datMICS[[toN]] <- datMICS[[fromN]]
}
if(!("Stratum" %in% names(datMICS))) datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

# Build / load inputs for this simulation (not timed)
inputsFile <- sprintf("savedOutput/simStudy1/inputs_sim%d_KMICS%03d_KDHS%02d_%02d.RData", simIndex, KMICS, KDHSu, KDHSr)
if(file.exists(inputsFile)) {
  cat("Loading precomputed inputs:", inputsFile, "\n")
  load(inputsFile) # loads inputsMDM
} else {
  cat("Creating inputs (this may take time; not included in model timings)\n")
  inputsMDM <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS, KMICS=KMICS,
                             KDHSurb=KDHSu, KDHSrur=KDHSr, saveNewIntPts=TRUE)
  save(inputsMDM, file=inputsFile)
  cat("Saved inputs to:", inputsFile, "\n")
}

# canonicalize DHS cov names (urb->urban, pop->normPop) if present as dimnames
fix_dhs_names <- function(mat) {
  if(is.null(mat)) return(mat)
  cn <- colnames(mat)
  if(is.null(cn)) {
    dn <- dimnames(mat)
    if(!is.null(dn) && length(dn) >= 2 && !is.null(dn[[2]])) cn <- dn[[2]]
  }
  if(!is.null(cn)) {
    cn <- gsub("^urb$", "urban", cn)
    cn <- gsub("^pop$", "normPop", cn)
    colnames(mat) <- cn
  }
  mat
}
inputsMDM$intPtsDHS$covsUrb <- fix_dhs_names(inputsMDM$intPtsDHS$covsUrb)
inputsMDM$intPtsDHS$covsRur <- fix_dhs_names(inputsMDM$intPtsDHS$covsRur)

# 1) Combined MICS + DHS FE+nug (map spatial terms out) -----------------
cat("\n=== Combined MICS+DHS: FE+nug (map spatial out) ===\n")
res_comb <- fitFEMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                    covariates=covs, KMICS=KMICS, Qgh=Qgh,
                    fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

pe_comb <- summary(res_comb$TMBsd, 'fixed')
alpha_comb <- pe_comb['alpha','Estimate']; alpha_comb_se <- pe_comb['alpha','Std. Error']
beta_comb <- pe_comb[grep('^beta', rownames(pe_comb)),'Estimate']
beta_comb_se <- pe_comb[grep('^beta', rownames(pe_comb)),'Std. Error']

# 2) MICS-only using fitFEM -----------------
cat('\n=== MICS-only: FE+nug via fitFEM ===\n')
res_mics <- fitFEM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                   covariates=covs, KMICS=KMICS, Qgh=Qgh,
                   fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

peM <- summary(res_mics$TMBsd, 'fixed')
mics_alpha <- peM['alpha','Estimate']; mics_alpha_se <- peM['alpha','Std. Error']
mics_beta <- peM[grep('^beta', rownames(peM)),'Estimate']
mics_beta_se <- peM[grep('^beta', rownames(peM)),'Std. Error']

# 3) DHS-only FE+nug using fitFED ---------------------------------
cat('\n=== DHS-only: FE+nug via fitFED ===\n')
res_dhs <- fitFED(datDHS=datDHS, inputsMDM=inputsMDM,
                  covariates=covs, Qgh=Qgh,
                  fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

peD <- summary(res_dhs$TMBsd, 'fixed')
alpha_dhs <- peD['alpha','Estimate']; alpha_dhs_se <- peD['alpha','Std. Error']
beta_dhs <- peD[grep('^beta', rownames(peD)),'Estimate']
beta_dhs_se <- peD[grep('^beta', rownames(peD)),'Std. Error']

# Summarize results into table
param_names <- c('alpha', paste0('beta[', covs, ']'))

est_comb_v <- c(alpha_comb, beta_comb)
se_comb_v <- c(alpha_comb_se, beta_comb_se)

est_mics <- c(mics_alpha, mics_beta)
se_mics_v <- c(mics_alpha_se, mics_beta_se)

est_dhs <- c(alpha_dhs, beta_dhs)
se_dhs_v <- c(alpha_dhs_se, beta_dhs_se)

res_df <- data.frame(Parameter=param_names,
                     Combined_Est=round(est_comb_v,4), Combined_SE=round(se_comb_v,4),
                     MICS_Est=round(est_mics,4), MICS_SE=round(se_mics_v,4),
                     DHS_Est=round(est_dhs,4), DHS_SE=round(se_dhs_v,4),
                     stringsAsFactors=FALSE)

cat('\n=== Parameter estimates and SEs ===\n')
print(res_df, row.names=FALSE)

cat('\n=== Timing (s) ===\n')
cat(sprintf('Combined: opt=%.2f  SE=%.2f  total=%.2f\n', res_comb$totalTime, res_comb$sdTime, res_comb$totalTime+res_comb$sdTime))
cat(sprintf('MICS-only: opt=%.2f  SE=%.2f  total=%.2f\n', res_mics$totalTime, res_mics$sdTime, res_mics$totalTime+res_mics$sdTime))
cat(sprintf('DHS-only: opt=%.2f  SE=%.2f  total=%.2f\n', res_dhs$totalTime, res_dhs$sdTime, res_dhs$totalTime+res_dhs$sdTime))

cat('\nDone.\n')
