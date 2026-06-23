# ============================================================================
# Commensurate-prior M_DM on the real application data, compared to the standard
# M_DM. Same data/inputs (makeInputsMDM -> inherits the area-mapping fix) and the
# same BYM2 core as modMDM_BYM2_GH_v2; the only addition is separate MICS effects
# (alpha_M, beta_M) with a commensurate prior (beta_M - beta) ~ N(0, sigma_comm^2).
# Predictions/parameters reported are the DHS-side alpha/beta (what predGrid pulls
# by name), matching the paper convention.
#
# Output: savedOutput/application/application_M_DM_comm.RData (+ fit checkpoint),
#         and a side-by-side comparison vs the standard M_DM parameter table.
# ============================================================================
if(file.exists("code/setup.R")) { source("code/setup.R"); source("code/validationCommon.R"); source("code/modFEMD_comm.R") } else { source("setup.R"); source("validationCommon.R"); source("modFEMD_comm.R") }
options(error = traceback)

KMICS <- 100; KDHSu <- 16; KDHSr <- 21; Qgh <- 10; NDRAWS <- 1000
covs  <- c("urban","access","elev","distRiversLakes","normPop")
outDir <- "savedOutput/application"; if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
fitFile <- sprintf("%s/fit_M_DM_comm.RData", outDir)

parRowNames <- c("(Int)","Urban","Healthcare Inaccessibility","Elevation",
                 "Dist. Rivers & Lakes","Population","sigma^2","phi","sigma_eps^2")
sumDraws <- function(v) c(Est=mean(v,na.rm=TRUE), L=unname(quantile(v,.025,na.rm=TRUE)), U=unname(quantile(v,.975,na.rm=TRUE)))
parTableFromGrid <- function(grid){
  M <- rbind(grid$alphaDraws, grid$betaDraws, grid$sigmaSqDraws, grid$phiDraws, grid$sigmaEpsSqDraws)
  tab <- t(apply(M, 1, sumDraws)); rownames(tab) <- parRowNames; round(tab, 3)
}

# ---- fit commensurate M_DM (BYM2 enabled), reusing the standard inputs ----
if(file.exists(fitFile)) {
  cat(sprintf("[%s] loading cached commensurate fit\n", format(Sys.time())))
  e <- new.env(); load(fitFile, envir=e); res <- e$res
} else {
  cat(sprintf("[%s] fitting commensurate M_DM (fixedEffectsOnly=FALSE) ...\n", format(Sys.time()))); flush.console()
  t0 <- proc.time()[3]
  res <- fitFEMD_comm(datDHS=ed, datMICS=edMICS, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                      Qgh=Qgh, fixedEffectsOnly=FALSE, getSDs=TRUE, verbose=FALSE)
  cat(sprintf("  fit %.1f min\n", (proc.time()[3]-t0)/60)); flush.console()
  save(res, file=fitFile)
}

# ---- sigma_comm (the key commensurate hyperparameter) ----
ssum <- tryCatch(summary(res$TMBsd, "fixed"), error=function(e) NULL)
sigmaCommStr <- "NA"
if(!is.null(ssum) && "log_sigma_comm" %in% rownames(ssum)) {
  lsc <- ssum["log_sigma_comm","Estimate"]; lscSE <- ssum["log_sigma_comm","Std. Error"]
  sc <- exp(lsc); scL <- exp(lsc-1.96*lscSE); scU <- exp(lsc+1.96*lscSE)
  sigmaCommStr <- sprintf("%.3f (%.3f, %.3f)", sc, scL, scU)
}
cat(sprintf("\nsigma_comm (commensurability; large => MICS & DHS effects free to differ) = %s\n", sigmaCommStr))

# ---- predictions + parameter table (DHS-side), reuse predGrid/predArea ----
cat(sprintf("[%s] predicting (Gaussian draws) ...\n", format(Sys.time()))); flush.console()
grid <- predGrid(res$TMBsd, popMat=popMatNGAThresh, nsim=NDRAWS, obj=res$TMBobj,
                 admLevel="stratMICS", useInla=FALSE)
parTab <- parTableFromGrid(grid)
a2 <- predArea(grid, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
a1 <- predArea(grid, areaVarName="area",    orderedAreas=adm1@data$NAME_1)
ciW <- function(x) diff(quantile(x, c(.025,.975), na.rm=TRUE))
preds <- list(pixel  = data.frame(mean=rowMeans(grid$gridDraws,na.rm=TRUE),
                                  ciWidth=apply(grid$gridDraws,1,ciW)),
              admin2 = data.frame(region=a2$aggregationResults$region,
                                  mean=rowMeans(as.matrix(a2$aggregationResults$p),na.rm=TRUE),
                                  ciWidth=apply(as.matrix(a2$aggregationResults$p),1,ciW)),
              admin1 = data.frame(region=a1$aggregationResults$region,
                                  mean=rowMeans(as.matrix(a1$aggregationResults$p),na.rm=TRUE),
                                  ciWidth=apply(as.matrix(a1$aggregationResults$p),1,ciW)))
save(parTab, preds, sigmaCommStr, file=sprintf("%s/application_M_DM_comm.RData", outDir))

# ---- compare to standard M_DM ----
cmp <- data.frame(comm_Est=parTab[,"Est"],
                  comm_CI=sprintf("(%.2f, %.2f)", parTab[,"L"], parTab[,"U"]))
stdFile <- sprintf("%s/application_M_DM.RData", outDir)
if(file.exists(stdFile)) {
  s <- new.env(); load(stdFile, envir=s); std <- s$parTab
  cmp$std_Est <- std[,"Est"]; cmp$std_CI <- sprintf("(%.2f, %.2f)", std[,"L"], std[,"U"])
}
rownames(cmp) <- parRowNames
cat("\n========= Commensurate M_DM vs standard M_DM (DHS-side params) =========\n")
print(cmp)
cat(sprintf("\nsigma_comm = %s\n", sigmaCommStr))
cat(sprintf("[%s] done.\n", format(Sys.time())))
