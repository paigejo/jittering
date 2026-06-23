# ============================================================================
# Application to women's secondary education (real 2018 NDHS + 2016 NMICS).
# Fits the three paper approaches M_d, M_D, M_DM (BYM2, NO commensurate prior)
# and produces, reusing the SAME functions as the simulation study:
#   - parameter table (Int, Urban, Healthcare(access), Elevation, Dist R&L,
#     Population, sigma^2, phi, sigma_eps^2) with posterior mean + 95% CI
#     -> matches Table tab:parSummary1
#   - pixel / Admin2 / Admin1 predictions: posterior mean + 95% CI width
# All raster-derived inputs are precomputed (intPtsDHS_16_21, intPtsMICS_100,
# predGridInputs), so no GDAL is needed at fit time.
#
# Reused machinery (no reimplementation):
#   buildFullMDMInputs()  -> canonical TMB inputs via makeInputsMDM
#   fitMd / fitMD / fitMDM -> the BYM2 fits used in the sim study
#   predGrid() / predArea() -> posterior param draws + grid/areal predictions
# ============================================================================
# path-robust sourcing: run from project root (local) or code/ (cluster)
if(file.exists("code/setup.R")) { source("code/setup.R"); source("code/validationCommon.R") } else { source("setup.R"); source("validationCommon.R") }
options(error = traceback)

KMICS <- 100; KDHSu <- 16; KDHSr <- 21; Qgh <- 10; NDRAWS <- 1000
covs  <- c("urban","access","elev","distRiversLakes","normPop")
outDir <- "savedOutput/application"; if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
# Model labels follow the repo convention: "Md" (no underscore) vs "M_D" so they
# differ by more than case -> safe on the case-INSENSITIVE Windows filesystem.
APPROACHES <- c("Md","M_D","M_DM")

# ---- build canonical full-data inputs ONCE (uses precomputed integration pts) ----
cat(sprintf("[%s] building full MDM inputs ...\n", format(Sys.time()))); flush.console()
fullInp <- buildFullMDMInputs(KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr)

# ---- parameter table rows + helper to summarize predGrid's parameter draws ----
parRowNames <- c("(Int)","Urban","Healthcare Inaccessibility","Elevation",
                 "Dist. Rivers & Lakes","Population",
                 "sigma^2","phi","sigma_eps^2")
sumDraws <- function(v) c(Est=mean(v,na.rm=TRUE), L=unname(quantile(v,.025,na.rm=TRUE)), U=unname(quantile(v,.975,na.rm=TRUE)))
ciW <- function(x) diff(quantile(x, c(.025,.975), na.rm=TRUE))
parTableFromGrid <- function(grid){
  betaD <- grid$betaDraws                      # 5 x nsim, order = covs
  M <- rbind(grid$alphaDraws,
             betaD,
             grid$sigmaSqDraws,
             grid$phiDraws,
             grid$sigmaEpsSqDraws)
  tab <- t(apply(M, 1, sumDraws))
  rownames(tab) <- parRowNames
  round(tab, 3)
}

# pixel / area predictions: posterior mean + 95% CI width (a2/a1 precomputed)
predLevels <- function(grid, a2, a1){
  nNaN <- sum(is.na(grid$gridDraws))
  if(nNaN) cat(sprintf("    note: %d NaN among %d pixel-draws (na.rm in summaries)\n",
                       nNaN, length(grid$gridDraws)))
  pixMean  <- rowMeans(grid$gridDraws, na.rm=TRUE)
  pixWidth <- apply(grid$gridDraws, 1, ciW)
  a2p <- as.matrix(a2$aggregationResults$p); a1p <- as.matrix(a1$aggregationResults$p)
  list(pixel   = data.frame(mean=pixMean, ciWidth=pixWidth),
       admin2  = data.frame(region=a2$aggregationResults$region,
                            mean=rowMeans(a2p,na.rm=TRUE), ciWidth=apply(a2p,1,ciW)),
       admin1  = data.frame(region=a1$aggregationResults$region,
                            mean=rowMeans(a1p,na.rm=TRUE), ciWidth=apply(a1p,1,ciW)))
}

fitOne <- function(tag){
  fitFile <- sprintf("%s/fit_%s.RData", outDir, tag)
  if(file.exists(fitFile)) {                                # reuse checkpointed fit -> re-predict only
    cat(sprintf("\n[%s] loading cached fit %s ...\n", format(Sys.time()), tag)); flush.console()
    e <- new.env(); load(fitFile, envir=e); res <- e$res
  } else {
    cat(sprintf("\n[%s] fitting %s ...\n", format(Sys.time()), tag)); flush.console()
    t0 <- proc.time()[3]
    res <- switch(tag,
      Md   = fitMd (datDHS=ed, datMICS=edMICS, inputsMDM=fullInp, KMICS=KMICS,
                    KDHSurb=KDHSu, KDHSrur=KDHSr, verbose=FALSE),
      M_D  = fitMD (datDHS=ed, datMICS=edMICS, inputsMDM=fullInp, KMICS=KMICS,
                    KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE),
      M_DM = fitMDM(datDHS=ed, datMICS=edMICS, inputsMDM=fullInp, KMICS=KMICS,
                    KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE))
    cat(sprintf("  fit %.1f min\n", (proc.time()[3]-t0)/60)); flush.console()
    save(res, file=fitFile)   # checkpoint fit so re-prediction never re-fits
    cat(sprintf("  checkpointed fit -> %s\n", fitFile)); flush.console()
  }
  cat(sprintf("  [%s] predicting ...\n", tag)); flush.console()
  # Gaussian posterior draws (useInla=FALSE): established equivalent to INLA-style, faster.
  grid <- predGrid(res$TMBsd, popMat=popMatNGAThresh, nsim=NDRAWS,
                   obj=res$TMBobj, admLevel="stratMICS", useInla=FALSE)
  a2 <- predArea(grid, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
  a1 <- predArea(grid, areaVarName="area",    orderedAreas=adm1@data$NAME_1)
  # also save in the format plotAppPreds() expects (savedOutput/ed/<tag>SepRepar.RData)
  edDir <- "savedOutput/ed"; if(!dir.exists(edDir)) dir.create(edDir, recursive=TRUE)
  gridPreds <- grid;        save(gridPreds,   file=sprintf("%s/gridPreds%sSepRepar.RData",   edDir, tag))
  admin2Preds <- a2;        save(admin2Preds, file=sprintf("%s/admin2Preds%sSepRepar.RData", edDir, tag))
  admin1Preds <- a1;        save(admin1Preds, file=sprintf("%s/admin1Preds%sSepRepar.RData", edDir, tag))
  parTab <- parTableFromGrid(grid)
  preds  <- predLevels(grid, a2, a1)
  cat(sprintf("  %s parameter table:\n", tag)); print(parTab)
  save(parTab, preds, file=sprintf("%s/application_%s.RData", outDir, tag))
  list(parTab=parTab, preds=preds)
}

results <- list()
for(tag in APPROACHES)
  results[[tag]] <- tryCatch(fitOne(tag),
                             error=function(e){cat("  ERROR:",conditionMessage(e),"\n"); NULL})

# ---- combined parameter table (paper layout: M_DM | M_D | Md) ----
combined <- do.call(cbind, lapply(c("M_DM","M_D","Md"), function(t){
  if(is.null(results[[t]])) return(NULL)
  m <- results[[t]]$parTab
  out <- data.frame(Est=m[,"Est"], CI=sprintf("(%.2f, %.2f)", m[,"L"], m[,"U"]))
  names(out) <- paste0(t, c("_Est","_CI")); out
}))
rownames(combined) <- parRowNames
cat("\n================ Application parameter table (vs paper tab:parSummary1) ================\n")
print(combined)
save(results, combined, file=sprintf("%s/applicationSummary.RData", outDir))

# ---- figures: pixel / admin2 / admin1 maps (mean + 95% CI width) ----
cat(sprintf("\n[%s] generating application maps via plotAppPreds() ...\n", format(Sys.time())))
if(!exists("plotAppPreds"))
  source(if(file.exists("code/plotGenerator.R")) "code/plotGenerator.R" else "plotGenerator.R")
if(!dir.exists("figures/ed")) dir.create("figures/ed", recursive=TRUE)
tryCatch(plotAppPreds(adm2Model=FALSE),
         error=function(e) cat("  plotAppPreds error:", conditionMessage(e), "\n"))

cat(sprintf("\n[%s] application done; saved to %s/ ; figures in figures/ed/\n", format(Sys.time()), outDir))
