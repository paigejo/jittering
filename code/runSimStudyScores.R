# ============================================================================
# Score the 6 FE/BYM2 sim-study models on each of NSIM SPDE-simulated populations.
# For every (sim, model) pair, compute three sets of scores via getScores():
#   - fixed effects   against the known DGP truth  (TRUE_FE)
#   - area level      against true area prevalences      (areaPops[, sim])
#   - subarea level   against true subarea prevalences   (subareaPops[, sim])
# Saves savedOutput/simStudy1/scores/scores_<model>_sim<i>.RData per pair.
# ============================================================================

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/modM_DSep.R")    # fitMD  (BYM2)
source("code/modM_MSep.R")    # fitMM  (BYM2)
source("code/modM_DMSep.R")   # fitMDM (BYM2)
source("code/modMdSep.R")     # fitMd  (BYM2, observed-coords-as-truth)
source("code/makeInputsTMB.R")
source("code/modBYM2.R")      # predGrid / predArea

# ---- config ------------------------------------------------------------------
NSIM   <- 10
Qgh    <- 10
KMICS  <- 100
KDHSu  <- 16
KDHSr  <- 21
COVS   <- c("urban", "access", "elev", "distRiversLakes", "normPop")
NDRAWS <- 5000

# True DGP parameters baked into simData1().
TRUE_FE <- c(alpha                = -1.25,
             beta_urban           =  1.0,
             beta_access          =  0.0,
             beta_elev            =  0.0,
             beta_distRiversLakes =  0.0,
             beta_normPop         =  0.5)

# ---- load shared data --------------------------------------------------------
# simPopsSurveys.RData provides:
#   surveysDHS, surveysMICS  : per-sim survey data
#   subareaPops (775 x NSIM) : true subarea prevalences, rows = adm2 subareas
#   areaPops    (37  x NSIM) : true area    prevalences, rows = adm1 areas
load("savedOutput/simStudy1/simPopsSurveys.RData")
load("savedOutput/global/intPtsMICS_100.RData")   # MICS int pts, shared across sims

outDir <- "savedOutput/simStudy1/scores"
if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

# ============================================================================
# helpers
# ============================================================================

# Load sim `simIdx` and build the canonical TMB input list. Avoids GDAL by
# passing pre-built integration points directly into makeInputsMDM().
prepSim <- function(simIdx) {
    datDHS  <- surveysDHS[[simIdx]]
    datMICS <- surveysMICS[[simIdx]]

    # normalize survey column names (N/Z aliases for n/ns/y/ys)
    for(p in list(c("N","n"), c("N","ns"), c("Z","y"), c("Z","ys"))) {
        if(!(p[2] %in% names(datDHS)))  datDHS[[p[2]]]  <- datDHS[[p[1]]]
        if(!(p[2] %in% names(datMICS))) datMICS[[p[2]]] <- datMICS[[p[1]]]
    }
    if(!("Stratum" %in% names(datMICS)))
        datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

    # this sim's K=16/21 DHS integration points (populates `intPtsDHS`)
    load(sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K16_21.RData", simIdx))

    inp <- makeInputsMDM(datDHS=datDHS, datMICS=datMICS,
                         intPtsDHS=intPtsDHS, intPtsMICS=intPtsMICS,
                         KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                         saveNewIntPts=FALSE)

    # canonicalize DHS covariate names (urb -> urban, pop -> normPop)
    canon <- function(m) {
        cn <- colnames(m)
        if(!is.null(cn)) {
            cn <- gsub("^urb$", "urban",   cn)
            cn <- gsub("^pop$", "normPop", cn)
            colnames(m) <- cn
        }
        m
    }
    inp$intPtsDHS$covsUrb <- canon(inp$intPtsDHS$covsUrb)
    inp$intPtsDHS$covsRur <- canon(inp$intPtsDHS$covsRur)

    list(datDHS=datDHS, datMICS=datMICS, inputsMDM=inp)
}

# Draw NDRAWS posterior samples of all parameters (fixed + random) as a
# (nPar x NDRAWS) matrix with row-names = parameter names. Falls back to
# cov.fixed when jointPrecision isn't available (pure FE-only fits).
postDraws <- function(SD) {
    if(!is.null(SD$jointPrecision)) {
        mode  <- c(SD$par.fixed, SD$par.random)
        L     <- Matrix::Cholesky(SD$jointPrecision, LDL=FALSE, super=NA)
        z     <- matrix(rnorm(length(mode) * NDRAWS), nrow=length(mode))
        draws <- as.matrix(Matrix::solve(L, z, system="Lt")) + mode
    } else {
        mode  <- SD$par.fixed
        Lc    <- chol(SD$cov.fixed)
        z     <- matrix(rnorm(length(mode) * NDRAWS), nrow=length(mode))
        draws <- mode + t(Lc) %*% z
    }
    rownames(draws) <- names(mode)
    draws
}

# Score the fixed effects (alpha + 5 betas) against TRUE_FE.
# fitMd uses repar=TRUE: there is no explicit `alpha` parameter — the intercept
# is absorbed into w_bym2Star. In that case we recover alpha per posterior draw
# as the column-mean of the w_bym2Star draws (size nArea x NDRAWS).
scoreFE <- function(res) {
    SD <- res$TMBsd
    if(!inherits(SD, "sdreport")) return(NULL)
    draws <- tryCatch(postDraws(SD), error=function(e) NULL)
    if(is.null(draws)) return(NULL)

    aRow  <- draws[rownames(draws) == "alpha", , drop=FALSE]
    if(nrow(aRow) == 0) {
        wb <- draws[rownames(draws) == "w_bym2Star", , drop=FALSE]
        if(nrow(wb) > 0) aRow <- matrix(colMeans(wb), nrow=1)
    }
    bRows   <- draws[rownames(draws) == "beta", , drop=FALSE]
    feDraws <- rbind(aRow, bRows)
    if(nrow(feDraws) != length(TRUE_FE)) return(NULL)

    getScores(truth=TRUE_FE, estMat=feDraws,
              significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
              getAverage=FALSE, na.rm=TRUE)
}

# Md collapses positional uncertainty by zeroing all integration-point weights
# except the first column (the cluster's reported center). Reuses the existing
# inputsMDM but degenerates the K=16/21 weighting scheme to K=1-style.
makeMdInputs <- function(inputsMDM) {
    inp <- inputsMDM
    inp$intPtsDHS$wUrban[, 1]  <- 1; inp$intPtsDHS$wUrban[, -1] <- 0
    inp$intPtsDHS$wRural[, 1]  <- 1; inp$intPtsDHS$wRural[, -1] <- 0
    inp
}

# Score area-level (admin1) or subarea-level (adm2) predictions against the
# true population-weighted prevalences for sim simIdx.
scoreSpatial <- function(res, simIdx, level=c("area","subarea")) {
    level <- match.arg(level)
    if(!inherits(res$TMBsd, "sdreport")) return(NULL)

    grid <- tryCatch(
        predGrid(res$TMBsd, popMat=popMatNGAThresh, nsim=NDRAWS,
                 obj=res$TMBobj, admLevel="stratMICS"),
        error=function(e) NULL)
    if(is.null(grid)) return(NULL)

    if(level == "area") {
        agg   <- predArea(grid, areaVarName="area", orderedAreas=adm1@data$NAME_1)
        truth <- areaPops[, simIdx]
    } else {
        agg   <- predArea(grid, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
        truth <- subareaPops[, simIdx]
    }

    estMat <- agg$aggregationResults$p
    rownames(estMat) <- agg$aggregationResults$region
    # align estMat rows to truth (defensive — both should be admX$NAME_* order)
    estMat <- estMat[match(names(truth), rownames(estMat)), , drop=FALSE]

    getScores(truth=truth, estMat=estMat,
              significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
              getAverage=TRUE, na.rm=TRUE)
}

# Fit one of the named models on the prepped sim inputs.
#   Md / Md_FE: collapse to K=1 (observed cluster center) via makeMdInputs.
fitOne <- function(modName, inp) {
    switch(modName,
        Md_FE     = fitFED (datDHS=inp$datDHS, inputsMDM=makeMdInputs(inp$inputsMDM),
                            covariates=COVS, Qgh=Qgh, fixedEffectsOnly=TRUE,
                            getSDs=TRUE, verbose=FALSE),
        Md        = fitMd  (datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                            verbose=FALSE),
        M_D_FE    = fitFED (datDHS=inp$datDHS, inputsMDM=inp$inputsMDM,
                            covariates=COVS, Qgh=Qgh, fixedEffectsOnly=TRUE,
                            getSDs=TRUE, verbose=FALSE),
        M_M_FE    = fitFEM (datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            covariates=COVS, KMICS=KMICS, Qgh=Qgh, fixedEffectsOnly=TRUE,
                            getSDs=TRUE, verbose=FALSE),
        M_DM_FE   = fitFEMD(datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            covariates=COVS, KMICS=KMICS, Qgh=Qgh, fixedEffectsOnly=TRUE,
                            getSDs=TRUE, verbose=FALSE),
        M_D_BYM2  = fitMD  (datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                            Qgh=Qgh, getSDs=TRUE, verbose=FALSE),
        M_M_BYM2  = fitMM  (datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE),
        M_DM_BYM2 = fitMDM (datDHS=inp$datDHS, datMICS=inp$datMICS, inputsMDM=inp$inputsMDM,
                            KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                            Qgh=Qgh, getSDs=TRUE, verbose=FALSE))
}

MODELS <- c("Md_FE",  "Md",
            "M_D_FE", "M_M_FE",  "M_DM_FE",
            "M_D_BYM2", "M_M_BYM2", "M_DM_BYM2")

# ============================================================================
# main loop
# ============================================================================
for(simIdx in 1:NSIM) {
    cat("\n##### sim ", simIdx, " #####\n", sep="")
    inp <- prepSim(simIdx)

    for(modName in MODELS) {
        outFile <- sprintf("%s/scores_%s_sim%d.RData", outDir, modName, simIdx)
        if(file.exists(outFile)) { cat("  [", modName, "] cached, skipping\n", sep=""); next }

        cat("  [", modName, "] fitting ... ", sep="")
        t0  <- proc.time()[3]
        res <- tryCatch(fitOne(modName, inp), error=function(e) NULL)
        if(is.null(res)) { cat("FAILED to fit\n"); next }
        cat(sprintf("fit %.1f min\n", (proc.time()[3] - t0)/60))

        scoresFE      <- scoreFE     (res)
        scoresArea    <- scoreSpatial(res, simIdx, "area")
        scoresSubarea <- scoreSpatial(res, simIdx, "subarea")

        save(scoresFE, scoresArea, scoresSubarea, file=outFile)
    }
}

cat("\nAll done.\n")
