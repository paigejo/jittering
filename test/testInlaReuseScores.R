# DECISIVE check: do the DOWNSTREAM SCORES (what we actually report) agree
# between the tape-reuse INLA walk and the rebuild-each-step walk? The raw
# w_bym2Free draws differ ~9% in the unidentified spatial direction, but the
# spatial effect is sum-to-zero and population-aggregated to areas, so the
# areal scores should be far more stable. If scoresArea/FE/Hyper agree to
# ~1%, the 7x-faster reuse is safe to adopt.
#
# We FORCE the INLA path everywhere (M_D_BYM2's jointPrecision is PD, so
# useInla="auto" would otherwise pick the Gaussian path and not exercise the
# code we're comparing): override .postDraws to force useInla=TRUE (drives
# .scoreFE/.scoreHyper) and call predGrid(useInla=TRUE) for the area score.

source("code/setup.R")          # NEW reuse inlaStyleDraws
options(warn = 1)

SEED   <- 31415
SIMIDX <- 17
simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)
truths   <- modelTruths("bym2")
areaPops <- simEnv$areaPops_smoothRisk

cat(sprintf("\n=== %s | fitting M_D_BYM2 ===\n", format(Sys.time())))
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))

# Force INLA in the FE/Hyper path (override the helper used by .scoreFE/.scoreHyper)
.postDraws <<- function(res, NDRAWS, useInla = "auto")
    tryCatch(posteriorDraws(res, NDRAWS = NDRAWS, useInla = TRUE),
             error = function(e) NULL)

# Area scorer with INLA forced (mirrors .scoreArea but useInla=TRUE)
scoreAreaForceInla <- function(res, areaPops, simIdx, NDRAWS) {
    grid <- predGrid(res$TMBsd, popMat = popMatNGAThresh, nsim = NDRAWS,
                     obj = res$TMBobj, admLevel = "stratMICS", res = res,
                     useInla = TRUE)
    agg  <- predArea(grid, areaVarName = "area", orderedAreas = adm1@data$NAME_1)
    truth  <- areaPops[, simIdx]
    estMat <- as.matrix(agg$aggregationResults$p)
    rownames(estMat) <- agg$aggregationResults$region
    estMat <- estMat[match(names(truth), rownames(estMat)), , drop = FALSE]
    getScores(truth = truth, estMat = estMat, significance = c(.5,.8,.9,.95),
              doFuzzyReject = TRUE, getAverage = TRUE, na.rm = TRUE)
}

runScores <- function(label) {
    cat(sprintf("\n=== %s | scoring (%s) ===\n", format(Sys.time()), label))
    set.seed(SEED); fe    <- .scoreFE(res, 1000)
    set.seed(SEED); hyper <- .scoreHyper(res, truths, 1000)
    set.seed(SEED); area  <- scoreAreaForceInla(res, areaPops, SIMIDX, 1000)
    list(FE = fe, Hyper = hyper, Area = area)
}

sc_new <- runScores("NEW reuse")

cat(sprintf("\n=== %s | sourcing OLD rebuild inlaStyleDraws ===\n", format(Sys.time())))
source("code/inlaStyleDraws_OLD.R")
.postDraws <<- function(res, NDRAWS, useInla = "auto")   # re-assert override
    tryCatch(posteriorDraws(res, NDRAWS = NDRAWS, useInla = TRUE),
             error = function(e) NULL)
sc_old <- runScores("OLD rebuild")

cmp <- function(a, b, nm) {
    if(is.null(a) || is.null(b)) { cat(sprintf("  %-11s NULL (new=%s old=%s)\n",
                                               nm, !is.null(a), !is.null(b))); return(invisible()) }
    A <- as.matrix(a); B <- as.matrix(b)
    if(!all(dim(A) == dim(B))) { cat(sprintf("  %-11s DIM MISMATCH\n", nm)); return(invisible()) }
    amax <- max(abs(A - B), na.rm = TRUE)
    rel  <- max(abs(A - B) / (abs(B) + 1e-8), na.rm = TRUE)
    cat(sprintf("  %-11s max abs diff %.3e   max rel diff %.3e   %s\n",
                nm, amax, rel, ifelse(rel < 1e-2, "PASS(<1%)",
                                      ifelse(rel < 5e-2, "OK(<5%)", "INVESTIGATE"))))
}
cat("\n================ DOWNSTREAM SCORES: reuse vs rebuild ================\n")
cmp(sc_new$FE,    sc_old$FE,    "scoresFE")
cmp(sc_new$Hyper, sc_old$Hyper, "scoresHyper")
cmp(sc_new$Area,  sc_old$Area,  "scoresArea")

cat("\n--- scoresArea side by side (new / old) ---\n")
if(!is.null(sc_new$Area) && !is.null(sc_old$Area)) {
    A <- as.matrix(sc_new$Area); B <- as.matrix(sc_old$Area)
    show <- intersect(colnames(A), c("CRPS","MSE","Coverage50","Coverage80","Coverage90","Coverage95","IntervalScore90"))
    for(cc in show) cat(sprintf("  %-16s %.5f / %.5f\n", cc, A[1, cc], B[1, cc]))
}
cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
