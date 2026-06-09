# Cross-machine reproducibility test: do (sim, model) pairs that the CLUSTER
# completed reproduce on the laptop with the current (synced) data?
#
# Usage:
#   - Place cluster scores at savedOutput/simStudy1/scores_full/{BYM2,SPDE}/
#   - Run this from the project root: Rscript --vanilla code/crossMachineTest.R
#   - It picks a CLUSTER-COMPLETED sim ID (defaults to scanning sims 1..100 for
#     M_D_BYM2 and grabbing the first one with a cluster score file), refits
#     it locally with current code, then compares scoresFE/scoresHyper.
#
# Expected outcome if data + code match: max|delta| < ~0.05 on scoresFE/Hyper.
# If max|delta| > 0.5, machines are NOT producing the same fit, which would
# invalidate all laptop-side debugging.

source("code/setup.R")
options(error = recover)

generative   <- "BYM2"    # which sim DGP
model        <- "M_D_BYM2"  # which scoring model — start with the failure-prone one
candidateIds <- 1:100

clusterDir <- sprintf("savedOutput/simStudy1/scores_full/%s", generative)
chosenSim  <- NA
for(i in candidateIds) {
    p <- sprintf("%s/scores_%s_sim%d.RData", clusterDir, model, i)
    if(file.exists(p)) { chosenSim <- i; break }
}
if(is.na(chosenSim))
    stop(sprintf("No cluster score file found for %s in %s for sims %d..%d",
                 model, clusterDir, min(candidateIds), max(candidateIds)))

cat(sprintf("\n[xmach] Using sim %d, model %s, generative %s\n",
            chosenSim, model, generative))
cat(sprintf("[xmach] Cluster file: %s/scores_%s_sim%d.RData\n",
            clusterDir, model, chosenSim))

# Re-fit on laptop using current code
simEnv <- new.env()
simulateSurveys(tolower(generative), nsim = chosenSim, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

areaPops_smooth <- simEnv$areaPops_smoothRisk
areaPops_prev   <- simEnv$areaPops
truths <- modelTruths(tolower(generative))

outDir <- "savedOutput/simStudy1/scores_xmach"
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
outFile <- sprintf("%s/scores_%s_sim%d.RData", outDir, model, chosenSim)

cat(sprintf("\n[xmach] Re-scoring locally (cluster wallclock ~few min for %s) ...\n", model))
.scoreOneFit(simIdx  = chosenSim,
             modName = model,
             outFile = outFile,
             surveysDHS  = simEnv$surveysDHS,
             surveysMICS = simEnv$surveysMICS,
             intPtsMICS  = micsEnv$intPtsMICS,
             model       = tolower(generative),
             truths      = truths,
             areaPops    = areaPops_smooth,
             KMICS = 100, KDHSu = 16, KDHSr = 21,
             Qgh = 10, NDRAWS = 1000,
             COVS = c("urban","access","elev","distRiversLakes","normPop"))

# Compare
clustE <- new.env(); load(sprintf("%s/scores_%s_sim%d.RData", clusterDir, model, chosenSim), envir = clustE)
locE   <- new.env(); load(outFile, envir = locE)

cmp <- function(field) {
    cat(sprintf("\n----- %s -----\n", field))
    cM <- if(exists(field, envir = clustE)) as.matrix(get(field, envir = clustE)) else NULL
    lM <- if(exists(field, envir = locE))   as.matrix(get(field, envir = locE))   else NULL
    if(is.null(cM) || is.null(lM)) {
        cat(sprintf("  (cluster=%s, local=%s)\n",
                    ifelse(is.null(cM),"NULL","present"),
                    ifelse(is.null(lM),"NULL","present")))
        return(invisible())
    }
    if(!all(dim(cM) == dim(lM))) {
        cat(sprintf("  DIM MISMATCH: cluster=%s  local=%s\n",
                    paste(dim(cM), collapse="x"), paste(dim(lM), collapse="x")))
        return(invisible())
    }
    d <- abs(cM - lM)
    cat(sprintf("  max abs diff overall: %.4g\n", max(d, na.rm = TRUE)))
    pc <- apply(d, 2, max, na.rm = TRUE)
    cat("  per-column max diff:\n"); print(round(pc, 4))
}
cmp("scoresFE")
cmp("scoresHyper")
cmp("scoresArea")
