# Cross-machine reproducibility test for Md (BYM2 generative), sim 17.
# - Cluster file is at savedOutput/simStudy1/scores_full/BYM2/scores_Md_sim17.RData
# - We re-fit locally with the same code and compare scoresFE / scoresHyper
#
# If sims are byte-identical between machines AND the fit is reproducible,
# scoresFE and scoresHyper should agree to within ~sqrt(2/NDRAWS) on
# Monte-Carlo-affected columns (Bias, sd, IntervalScore, etc.) and EXACTLY
# on deterministic columns (est, truth, ...) provided the fit converges to
# the same MLE.

source("code/setup.R")
options(warn = 1)

cat(sprintf("\n=== %s | xmach test: Md sim 17 (BYM2 generative) ===\n",
            format(Sys.time())))

CLUSTER_FILE <- "savedOutput/simStudy1/scores_full/BYM2/scores_Md_sim17.RData"
LOCAL_OUT    <- "savedOutput/simStudy1/scores_xmach/scores_Md_sim17.RData"

# Sanity check that the cluster file is there.
if(!file.exists(CLUSTER_FILE)) stop("Missing cluster file: ", CLUSTER_FILE)

# Load sim 17 surveys from the BYM2 simPops file. Note: this is the LOCAL
# version. If the user synced their NEW (post-sum-to-zero patch) sims to the
# cluster, both machines used the same surveys for sim 17 -> apples-to-apples.
simEnv <- new.env()
simulateSurveys("bym2", nsim = 17, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

areaPops <- simEnv$areaPops_smoothRisk           # default truthQuantity = smoothRisk
truths   <- modelTruths("bym2")

dir.create(dirname(LOCAL_OUT), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\n=== %s | refitting Md sim 17 locally ===\n", format(Sys.time())))
.scoreOneFit(simIdx  = 17,
             modName = "Md",
             outFile = LOCAL_OUT,
             surveysDHS  = simEnv$surveysDHS,
             surveysMICS = simEnv$surveysMICS,
             intPtsMICS  = micsEnv$intPtsMICS,
             model       = "bym2",
             truths      = truths,
             areaPops    = areaPops,
             KMICS = 100, KDHSu = 16, KDHSr = 21,
             Qgh = 10, NDRAWS = 1000,
             COVS = c("urban","access","elev","distRiversLakes","normPop"))

cat(sprintf("\n=== %s | comparison ===\n", format(Sys.time())))

cE <- new.env(); load(CLUSTER_FILE, envir = cE)
lE <- new.env(); load(LOCAL_OUT,    envir = lE)

cmp <- function(field) {
    cat(sprintf("\n----- %s -----\n", field))
    cM <- if(exists(field, envir = cE)) as.matrix(get(field, envir = cE)) else NULL
    lM <- if(exists(field, envir = lE)) as.matrix(get(field, envir = lE)) else NULL
    if(is.null(cM) || is.null(lM)) {
        cat(sprintf("  cluster=%s  local=%s\n",
                    ifelse(is.null(cM),"NULL","present"),
                    ifelse(is.null(lM),"NULL","present")))
        return(invisible())
    }
    if(!all(dim(cM) == dim(lM))) {
        cat(sprintf("  DIM MISMATCH: cluster %s vs local %s\n",
                    paste(dim(cM), collapse="x"), paste(dim(lM), collapse="x")))
        return(invisible())
    }
    cat(sprintf("  dim: %s\n", paste(dim(cM), collapse=" x ")))
    cat("  cluster:\n"); print(round(cM, 4))
    cat("  local:\n");   print(round(lM, 4))
    d <- abs(cM - lM)
    cat(sprintf("  max abs diff overall: %.4g\n", max(d, na.rm = TRUE)))
    cat("  per-column max abs diff:\n"); print(round(apply(d, 2, max, na.rm = TRUE), 4))
}

cmp("scoresFE")
cmp("scoresHyper")
cmp("scoresArea")

# Also show fit/score wallclock for reference
cat("\n----- wallclock -----\n")
cat(sprintf("  cluster fit=%.2f score=%.2f total=%.2f min\n",
            cE$fitTime, cE$scoreTime, cE$fitTime + cE$scoreTime))
cat(sprintf("  local   fit=%.2f score=%.2f total=%.2f min\n",
            lE$fitTime, lE$scoreTime, lE$fitTime + lE$scoreTime))

cat(sprintf("\n=== %s | DONE ===\n", format(Sys.time())))
