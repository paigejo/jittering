# Validation re-score: M_DM_BYM2 sim 1 (BYM2 generative) using the CURRENT
# fixed code. Writes to a sandbox dir, then prints side-by-side comparison
# against the cluster's June-07 file.
#
# Hypothesis: the .postDraws Pt-step bug (May-21 code, fixed since)
# fully explains the local-vs-cluster discrepancy. If true, this local
# re-score should agree with the cluster file to within RNG noise on the
# 1000 posterior draws.

source("code/setup.R")
options(error = recover)

cat(sprintf("\n=== %s | START re-score M_DM_BYM2 sim1 ===\n", format(Sys.time())))

simEnv <- new.env()
simulateSurveys("bym2", nsim = 1, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

# Match cluster: fine-scale prevalence (smoothRisk field doesn't exist in the
# cluster's simPops file because that file pre-dates the simulator patch).
areaPops <- simEnv$areaPops
truths   <- modelTruths("bym2")

outDir  <- "savedOutput/simStudy1/scores_validate"
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
outFile <- file.path(outDir, "scores_M_DM_BYM2_sim1.RData")

cat(sprintf("\n=== %s | fitting + scoring (cluster wallclock ~36 min) ===\n",
            format(Sys.time())))

.scoreOneFit(simIdx  = 1,
             modName = "M_DM_BYM2",
             outFile = outFile,
             surveysDHS  = simEnv$surveysDHS,
             surveysMICS = simEnv$surveysMICS,
             intPtsMICS  = micsEnv$intPtsMICS,
             model       = "bym2",
             truths      = truths,
             areaPops    = areaPops,
             KMICS = 100, KDHSu = 16, KDHSr = 21,
             Qgh = 10, NDRAWS = 1000,
             COVS = c("urban","access","elev","distRiversLakes","normPop"))

cat(sprintf("\n=== %s | scoring done; running comparison ===\n", format(Sys.time())))

clustFile <- "c:/Users/jpaige/OneDrive - Norsk Regnesentral/Projects/jittering/jitterScores/scores_full_cluster/BYM2/scores_M_DM_BYM2_sim1.RData"
lE <- new.env(); load(outFile,   envir = lE)
cE <- new.env(); load(clustFile, envir = cE)

show <- function(field) {
    cat(sprintf("\n--- %s : cluster ---\n", field))
    if(exists(field, envir = cE)) print(round(get(field, envir = cE), 4)) else cat("(NULL)\n")
    cat(sprintf("\n--- %s : local re-score ---\n", field))
    if(exists(field, envir = lE)) print(round(get(field, envir = lE), 4)) else cat("(NULL)\n")
    if(exists(field, envir = cE) && exists(field, envir = lE)) {
        cM <- as.matrix(get(field, envir = cE))
        lM <- as.matrix(get(field, envir = lE))
        if(all(dim(cM) == dim(lM))) {
            cat(sprintf("\n--- %s : max|cluster - local| per col ---\n", field))
            print(round(apply(abs(cM - lM), 2, max, na.rm = TRUE), 4))
        } else {
            cat(sprintf("\n--- %s : DIM MISMATCH cluster=%s local=%s ---\n",
                        field,
                        paste(dim(cM), collapse="x"), paste(dim(lM), collapse="x")))
        }
    }
}

show("scoresFE")
show("scoresHyper")
show("scoresArea")

cat(sprintf("\n=== %s | DONE ===\n", format(Sys.time())))
