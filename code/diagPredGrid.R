# predGrid decomposition on a real Md fit (BYM2 generative, sim 17).
#
# Goal: localize the +0.11 areal bias. With TRUE parameters the nugget-
# integrated, target-pop-weighted prediction is ~0.30 national while the
# stored pSmoothRisk truth is ~0.34 — yet scored predictions sit ~0.45
# (prevalence level). Suspects:
#   S1 aggregation weights (predArea: total pop vs truth: target pop)
#   S2 sigma^2-vs-sigma slip in predGrid's nugget integration
#   S3 fitted latent level (intercept + BYM2) sitting too high
#
# Outputs:
#   - fitted params vs truth (alpha, sigmaBYM2, phi, sigmaEps)
#   - pixel-level prediction mean under total-pop and target-pop weighting
#   - per-area bias vs pSmoothRisk truth and vs prevalence truth, using the
#     exact .scoreArea aggregation path (predArea) AND a manual target-pop
#     aggregation
#   - mean bias reproduced (~+0.11?) and which change removes it

source("code/setup.R")
options(warn = 1)

SIMIDX <- 17
cat(sprintf("\n=== %s | predGrid decomposition: BYM2 sim %d ===\n",
            format(Sys.time()), SIMIDX))

simEnv <- new.env()
simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)

truthSR <- simEnv$areaPops_smoothRisk[, SIMIDX]
truthPv <- simEnv$areaPops[, SIMIDX]
truths  <- modelTruths("bym2")
wTot    <- popMatNGAThresh$pop
wTar    <- popMatNGAedThresh$pop
areaFac <- popMatNGAThresh$area

decompose <- function(modName) {
    cat(sprintf("\n############ %s ############\n", modName))
    t0  <- proc.time()[3]
    res <- .fitOne(modName, ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
                   COVS = c("urban","access","elev","distRiversLakes","normPop"))
    cat(sprintf("fit time: %.1f min\n", (proc.time()[3] - t0)/60))

    pf <- res$TMBsd$par.fixed
    cat("--- fitted par.fixed (note position 1!) ---\n"); print(round(pf, 4))
    aHat <- if("alpha" %in% names(pf)) pf["alpha"]
            else if("alpha" %in% names(res$TMBsd$par.random))
                res$TMBsd$par.random[names(res$TMBsd$par.random) == "alpha"]
            else NA
    cat(sprintf("alpha-hat (by name) = %.4f (truth %.4f);  par.fixed[1] = %.4f <- what old predGrid used\n",
                aHat, -1.25, pf[1]))
    if("log_tauEps" %in% names(pf))
        cat(sprintf("sigmaEps-hat = %.4f (truth %.4f)\n",
                    exp(-0.5*pf["log_tauEps"]), truths$sigmaEps))

    grid <- tryCatch(
        predGrid(res$TMBsd, popMat = popMatNGAThresh, nsim = 1000,
                 obj = res$TMBobj, admLevel = "stratMICS", res = res),
        error = function(e) { cat("predGrid ERROR: ", conditionMessage(e), "\n"); NULL })
    if(is.null(grid)) return(invisible(NULL))

    predsPix <- grid$preds
    cat(sprintf("\npixel preds: mean=%.4f  TOTAL-pop wtd=%.4f  TARGET-pop wtd=%.4f\n",
                mean(predsPix, na.rm = TRUE),
                sum(wTot*predsPix, na.rm = TRUE)/sum(wTot[!is.na(predsPix)]),
                sum(wTar*predsPix, na.rm = TRUE)/sum(wTar[!is.na(predsPix)])))

    agg    <- tryCatch(predArea(grid, areaVarName = "area",
                                orderedAreas = adm1@data$NAME_1),
                       error = function(e) { cat("predArea ERROR: ",
                                                 conditionMessage(e), "\n"); NULL })
    estVec <- if(!is.null(agg)) {
        v <- agg$aggregationResults$p
        names(v) <- agg$aggregationResults$region
        v[match(names(truthSR), names(v))]
    } else NULL

    manTar <- tapply(seq_along(predsPix), areaFac, function(ii)
        sum(wTar[ii]*predsPix[ii], na.rm = TRUE)/sum(wTar[ii][!is.na(predsPix[ii])]))
    manTar <- manTar[match(names(truthSR), names(manTar))]
    manTot <- tapply(seq_along(predsPix), areaFac, function(ii)
        sum(wTot[ii]*predsPix[ii], na.rm = TRUE)/sum(wTot[ii][!is.na(predsPix[ii])]))
    manTot <- manTot[match(names(truthSR), names(manTot))]

    cat("\n---------------- AREA-LEVEL BIAS SUMMARY ----------------\n")
    cat(sprintf("mean truth pSmoothRisk        : %.4f\n", mean(truthSR, na.rm = TRUE)))
    cat(sprintf("mean truth prevalence         : %.4f\n", mean(truthPv, na.rm = TRUE)))
    if(!is.null(estVec))
        cat(sprintf("mean pred (predArea path)     : %.4f  bias vs SR: %+.4f  vs Prev: %+.4f\n",
                    mean(estVec, na.rm = TRUE),
                    mean(estVec - truthSR, na.rm = TRUE),
                    mean(estVec - truthPv, na.rm = TRUE)))
    cat(sprintf("mean pred (manual TOTAL-pop)  : %.4f  bias vs SR: %+.4f\n",
                mean(manTot, na.rm = TRUE), mean(manTot - truthSR, na.rm = TRUE)))
    cat(sprintf("mean pred (manual TARGET-pop) : %.4f  bias vs SR: %+.4f\n",
                mean(manTar, na.rm = TRUE), mean(manTar - truthSR, na.rm = TRUE)))
    invisible(NULL)
}

# Md_FE is where the +0.11 was observed; Md exercises the (previously
# crashing) Star path; M_D_BYM2 exercises the w_bym2Free GH path shared
# with M_M_BYM2 / M_DM_BYM2.
decompose("Md_FE")
decompose("Md")
decompose("M_D_BYM2")

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
