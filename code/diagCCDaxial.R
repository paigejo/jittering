# Does warm-start actually drive the inner mode at an EXTREME axial CCD point,
# and does my reuse warm-start write take effect? At a theta with logit_phi
# pushed toward the boundary (near-singular BYM2 inner problem), compute the
# inner mode three ways and compare:
#   (R)  rebuild, innerWarm = centre mode      (.makeFixedHyperObj path)
#   (Wc) reuse,   innerStart = centre mode      (my .evalCCDreuse write)
#   (Wz) reuse,   innerStart = zeros            (deliberately different init)
# If R == Wc  != Wz  -> warm-start drives the mode AND my write works; the
#   walk difference is then because the reuse WALK isn't really centre-warming.
# If Wc == Wz         -> my write is ignored (wrong TMB field) OR mode unique.
# If R != Wc          -> rebuild vs reuse differ even with identical intended
#   init -> deeper structural difference.

source("code/setup.R")
options(warn = 1)

SIMIDX <- 17
simEnv <- new.env(); simulateSurveys("bym2", nsim = SIMIDX, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env(); load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
ip <- buildInputsForSim(SIMIDX, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2", KMICS = 100, KDHSu = 16, KDHSr = 21)
res <- .fitOne("M_D_BYM2", ip, KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
               COVS = c("urban","access","elev","distRiversLakes","normPop"))
SD  <- res$TMBsd
obj <- res$TMBobj
dataList <- obj$env$data; paramsTemplate <- obj$env$parameters; DLL <- obj$env$DLL
hyperNames <- c("log_tau","logit_phi","log_tauEps")

# centre mode (inner) from a fresh mapped obj at the hyper modes
thetaC <- SD$par.fixed[hyperNames]
objC <- .makeFixedHyperObj(dataList, paramsTemplate, thetaC, DLL, innerWarm = NULL)
invisible(objC$fn(numeric())); centreInner <- objC$env$last.par.best
cat(sprintf("centre inner mode captured (n=%d)\n", length(centreInner)))

# EXTREME axial theta: push logit_phi well toward +boundary (phi -> ~0.95+)
thetaX <- thetaC; thetaX["logit_phi"] <- thetaC["logit_phi"] + 3.0
cat("extreme theta:\n"); print(thetaX)

# (R) rebuild with innerWarm = centre
objR <- .makeFixedHyperObj(dataList, paramsTemplate, thetaX, DLL,
                           innerWarm = list(
                             alpha      = as.numeric(centreInner[names(centreInner)=="alpha"]),
                             beta       = as.numeric(centreInner[names(centreInner)=="beta"]),
                             w_bym2Free = as.numeric(centreInner[names(centreInner)=="w_bym2Free"]),
                             u_bym2Free = as.numeric(centreInner[names(centreInner)=="u_bym2Free"])))
invisible(objR$fn(numeric())); muR <- objR$env$last.par.best

# reuse obj (shared), warm-start via last.par write
objW <- MakeADFun(data = dataList, parameters = paramsTemplate,
                  random = c("alpha","beta","w_bym2Free","u_bym2Free"),
                  DLL = DLL, silent = TRUE)
randIdx <- objW$env$random
setInner <- function(o, vals) {
    lp <- o$env$last.par; lp[o$env$random] <- vals; o$env$last.par <- lp
    if(!is.null(o$env$last.par.best)) { lpb <- o$env$last.par.best; lpb[o$env$random] <- vals; o$env$last.par.best <- lpb }
}
thX <- as.numeric(thetaX[names(objW$par)])

setInner(objW, centreInner); invisible(objW$fn(thX)); muWc <- objW$env$last.par[randIdx]
setInner(objW, rep(0, length(randIdx))); invisible(objW$fn(thX)); muWz <- objW$env$last.par[randIdx]

cat(sprintf("\nmax|muR - muWc| (rebuild-centre vs reuse-centre): %.3e\n", max(abs(muR - muWc))))
cat(sprintf("max|muWc - muWz| (reuse-centre vs reuse-zeros):    %.3e\n", max(abs(muWc - muWz))))
cat(sprintf("max|muR  - muWz| (rebuild-centre vs reuse-zeros):  %.3e\n", max(abs(muR  - muWz))))
cat("\nper-group max|muR - muWc|:\n")
for(nm in unique(names(muR)))
    cat(sprintf("  %-12s %.3e\n", nm, max(abs(muR[names(muR)==nm] - muWc[names(muWc)==nm]))))
cat("\nper-group max|muWc - muWz| (does init change the mode here?):\n")
for(nm in unique(names(muWc)))
    cat(sprintf("  %-12s %.3e\n", nm, max(abs(muWc[names(muWc)==nm] - muWz[names(muWz)==nm]))))
cat("\n=== done ===\n")
