# Did the SPDE surveysDHS east/north/urban change after the doSmoothRisk patch?
# If not, the cached intPtsDHS_SPDE_*.RData files are still valid.

oldEnv <- new.env()
load("savedOutput/simStudy1/simPopsSurveys_SPDE.RData.bak", envir = oldEnv)
newEnv <- new.env()
load("savedOutput/simStudy1/simPopsSurveys_SPDE.RData",     envir = newEnv)

for(i in 1:3) {
    o <- oldEnv$surveysDHS[[i]]; n <- newEnv$surveysDHS[[i]]
    coordSame <- isTRUE(all.equal(o$east,  n$east )) &&
                 isTRUE(all.equal(o$north, n$north)) &&
                 isTRUE(all.equal(o$urban, n$urban))
    zSame     <- isTRUE(all.equal(o$Z, n$Z))
    nSame     <- isTRUE(all.equal(o$N, n$N))
    cat(sprintf("sim %d: nrow %d→%d  coords=%s  Z=%s  N=%s\n",
                i, nrow(o), nrow(n),
                ifelse(coordSame, "SAME", "DIFFER"),
                ifelse(zSame,     "SAME", "DIFFER"),
                ifelse(nSame,     "SAME", "DIFFER")))
}
