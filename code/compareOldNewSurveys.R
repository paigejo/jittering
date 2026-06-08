# Per-survey-type comparison: does BYM2 have identical DHS+MICS to May-21,
# and SPDE have changed DHS+MICS? Or is there an asymmetry between DHS and
# MICS within either generative model?

for(tag in c("BYM2","SPDE")) {
    cat(sprintf("\n========== %s ==========\n", tag))
    oldE <- new.env()
    load(sprintf("savedOutput/simStudy1/simPopsSurveys_%s.RData.bak", tag), envir = oldE)
    newE <- new.env()
    load(sprintf("savedOutput/simStudy1/simPopsSurveys_%s.RData",     tag), envir = newE)
    for(svName in c("surveysDHS","surveysMICS")) {
        cat(sprintf("\n  %s:\n", svName))
        nMatch    <- 0
        coordsAll <- 0
        zAll      <- 0
        for(i in 1:10) {
            o <- get(svName, envir = oldE)[[i]]
            n <- get(svName, envir = newE)[[i]]
            coordEq <- isTRUE(all.equal(o$east,  n$east )) &&
                       isTRUE(all.equal(o$north, n$north)) &&
                       isTRUE(all.equal(o$urban, n$urban))
            zEq     <- isTRUE(all.equal(o$Z, n$Z))
            nEq     <- isTRUE(all.equal(o$N, n$N))
            if(coordEq) coordsAll <- coordsAll + 1
            if(zEq)     zAll      <- zAll + 1
            if(coordEq && zEq && nEq) nMatch <- nMatch + 1
        }
        cat(sprintf("    sims 1-10: coords match in %d/10  Z match in %d/10  all match in %d/10\n",
                    coordsAll, zAll, nMatch))
    }
}
