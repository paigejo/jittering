# Identify which BYM2 sim (1-10) differs between old/new and characterize.
oldE <- new.env(); load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData.bak", envir = oldE)
newE <- new.env(); load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData",     envir = newE)
for(i in 1:10) {
    od <- oldE$surveysDHS[[i]];  nd <- newE$surveysDHS[[i]]
    om <- oldE$surveysMICS[[i]]; nm <- newE$surveysMICS[[i]]
    dDhs  <- !isTRUE(all.equal(od$Z, nd$Z))
    dMics <- !isTRUE(all.equal(om$Z, nm$Z))
    coordDhsDiff  <- !isTRUE(all.equal(od$east, nd$east))
    coordMicsDiff <- !isTRUE(all.equal(om$east, nm$east))
    if(dDhs || dMics) {
        cat(sprintf("sim %d:  DHS Z differ=%s  DHS coords differ=%s  MICS Z differ=%s  MICS coords differ=%s\n",
                    i, dDhs, coordDhsDiff, dMics, coordMicsDiff))
        cat(sprintf("   DHS: oldZ=%d newZ=%d  oldN=%d newN=%d\n",
                    sum(od$Z), sum(nd$Z), sum(od$N), sum(nd$N)))
        cat(sprintf("   MICS: oldZ=%d newZ=%d  oldN=%d newN=%d\n",
                    sum(om$Z), sum(nm$Z), sum(om$N), sum(nm$N)))
    }
}
