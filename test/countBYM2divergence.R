# How many of the 100 BYM2 sims actually match May-21? Once RNG drifts, all
# subsequent sims diverge — so we expect a clear cutover.
oldE <- new.env(); load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData.bak", envir = oldE)
newE <- new.env(); load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData",     envir = newE)
match <- logical(100)
for(i in 1:100) {
    od <- oldE$surveysDHS[[i]]; nd <- newE$surveysDHS[[i]]
    om <- oldE$surveysMICS[[i]]; nm <- newE$surveysMICS[[i]]
    match[i] <- isTRUE(all.equal(od$Z, nd$Z)) &&
                isTRUE(all.equal(od$east, nd$east)) &&
                isTRUE(all.equal(om$Z, nm$Z)) &&
                isTRUE(all.equal(om$east, nm$east))
}
cat(sprintf("BYM2 sims matching May-21: %d/100\n", sum(match)))
cat(sprintf("First non-matching: %d\n", which(!match)[1]))
cat(sprintf("Matching idxs: %s\n", paste(which(match), collapse=",")))
