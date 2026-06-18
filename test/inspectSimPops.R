# Inventory a simPopsSurveys_<MODEL>.RData file: list each top-level object,
# its class/length, and (for surveysDHS/surveysMICS) the first sim's record
# counts. Run on BOTH machines and compare to localise the simPops divergence.

for(modelTag in c("BYM2","SPDE")) {
    p <- sprintf("savedOutput/simStudy1/simPopsSurveys_%s.RData", modelTag)
    cat(sprintf("\n##### %s -- %d bytes #####\n", p, file.info(p)$size))
    e <- new.env(); load(p, envir = e)
    for(o in ls(e)) {
        x <- get(o, envir = e)
        cls <- paste(class(x), collapse = "/")
        dims <- if(is.matrix(x) || is.data.frame(x)) paste(dim(x), collapse = "x")
                else if(is.list(x))                  paste0("len=", length(x))
                else                                  paste0("len=", length(x))
        cat(sprintf("  %-30s  %-20s  %s\n", o, cls, dims))
    }
    # Sim 1 of surveysDHS/surveysMICS — record-counts + Z totals
    for(svName in c("surveysDHS","surveysMICS")) {
        if(!exists(svName, envir = e)) next
        sv <- get(svName, envir = e)
        if(length(sv) < 1) next
        s1 <- sv[[1]]
        zSum  <- if("Z"  %in% names(s1)) sum(s1$Z)  else NA
        nSum  <- if("N"  %in% names(s1)) sum(s1$N)  else NA
        cat(sprintf("  %s[[1]]: nrow=%d  sum(Z)=%s  sum(N)=%s\n",
                    svName, nrow(s1), zSum, nSum))
    }
    # areaPops/subareaPops/stratumPops — sims 1 column-sums for fingerprinting
    for(popName in c("subareaPops","areaPops","stratumPops",
                     "subareaPops_smoothRisk","areaPops_smoothRisk","stratumPops_smoothRisk")) {
        if(!exists(popName, envir = e)) next
        m <- get(popName, envir = e)
        if(is.null(m)) next
        if(ncol(m) < 1) next
        cat(sprintf("  %-30s  dim=%s  sum(col1)=%.6f\n",
                    popName, paste(dim(m), collapse="x"), sum(m[, 1])))
    }
}
