# Shared validation utilities for the FE+nug+GH (and later BYM2) family of models.
#
# Architecture:
#   - buildFullMDMInputs(): build the canonical TMB input list once on the full DHS+MICS
#     data, using the existing makeInputsMDM machinery.
#   - subsetMDMInputs(fullInp, dhsKeep, micsKeep): subset all DHS/MICS-row pieces by
#     two logical vectors. The output is the same shape as makeInputsMDM's output, so it
#     can be passed directly as `inputsMDM=...` to fitFEM, fitFED, fitFEMD, fitFEMD_comm.
#   - predFEclusters(): produce posterior draws of held-out cluster prevalences for
#     FE+nug+GH model output.
#   - formatPredsForScoring(): wrap probDraws into the list shape that scoreValidationPreds
#     expects, save under savedOutput/validation/folds/preds<Tag>_fold<N>.RData.
#
# Fold assignment uses the existing stratKmics() / stratKdhs() functions in validation.R
# (sourced separately). DHS folds = 1..10, MICS folds = 11..20.

# ---- Stratified K-fold assignment (mirrors stratKdhs / stratKmics) -------------
# strat: vector of stratum labels (length n)
# urb:   logical vector of urbanicity (length n)
# Returns a length-n integer fold vector with values in 1..K.
assignKfoldStratified = function(strat, urb, K=10, seed=1234) {
    set.seed(seed)
    n = length(strat)
    foldVec = integer(n)
    inds    = 1:n
    uniqueStrata = sort(unique(strat))

    drawFold = function(thisStrat, thisUrb) {
        thisInds = inds[(strat == thisStrat) & (urb == thisUrb)]
        thisN = length(thisInds)
        if(thisN == 0) return(invisible(NULL))
        nBig = ((thisN - 1) %% K) + 1
        bigSize = ceiling(thisN / K)
        smallSize = bigSize - 1
        bigInds = sample(1:K, nBig)
        for(k in 1:K) {
            thisFoldSize = if(k %in% bigInds) bigSize else smallSize
            if(thisFoldSize == 0) next
            if(length(thisInds) == 1) {
                if(thisFoldSize == 1) tempInds = thisInds
                else stop("singleton stratum but thisFoldSize > 1")
            } else {
                tempInds = sample(thisInds, thisFoldSize)
            }
            thisInds = thisInds[-match(tempInds, thisInds)]
            foldVec[tempInds] <<- k
        }
    }

    for(s in uniqueStrata) {
        drawFold(s, TRUE)
        drawFold(s, FALSE)
    }
    foldVec
}

# ---- Build full MDM inputs once -------------------------------------------------
buildFullMDMInputs = function(KMICS=100, KDHSurb=16, KDHSrur=21,
                              JInnerUrban=4, JInnerRural=4, JOuterRural=1,
                              datDHS=NULL, datMICS=NULL,
                              admMICS=admFinal, adm2DHS=adm2Full,
                              adm2AsCovariate=FALSE,
                              K=10, seed=1234) {
    if(is.null(datDHS))  { load("savedOutput/global/ed.RData");     datDHS  = ed }
    if(is.null(datMICS)) { load("savedOutput/global/edMICS.RData"); datMICS = edMICS }
    if(!("Stratum" %in% names(datMICS))) {
        datMICS$Stratum = adm2ToStratumMICS(datMICS$subarea)
    }

    # Assign folds BEFORE any reordering, using the same stratification keys as
    # stratKdhs / stratKmics in code/validation.R. The fold column then travels with
    # each row through subsequent reorders, so we never have an alignment problem.
    datDHS$fold  = assignKfoldStratified(datDHS$admin1,    datDHS$urban,  K=K, seed=seed)
    datMICS$fold = assignKfoldStratified(datMICS$StratumID, datMICS$urban, K=K, seed=seed)

    datMICS = sortByCol(datMICS, "Stratum", admMICS$NAME_FINAL)

    # Match makeInputsMDM's "N -> ns/n, Z -> ys/y" column normalization.
    nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
    for(i in 1:nrow(nameTab)) {
        fromN = nameTab[i,1]; toN = nameTab[i,2]
        if(!(toN %in% names(datDHS)))  datDHS[[toN]]  = datDHS[[fromN]]
        if(!(toN %in% names(datMICS))) datMICS[[toN]] = datMICS[[fromN]]
    }

    inp = makeInputsMDM(datDHS, datMICS,
                        intPtsMICS=NULL, intPtsDHS=NULL,
                        KMICS=KMICS,
                        KDHSurb=KDHSurb, JInnerUrban=JInnerUrban,
                        KDHSrur=KDHSrur, JInnerRural=JInnerRural,
                        JOuterRural=JOuterRural,
                        admMICS=admMICS, adm2DHS=adm2DHS,
                        adm2AsCovariate=adm2AsCovariate)

    # Stash the cluster-row identifiers so subsetMDMInputs knows which DHS / MICS
    # cluster each row corresponds to. makeInputsMDM groups by urban/rural; we
    # need the urban-flagged DHS clusters and urban-flagged MICS clusters in their
    # internal ordering.
    dhsUrb  = datDHS$urban
    micsUrb = datMICS$urban

    # Indices of urb/rur clusters in the original ed and edMICS row ordering.
    inp$dhsUrbIdx  = which(dhsUrb)
    inp$dhsRurIdx  = which(!dhsUrb)
    inp$micsUrbIdx = which(micsUrb)
    inp$micsRurIdx = which(!micsUrb)

    inp$nDHSurb  = length(inp$dhsUrbIdx)
    inp$nDHSrur  = length(inp$dhsRurIdx)
    inp$nMICSurb = length(inp$micsUrbIdx)
    inp$nMICSrur = length(inp$micsRurIdx)

    inp$datDHS   = datDHS    # store originals for later out-of-sample lookup
    inp$datMICS  = datMICS
    inp$KMICS    = KMICS
    inp$KDHSurb  = KDHSurb
    inp$KDHSrur  = KDHSrur

    inp
}

# ---- Subset full MDM inputs by per-cluster keep-flags ---------------------------
# dhsKeep:  logical vector, length = nrow(datDHS),  TRUE for clusters to retain
# micsKeep: logical vector, length = nrow(datMICS), TRUE for clusters to retain
subsetMDMInputs = function(fullInp, dhsKeep, micsKeep) {
    out = fullInp
    nUobs = fullInp$nDHSurb;  nRobs = fullInp$nDHSrur
    nUmic = fullInp$nMICSurb; nRmic = fullInp$nMICSrur

    keepUDHS  = dhsKeep[fullInp$dhsUrbIdx]
    keepRDHS  = dhsKeep[fullInp$dhsRurIdx]
    keepUMICS = micsKeep[fullInp$micsUrbIdx]
    keepRMICS = micsKeep[fullInp$micsRurIdx]

    # Subset DHS pieces ----------------------------------------------------------
    out$ysUrbDHS  = fullInp$ysUrbDHS[keepUDHS]
    out$nsUrbDHS  = fullInp$nsUrbDHS[keepUDHS]
    out$ysRurDHS  = fullInp$ysRurDHS[keepRDHS]
    out$nsRurDHS  = fullInp$nsRurDHS[keepRDHS]
    out$areaidxlocUrbanDHS = fullInp$areaidxlocUrbanDHS[keepUDHS]
    out$areaidxlocRuralDHS = fullInp$areaidxlocRuralDHS[keepRDHS]

    # intPtsDHS: design matrices are stored as nObs * Kint rows, indexed (mod nObs) over k.
    intD = fullInp$intPtsDHS
    KU = nrow(intD$covsUrb) / nUobs
    KR = nrow(intD$covsRur) / nRobs

    # Row indices to keep for the [nObs * K] x p design matrix:
    rowsU = as.vector(outer(which(keepUDHS), (0:(KU-1))*nUobs, "+"))
    rowsR = as.vector(outer(which(keepRDHS), (0:(KR-1))*nRobs, "+"))

    out$intPtsDHS = intD
    out$intPtsDHS$covsUrb = intD$covsUrb[rowsU, , drop=FALSE]
    out$intPtsDHS$covsRur = intD$covsRur[rowsR, , drop=FALSE]
    out$intPtsDHS$wUrban  = intD$wUrban[keepUDHS, , drop=FALSE]
    out$intPtsDHS$wRural  = intD$wRural[keepRDHS, , drop=FALSE]

    # Subset MICS pieces ---------------------------------------------------------
    out$ysUrbMICS = fullInp$ysUrbMICS[keepUMICS]
    out$nsUrbMICS = fullInp$nsUrbMICS[keepUMICS]
    out$ysRurMICS = fullInp$ysRurMICS[keepRMICS]
    out$nsRurMICS = fullInp$nsRurMICS[keepRMICS]
    out$areaidxlocUrbanMICS = fullInp$areaidxlocUrbanMICS[keepUMICS]
    out$areaidxlocRuralMICS = fullInp$areaidxlocRuralMICS[keepRMICS]

    intM = fullInp$intPtsMICS
    KMu = nrow(intM$XUrb) / nUmic
    KMr = nrow(intM$XRur) / nRmic
    rowsMU = as.vector(outer(which(keepUMICS), (0:(KMu-1))*nUmic, "+"))
    rowsMR = as.vector(outer(which(keepRMICS), (0:(KMr-1))*nRmic, "+"))

    out$intPtsMICS = intM
    out$intPtsMICS$XUrb   = intM$XUrb[rowsMU, , drop=FALSE]
    out$intPtsMICS$XRur   = intM$XRur[rowsMR, , drop=FALSE]
    out$intPtsMICS$wUrban = intM$wUrban[keepUMICS, , drop=FALSE]
    out$intPtsMICS$wRural = intM$wRural[keepRMICS, , drop=FALSE]

    # Refresh the bookkeeping so a downstream call (e.g. another subset) is correct.
    out$dhsUrbIdx  = fullInp$dhsUrbIdx[keepUDHS]
    out$dhsRurIdx  = fullInp$dhsRurIdx[keepRDHS]
    out$micsUrbIdx = fullInp$micsUrbIdx[keepUMICS]
    out$micsRurIdx = fullInp$micsRurIdx[keepRMICS]
    out$nDHSurb    = sum(keepUDHS)
    out$nDHSrur    = sum(keepRDHS)
    out$nMICSurb   = sum(keepUMICS)
    out$nMICSrur   = sum(keepRMICS)

    out
}

# ---- FE+nug+GH prediction for held-out clusters ---------------------------------
# Inputs:
#   fitRes:        result list from fitFEM / fitFED / fitFEMD / fitFEMD_comm (TMB path)
#   outOfSample:   list with components covsUrb (matrix, nObs*K x p), covsRur, wUrban,
#                  wRural, ys/nsUrb, ys/nsRur — same structure as the corresponding
#                  pieces of buildFullMDMInputs / subsetMDMInputs output.
#   nObsUrb, nObsRur: number of held-out urban/rural clusters
#   nsim:          number of posterior draws
#   side:          "DHS" (default) or "MICS" — which side's alpha/beta to use for the
#                  commensurate fusion model (for the standard models there is only one).
predFEclusters = function(fitRes, outOfSample,
                          nObsUrb, nObsRur, nsim=10000,
                          side=c("DHS", "MICS"), areaidx=NULL) {
    side = match.arg(side)
    SD = fitRes$TMBsd
    if(!inherits(SD, "sdreport")) {
        stop("fitRes$TMBsd is not a TMB sdreport — cannot draw from posterior")
    }

    # Posterior parameter draws (params x nsim, rownames = parameter names),
    # reusing the reviewed draw machinery in inlaStyleDraws.R (sourced via setup.R):
    #   - jointPrecision present  -> posteriorDraws(useInla="auto"): joint Gaussian
    #     draws when the jointPrecision Cholesky succeeds, else automatic INLA-style
    #     fall-back when the joint Hessian is non-PD (e.g. the commensurate model
    #     near a sigma_comm boundary) -- same robustness as the BYM2/predGrid path.
    #   - no jointPrecision but the model HAS random effects (e.g. commensurate with
    #     marginalised FEs, when sdreport only succeeded at getJointPrecision=FALSE)
    #     -> INLA-style draws (cov.fixed alone can't supply the random alpha/beta).
    #   - no jointPrecision and no random effects (pure FE) -> cov.fixed Gaussian
    #     draws (the historical FE path).
    if(!is.null(SD$jointPrecision)) {
        parDraws = posteriorDraws(fitRes, NDRAWS=nsim, useInla="auto")
    } else if(length(SD$par.random) > 0) {
        parDraws = inlaStyleDraws(fitRes, NDRAWS=nsim)
    } else {
        parDraws = .gaussianPosteriorDraws(fitRes, NDRAWS=nsim, requireJoint=FALSE)
    }
    nm = rownames(parDraws)

    # Pick the right rows for alpha and beta by name. The commensurate model has
    # both "alpha"/"alpha_M" and "beta"/"beta_M"; standard models have just "alpha"/"beta".
    if(side == "DHS" || !("alpha_M" %in% nm)) {
        alphaName = "alpha";   betaName = "beta"
    } else {
        alphaName = "alpha_M"; betaName = "beta_M"
    }
    alphaDraws    = parDraws[which(nm == alphaName), , drop=TRUE]
    betaRows      = which(nm == betaName)
    betaDraws     = parDraws[betaRows, , drop=FALSE]            # p x nsim
    logTauEpsRow  = which(nm == "log_tauEps")
    logTauEpsDraw = parDraws[logTauEpsRow, , drop=TRUE]
    sigmaEpsDraw  = exp(-0.5 * logTauEpsDraw)

    # BYM2 spatial field, if this is a BYM2 fit (w_bym2Free in the draws) and the
    # held-out clusters' area indices were supplied. The GH templates apply the
    # field PER CLUSTER at its nominal area (base_eta = alpha + Xbeta +
    # w_bym2[areaidx(i)]; constant over a cluster's jitter points), so we add a
    # per-cluster field term alongside alpha. Field is reconstructed exactly as
    # predGrid does: w_bym2 = c(w_bym2Free, -sum(w_bym2Free)) over the nAreas areas.
    # areaidx is 0-based (C++ index), so +1 for R. FE fits (no w_bym2Free) skip this.
    fieldFor = function(idx0) NULL
    if("w_bym2Free" %in% nm && !is.null(areaidx)) {
        wFree = parDraws[which(nm == "w_bym2Free"), , drop=FALSE]   # (nAreas-1) x nsim
        bym2Field = rbind(wFree, -colSums(wFree))                   # nAreas x nsim
        fieldFor = function(idx0) {
            if(length(idx0) == 0) return(matrix(0, nrow=0, ncol=nsim))
            bym2Field[idx0 + 1L, , drop=FALSE]                      # nObs x nsim (per cluster)
        }
    }

    XUrb = as.matrix(outOfSample$covsUrb)   # nObsUrb*KU x p
    XRur = as.matrix(outOfSample$covsRur)   # nObsRur*KR x p
    wUrb = as.matrix(outOfSample$wUrban)    # nObsUrb x KU
    wRur = as.matrix(outOfSample$wRural)    # nObsRur x KR
    KU = ncol(wUrb); KR = ncol(wRur)

    drawClusterPreds = function(X, w, nObs, K, fieldMat) {
        if(nObs == 0) return(matrix(0, nrow=0, ncol=nsim))
        # X is nObs*K x p with X-row index = nObs*(k-1) + i  (1-indexed k, i)
        # X %*% beta yields a length nObs*K vector per draw.
        XBeta = X %*% betaDraws                              # (nObs*K) x nsim
        XBetaArr = array(as.matrix(XBeta), dim=c(nObs, K, nsim))

        eps = matrix(rnorm(nObs * nsim), nrow=nObs, ncol=nsim) *
              matrix(rep(sigmaEpsDraw, each=nObs), nrow=nObs, ncol=nsim)
        alphaMat = matrix(alphaDraws, nrow=nObs, ncol=nsim, byrow=TRUE)
        if(!is.null(fieldMat) && nrow(fieldMat) == nObs) alphaMat = alphaMat + fieldMat  # + per-cluster BYM2 field

        out = matrix(0, nrow=nObs, ncol=nsim)
        for(k in 1:K) {
            etaK = alphaMat + XBetaArr[, k, ] + eps          # nObs x nsim
            pK   = 1 / (1 + exp(-etaK))
            out  = out + w[, k] * pK
        }
        out
    }

    idxUrb = if(is.null(areaidx)) integer(0) else areaidx$urb
    idxRur = if(is.null(areaidx)) integer(0) else areaidx$rur
    probDrawsUrb = drawClusterPreds(XUrb, wUrb, nObsUrb, KU, fieldFor(idxUrb))
    probDrawsRur = drawClusterPreds(XRur, wRur, nObsRur, KR, fieldFor(idxRur))

    list(probDrawsUrb=probDrawsUrb, probDrawsRur=probDrawsRur,
         predsUrb=if(nObsUrb > 0) rowMeans(probDrawsUrb) else numeric(0),
         predsRur=if(nObsRur > 0) rowMeans(probDrawsRur) else numeric(0))
}

# ---- Wrap a prediction object into the structure scoreValidationPreds expects ----
# Saves to savedOutput/validation/folds/preds<tag>_fold<fold>.RData
formatAndSavePreds = function(predObj, ysUrb, nsUrb, ysRur, nsRur,
                              fold, tag, totalTime,
                              outDir="savedOutput/validation/folds") {
    if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    preds = list(
        probDrawsUrb = predObj$probDrawsUrb,
        probDrawsRur = predObj$probDrawsRur,
        predsUrb     = predObj$predsUrb,
        predsRur     = predObj$predsRur,
        quantsUrb    = if(length(predObj$predsUrb) > 0)
                          t(apply(predObj$probDrawsUrb, 1, quantile,
                                  probs=c(0.025, 0.1, 0.5, 0.9, 0.975))) else NULL,
        quantsRur    = if(length(predObj$predsRur) > 0)
                          t(apply(predObj$probDrawsRur, 1, quantile,
                                  probs=c(0.025, 0.1, 0.5, 0.9, 0.975))) else NULL,
        yUrb         = ysUrb,
        nUrb         = nsUrb,
        yRur         = ysRur,
        nRur         = nsRur
    )
    fname = file.path(outDir, paste0("preds", tag, "_fold", fold, ".RData"))
    save(preds, totalTime, file=fname)
    invisible(fname)
}

# ---- Score one tag's preds file using the existing getScores() ---------------
# Reads savedOutput/validation/folds/preds<tag>_fold<N>.RData (saved by
# formatAndSavePreds) and computes urban/rural/combined scores via getScores from
# code/scores.R. Saves scores<tag>_fold<N>.RData with scoresUrb / scoresRur / scores.
scoreTagFold = function(tag, fold,
                        outDir="savedOutput/validation/folds",
                        useSampleWt=FALSE) {
    predFile  = file.path(outDir, paste0("preds",  tag, "_fold", fold, ".RData"))
    scoreFile = file.path(outDir, paste0("scores", tag, "_fold", fold, ".RData"))
    if(!file.exists(predFile)) return(invisible(NULL))

    e = new.env(); load(predFile, envir=e)
    preds      = e$preds
    totalTime  = e$totalTime

    yUrb = preds$yUrb; nUrb = preds$nUrb
    yRur = preds$yRur; nRur = preds$nRur
    probDrawsUrb = preds$probDrawsUrb
    probDrawsRur = preds$probDrawsRur

    # If the held-out side had no clusters of a given urbanicity, getScores would
    # error on length-0 input — handle each piece guardedly.
    safeScore = function(y, n, draws) {
        if(length(y) == 0) return(NULL)
        getScores(truth=y/n, estMat=draws, weights=n,
                  significance=c(.5, .8, .9, .95), doFuzzyReject=TRUE,
                  getAverage=TRUE, na.rm=TRUE, ns=n)
    }

    scoresUrb = safeScore(yUrb, nUrb, probDrawsUrb)
    scoresRur = safeScore(yRur, nRur, probDrawsRur)

    yAll = c(yUrb, yRur); nAll = c(nUrb, nRur)
    drawsAll = rbind(probDrawsUrb, probDrawsRur)
    scores = safeScore(yAll, nAll, drawsAll)

    if(!is.null(scoresUrb)) scoresUrb = cbind(scoresUrb, Time=totalTime/60)
    if(!is.null(scoresRur)) scoresRur = cbind(scoresRur, Time=totalTime/60)
    if(!is.null(scores))    scores    = cbind(scores,    Time=totalTime/60)

    save(scoresUrb, scoresRur, scores, preds, file=scoreFile)
    invisible(scoreFile)
}

# Score every per-fold preds file produced by runValidationFE.R for a given model.
# model must be one of "M_D", "M_M", "M_DM", "M_DM_comm".
scoreAllFoldsFE = function(model, outDir="savedOutput/validation/folds") {
    for(j in 1:20) {
        side = if(j <= 10) "DHS" else "MICS"
        # Only valid (model, side) combos exist; others were skipped at fit time.
        valid = switch(model,
                       "M_D"       = side == "DHS",
                       "M_M"       = side == "MICS",
                       "M_DM"      = TRUE,
                       "M_DM_comm" = TRUE,
                       FALSE)
        if(!valid) next
        tag = paste0(model, ifelse(side == "DHS", "_dhsFold", "_micsFold"),
                     ifelse(side == "DHS", j, j - 10))
        scoreTagFold(tag, fold=j)
    }
}

# ---- Decode (model, fold) from a single SLURM array index ----------------------
# Returns a list with: model (character), fold (1..20), side ("DHS" or "MICS"),
# isValid (whether this combination is meaningful for the given model).
# Models for FE family:  i = 1..4 = M_D, M_M, M_DM, M_DM_comm
decodeFEjob = function(index, maxJ=20) {
    jobInds = getJobIndices(index, maxJ=maxJ, rev=TRUE)
    i = jobInds[1]; j = jobInds[2]
    model = c("M_D", "M_M", "M_DM", "M_DM_comm")[i]
    side  = if(j <= 10) "DHS" else "MICS"
    foldDHS  = if(side == "DHS")  j         else NA_integer_
    foldMICS = if(side == "MICS") j - 10    else NA_integer_

    isValid = switch(model,
                     "M_D"        = side == "DHS",
                     "M_M"        = side == "MICS",
                     "M_DM"       = TRUE,
                     "M_DM_comm"  = TRUE,
                     FALSE)

    list(model=model, fold=j, side=side,
         foldDHS=foldDHS, foldMICS=foldMICS, isValid=isValid,
         tag=paste0(model, ifelse(side == "DHS", "_dhsFold", "_micsFold"),
                    ifelse(side == "DHS", foldDHS, foldMICS)))
}
