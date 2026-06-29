# Validation model set, FE fit dispatch, and per-fold cluster runner.
# Functions (not scripts) so both the local FE test and the cluster array reuse
# the SAME code. Reuses: subsetMDMInputs, makeMdInputs, makeInputsMdm,
# fitFED/fitFEM/fitFEMD, predFEclusters, formatAndSavePreds (all in
# validationCommon.R / model files), so nothing here is reimplemented.
#
# The 5 approaches (see project memory project-mdm-naming-convention):
#   Md    DHS only, observed coords as truth (weights collapsed to centre)
#   M_D   DHS only, jitter integrated
#   M_dm  DHS observed + MICS single point ~ integration weight (no jitter)
#   M_DM  DHS+MICS, jitter integrated (standard shared effects)
#   M_M   MICS only
# (The commensurate and fully-non-commensurate/separate-effects M_DM variants were
#  dropped: the DHS and MICS effect estimates agree in the application, so there is
#  no need to model the MICS effects as either commensurate-shrunk or separate.)

VALMODELS <- c("Md","M_D","M_dm","M_DM","M_M")

# Which survey(s) a model is FIT on (used by the areal fit and the combined-fit
# subset): DHS-only models drop MICS, MICS-only drops DHS. (Every model now
# PREDICTS both surveys' held-out clusters; this is only about what it's fit on.)
valModelValid <- function(model, side) {
  switch(model,
         Md   = side == "DHS",
         M_D  = side == "DHS",
         M_M  = side == "MICS",
         TRUE)            # M_dm, M_DM use both
}

# Enumerate every validation job. The CV is 20-fold (10 DHS + 10 MICS held-out
# folds) for all models, collapsed for single-dataset models: a survey NOT used in
# fitting contributes ONE job (fold = NA) predicting ALL of it, since its 10 folds
# share the same fit. Per model: predict-DHS jobs (10 folds if DHS-fit, else 1
# all-DHS job) + predict-MICS jobs (10 if MICS-fit, else 1). Then one areal job
# per model. Counts: Md/M_D/M_M = 11 cluster jobs each, M_dm/M_DM = 20 each
# (73 cluster + 5 areal = 78 for the 5 VALMODELS).
.buildValJobList <- function(models = VALMODELS) {
  jobs <- list()
  add  <- function(...) jobs[[length(jobs) + 1]] <<- list(...)
  for(m in models) {
    if(.usesDHS(m))  for(k in 1:10) add(type="cluster", model=m, predict="DHS",  fold=k)
    else                            add(type="cluster", model=m, predict="DHS",  fold=NA_integer_)
    if(.usesMICS(m)) for(k in 1:10) add(type="cluster", model=m, predict="MICS", fold=k)
    else                            add(type="cluster", model=m, predict="MICS", fold=NA_integer_)
  }
  for(m in models) add(type="areal", model=m, predict=NA_character_, fold=NA_integer_)
  jobs
}

# Decode a SLURM array index into one job from .buildValJobList().
decodeValJob <- function(index, models = VALMODELS) {
  jobs <- .buildValJobList(models)
  if(index < 1 || index > length(jobs))
    return(list(type="invalid", isValid=FALSE, tag=paste0("invalid_", index)))
  j <- jobs[[index]]; j$isValid <- TRUE
  j$tag <- if(j$type == "areal") paste0(j$model, "_areal")
           else paste0(j$model, "_pred", j$predict,
                       if(is.na(j$fold)) "_all" else paste0("_fold", j$fold))
  j
}

# Apply a model's positional-uncertainty treatment to an inputs object. Md/M_dm
# collapse DHS to the observed centre; M_dm additionally draws one MICS point per
# cluster proportional to its integration weight (makeInputsMdm). All others keep
# the full jitter integration. Used identically for in-sample inputs AND held-out
# inputs (so the held-out prediction matches the model's assumption).
applyPositionalTreatment <- function(inp, model, seed = 1) {
  if(model == "Md") {
    # observed coords as truth: all DHS weight on the centre (integration pt 1).
    # (defined inline so we don't depend on the script-local makeMdInputs)
    inp$intPtsDHS$wUrban[, -1] <- 0; inp$intPtsDHS$wUrban[, 1] <- 1
    inp$intPtsDHS$wRural[, -1] <- 0; inp$intPtsDHS$wRural[, 1] <- 1
    inp
  } else if(model == "M_dm") {
    makeInputsMdm(inp, seed = seed)   # DHS centre + MICS single point ~ weight
  } else {
    inp
  }
}

# Fit one FE-family model (FE+nug, no spatial field) on already-subset inputs.
fitValModelFE <- function(model, fullInp, inSampleInp,
                          KMICS = fullInp$KMICS, KDHSu = fullInp$KDHSurb,
                          KDHSr = fullInp$KDHSrur, Qgh = 10, verbose = FALSE) {
  dD <- fullInp$datDHS; dM <- fullInp$datMICS
  inp <- applyPositionalTreatment(inSampleInp, model)
  switch(model,
    Md   = fitFED (datDHS=dD,            inputsMDM=inp, KDHSurb=KDHSu, KDHSrur=KDHSr,
                   Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_D  = fitFED (datDHS=dD,            inputsMDM=inp, KDHSurb=KDHSu, KDHSrur=KDHSr,
                   Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_M  = fitFEM (datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS,
                   Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_dm = fitFEMD(datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                   Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_DM = fitFEMD(datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                   Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    stop("fitValModelFE: unknown model ", model))
}

# Fit one BYM2-family model (spatial field + nugget) on already-subset inputs --
# the spatial analogue of fitValModelFE, dispatching to the application's BYM2
# fitters (fitMd/fitMD/fitMDM/fitMM). Same positional treatment as the FE path
# (Md observed coords; M_dm DHS centre + MICS single point). fitMd has no jitter
# integration (observed coords) so takes no Qgh; the others integrate over jitter.
fitValModelBYM2 <- function(model, fullInp, inSampleInp,
                            KMICS = fullInp$KMICS, KDHSu = fullInp$KDHSurb,
                            KDHSr = fullInp$KDHSrur, Qgh = 10, verbose = FALSE) {
  dD <- fullInp$datDHS; dM <- fullInp$datMICS
  inp <- applyPositionalTreatment(inSampleInp, model)
  switch(model,
    Md   = fitMd (datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                  getSDs=TRUE, verbose=verbose),
    M_D  = fitMD (datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                  Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_dm = fitMDM(datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                  Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_DM = fitMDM(datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                  Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    M_M  = fitMM (datDHS=dD, datMICS=dM, inputsMDM=inp, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                  Qgh=Qgh, getSDs=TRUE, verbose=verbose),
    stop("fitValModelBYM2: unknown model ", model))
}

# Pick the fitter for a family ("FE" or "BYM2").
fitValModelFamily <- function(model, fullInp, inSampleInp, family = c("FE","BYM2"), ...) {
  family <- match.arg(family)
  if(family == "FE") fitValModelFE(model, fullInp, inSampleInp, ...)
  else               fitValModelBYM2(model, fullInp, inSampleInp, ...)
}

# Per-fold CLUSTER-level validation for one (model, fold, side). Subsets to
# in-sample (fold held out), fits, predicts the held-out clusters conditional on
# their observed locations, scores. fold: 1..10. side: "DHS" or "MICS".
# Returns the tag; writes fit/preds/scores under outDir.
runValFoldFE <- function(model, fold, side, fullInp,
                         outDir = "savedOutput/validation/folds",
                         nsim = 10000, KMICS = fullInp$KMICS,
                         KDHSu = fullInp$KDHSurb, KDHSr = fullInp$KDHSrur, Qgh = 10) {
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  tag <- paste0(model, if(side=="DHS") "_dhsFold" else "_micsFold", fold)
  nDHS <- nrow(fullInp$datDHS); nMICS <- nrow(fullInp$datMICS)
  if(side == "DHS") {
    dhsKeep  <- fullInp$datDHS$fold != fold;  micsKeep <- rep(TRUE, nMICS)
    dhsHeld  <- !dhsKeep;                      micsHeld <- rep(FALSE, nMICS)
  } else {
    dhsKeep  <- rep(TRUE, nDHS);               micsKeep <- fullInp$datMICS$fold != fold
    dhsHeld  <- rep(FALSE, nDHS);              micsHeld <- !micsKeep
  }
  inSamp <- subsetMDMInputs(fullInp, dhsKeep, micsKeep)
  held   <- subsetMDMInputs(fullInp, dhsHeld, micsHeld)

  t0 <- proc.time()[3]
  fitRes <- fitValModelFE(model, fullInp, inSamp, KMICS=KMICS, KDHSu=KDHSu, KDHSr=KDHSr, Qgh=Qgh)
  totalTime <- proc.time()[3] - t0
  save(fitRes, totalTime, file = file.path(outDir, paste0("fit", tag, ".RData")))

  # held-out prediction uses the model's positional treatment too
  heldT <- applyPositionalTreatment(held, model)
  if(side == "DHS") {
    nUrb <- heldT$nDHSurb; nRur <- heldT$nDHSrur
    oos  <- list(covsUrb=heldT$intPtsDHS$covsUrb, covsRur=heldT$intPtsDHS$covsRur,
                 wUrban =heldT$intPtsDHS$wUrban,  wRural =heldT$intPtsDHS$wRural)
    ysU<-heldT$ysUrbDHS; nsU<-heldT$nsUrbDHS; ysR<-heldT$ysRurDHS; nsR<-heldT$nsRurDHS
  } else {
    nUrb <- heldT$nMICSurb; nRur <- heldT$nMICSrur
    oos  <- list(covsUrb=heldT$intPtsMICS$XUrb, covsRur=heldT$intPtsMICS$XRur,
                 wUrban =heldT$intPtsMICS$wUrban, wRural =heldT$intPtsMICS$wRural)
    ysU<-heldT$ysUrbMICS; nsU<-heldT$nsUrbMICS; ysR<-heldT$ysRurMICS; nsR<-heldT$nsRurMICS
  }
  predObj <- predFEclusters(fitRes, oos, nUrb, nRur, nsim=nsim, side="DHS")
  formatAndSavePreds(predObj, ysU, nsU, ysR, nsR, fold=fold, tag=tag, totalTime=totalTime, outDir=outDir)
  scoreTagFold(tag, fold, outDir=outDir)
  invisible(tag)
}

# Out-of-sample prediction pieces (covariates + integration weights + held-out
# y/n) for one survey, pulled from a (treated) held-out inputs object.
.heldOOS <- function(heldT, survey) {
  if(survey == "DHS")
    list(oos = list(covsUrb=heldT$intPtsDHS$covsUrb, covsRur=heldT$intPtsDHS$covsRur,
                    wUrban =heldT$intPtsDHS$wUrban,  wRural =heldT$intPtsDHS$wRural),
         nUrb=heldT$nDHSurb, nRur=heldT$nDHSrur,
         yU=heldT$ysUrbDHS, nU=heldT$nsUrbDHS, yR=heldT$ysRurDHS, nR=heldT$nsRurDHS)
  else
    list(oos = list(covsUrb=heldT$intPtsMICS$XUrb, covsRur=heldT$intPtsMICS$XRur,
                    wUrban =heldT$intPtsMICS$wUrban, wRural =heldT$intPtsMICS$wRural),
         nUrb=heldT$nMICSurb, nRur=heldT$nMICSrur,
         yU=heldT$ysUrbMICS, nU=heldT$nsUrbMICS, yR=heldT$ysRurMICS, nR=heldT$nsRurMICS)
}

# Which surveys a model is FIT on.
.usesDHS  <- function(model) model %in% c("Md","M_D","M_dm","M_DM")
.usesMICS <- function(model) model %in% c("M_dm","M_DM","M_M")

# ONE cluster-level CV job: fit the model, then predict the held-out clusters of
# ONE survey. The full design is a 20-fold CV (10 DHS held-out folds + 10 MICS
# held-out folds) for EVERY model, but it collapses for the single-dataset models:
# when the predicted survey is NOT used in fitting, leaving out one of its folds
# changes neither the fit nor the predictions, so all 10 folds share one fit on
# the full other survey and jointly predict the entire survey -- i.e. one job
# (fold = NA) predicting ALL of that survey ("the 11th fold"). So:
#   * predictSurvey IS fit-used  -> leave-one-fold-out: hold out fold k, fit on
#     (that survey \ fold k) + the other fit-used survey in FULL, predict fold k.
#   * predictSurvey NOT fit-used -> fold = NA: fit on the full fit-used survey(s),
#     predict ALL clusters of predictSurvey (cross-survey prediction).
# Held-out positional treatment mirrors the original validation.R mapping:
# observed/single-point for lowercase models (Md -> DHS observed + MICS single
# point), full jitter/geomask integration for uppercase (M_D / M_DM / M_M).
# Saves the held-out probability draws + y/n so .scoreModelCluster can pool every
# fold of a model into per-survey and denominator-weighted combined scores.
runValJobCluster <- function(model, predictSurvey, fold, fullInp,
                             outDir = "savedOutput/validation/folds",
                             nsim = 10000, family = c("FE","BYM2"), regenerate = FALSE,
                             KMICS = fullInp$KMICS,
                             KDHSu = fullInp$KDHSurb, KDHSr = fullInp$KDHSrur, Qgh = 10) {
  family <- match.arg(family)
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  tag <- paste0(model, "_pred", predictSurvey, if(is.na(fold)) "_all" else paste0("_fold", fold))
  outFile <- file.path(outDir, paste0("predsJob_", tag, ".RData"))
  if(!regenerate && file.exists(outFile)) return(invisible(tag))   # skip completed job (re-collate cheaply)
  nDHS <- nrow(fullInp$datDHS); nMICS <- nrow(fullInp$datMICS)
  uD <- .usesDHS(model); uM <- .usesMICS(model)
  predIsFit <- (predictSurvey == "DHS" && uD) || (predictSurvey == "MICS" && uM)

  # in-sample: full fit-used survey(s); if predicting a fit-used survey, hold out
  # its fold k (the OTHER fit-used survey stays in full).
  dhsKeep  <- rep(uD, nDHS); micsKeep <- rep(uM, nMICS)
  if(predIsFit) {
    if(predictSurvey == "DHS") dhsKeep  <- fullInp$datDHS$fold  != fold
    else                       micsKeep <- fullInp$datMICS$fold != fold
  }
  inSamp <- subsetMDMInputs(fullInp, dhsKeep, micsKeep)
  t0 <- proc.time()[3]
  fitRes <- fitValModelFamily(model, fullInp, inSamp, family=family, KMICS=KMICS, KDHSu=KDHSu, KDHSr=KDHSr, Qgh=Qgh)
  totalTime <- proc.time()[3] - t0

  # held-out set: fold k of predictSurvey if fit-used, else ALL of predictSurvey.
  if(predictSurvey == "DHS")
    held <- subsetMDMInputs(fullInp,
                            if(predIsFit) fullInp$datDHS$fold == fold else rep(TRUE, nDHS),
                            rep(FALSE, nMICS))
  else
    held <- subsetMDMInputs(fullInp, rep(FALSE, nDHS),
                            if(predIsFit) fullInp$datMICS$fold == fold else rep(TRUE, nMICS))

  treatAs <- if(predictSurvey == "DHS") (if(model %in% c("Md","M_dm")) "Md" else "M_D")
             else                       (if(model %in% c("Md","M_dm")) "M_dm" else "M_DM")
  heldT <- applyPositionalTreatment(held, treatAs)
  z  <- .heldOOS(heldT, predictSurvey)
  # held-out clusters' nominal area indices (0-based) for the BYM2 field term.
  areaidx <- if(predictSurvey == "DHS")
    list(urb=heldT$areaidxlocUrbanDHS,  rur=heldT$areaidxlocRuralDHS)
  else
    list(urb=heldT$areaidxlocUrbanMICS, rur=heldT$areaidxlocRuralMICS)
  pd <- predFEclusters(fitRes, z$oos, z$nUrb, z$nRur, nsim = nsim, side = "DHS", areaidx = areaidx)
  preds <- list(probDrawsUrb=pd$probDrawsUrb, probDrawsRur=pd$probDrawsRur,
                yUrb=z$yU, nUrb=z$nU, yRur=z$yR, nRur=z$nR,
                predictSurvey=predictSurvey)
  save(preds, model, predictSurvey, fold, totalTime, file = outFile)   # tag/outFile set at top
  invisible(tag)
}

# Pool ALL of one model's saved cluster jobs (predsJob_<model>_*) into per-survey
# (DHS, MICS) and a denominator-weighted COMBINED score. Each held-out cluster is
# predicted exactly once across the jobs (DHS folds 1..10 + MICS folds 1..10, or
# the collapsed all-survey job). Weighting is weights = n (cluster denominator /
# size, the original validation's getScores weight); pooling DHS+MICS clusters by
# n weights each survey relative to the other by its TOTAL denominator, so the
# non-comparable per-survey sampling weights never enter. Returns a data.frame
# with DHS / MICS / combined rows, or NULL if no job files are present.
.scoreModelCluster <- function(model, foldDir = "savedOutput/validation/folds",
                               significance = c(.5,.8,.9,.95)) {
  files <- list.files(foldDir, pattern = paste0("^predsJob_", model, "_pred.*\\.RData$"),
                      full.names = TRUE)
  if(length(files) == 0) return(NULL)
  acc <- list(DHS=list(truth=NULL, draws=NULL, n=NULL),
              MICS=list(truth=NULL, draws=NULL, n=NULL))
  tsec <- numeric(0)                                    # per-job fit times (seconds)
  for(f in files) {
    e <- new.env(); load(f, envir = e); p <- e$preds; s <- p$predictSurvey
    truth <- c(p$yUrb, p$yRur) / c(p$nUrb, p$nRur)
    draws <- rbind(p$probDrawsUrb, p$probDrawsRur)
    n     <- c(p$nUrb, p$nRur)
    acc[[s]]$truth <- c(acc[[s]]$truth, truth)
    acc[[s]]$draws <- rbind(acc[[s]]$draws, draws)
    acc[[s]]$n     <- c(acc[[s]]$n, n)
    if(!is.null(e$totalTime)) tsec <- c(tsec, e$totalTime)
  }
  tmin <- if(length(tsec)) mean(tsec, na.rm=TRUE)/60 else NA_real_   # mean fit time (minutes)
  scoreSet <- function(a) {
    if(is.null(a$truth) || length(a$truth) == 0) return(NULL)
    getScores(truth=a$truth, estMat=a$draws, weights=a$n, significance=significance,
              doFuzzyReject=TRUE, getAverage=TRUE, na.rm=TRUE, ns=a$n, addBinomVar=TRUE)
  }
  mk <- function(s, side) if(is.null(s)) NULL else
    data.frame(model=model, side=side, as.data.frame(as.list(s)), stringsAsFactors=FALSE)
  rowD <- mk(scoreSet(acc$DHS),  "DHS")
  rowM <- mk(scoreSet(acc$MICS), "MICS")
  # Combined = denominator-weighted average of the DHS and MICS score rows
  # (weights = total cluster denominator per survey). Computed directly from the
  # two per-survey rows -- NOT by re-pooling and re-scoring all clusters -- so the
  # combined score is GUARANTEED to lie between the DHS and MICS scores (never
  # above the max). Re-pooling could exceed the max because addBinomVar and the
  # fuzzy-reject coverage draw fresh randomness independently of the per-survey
  # calls, breaking the exact weighted-average identity.
  rowC <- NULL
  if(!is.null(rowD) && !is.null(rowM)) {
    wD <- sum(acc$DHS$n); wM <- sum(acc$MICS$n)
    rowC <- rowD; rowC$side <- "combined"
    numCols <- setdiff(names(rowD), c("model", "side"))
    rowC[numCols] <- (wD * rowD[numCols] + wM * rowM[numCols]) / (wD + wM)
  } else if(!is.null(rowD)) { rowC <- rowD; rowC$side <- "combined" }
    else if(!is.null(rowM)) { rowC <- rowM; rowC$side <- "combined" }
  # mean per-fit fit time (minutes), a model-level quantity (same across sides)
  for(nm in c("rowD","rowM","rowC")) if(!is.null(get(nm))) assign(nm, `[[<-`(get(nm), "Time", tmin))
  do.call(rbind, list(rowD, rowM, rowC))
}

# Admin1 svyglm Horvitz-Thompson direct estimates from the HELD-OUT folds, per
# survey and combined. Reuses getDirectEsts (the workhorse svyglm estimator) on
# the fullInp data (so the held-out folds match the model fit's in-sample folds),
# then forms the precision-weighted combined estimate on the LOGIT scale with NA
# handling (an area un-estimable in one survey falls back to the other; NA in both
# stays NA and is dropped at scoring). Returns list(DHS, MICS, combined), each a
# data.frame with columns area, logit.est, logit.var.
arealDirectEsts <- function(fullInp, leftOutFolds = 6:10,
                            signifs = c(.5, .8, .9, .95)) {
  dD <- fullInp$datDHS[fullInp$datDHS$fold %in% leftOutFolds, , drop=FALSE]
  names(dD)[names(dD) == "clusterID"] <- "clustID"
  estsDHS <- getDirectEsts(dD, signifs=signifs)

  dM <- fullInp$datMICS[fullInp$datMICS$fold %in% leftOutFolds, , drop=FALSE]
  names(dM)[names(dM) == "ys"]   <- "y"
  names(dM)[names(dM) == "ns"]   <- "n"
  names(dM)[names(dM) == "Area"] <- "area"
  estsMICS <- getDirectEsts(dM, signifs=signifs, customStratVarName="Stratum")

  # combine on the logit scale, aligning by area name
  areas <- sort(unique(c(estsDHS$area, estsMICS$area)))
  d <- estsDHS[match(areas, estsDHS$area), c("logit.est","logit.var")]
  m <- estsMICS[match(areas, estsMICS$area), c("logit.est","logit.var")]
  precD <- 1/d$logit.var; precM <- 1/m$logit.var
  precD[!is.finite(precD)] <- 0; precM[!is.finite(precM)] <- 0
  wsum <- precD + precM
  cl.est <- (ifelse(precD>0, d$logit.est, 0)*precD + ifelse(precM>0, m$logit.est, 0)*precM) / wsum
  cl.var <- 1 / wsum                       # precision-weighted combine, logit scale
  cl.est[wsum == 0] <- NA; cl.var[!is.finite(cl.var)] <- NA
  combined <- data.frame(area=areas, logit.est=cl.est, logit.var=cl.var,
                         stringsAsFactors=FALSE)
  list(DHS=estsDHS, MICS=estsMICS, combined=combined)
}

# AREAL validation for one model: PER-AREA leave-half-out, matching the original
# validation.R getValidationData*(areal=TRUE) (inSample = (area != foldArea) |
# (fold <= 5); outOfSample = (area == foldArea) & (fold > 5)). For EACH admin1
# area we hold out that area's leftOutFolds (6:10) clusters, fit the model on
# everything else (ALL other areas in full + this area's OTHER half), and predict
# just that area (predGrid predAtArea -> cheap, one area's pixels). Each area is
# then scored against the svyglm Horvitz-Thompson direct estimate built from its
# OWN held-out clusters (DHS, MICS, and precision-weighted combined). Scoring is
# Gascoigne-style: getScoresDirectEstimates simulates the direct estimate's
# sampling distribution and scores (Yhat - Y) against 0 over the areas.
#
# The direct estimates are unchanged by going per-area: each area's direct est
# depends only on its own held-out clusters, so arealDirectEsts(leftOutFolds)
# already gives the correct per-area "truths". What changes vs a single global
# fold-1:5 fit is the MODEL prediction -- now a proper per-area CV (the model
# never sees the held-out half of the area it predicts, but DOES use all other
# areas + that area's in-sample half, as in real small-area prediction).
# Returns a scores data.frame (combined / DHS / MICS direct-est rows).
runValArealFE <- function(model, fullInp,
                          outDir = "savedOutput/validation/areal",
                          leftOutFolds = 6:10, nsim = 10000,
                          popMat = popMatNGAThresh, adm1obj = adm1,
                          KMICS = fullInp$KMICS, KDHSu = fullInp$KDHSurb,
                          KDHSr = fullInp$KDHSrur, Qgh = 10, verbose = FALSE) {
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  tag <- paste0(model, "_areal_fold", paste(range(leftOutFolds), collapse="-"))

  # per-area direct-estimate "truths" (DHS / MICS / combined), unchanged by the
  # per-area fitting since each area's estimate uses only its own held-out clusters.
  de <- arealDirectEsts(fullInp, leftOutFolds = leftOutFolds)
  areas <- de$combined$area                      # admin1 areas (sorted), the scoring order

  nDHS <- nrow(fullInp$datDHS); nMICS <- nrow(fullInp$datMICS)
  uD <- .usesDHS(model); uM <- .usesMICS(model)
  estMat   <- matrix(NA_real_, nrow = length(areas), ncol = nsim)   # per-area model pred draws
  fitTimes <- numeric(length(areas))
  t0all <- proc.time()[3]
  for(ai in seq_along(areas)) {
    a <- areas[ai]
    # in-sample = everything EXCEPT this area's held-out folds (other areas full +
    # this area's in-sample half); fit-used surveys only.
    dhsKeep  <- if(uD) !((fullInp$datDHS$area  == a) & (fullInp$datDHS$fold  %in% leftOutFolds)) else rep(FALSE, nDHS)
    micsKeep <- if(uM) !((fullInp$datMICS$Area == a) & (fullInp$datMICS$fold %in% leftOutFolds)) else rep(FALSE, nMICS)
    inSamp <- subsetMDMInputs(fullInp, dhsKeep, micsKeep)
    t0 <- proc.time()[3]
    fitRes <- fitValModelFE(model, fullInp, inSamp, KMICS=KMICS, KDHSu=KDHSu, KDHSr=KDHSr, Qgh=Qgh, verbose=FALSE)
    fitTimes[ai] <- proc.time()[3] - t0
    grid <- predGrid(fitRes$TMBsd, popMat=popMat, nsim=nsim, obj=fitRes$TMBobj,
                     admLevel="stratMICS", useInla=FALSE, predAtArea=a)
    a1 <- predArea(grid, areaVarName="area", orderedAreas=a)
    p  <- as.matrix(a1$aggregationResults$p)
    if(nrow(p) >= 1) estMat[ai, ] <- p[1, ]
    if(verbose) cat(sprintf("  [areal %s] %2d/%d %-22s (%.2f min)\n",
                            model, ai, length(areas), a, fitTimes[ai]/60))
  }
  totalTime <- proc.time()[3] - t0all

  # score per-area model preds vs each direct-est truth (estMat rows aligned to `areas`)
  finiteRow <- apply(estMat, 1, function(r) all(is.finite(r)))
  scoreOne <- function(dEst, label) {
    le <- dEst$logit.est[match(areas, dEst$area)]
    lv <- dEst$logit.var[match(areas, dEst$area)]
    ok <- is.finite(le) & is.finite(lv) & finiteRow
    if(sum(ok) == 0) return(NULL)
    sc <- getScoresDirectEstimates(le[ok], lv[ok], estMat = estMat[ok, , drop=FALSE], na.rm=TRUE)
    data.frame(model=model, directEst=label, nAreas=sum(ok),
               as.data.frame(as.list(sc)), stringsAsFactors=FALSE)
  }
  scores <- do.call(rbind, list(scoreOne(de$combined, "combined"),
                                scoreOne(de$DHS,      "DHS"),
                                scoreOne(de$MICS,     "MICS")))
  save(de, estMat, areas, scores, fitTimes, totalTime, tag,
       file = file.path(outDir, paste0("areal_", tag, ".RData")))
  scores
}

# ONE area of the per-area leave-half-out areal validation (the parallelizable
# unit, so the ~37 areas of a model spread across workers instead of one long
# looping job). Holds out THIS area's leftOutFolds clusters, fits on everything
# else, predicts just this area (predGrid predAtArea), and saves the area's
# posterior prevalence draws. areaIdx indexes the canonical area order
# (arealDirectEsts()$combined$area) so scoreModelAreal can reassemble estMat.
runValArealOneArea <- function(model, areaIdx, area, fullInp,
                               outDir = "savedOutput/validation/areal",
                               leftOutFolds = 6:10, nsim = 2000, popMat = popMatNGAThresh,
                               family = c("FE","BYM2"), useInla = NULL, regenerate = FALSE,
                               KMICS = fullInp$KMICS, KDHSu = fullInp$KDHSurb,
                               KDHSr = fullInp$KDHSrur, Qgh = 10) {
  family <- match.arg(family)
  if(is.null(useInla)) useInla <- if(family == "BYM2") "auto" else FALSE  # Gaussian w/ INLA fallback
  if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
  tag <- sprintf("%s_area%03d", model, areaIdx)
  outFile <- file.path(outDir, paste0("arealPred_", tag, ".RData"))
  if(!regenerate && file.exists(outFile)) return(invisible(tag))   # skip completed area
  nDHS <- nrow(fullInp$datDHS); nMICS <- nrow(fullInp$datMICS)
  uD <- .usesDHS(model); uM <- .usesMICS(model)
  dhsKeep  <- if(uD) !((fullInp$datDHS$area  == area) & (fullInp$datDHS$fold  %in% leftOutFolds)) else rep(FALSE, nDHS)
  micsKeep <- if(uM) !((fullInp$datMICS$Area == area) & (fullInp$datMICS$fold %in% leftOutFolds)) else rep(FALSE, nMICS)
  inSamp <- subsetMDMInputs(fullInp, dhsKeep, micsKeep)
  t0 <- proc.time()[3]
  fitRes <- fitValModelFamily(model, fullInp, inSamp, family=family, KMICS=KMICS, KDHSu=KDHSu, KDHSr=KDHSr, Qgh=Qgh)
  # res=fitRes is required for useInla="auto"/"yes" (posteriorDraws needs the full
  # fit: TMBobj + TMBsd). Harmless for useInla=FALSE (the FE path).
  grid <- predGrid(fitRes$TMBsd, popMat=popMat, nsim=nsim, obj=fitRes$TMBobj, res=fitRes,
                   admLevel="stratMICS", useInla=useInla, predAtArea=area)
  a1 <- predArea(grid, areaVarName="area", orderedAreas=area)
  draws <- as.matrix(a1$aggregationResults$p)[1, ]
  fitTime <- proc.time()[3] - t0
  save(draws, model, area, areaIdx, fitTime, file = outFile)   # tag/outFile set at top
  invisible(tag)
}

# Collate one model's per-area predictions (arealPred_<model>_areaNNN) into
# DHS/MICS/combined areal scores. de (the per-area direct estimates) is the same
# for every model, so pass it in once; recomputed if NULL. Areas are ordered by
# de$combined$area to match runValArealOneArea's areaIdx.
scoreModelAreal <- function(model, fullInp, arealDir = "savedOutput/validation/areal",
                            leftOutFolds = 6:10, de = NULL) {
  if(is.null(de)) de <- arealDirectEsts(fullInp, leftOutFolds = leftOutFolds)
  areas <- de$combined$area
  files <- list.files(arealDir, pattern = paste0("^arealPred_", model, "_area[0-9]+\\.RData$"),
                      full.names = TRUE)
  if(length(files) == 0) return(NULL)
  drawsList <- vector("list", length(areas))
  tsec <- numeric(0)                                    # per-area fit times (seconds)
  for(f in files) { e <- new.env(); load(f, envir = e); drawsList[[e$areaIdx]] <- e$draws
                    if(!is.null(e$fitTime)) tsec <- c(tsec, e$fitTime) }
  tmin <- if(length(tsec)) mean(tsec, na.rm=TRUE)/60 else NA_real_   # mean fit time (minutes)
  ncols  <- max(sapply(drawsList, function(d) if(is.null(d)) 0L else length(d)))
  estMat <- matrix(NA_real_, nrow = length(areas), ncol = ncols)
  for(i in seq_along(drawsList)) if(!is.null(drawsList[[i]])) estMat[i, ] <- drawsList[[i]]
  finiteRow <- apply(estMat, 1, function(r) all(is.finite(r)))
  # Score DHS / MICS / combined on a COMMON area set (areas where ALL THREE direct
  # estimates are finite AND the model prediction is finite). Otherwise each is on
  # a different support -- the combined direct est is finite wherever EITHER survey
  # has data, so it covers more areas than DHS- or MICS-only, making its coverage
  # incomparable (and able to exceed both). A common set makes them comparable, so
  # the combined coverage falls between the two as expected.
  matchDE <- function(dEst, col) dEst[[col]][match(areas, dEst$area)]
  finOK   <- function(dEst) is.finite(matchDE(dEst, "logit.est")) & is.finite(matchDE(dEst, "logit.var"))
  okCommon <- finiteRow & finOK(de$DHS) & finOK(de$MICS) & finOK(de$combined)
  scoreOne <- function(dEst, label) {
    if(sum(okCommon) == 0) return(NULL)
    le <- matchDE(dEst, "logit.est")[okCommon]; lv <- matchDE(dEst, "logit.var")[okCommon]
    sc <- getScoresDirectEstimates(le, lv, estMat = estMat[okCommon, , drop=FALSE], na.rm=TRUE)
    data.frame(model=model, directEst=label, nAreas=sum(okCommon),
               as.data.frame(as.list(sc)), Time=tmin, stringsAsFactors=FALSE)
  }
  do.call(rbind, list(scoreOne(de$combined, "combined"),
                      scoreOne(de$DHS, "DHS"), scoreOne(de$MICS, "MICS")))
}
