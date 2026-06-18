# simulate data from an SPDE model for simulation study 1

# Add DHS-style positional jittering to a sampled cluster survey (output of
# SUMMER::sampleClusterSurveys). Each cluster's PUBLISHED location is displaced
# from its TRUE location by a uniform direction (0..2pi) and uniform radius
# (0..maxDist), with maxDist = 2km (urban) or 5km w.p. 0.99 / 10km w.p. 0.01
# (rural), per the DHS geomasking scheme (see dataFusionGeomasking.tex). Jitter
# must NOT cross the boundary of the area in `adminMap` containing the true
# location (default Admin2); this is enforced by rejection sampling. Coordinates
# are easting/northing in km, while `adminMap` is lon/lat, so a projection
# function `proj` (default projNigeria) maps between them.
#
# Efficiency: only true locations within maxDist of their own area boundary can
# ever produce a rejected jitter (cf. the DHS integration-point weight code,
# which zeroes cross-border points). We compute distance-to-boundary once
# (dist2Line on the point's own polygon) and rejection-sample ONLY those points;
# all others are jittered in a single accepted draw.
#
# Returns the survey with east/north/lon/lat replaced by the jittered (published)
# values and eastTrue/northTrue/lonTrue/latTrue holding the originals. The
# responses (Z, N, pFineScale*) are unchanged — they were generated at the true
# location; only the observed coordinates are anonymized.
addJitterToDHS = function(survey, adminMap=adm2Full, proj=projNigeria,
                          areaNameVar="NAME_2", areaCol="subarea", urbanCol="urban",
                          eastCol="east", northCol="north",
                          gridRes=5, maxItr=1000, verbose=FALSE) {
  # Jitter DHS published coords away from the true (pixel-center) location while
  # keeping them inside the cluster's NOMINAL area (the claimed `areaCol`, an
  # adm2/NAME_2 label) -- the SAME area the integration points and the model's
  # areal effect use. We deliberately do NOT re-derive the area from over(truePt):
  # a true center can fall in a different polygon than its claimed label (snapping
  # near borders), and constraining to that geographic neighbor would associate
  # the cluster with the wrong areal effect.
  #   - true point INSIDE its nominal polygon (the rule): jitter uniformly in
  #     angle/radius (urban 2km; rural 5km w.p. .99, 10km w.p. .01), rejection-
  #     sampling so the draw stays in the nominal polygon.
  #   - true point OUTSIDE its nominal polygon (snapping artifact): OPTION 2 --
  #     set the published location to the closest point inside the nominal area
  #     (project onto its border, nudge just inside). Guarantees in-area without
  #     ever assigning a wrong area, and bounded by the snapping error.
  # gridRes: popMat cell size (km); only caps the inward search for option 2.
  require(geosphere)
  n = nrow(survey)
  east = survey[[eastCol]]; north = survey[[northCol]]
  urban = as.logical(survey[[urbanCol]])
  nomArea = as.character(survey[[areaCol]])
  nomID = match(nomArea, as.character(adminMap[[areaNameVar]]))
  if(anyNA(nomID))
    stop(sprintf("addJitterToDHS: %d clusters have a nominal %s not found in adminMap$%s",
                 sum(is.na(nomID)), areaCol, areaNameVar))

  # per-cluster max jitter distance; sample the rural 99/1 (5km/10km) up front
  maxDist = ifelse(urban, 2, ifelse(runif(n) < 0.99, 5, 10))

  trueLL = proj(cbind(east, north), inverse=TRUE)
  spLL   = SpatialPoints(trueLL, adminMap@proj4string)
  adminPolys = as.SpatialPolygons.PolygonsList(adminMap@polygons, adminMap@proj4string)

  # is each true point inside its OWN nominal polygon?
  ovName = as.character(over(spLL, adminMap, returnList=FALSE)[[areaNameVar]])
  inNom  = !is.na(ovName) & (ovName == nomArea)
  # distance (km) from each true point to its nominal polygon (border)
  distNom = sapply(seq_len(n), function(i) dist2Line(spLL[i], adminPolys[nomID[i]])[1]) / 1000

  jE = east; jN = north

  # ---- points INSIDE their nominal polygon ----
  # "free" = at least maxDist inside the border: a single draw can never cross.
  freeMask = inNom & (distNom >= maxDist)
  free = which(freeMask)
  if(length(free)) {
    th = runif(length(free), 0, 2*pi); r = runif(length(free), 0, maxDist[free])
    jE[free] = east[free] + r*cos(th); jN[free] = north[free] + r*sin(th)
  }
  # inside-but-near-border: rejection-sample within the nominal polygon (always
  # feasible -- the true point itself is inside, so small radii are accepted).
  need = inNom & !freeMask
  for(it in seq_len(maxItr)) {
    idx = which(need); if(!length(idx)) break
    th = runif(length(idx), 0, 2*pi); r = runif(length(idx), 0, maxDist[idx])
    ce = east[idx] + r*cos(th); cn = north[idx] + r*sin(th)
    a  = as.character(over(SpatialPoints(proj(cbind(ce, cn), inverse=TRUE),
                                         adminMap@proj4string), adminMap)[[areaNameVar]])
    ok = !is.na(a) & (a == nomArea[idx])
    jE[idx[ok]] = ce[ok]; jN[idx[ok]] = cn[ok]; need[idx[ok]] = FALSE
    if(verbose) cat(sprintf("  [jitter] iter %d: %d near-border points remaining\n", it, sum(need)))
  }
  if(any(need))
    warning(sprintf("addJitterToDHS: %d in-area clusters not accepted within %d iters; kept true location",
                    sum(need), maxItr))

  # ---- points OUTSIDE their nominal polygon: OPTION 2 (closest in-area point) ----
  outIdx = which(!inNom)
  nKeptTrue = 0L
  for(i in outIdx) {
    d2l = dist2Line(spLL[i], adminPolys[nomID[i]])   # [1]=dist(m) [2]=lon [3]=lat of nearest border pt
    bEN = proj(matrix(c(d2l[2], d2l[3]), nrow=1))    # nearest border point in east/north
    bE = bEN[1]; bN = bEN[2]
    # inward direction = from the (outside) true point toward the border point;
    # continuing past the border point steps INTO the nominal polygon.
    dE = bE - east[i]; dN = bN - north[i]; len = sqrt(dE^2 + dN^2)
    placed = FALSE
    if(len > 1e-9) {
      ux = dE/len; uy = dN/len
      step = 0.05; cap = maxDist[i] + sqrt(2)*gridRes
      while(step <= cap) {
        ce = bE + step*ux; cn = bN + step*uy
        a = as.character(over(SpatialPoints(proj(cbind(ce, cn), inverse=TRUE),
                                            adminMap@proj4string), adminMap)[[areaNameVar]])
        if(!is.na(a) && a == nomArea[i]) { jE[i] = ce; jN[i] = cn; placed = TRUE; break }
        step = step * 2
      }
    }
    if(!placed) { jE[i] = east[i]; jN[i] = north[i]; nKeptTrue = nKeptTrue + 1L }  # fallback: keep true loc
  }
  if(nKeptTrue)
    warning(sprintf("addJitterToDHS: %d outside-nominal clusters could not be projected in-area; kept true location",
                    nKeptTrue))

  survey$eastTrue = east; survey$northTrue = north
  survey[[eastCol]] = jE; survey[[northCol]] = jN
  if(all(c("lon","lat") %in% names(survey))) {
    newLL = proj(cbind(jE, jN), inverse=TRUE)
    survey$lonTrue = survey$lon; survey$latTrue = survey$lat
    survey$lon = newLL[,1]; survey$lat = newLL[,2]
  }
  attr(survey, "nOutsideNominal") = length(outIdx)
  attr(survey, "nKeptTrue") = nKeptTrue
  survey
}

# Replace each sampled cluster's pixel-CENTER coordinate with a uniform draw
# WITHIN its pixel, so the generated "true" locations are continuous (matching
# the continuous uniform-true-location assumption the jitter/integration model
# makes) rather than snapped to 5km grid centers. The only consistency guard is
# that each draw stays attributed to the SAME pixel it was sampled from (its
# nearest popMat pixel == survey$pixelIs); that pixel's claimed area/field are
# what generated the response, so they remain the ones the location maps to. We
# do NOT require the draw (or center) to fall in the pixel's claimed Admin2
# polygon -- the claimed label stands regardless of the geographic over(). The
# nearest-pixel test uses the grid: pre-filter popMat to pixels within one grid
# step in east AND north, then take the min-distance pixel among those.
# Coordinates are east/north km. Run BEFORE addJitterToDHS.
sampleTrueWithinPixel = function(survey, popMat, proj=projNigeria, maxItr=1000, verbose=FALSE) {
  if(is.null(survey$pixelIs)) stop("sampleTrueWithinPixel: survey lacks $pixelIs")
  if(is.null(popMat$lon) || is.null(popMat$lat))
    stop("sampleTrueWithinPixel: popMat needs lon/lat (the grid is regular in lon/lat, not east/north)")
  plon = popMat$lon; plat = popMat$lat; pe = popMat$east; pn = popMat$north
  pixI = survey$pixelIs
  lon0 = plon[pixI]; lat0 = plat[pixI]
  # cell half-widths in the NATIVE lon/lat grid (median consecutive spacing is
  # robust to projection irregularity and to grid gaps; min() is not).
  dlon = median(diff(sort(unique(round(plon, 6)))))
  dlat = median(diff(sort(unique(round(plat, 6)))))
  hLon = dlon/2; hLat = dlat/2
  n = nrow(survey); tLon = lon0; tLat = lat0; need = rep(TRUE, n)
  for(it in seq_len(maxItr)) {
    idx = which(need); if(!length(idx)) break
    clon = lon0[idx] + runif(length(idx), -hLon, hLon)
    clat = lat0[idx] + runif(length(idx), -hLat, hLat)
    en = proj(cbind(clon, clat))                 # lon/lat -> east/north (km)
    ce = en[,1]; cn = en[,2]
    # nearest popMat pixel (by east/north km distance) == attributed pixel;
    # pre-filter to one cell step in lon AND lat so the min-dist scan is tiny.
    nearOK = vapply(seq_along(idx), function(k) {
      nb = which(abs(plon - clon[k]) < dlon & abs(plat - clat[k]) < dlat)
      nb[which.min((pe[nb]-ce[k])^2 + (pn[nb]-cn[k])^2)] == pixI[idx[k]]
    }, logical(1))
    tLon[idx[nearOK]] = clon[nearOK]; tLat[idx[nearOK]] = clat[nearOK]
    need[idx[nearOK]] = FALSE
    if(verbose) cat(sprintf("  [withinPixel] iter %d: %d remaining\n", it, sum(need)))
  }
  if(any(need))
    warning(sprintf("sampleTrueWithinPixel: %d clusters kept at pixel center (rejection cap)", sum(need)))
  en = proj(cbind(tLon, tLat))                   # final true locations -> east/north
  survey[["east"]] = en[,1]; survey[["north"]] = en[,2]
  survey$lon = tLon; survey$lat = tLat
  survey
}

simData1 = function(nsim=100, margVar=.5, effRange=200, sigmaEpsilon=sqrt(1.5),
                    beta0=-1.25, gamma=1, betaRest=c(0, 0, 0, .5), 
                    mesh=getSPDEMesh(), easpaDat=easpaNGAed, 
                    popMat=popMatNGAThresh, targetPopMat=popMatNGAedThresh, 
                    poppsub=poppsubNGAThresh, nHHMICS=16, nHHDHS=25, seed=123, 
                    useThreshPopMat=TRUE, fixPopPerHH=NULL, 
                    eaSampleStrat="pps", regenPop=FALSE) {
  set.seed(seed)

  # make sure everything is ordered nicely. CRITICAL: targetPopMat must be
  # permuted with the SAME row order as popMat — SUMMER's simPop machinery
  # assumes the two are row-parallel (it weights pixel-level smooth risk by
  # targetPopMat$pop). Sorting only popMat misaligned 99.6% of rows and
  # silently flattened the smooth-risk truth weights to ~equal weighting
  # (~0.30 instead of ~0.47 nationally) — the source of the +0.11 areal
  # "bias" of (correct) predictions against the (corrupted) truth.
  stopifnot(nrow(targetPopMat) == nrow(popMat))
  .perm = order(popMat$subarea)
  popMat = popMat[.perm,]
  targetPopMat = targetPopMat[.perm,]
  stopifnot(all(popMat$subarea == targetPopMat$subarea))
  poppsub = poppsub[order(poppsub$subarea),]

  # construct logit offset vector based on covariates in betaRest
  # first get the design matrix
  print("Constructing offset based on covariates...")
  LLcoords = cbind(popMat$lon, popMat$lat)
  tempDesMat = getDesignMat(LLcoords, normalized=TRUE, useThreshPopMat=useThreshPopMat)
  # cbind(int=1, pop=popVals, urb=urbVals, access=accessVals, elev=elevVals, 
  #       distRiversLakes=distVals, urbanicity=urbanicityVals)
  
  # normalized population density is calculated differently
  load("savedOutput/global/covariatesNorm.RData")
  popVals = extract(pop, LLcoords, method="bilinear")
  
  load("savedOutput/global/popMeanSDCal.RData")
  popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
  popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
  normPop=(log1p(popVals)-popMeanCal)/popSDCal
  normPop[is.na(normPop)] = min(normPop, na.rm=TRUE)
  
  # get final design matrix
  covRestMat = tempDesMat[,-c(1:3, 7)] # remove int, pop, urb, urbanicity
  covRestMat = cbind(covRestMat, normPop=normPop) # add in normalized population density
  
  # calculate offset
  offset = covRestMat %*% betaRest
  
  # get aggregation info from admin2 areas to MICS strata
  tempAreasFrom = popMat$subarea
  tempAreasTo = popMat$stratumMICS
  areasFrom = sort(unique(tempAreasFrom))
  areasToI = match(areasFrom, tempAreasFrom)
  areasTo = tempAreasTo[areasToI]
  
  # simulate populations and surveys
  outFile = "savedOutput/simStudy1/simPopsSurveys_SPDE.RData"

  # Default initial state (fresh run); may be overridden by checkpoint resume.
  surveysDHS  = list()
  surveysMICS = list()
  subareaPops = NULL
  areaPops    = NULL
  stratumPops = NULL
  subareaPops_smoothRisk = NULL
  areaPops_smoothRisk    = NULL
  stratumPops_smoothRisk = NULL
  startI = 1

  # Try to resume from a compatible checkpoint (same seed, <= nsim sims saved).
  if(file.exists(outFile)) {
    chk = new.env()
    try(load(outFile, envir = chk), silent = TRUE)
    canResume <- exists("surveysDHS",      envir = chk) &&
                 exists(".rngStateAtStart", envir = chk) &&
                 exists(".seedUsed",       envir = chk) &&
                 isTRUE(chk$.seedUsed == seed) &&
                 length(chk$surveysDHS) <= nsim
    if(canResume) {
      surveysDHS  <- chk$surveysDHS
      surveysMICS <- chk$surveysMICS
      subareaPops <- chk$subareaPops
      areaPops    <- chk$areaPops
      stratumPops <- chk$stratumPops
      # Smooth-risk fields may not exist in older checkpoints; restore if present.
      if(exists("subareaPops_smoothRisk", envir = chk)) subareaPops_smoothRisk <- chk$subareaPops_smoothRisk
      if(exists("areaPops_smoothRisk",    envir = chk)) areaPops_smoothRisk    <- chk$areaPops_smoothRisk
      if(exists("stratumPops_smoothRisk", envir = chk)) stratumPops_smoothRisk <- chk$stratumPops_smoothRisk
      startI      <- length(surveysDHS) + 1
      if(startI > nsim) {
        print(paste0("All ", nsim, " sims already complete in ", outFile, "; skipping."))
        return(invisible(NULL))
      }
      assign(".Random.seed", chk$.rngStateAtStart, envir = .GlobalEnv)
      print(paste0("Resuming SPDE sims from ", startI, "/", nsim,
                   " (", startI - 1, " already in ", outFile, ")"))
    } else {
      print("Existing checkpoint is incompatible (different seed or no metadata); starting fresh.")
    }
  }

  print("Simulating populations and surveys...")
  startT = proc.time()[3]
  for(i in startI:nsim) {
    # simulate population at pixel, EA levels 
    print(paste0("Simulating population ", i, "/", nsim))

    simPop =
      SUMMER::simPopSPDE(nsim=1, easpa=easpaDat, popMat=popMat, targetPopMat=targetPopMat,
                         poppsub=poppsub, spdeMesh=mesh,
                         margVar=margVar, sigmaEpsilon=sigmaEpsilon, effRange=effRange,
                         gamma=gamma, beta0=beta0, seed=NULL, nHHSampled=nHHSampled,
                         stratifyByUrban=TRUE, subareaLevel=TRUE, offset=offset,
                         doFineScaleRisk=FALSE, doSmoothRisk=TRUE, min1PerSubarea=TRUE,
                         verbose=FALSE
      )

    # calculate stratum level population information
    stratPop = SUMMER::areaPopToArea(areaLevelPop=simPop$subareaPop,
                                     areasFrom=areasFrom,
                                     areasTo=areasTo,
                                     stratifyByUrban=TRUE, doFineScaleRisk=FALSE, doSmoothRisk=TRUE)

    # append population information.  We store BOTH the fine-scale prevalence
    # (binomial-realised Z_EA / N_EA, the historical truth used for scoring)
    # AND the smooth (pixel-level expected) risk, since predGrid() produces
    # smooth risk predictions and the consistent comparison is at that level.
    if(is.null(subareaPops)) {
      subareaPops      = simPop$subareaPop$aggregationResults$pFineScalePrevalence
      areaPops         = simPop$areaPop$aggregationResults$pFineScalePrevalence
      stratumPops      = stratPop$aggregationResults$pFineScalePrevalence
      subareaPops_smoothRisk = simPop$subareaPop$aggregationResults$pSmoothRisk
      areaPops_smoothRisk    = simPop$areaPop$aggregationResults$pSmoothRisk
      stratumPops_smoothRisk = stratPop$aggregationResults$pSmoothRisk
    } else {
      # cbind the new pop info to the full set of populations
      subareaPops = cbind(subareaPops,
                          simPop$subareaPop$aggregationResults$pFineScalePrevalence)
      areaPops    = cbind(areaPops,
                          simPop$areaPop$aggregationResults$pFineScalePrevalence)
      stratumPops = cbind(stratumPops,
                          stratPop$aggregationResults$pFineScalePrevalence)
      subareaPops_smoothRisk = cbind(subareaPops_smoothRisk,
                                     simPop$subareaPop$aggregationResults$pSmoothRisk)
      areaPops_smoothRisk    = cbind(areaPops_smoothRisk,
                                     simPop$areaPop$aggregationResults$pSmoothRisk)
      stratumPops_smoothRisk = cbind(stratumPops_smoothRisk,
                                     stratPop$aggregationResults$pSmoothRisk)
    }
    
    # generate surveys
    print(paste0("Generating surveys for population ", i, "/", nsim))
    # get EA level population information for population i
    thisEApop = simPop$eaPop$eaDatList[1]
    
    # get associated HH level population information
    thisHHpop = SUMMER::getHHpop(thisEApop, fixPopPerHH=fixPopPerHH)
    
    # sample DHS survey for this population
    survDHS = SUMMER::sampleClusterSurveys(1, thisHHpop, HHperClust=nHHDHS, clustpaList=list(clustpaDHSed))
    # Anonymize DHS cluster coordinates by jittering (true -> published) within
    # Admin2; responses (generated at the true location) are unchanged.
    # (sampleTrueWithinPixel for continuous within-pixel true locations is a
    # PARKED long-term faithfulness item -- not applied here; see memory.)
    survDHS[[1]] = addJitterToDHS(survDHS[[1]])
    
    # now sample the MICS survey. Do some gymnastics to make sure it works for MICS strata
    tempClustpa = clustpaMICSed
    names(tempClustpa)[1] = "area"
    
    thisHHpop[[1]]$area = adm2ToStratumMICS(thisHHpop[[1]]$subarea)
    
    survMICS = SUMMER::sampleClusterSurveys(1, thisHHpop, HHperClust=nHHMICS, clustpaList=list(tempClustpa))
    
    if(FALSE) {
      pixelIs = c(survDHS[[1]]$pixelIs, survMICS[[1]]$pixelIs)
      squilt(c(survDHS[[1]]$east, survMICS[[1]]$east), 
             c(survDHS[[1]]$north, survMICS[[1]]$north), 
             c(survDHS[[1]]$pFineScalePrevalence, survMICS[[1]]$pFineScalePrevalence), 
             nx=70, ny=70)
      
      plot(normPop[pixelIs], 
           c(survDHS[[1]]$pFineScalePrevalence, survMICS[[1]]$pFineScalePrevalence), 
           pch=19, col="blue", cex=.3)
      
      ps = c(survDHS[[1]]$pFineScalePrevalence, survMICS[[1]]$pFineScalePrevalence)
      pPops = normPop[pixelIs]
      urbs = c(survDHS[[1]]$urban, survMICS[[1]]$urban)
      summary(lm(ps ~ urbs + pPops))
      
      binMat = cbind(c(survDHS[[1]]$Z, survMICS[[1]]$Z), 
                     c(survDHS[[1]]$N, survMICS[[1]]$N) - c(survDHS[[1]]$Z, survMICS[[1]]$Z))
      summary(glm(binMat ~ urbs + pPops, family=binomial()))
      
      summary(lm(simPop$logitRiskDraws ~ popMat$urban + normPop))
    }
    
    # concatenate results
    surveysDHS = c(surveysDHS, survDHS)
    surveysMICS = c(surveysMICS, survMICS)

    # Save checkpoint atomically (tmpfile + rename) so a crash mid-save can't
    # leave a corrupt outFile. Stored metadata lets the next call resume.
    .rngStateAtStart <- get(".Random.seed", envir = .GlobalEnv)
    .seedUsed        <- seed
    .nsimRequested   <- nsim
    tmpFile <- paste0(outFile, ".tmp")
    save(subareaPops, areaPops, stratumPops,
         subareaPops_smoothRisk, areaPops_smoothRisk, stratumPops_smoothRisk,
         surveysDHS, surveysMICS,
         .rngStateAtStart, .seedUsed, .nsimRequested,
         file = tmpFile)
    file.rename(tmpFile, outFile)

    # estimate time left and print
    thisT = proc.time()[3]
    timePerIter = (thisT - startT)/(i - startI + 1)
    timeLeft = timePerIter * (nsim - i)

    print(paste0("estimated time remaining: ", (timeLeft/60)/24, " hours"))
  }

  invisible(NULL)
}

# simulates household level population data given EA level population data
# 
# 
getHHpop = function(popSim, fixPopPerHH=NULL, verbose=TRUE) {
  # first get ea level data from popSim
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      eaPopDat = popSim
    } else {
      stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
    }
  } else {
    stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  if((length(eaPopDat) > 1) && verbose) {
    warning("length(eaPopDat) > 1, so there may be a lot of household level data, and memory issues accordingly...")
  }
  
  HHdat = list()
  for(i in 1:length(eaPopDat)) {
    eaDat = eaPopDat[[i]]
    
    # first get the number of households
    numHouseholds = eaDat$nHH
    
    # now expand the eaDat table to be in the long format, where each row is a house
    rowsLong = rep(1:nrow(eaDat), numHouseholds)
    thisHHdat = eaDat[rowsLong, ]
    thisHHdat$eaIs = rowsLong
    urbanHH = thisHHdat$urban
    
    # function for randomly spreading people among households in long form data:
    extendThisDat = function(xs, nHH) {
      revenMultinom = function(sizeK) {
        size = sizeK[1]
        k = sizeK[2]
        prob = rep(1/k, k)
        rmultinom(1, size, prob)
      }
      unlist(apply(cbind(xs, nHH), 1, revenMultinom))
    }
    
    # generate how many of the target population are in each cluster
    if(is.null(fixPopPerHH)) {
      lived = extendThisDat(eaDat$N - eaDat$Z, numHouseholds)
      died = extendThisDat(eaDat$Z, numHouseholds)
    } else if(fixPopPerHH == 1) {
      extendDatEven = function(Ns, Zs) {
        
        spreadAmongHHs = function(thisRow) {
          thisN = thisRow[1]
          thisZ = thisRow[2]
          
          # spread population evenly among households
          # nHH = thisRow$nHH
          # hhI = sample(1:nHH, nHH, replace=FALSE)
          c(rep(1, thisZ), rep(0, thisN - thisZ))
        }
        c(unlist(apply(cbind(Ns, Zs), 1, spreadAmongHHs)))
      }
      died = extendDatEven(eaDat$N, eaDat$Z)
      lived = 1 - died
    } else {
      stop("If fixPopPerHH is not NULL it must be 1. Other values not currently supported")
    }
    
    thisHHdat$N = died + lived
    thisHHdat$Z = died
    thisHHdat$nHH = 1
    thisHHdat$pFineScalePrevalence = thisHHdat$Z / thisHHdat$N
    thisHHdat$pFineScalePrevalence[thisHHdat$N == 0] = 0 # NaN otherwise
    
    HHdat = c(HHdat, list(thisHHdat))
  }
  
  HHdat
}

# function for sampling clusters from simulated population
sampleClusterSurveys = function(n=NULL, popSim=NULL, HHperClust=25, fixPopPerHH=NULL, 
                                eaSampleStrat=c("pps", "srs"), clustpaList=list(clustpaDHSed), 
                                seed=NULL, verbose=FALSE) {
  
  eaSampleStrat = match.arg(eaSampleStrat)
  
  # set random seed if supplied
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # initialize EA and HH level population info to NULL
  hhPopDat = NULL
  eaPopDat = NULL
  
  # first get ea level data from popSim if need be
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      if(all(popSim[[1]]$nHH == 1)) {
        hhPopDat = popSim
      } else {
        eaPopDat = popSim
      }
    }
  } else {
    stop("popSim has no EA or HH level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  
  # set n if unset
  if(is.null(n)) {
    if(!is.null(hhPopDat)) {
      n = length(hhDatList)
    } else {
      n = length(eaDatList)
    }
  }
  
  if((length(clustpaList) > 1) && (length(clustpaList) != n)) {
    stop(paste0("length mismatch between n (", n, ") and length(clustpaList) (", length(clustpaList), ")"))
  }
  
  # if popDat has length 1 and n > 1, sample multiple surveys from the same population
  if(!is.null(eaPopDat)) {
    resamplePop = ifelse((length(eaPopDat) == 1) && (n > 1), TRUE, FALSE)
  } else if(!is.null(hhPopDat)) {
    resamplePop = ifelse((length(hhPopDat) == 1) && (n > 1), TRUE, FALSE)
  }
  
  if(verbose) {
    print("Simulating surveys")
  }
  surveys = list()
  for(i in 1:n) {
    if(verbose) {
      print(paste0("Sampling survey ", i, "/", n))
    }
    
    # get number of clusters to sample per area
    thisClustpa = clustpaList[[min(c(length(clustpaList), i))]]
    
    # get HH level population information
    if(!is.null(eaPopDat)) {
      # if no HH info available, first get the EA level info, then use to get HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        eaDat = eaPopDat[[1]]
      } else {
        eaDat = eaPopDat[[i]]
      }
      
      hhDat = SUMMER::getHHpop(list(eaDatList = list(eaDat)), fixPopPerHH = fixPopPerHH)[[1]]
      
    } else if(!is.null(hhPopDat)) {
      # we already have HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        hhDat = hhPopDat[[1]]
      } else {
        hhDat = hhPopDat[[i]]
      }
      eaDat = NULL
    }
    
    # obtain info about the EAs
    uniqueEAIs = sort(unique(hhDat$eaIs))
    eaUrbs = hhDat$urban[match(uniqueEAIs, hhDat$eaIs)]
    eaAreas = hhDat$area[match(uniqueEAIs, hhDat$eaIs)]
    
    if(eaSampleStrat == "pps") {
      # calculate number of HHs per EA
      aggHHs = aggregate(hhDat$nHH, by=list(hhDat$eaIs), FUN=sum)
      eaHHs = aggHHs$x
    }
    
    # sample EAs
    sampledEAIs = numeric(0)
    inclusionProbs = numeric(0)
    for(j in 1:nrow(thisClustpa)) {
      
      thisArea = thisClustpa$area[j]
      nUrbEA = thisClustpa$EAUrb[j]
      nRurEA = thisClustpa$EARur[j]
      thisEAIs = uniqueEAIs[eaAreas == thisArea]
      thisEAUrbs = eaUrbs[eaAreas == thisArea]
      thisEAIsUrb = thisEAIs[thisEAUrbs]
      thisEAIsRur = thisEAIs[!thisEAUrbs]
      
      if(eaSampleStrat == "pps") {
        thisEAhhs = eaHHs[eaAreas == thisArea]
        thisEAhhsUrb = thisEAhhs[thisEAUrbs]
        thisEAhhsRur = thisEAhhs[!thisEAUrbs]
      }
      
      # sample EAs and HHs for this area
      sampUrbEAIs = numeric(0)
      sampRurEAIs = numeric(0)
      inclusionProbsUrb = numeric(0)
      inclusionProbsRur = numeric(0)
      if(eaSampleStrat == "srs") {
        if(nUrbEA != 0) {
          sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F)
          inclusionProbsUrb = rep(nUrbEA/length(thisEAIsUrb), nUrbEA)
        }
        if(nRurEA != 0) {
          sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F)
          inclusionProbsRur = rep(nRurEA/length(thisEAIsRur), nRurEA)
        }
      } else if(eaSampleStrat == "pps") {
        require(sampling)
        if(nUrbEA != 0) {
          # sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F, prob=thisEAhhsUrb/sum(thisEAhhsUrb))
          inclusionProbsUrb = nUrbEA * thisEAhhsUrb/sum(thisEAhhsUrb)
          sampUrbEAIs = thisEAIsUrb[as.logical(UPmidzuno(inclusionProbsUrb))]
        }
        if(nRurEA != 0) {
          # sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F, prob=thisEAhhsRur/sum(thisEAhhsRur))
          inclusionProbsRur = nRurEA * thisEAhhsRur/sum(thisEAhhsRur)
          sampRurEAIs = thisEAIsRur[as.logical(UPmidzuno(inclusionProbsRur))]
        }
      } else {
        stop(paste0("eaSampleStrat '", eaSampleStrat, "' not supported"))
      }
      
      # concatenate urban and rural EAs sampled from this area
      thisSampEAIs = c(sampUrbEAIs, sampRurEAIs)
      thisInclusionProbs = c(inclusionProbsUrb, inclusionProbsRur)
      
      # concatenate to vector of all EAs sampled from all areas
      sampledEAIs = c(sampledEAIs, thisSampEAIs)
      inclusionProbs = c(inclusionProbs, thisInclusionProbs)
    }
    
    # subset hhDat to only EAs sampled
    hhSubdat = hhDat[hhDat$eaIs %in% sampledEAIs,]
    
    # sample HHs within chosen EAs
    hhIsTab = aggregate(1:nrow(hhSubdat), by=list(eaIs=hhSubdat$eaIs), function(x) {
      sample(x, HHperClust, replace=FALSE)
    })
    hhIs = sort(c(as.matrix(hhIsTab[,-1])))
    
    hhDatSample = hhSubdat[hhIs,]
    
    # aggregate HH level data for only the EAs sampled
    aggTab = lapply(1:ncol(hhSubdat), function(j) {
      varName = names(hhSubdat)[j]
      aggregate(hhSubdat[,j], by=list(hhSubdat$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhSubdat)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    eaDatSampled = aggTab
    
    # do the same for the actual sampled HHs
    aggTab = lapply(1:ncol(hhDatSample), function(j) {
      varName = names(hhDatSample)[j]
      aggregate(hhDatSample[,j], by=list(hhDatSample$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhDatSample)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    surveyDat = aggTab
    
    # calculate final sampling weight and number HHs in the full EA
    surveyDat$nHHsEA = eaDatSampled$nHH[match(surveyDat$eaIs, eaDatSampled$eaIs)]
    surveyDat$includeProbHH = surveyDat$nHH / surveyDat$nHHsEA
    surveyDat$samplingWeight = surveyDat$N / (surveyDat$includeProbHH * surveyDat$includeProbEA)
    
    
    # concatenate results
    surveys = c(surveys, list(surveyDat))
  }
  
  surveys
}

# stratOrder: if provided, will sort resulting strata to be in the same order as in 
#             stratOrder Otherwise, sorted alphabetically
getClustpaFromSurvey = function(survDat=ed, stratOrder=easpaNGA$area, stratName="area", HHperEA=25) {
  
  if("ns" %in% names(survDat)) {
    survDat$n = survDat$ns
  }
  
  nAreas = length(unique(survDat[[stratName]]))
  
  # first the number of sampled clusters
  nEAtabTmp = aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=length, drop=FALSE)
  nEAtabTmp[is.na(nEAtabTmp[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  nEAtab = data.frame(nEAtabTmp[1:nAreas, 1], EAUrb=nEAtabTmp[(nAreas+1):(2*nAreas), 3], EARur=nEAtabTmp[1:nAreas, 3])
  names(nEAtab)[1] = stratName
  nEAtab$EATotal = nEAtab$EAUrb + nEAtab$EARur
  
  # initialize clustpa
  clustpa = nEAtab
  
  # second the number of households
  clustpa$HHUrb = clustpa$EAUrb * HHperEA
  clustpa$HHRur = clustpa$EARur * HHperEA
  clustpa$HHTotal = clustpa$EATotal * HHperEA
  
  # third the number of people
  popTabTmp = aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=sum, drop=FALSE)
  popTabTmp[is.na(popTabTmp[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  popTab = data.frame(popTabTmp[1:nAreas, 1], popUrb=popTabTmp[(nAreas+1):(2*nAreas), 3], popRur=popTabTmp[1:nAreas, 3])
  names(popTab)[1] = stratName
  popTab$popTotal = popTab$popUrb + popTab$popRur
  
  # concatenate cluster level denominator info to clustpa info
  clustpa$popUrb = popTab$popUrb
  clustpa$popRur = popTab$popRur
  clustpa$popTotal = popTab$popTotal
  
  # sort if need be
  if(!is.null(stratOrder)) {
    ordering = order(stratOrder)
    clustpa = clustpa[ordering,]
  }
  
  clustpa
}

makeIntegrationPointsSim1 = function(regen=FALSE, simStr="") {
  if(! (simStr %in% c("", "_BYM2"))) {
    stop(paste0("invalid simStr: ", simStr))
  }
  
  KDHSurb = 16 # 1 + 3*5: center + 3 inner rings of 5 each
  JInnerUrban = 4
  KDHSrur = 21 # 1 + 3*5 + 1*5: center + 3 inner + 1 outer rings
  JInnerRural = 4
  JOuterRural = 1
  
  out = load(paste0("savedOutput/simStudy1/simPopsSurveys", simStr, ".RData"))
  
  for(i in 1:length(surveysDHS)) {
    if(regen || !file.exists(paste0("savedOutput/simStudy1/intPtsDHS_simStudy1_", i, simStr, ".RData"))) {
      print(paste0("Making integration points survey ", i, "/", length(surveysDHS)))
      thisEdDHS = surveysDHS[[i]]
      
      intPtsDHS = makeAllIntegrationPointsDHS(cbind(thisEdDHS$east, thisEdDHS$north), thisEdDHS$urban, 
                                              areaNames=thisEdDHS$subarea, popPrior=TRUE, 
                                              numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                              JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                              JOuterRural=JOuterRural, adminMap=adm2Full, saveOutput=FALSE)
      
      save(intPtsDHS, file=paste0("savedOutput/simStudy1/intPtsDHS_simStudy1_", i, simStr, ".RData"))
    }
  }
  
  invisible(NULL)
}

makeInputsSim1 = function() {
  
  out = load("savedOutput/simStudy1/simPopsSurveys.RData")
  
  browser()
  
  for(i in 1:length(surveysDHS)) {
    thisEdDHS = surveysDHS[[i]]
    thisEdMICS = surveysMICS[[i]]
    
    out = load(paste0("intPtsDHS_simStudy1_", i, ".RData"))
    
    
    
  }
  
  
}


# simulates data from a BYM2 model instead of SPDE
# BYM2 = spatially structured (ICAR) + unstructured spatial effects at subarea level
# Otherwise matches simPopSPDE() as closely as possible
#
# BYM2 covariance (marginalized form):
#   Cov = (1/tau) * [(1-phi)*I + phi * Q_besag^{-1}_scaled]
# Using eigendecomposition Q_scaled = V diag(gammas) V^T:
#   Cov = V * diag((1/tau) * (1 + phi * gammaTildesm1)) * V^T
# where gammaTildesm1 = 1/gammas - 1 (with 0 for zero eigenvalues)
#
# To simulate: Epsilon_bym2 = V * diag(sqrt(variances)) * z, z ~ N(0,I)

simPopBYM2 = function(nsim=1, easpa, popMat, targetPopMat, poppsub,
                      graphObj,                     # adjacency graph for BYM2 (e.g. admFinalMat)
                      areaCol="area",                # column of popMat specifying which areas match the graph
                      sigmaBYM2=sqrt(0.243),        # BYM2 marginal SD (= 1/sqrt(tau))
                      phi=0.8,                       # mixing: phi*structured + (1-phi)*IID
                      sigmaEpsilon=sqrt(0.463), gamma=0.009,
                      beta0=-3.922, seed=NULL, inla.seed=-1L,
                      nHHSampled=25, stratifyByUrban=TRUE, subareaLevel=TRUE,
                      doFineScaleRisk=FALSE, doSmoothRisk=FALSE,
                      doSmoothRiskLogisticApprox=FALSE,
                      min1PerSubarea=TRUE, offset=NULL, verbose=TRUE,
                      constr=TRUE, scale.model=TRUE,
                      bym2ArgsTMB=NULL) {
  # bym2ArgsTMB: optional list with $Q, $V, $gammaTildesm1 (as returned by
  # prepareBYM2argumentsForTMB). When supplied, we re-use it instead of
  # calling makeQBesag + eigen.spam from scratch, guaranteeing that the
  # simulator's Q, V, and eigenvalues are bit-identical to what TMB consumes.
  # This eliminates a class of cross-platform discrepancies (different LAPACK
  # versions can produce eigendecompositions that differ in sign / ordering /
  # last-bit precision, and QinvSumsNorm computed from V * diag * V^T is
  # particularly sensitive to such drift).

  if (!is.null(seed)) {
    set.seed(seed)
    if (inla.seed < 0) {
      stop("seed specified, but not inla.seed. Set inla.seed to a positive integer to ensure reproducibility")
    }
  }

  pixelCoords = cbind(popMat$east, popMat$north)

  if (is.null(bym2ArgsTMB)) {
    # Build scaled Besag precision matrix and eigendecomposition
    if (verbose) print("Building BYM2 precision matrix and eigendecomposition...")
    Q = makeQBesag(graphObj, constr=constr, scale.model=scale.model, matrixType="spam")
    nGraphAreas = nrow(Q)

    # Eigendecomposition of Q
    eigQ = eigen.spam(Q, symmetric=TRUE)
    gammas = eigQ$values
    V = eigQ$vectors  # nAreas x nAreas eigenvector matrix

    # Compute gammaTildes and gammaTildesm1
    tol = 1e-8
    gammaTildes = 1/gammas
    gammaTildes[abs(gammas) < tol] = 0
    gammaTildesm1 = gammaTildes - 1
  } else {
    if (verbose) print("Re-using supplied bym2ArgsTMB (Q, V, gammaTildesm1)...")
    Q             = bym2ArgsTMB$Q
    V             = bym2ArgsTMB$V
    gammaTildesm1 = bym2ArgsTMB$gammaTildesm1
    gammaTildes   = gammaTildesm1 + 1
    nGraphAreas   = nrow(Q)
    tol           = 1e-8
    # Reconstruct gammas (eigenvalues of Q) from gammaTildes = 1/gammas (with
    # 0 for the null-space mode). Non-finite -> 0 means zero eigenvalue here.
    gammas        = ifelse(abs(gammaTildes) > tol, 1 / gammaTildes, 0)
  }
  
  # Compute marginal variances in eigenbasis for BYM2
  # Var_i = (1/tau) * (1 + phi * gammaTildesm1_i)
  # where tau = 1/sigmaBYM2^2
  tau = 1 / sigmaBYM2^2
  eigenVars = (1/tau) * (1 + phi * gammaTildesm1)
  # For the zero eigenvalue(s), gammaTildesm1 = -1, so eigenVars = (1/tau)*(1-phi)
  eigenVars[eigenVars < 0] = 0  # safety: shouldn't happen for valid phi in [0,1]
  # ENFORCE SUM-TO-ZERO on the BYM2 effect to match the TMB template:
  # TMB samples w_bym2 / u_bym2 on the (n-1)-dim free subspace orthogonal
  # to (1,...,1); see modMDM_BYM2_GH_v2.cpp lines 122-133. The simulator
  # had previously included the constant-mode contribution at variance
  # (1/tau)*(1-phi), which made simulated Epsilon_bym2 lie OUTSIDE the
  # sum-to-zero subspace TMB lives in. The constant got absorbed into alpha
  # and introduced sim-to-sim FE bias of SD ~ sqrt((1-phi)/(tau*n)).
  # Zeroing out the zero-mode eigenvariance below puts the simulator in
  # exactly the same subspace as TMB.
  zeroModeIdx = which(abs(gammas) < tol)
  if (length(zeroModeIdx) > 0) {
    eigenVars[zeroModeIdx] = 0
  }
  eigenSDs = sqrt(eigenVars)
  
  if (verbose) {
    cat(sprintf("  nGraphAreas = %d, tau = %.4f, phi = %.4f\n", nGraphAreas, tau, phi))
    cat(sprintf("  BYM2 marginal SD = %.4f, eigenVar range = [%.4f, %.4f]\n",
                sigmaBYM2, min(eigenVars), max(eigenVars)))
  }
  
  # Map pixels to the areas defined by areaCol
  if (!(areaCol %in% names(popMat))) {
    stop(paste0("Column '", areaCol, "' not found in popMat. Available columns: ",
                paste(names(popMat), collapse=", ")))
  }
  pixelAreas = popMat[[areaCol]]
  # CRITICAL: Epsilon_bym2 = V %*% (...) is in the GRAPH's vertex order (the
  # row/col order of graphObj / Q, preserved by V's eigendecomposition). Pixels
  # MUST be mapped to that same order. Using sort(unique(pixelAreas))
  # (alphabetical) silently SCRAMBLED the field whenever the graph order is not
  # alphabetical (it isn't -- NAME_FINAL is in polygon order): each area then
  # received a DIFFERENT area's effect, destroying the spatial structure
  # (geographic Moran ~0.29 -> ~0.06) and collapsing fitted phi. Map in graph order.
  graphOrder = colnames(graphObj)
  if (is.null(graphOrder)) graphOrder = colnames(as.matrix(Q))
  if (is.null(graphOrder)) graphOrder = rownames(as.matrix(Q))
  if (is.null(graphOrder))
    stop("simPopBYM2: cannot determine graph vertex order (graphObj/Q unnamed); ",
         "field-to-area alignment would be ambiguous")
  nUniqueGraphAreas = length(unique(pixelAreas))
  if (nUniqueGraphAreas != nGraphAreas) {
    stop(paste0("Number of unique values in popMat$", areaCol, " (", nUniqueGraphAreas,
                ") doesn't match graph dimension (", nGraphAreas, ")"))
  }
  if (verbose) print(paste0("  BYM2 at level '", areaCol, "': ", nGraphAreas, " areas (graph order)"))
  pixelAreaIdx = match(pixelAreas, graphOrder)
  if (anyNA(pixelAreaIdx))
    stop("simPopBYM2: some popMat[[areaCol]] values are absent from the graph's ",
         "area names; cannot align the BYM2 field to pixels")
  
  # Simulate BYM2 spatial effects for each simulation
  # For each sim: z ~ N(0,I), then Epsilon_bym2 = V %*% diag(eigenSDs) %*% z
  simVals = matrix(NA, nrow=nrow(pixelCoords), ncol=nsim)
  
  for (i in 1:nsim) {
    z = rnorm(nGraphAreas)
    Epsilon_bym2 = as.numeric(V %*% (eigenSDs * z))  # nGraphAreas-vector
    
    # Expand to pixel level via area/subarea matching
    simVals[, i] = Epsilon_bym2[pixelAreaIdx]
  }
  
  # Add intercept and urban effect
  simVals = simVals + beta0
  simVals = sweep(simVals, 1, gamma * popMat$urban, "+")
  
  # EA-level nugget draws
  totalEAs = sum(easpa$EATotal)
  epsc = matrix(stats::rnorm(totalEAs * nsim, sd=sigmaEpsilon), ncol=nsim)
  
  probsNoNug = expit(simVals)
  logitRiskDraws = simVals
  
  if (!is.null(offset)) {
    logitRiskDraws = sweep(logitRiskDraws, 1, offset, "+")
  }
  
  sigmaEpsilonDraws = rep(sigmaEpsilon, nsim)
  
  if (verbose) print("Using BYM2 model to simulate EA and pixel level populations")
  
  outPixelLevel = simPopCustom(logitRiskDraws=logitRiskDraws,
                               sigmaEpsilonDraws=sigmaEpsilonDraws,
                               easpa=easpa,
                               popMat=popMat,
                               targetPopMat=targetPopMat,
                               stratifyByUrban=stratifyByUrban,
                               doSmoothRisk=doSmoothRisk,
                               doSmoothRiskLogisticApprox=doSmoothRiskLogisticApprox,
                               doFineScaleRisk=doFineScaleRisk,
                               poppsub=poppsub,
                               subareaLevel=subareaLevel,
                               min1PerSubarea=min1PerSubarea,
                               returnEAinfo=TRUE,
                               epsc=epsc,
                               verbose=verbose)
  
  eaPop = list(eaDatList=outPixelLevel$eaDatList, eaSamples=outPixelLevel$eaSamples)
  outPixelLevel$eaDatList = NULL
  outPixelLevel$eaSamples = NULL
  
  if (subareaLevel) {
    if (verbose) print("aggregating from pixel level to subarea level")
    outSubareaLevel = pixelPopToArea(pixelLevelPop=outPixelLevel,
                                     eaSamples=eaPop$eaSamples,
                                     areas=popMat$subarea,
                                     stratifyByUrban=stratifyByUrban,
                                     targetPopMat=targetPopMat,
                                     doFineScaleRisk=doFineScaleRisk,
                                     doSmoothRisk=doSmoothRisk)
    if (verbose) print("aggregating from subarea level to area level")
    tempAreasFrom = popMat$subarea
    tempAreasTo = popMat$area
    areasFrom = sort(unique(tempAreasFrom))
    areasToI = match(areasFrom, tempAreasFrom)
    areasTo = tempAreasTo[areasToI]
    outAreaLevel = areaPopToArea(areaLevelPop=outSubareaLevel,
                                 areasFrom=areasFrom,
                                 areasTo=areasTo,
                                 stratifyByUrban=stratifyByUrban,
                                 doFineScaleRisk=doFineScaleRisk,
                                 doSmoothRisk=doSmoothRisk)
  } else {
    outSubareaLevel = NULL
    if (verbose) print("aggregating from pixel level to area level")
    outAreaLevel = pixelPopToArea(pixelLevelPop=outPixelLevel,
                                  eaSamples=eaPop$eaSamples,
                                  areas=popMat$area,
                                  stratifyByUrban=stratifyByUrban,
                                  doFineScaleRisk=doFineScaleRisk,
                                  doSmoothRisk=doSmoothRisk)
  }
  
  list(eaPop=eaPop, pixelPop=outPixelLevel, subareaPop=outSubareaLevel,
       areaPop=outAreaLevel, logitRiskDraws=logitRiskDraws)
}

# Data simulation for BYM2 model (replaces simData1 with BYM2 spatial effects)
# By default uses admFinalMat (41 MICS stratum level areas) and popMat$stratumMICS.
# Pass graphObj and areaCol to use a different graph/area level.
simData1BYM2 = function(nsim=100, sigmaBYM2=sqrt(0.5), phi=0.8,
                        sigmaEpsilon=sqrt(1.5),
                        beta0=-1.25, gamma=1, betaRest=c(0, 0, 0, .5),
                        easpaDat=easpaNGAed,
                        popMat=popMatNGAThresh,
                        targetPopMat=popMatNGAedThresh,
                        poppsub=poppsubNGAThresh,
                        graphObj=NULL, areaCol="stratumMICS",
                        nHHMICS=16, nHHDHS=25, seed=123,
                        useThreshPopMat=TRUE, fixPopPerHH=NULL,
                        eaSampleStrat="pps", regenPop=FALSE,
                        bym2ArgsTMB=NULL) {
  
  set.seed(seed)
  
  # make sure everything is ordered nicely. CRITICAL: targetPopMat must be
  # permuted with the SAME row order as popMat — see the matching comment in
  # simData1. Sorting only popMat misaligned the smooth-risk truth weights.
  stopifnot(nrow(targetPopMat) == nrow(popMat))
  .perm = order(popMat$subarea)
  popMat = popMat[.perm,]
  targetPopMat = targetPopMat[.perm,]
  stopifnot(all(popMat$subarea == targetPopMat$subarea))
  poppsub = poppsub[order(poppsub$subarea),]

  # construct logit offset vector based on covariates in betaRest
  print("Constructing offset based on covariates...")
  LLcoords = cbind(popMat$lon, popMat$lat)
  tempDesMat = getDesignMat(LLcoords, normalized=TRUE, useThreshPopMat=useThreshPopMat)
  
  load("savedOutput/global/covariatesNorm.RData")
  popVals = extract(pop, LLcoords, method="bilinear")
  
  load("savedOutput/global/popMeanSDCal.RData")
  popMean = ifelse(useThreshPopMat, popMeanCalThresh, popMeanCal)
  popSD = ifelse(useThreshPopMat, popSDCalThresh, popSDCal)
  normPop=(log1p(popVals)-popMeanCal)/popSDCal
  normPop[is.na(normPop)] = min(normPop, na.rm=TRUE)
  
  # get final design matrix
  covRestMat = tempDesMat[,-c(1:3, 7)] # remove int, pop, urb, urbanicity
  covRestMat = cbind(covRestMat, normPop=normPop) # add in normalized population density
  
  # calculate offset
  offset = covRestMat %*% betaRest
  
  # get aggregation info from admin2 areas to MICS strata
  tempAreasFrom = popMat$subarea
  tempAreasTo = popMat$stratumMICS
  areasFrom = sort(unique(tempAreasFrom))
  areasToI = match(areasFrom, tempAreasFrom)
  areasTo = tempAreasTo[areasToI]
  
  # load adjacency graph for BYM2 (default: admFinalMat = 41 MICS strata)
  if (is.null(graphObj)) {
    out = load("savedOutput/global/admFinalMat.RData")
    graphObj = admFinalMat
  }

  # By default, build the BYM2 arguments via the SAME constructor the fitting
  # functions use (prepareBYM2argumentsForTMB with u=0.5, alpha=1/3,
  # constr=TRUE, scale.model=TRUE — see modM_*Sep.R). This makes the
  # simulator's Q / eigenvectors / gammaTildes structurally identical to what
  # TMB consumes, rather than rebuilt independently via makeQBesag +
  # eigen.spam (which can differ in eigenvector sign/order across LAPACK
  # versions and is a known cross-platform reproducibility hazard, especially
  # for QinvSumsNorm). Callers can still override by passing bym2ArgsTMB.
  if (is.null(bym2ArgsTMB)) {
    bym2ArgsTMB = prepareBYM2argumentsForTMB(graphObj, u = 0.5, alpha = 1/3,
                                             constr = TRUE, scale.model = TRUE,
                                             matrixType = "TsparseMatrix")
  }

  # ensure areaCol exists in popMat; add stratumMICS if needed
  if (!(areaCol %in% names(popMat))) {
    if (areaCol == "stratumMICS") {
      print("Adding stratumMICS column to popMat via adm2ToStratumMICS()...")
      popMat$stratumMICS = adm2ToStratumMICS(popMat$subarea)
    } else {
      stop(paste0("Column '", areaCol, "' not found in popMat"))
    }
  }
  
  # Call simPopBYM2 one simulation at a time (matching simData1 structure)
  outFile = "savedOutput/simStudy1/simPopsSurveys_BYM2.RData"

  # Default initial state (fresh run); may be overridden by checkpoint resume.
  surveysDHS  = list()
  surveysMICS = list()
  subareaPops = NULL
  areaPops    = NULL
  stratumPops = NULL
  subareaPops_smoothRisk = NULL
  areaPops_smoothRisk    = NULL
  stratumPops_smoothRisk = NULL
  startI = 1

  # Try to resume from a compatible checkpoint (same seed, <= nsim sims saved).
  if(file.exists(outFile)) {
    chk = new.env()
    try(load(outFile, envir = chk), silent = TRUE)
    canResume <- exists("surveysDHS",      envir = chk) &&
                 exists(".rngStateAtStart", envir = chk) &&
                 exists(".seedUsed",       envir = chk) &&
                 isTRUE(chk$.seedUsed == seed) &&
                 length(chk$surveysDHS) <= nsim
    if(canResume) {
      surveysDHS  <- chk$surveysDHS
      surveysMICS <- chk$surveysMICS
      subareaPops <- chk$subareaPops
      areaPops    <- chk$areaPops
      stratumPops <- chk$stratumPops
      # Smooth-risk fields may not exist in older checkpoints; restore if present.
      if(exists("subareaPops_smoothRisk", envir = chk)) subareaPops_smoothRisk <- chk$subareaPops_smoothRisk
      if(exists("areaPops_smoothRisk",    envir = chk)) areaPops_smoothRisk    <- chk$areaPops_smoothRisk
      if(exists("stratumPops_smoothRisk", envir = chk)) stratumPops_smoothRisk <- chk$stratumPops_smoothRisk
      startI      <- length(surveysDHS) + 1
      if(startI > nsim) {
        print(paste0("All ", nsim, " sims already complete in ", outFile, "; skipping."))
        return(invisible(NULL))
      }
      assign(".Random.seed", chk$.rngStateAtStart, envir = .GlobalEnv)
      print(paste0("Resuming BYM2 sims from ", startI, "/", nsim,
                   " (", startI - 1, " already in ", outFile, ")"))
    } else {
      print("Existing checkpoint is incompatible (different seed or no metadata); starting fresh.")
    }
  }

  print("Simulating populations and surveys...")
  startT = proc.time()[3]
  for(i in startI:nsim) {
    # simulate population at pixel, EA levels 
    print(paste0("Simulating population ", i, "/", nsim))
    
    simPop = simPopBYM2(nsim=1, easpa=easpaDat, popMat=popMat,
                        targetPopMat=targetPopMat, poppsub=poppsub,
                        graphObj=graphObj, areaCol=areaCol,
                        sigmaBYM2=sigmaBYM2, phi=phi,
                        sigmaEpsilon=sigmaEpsilon, gamma=gamma,
                        beta0=beta0, seed=NULL,
                        nHHSampled=nHHDHS, stratifyByUrban=TRUE,
                        subareaLevel=TRUE, doFineScaleRisk=FALSE,
                        doSmoothRisk=TRUE, doSmoothRiskLogisticApprox=FALSE,
                        min1PerSubarea=TRUE, offset=offset, verbose=FALSE,
                        bym2ArgsTMB=bym2ArgsTMB)

    # calculate stratum level population information
    stratPop = SUMMER::areaPopToArea(areaLevelPop=simPop$subareaPop,
                                     areasFrom=areasFrom,
                                     areasTo=areasTo,
                                     stratifyByUrban=TRUE, doFineScaleRisk=FALSE, doSmoothRisk=TRUE)

    # append population information.  We store BOTH the fine-scale prevalence
    # (binomial-realised Z_EA / N_EA, the historical truth used for scoring)
    # AND the smooth (pixel-level expected) risk, since predGrid() produces
    # smooth risk predictions and the consistent comparison is at that level.
    if(is.null(subareaPops)) {
      subareaPops      = simPop$subareaPop$aggregationResults$pFineScalePrevalence
      areaPops         = simPop$areaPop$aggregationResults$pFineScalePrevalence
      stratumPops      = stratPop$aggregationResults$pFineScalePrevalence
      subareaPops_smoothRisk = simPop$subareaPop$aggregationResults$pSmoothRisk
      areaPops_smoothRisk    = simPop$areaPop$aggregationResults$pSmoothRisk
      stratumPops_smoothRisk = stratPop$aggregationResults$pSmoothRisk
    } else {
      # cbind the new pop info to the full set of populations
      subareaPops = cbind(subareaPops,
                          simPop$subareaPop$aggregationResults$pFineScalePrevalence)
      areaPops    = cbind(areaPops,
                          simPop$areaPop$aggregationResults$pFineScalePrevalence)
      stratumPops = cbind(stratumPops,
                          stratPop$aggregationResults$pFineScalePrevalence)
      subareaPops_smoothRisk = cbind(subareaPops_smoothRisk,
                                     simPop$subareaPop$aggregationResults$pSmoothRisk)
      areaPops_smoothRisk    = cbind(areaPops_smoothRisk,
                                     simPop$areaPop$aggregationResults$pSmoothRisk)
      stratumPops_smoothRisk = cbind(stratumPops_smoothRisk,
                                     stratPop$aggregationResults$pSmoothRisk)
    }

    # generate surveys
    print(paste0("Generating surveys for population ", i, "/", nsim))
    # get EA level population information for population i
    thisEApop = simPop$eaPop$eaDatList[1]

    # get associated HH level population information
    thisHHpop = SUMMER::getHHpop(thisEApop, fixPopPerHH=fixPopPerHH)

    # sample DHS survey for this population
    survDHS = SUMMER::sampleClusterSurveys(1, thisHHpop, HHperClust=25, clustpaList=list(clustpaDHSed))
    # Anonymize DHS cluster coordinates by jittering (true -> published) within
    # Admin2; responses (generated at the true location) are unchanged.
    # (sampleTrueWithinPixel for continuous within-pixel true locations is a
    # PARKED long-term faithfulness item -- not applied here; see memory.)
    survDHS[[1]] = addJitterToDHS(survDHS[[1]])

    # now sample the MICS survey. Do some gymnastics to make sure it works for MICS strata
    tempClustpa = clustpaMICSed
    names(tempClustpa)[1] = "area"

    thisHHpop[[1]]$area = adm2ToStratumMICS(thisHHpop[[1]]$subarea)

    survMICS = SUMMER::sampleClusterSurveys(1, thisHHpop, HHperClust=16, clustpaList=list(tempClustpa))

    # concatenate results
    surveysDHS = c(surveysDHS, survDHS)
    surveysMICS = c(surveysMICS, survMICS)

    # Save checkpoint atomically (tmpfile + rename) so a crash mid-save can't
    # leave a corrupt outFile. Stored metadata lets the next call resume.
    .rngStateAtStart <- get(".Random.seed", envir = .GlobalEnv)
    .seedUsed        <- seed
    .nsimRequested   <- nsim
    tmpFile <- paste0(outFile, ".tmp")
    save(subareaPops, areaPops, stratumPops,
         subareaPops_smoothRisk, areaPops_smoothRisk, stratumPops_smoothRisk,
         surveysDHS, surveysMICS,
         .rngStateAtStart, .seedUsed, .nsimRequested,
         file = tmpFile)
    file.rename(tmpFile, outFile)

    # estimate time left and print
    thisT = proc.time()[3]
    timePerIter = (thisT - startT)/(i - startI + 1)
    timeLeft = timePerIter * (nsim - i)

    print(paste0("estimated time remaining: ", (timeLeft/60)/24, " hours"))
  }

  invisible(NULL)
}


