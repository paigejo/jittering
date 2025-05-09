# script for women's secondary education in Nigeria application

# load datasets ----
out = load("savedOutput/global/ed.RData")
out = load("savedOutput/global/edMICS.RData")
edMICS = sortByCol(edMICS, "Stratum", admFinal$NAME_FINAL)

# set parameters ----
# Umut settings: 
#   5 urban rings of 15 each (61 points total)
#   10 rural rings of 15 each (136 points total)
# KMICS=100
# KDHSurb = 31 # 4 rings of 10 each
# JInnerUrban = 4
# KDHSrur = 71 # 4 inner + 4 outer rings of 10 each
# JInnerRural = 4
# JOuterRural = 4
KMICS=100
KDHSurb = 11 # 3 rings of 5 each
JInnerUrban = 3
KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
JInnerRural = 3
JOuterRural = 1

if(FALSE) {
  # do some precomputation ----
  
  # make integration points if necessary
  intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
                                            numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE, 
  #                                         numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
  #                                         JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
  #                                         JOuterRural=JOuterRural)
  intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, 
                                          areaNames=ed$subarea, popPrior=TRUE, 
                                          numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                          JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                          JOuterRural=JOuterRural, adminMap=adm2Full)
  
  load("savedOutput/global/intPtsDHS.RData")
  load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  
  AUrbDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasUrban), admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  ARurDHS = makeApointToArea(adm2ToStratumMICS(intPtsDHS$areasRural), admFinal$NAME_FINAL) # 41 x 810
  
  # AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
  # ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICS$Stratum, by=list(strat=edMICS$Stratum, urb=edMICS$urban), FUN=length, drop=FALSE)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  # TODO: EXTEND AMICS TO BE LARGER, INCLUDE DIFFERENT ROW FOR EACH INTEGRATION POINT AND OBSERVATION
  
  # numPerStratUrb = table(edMICS$Stratum[edMICS$urban])
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # startInds = seq(1, KMICS*length(admFinal@data$NAME_FINAL), by=KMICS)
  # getInds = function(intPtI = 1, numPerStrat) {
  #   unlist(mapply(rep, startInds+intPtI-1, each=numPerStrat))
  # }
  # actualIndexUrb = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratUrb))
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  # 41 strata, so rows 1:41 are first int pt, 42:82 second, etc.
  # > head(intPtsMICS$XUrb)
  #        east     north       lon      lat        pop urban      area       subarea     strat
  # 1  547.0015 118.20927  7.351590 5.063059 27186.7228     1      Abia     Aba South      Abia
  # 2 1063.9174 494.24425 12.045797 8.397973  4943.1894     1   Adamawa         Ganye   Adamawa
  # 3  591.5783  93.82219  7.752028 4.840914   158.1402     1 Akwa Ibom     Oruk Anam Akwa Ibom
  # 4  494.2189 197.59492  6.878789 5.782479 12464.0931     1   Anambra        Ihiala   Anambra
  # 5  786.9519 637.69014  9.563900 9.730510  4567.1614     1    Bauchi Tafawa-Balewa    Bauchi
  # 6  395.3649  34.94804  5.983221 4.314700  5334.2306     1   Bayelsa         Brass   Bayelsa
  #   popMatIs int     access       elev distRiversLakes urbanicity  normPop
  # 1      922   1 -2.8759693 -0.5412159       0.5633840  2.6536520 2.559150
  # 2    12438   1 -0.4855544  1.4601349       1.1318496  0.0000000 1.286019
  # 3      574   1 -1.6545579 -0.8598766       1.3266178  0.0000000 1.513281
  # 4     2179   1 -2.5005279 -1.0000544      -0.7473778  0.2977891 1.856981
  # 5    18645   1 -2.3125340  0.7973164      -0.9779331  0.0000000 1.427486
  # 6        1   1  0.9488462  0.0000000      -1.0094755  0.0000000 1.197035
  startStratInds = which(XUrb$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XUrb)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratUrb[x])})) # length nUrb, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nUrb*KMICS, range = 1:41. gives area index for each integration point of each observation
  nUrb = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nUrb) # length nUrb*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIUrb = allAreaIs + (allIntIs-1)*nAreas
  XUrb = XUrb[transformIUrb,] # now XUrb is [K * nObsUrb] x nVar
  
  # cbind(numPerStratUrb, table(test)) # (matches up exactly)
  # cbind(numPerStratUrb, table(transformI)) # (matches up exactly)
  
  # AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  # numPerStratRur = table(edMICS$Stratum[!edMICS$urban])
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # actualIndexRur = c(sapply(1:KMICS, getInds, numPerStrat=numPerStratRur))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  # ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  startStratInds = which(XRur$strat == "Abia") # 1, 42, 83, .... Add 1 to this to get Adamawa inds
  nAreas = nrow(XRur)/KMICS
  areaI = unlist(sapply(1:nAreas, function(x) {rep(x, each=numPerStratRur[x])})) # length nRur, range = 1:41. gives area index for each obs
  allAreaIs = rep(areaI, KMICS) # length nRur*KMICS, range = 1:41. gives area index for each integration point of each observation
  nRur = length(allAreaIs)/KMICS
  allIntIs = rep(1:KMICS, each=nRur) # length nRur*KMICS, range = 1:KMICS. gives int point index for each integration point of each observation
  transformIRur = allAreaIs + (allIntIs-1)*nAreas
  XRur = XRur[transformIRur,]
  
  XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # w matrices are nStrata x K. They should be nObs x K
  wUrban = intPtsMICS$wUrban
  stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  wUrban = wUrban[stratIndexUrbW,]
  
  wRural = intPtsMICS$wRural
  stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICS$Stratum, admFinal$NAME_FINAL)
  # edMICS = edMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  ysUrbMICS = edMICS[edMICS$urban,]$ys
  nsUrbMICS = edMICS[edMICS$urban,]$ns
  ysRurMICS = edMICS[!edMICS$urban,]$ys
  nsRurMICS = edMICS[!edMICS$urban,]$ns
  
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  AUrbMICS = t(AUrbMICS)
  ARurMICS = t(ARurMICS)
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  mode(AUrbMICS) = "numeric"
  mode(ARurMICS) = "numeric"
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  intPtsMICS$XRur = XRur[,-(2:3)]
  intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  intPtsMICS$wUrban = wUrban
  intPtsMICS$wRural = wRural
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  
  # convert A matrices to sparse matrices
  # AUrbMICS = as(AUrbMICS, "sparseMatrix")
  # ARurMICS = as(ARurMICS, "sparseMatrix")
  # AUrbDHS = as(AUrbDHS, "sparseMatrix")
  # ARurDHS = as(ARurDHS, "sparseMatrix")
  AUrbMICS = as.matrix(AUrbMICS)
  ARurMICS = as.matrix(ARurMICS)
  AUrbDHS = as.matrix(AUrbDHS)
  ARurDHS = as.matrix(ARurDHS)
  
  areaidxlocUrbanMICS = apply(AUrbMICS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
  areaidxlocRuralMICS = apply(ARurMICS, 1, function(x) {match(1, x)}) - 1
  areaidxlocUrbanMICS = as.integer(areaidxlocUrbanMICS)
  areaidxlocRuralMICS = as.integer(areaidxlocRuralMICS)
  areaidxlocUrbanDHS = apply(AUrbDHS, 1, function(x) {match(1, x)}) - 1 # TMB indices start from 0
  areaidxlocRuralDHS = apply(ARurDHS, 1, function(x) {match(1, x)}) - 1
  areaidxlocUrbanDHS = as.integer(areaidxlocUrbanDHS)
  areaidxlocRuralDHS = as.integer(areaidxlocRuralDHS)
  
  save(AUrbMICS, ARurMICS, AUrbDHS, ARurDHS, intPtsDHS, intPtsMICS, 
       areaidxlocUrbanMICS, areaidxlocRuralMICS, 
       areaidxlocUrbanDHS, areaidxlocRuralDHS, 
       ysUrbMICS, nsUrbMICS, ysRurMICS, nsRurMICS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
       file="savedOutput/global/edM_DMInputs.RData")
  
  # compile model ----
  # dyn.unload( dynlib("code/modM_MSepsparse"))
  # compile( "code/modM_MSepsparse.cpp", framework="TMBad", safebounds=FALSE)
  
  dyn.unload( dynlib("code/modM_DSep"))
  compile( "code/modM_DSep.cpp", 
           framework="TMBad", safebounds=FALSE)
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/TMB/include" -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppEigen/include"  -DTMB_SAFEBOUNDS -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modM_DSep  -DTMB_LIB_INIT=R_init_modM_DSep  -DTMBAD_FRAMEWORK  -I/usr/local/include   -fPIC  -Wall -g -O2  -c code/modM_DSep.cpp -o code/modM_DSep.o
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o code/modM_DSep.so code/modM_DSep.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
  
  # on Idun:
  # g++ -std=gnu++14 -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/include" -DNDEBUG -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/TMB/include" -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/RcppEigen/include"   -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modM_DSep  -DTMB_LIB_INIT=R_init_modM_DSep  -DTMBAD_FRAMEWORK  -I/cluster/apps/eb/software/OpenSSL/1.1/include -I/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/include -I/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/include -I/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/include -I/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/include -I/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/include -I/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Java/11.0.2/include -I/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/include -I/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/include -I/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include/flexiblas   -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c code/modM_DSep.cpp -o code/modM_DSep.o
  # g++ -std=gnu++14 -shared -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -L/cluster/apps/eb/software/OpenSSL/1.1/lib64 -L/cluster/apps/eb/software/OpenSSL/1.1/lib -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib64 -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib64 -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Java/11.0.2/lib64 -L/cluster/apps/eb/software/Java/11.0.2/lib -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib64 -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib -L/cluster/apps/eb/software/GCCcore/11.3.0/lib64 -L/cluster/apps/eb/software/GCCcore/11.3.0/lib -o code/modM_DSep.so code/modM_DSep.o -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -lR
}

if(FALSE) {
  # test MICS integration points
  image(cbind(c(data_full$wUrbanMICS)[1:734]*100-2, 4*match(edMICS$Stratum[edMICS$urban], sort(unique(edMICS$Stratum)))/41-2, data_full$X_betaUrbanMICS[1:734,-1]), col = hcl.colors(200, "YlOrRd", rev = TRUE))
  image(cbind(4*match(edMICS$Stratum[!edMICS$urban], sort(unique(edMICS$Stratum)))/41-2, data_full$X_betaRuralMICS[1:1449,-1]), col = hcl.colors(200, "YlOrRd", rev = TRUE))
  
  out = load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
  tempXUrb = intPtsMICS$XRur[transformIUrb,]
  tempXRur = intPtsMICS$XRur[transformIRur,]
  head(data_full$X_betaRuralMICS)
  plotWithColor(tempXRur$lon, tempXRur$lat, data_full$X_betaRuralMICS[,2], pch=19, cex=.3)
  plotMapDat(admFinal, new=FALSE)
  plotWithColor(tempXRur$lon, tempXRur$lat, data_full$X_betaRuralMICS[,3], pch=19, cex=.3)
  plotMapDat(admFinal, new=FALSE)
  plotWithColor(tempXRur$lon, tempXRur$lat, data_full$X_betaRuralMICS[,4], pch=19, cex=.3)
  plotMapDat(admFinal, new=FALSE)
  plotWithColor(tempXRur$lon, tempXRur$lat, data_full$X_betaRuralMICS[,5], pch=19, cex=.3)
  plotMapDat(admFinal, new=FALSE)
  plotWithColor(tempXRur$lon, tempXRur$lat, log(c(data_full$wRuralMICS)), pch=19, cex=.3)
  plotMapDat(admFinal, new=FALSE)
  plot(log(c(data_full$wRuralMICS)), data_full$X_betaRuralMICS[,5])
  cor(log(c(data_full$wRuralMICS)), data_full$X_betaRuralMICS[,5])
  
  allNumPerStrat = aggregate(edMICS$Stratum, by=list(strat=edMICS$Stratum, urb=edMICS$urban), FUN=length, drop=FALSE)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  names(numPerStratUrb) = allNumPerStrat[1:41,1]
  names(numPerStratRur) = allNumPerStrat[1:41,1]
  
  # get average covariate value for each stratum
  
  # check weights
  image(data_full$wUrbanMICS)
  image(data_full$wRuralMICS)
  
  out = load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
  nUrb = nrow(XUrb)/KMICS
  nRur = nrow(XRur)/KMICS
  head(XUrb)
  head(intPtsMICS$XUrb)
  XUrb[1:20,]
  
}

# load in TMB function inputs
out = load("savedOutput/global/edM_DMInputs.RData")

# set priors ----
alpha_pri = c(0, 100^2)
beta_pri = c(0, sqrt(1000))

out = load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, 
                                         constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
lambdaTau = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for bym2 precision
lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision

# conditioning by Kriging from Eq (2.30) in Rue Held:
# Ax = e (for A = (0^T 1^T), e = 0), x = (w^T u^T)^T
# x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x - e)
# x* = x - Q_x^-1 A^T (A Q_x^-1 A^T)^-1 (A x)
# x* = x - (Q_x^-1 A^T A x) / sum(Q^+)
# x* = x - (sqrt(phi/tau) Q_{+:}^+ \\ Q_{+:}^+) * sum(u) / sum(Q^+)
# for Q_{+:}^+ = rowSums(Q^+), where * denotes the constrained version of the effect
# Hence, we need Q_{+:}^+ / sum(Q^+):
Qinv = bym2ArgsTMB$V %*% bym2ArgsTMB$Q %*% t(bym2ArgsTMB$V)
QinvSumsNorm = rowSums(Qinv)/sum(Qinv)

# Specify inputs for TMB ----

## specify random effects
rand_effs <- c('alpha', 'beta', 'w_bym2Star', 'u_bym2Star', 
               'nuggetUrbDHS', 'nuggetRurDHS')

# collect input data

data_full = list(
  y_iUrbanDHS=ysUrbDHS, # observed binomial experiment at point i (clust)
  y_iRuralDHS=ysRurDHS, # 
  n_iUrbanDHS=nsUrbDHS, # number binomial trials
  n_iRuralDHS=nsRurDHS, # 
  # AprojUrbanMICS=AUrbMICS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  # AprojRuralMICS=ARurMICS, # 
  areaidxlocUrban=areaidxlocUrbanDHS, # [nIntegrationPointsUrban * nObsUrban] length vector of areal indices associated with each observation
  areaidxlocRural=areaidxlocRuralDHS, # [nIntegrationPointsRural * nObsRural] length vector of areal indices associated with each observation
  X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralDHS=intPtsDHS$covsRur, # 
  wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
  wRuralDHS=intPtsDHS$wRural, # 
  
  Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
  # V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
  alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
  beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
  tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
  gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
  QinvSumsNorm=QinvSumsNorm, 
  lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
  lambdaTau=lambdaTau, # determines PC prior for tau
  lambdaTauEps=lambdaTauEps, 
  options=0 # 1 for adreport of log tau and logit phi
)

anyna = function(x) {any(is.na(x))}
myDim = function(x) {
  dims = dim(x)
  if(is.null(dims)) {
    length(x)
  }
  else {
    dims
  }
}

if(FALSE) {
  # for testing purposes
  any(sapply(data_full, anyna))
  sapply(data_full, anyna)
  sapply(data_full, myDim)
  sapply(data_full, class)
  hist(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS, breaks=50)
  hist(data_full$y_iRuralDHS/data_full$n_iRuralDHS, breaks=50)
  hist(data_full$n_iUrbanDHS-data_full$y_iUrbanDHS, breaks=seq(-0.5, 17.5, by=1))
  hist(data_full$n_iRuralDHS-data_full$y_iRuralDHS, breaks=seq(-0.5, 24.5, by=1))
  hist(data_full$y_iRuralDHS, breaks=seq(-0.5, 16.5, by=1))
  mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1)
  mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0)
  mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1) + mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 0)
  mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0) + mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 1)
  
  hist(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, breaks=50)
  hist(data_full$y_iRuralMICS/data_full$n_iRuralMICS, breaks=50)
  hist(data_full$n_iUrbanMICS-data_full$y_iUrbanMICS, breaks=seq(-0.5, 17.5, by=1))
  hist(data_full$y_iRuralMICS, breaks=seq(-0.5, 16.5, by=1))
  mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1)
  mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0)
  mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1) + mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 0)
  mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0) + mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 1)
}

# initial parameters
initUrbP = sum(c(data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanDHS))
initRurP = sum(c(data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

tmb_params <- list(log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   alpha = initAlpha, # intercept
                   beta = c(initBeta1, rep(0, ncol(data_full$X_betaUrbanDHS)-1)), 
                   w_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                   u_bym2Star = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                   nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                   nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
)

# make TMB fun and grad ----
# dyn.load( dynlib("code/modM_DSepsparse"))
dyn.load( dynlib("code/modM_DSep"))
TMB::config(tmbad.sparse_hessian_compress = 1)
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modM_DSep')
# objFull <- MakeADFun(data=data_full,
#                      parameters=tmb_params,
#                      hessian=TRUE,
#                      DLL='modM_DSep')

lower = rep(-10, length(obj[['par']]))
upper = rep( 10, length(obj[['par']]))

# make wrapper functions that print out parameters and function values
funWrapper = function(par, badParVal=1e10) {
  if(any(par < lower) || any(par > upper)) {
    return(badParVal)
  }
  print(par)
  objVal = testObj[['fn']](par)
  parNames = names(par)
  parVals = par
  parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
  parStr = paste(parStrs, collapse=", ")
  print(paste0("objective: ", objVal, " for parameters, ", parStr))
  objVal
}

grWrapper = function(par, badParVal=1e10) {
  if(any(par < lower) || any(par > upper)) {
    return(rep(badParVal, length(par)))
  }
  print(par)
  grVal = testObj[['gr']](par)
  parNames = names(par)
  parVals = par
  parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
  parStr = paste(parStrs, collapse=", ")
  grStr = paste(grVal, collapse=", ")
  print(paste0("gradient: ", grStr, " for parameters, ", parStr))
  grVal
}

if(FALSE) {
  # testing before running optimization
  initAlphaBeta = rep(1, 6)
  initParFull = objFull$par
  initParFull[1:6] = initAlphaBeta
  objFull$fn(initParFull)
  testRep = obj$report(initParFull)
  
  system.time(test <- obj$fn(obj$par))
  # 363.236  for non-sparse
  # 355.211 for sparse
  
  system.time(test <- obj$gr(obj$par))
  # 54.5 for non-sparse
  # 55.3 for sparse
  
  testRep = obj$report()
  
  range(testRep$latentFieldUrbMICS)
  range(testRep$latentFieldRurMICS)
  range(testRep$latentFieldUrbDHS)
  range(testRep$latentFieldRurDHS)
  
  range(testRep$fe_iUrbanMICS)
  range(testRep$fe_iRuralMICS)
  range(testRep$fe_iUrbanDHS)
  range(testRep$fe_iRuralDHS)
  
  range(testRep$projepsilon_iUrbanMICS)
  range(testRep$projepsilon_iRuralMICS)
  range(testRep$projepsilon_iUrbanDHS)
  range(testRep$projepsilon_iRuralDHS)
  
  range(testRep$liksUrbMICS)
  range(testRep$liksRurMICS)
  range(testRep$liksUrbDHS)
  range(testRep$liksRurDHS)
  
  range(testRep$logDet)
  sapply(testRep, range)
  
  any(apply(data_full$wUrbanDHS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wRuralDHS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wUrbanMICS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wRuralMICS, 1, function(x) {all(x == 0)}))
}

# * Run TMB ----
# set start of optimization
# obj$par = optParINLA
# obj$par = optPar

{
  # tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
  tolSeq = 1e-06
  testObj = obj
  optPar = testObj$par
  startTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = obj
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
    testObj$env$tracepar = TRUE
    print(paste0("optimizing for tol = ", thisTol, "."))
    opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                  method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
    optPar = opt1$par
    if(!is.null(opt1$message)) {
      print(paste0("error for tol = ", thisTol, ". Message:"))
      print(opt1$message)
      next
    }
    else {
      print(paste0("completed optimization for tol = ", thisTol, ""))
      
      ## Get standard errors
      print("getting standard errors...")
      sdTime = system.time(
        SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                             bias.correct = TRUE,
                             bias.correct.control = list(sd = TRUE))
      )[3]
      # SD0
      print(sdTime/60)
      # 3.9447 minutes for intern=FALSE
      
      if(SD0$pdHess) {
        print("Optimization and PD hess calculation done!")
        break
      }
      else {
        print("Hessan not PD. Rerunning optimization with stricter tol...")
        
      }
    }
    
    
    
  }
  endTime = proc.time()[3]
  sdTime/60
  totalTime = endTime - startTime
  print(paste0("optimization took ", totalTime/60, " minutes"))
  # optimization took 2.36478333333313 minutes (for intern=FALSE)
}

if(FALSE) {
  badPar = c(SD0$par.fixed, SD0$par.random)
  
  tempGrad = obj$env$f(badPar,order=1)
  tempHess = obj$env$spHess(badPar,random=TRUE)
  range(tempHess)
  
  tempEig = eigen(tempHess)
  range(tempEig$values)
}

# opt0 <- nlminb(start       =    obj[['par']],
#                objective   =    obj[['fn']],
#                gradient    =    obj[['gr']],
#                lower = rep(-10, length(obj[['par']])),
#                upper = rep( 10, length(obj[['par']])),
#                control     =    list(trace=1))

# * TMB Posterior Sampling ----

if(FALSE) {
  funWrapper(optPar)
  testRep = obj$report(obj$env$last.par)
  testPar = obj$env$last.par
  testPar[6:9]= 1
  testRep = obj$report(testPar)
  feUrb = testRep$fe_iUrbanMICS
  feRur = testRep$fe_iRuralMICS
  latentUrb = testRep$latentFieldUrbMICS
  latentRur = testRep$latentFieldRurMICS
  
  out = load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
  intPtsMICS = straightenMICS(intPtsMICS)
  tempXUrb = intPtsMICS$XUrb[transformIUrb,]
  tempXRur = intPtsMICS$XRur[transformIRur,]
  
  pdf("figures/test/edM_Dlatent.pdf", width=6, height=6)
  zlim = range(c(latentUrb, latentRur))
  plotWithColor(tempXRur$lon, tempXRur$lat, latentRur, pch=19, cex=.3, zlim = zlim)
  plotWithColor(tempXUrb$lon, tempXUrb$lat, latentUrb, pch=19, cex=.3, add=TRUE)
  plotMapDat(admFinal, new=FALSE)
  dev.off()
  
  pdf("figures/test/edM_DpredUrban.pdf", width=6, height=6)
  plotWithColor(tempXRur$lon, tempXRur$lat, tempXRur$urban, pch=19, cex=.3, zlim=c(0,1))
  plotWithColor(tempXUrb$lon, tempXUrb$lat, tempXUrb$urban, pch=19, cex=.3, add=TRUE, zlim=c(0,1))
  plotMapDat(admFinal, new=FALSE)
  dev.off()
  
  pdf("figures/test/edM_DempiricalProp.pdf", width=6, height=6)
  plotWithColor(jitter(tempXRur$lon, amount=.08), jitter(tempXRur$lat, amount=.08), rep(data_full$y_iRuralMICS/data_full$n_iRuralMICS, KMICS), pch=19, cex=c(data_full$wRuralMICS), zlim=c(0,1))
  plotWithColor(jitter(tempXUrb$lon, amount=.08), jitter(tempXUrb$lat, amount=.08), rep(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, KMICS), pch=19, cex=c(data_full$wUrbanMICS), add=TRUE, zlim=c(0,1))
  plotMapDat(admFinal, new=FALSE)
  dev.off()
  
  pdf("figures/test/edM_DtrueUrban.pdf", width=6, height=6)
  quilt.plot(popMatNGAThresh$lon, popMatNGAThresh$lat, popMatNGAThresh$urban, nx=259, ny=212)
  plotMapDat(admFinal, new=FALSE)
  dev.off()
}



## summary(SD0, 'report')
## summary(SD0, 'fixed')

save(SD0, obj, totalTime, sdTime, file="savedOutput/ed/fitM_DSep.RData")
out = load("savedOutput/ed/fitM_DSep.RData")

gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=1000, admLevel="stratMICS", 
                     quantiles=c(0.025, 0.1, 0.9, 0.975), sep=TRUE)
# new results:
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\ 
# \hline
# (Int) & -1.74 & -1.95 & -1.88 & -1.60 & -1.54 \\ 
# urb & 0.48 & 0.26 & 0.34 & 0.64 & 0.71 \\ 
# access & -0.14 & -0.24 & -0.20 & -0.08 & -0.05 \\ 
# elev & 0.15 & 0.02 & 0.06 & 0.23 & 0.28 \\ 
# distRiversLakes & 0.09 & -0.03 & 0.01 & 0.17 & 0.21 \\ 
# popValsNorm & 0.81 & 0.60 & 0.68 & 0.95 & 1.01 \\ 
# sigmaSq & 0.61 & 0.36 & 0.43 & 0.81 & 0.99 \\ 
# phi & 0.84 & 0.45 & 0.66 & 0.97 & 0.99 \\ 
# sigmaEpsSq & 1.34 & 1.15 & 1.21 & 1.47 & 1.54 \\ 
# \hline
# \end{tabular}
# \end{table}

# old results:
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\ 
# \hline
# (Int) & -1.76 & -1.97 & -1.90 & -1.63 & -1.56 \\ 
# urb & 0.32 & 0.08 & 0.17 & 0.48 & 0.56 \\ 
# access & -0.19 & -0.30 & -0.26 & -0.12 & -0.09 \\ 
# elev & 0.14 & 0.01 & 0.06 & 0.23 & 0.27 \\ 
# distRiversLakes & 0.09 & -0.04 & 0.01 & 0.17 & 0.21 \\ 
# popValsNorm & 0.82 & 0.62 & 0.69 & 0.95 & 1.02 \\ 
# sigmaSq & 0.63 & 0.37 & 0.44 & 0.84 & 1.00 \\ 
# phi & 0.85 & 0.48 & 0.65 & 0.97 & 0.99 \\ 
# sigmaEpsSq & 1.42 & 1.23 & 1.29 & 1.56 & 1.62 \\  
# \hline
# \end{tabular}
# \end{table}
save(gridPreds, file="savedOutput/ed/gridPredsM_DSep.RData")
out = load("savedOutput/ed/gridPredsM_DSep.RData")

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPreds, file="savedOutput/ed/stratPredsM_DSep.RData")
save(admin1Preds, file="savedOutput/ed/admin1PredsM_DSep.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsM_DSep.RData")
out = load("savedOutput/ed/stratPredsM_DSep.RData")
out = load("savedOutput/ed/admin1PredsM_DSep.RData")
out = load("savedOutput/ed/admin2PredsM_DSep.RData")

summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, 
               gridPreds=gridPreds)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrr}
# \hline
# & Est & Q0.025 & Q0.975 \\ 
# \hline
# X.Int. & -1.80 & -1.96 & -1.91 & -1.68 & -1.61 \\ 
# beta & 0.94 & 0.77 & 0.83 & 1.05 & 1.11 \\ 
# beta.1 & -0.04 & -0.13 & -0.10 & 0.01 & 0.05 \\ 
# beta.2 & 0.16 & 0.01 & 0.06 & 0.26 & 0.31 \\ 
# beta.3 & -0.01 & -0.17 & -0.11 & 0.10 & 0.16 \\ 
# beta.4 & 0.83 & 0.63 & 0.71 & 0.96 & 1.03 \\ 
# sigmaSq & 1.36 & 1.00 & 1.10 & 1.63 & 1.81 \\ 
# phi & 0.90 & 0.69 & 0.80 & 0.98 & 0.99 \\ 
# sigmaEpsSq & 0.50 & 0.42 & 0.45 & 0.56 & 0.59 \\ 
# \hline
# \end{tabular}
# \end{table}
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=NULL, 
          plotNameRoot="edFusionM_DSep")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edFusionM_DSep", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edFusionM_DSep", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edFusionM_DSep", plotNameRootAreal="Admin2")


