# edSimpleINLA

# script for women's secondary education in Nigeria application

# load datasets ----
out = load("savedOutput/global/ed.RData")
out = load("savedOutput/global/edMICS.RData")

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
KMICS=300
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
  
  # AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  # ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = makeApointToArea(ed$subarea[ed$urban], adm2$NAME_2) # 775 x 569 nArea x nObsUrb
  ARurDHS = makeApointToArea(ed$subarea[!ed$urban], adm2$NAME_2) # 775 x 810
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICS$Stratum, by=list(strat=edMICS$Stratum, urb=edMICS$urban), FUN=length, drop=FALSE)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  # XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  # AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  # TODO: EXTEND AMICS TO BE LARGER, INCLUDE DIFFERENT ROW FOR EACH INTEGRATION POINT AND OBSERVATION
  
  # numPerStratUrb = table(edMICS$Stratum[edMICS$urban])
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  # AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  # XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  # ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  # numPerStratRur = table(edMICS$Stratum[!edMICS$urban])
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  # ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  # XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # w matrices are nStrata x K. They should be nObs x K
  # wUrban = intPtsMICS$wUrban
  # stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  # wUrban = wUrban[stratIndexUrbW,]
  # 
  # wRural = intPtsMICS$wRural
  # stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  # wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICS$Stratum, admFinal$NAME_FINAL)
  # edMICS = edMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  # ysUrbMICS = edMICS[edMICS$urban,]$ys
  # nsUrbMICS = edMICS[edMICS$urban,]$ns
  # ysRurMICS = edMICS[!edMICS$urban,]$ys
  # nsRurMICS = edMICS[!edMICS$urban,]$ns
  
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  # AUrbMICS = t(AUrbMICS)
  # ARurMICS = t(ARurMICS)
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  # mode(AUrbMICS) = "numeric"
  # mode(ARurMICS) = "numeric"
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  # intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  # intPtsMICS$XRur = XRur[,-(2:3)]
  # intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  # intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  # intPtsMICS$wUrban = wUrban
  # intPtsMICS$wRural = wRural
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[1:nrow(AUrbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[1:nrow(ARurDHS),-1]
  
  # convert A matrices to sparse matrices
  # AUrbMICS = as(AUrbMICS, "sparseMatrix")
  # ARurMICS = as(ARurMICS, "sparseMatrix")
  AUrbDHS = as(AUrbDHS, "sparseMatrix")
  ARurDHS = as(ARurDHS, "sparseMatrix")
  # AUrbMICS = as.matrix(AUrbMICS)
  # ARurMICS = as.matrix(ARurMICS)
  # AUrbDHS = as.matrix(AUrbDHS)
  # ARurDHS = as.matrix(ARurDHS)
  
  save(AUrbDHS, ARurDHS, intPtsDHS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
       file="savedOutput/global/edSimpleInputs.RData")
  
  # compile model ----
  # dyn.unload( dynlib("code/modBYM2JitterDHS2sparse"))
  # compile( "code/modBYM2JitterDHS2sparse.cpp", framework="TMBad", safebounds=FALSE)
  
  dyn.unload( dynlib("code/modBYM2simple"))
  compile( "code/modBYM2simple.cpp", 
           framework="TMBad", safebounds=FALSE)
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/TMB/include" -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppEigen/include"  -DTMB_SAFEBOUNDS -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/usr/local/include   -fPIC  -Wall -g -O2  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
  
  # on Idun:
  # g++ -std=gnu++14 -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/include" -DNDEBUG -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/TMB/include" -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/RcppEigen/include"   -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/cluster/apps/eb/software/OpenSSL/1.1/include -I/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/include -I/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/include -I/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/include -I/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/include -I/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/include -I/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Java/11.0.2/include -I/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/include -I/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/include -I/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include/flexiblas   -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # g++ -std=gnu++14 -shared -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -L/cluster/apps/eb/software/OpenSSL/1.1/lib64 -L/cluster/apps/eb/software/OpenSSL/1.1/lib -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib64 -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib64 -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Java/11.0.2/lib64 -L/cluster/apps/eb/software/Java/11.0.2/lib -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib64 -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib -L/cluster/apps/eb/software/GCCcore/11.3.0/lib64 -L/cluster/apps/eb/software/GCCcore/11.3.0/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -lR
}

# load in TMB function inputs
out = load("savedOutput/global/edSimpleIntInputs.RData")

urbI = ed$urban
ordEd = rbind(ed[urbI,], 
              ed[!urbI,])
covMat = rbind(intPtsDHS$covsUrb, 
               intPtsDHS$covsRur)

system.time(mod <- runBYM2simple(dat=ed, strategy="eb", fast=FALSE, areaOrder=adm2$NAME_2))
summary(mod)
optParINLA = mod$summary.hyperpar[,1] # tau, phi, tauEps
optParINLA = c(log(optParINLA[1]), logit(optParINLA[2]), log(optParINLA[3]))



predCols = makeBlueGreenYellowSequentialColors(64)
quantCols = makePurpleYellowSequentialColors(64)

plotMapDat(adm2, expit(mod$summary.random$idx[1:length(adm2),2] + mod$summary.fixed[1,1]), varAreas=adm2$NAME_2, 
           regionNames=adm2$NAME_2, border=rgb(.6, .6, .6), cols=predCols, lwd=.6)
plotMapDat(adm1, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, new=FALSE)

plotMapDat(adm2, mod$summary.random$idx[1:length(adm2),2] + mod$summary.fixed[1,1], varAreas=adm2$NAME_2, 
           regionNames=adm2$NAME_2, border=rgb(.4, .4, .4))
plotMapDat(adm1, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, new=FALSE)

CI80 = 
  expit(mod$summary.random$idx[1:length(adm2),6] + mod$summary.fixed[1,2]) - 
  expit(mod$summary.random$idx[1:length(adm2),4] + mod$summary.fixed[1,2])
plotMapDat(adm2, CI80, varAreas=adm2$NAME_2, asp=1, 
           regionNames=adm2$NAME_2, border="grey")
plotMapDat(adm1, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, asp=1, new=FALSE)




gridPreds = predGridINLA(mod, popMat=popMatNGAThresh, nsim=1000, admLevel="adm2", 
                         quantiles=c(0.025, 0.1, 0.9, 0.975))
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\ 
# \hline
# (Int) & -0.48 & -0.59 & -0.54 & -0.41 & -0.38 \\ 
# sigmaSq & 1.97 & 1.56 & 1.69 & 2.24 & 2.40 \\ 
# phi & 0.94 & 0.88 & 0.88 & 0.98 & 0.99 \\ 
# sigmaEpsSq & 1.40 & 1.22 & 1.26 & 1.57 & 1.63 \\ 
# \hline
# \end{tabular}
# \end{table}

save(gridPreds, file="savedOutput/ed/gridPredsSimpleIntINLA.RData")
out = load("savedOutput/ed/gridPredsSimpleIntINLA.RData")

# stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
# admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
# save(stratPreds, file="savedOutput/ed/stratPredsSimpleIntINLA.RData")
# save(admin1Preds, file="savedOutput/ed/admin1PredsSimpleIntINLA.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsSimpleIntINLA.RData")
# out = load("savedOutput/ed/stratPredsSimpleIntINLA.RData")
# out = load("savedOutput/ed/admin1PredsSimpleIntINLA.RData")
out = load("savedOutput/ed/admin2PredsSimpleIntINLA.RData")

# plotPreds(SD0, obj, popMat=popMatNGAThresh, 
#           gridPreds=gridPreds, arealPreds=NULL, 
#           plotNameRoot="edSimpleIntINLA")
# plotPreds(SD0, obj, popMat=popMatNGAThresh, 
#           gridPreds=gridPreds, arealPreds=stratPreds, 
#           plotNameRoot="edSimpleIntINLA", plotNameRootAreal="Strat")
# plotPreds(SD0, obj, popMat=popMatNGAThresh, 
#           gridPreds=gridPreds, arealPreds=admin1Preds, 
#           plotNameRoot="edSimpleIntINLA", plotNameRootAreal="Admin1")
plotPreds(SD0=NULL, tmbObj=NULL, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edSimpleIntINLA", plotNameRootAreal="Admin2")



