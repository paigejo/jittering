# some tests
out = load("savedOutput/validation/edVal.RData")
out = load("savedOutput/validation/edMICSval.RData")

nUrban = length(ysUrban)
is = 1:nUrban
par(mfrow=c(2,2))
plot(coordsUrban[is,], pch=".")
plot(coordsUrban[is+nUrban,], pch=".")
plot(coordsUrban[is+nUrban*2,], pch=".")
# plot(coordsUrban[is+nUrban*3,], pch=".")

is = apply(AUrban!= 0, 1, which)
coords = mesh.s$loc[is,1:2]
plot(coords, pch=".", col="blue")

nRural = length(ysRural)
is = 1:nRural
par(mfrow=c(2,2))
plot(coordsRural[is,], pch=".")
plot(coordsRural[is+nRural,], pch=".")
plot(coordsRural[is+nRural*2,], pch=".")
# plot(coordsRural[is+nRural*3,], pch=".")

is = apply(ARural!= 0, 1, which)
coords = mesh.s$loc[is,1:2]
plot(coords, pch=".", col="blue")


wUrbanTest = matrix(c(rep(1, nrow(wUrban)), rep(0, nrow(wUrban) * (ncol(wUrban) - 1))), ncol=ncol(wUrban))
wRuralTest = matrix(c(rep(1, nrow(wRural)), rep(0, nrow(wRural) * (ncol(wRural) - 1))), ncol=ncol(wRural))
image.plot(wRural)
image.plot(wRuralTest)

data_full <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
                  num_iRural = length(ysRural),  # Total number of rural observations
                  num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                  y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
                  y_iRural   = ysRural, # num. of pos. rural obs in the cluster
                  n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
                  n_iRural   = nsRural,  # num. of rural exposures in the cluster
                  n_integrationPointsUrban = n_integrationPointsUrban, 
                  n_integrationPointsRural = n_integrationPointsRural, 
                  wUrban = wUrbanTest, 
                  wRural = wRuralTest, 
                  X_alphaUrban  = matrix(1, nrow = length(ysUrban), ncol = 1),# des.mat for urban observations
                  X_alphaRural  = matrix(1, nrow = length(ysRural), ncol = 1),# des.mat for rural observations
                  M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                  M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                  M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                  AprojUrban = AUrban,             # Projection matrix (urban)
                  AprojRural = ARural,             # Projection matrix (rural)
                  options = c(1, ## if 1, use normalization trick
                              1), ## if 1, run adreport
                  # normalization flag.
                  flag = 1,
                  alpha_pri = alpha.pri, ## normal vector of mean, sd for intercept
                  matern_pri = matern.pri ## 
)

tmb_params <- list(alpha = 0.0, # intercept
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   log_kappa = 0, # SPDE parameter related to the range
                   Epsilon_s = rep(0, mesh.s[['n']]) # RE on mesh vertices
)

## make a list of things that are random effects
rand_effs <- c('Epsilon_s')

## make the autodiff generated likelihood func & gradient
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modSPDEJitter')

## we can normalize the GMRF outside of the nested optimization,
## avoiding unnecessary and expensive cholesky operations.
obj <- normalize(obj, flag="flag", value = 0)

obj[['par']]
objSimple[['par']]
obj[['fn']](obj[['par']])
objSimple[['fn']](objSimple[['par']])

# urbanInds = 1 + c(outer(which(dat$urban)-1, (0:(length(ysUrban)-1))*length(ysUrban), "+"))
A.projUrban = A.proj[dat$urban,]
testA.projUrban = AUrban[1:length(ysUrban),]
image.plot(as.matrix(A.projUrban))
image.plot(as.matrix(testA.projUrban))
all.equal(as.matrix(A.projUrban), as.matrix(testA.projUrban))

# Check which DHS points have misclassified urbanicity under our 
# classification (assuming DHS coordinates are exact)
LLcoords = cbind(ed$lon, ed$lat)
Xmat = getDesignMat(LLcoords, normalized=TRUE)
head(Xmat)
classUrban = Xmat[,3]
pdf("figures/test/urbanMisclassPts.pdf", width=6, height=6)
plot(LLcoords[,1], LLcoords[,2], col=c("black", "red")[(ed$urban != classUrban) + 1], pch=19, cex=.3, asp=1)
plotMapDat(adm1, new=FALSE)
dev.off()

predCols = makeBlueGreenYellowSequentialColors(64)
quantCols = makePurpleYellowSequentialColors(64)

# calculate percentage of points in each adm1 area with misclassified urbanicity
out = aggregate(ed$urban != classUrban, by=list(ed$area), FUN=mean, drop=FALSE)
out = out[match(out$Group.1, adm1$NAME_1),]
pdf("figures/test/urbanMisclassAdm1.pdf", width=6, height=6)
plotMapDat(adm1, out$x, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, asp=1, cols=predCols)
dev.off()

# calculate percentage of points in each adm1 area with misclassified urbanicity
out = aggregate(ed$urban != classUrban, by=list(factor(ed$subarea, levels=adm2$NAME_2)), FUN=mean, drop=FALSE)
out = out[match(out$Group.1, adm2$NAME_2),]
pdf("figures/test/urbanMisclassAdm2.pdf", width=6, height=6)
plotMapDat(adm2, out$x, varAreas=adm2$NAME_2, regionNames=adm2$NAME_2, asp=1, border=rgb(.6, .6, .6), lwd=.6, col=rev(quantCols))
plotMapDat(adm1, lwd=1, new=FALSE)
dev.off()



urbCols = makeBlueSequentialColors(64)
urbColsDiverging = ThreeColorDivergingScale(64, valRange=range(poppaNGA$pctUrb), center=median(poppaNGA$pctUrb), 
                                            col1="green", col2=rgb(.95,.95,.95), col3="blue")
pdf("figures/test/pctUrbanAdm1.pdf", width=6, height=6)
plotMapDat(adm1, poppaNGA$pctUrb, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, asp=1, cols=urbColsDiverging)
dev.off()

out = aggregate(edMICS$ns, by=list(edMICS$Area), FUN=sum, drop=FALSE)
out = out[match(out$Group.1, adm1$NAME_1),]
pdf("figures/test/nMICSadm1.pdf", width=6, height=6)
plotMapDat(adm1, out$x, varAreas=adm1$NAME_1, regionNames=adm1$NAME_1, asp=1, cols=predCols)
dev.off()



# check out of sample Xurb for MICS dataset, fold 1
out = load("~/git/jittering/savedOutput/validation/datM_DM.RData")
out
dat = datMDM[[11]]
thisEdMICS = dat$edMICSOutOfSample
thisEdMICSUrb = thisEdMICS[thisEdMICS$urban,]
Xurb = dat$dataOutOfSample$X_betaUrbanMICS
nUrb = nrow(thisEdMICSUrb)
ys = dat$edMICSOutOfSample$ys[dat$edMICSOutOfSample$urban]
ns = dat$edMICSOutOfSample$ns[dat$edMICSOutOfSample$urban]

head(Xurb)
res=100
out=load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
head(intPtsMICS$XUrb)
head(intPtsMICS$XUrb[1:41,c(3,4,7,9,10,12:14,16)])
trueXurb = intPtsMICS$XUrb[,c(3,4,7,9,10,12:14,16)]

allNumPerStrat = aggregate(thisEdMICS$Stratum, by=list(strat=thisEdMICS$Stratum, urb=thisEdMICS$urban), FUN=length, drop=FALSE)
numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
numPerStratUrb[is.na(numPerStratUrb)] = 0
numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
numPerStratRur[is.na(numPerStratRur)] = 0

thisI = unlist(sapply(1:41, function(i) {rep(i, each=numPerStratUrb[i])}))
head(trueXurb[thisI,], 10)
head(Xurb, 10)

head(trueXurb[thisI+41,], 10)
Xurb[(nUrb+1):(nUrb+10),]


# Now check weights
wUrbTrue = intPtsMICS$wUrban
wUrb = dat$dataOutOfSample$wUrbanMICS

wUrb[1:5,1:5]
wUrbTrue[thisI,][1:5,1:5]
all.equal(wUrb, wUrbTrue[thisI,])

Aurb = dat$dataOutOfSample$AprojUrbanMICS
stratumI = apply(Aurb, 1, function(x) {which(x == 1)})
stratum = admFinal$NAME_FINAL[stratumI]

stratumTrue = trueXurb[thisI,]$strat
all.equal(stratum, stratumTrue)

# check ys, ns
sum(edMICSval$ys[(edMICSval$fold == 1) & (edMICSval$urban)])/sum(edMICSval$ns[(edMICSval$fold == 1) & (edMICSval$urban)])
sum(ys)/sum(ns)

# Now check in sample covariates

thisEdMICS = dat$edMICSInSample
thisEdMICSUrb = thisEdMICS[thisEdMICS$urban,]
Xurb = dat$MakeADFunInputs$data$X_betaUrbanMICS
nUrb = nrow(thisEdMICSUrb)
ys = thisEdMICS$ys[thisEdMICS$urban]
ns = thisEdMICS$ns[thisEdMICS$urban]

head(Xurb)
res=100
out=load(paste0("savedOutput/global/intPtsMICS_", res, "_adm2Cov.RData"))
head(intPtsMICS$XUrb)
head(intPtsMICS$XUrb[1:41,c(3,4,7,9,10,12:14,16)])
trueXurb = intPtsMICS$XUrb[,c(3,4,7,9,10,12:14,16)]

allNumPerStrat = aggregate(thisEdMICS$Stratum, by=list(strat=thisEdMICS$Stratum, urb=thisEdMICS$urban), FUN=length, drop=FALSE)
numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
numPerStratUrb[is.na(numPerStratUrb)] = 0
numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
numPerStratRur[is.na(numPerStratRur)] = 0

thisI = unlist(sapply(1:41, function(i) {rep(i, each=numPerStratUrb[i])}))
head(trueXurb[thisI,], 10)
head(Xurb, 10)

head(trueXurb[thisI+41,], 10)
Xurb[(nUrb+1):(nUrb+10),]

all.equal(as.matrix(Xurb[1:nUrb,2:5]), as.matrix(trueXurb[thisI,c(6:9)]))
all.equal(as.matrix(Xurb[1:nUrb + nUrb,2:5]), as.matrix(trueXurb[thisI+41,c(6:9)]))
all.equal(as.matrix(Xurb[1:nUrb + nUrb*(res-1),2:5]), as.matrix(trueXurb[thisI+41*(res-1),c(6:9)]))

# check weights

wUrb = dat$MakeADFunInputs$data$wUrbanMICS

wUrb[1:5,1:5]
wUrbTrue[thisI,][1:5,1:5]
all.equal(wUrb, wUrbTrue[thisI,])

Aurb = dat$MakeADFunInputs$data$AprojUrbanMICS
stratumI = apply(Aurb, 1, function(x) {which(x == 1)})
stratum = admFinal$NAME_FINAL[stratumI]

stratumTrue = trueXurb[thisI,]$strat
all.equal(stratum, stratumTrue)

# check ys, ns
sum(edMICSval$ys[edMICSval$fold == 1])/sum(edMICSval$ns[edMICSval$fold == 1])

# check all covariates

trueXurb = intPtsMICS$XUrb[,c(3,4,7,9,10,12:14,16)]
head(trueXurb)
predScale = makeBlueGreenYellowSequentialColors(64)
plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$access, pch=19, cex=.3, colScale=predScale)
plotMapDat(admFinal, new=FALSE)
plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$access, pch=19, cex=.3, colScale=predScale, add=TRUE)

plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$elev, pch=19, cex=.3, colScale=predScale)
plotMapDat(admFinal, new=FALSE)
plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$elev, pch=19, cex=.3, colScale=predScale, add=TRUE)

plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$normPop, pch=19, cex=.3, colScale=predScale)
plotMapDat(admFinal, new=FALSE)
plotWithColor(trueXurb$lon, trueXurb$lat, trueXurb$normPop, pch=19, cex=.3, colScale=predScale, add=TRUE)

edPixelI = numeric(nrow(edVal))
for(i in 1:nrow(edVal)) {
  print(paste0("observation ", i, "/", nrow(edVal)))
  thisAdm2 = edVal[i,]$subarea
  popMatSubset = which(popMatNGAThresh$subarea == thisAdm2)
  thisPopMat = popMatNGAThresh[popMatSubset,]
  dists = rdist(cbind(edVal[i,]$east, edVal[i,]$north), cbind(thisPopMat$east, thisPopMat$north))
  edPixelI[i] = popMatSubset[which.min(dists)]
}
pixelUrb = popMatNGAThresh[edPixelI,]$urban
equalClassification = pixelUrb == edVal$urban
mean(equalClassification) # 0.8462654
mean(equalClassification[edVal$urban]) # 0.8260105
mean(equalClassification[!edVal$urban]) # 0.8604938
mean(equalClassification[pixelUrb]) # 0.806175
mean(equalClassification[!pixelUrb]) # 0.8756281


plotWithColor(edVal$lon, edVal$lat, equalClassification, colScale=c("red", "green"), pch=19, cex=.4, 
              xlab="longitude", ylab="latitude")
plotMapDat(admFinal, new=FALSE)
plotWithColor(edVal$lon, edVal$lat, equalClassification, colScale=c("red", "green"), pch=19, cex=.4, add=TRUE)

plotWithColor(edVal$lon[edVal$urban], edVal$lat[edVal$urban], equalClassification[edVal$urban], colScale=c("red", "green"), pch=19, cex=.4, 
              xlab="longitude", ylab="latitude")
plotMapDat(admFinal, new=FALSE)
plotWithColor(edVal$lon[edVal$urban], edVal$lat[edVal$urban], equalClassification[edVal$urban], colScale=c("red", "green"), pch=19, cex=.4, add=TRUE)

# load covariates at prediction locations
LLcoords = cbind(popMatNGAThresh[edPixelI,]$lon, popMatNGAThresh[edPixelI,]$lat)
Xmat = getDesignMat(LLcoords, normalized = TRUE)

out = load("savedOutput/global/popMeanSDCal.RData")
popMean = popMeanCalThresh
popSD = popSDCalThresh
popValsNorm = (log1p(Xmat[,2]) - popMean) * (1/popSD)

plot(popValsNorm, equalClassification)
plot(popValsNorm, edVal$urban)
plot(Xmat[,4], equalClassification)
plot(Xmat[,4], edVal$urban)
plot(Xmat[,5], equalClassification)
plot(Xmat[,5], edVal$urban)
plot(Xmat[,6], equalClassification)
plot(Xmat[,6], edVal$urban)
tempTab = data.frame(pop=popValsNorm, access=Xmat[,4], elev=Xmat[,5], distRiversLakes=Xmat[,6], 
                     equalClass=equalClassification, pixelUrb=pixelUrb)
pairs(tempTab, pch=19, cex=.6, col=rgb(0,0,0,.1))
glm(edVal$urban ~ popValsNorm, family=binomial)

out = load("savedOutput/global/intPtsMICS_100_adm2Cov.RData")
test = straightenMICS(intPtsMICS)
test$strataMICS
intPtsMICS$strataMICS
rowSums(intPtsMICS$wRural)
rowSums(test$wRural)

# check urban/rural data averages both weighted and unweighted
# [1] "mean predicted urban prob: 0.654948768863194"
# [1] "mean predicted rural prob: 0.308878898550139"
sum(ed$y[ed$urban]/ed$n[ed$urban] * ed$samplingWeight[ed$urban]) / sum(ed$samplingWeight[ed$urban])
# 0.6872929
sum(ed$y[!ed$urban]/ed$n[!ed$urban] * ed$samplingWeight[!ed$urban]) / sum(ed$samplingWeight[!ed$urban])
# 0.3051247
sum(edMICS$ys[edMICS$urban]/edMICS$ns[edMICS$urban] * edMICS$samplingWeight[edMICS$urban]) / sum(edMICS$samplingWeight[edMICS$urban])
# 0.659839
sum(edMICS$ys[!edMICS$urban]/edMICS$ns[!edMICS$urban] * edMICS$samplingWeight[!edMICS$urban]) / sum(edMICS$samplingWeight[!edMICS$urban])
# 0.2446335

sum(ed$y[ed$urban])/sum(ed$n[ed$urban])
# 0.6480836
sum(ed$y[!ed$urban])/sum(ed$n[!ed$urban])
# 0.2788789
sum(edMICS$ys[edMICS$urban])/sum(edMICS$ns[edMICS$urban])
# 0.700167
sum(edMICS$ys[!edMICS$urban])/sum(edMICS$ns[!edMICS$urban])
# 0.285413

tempUrb = 0.5*(0.6480836 - 0.2788789 + 0.700167 - 0.285413)
# 0.3919794
tempInt = mean(c(0.6480836, 0.2788789, 0.700167, 0.285413)) - 0.5*tempUrb
# 0.2821459
urbVal = tempInt + tempUrb
# 0.6741253


# check validation predictive distribution ----
zlimsStrat = NULL
zlimsAdmin1 = NULL
zlimsAdmin2 = NULL
zlimsGrid = NULL
gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=1000, admLevel="stratMICS", 
                     quantiles=c(0.025, 0.1, 0.9, 0.975), sep=TRUE)

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)

plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edFusionM_DMSepValTest2", plotNameRootAreal="Strat", CIwidthLims=zlimsStrat)
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edFusionM_DMSepValTest2", plotNameRootAreal="Admin1", CIwidthLims=zlimsAdmin1)
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edFusionM_DMSepValTest2", plotNameRootAreal="Admin2", CIwidthLims=zlimsAdmin2)
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=NULL, 
          plotNameRoot="edFusionM_DMSepValTest2", CIwidthLims=zlimsGrid)
