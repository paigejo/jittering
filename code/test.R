# some tests


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