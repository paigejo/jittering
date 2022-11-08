# preprocessing script to produce RData files to load
setwd("~/git/jittering/")

# Load and compress GADM shapefiles ----
library(rgdal)
adm0 = readOGR(dsn = "data/gadm41_NGA_shp/", layer = "gadm41_NGA_0")
adm1 = readOGR(dsn = "data/gadm41_NGA_shp/", layer = "gadm41_NGA_1")
adm2 = readOGR(dsn = "data/gadm41_NGA_shp/", layer = "gadm41_NGA_2")

# compress maps for easy plotting
require(rgeos)
adm0compressed = gSimplify(adm0, tol=.01, topologyPreserve = TRUE)
adm1compressed = gSimplify(adm1, tol=.01, topologyPreserve = TRUE)
adm2compressed = gSimplify(adm2, tol=.01, topologyPreserve = TRUE)
adm0compressed = SpatialPolygonsDataFrame(adm0compressed, adm0@data, match.ID=FALSE)
adm1compressed = SpatialPolygonsDataFrame(adm1compressed, adm1@data, match.ID=FALSE)
adm2compressed = SpatialPolygonsDataFrame(adm2compressed, adm2@data, match.ID=FALSE)

# Construct senatorial districts in Lagos and Kano ----
# we are interested in lagos and kano
head(adm1compressed@data)
adm1compressed@data$NAME_1
lkan = which(adm1compressed@data$NAME_1 %in% c("Lagos", "Kano"))
lkan
plot(adm2compressed, border="lightgray", lwd=.5)
plot(adm1compressed, add=TRUE)
polygon(adm1compressed@polygons[20][[1]]@Polygons[[1]]@coords, border="blue")
polygon(adm1compressed@polygons[25][[1]]@Polygons[[1]]@coords, border="purple")

# senatorial districts in Lagos:
# LaE, LaC LaW (Lagos East, Lagos Central, Lagos West)

# senatorial districts in Kano:
# KaS, KaC, KaN (Kano South, Kano Central, Kano North)

adm2ToSen = read.csv2("data/admin2ToSen.csv")
kanoSens = c("Kano South", "Kano Central", "Kano North")
lagosSens = c("Lagos West", "Lagos Central", "Lagos East")
allSens = c(kanoSens, lagosSens)

## make sure names in the senatorial district list match with gadm ----
# check which admin2 areas are unnaccounted for:
accountedForNamesKano = c()
for(i in 1:length(kanoSens)) {
  thisSen = kanoSens[i]
  
  # get Admin2 areas corresponding to this senatorial district and their indices
  accountedForNamesKano = c(accountedForNamesKano, adm2ToSen$admin2Name_en[adm2ToSen$SenDist_en == thisSen])
}

accountedForNamesLagos = c()
for(i in 1:length(lagosSens)) {
  thisSen = lagosSens[i]
  
  # get Admin2 areas corresponding to this senatorial district and their indices
  accountedForNamesLagos = c(accountedForNamesLagos, adm2ToSen$admin2Name_en[adm2ToSen$SenDist_en == thisSen])
}

senAdm2NamesKano = accountedForNamesKano
senAdm2NamesLagos = accountedForNamesLagos

gadmNamesKano = adm2compressed@data$NAME_2[adm2compressed@data$NAME_1 == "Kano"]
gadmNamesLagos = adm2compressed@data$NAME_2[adm2compressed@data$NAME_1 == "Lagos"]

any(sort(senAdm2NamesKano) != sort(gadmNamesKano))
badNamesKano = which(sort(senAdm2NamesKano) != sort(gadmNamesKano))
cbind(sort(senAdm2NamesKano), sort(gadmNamesKano))[badNamesKano,]

any(sort(senAdm2NamesLagos) != sort(gadmNamesLagos))
badNamesLagos = which(sort(senAdm2NamesLagos) != sort(gadmNamesLagos))
cbind(sort(senAdm2NamesLagos), sort(gadmNamesLagos))[badNamesLagos,]

adm2ToSen$admin2Name_en[adm2ToSen$admin1Name_en == "Kano"] = sort(gadmNamesKano)
adm2ToSen$admin2Name_en[adm2ToSen$admin1Name_en == "Lagos"] = sort(gadmNamesLagos)

## make a new SpatialPolygonsDataFrame from the 6 senatorial districts ----
sen = adm2compressed[-(1:775),]
for(i in 1:length(allSens)) {
  thisSen = allSens[i]
  thisAdm1 = ifelse(i <= 3, "Kano", "Lagos")
  
  # get Admin2 areas corresponding to this senatorial district and their indices
  theseAdmin2s = adm2ToSen$admin2Name_en[adm2ToSen$SenDist_en == thisSen]
  theseAdmin2Is = which((adm2@data$NAME_2 %in% theseAdmin2s) & (adm2@data$NAME_1 == thisAdm1))
  frameToAdd = adm2[theseAdmin2Is,]
  
  if(FALSE) {
    # plot the Admin2 areas one by one for debugging
    for(j in 1:length(theseAdmin2Is)) {
      plot(frameToAdd[j,], lwd=1, border="orange", add=TRUE)
    }
  }
  
  # union all the Admin2 area polygons to get the senatorial district polygon
  frameToAdd$NAME_2 = thisSen
  thisSenPoly = unionSpatialPolygons(frameToAdd, rep(thisSen, length(theseAdmin2Is)))
  
  # construct associated data.frame based on the admin2 data.frame
  tempData = frameToAdd@data[1,]
  row.names(tempData) = thisSen
  
  # append resulting SpatialPolygonsDataFrame row to sen
  thisSenFrame = SpatialPolygonsDataFrame(thisSenPoly, tempData)
  sen = rbind(sen, thisSenFrame)
}

sencompressed = gSimplify(sen, tol=.01, topologyPreserve = TRUE)
sencompressed = SpatialPolygonsDataFrame(sencompressed, sen@data, match.ID=FALSE)

# plot results to test them
plot(adm2compressed, border="lightgray", lwd=.5)
plot(adm1compressed, add=TRUE, lwd=2)
polygon(adm1compressed@polygons[20][[1]]@Polygons[[1]]@coords, border="blue", lwd=1.5)
polygon(adm1compressed@polygons[25][[1]]@Polygons[[1]]@coords, border="purple", lwd=1.5)
plot(sencompressed, lwd=.5, border=rainbow(nrow(sencompressed)), add=TRUE, lwd=1.2)

## replace Lagos and Kano with senatorial districts in final modified 'adm1' map ----

# start by subsetting the sen data attribute so it's compatible with adm1
keepColsSen = which(names(sen@data) %in% names(adm1@data))
tempSen = sen
tempSen@data = tempSen@data[, keepColsSen]
tempSencompressed = sencompressed
tempSencompressed@data = tempSencompressed@data[, keepColsSen]

# do the same for adm1
keepColsAdm1 = which(names(adm1@data) %in% names(tempSen@data))
admFinal = adm1[-lkan, keepColsAdm1]
admFinalcompressed = adm1compressed[-lkan, keepColsAdm1]

# combine senatorial districts and admin1 districts aside from Lagos and Kano
admFinal = rbind(admFinal, tempSen)
admFinal@data$NAME_FINAL = admFinal@data$NAME_1
admFinal@data$NAME_FINAL[(nrow(admFinal@data)-5):nrow(admFinal@data)] = 
  row.names(admFinal@data[(nrow(admFinal@data)-5):nrow(admFinal@data),])
admFinalcompressed = rbind(admFinalcompressed, tempSencompressed)
admFinalcompressed@data$NAME_FINAL = admFinalcompressed@data$NAME_1
admFinalcompressed@data$NAME_FINAL[(nrow(admFinalcompressed@data)-5):nrow(admFinalcompressed@data)] = 
  row.names(admFinalcompressed@data[(nrow(admFinalcompressed@data)-5):nrow(admFinalcompressed@data),])

# plot it to check if everything worked
plot(adm2compressed, border="lightgray", lwd=.5)
plot(adm1compressed, add=TRUE, lwd=2)
polygon(adm1compressed@polygons[20][[1]]@Polygons[[1]]@coords, border="blue", lwd=1.5)
polygon(adm1compressed@polygons[25][[1]]@Polygons[[1]]@coords, border="purple", lwd=1.5)
plot(sencompressed, lwd=.5, border=rainbow(nrow(sencompressed)), add=TRUE, lwd=1.2)

plot(adm2compressed, border="lightgray", lwd=.5)
plot(admFinalcompressed, add=TRUE, lwd=2)
polygon(adm1compressed@polygons[20][[1]]@Polygons[[1]]@coords, border="blue", lwd=1.5)
polygon(adm1compressed@polygons[25][[1]]@Polygons[[1]]@coords, border="purple", lwd=1.5)

# rename so that adm maps are by default compressed (this will prevent accidental R crashes)
adm0Full = adm0
adm1Full = adm1
adm2Full = adm2
admFinalFull = admFinal
senFull = sen
adm0 = adm0compressed
adm1 = adm1compressed
adm2 = adm2compressed
sen = sencompressed
admFinal = admFinalcompressed

## construct domain polygon ----
areas = sapply(adm0Full@polygons[[1]]@Polygons, function(x) {x@area})
maxI = which.max(areas)
adm0Poly = adm0Full@polygons[[1]]@Polygons[[maxI]]@coords

# save results
save(adm0, adm1, adm2, adm0Full, adm1Full, adm2Full, 
     sen, senFull, admFinal, admFinalFull, adm0Poly, 
     file="savedOutput/global/NigeriaMapData.RData")

# Read and clean MICS data ----
# https://www.unicef.org/nigeria/media/1406/file/Nigeria-MICS-2016-17.pdf.pdf
library(dplyr)
library(haven)
dat = read_sav("data/NigeriaMICS5/NigeriaMICS2016-17SPSS/wm.sav")
dim(dat)
# 36176 women elligible for interview, and 34,376 responded (~95% response rate)
# WM7: 1 if they were fully interviewed, 4 if partly interviewed
# WM1, HH1: Cluster number
# WM2: Household number
# WB2: age in years
# WB3: have you ever attended school or preschool (1 is yes, no -> WB7)
# WB4: Highest level of school attended (0: preschool. 1: primary. 2: Secondary. 3: Higher. 0 -> WB7)
# WB5: Highest grade achieved at that level ("00" means first grade at this level not completed)
# WB6: attended secondary school
# WB7: literacy (1: cannot read. 2: read part of sentence. 3: read whole sentence. 4: wrong language. 5: blind)
# NOTE on Nigerian education system: 1 year pre-primary, 6 years primary, 3 years 
#      junior secondary, 3 years senior secondary and 4 years tertiary education. 
#      Hence, new interpretation of WB5: first number is WB4 and second is number 
#      of years in that part of school. This implies 26 or higher means completed 
#      secondary education as evidenced by 4-9, 17-19 and 27-29 not existing. More 
#      info is available at end of MICS report in questionnaires Sec.
# stratum: 10, 20, ..., 180, 191, 192, 193, 200, 210, 220, 230, 241, 242, 243, 250, 260, ..., 370
#          each represent a state, or, in the case of Lagos and Kano, their 3 
#          senatorial districts. It's in alphabetical order by state except FCT is at end.
clustID = dat$WM1
houseID = dat$WM2
ageMICS = dat$WB2
wb3orig = dat$WB3
wb4orig = dat$WB4
wb5orig = dat$WB5
# wb6orig = dat$WB6
literacyOrig = dat$WB7
wb3 = as.numeric(as.character(dat$WB3))
wb4 = as.numeric(as.character(dat$WB4))
wb5 = as.numeric(as.character(dat$WB5))
# wb6 = as.numeric(as.character(dat$WB6))
literacyMICS = as.numeric(as.character(dat$WB7))
fullyInterviewed = dat$WM7

# histogram and summary statistics about literacy of people without formal education
hist(literacyMICS[wb4 == 4], breaks=c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 9.5))
mean(literacyMICS[wb4 == 4] == 1, na.rm=TRUE) # 78% can't read at all
mean(literacyMICS[wb4 == 4] == 2, na.rm=TRUE) # 9% can only read part of sentence
mean(literacyMICS[wb4 == 4] == 3, na.rm=TRUE) # about 1% can read whole sentence
mean(literacyMICS[wb4 == 4] == 4, na.rm=TRUE) # 11% spoke a different language
mean(literacyMICS[wb4 == 4] == 5, na.rm=TRUE) # none are visually impaired
mean(literacyMICS[wb4 == 4] == 9, na.rm=TRUE) # about 1.6% missing but not NA

# compare above summary statistics to those of primary school (above numbers are worse):
mean(literacyMICS[wb4 == 1] == 1, na.rm=TRUE) # about 55% of primary school students can't read at all
mean(literacyMICS[wb4 == 1] == 2, na.rm=TRUE) # about 32% of primary school students can only read part of sentence
mean(literacyMICS[wb4 == 1] == 3, na.rm=TRUE) # about 9% of primary school students can read whole sentence
mean(literacyMICS[wb4 == 1] == 4, na.rm=TRUE) # about 3% of primary school students spoke a different language
mean(literacyMICS[wb4 == 1] == 5, na.rm=TRUE) # about .2% of primary school students are visually impaired
mean(literacyMICS[wb4 == 1] == 9, na.rm=TRUE) # about .2% of primary school students data is missing but not NA

# determine if the women completed their secondary education. Use other information 
# than just wb5 to help in this determination. If education is informal, they 
# have not completed a formal secondary education (and, based on literacy rates, 
# they haven't completed their primary education either...)
edMICS = wb5 >= 26
wb3[is.na(wb3)] = 9
wb4[is.na(wb4)] = 9
edMICS[wb3 == 2] = FALSE
edMICS[(wb4 < 2) | (wb4 == 4)] = FALSE
edMICS[(wb4 > 2) & (wb4 < 4)] = TRUE
mean(is.na(edMICS)) # Nice! We match the missing data rate of 5%

# fill in literacy NAs in a similar way
literacyMICS[(wb4 == 2) | (wb4 == 3)] = TRUE
mean(is.na(literacyMICS)) # Nice! We match the missing data rate of 5%

datMICSwm = data.frame(clustID=clustID, hhID=houseID, age=ageMICS, 
                       secondaryEd=edMICS, literacy=literacyMICS)

# now read in the household data because it contains stratum information
dat = read_sav("data/NigeriaMICS5/NigeriaMICS2016-17SPSS/hh.sav")
stratum = dat$stratum
clustIDhh = dat$HH1
houseIDhh = dat$HH2
urbanhh = dat$HH6 # 1: urban, 2: rural

# determine stratum name from stratum ID. Senatorial districts within 
# Kano and Lagos are in alphabetical order, FCT is at end.
stateNames = sort(unique(admFinal@data$NAME_1))
stratumNames = sort(unique(admFinal@data$NAME_FINAL))
stratumIDs = sort(unique(stratum))
stateStrata = stratumIDs
stateStrata[(stateStrata < 200) & (stateStrata > 190)] = 190
stateStrata[(stateStrata < 250) & (stateStrata > 240)] = 240
stateNames = stateNames[c(1:14, 16:20, 20, 20:25, 25, 25:37, 15)]
stratumNames = stratumNames[c(1:14, 16:41, 15)]
stratMatchTable = data.frame(Area=stateNames, Stratum=stratumNames, StratumID = stratumIDs)
datMICShh = stratMatchTable[match(stratum, stratMatchTable$StratumID),]

# merge the stratum/area information into the women's education data.frame
fullHouseIDhh = paste(as.character(clustIDhh), as.character(houseIDhh), sep=", ")
fullHouseIDwm = paste(as.character(clustID), as.character(houseID), sep=", ")
datMICSwm = cbind(datMICSwm, fullHouseID=fullHouseIDwm)
datMICShh = cbind(datMICShh, fullHouseID=fullHouseIDhh, urban=urbanhh==1)

datMICS = merge(datMICShh, datMICSwm, by="fullHouseID")
datMICS = datMICS[,c(2:ncol(datMICS), 1)]

save(datMICS, file="savedOutput/global/datMICS.RData")

# subset by 20-29 year age
lowAge = 20
highAge = 29
correctAge = (ageMICS >= lowAge) & (ageMICS <= highAge)
correctAge[is.na(correctAge)] = FALSE

datMICS = datMICS[correctAge,]
edMICS = datMICS[!is.na(datMICS$secondaryEd),]
save(edMICS, file="savedOutput/global/edMICS.RData")

# Read DHS data ----
library(haven)
library(fields)
library(zoo)
library(latex2exp)
library(maptools)
library(data.table)

# lines represent births
# data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))
dat <- data.frame(read_dta("data/NG_2018_DHS_09022022_150_113908/NGIR7BDT/NGIR7BFL.DTA"))

## variables from the birth recode
# extract the columns of interest
# b1 - month of birth of child,
# b2 - year of birth of child, 
# b4 - gender of child, # TODO: important?
# b5 - child alive at time of interview, 
# b7 - age of child at death in months completed, 
# v024 - region of residence
# v025 - type of place of residence urban rural (32.3% of pop. is urban in 2009 and 2014, see FR p.30/2)
# v001 - cluster number (in theory there are 1,612 total.  In actuality there are 1,593)
# v002 - household number
# v005 - sample weight out of 1000000. normalized weights sum to number of households
# v006 - month of interview
# V007 - year of interview
# V136 - total number of household members
# V138 - number of women in the household aged 15 to 49
# V137 - number of children in the household aged 5 or under

## variables from this (individual/women) recode
# v024 - region of residence
# v025 - type of place of residence urban rural (32.3% of pop. is urban in 2009 and 2014, see FR p.30/2)
# v001 - cluster number (in theory there are 1,612 total.  In actuality there are 1,593)
# v002 - household number
# v005 - sample weight out of 1000000. normalized weights sum to number of households
# v006 - month of interview
# V007 - year of interview
# v012 - current age in completed years
# v106 - highest education level attended in order of: no educ, primary, secondary, higher
# v107 - highest year comleted (years completed at level given in v106)
# v149 - educational achievement: none (0), incomplete primary (1), completed primary (2), incomplete secondary (3), 
#        complete secondary (4), higher education (5)
subdata <- data.frame(dat[,c('v024', 'v025', 'v001', 'v002', 'v005', 'v012', 'v106', 'v107', 'v149')])

# extract women with the given age range
lowAge = 20
highAge = 29
subdata <- subdata[(subdata[,'v012'] >= lowAge & subdata[,'v012'] <= highAge),]

# add a column for the stratification variable as an interaction between
# the urban/rural indicator 'v025' (1: urban, 2:rural) and the region indicator 'v024'
subdata$regionUral <- with(subdata, interaction(v024, v025), drop=TRUE)

# add a column for the unique households with interaction between
# the household indicator 'v002' and the cluster indicator 'v001'
subdata$hhold <- with(subdata, interaction(v001, v002), drop=TRUE)

# find for each cluster the regionUral indicator
clStrat = subdata[,c("v001", "regionUral", "v024", "v025", "v005")]
clStrat = clStrat[!duplicated(clStrat), ]
colnames(clStrat) = c("clusterID", "regionRural", "region", "urban", "samplingWeight")
clStrat$urban =  clStrat$urban == 1

# determine whether each woman completed their secondary education
completed = subdata$v149 >= 4
completed[is.na(completed)] = TRUE

# get the number of women by cluster
n <- table(subdata[,'v001'])
clusterid <- dimnames(n)[[1]]
n.data = data.frame(clusterID=clusterid, n=as.vector(n))

# get the number of women who completed their secondary education by cluster
# y <- table(subdata[,'v001'])
y = aggregate(completed, list(subdata[,'v001']), sum)
n.data$y = y$x

# add in strata
ed <- merge(n.data, clStrat, by='clusterID', all=TRUE, sort=TRUE)

# Read DHS geographical information
library(rgdal)
spObj = readOGR(dsn = "~/git/jittering/data/NG_2018_DHS_09022022_150_113908/NGGE7BFL/", layer = "KEGE71FL")

# Extract (lon, lat) coordinates of all clusters
geoObj = data.frame(cId = spObj$DHSCLUST, lon = spObj$LONGNUM, lat = spObj$LATNUM)

# Extract coordinates of clusters with data
idx = match(ed$clusterID, geoObj$cId)
ed$lon = geoObj$lon[idx]
ed$lat = geoObj$lat[idx]

# Missing geographical information is assigned value (0,0)
# Remove these
missIdx = which(ed$lon == 0)
ed = ed[-missIdx,]

library(SUMMER)
library(foreign)

gpsDat = readOGR(dsn = "data/NG_2018_DHS_09022022_150_113908/NGGE7BFL/", layer = "NGGE7BFL")
# gpsDat = readShapePoints("Kenya2014gps/KEGE71FL.shp")
coords = attr(gpsDat, "coords")
plot(coords)
world(add=TRUE)
test = coords[,1] < 2 #  these are the observations  whose source is missing.  remove these
sum(test)

#  remove observations with unknown locations and set unknown data to NA
gpsDat=gpsDat[!test,]
names(gpsDat)=c("ID","countryID", "year", "clustID", "countryIDFIPS", "countryIDAdminID", 
                "admin1FIPS","admin1IDSALB", 
                "admin1SALB", "admin1ID", "admin1", "regionID", "region", 
                "source", "urban", "lat", "lon", "altGPS", "altRadar", "coordRef")
gpsDat$altGPS[gpsDat$altGPS == 9999] = NA
gpsDat$altGPS[gpsDat$altRadar == 9999] = NA
gpsDat$urban = gpsDat$urban == "U"
gpsDat$countryIDFIPS[gpsDat$countryIDFIPS == "NULL"] = NA
gpsDat$admin1FIPS[gpsDat$admin1FIPS == "NULL"] = NA
gpsDat$admin1IDSALB[gpsDat$admin1IDSALB == "NULL"] = NA
gpsDat$admin1SALB[gpsDat$admin1SALB == "NULL"] = NA
gpsDat$countryIDAdminID[gpsDat$countryIDAdminID == "NULL"] = NA

save(gpsDat, file="savedOutput/global/gpsDatEd.RData")

# get region and admin data from gps data, add to clusters in ed dataset
gpsI = match(data.frame(rbind(ed$lon, ed$lat)), data.frame(rbind(gpsDat$lon, gpsDat$lat)))
ed$admin1 = gpsDat$admin1[gpsI]
ed$region = gpsDat$region[gpsI]

# get easting and northing using projection
tmp = projNigeria(ed$lon, ed$lat)
ed$east = tmp[,1]
ed$north = tmp[,2]
save(ed, file="savedOutput/global/ed.RData")

# Cleaning urban proportions ----
urbPropsRaw = read.csv2("data/urbProp2017.csv", header=TRUE)
urbProps = urbPropsRaw

# update state names (fill in gaps)
lastStateName = ""
for(i in 1:nrow(urbProps)) {
  thisState = urbProps[i,1]
  
  if(thisState != "") {
    lastStateName = thisState
  } else {
    urbProps[i,1] = lastStateName
  }
}

# to title case
library(stringr)
urbProps$State = str_to_title(urbProps$State)
urbProps$LGA = str_to_title(urbProps$LGA)
urbProps$State[urbProps$State == "Fct Abuja"] = "FCT Abuja"

# extract rows with "Total" in State name. Make separate State level props table
statePropIs = grepl("Total", urbProps$State)
stateProps = urbProps[statePropIs,]
urbProps = urbProps[!statePropIs,]

stateProps = stateProps[-nrow(stateProps),]
stateProps$LGA = NULL

# adjust state prop names to remove "Total"
newStates = gsub(" Total", "", stateProps$State)
stateProps$State = newStates
stateProps$State[stateProps$State == "Fct Abuja"] = "FCT Abuja"

# make sure the state and LGA names match GADM names
gadmStateNames = sort(adm1@data$NAME_1)
propStateNames = sort(stateProps$State)
unmatchedNames = gadmStateNames != propStateNames
cbind(gadm=gadmStateNames, prop=propStateNames)[unmatchedNames,]

gadmNames = sort(adm2@data$NAME_2)
propNames = sort(urbProps$LGA)
unmatchedNames = gadmNames != propNames
cbind(gadm=gadmNames, prop=propNames)[unmatchedNames,]

gadmNames[grepl("Danko", gadmNames)]
urbProps$LGA[urbProps$LGA == "Wasagu/Danko"] = "Danko Wasagu"

propNames = sort(propNames)
unmatchedNames = gadmNames != propNames
cbind(gadm=gadmNames, prop=propNames)

propNames[grepl("Egbado", propNames)]
urbProps$LGA = gsub("Egbado", "Yewa", urbProps$LGA)

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

urbProps$LGA = gsub("  ", " ", urbProps$LGA)

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

urbProps$LGA[urbProps$LGA == "Onuimo"] = "Unuimo"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

temp = adm2@data[adm2@data$NAME_1 == "Adamawa",]
temp[order(temp$NAME_2),]
urbProps$LGA[urbProps$LGA == "Toungo"] = "Teungo"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

adm2@data[adm2@data$NAME_2 == "Mainland",]
urbProps$LGA[urbProps$LGA == "Lagos Mainland"] = "Mainland"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

adm2@data[adm2@data$NAME_2 == "Lake Chad",]
urbProps[urbProps$State == "Borno",]

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames, prop=propNames)

urbProps$LGA[urbProps$LGA == "Kano Municipal"] = "Kano"

cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Aiyekire(Gbonyin)"] = "Gboyin"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Birnin-Magaji"] = "Birnin-Magaji/Kiyaw"
urbProps$LGA[urbProps$LGA == "Birnin Kebbi"] = "Birnin Kudu"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

adm2@data[adm2@data$NAME_2 == "Abuja",]
urbProps[urbProps$State == "FCT Abuja",]
adm2@data[adm2@data$NAME_1 == "Federal Capital Territory",]
urbProps$LGA[urbProps$LGA == "Amac"] = "Abuja"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Munya"] = "Muya"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Dan Musa"] = "Danmusa"
propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Buruku"] = "Bukuru"
propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Bursari"] = "Borsari"
propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Birninwa"] = "Biriniwa"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps[urbProps$LGA == "Ayedade",]
cbind(adm2@data$NAME_2[adm2@data$NAME_1 == "Osun"], 
      urbProps$LGA[urbProps$State == "Osun"])
urbProps$LGA[urbProps$LGA == "Ayedade"] = "Aiyedade"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Akoko Edo"] = "Akoko-Edo"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Umunneochi"] = "Umu Nneochi"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Sagamu"] = "Shagamu"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps$LGA[urbProps$LGA == "Isiukwuato"] = "Isuikwuato"

propNames = sort(urbProps$LGA)
cbind(gadm=gadmNames[!(gadmNames == "Lake Chad")], prop=propNames[-1])

urbProps = urbProps[order(urbProps$LGA),]
urbProps = urbProps[urbProps$LGA != "* Disputed Areas",]

all.equal(sort(stateProps$State), sort(unique(urbProps$State)))

# set all Lagos LGAs to be entirely urban

save(urbProps, stateProps, file="savedOutput/global/urbProps.RData")

# Cleaning population totals ----
popTabRaw = read.csv2("data/popTabNGA_DHS2017.csv", header=TRUE)
popTab = popTabRaw

# make sure state names align with GADM
gadmNames = sort(unique(adm1@data$NAME_1))
tabNames = sort(unique(popTab$state))

sum(gadmNames != tabNames)
which(gadmNames != tabNames)
cbind(GADM=gadmNames, tabNames=tabNames)

popTab = popTab[order(popTab$state),]
popTab$state = gadmNames

# project population forward to 2017 based on growth rates
growth = (1 + popTab$growthRate/100)^11
propUrb = popTab$pctUrban2017/100
popTab$popTotal = round(popTab$popTotal * growth)
popTab$popRural = round(popTab$popTotal * (1 - propUrb))
popTab$popUrban = popTab$popTotal - popTab$popRural

poppa = popTab
names(poppa) = c("area", "growthRateSince2006", "popUrb", "popRur", "popTotal", "pctTotal", "pctUrb")
poppa = poppa[,c(1, 3:7, 2)]
poppaNGA = poppa
save(poppaNGA, file="savedOutput/global/poppaNGA.RData")

# Covariates ----
## Load covariates ----
pop = raster("data/covariates/WorldPopDataNigeria/worldpop/nga_ppp_2020_constrained.tif")
pop = crop(pop, adm0)
pop = mask(pop, adm0)
writeRaster(pop, file="savedOutput/global/pop.tif", format="GTiff")
urb = raster("data/covariates/Urbanization/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0.tif")
urb = projectRaster(urb, crs=CRS(SRS_string="EPSG:4326"))
writeRaster(urb, file="savedOutput/global/urbTemp.tif", format="GTiff", overwrite=TRUE)
urb = crop(urb, adm0)
urb = mask(urb, adm0)
writeRaster(urb, file="savedOutput/global/urb.tif", format="GTiff", overwrite=TRUE)
access = raster("data/covariates/2015_accessibility_to_cities_v1/2015_accessibility_to_cities_v1.0.tif")
access = crop(access, adm0)
access = mask(access, adm0)
writeRaster(access, file="savedOutput/global/access.tif", format="GTiff", overwrite=TRUE)

# alternative try and getting urbanicity, since above is entirely NAs...
urb = raster("data/covariates/Urbanization/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0.tif")
# urb = projectRaster(urb, crs=CRS(SRS_string="EPSG:4326"))
adm0FullProjUrb = spTransform(adm0Full, CRS(urb@crs@projargs))
urb = crop(urb, adm0FullProjUrb)
urb = projectRaster(urb, crs=CRS(SRS_string="EPSG:4326"), method="ngb")
urb = mask(urb, adm0)
writeRaster(urb, file="savedOutput/global/urb.tif", format="GTiff", overwrite=TRUE)

# elevation requires a bit of extra work
elev = readBin(con = "data/covariates/Elevation/g10g" ,what = "integer",n= 129600000, size = 2, signed = TRUE,
               endian = "little")

elev <- matrix(elev, nrow = 6000, ncol = 10800, byrow = TRUE)
elev <- raster(elev)

# Give it lat/lon coords
extent(elev) <- c(0,90,0,50)

# Assign a projection
projection(elev) <- CRS("+proj=longlat +datum=WGS84")

## distance from rivers/lakes
# rivers and lakes centerlines
library(rgdal)
centerlines = readOGR(dsn = "data/covariates/DistanceFromRiverOrLakes/RiversLakesCenterlines/ne_10m_rivers_lake_centerlines/", layer = "ne_10m_rivers_lake_centerlines")

# crop it to bounding box of Nigeria (plus a 2 degree buffer)
ext = extent(adm0@bbox + rbind(c(-2, 2), c(-2, 2)))
centerlines = crop(centerlines, ext)

# project centerlines to easting/northing in meters
centerlinesEN = spTransform(centerlines, sp::CRS("+init=epsg:26391 +units=m"))

ext <- elev@extent
thisCRS <- elev@crs
res <- res(elev)
minDistRiverLakes <- raster(crs =thisCRS, ext = ext, resolution = res)

ext = extent(adm0@bbox + rbind(c(-2, 2), c(-2, 2)))
minDistRiverLakes = crop(minDistRiverLakes, ext)

# projectRaster(urban, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method="ngb")
minDistRiverLakesEN = projectRaster(minDistRiverLakes, crs="+init=epsg:26391 +units=m")

# get the center coordinates of the cells
centers <- xyFromCell(minDistRiverLakesEN, 1:ncell(minDistRiverLakesEN))

# Compute the minimum distances between each cell center and the components of
# the SpatialLinesDataFrame which is called "centerlines"
require(fields)
i=2
minDistRiverLake = c()
require(rgeos)
startTime = proc.time()[3]
for (i in 1:(nrow(centers))) {
  if((i %% 1000) == 0) {
    currentTime = proc.time()[3]
    propDone = i/nrow(centers)
    propLeft = 1-propDone
    timeTaken = (currentTime - startTime) / 60
    totalTimeEst = timeTaken / propDone
    timeLeftEst = totalTimeEst * propLeft
    print(paste0("point ", i, "/", nrow(centers)))
    print(paste0("Time taken (minutes): ", timeTaken))
    print(paste0("Estimated time left (minutes): ", timeLeftEst))
  }
  # point = SpatialPoints(cbind(centers[i,1], centers[i,2]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), bbox = NULL)
  # m <- gDistance(centerlines, point, byid=TRUE)
  
  point = SpatialPoints(cbind(centers[i,1], centers[i,2]), proj4string = CRS("+init=epsg:26391 +units=m"), bbox = NULL)
  m <- gDistance(centerlinesEN, point, byid=TRUE)
  m<-apply(m, 1, min)
  m <- unname(m)
  m[is.nan(m) | is.na(m)] = 0
  
  minDistRiverLake[i] = m
}

# convert from m to km
minDistRiverLake = minDistRiverLake * (1/1000)

# set the values in the raster
test = minDistRiverLakesEN
test = setValues(test, minDistRiverLake)
test = projectRaster(test, crs=CRS(SRS_string="EPSG:4326"))

plot(test)
plot(adm0compressed, add=TRUE)

minDistRiverLakes = test

# save temporarily
save(minDistRiverLakes, file="savedOutput/global/minDistRiverLakes.RData")
load("savedOutput/global/minDistRiverLakes.RData")

# crop and mask rasters to adm0
# pop, urb, access, elev
pop = crop(pop, adm0)
pop = mask(pop, adm0)
plot(pop)

urb = crop(urb, adm0)
urb = mask(urb, adm0)
plot(urb)

access = crop(access, adm0)
access = mask(access, adm0)
plot(access)

elev = crop(elev, adm0)
elev = mask(elev, adm0)
plot(elev)

minDistRiverLakes = crop(minDistRiverLakes, adm0)
minDistRiverLakes = mask(minDistRiverLakes, adm0)
plot(minDistRiverLakes)

## save results ----
writeRaster(pop, file="savedOutput/global/pop.tif", format="GTiff")
writeRaster(urb, file="savedOutput/global/urb.tif", format="GTiff")
writeRaster(access, file="savedOutput/global/access.tif", format="GTiff")
writeRaster(elev, file="savedOutput/global/elev.tif", format="GTiff", overwrite=TRUE)
writeRaster(minDistRiverLakes, file="savedOutput/global/minDistRiverLakes.tif", format="GTiff", overwrite=TRUE)
pop = raster("savedOutput/global/pop.tif")
urb = raster("savedOutput/global/urb.tif")
access = raster("savedOutput/global/access.tif")
elev = raster("savedOutput/global/elev.tif")
minDistRiverLakes = raster("savedOutput/global/minDistRiverLakes.tif")

if(FALSE) {
  # file paths are different on the cluster. Adjust for laptop accordingly:
  pop@file@name = "~/git/jittering/savedOutput/global/pop.tif"
  urb@file@name = "~/git/jittering/savedOutput/global/urb.tif"
  access@file@name = "~/git/jittering/savedOutput/global/access.tif"
  elev@file@name = "~/git/jittering/savedOutput/global/elev.tif"
  minDistRiverLakes@file@name = "~/git/jittering/savedOutput/global/minDistRiverLakes.tif"
}

save(pop, urb, access, elev, minDistRiverLakes, file="savedOutput/global/covariates.RData")

## stratify population ----
# start by getting Admin2/LGA areas of each point in population raster
popMat = rasterToPoints(pop) # NAs are not converted
popCoords = popMat[,1:2]
lgaNamesFull = getRegion2(popCoords, project=FALSE, mapDat=adm2Full, nameVar="NAME_2")
lgaNames = lgaNamesFull$regionNames
mean(lgaNamesFull$multipleRegs) # 0.001441374 x 100% of points are in multiple LGAs due to overlapping LGA geometries

# get urbanicities associated with population density at each point on 0-1 scale
urbVals = extract(urb, SpatialPoints(popCoords, CRS(SRS_string=pop@crs@projargs)))/100

# get urban proportions according to 

# load urban props for each LGA
out = load("savedOutput/global/urbProps.RData")


## transform rasters ----
popVals = getValues(pop)
urbVals = getValues(urb)
accessVals = getValues(access)
elevVals = getValues(elev)
distRiverLakesVals = getValues(minDistRiverLakes)

popValsNorm = log1p(popVals)
popValsNorm = (popValsNorm - mean(popValsNorm, na.rm=TRUE))/sd(popValsNorm, na.rm=TRUE)

urbValsNorm = log1p(urbVals)
urbValsNorm = (urbValsNorm - mean(urbValsNorm, na.rm=TRUE))/sd(urbValsNorm, na.rm=TRUE)

accessValsNorm = log1p(accessVals)
accessValsNorm = (accessValsNorm - mean(accessValsNorm, na.rm=TRUE))/sd(accessValsNorm, na.rm=TRUE)

elevValsNorm = log1p(elevVals)
elevValsNorm = (elevValsNorm - mean(elevValsNorm, na.rm=TRUE))/sd(elevValsNorm, na.rm=TRUE)

distRiverLakesValsNorm = log1p(distRiverLakesVals)
distRiverLakesValsNorm = (distRiversLakesValsNorm - mean(distRiversLakesValsNorm, na.rm=TRUE))/sd(distRiversLakesValsNorm, na.rm=TRUE)

popNorm = pop
values(popNorm) = popValsNorm
urbNorm = urb
values(urbNorm) = urbValsNorm
accessNorm = access
values(accessNorm) = accessValsNorm
elevNorm = elev
values(elevNorm) = elevValsNorm
minDistRiverLakesNorm = minDistRiverLakes
values(minDistRiverLakesNorm) = distRiverLakesValsNorm

writeRaster(popNorm, file="savedOutput/global/popNorm.tif", format="GTiff", overwrite=TRUE)
writeRaster(urbNorm, file="savedOutput/global/urbNorm.tif", format="GTiff", overwrite=TRUE)
writeRaster(accessNorm, file="savedOutput/global/accessNorm.tif", format="GTiff", overwrite=TRUE)
writeRaster(elevNorm, file="savedOutput/global/elevNorm.tif", format="GTiff", overwrite=TRUE)
writeRaster(minDistRiverLakesNorm, file="savedOutput/global/minDistRiverLakesNorm.tif", format="GTiff", overwrite=TRUE)
popNorm = raster("savedOutput/global/popNorm.tif")
urbNorm = raster("savedOutput/global/urbNorm.tif")
accessNorm = raster("savedOutput/global/accessNorm.tif")
elevNorm = raster("savedOutput/global/elevNorm.tif")
minDistRiverLakesNorm = raster("savedOutput/global/minDistRiverLakesNorm.tif")
save(popNorm, urbNorm, accessNorm, elevNorm, minDistRiverLakesNorm, file="savedOutput/global/covariatesNorm.RData")


# Matching urban proportions ----

# first calculate population per LGA
out = load("savedOutput/global/covariates.RData")
adm2FullProj = projNigeriaArea(adm2Full)
adm1FullProj = projNigeriaArea(adm1Full)
lgaArea = getArea(adm2FullProj, "NAME_2", TRUE) / 1000^2
stateArea = getArea(adm1FullProj, "NAME_1", TRUE) / 1000^2
lgaAreaList = list(subarea = names(lgaArea), spatialArea = lgaArea)
stateAreaList = list(area = names(stateArea), spatialArea = stateArea)
poppsubNGA = SUMMER::getPoppsub(kmRes=1, pop=pop, domainPoly=adm0Poly, eastLim=eastLimNGA, northLim=northLimNGA, mapProjection=projNigeria, 
                            poppa=poppaNGA, areapa=stateAreaList, areapsub=lgaAreaList, subareaMapDat=adm2Full, subareaNameVar="NAME_2", 
                            stratifyByUrban=TRUE, areaMapDat=adm1Full, areaNameVar="NAME_1", 
                            areaPolygonSubsetI=NULL, subareaPolygonSubsetI=NULL, 
                            mean.neighbor=50, delta=.1)
SUMMER::poppRegionFromPopMat()

# Use SUMMER's built in urban matching functions

poppsub = SUMMER::setThresholdsByRegion
