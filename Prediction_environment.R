## Creating a dataset for prediction with values every degree
## Created: 16 / 4 / 15
## Last edited: 3 / 6 / 15
## Isabel Fenton
##
## Based on the code from Analysis_MARGO.R
##
## Previous file: 1311 LDGPaper/Reanalysis/Analysis_MARGO.R
## Next file: 1311 LDGPaper/Reanalysis/Analysis_MARGO.R

## Inputs ------------------------------------------------------------------
# requires ldg.margo.mod

## Outputs -----------------------------------------------------------------
# 150421_OceanBoundaries.RData - points for defining the ocean boundaries

## Libraries ---------------------------------------------------------------
library(fields) # map
data(world.dat) # map
library(sp) # point.in.polygon
library("rhdf5") # for processing h5 files
library(colorRamps) # for matlab.like colours
source("../../../Code/maps.R") # for maps
setwd("C:/Documents/Science/PhD/Work/1311 LDGPaper/Reanalysis/")


## 1. Generate dataframe for predicting ------------------------------------
# meanSST.1deg, sdSST.1deg, mean.mld.t, depth10deg, logProd.mn.ann, meanSal.0m, sdSal.0m, prop2.oxy, Ocean2, delta_carb_ion
ldg.p.margo <- data.frame(Longitude = rep(-179.5:179.5, 180), Latitude = rep(-89.5:89.5, each = 360))
names(rsr.margo.mod)


## 2. Ocean ---------------------------------------------------------------
with(rsr.margo.mod, distrib.map(Longitude, Latitude, Ocean2, key = "FALSE"))

## 2i. Create polygons of the ocean boundaries -----------------------------
# Each ocean is defined using the locator() function. Save these results
# save(pacific, atlantic, atlantic.2, indian, pacific.2, file = "Outputs/150421_OceanBoundaries.RData")
load("Outputs/150421_OceanBoundaries.RData")
# with these boundaries, then define polygons of the oceans. Exclude the Mediterranean and high northern latitudes where I have no data
pacific.1 <- pacific
pacific.1$x <- c(pacific$x[length(pacific$x)], -180, -180, pacific$x[2:length(pacific$x)])
pacific.1$y <- c(pacific$y[length(pacific$y)], -90, pacific$y[2], pacific$y[2:length(pacific$y)])
points(pacific.1$x, pacific.1$y, type = "l", col = "blue")

# without the gulf of mexico
atlantic.1 <- atlantic
atlantic.1$x <- c(atlantic$x, pacific$x[9:length(pacific$x)], rev(atlantic.2$x), atlantic$x[1])
atlantic.1$y <- c(atlantic$y, pacific$y[9:length(pacific$x)], rev(atlantic.2$y), atlantic$y[1])
points(atlantic.1$x, atlantic.1$y, type = "l", col = "yellow")

# if I want to include the gulf of mexico
atlantic.3 <- atlantic
atlantic.3$x <- c(atlantic$x[1:9], pacific$x[4:length(pacific$x)], rev(atlantic.2$x), atlantic$x[1])
atlantic.3$y <- c(atlantic$y[1:9], pacific$y[4:length(pacific$x)], rev(atlantic.2$y), atlantic$y[1])
points(atlantic.3$x, atlantic.3$y, type = "l", col = "black")

indian.1 <- indian
indian.1$x <- c(indian$x[length(indian$x)], rev(atlantic.2$x)[1:5], indian$x)
indian.1$y <- c(indian$y[length(indian$y)],rev(atlantic.2$y)[1:5], indian$y)
points(indian.1$x, indian.1$y, type = "l", col = "purple")

pacific.3 <- pacific.2
pacific.3$x <- c(180, rev(indian$x)[1:(length(indian$x) - 7)], rev(pacific.2$x), 180)
pacific.3$y <- c(-90, rev(indian$y)[1:(length(indian$y) - 7)], rev(pacific.2$y), -90)
points(pacific.3$x, pacific.3$y, type = "l", col = "green")

save(pacific.1, pacific.3, atlantic.1, indian.1, file = "Outputs/150509_OceanPolygons.RData")
save(pacific.1, pacific.3, atlantic.3, indian.1, file = "Outputs/150603_Gulf_OceanPolygons.RData")
rm(pacific, atlantic, atlantic.2, indian, pacific.2, atlantic.1)

## 2ii. Add column to ldg.p.margo ------------------------------------------
# calculate the ocean for each point
ldg.p.margo$Ocean2 <- NA
ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, pacific.1$x, pacific.1$y) == 1)] <- "Pacific"
ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, pacific.3$x, pacific.3$y) == 1)] <- "Pacific"
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Pacific", ], points(Longitude, Latitude, col = "yellow"))

ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, atlantic.3$x, atlantic.3$y) == 1)] <- "Atlantic"
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Atlantic", ], points(Longitude, Latitude, col = "red"))

ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, indian.1$x, indian.1$y) == 1)] <- "Indian"
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Indian", ], points(Longitude, Latitude, col = "orange"))

# set the sites on land back to NA
# n.b. necessary to add points for the whole of antarctica.
tmp <- which(world.dat$x == 180)[2]
land <- world.dat
land$x <- c(world.dat$x[1:tmp], 180, -180, -180, world.dat$x[(tmp + 1):length(world.dat$x)])
land$y <- c(world.dat$y[1:tmp], -90, -90, world.dat$y[which(world.dat$x == -180)], world.dat$y[(tmp + 1):length(world.dat$x)])
ldg.p.margo$Ocean2 <- as.factor(ldg.p.margo$Ocean2)

points(land$x, land$y, type = "l", col = "green")

ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, land$x, land$y) == 1)] <- NA

## 2iii. Check and tidy up -------------------------------------------------
with(rsr.margo.mod, distrib.map(Longitude, Latitude, Ocean2, key = "FALSE"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Atlantic", ], points(Longitude, Latitude, col = "red"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Indian", ], points(Longitude, Latitude, col = "orange"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Pacific", ], points(Longitude, Latitude, col = "yellow"))

rm(tmp, atlantic.3, indian.1, pacific.1, pacific.3)


## 3. SST.4km -------------------------------------------------------------
setwd("../../../Project/BFD/Environmental/SST_4km/")

## 3i. Extract data for lat / long -----------------------------------------
# extract the lat/long values (this is the same for all months)
lon <- h5read("month01_combined.h5", "Longitude")
lat <- h5read("month01_combined.h5", "Latitude")

# calculate the closest coordinates for each site
long.margo <- sapply(ldg.p.margo$Longitude, match.4km, lon)
lat.margo <- sapply(ldg.p.margo$Latitude, match.4km, lat)

# set up a dataframe to hold the values
SST.4km.p <- ldg.p.margo

# for each month
for (i in 1:12) {
  # read in the file
  if (i < 10) {
    file.nm <- paste("month0", i, "_combined.h5", sep = "")
  } else {
    file.nm <- paste("month", i, "_combined.h5", sep = "")
  }
  sst <- as.vector(h5read(file.nm, "Clim_SST_Filled"))
  # extract the data for the margo sites
  SST.4km.p <- cbind(SST.4km.p, sst[long.margo + length(lon) * (lat.margo - 1)])
  colnames(SST.4km.p)[ncol(SST.4km.p)] <- paste("sst.", i, sep = "")
}
rm(i, file.nm, sst, lat, lon, lat.margo, long.margo)
# save this out so that I can come back to it if necessary
save(SST.4km.p, file = "150423_SST.4km.p_working.RData")

## convert to SST
# 0 is missing data, 1 is land. Set both of these to NA
for (i in 1:12) {
  SST.4km.p[which(SST.4km.p[, paste("sst.", i, sep = "")] == 0), paste("sst.", i, sep = "")] <- NA
  SST.4km.p[which(SST.4km.p[, paste("sst.", i, sep = "")] == 1), paste("sst.", i, sep = "")] <- NA
}
rm(i)

# check this has worked
head(SST.4km.p)

# convert to actual values
# scale * pixel_Value + offset
# scale = 0.075 and offset = - 3.0
SST.4km.p[, grep("sst", names(SST.4km.p))] <- SST.4km.p[, grep("sst", names(SST.4km.p))] * 0.075 - 3

## 3ii. Calculate mean and SD ----------------------------------------------
SST.4km.p$meanSST <- rowMeans(SST.4km.p[, grep("sst", names(SST.4km.p))], na.rm = T)
SST.4km.p$sdSST <- apply(SST.4km.p[, grep("sst", names(SST.4km.p))], 1, sd, na.rm = T)
# check these
head(SST.4km.p)

## 3iii. Images ------------------------------------------------------------
# create maps of monthly sst
for (i in 1:12) 
{
  png(paste("SST_4km_", i, ".png", sep = ""), 800, 500)
  distrib.map(SST.4km.p$Longitude, SST.4km.p$Latitude, SST.4km.p[, paste("sst.", i, sep = "")], palette = "matlab.like", col.land = "black", col.water = "white", pch = 15, cex = 0.4, max.col = 36)
  dev.off()
}
rm(i)

# map of mean SST
png(file = "meanSST_4km.png", width = 800, height = 500)
with(SST.4km.p, distrib.map(Longitude, Latitude, meanSST, palette = "matlab.like", col.land = "black", col.water = "white", pch = 15, cex = 0.4))
dev.off()

# map of sd SST
png(file = "sdSST_4km.png", width = 800, height = 500)
with(SST.4km.p, distrib.map(Longitude, Latitude, sdSST, palette = "matlab.like", col.land = "black", col.water = "white", pch = 15, cex = 0.4))
dev.off()

## 3iv. Add data to ldg.p.margo --------------------------------------------
ldg.p.margo$meanSST.4km <- SST.4km.p$meanSST
ldg.p.margo$sdSST.4km <- SST.4km.p$sdSST

save(SST.4km.p, file = "150423_SST_4km_p.RData")
rm(SST.4km.p)


## 4. SST.1deg ---------------------------------------------------------
# load in the data
setwd("../SST_1deg_SD")
load("150414_Temperature_1deg.RData")

# check that the lat / long match in length
sum(mean.t.depth$Long[order(mean.t.depth$Lat, mean.t.depth$Long)] != ldg.p.margo$Long)
sum(mean.t.depth$Lat[order(mean.t.depth$Lat, mean.t.depth$Long)] != ldg.p.margo$Lat)

# add a column for mean SST
ldg.p.margo$meanSST.1deg <- mean.t.depth$depth0m[order(mean.t.depth$Lat, mean.t.depth$Long)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, meanSST.1deg))

# add a column for sd SST
ldg.p.margo$sdSST.1deg <- sd.t.depth$depth0m[order(sd.t.depth$Lat, sd.t.depth$Long)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, sdSST.1deg))

rm(mean.t.depth, sd.t.depth, temp.margo)


## 5. Mixed layer depth ------------------------------------------------------------
setwd("../MLD_montegut/")
load("150414_mld.RData")

# identify lat / long
mld.long <- sapply(ldg.p.margo$Longitude, match.2deg, unique(mld.2deg$Long))
mld.lat <- sapply(ldg.p.margo$Latitude, match.2deg, unique(mld.2deg$Lat))

# add columns
ldg.p.margo <- cbind(ldg.p.margo, mld.2deg[mld.long[1, ] + length(unique(mld.2deg$Long)) * (mld.lat[1, ] - 1), grep("\\.mld", names(mld.2deg))])

with(ldg.p.margo, distrib.map(Longitude, Latitude, mean.mld.t))

rm(mld.long, mld.lat, mld.2deg, mld.margo)


## 6. depth10deg ----------------------------------------------------------
setwd("../SST_1deg_Mean/")
load("150414_Strat.RData")

# check that the lat / long match in length
sum(strat.1deg$Long[order(strat.1deg$Lat, strat.1deg$Long)] != ldg.p.margo$Long)
sum(strat.1deg$Lat[order(strat.1deg$Lat, strat.1deg$Long)] != ldg.p.margo$Lat)

# add a column for mean SST
ldg.p.margo$depth10deg <- strat.1deg$depth10deg[order(strat.1deg$Lat, strat.1deg$Long)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, depth10deg))

# set the ocean points depth10deg NAs to 0
with(ldg.p.margo[is.na(ldg.p.margo$depth10deg),], distrib.map(Longitude, Latitude, Longitude))
ldg.p.margo$depth10deg[is.na(ldg.p.margo$depth10deg)] <- 0

# identify points on land
ldg.p.margo$depth10deg[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, land$x, land$y) == 1)] <- NA

with(ldg.p.margo[is.na(ldg.p.margo$depth10deg),], distrib.map(Longitude, Latitude, Longitude))

rm(strat.1deg)


## 7. Productivity ---------------------------------------------------------
setwd("../Productivity/")

## 7i. Extract data for each site ---------------------------------------
lon <- seq(-180, 179 + 5/6, 1/6)
lat <- seq(90, -89 - 5/6, -1/6)

# calculate the closest coordinates for each p site
long.p <- sapply(ldg.p.margo$Longitude, match.4km, lon)
lat.p <- sapply(ldg.p.margo$Latitude, match.4km, lat)

# set up a dataframe to hold the values
prod.p <- ldg.p.margo[,1:2]
# for 2002, we don't have complete data, therefore
j <- 2
for (i in 7:12) {
  # read in the file
  if (i < 10) {
    file.nm <- paste("data_by_month/npp.0", j, "0", i, ".h5", sep = "")
  } else {
    file.nm <- paste("data_by_month/npp.0", j, i, ".h5", sep = "")
  }
  prod <- as.vector(h5read(file.nm, "npp"))
  # extract the data for the p sites
  prod.p <- cbind(prod.p, prod[long.p + length(lon) * (lat.p - 1)])
  colnames(prod.p)[ncol(prod.p)] <- paste("prod.", j, ".", i, sep = "")
}

# for each year between
for (j in 3:14) {
  # for each month
  for (i in 1:12) {
    # read in the file
    if (i < 10) {
      file.nm <- paste("data_by_month/npp.0", j, "0", i, ".h5", sep = "")
    } else {
      file.nm <- paste("data_by_month/npp.0", j, i, ".h5", sep = "")
    }
    prod <- as.vector(h5read(file.nm, "npp"))
    # extract the data for the p sites
    prod.p <- cbind(prod.p, prod[long.p + length(lon) * (lat.p - 1)])
    colnames(prod.p)[ncol(prod.p)] <- paste("prod.", j, ".", i, sep = "")
  }
}

# for 2015
j <- 15
i <- 1
file.nm <- paste("data_by_month/npp.0", j, "0", i, ".h5", sep = "")
prod <- as.vector(h5read(file.nm, "npp"))
# extract the data for the p sites
prod.p <- cbind(prod.p, prod[long.p + length(lon) * (lat.p - 1)])
colnames(prod.p)[ncol(prod.p)] <- paste("prod.", j, ".", i, sep = "")

rm(i, j, file.nm, prod, lat, lon, lat.p, long.p)
# save this out so that I can come back to it if necessary
save(prod.p, file = "150423_prod.p_working.RData")

## convert to SST
# -9999 is missing data, set to NA
for (i in 3:ncol(prod.p)) {
  prod.p[which(prod.p[, i] == -9999), i] <- NA
}
rm(i)
# check this has worked
head(prod.p)

## 7ii. Calculate monthly averages across all years ------------------------
prod.month.p <- ldg.p.margo[,1:2]

for (i in 1:12) {
  prod.month.p <- cbind(prod.month.p, NA)
  prod.month.p[ncol(prod.month.p)] <- rowMeans(prod.p[, grep(paste("\\.", i, "$", sep = ""), names(prod.p))], na.rm = TRUE)
  names(prod.month.p)[ncol(prod.month.p)] <- paste("prod.mn.", i, sep = "")
}
rm(i)

for (i in 1:12) {
  prod.month.p <- cbind(prod.month.p, NA)
  prod.month.p[ncol(prod.month.p)] <- apply(prod.p[, grep(paste("\\.", i, "$", sep = ""), names(prod.p))], 1, sd, na.rm = T)
  names(prod.month.p)[ncol(prod.month.p)] <- paste("prod.sd.", i, sep = "")
}
rm(i)

## 7iii. Plot these productivity data ----------------------------------------
# the data follows an exponential distribution, so plot < 1000 to see the major patterns
# plot these up
for (i in 1:12) 
{
  col <- grep(paste(".mn.", i, "$", sep = ""), names(prod.month.p))
  sites1000 <- prod.month.p[, col] < 1000
  png(paste("prodMn_", i, "_1000.png", sep = ""), 800, 500)
  with(prod.month.p[sites1000, ], distrib.map(Longitude, Latitude, prod.month.p[sites1000, col], palette = "matlab.like", col.land = "black", col.water = "white"))
  dev.off()
  col <- grep(paste(".sd.", i, "$", sep = ""), names(prod.month.p))
  sites1000 <- prod.month.p[, col] < 1000
  png(paste("prodSd_", i, "_1000.png", sep = ""), 800, 500)
  with(prod.month.p[sites1000, ], distrib.map(Longitude, Latitude, prod.month.p[sites1000, col], palette = "matlab.like", col.land = "black", col.water = "white"))
  dev.off()
}
rm(i, col, sites1000)

## 7iv. Calculate mean and sd for these ------------------------------------
prod.month.p$prod.mn.ann <- rowMeans(prod.month.p[, grep("mn", names(prod.month.p))], na.rm = T)
prod.month.p$prod.sd.ann <- rowMeans(prod.month.p[, grep("sd", names(prod.month.p))], na.rm = T)

# due to the spread of productivity, also calculate logs
prod.month.p$logProd.mn.ann <- rowMeans(log(prod.month.p[, grep("mn", names(prod.month.p))][-13]), na.rm = T)
prod.month.p$logProd.sd.ann <- rowMeans(log(prod.month.p[, grep("sd", names(prod.month.p))][-13]), na.rm = T)

# plot these
for (i in grep("ann", names(prod.month.p))) {
  png(paste(names(prod.month.p)[i], ".png", sep = ""), width = 800, height = 500)
  with(prod.month.p[prod.month.p$prod.mn.ann < 1000, ], distrib.map(Longitude, Latitude, prod.month.p[prod.month.p$prod.mn.ann < 1000, i], palette = "matlab.like", col.land = "black", col.water = "white"))
  dev.off()
}

## 7v. add these to ldg.p.margo -----------------------------------------------
ldg.p.margo$prod.mn.ann <- prod.month.p$prod.mn.ann
ldg.p.margo$prod.sd.ann <- prod.month.p$prod.sd.ann
ldg.p.margo$logProd.mn.ann <- prod.month.p$logProd.mn.ann
ldg.p.margo$logProd.sd.ann <- prod.month.p$logProd.sd.ann

# tidy up
save(prod.p, prod.month.p, file = "150423_productivity_p.RData")
rm(prod.p, prod.month.p)


## 8. Salinity ----------------------------------------------------------
setwd("../Salinity")
load("150414_salinity.RData")

## 8i. mean Salinity -------------------------------------------------------
head(sal.mean.depth)
ldg.p.margo$meanSal.0m <- NA
for (i in 1:nrow(sal.mean.depth)) {
  ldg.p.margo$meanSal.0m[sal.mean.depth$Longitude[i] == ldg.p.margo$Longitude & sal.mean.depth$Latitude[i] == ldg.p.margo$Latitude] <- sal.mean.depth$Depth0m[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, meanSal.0m))

## 8ii. sd Salinity -------------------------------------------------------
ldg.p.margo$sdSal.0m <- NA
for (i in 1:nrow(sal.sd.depth)) {
  ldg.p.margo$sdSal.0m[sal.sd.depth$Longitude[i] == ldg.p.margo$Longitude & sal.sd.depth$Latitude[i] == ldg.p.margo$Latitude] <- sal.sd.depth$Depth0m[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, sdSal.0m))

rm(sal.mean.depth, sal.sd.depth, sal.margo)


## 9. Oxygen stress --------------------------------------------------------
setwd("../Oxygen stress/")
load("150506_oxygen.RData")

# mean oxygen
ldg.p.margo$meanOxy <- NA
for (i in 1:nrow(oxy.1deg)) {
  ldg.p.margo$meanOxy[oxy.1deg$Longitude[i] == ldg.p.margo$Longitude & oxy.1deg$Latitude[i] == ldg.p.margo$Latitude] <- oxy.1deg$meanOxy[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, meanOxy))

# sd oxygen
ldg.p.margo$sdOxy <- NA
for (i in 1:nrow(oxy.1deg)) {
  ldg.p.margo$sdOxy[oxy.1deg$Longitude[i] == ldg.p.margo$Longitude & oxy.1deg$Latitude[i] == ldg.p.margo$Latitude] <- oxy.1deg$sdOxy[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, sdOxy))

# prop oxygen
ldg.p.margo$prop2.oxy <- NA
for (i in 1:nrow(oxy.1deg)) {
  ldg.p.margo$prop2.oxy[oxy.1deg$Longitude[i] == ldg.p.margo$Longitude & oxy.1deg$Latitude[i] == ldg.p.margo$Latitude] <- oxy.1deg$prop2.oxy[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, prop2.oxy))

rm(oxy.1deg, oxy.margo)


## 10. delta_carb_ion ------------------------------------------------------------
# carb_ion
ldg.p.margo$delta_carb_ion <- 0
# pick this as high values of carbonate ion saturation indicate no dissolution

ldg.p.margo$delta_carb_ion[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, land$x, land$y) == 1)] <- NA

with(ldg.p.margo[!is.na(ldg.p.margo$delta_carb_ion), ], distrib.map(Longitude, Latitude, Longitude))


## 11. Save the data -------------------------------------------------------
save(ldg.p.margo, file = "../../../../Work/1311 LDGPaper/Reanalysis/Outputs/ldg_p_margo.RData")

rm(land)
