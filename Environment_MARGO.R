## 26 / 3 / 2015
## Isabel Fenton
## Collecting the Environmental data from the ldg analysis for the MARGO points
##
## Based on 1311 LDGPaper/Code/Environment.R
##
## Input: 
## 150326_Environment.RData - the environmental data already downloaded for the bfd

## BFD_SST.txt (BFD/Environmental/SST_4km) temperature at 4km res
## files of mean and sd of 1deg temperature at various depths (BFD/Environmental/SST_1deg_Mean & SST_1deg_SD)
## MLD_dec_ed.txt (BFD/Environmental/MLD) mixed layer depth yearly values
## monthly_chl.txt (BFD/Environmental/Chl_4km) chlorophyll data at 4km resolution
## files of Chla_1deg_means.csv (BFD/Environmental/Chl_1deg) chlorophyll data at 1 degree resolution
## s00mn01.csv (BFD/Environment/Salinity) mean salinity with depth data at 1 degree resolution
## s00sd01.csv (BFD/Environment/Salinity) sd salinity with depth data at 1 degree resolution
##
## Output: 


## distrib.maps of all files 
## temperature sections with depth for range of stratification measures
## 140522_ldg_env.Rdata - environmental data for BFD sites
## 140522_Environment_ws.Rdata - Complete workspace
##

source("C:/Documents/Science/PhD/Code/maps.R")
source("C:/Documents/Science/PhD/Code/palettes.R")
source("C:/Documents/Science/PhD/Code/sp_mat_2_df.R") # for section 6i. 
library(sp) # point in polygon, section 6ii. 
load("C:/Documents/Science/PhD/Work/1311\ LDGPaper/Output/150326_Environment.RData")

## 1. SST at 4km resolution ------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/SST_4km/")

## 1i. load in the data ----------------------------------------------------
# This data has been obtained from NOAA. See temperature.docx and SSTmeta.txt for more info
SST.4km.bfd <- read.table("BFD_SST.txt", header = T)

## convert to SST
# scale*pixel_Value + offset
# scale = 0.075 and offset = - 3.0
SST.4km.bfd[, 4:15] <- SST.4km.bfd[, 4:15] * 0.075 - 3

## 1ii. calculate mean and SD ----------------------------------------------
# add columns of mean and SD
SST.4km.bfd$meanSST <- rowMeans(SST.4km.bfd[, 4:15])
SST.4km.bfd$sdSST <- apply(SST.4km.bfd[, 4:15], 1, sd)

## 1iii. output images ------------------------------------------------------
# notice there is one point in the tropical Atlantic which is coming out very low T
order(SST.4km.bfd$meanSST)[1]
(sort(SST.4km.bfd$meanSST)[1] + 3) / 0.075
which(SST.4km.bfd$sdSST == 0)

# for some reason this hasn't worked, so set it to NA
SST.4km.bfd$meanSST[SST.4km.bfd$meanSST == min(SST.4km.bfd$meanSST)] <- NA
SST.4km.bfd$sdSST[SST.4km.bfd$sdSST == 0] <- NA

png(file = "meanSST_4km.png", width = 800, height = 500)
with(SST.4km.bfd[!is.na(SST.4km.bfd$meanSST), ], distrib.filled(interp(LONG, LAT, meanSST + 0.5, duplicate = "mean"), nlevels = 100))
dev.off()

png(file = "sdSST_4km.png", width = 800, height = 500)
with(SST.4km.bfd[!is.na(SST.4km.bfd$sdSST), ], distrib.filled(interp(LONG, LAT, sdSST, duplicate = "mean"), nlevels = 100))
dev.off()

## 1iv. put values for mean and sd in a dataframe ----------------------------------------
ldg.env <- ldg.data[, 1:6]

# add that data to ldg.data
ldg.env$meanSST.4km <- SST.4km.bfd$meanSST
ldg.env$sdSST.4km <- SST.4km.bfd$sdSST


## 2. Temperature with depth at 1deg resolution ----------------------------

## 2i. Create a function for loading data and one for images -------------------
# actual depths
WOA.depths <- c(0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)

## function to load the data in for mean / sd and to produce a dataframe
load.depths <- function(value){
  # value - either "mn" or "sd"
  # for each depth
  for (i in 1:33) {
    # load in the file
    if (i < 10) {
      tmp <- read.table(paste("t00", value, "1.00", i, sep = ""), header = T)
    } else {
      tmp <- read.table(paste("t00", value, "1.0", i, sep = ""), header = T)
    }
    # for the first load the entire file, for the rest, only add the value column
    if (i == 1) {
      t1 <- tmp
      colnames(t1)[1:2] <- c("Lat", "Long")
    } else {
      t1 <- cbind(t1, tmp$WOA_VALUE)
    }
    # change column names to relevant depth
    colnames(t1)[i + 2] <- paste("depth", WOA.depths[i], "m", sep = "")
    # change missing values to NA
    t1[,i + 2][t1[,i + 2] == -99.999] <- NA  
  }
  return(t1)
}

# function to produce image files for each depth
data.img <- function(data, name) {
  # data - dataframe of temp data with depth
  # name - either "mn" or "sd"
  for (i in 1:33) {
    # produce .png images in the same file  
    png(paste("t", name, i,  ".png", sep = ""), width = 800, height = 500)
    # n.b. the scale is offset by 3 (so currently no key)
    with(data, distrib.map(Long, Lat, data[,i + 2], pch = 15, maintitle = colnames(data)[i + 2], cex = 0.4))
    dev.off()
  }
}

## 2ii. Run function for mean T ----------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/SST_1deg_Mean")
mean.t.depth <- load.depths("mn")
data.img(mean.t.depth, "mn")

## 2iii. Run function for SD --------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/SST_1deg_SD")
sd.t.depth <- load.depths("sd")
data.img(sd.t.depth, "sd")

## 2iv. Identifying environmental conditions for BFD points ----------------------
# write a function that will resolve a value to a given accuracy
matching <- function(coord, to.match) {
  # coord - the latitude or longitude of the site (single number)
  # to.match - vector of the new resolutions
  # give it data to a certain resolution, return data to a new resolution
  # what is the new resolution?
  new.res <- unique(to.match %% 1)
  round.coord <- coord %% 1
  # set up the resolution comparisons
  new.res <- c(max(new.res) - 1, new.res, min(new.res) + 1)
  # return the new resolution ending
  new.coord <- new.res[which(abs(new.res - round.coord) == min(abs(new.res - round.coord)))][1]
  new.coord + coord %/% 1
}

# a function that will return the temperature data for a given set of coordinates
WOA.temp <- function(bfd.row, T.data, T.cols) {
  # bfd.row - a row from a sites dataframe, with columns 'Lat' and 'Long
  # T.data - an environmental dataframe, with columns 'Lat' and 'Long'
  # T.data - the columns of the environmental dataframe which relate to the environmental variable
  # calculate the possible lats and longs
  long <- matching(bfd.row$Long, T.data$Long)
  lat <- matching(bfd.row$Lat, T.data$Lat)
  # which row has both 
  row <- which(T.data$Long == long & T.data$Lat == lat)
  # if there is no data throw an error
  # return either a set of values or single value
  if (!is.null(dim(T.cols))) {
    if (sum(!is.na(T.cols[row,])) == 0) stop ("All data are NAs")
    return(T.cols[row,])
  } else {
    if (is.na(T.cols[row])) stop ("Data is NA")
    return(T.cols[row])
  }
}

# generate a dataframe containing the coordinates to calculate temperature at those points. 
temp.bfd <- ldg.data[,1:6]
# create columns for mean t
for (i in 3:ncol(mean.t.depth)){
  temp.bfd <- cbind(temp.bfd,NA)
  colnames(temp.bfd)[4 + i] <- as.character(paste("mean", colnames(mean.t.depth[i]), sep = "."))
}
# and for sd
for (i in 3:ncol(sd.t.depth)){
  temp.bfd <- cbind(temp.bfd,NA)
  colnames(temp.bfd)[37 + i] <- as.character(paste("sd", colnames(sd.t.depth[i]), sep = "."))
}

# add the temperature data (mean and sd) for the whole of the dataset
for(i in 1:nrow(temp.bfd))
{
  tryCatch(temp.bfd[i, grep("mean", colnames(temp.bfd))] <- WOA.temp(temp.bfd[i,], mean.t.depth, mean.t.depth[, 3:ncol(mean.t.depth)]), error = function(e) print(i))
  tryCatch(temp.bfd[i, grep("sd", colnames(temp.bfd))] <- WOA.temp(temp.bfd[i,], sd.t.depth, sd.t.depth[, 3:ncol(sd.t.depth)]), error = function(e) print(i))
}

# identify where there is no data
which(is.na(temp.bfd$mean.depth0m))
which(is.na(temp.bfd$sd.depth0m))

# for those sites average the neighbouring sites
WOA.na.temp <- function(bfd.row, T.data, T.cols) {
  # bfd.row - a row from a sites dataframe, with columns 'Lat' and 'Long
  # T.data - an environmental dataframe, with columns 'Lat' and 'Long'
  # T.data - the columns of the environmental dataframe which relate to the environmental variable
  # calculate the possible lats and longs
  long <- matching(bfd.row$Long, T.data$Long)
  lat <- matching(bfd.row$Lat, T.data$Lat)
  # generate a data frame of nearest points
  nearest <- data.frame(long = c(long, long + 1, long - 1, long, long), lat = c(rep(lat, 3), lat + 1, lat - 1))
  # cycle through those points to look for data
  tmp.values <- NULL
  if (!is.null(dim(T.cols))) {
    for (i in 1:nrow(nearest)) {
      row <- which(T.data$Long == nearest$long[i] & T.data$Lat == nearest$lat[i])
      tmp.values <- rbind(tmp.values, T.cols[row, ])
    }
    # calculate the mean values
    avg.value <- apply(tmp.values, 2, mean, na.rm = T)
    avg.value[is.na(avg.value)] <- NA
    # if all data is NA, throw an error
    if (sum(!is.na(avg.value)) == 0) stop ("All data are NAs")
    return(avg.value)
  } else {
    for (i in 1:nrow(nearest)) {
      row <- which(T.data$Long == nearest$long[i] & T.data$Lat == nearest$lat[i])
      tmp.values <- c(tmp.values, T.cols[row])
    }
    return(mean(tmp.values, na.rm = T))
  }
}

na.sites <- sort(unique(c(which(is.na(temp.bfd$mean.depth0m)), which(is.na(temp.bfd$sd.depth0m)))))

# add a column to the dataframe showing which are exact
temp.bfd$exact <- "Y"
temp.bfd$exact[na.sites] <- "N"

# add the data (mean and sd) for those missing values by averaging surrounding sites
for(i in na.sites)
{
  tryCatch(temp.bfd[i, grep("mean", colnames(temp.bfd))] <- WOA.na.temp(temp.bfd[i,], mean.t.depth, mean.t.depth[, 3:ncol(mean.t.depth)]), error = function(e) print(i))
  
  tryCatch(temp.bfd[i, grep("sd", colnames(temp.bfd))] <- WOA.na.temp(temp.bfd[i,], sd.t.depth, sd.t.depth[, 3:ncol(sd.t.depth)]), error = function(e) print(i))
}

rm(na.sites)

# check whether this has worked
which(is.na(temp.bfd$mean.depth0m)) # there now aren't any
which(is.na(temp.bfd$sd.depth0m)) # or of these

# add the data to the ldg.env dataframe
ldg.env$meanSST.1deg <- temp.bfd$mean.depth0m
ldg.env$sdSST.1deg <- temp.bfd$sd.depth0m
ldg.env$SST.1deg.exact <- temp.bfd$exact


## 3. Calculating MLD data -------------------------------------------------

## 3i. Load in the MLD data ------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/MLD/")
mld.2 <- read.table("MLD_dec_ed.txt", header = T)

## 3ii. Calculate the mean and SD ------------------------------------------
# generate a dataframe with a row for each site (MLD_Apr is arbitrary choice)
names(mld.2)
mld.ann <- mld.2[mld.2$Cruise == "MLD_Apr", 2:8]
dim(mld.ann)
dim(mld.2)

# rename colnames in mld.ann to "Lat" and "Long"
names(mld.ann)
names(mld.ann) <- gsub("Latitude", "Lat", names(mld.ann))
names(mld.ann) <- gsub("Longitude", "Long", names(mld.ann))
names(mld.ann)

# add mean values for the year for potential temperature
# n.b. Station is unique for each site
mld.ann$mean.pt <- unlist(lapply(mld.ann$Station, function (x) mean(mld.2$MLD.pt.m.[mld.2$Station == x])))
# there are 110 NAs, but these are from sites that have no MLD data
sum(is.na(mld.2$MLD.pt.m.))/12

# add sd values for the year
mld.ann$sd.pt <- unlist(lapply(mld.ann$Station, function (x) sd(mld.2$MLD.pt.m.[mld.2$Station == x])))

# add for the other measures of mld (potential density and variable density)
mld.ann$mean.pd <- unlist(lapply(mld.ann$Station, function (x) mean(mld.2$MLD.pd.m.[mld.2$Station == x])))
mld.ann$sd.pd <- unlist(lapply(mld.ann$Station, function (x) sd(mld.2$MLD.pd.m.[mld.2$Station == x])))
mld.ann$mean.vd <- unlist(lapply(mld.ann$Station, function (x) mean(mld.2$MLD.vd.m.[mld.2$Station == x])))
mld.ann$sd.vd <- unlist(lapply(mld.ann$Station, function (x) sd(mld.2$MLD.vd.m.[mld.2$Station == x])))

# check these are giving what I think they do
head(mld.ann)
mld.2[mld.2$Station == 1,]
mean(mld.2$MLD.pt.m.[mld.2$Station == 1]) == mld.ann$mean.pt[mld.ann$Station == 1]
sd(mld.2$MLD.pt.m.[mld.2$Station == 1]) == mld.ann$sd.pt[mld.ann$Station == 1]
mean(mld.2$MLD.pd.m.[mld.2$Station == 1]) == mld.ann$mean.pd[mld.ann$Station == 1]
sd(mld.2$MLD.pd.m.[mld.2$Station == 1]) == mld.ann$sd.pd[mld.ann$Station == 1]
mean(mld.2$MLD.vd.m.[mld.2$Station == 1]) == mld.ann$mean.vd[mld.ann$Station == 1]
sd(mld.2$MLD.vd.m.[mld.2$Station == 1]) == mld.ann$sd.vd[mld.ann$Station == 1]

# these are all true, so that has worked

## 3iii. Plot the data -----------------------------------------------------
tmp <- c("mean.pt", "sd.pt", "mean.pd", "sd.pd", "mean.vd", "sd.vd")
# plot the environmental data
for (i in tmp) {
  png(file = paste(i, ".png", sep = ""), width = 800, height = 500)
  with(mld.ann[!is.na(mld.ann[,i]), ], distrib.map(Long, Lat, mld.ann[,i][!is.na(mld.ann[,i]) ], pch = 15, cex = 0.4))
  dev.off()
}

rm(i, tmp)

# but note that the max values are not the same as the max on the scale
max(mld.ann$mean.pt, na.rm = T)
plot(sort(mld.ann$mean.pt))
plot(sort(log(mld.ann$mean.pt)))

# so log the data and then replot
for (i in tmp) {
  png(file = paste(i, "_log.png", sep = ""), width = 800, height = 500)
  with(mld.ann[!is.na(mld.ann[,i]), ], distrib.map(Long, Lat, log(mld.ann[,i][!is.na(mld.ann[,i]) ] + 1), pch = 15, cex = 0.4))
  dev.off()
}

rm(i, tmp)

## 3iv. Obtain values for BFD points ---------------------------------------
# create a dataframe for the data
mld.bfd <- ldg.data[,1:6]

# create columns for mld data
for (i in 8:ncol(mld.ann)){
  mld.bfd <- cbind(mld.bfd,NA)
  colnames(mld.bfd)[i - 1] <- colnames(mld.ann[i])
}

head(mld.bfd)

# not all rows are present, so obtain data for those that do and identify those that don't 
for(i in 1:nrow(mld.bfd))
{
  tryCatch(mld.bfd[i, 7:ncol(mld.bfd)] <- WOA.temp(mld.bfd[i,], mld.ann, mld.ann[, 8:ncol(mld.ann)]), error = function(e) print(i))
}

# check that this function works for the first one
mld.bfd[1,]
mld.ann[mld.ann$Lat == 27.5 & mld.ann$Long == -38.5, ]

# seems only three don't have data
which(is.na(mld.bfd$mean.pt))

# create a list of those sites where it hasn't worked
na.sites <- which(is.na(mld.bfd$mean.pt))

# add a column to the dataframe showing which are exact
mld.bfd$exact <- "Y"
mld.bfd$exact[na.sites] <- "N"

table(mld.bfd$exact)

# add the salinity data (mean and sd) for the whole of the dataset
for(i in na.sites)
{
  tryCatch(mld.bfd[i, 7:(ncol(mld.bfd) - 1)] <- WOA.na.temp(mld.bfd[i,], mld.ann, mld.ann[, 8:ncol(mld.ann)]), error = function(e) print(i))
}

rm(na.sites)

# check there are now no na's
which(is.na(mld.bfd$mean.pt))

# compare this with the mld.ann to see if it looks reasonable
with(mld.bfd, distrib.filled(interp(Long, Lat, mean.pt, duplicate = "mean")))
with(mld.ann[!is.na(mld.ann$mean.pt), ], distrib.filled(interp(Long, Lat, mean.pt, duplicate = "mean")))

ldg.env <- cbind(ldg.env, mld.bfd[,7:ncol(mld.bfd)])
names(ldg.env)[names(ldg.env) == "exact"] <- "mld.exact"


## 4. Stratification -------------------------------------------------------

## 4i. A function for plotting depth gradients ------------------------------
depth.grad <- function(env.data, depth.data, x.axis = "", main = "", ...)
{
  # env.data - data.frame row, the environmental data
  # depth.data - vector, the depths
  # x.axis - string
  # main - string
  # plot a temperature by depth graph
  plot(as.numeric(env.data), -depth.data, bty = "l", xlab = x.axis, ...)
}

depth.grad(temp.bfd[1, grep("mean", names(temp.bfd))], WOA.depths)

## 4ii. plotting latitude depth profiles  -----------------------------------
# currently plotting at 0.5 deg Long
filled.contour(mean.t.depth$Lat[1:180], WOA.depths, as.matrix(mean.t.depth[1:180, 3:ncol(mean.t.depth)]), ylim =  rev(range(WOA.depths)), zlim = range(mean.t.depth[1:180,3:ncol(mean.t.depth)], na.rm = T), color.palette = log.heat)
# for the top 1500m
filled.contour(mean.t.depth$Lat[1:180], WOA.depths, as.matrix(mean.t.depth[1:180, 3:ncol(mean.t.depth)]), ylim = c(1500, 0), zlim = range(mean.t.depth[1:180,3:ncol(mean.t.depth)], na.rm = T), color.palette = log.heat)

## 4iii. Calculate 10deg contour for temp.bfd -------------------------------

# create a function for identifying the 10 degree contour for a given longitudinal transect
calc.10deg <- function(row, data, cols) {
  # row - the row of the dataset
  # data - dataframe of temperature with depth
  # cols - a list of the columns containing the environmental data
  WOA.depths <- c(0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)
  tmp <- data[row, cols]
  # check row has data
  if(sum(!is.na(tmp)) > 0) {
    # check it includes the 10degree contour
    if(max(tmp, na.rm = T) > 10 & min(tmp, na.rm = T) < 10) {
      # if it does calculate the depth closest to 10
      return(WOA.depths[which(abs(tmp - 10) == min(abs(tmp - 10), na.rm = T))[1]])
    }
  }
  return(NA)
}

# test this for Long = 0.5
tmp.depth <- rep(NA, 180)
tmp.lat <- rep(NA, 180)
for(i in 1:180) {
  tmp.depth[i] <- calc.10deg(i, mean.t.depth, 3:ncol(mean.t.depth))
  tmp.lat[i] <- mean.t.depth$Lat[i]
}
plot(tmp.lat, tmp.depth, xlim = c(-90, 90), ylim = c(2000, 0), pch = 20)
rm(tmp.depth, tmp.lat)

# generate a dataframe for adding stratification variables
strat <- mean.t.depth[, 1:2]

# calculate the 10 deg depth for all locations
strat$depth10deg <- NA
# for each line of longitude
for (i in 1:length(unique(mean.t.depth$Long))) {
  # for each degree of latitude
  for(j in 1:180) {
    # check that the T range includes 10
    strat$depth10deg[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- calc.10deg(j + 180 * (i - 1), mean.t.depth, 3:ncol(mean.t.depth))
  }
}

# check this has worked
plot(strat$Lat[1:180], strat$depth10deg[1:180], xlim = c(-90, 90), ylim = c(2000, 0), pch = 20)

# cycle through these with the contour plot
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, depth10deg, pch = 20) }))
}

rm(i, tmp)

# Look at the spatial variation of this
png("depth10deg.png", width = 800, height = 500)
with(strat, distrib.map(Long, Lat, depth10deg, palette = "rev.log.heat", pch = 15, cex = 0.4))
dev.off()

# calculate the BFD values
temp.bfd$depth10deg <- rep(NA, nrow(temp.bfd))
for (i in 1:nrow(temp.bfd)) {
  tryCatch(temp.bfd$depth10deg[i] <- WOA.temp(temp.bfd[i,], strat, strat$depth10deg), error = function(e) print(i))
}

summary(temp.bfd$depth10deg)
with(temp.bfd,distrib.map(Long, Lat, depth10deg))
with(temp.bfd,distrib.map(Long, Lat, depth10deg, palette = "rev.log.heat"))

# where are the NAs
with(temp.bfd[is.na(temp.bfd$depth10deg),], distrib.map(Long, Lat, 2))

ldg.env$depth10deg <- temp.bfd$depth10deg

## 4iv. Thermocline as fastest rate of change below MLD ------------------------------------
# create a function to calculate this
calc.thermo.rate <- function(num) {
  WOA.depths <- c(0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)
  lat <- mean.t.depth$Lat[num]
  long <- mean.t.depth$Long[num]
  env <- mean.t.depth[num, 3:ncol(mean.t.depth)]
  
  # calculate thermocline depth for all these
  # thermocline is maximum rate of change of Temp with depth below the mld
  depth <- NA
  # remove depths shallower than MLD
  mld.depth <- mld.ann$mean.pt[mld.ann$Lat == lat & mld.ann$Long == long]
  # if no mld data for that site
  if (length(mld.depth) == 0 || is.na(mld.depth))
    return(depth)
  env.diff <- diff(as.numeric(env[WOA.depths > mld.depth]))
  short.depths <- WOA.depths[WOA.depths > mld.depth]
  # calculate thermocline depth
  if (sum(!is.na(env.diff)) > 0) {
    # sometimes get two depths with the same rate of change - select the first
    num <- which(env.diff / diff(short.depths) == min(env.diff / diff(short.depths), na.rm = T))[1]
    # as we are calculating diff, should return an average of the two depths
    depth <- mean(c(short.depths[num], short.depths[num + 1]))
  }
  return(depth)
}

# calculate this for all points
strat$thermo.rate <- NA
# for each line of longitude
for (i in 1:length(unique(mean.t.depth$Long))) {
  # for each degree of latitude
  for(j in 1:180) {
    # check that the T range includes 10
    strat$thermo.rate[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- calc.thermo.rate(j + 180 * (i - 1))
  }
}

# check this
summary(strat$thermo.rate)

# cycle through these with the contour plot
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thermo.rate, pch = 20) }))
}

rm(i, tmp)

# Look at the spatial variation of this
with(strat, distrib.map(Long, Lat, thermo.rate, palette = "rev.log.heat"))

# calculate the BFD values
temp.bfd$thermo.rate <- rep(NA, nrow(temp.bfd))
for (i in 1:nrow(temp.bfd)) {
  tryCatch(temp.bfd$thermo.rate[i] <- WOA.temp(temp.bfd[i,], strat, strat$thermo.rate), error = function(e) print(i))
}

# which have no values
which(is.na(temp.bfd$thermo.rate))

# create a list of those sites
na.sites <- which(is.na(temp.bfd$thermo.rate))

# add a column to the dataframe showing which are exact
temp.bfd$thra.exact <- "Y"
temp.bfd$thra.exact[na.sites] <- "N"

# calculate means for sites where there are NA's. 
for(i in na.sites)
{
  tryCatch(temp.bfd$thermo.rate[i] <- WOA.na.temp(temp.bfd[i,], strat, strat$thermo.rate), error = function(e) print(i))
}

rm(na.sites)

# look at this
summary(temp.bfd$thermo.rate)
with(temp.bfd,distrib.map(Long, Lat, temp.bfd$thermo.rate))
with(temp.bfd,distrib.map(Long, Lat, temp.bfd$thermo.rate, palette = "rev.log.heat"))

# add this data to the bfd dataframe
ldg.env$thermo.rate <- temp.bfd$thermo.rate
ldg.env$thra.exact <- temp.bfd$thra.exact

# points should be below MLD (as forced), but check
plot(temp.bfd$thermo.rate, ldg.env$mean.pt)
abline(1,1)
# also how where it is most different from MLD
with(ldg.env, distrib.map(Long, Lat, thermo.rate - mean.pt, palette = "rev.log.heat"))

# look at those which are deepest
rev(order(temp.bfd$thermo.rate))[1:10]
plot(-WOA.depths ~ as.numeric(temp.bfd[954, 7:39]), type = "l")
abline(h = -ldg.env$mean.pt[954])
abline(h = -temp.bfd$thermo.rate[954])

# can check multiple of these plots
par(ask = T)
for (i in 1:nrow(temp.bfd)) {
  plot(-WOA.depths ~ as.numeric(temp.bfd[i, 7:39]), type = "l", main = i)
  abline(h = -ldg.env$mean.pt[i])
  abline(h = -temp.bfd$thermo.rate[i])
}
par(ask = F)

## 4v. Depth of temp halfway between MLD and zero (mid-thermocline) -------
# create a function to calculate this
calc.thermo.depth <- function(num) {
  WOA.depths <- c(0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)
  lat <- mean.t.depth$Lat[num]
  long <- mean.t.depth$Long[num]
  env <- mean.t.depth[num, 3:ncol(mean.t.depth)]
  
  # calculate thermocline depth for all these
  # thermocline is maximum rate of change of Temp with depth below the mld
  depth <- NA
  # remove depths shallower than MLD
  mld.depth <- mld.ann$mean.pt[mld.ann$Lat == lat & mld.ann$Long == long]
  # if no mld data for that site
  if (length(mld.depth) == 0 || is.na(mld.depth))
    return(depth)
  env.mld <- as.numeric(env[WOA.depths > mld.depth])
  short.depths <- WOA.depths[WOA.depths > mld.depth]
  # calculate thermocline depth
  if (sum(!is.na(env.mld)) > 0 & !is.na(env.mld[1])) {
    # identify the temperature at the mld
    # assume bottom water at 0 degrees, so calculate depth of halfway between mld and 0
    mid.temp <- env.mld[1]/2  
    # stop if mid.temp isn't within the temp range
    if (sum(env.mld - mid.temp < 0, na.rm = T) != 0) {
      # find midpoint of thermocline
      mid.val <- which(abs(env.mld - mid.temp) == min(abs(env.mld - mid.temp), na.rm = T))
      # calc mean so that if there are multiple (i.e. mid is exactly between 2), then get the average depth
      depth <- mean(short.depths[mid.val])
    }
  }
  return(depth)
}

# calculate this for all points
strat$thermo.depth <- NA
# for each line of longitude
for (i in 1:length(unique(mean.t.depth$Long))) {
  # for each degree of latitude
  for(j in 1:180) {
    # check that the T range includes 10
    strat$thermo.depth[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- calc.thermo.depth(j + 180 * (i - 1))
  }
}

# check this
summary(strat$thermo.depth)

# cycle through these with the contour plot
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thermo.depth, pch = 20) }))
}

rm(i, tmp)

# Look at the spatial variation of this
with(strat, distrib.map(Long, Lat, thermo.depth, palette = "rev.log.heat"))

# calculate the BFD values
temp.bfd$thermo.depth <- rep(NA, nrow(temp.bfd))
for (i in 1:nrow(temp.bfd)) {
  tryCatch(temp.bfd$thermo.depth[i] <- WOA.temp(temp.bfd[i,], strat, strat$thermo.depth), error = function(e) print(i))
}

# which have no values
which(is.na(temp.bfd$thermo.depth))

# create a list of those sites
na.sites <- which(is.na(temp.bfd$thermo.depth))

# add a column to the dataframe showing which are exact
temp.bfd$thde.exact <- "Y"
temp.bfd$thde.exact[na.sites] <- "N"

# calculate the mean values for those sites which are NAs
for(i in na.sites)
{
  tryCatch(temp.bfd$thermo.depth[i] <- WOA.na.temp(temp.bfd[i,], strat, strat$thermo.depth), error = function(e) print(i))
}

rm(na.sites)

# look at this
summary(temp.bfd$thermo.depth)
with(temp.bfd,distrib.map(Long, Lat, temp.bfd$thermo.depth))
with(temp.bfd,distrib.map(Long, Lat, temp.bfd$thermo.depth, palette = "rev.log.heat"))

# add these rows to the bfd dataset
ldg.env$thermo.depth <- temp.bfd$thermo.depth
ldg.env$thde.exact <- temp.bfd$thde.exact

# points should be below MLD (as forced), but check
plot(temp.bfd$thermo.depth, ldg.env$mean.pt)
abline(1,1)
# also how where it is most different from MLD
with(ldg.env, distrib.map(Long, Lat, thermo.depth - mean.pt + 10, palette = "rev.log.heat"))

## 4vi. Mid, width, top and bottom of thermocline -------------------------
# create a function to calculate this
calc.thermocline.val <- function(num) {
  WOA.depths <- c(0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500)
  lat <- mean.t.depth$Lat[num]
  long <- mean.t.depth$Long[num]
  env <- mean.t.depth[num, 3:ncol(mean.t.depth)]
  
  # calculate thermocline depth for all these
  # thermocline is maximum rate of change of Temp with depth below the mld
  mid.depth <- NA
  width.thermo <- NA
  top.depth <- NA
  bot.depth <- NA
  # remove depths shallower than MLD
  mld.depth <- mld.ann$mean.pt[mld.ann$Lat == lat & mld.ann$Long == long]
  # if no mld data for that site
  if (length(mld.depth) == 0 || is.na(mld.depth))
    return(c(mid.depth, width.thermo, top.depth, bot.depth))
  env.mld <- as.numeric(env[WOA.depths > mld.depth])
  short.depths <- WOA.depths[WOA.depths > mld.depth]
  # calculate thermocline depth
  if (sum(!is.na(env.mld)) > 0) {
    # identify the maximum temperature
    top.temp <- max(env.mld, na.rm = T)
    # identify the deepest temperature
    deep.temp <- as.numeric(rev(env.mld[!is.na(env.mld)])[1])
    # stop if data too incomplete
    if (top.temp != deep.temp & deep.temp < 10) {
      # calculate the mid temperature
      mid.temp <- (top.temp + deep.temp)/2  
      # calculate the temp range
      t.range <- top.temp - deep.temp
      
      # top of thermocline is MLD
      top.depth <- short.depths[1]
      
      # find bottom of thermocline (10% warmer than min)
      bot.val <- max(which(env.mld - deep.temp > t.range/10), na.rm = T)
      bot.depth <- short.depths[bot.val]
      
      # width of thermocline
      width.thermo <- bot.depth - top.depth
      
      # find midpoint of thermocline
      mid.val <- which(abs(env.mld - mid.temp) == min(abs(env.mld - mid.temp), na.rm = T))
      mid.depth <- short.depths[mid.val][1]
    }
  }
  return(c(mid.depth, width.thermo, top.depth, bot.depth))
}

# check this for some random numbers
calc.thermocline.val(runif(1, 1, 360*180))

# add columns to the dataset with mid and width of the thermocline as well as top and bottom temperatures
strat$thcl.mid <- NA
strat$thcl.width <- NA
strat$thcl.top <- NA
strat$thcl.bottom <- NA

# for each line of longitude
for (i in 1:length(unique(mean.t.depth$Long))) {
  # for each degree of latitude
  for(j in 1:180) {
    tmp <- calc.thermocline.val(j + 180 * (i - 1))
    strat$thcl.mid[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- tmp[1]
    strat$thcl.width[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- tmp[2]
    strat$thcl.top[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- tmp[3]
    strat$thcl.bottom[strat$Long == unique(mean.t.depth$Long)[i] & strat$Lat == mean.t.depth$Lat[j]] <- tmp[4]
  }
}

rm(i, j, tmp)

# check for NAs and overall spread of data
summary(strat$thcl.mid)
summary(strat$thcl.width)
summary(strat$thcl.top)
summary(strat$thcl.bottom)

# cycle through these with the contour plot
# mid
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Mid : Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thcl.mid, pch = 20) }))
}

rm(i, tmp)

# width
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Width : Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thcl.width, pch = 20) }))
}

rm(i, tmp)

# MLD 
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("MLD : Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thcl.top, pch = 20) }))
}

rm(i, tmp)

# bottom
for (i in 1:length(unique(strat$Long))) {
  # matrix for the contour background
  tmp <- as.matrix(mean.t.depth[mean.t.depth$Long == unique(mean.t.depth$Long)[i], 3:ncol(mean.t.depth)])
  with(strat[strat$Long == unique(strat$Long)[i], ], filled.contour(Lat, WOA.depths, tmp, ylim = c(1500, 0), zlim = range(tmp, na.rm = T), color.palette = log.heat, main.title = title(paste("Bottom : Longitude =", unique(strat$Long)[i])), plot.axes = { axis(1); axis(2); points(Lat, thcl.bottom, pch = 20) }))
}

rm(i, tmp)

# Look at the spatial variation of this
with(strat, distrib.map(Long, Lat, thcl.mid, palette = "rev.log.heat"))
with(strat, distrib.map(Long, Lat, thcl.width, palette = "rev.log.heat"))
with(strat, distrib.map(Long, Lat, thcl.top, palette = "rev.log.heat"))
with(strat, distrib.map(Long, Lat, thcl.bottom, palette = "rev.log.heat"))

temp.bfd$thcl.mid <- NA
temp.bfd$thcl.width <- NA
temp.bfd$thcl.top <- NA
temp.bfd$thcl.bottom <- NA

# test this
WOA.temp(temp.bfd[runif(1, 1, 1265),], strat, strat[, grep("thcl", names(strat))])

# identify values for bfd
for (i in 1:nrow(temp.bfd)) {
  tryCatch(temp.bfd[i, grep("thcl", names(temp.bfd))] <- WOA.temp(temp.bfd[i,], strat, strat[, grep("thcl", names(strat))]), error = function(e) print(i))
}

# which have no values
which(is.na(temp.bfd$thcl.mid))
which(is.na(temp.bfd$thcl.width))
which(is.na(temp.bfd$thcl.top))
which(is.na(temp.bfd$thcl.bottom))

# create a list of those sites
na.sites <- unique(c(which(is.na(temp.bfd$thcl.mid)), which(is.na(temp.bfd$thcl.width)), which(is.na(temp.bfd$thcl.top)), which(is.na(temp.bfd$thcl.bottom))))

# add a column to the dataframe showing which are exact
temp.bfd$thercl.exact <- "Y"
temp.bfd$thercl.exact[na.sites] <- "N"

# calculate the mean values for those sites which are NAs
for(i in na.sites) {
  tryCatch(temp.bfd[i, grep("thcl", names(temp.bfd))] <- WOA.na.temp(temp.bfd[i,], strat, strat[, grep("thcl", names(strat))]), error = function(e) print(i))
}

rm(na.sites)

# look at these
summary(temp.bfd$thcl.mid)
summary(temp.bfd$thcl.width)
summary(temp.bfd$thcl.top)
summary(temp.bfd$thcl.bottom)
with(temp.bfd,distrib.map(Long, Lat, thcl.mid))
with(temp.bfd,distrib.map(Long, Lat, thcl.mid, palette = "rev.log.heat"))
with(temp.bfd,distrib.map(Long, Lat, thcl.width))
with(temp.bfd,distrib.map(Long, Lat, thcl.width, palette = "rev.log.heat"))
with(temp.bfd,distrib.map(Long, Lat, thcl.top))
with(temp.bfd,distrib.map(Long, Lat, thcl.top, palette = "rev.log.heat"))
with(temp.bfd,distrib.map(Long, Lat, thcl.bottom))
with(temp.bfd,distrib.map(Long, Lat, thcl.bottom, palette = "rev.log.heat"))

# add these rows to the bfd dataset
ldg.env$mid <- temp.bfd$thcl.mid
ldg.env$width <- temp.bfd$thcl.width
ldg.env$top <- temp.bfd$thcl.top
ldg.env$bottom <- temp.bfd$thcl.bottom
ldg.env$thcl.exact <- temp.bfd$thercl.exact

# points should be below MLD (as forced), but check, n.b. this doesn't apply to width as it is relative
with(ldg.env, plot(top, mean.pt, xlim = c(0, 1000), pch = 20, col = 1))
with(ldg.env, points(mid, mean.pt, pch = 20, col = 2))
with(ldg.env, points(bottom, mean.pt, pch = 20, col = 3))
abline(1,1)

## 4vii. plotting images for stratification --------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/Stratification")
tmp.names <- names(strat)[3:ncol(strat)]
# plot distrib.maps for all the possible strat variables
for (i in tmp.names) {
  png(file = paste(i, ".png", sep = ""), width = 800, height = 500)
  with(strat, distrib.map(Long, Lat, strat[,i], pch = 15, maintitle = i, cex = 0.4, palette = "rev.log.heat"))
  dev.off()
}

# plot contour plots at Longitude = -179.5, -33.5 for all the different strat variables
for (i in tmp.names) {
  # -179.5
  png(file = paste(i, "_cont_179.5W.png", sep = ""), width = 800, height = 500)
  # matrix for the contour background
  tmp.mat <- as.matrix(mean.t.depth[mean.t.depth$Long == -179.5, 3:ncol(mean.t.depth)])
  with(strat[strat$Long == -179.5, ], filled.contour(Lat, WOA.depths, tmp.mat, ylim = c(1500, 0), zlim = range(tmp.mat, na.rm = T), color.palette = log.heat, main.title = title(paste(i, ": Longitude =", Long[1])), xlab = "Latitude / degrees", ylab = "Depth / m", plot.axes = { axis(1); axis(2); points(Lat, strat[strat$Long == -179.5, i], pch = 20) }))
  dev.off()
  # -33.5
  png(file = paste(i, "_cont_33.5W.png", sep = ""), width = 800, height = 500)
  tmp.mat <- as.matrix(mean.t.depth[mean.t.depth$Long == -33.5, 3:ncol(mean.t.depth)])
  with(strat[strat$Long == -33.5, ], filled.contour(Lat, WOA.depths, tmp.mat, ylim = c(1500, 0), zlim = range(tmp.mat, na.rm = T), color.palette = log.heat, main.title = title(paste(i, ": Longitude =", Long[1])), xlab = "Latitude / degrees", ylab = "Depth / m", plot.axes = { axis(1); axis(2); points(Lat, strat[strat$Long == -33.5, i], pch = 20) }))
  dev.off()
}

# plot contour plots at Longitude = -179.5, -33.5 for all the different strat variables
for (i in tmp.names) {
  # -179.5
  png(file = paste(i, "_cont_179.5W_bw.png", sep = ""), width = 800, height = 500)
  # matrix for the contour background
  tmp.mat <- as.matrix(mean.t.depth[mean.t.depth$Long == -179.5, 3:ncol(mean.t.depth)])
  with(strat[strat$Long == -179.5, ], filled.contour(Lat, WOA.depths, tmp.mat, ylim = c(1500, 0), zlim = range(tmp.mat, na.rm = T), color.palette = gray.palette, main.title = title(paste(i, ": Longitude =", Long[1])), xlab = "Latitude / degrees", ylab = "Depth / m", plot.axes = { axis(1); axis(2); points(Lat, strat[strat$Long == -179.5, i], pch = 20) }))
  dev.off()
  # -33.5
  png(file = paste(i, "_cont_33.5W_bw.png", sep = ""), width = 800, height = 500)
  tmp.mat <- as.matrix(mean.t.depth[mean.t.depth$Long == -33.5, 3:ncol(mean.t.depth)])
  with(strat[strat$Long == -33.5, ], filled.contour(Lat, WOA.depths, tmp.mat, ylim = c(1500, 0), zlim = range(tmp.mat, na.rm = T), color.palette = gray.palette, main.title = title(paste(i, ": Longitude =", Long[1])), xlab = "Latitude / degrees", ylab = "Depth / m", plot.axes = { axis(1); axis(2); points(Lat, strat[strat$Long == -33.5, i], pch = 20) }))
  dev.off()
}

rm(tmp.names, tmp.mat)


## 5. Calculating Chla (4km) ------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/Chl_4km")

## 5i. Load in the data ----------------------------------------------------
chla <- read.table("monthly_chl.txt",header = T)
head(chla)

# have a look at the data
summary(chla)

## set the -32767 to NA for each month
for (i in 4:15){
  chla[which(chla[,i] == -32767),i] <- NA
}
summary(chla)

## 5ii. Save images of monthly chla ----------------------------------------
tmp <- names(chla)[4:15]
for (i in tmp) {
  png(file = paste(i, ".png", sep = ""), width = 800, height = 500)
  with(chla[!is.na(chla[,i]), ], distrib.filled(interp(LONG, LAT, chla[,i][!is.na(chla[,i]) ], duplicate = "mean"), nlevels = 100, palette = "rev.log.heat", maintitle = i))
  dev.off()
}

rm(i, tmp)

## 5iii. calculate mean and sd ---------------------------------------------
chla$meanChla <- rowMeans(chla[, 4:15], na.rm = TRUE)
chla$sdChla <- apply(chla[, 4:15], 1, sd, na.rm = TRUE)

# need na.rm = TRUE or many points are missed
summary(chla$meanChla)
summary(chla$sdChla)

# there is one site which has no data
which(is.na(chla$sdChla))
chla[262, ]
with(chla[262, ], distrib.map(LONG, LAT, 2))
# this is in the Indian Ocean - no particularly good reason for it

## 5iv. plot the mean and sd -----------------------------------------------
png(file = "meanChla.png", width = 800, height = 500)
with(chla[!is.na(chla$meanChla), ], distrib.filled(interp(LONG, LAT, meanChla, duplicate = "mean"), nlevels = 100, palette = "rev.log.heat"))
dev.off()

png(file = "sdChla.png", width = 800, height = 500)
with(chla[!is.na(chla$sdChla), ], distrib.filled(interp(LONG, LAT, sdChla, duplicate = "mean"), nlevels = 100, palette = "rev.log.heat"))
dev.off()

## 5v. add the columns to ldg.env ------------------------------------------
ldg.env$meanChla.4km <- chla$meanChla
ldg.env$sdChla.4km <- chla$sdChla


# 6. Chla for 1deg --------------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/Chl_1deg")

# 6i. Read in data for chla months ----------------------------------------
chla.1deg <- as.data.frame(cbind(rep(-179.5:179.5, 180), rep(89.5:-89.5, each = 360)))
names(chla.1deg) <- c("Long", "Lat")

for (i in 1:12) {
  if(i < 10){
    tmp <- read.csv(paste("Chla_1deg_means_0",i, ".csv", sep = ""), header = FALSE)
  } else {
    tmp <- read.csv(paste("Chla_1deg_means_",i, ".csv", sep = ""), header = FALSE)
  }
  chla.1deg <- cbind(chla.1deg, sp.mat.2.df(-179.5:179.5, 89.5:-89.5, tmp)[3])
  names(chla.1deg)[ncol(chla.1deg)] <- paste("Chla", i, sep = ".")
  print(paste("Completed", i))
}

rm(i, tmp)

head(chla.1deg)
tail(chla.1deg)
dim(chla.1deg)

# 6ii. Plot these chlorophyll data ----------------------------------------
for (i in 1:12) {
  # produce .png images in the same file  
  if (i < 10) {
    png(paste("log_chla_0", i,  ".png", sep = ""), width = 800, height = 500)
  } else {
    png(paste("log_chla_", i,  ".png", sep = ""), width = 800, height = 500)
  }
  with(chla.1deg[chla.1deg[, i + 2] < 15, ], distrib.map(Long, Lat, log(chla.1deg[chla.1deg[, i + 2] < 15, i + 2] + 1), pch = 15, maintitle = colnames(chla.1deg)[i + 2], cex = 0.4))
  dev.off()
}

# 6iii. Calculate mean and sd for these ------------------------------------
chla.1deg$meanChl <- rowMeans(chla.1deg[, 3:14], na.rm = T)
chla.1deg$sdChl <- apply(chla.1deg[, 3:14], 1, sd, na.rm = T)
# due to the spread of chl, also calculate logs
chla.1deg$mean.logChl <- rowMeans(log(chla.1deg[, 3:14] + 1), na.rm = T)
chla.1deg$sd.logChl <- apply(log(chla.1deg[, 3:14] + 1), 1, sd, na.rm = T)

# plot these
png(paste("meanChl_1deg.png", sep = ""), width = 800, height = 500)
with(chla.1deg[chla.1deg$meanChl < 10, ], distrib.map(Long, Lat, meanChl, pch = 15, cex = 0.4, palette = "rev.log.heat"))
dev.off()

summary(chla.1deg$sdChl)
png(paste("sdChl_1deg.png", sep = ""), width = 800, height = 500)
na.sites <- chla.1deg[!is.na(chla.1deg$sdChl), ]
with(na.sites[na.sites$meanChl < 10, ], distrib.map(Long, Lat, sdChl, pch = 15, cex = 0.4, palette = "rev.log.heat"))
dev.off()

rm(na.sites)

png(paste("mean_logChl_1deg.png", sep = ""), width = 800, height = 500)
with(chla.1deg[chla.1deg$meanChl < 10, ], distrib.map(Long, Lat, mean.logChl, pch = 15, cex = 0.4, palette = "rev.log.heat"))
dev.off()

png(paste("sd_logChl_1deg.png", sep = ""), width = 800, height = 500)
with(chla.1deg[chla.1deg$meanChl < 10, ], distrib.map(Long, Lat, sd.logChl, pch = 15, cex = 0.4, palette = "rev.log.heat"))
dev.off()

# 6iii. Calculate values for bfd ------------------------------------------
# generate a dataframe containing the coordinates of bfd 
chla.1deg.bfd <- ldg.data[,1:6]
# create columns
chla.1deg.bfd$meanChl <- NA
chla.1deg.bfd$sdChl <- NA
chla.1deg.bfd$mean.logChl <- NA
chla.1deg.bfd$sd.logChl <- NA

# try adding points
for(i in 1:nrow(chla.1deg.bfd))
{
  tryCatch(chla.1deg.bfd[i, 7:ncol(chla.1deg.bfd)] <- WOA.temp(chla.1deg.bfd[i,], chla.1deg, chla.1deg[, 15:ncol(chla.1deg)]), error = function(e) print(i))
}

# identify where it hasn't worked
na.sites <- which(is.na(chla.1deg.bfd$meanChl))
which(is.na(chla.1deg.bfd$sdChl))
which(is.na(chla.1deg.bfd$mean.logChl))
which(is.na(chla.1deg.bfd$sd.logChl))
# n.b. these are identical

# add a column to the dataframe showing which are exact
chla.1deg.bfd$exact <- "Y"
chla.1deg.bfd$exact[na.sites] <- "N"

# for those sites where it hasn't worked, calculate an average
for(i in na.sites)
{
  tryCatch(chla.1deg.bfd[i, 7:(ncol(chla.1deg.bfd) - 1)] <- WOA.na.temp(chla.1deg.bfd[i,], chla.1deg, chla.1deg[, 15:ncol(chla.1deg)]), error = function(e) print(i))
}

rm(na.sites)

# check whether this has worked
which(is.na(chla.1deg.bfd$meanChl)) # there now aren't any
which(is.na(chla.1deg.bfd$sdChl)) # or of these
which(is.na(chla.1deg.bfd$mean.logChl)) # there now aren't any
which(is.na(chla.1deg.bfd$sd.logChl)) # or of these

# plot the data to check it looks sensible when compared with the global plots
with(chla.1deg.bfd[chla.1deg.bfd$meanChl < 10, ], distrib.map(Long, Lat, meanChl, palette = "rev.log.heat"))
with(chla.1deg.bfd[chla.1deg.bfd$meanChl < 10, ], distrib.map(Long, Lat, sdChl, palette = "rev.log.heat"))
with(chla.1deg.bfd[chla.1deg.bfd$meanChl < 10, ], distrib.map(Long, Lat, meanChl, palette = "rev.log.heat"))
with(chla.1deg.bfd[chla.1deg.bfd$meanChl < 10, ], distrib.map(Long, Lat, meanChl, palette = "rev.log.heat"))

# 6iv. add these to ldg.env -----------------------------------------------
ldg.env$meanChla.1deg <- chla.1deg.bfd$meanChl
ldg.env$sdChla.1deg <- chla.1deg.bfd$sdChl
ldg.env$mean.logChla.1deg <- chla.1deg.bfd$mean.logChl
ldg.env$sd.logChla.1deg <- chla.1deg.bfd$sd.logChl
ldg.env$chl.exact <- chla.1deg.bfd$exact


## 7. Calculating salinity -------------------------------------------------
setwd("C:/Documents/Science/PhD/Project/BFD/Environmental/Salinity")
## Originally I tried extracting this at 0.25 degree resolution, but there was insufficient data (see Environment_test.R)

## 7i. Load in the data (1 degree resolution) --------------------------------------------
sal.mean.depth <- read.csv("s00mn01.csv")
sal.sd.depth <- read.csv("s00sd01.csv")

head(sal.mean.depth)
head(sal.sd.depth)

## 7ii. Save images of salinity --------------------------------------------

# generate plots for mean (n.b. not worth making a function as I want to change parameters)
for (i in 1:33) {
  # produce .png images in the same file  
  data.na <- sal.mean.depth[!is.na(sal.mean.depth[, i + 2]),  ]
  png(paste("sal_mn_", i,  ".png", sep = ""), width = 800, height = 500)
  with(data.na[data.na[, i + 2] > 20, ], distrib.map(Longitude, Latitude, data.na[data.na[, i + 2] > 20, i + 2], pch = 15, maintitle = colnames(data)[i + 2], cex = 0.4))
  dev.off()
}

rm(data.na)

# repeat for sd
for (i in 1:33) {
  # produce .png images in the same file  
  data.na <- sal.sd.depth[!is.na(sal.sd.depth[, i + 2]),  ]
  png(paste("sal_sd_", i,  ".png", sep = ""), width = 800, height = 500)
  with(data.na[data.na[, i + 2] < 3.5, ], distrib.map(Longitude, Latitude, data.na[data.na[, i + 2] < 3.5, i + 2], pch = 15, maintitle = colnames(data)[i + 2], cex = 0.4))
  dev.off()
}

rm(data.na)

## 7iii. calculate the bfd sites ------------------------------------------
# only add the surface salinity at the moment - not clear whether it should only be surface salinity

# generate a dataframe containing the coordinates to calculate temperature at those points. 
sal.bfd <- ldg.data[,1:6]
# create columns for mean salinity
for (i in 3:ncol(sal.mean.depth)){
  sal.bfd <- cbind(sal.bfd,NA)
  colnames(sal.bfd)[4 + i] <- as.character(paste("mean", colnames(sal.mean.depth[i]), sep = "."))
}
# and for sd
for (i in 3:ncol(sal.sd.depth)){
  sal.bfd <- cbind(sal.bfd,NA)
  colnames(sal.bfd)[106 + i] <- as.character(paste("sd", colnames(sal.sd.depth[i]), sep = "."))
}

# add the salinity data (mean and sd) for the whole of the dataset
for(i in 1:nrow(sal.bfd))
{
  tryCatch(sal.bfd[i, grep("mean", colnames(sal.bfd))] <- WOA.temp(sal.bfd[i,], sal.mean.depth, sal.mean.depth[, 3:ncol(sal.mean.depth)]), error = function(e) print(i))
  
  tryCatch(sal.bfd[i, grep("sd", colnames(sal.bfd))] <- WOA.temp(sal.bfd[i,], sal.sd.depth, sal.sd.depth[, 3:ncol(sal.sd.depth)]), error = function(e) print(i))
}

# identify where it hasn't worked (only really interested in surface atm)
which(is.na(sal.bfd$mean.Depth0m))
sum(is.na(sal.bfd$mean.Depth0m) != is.na(sal.bfd$sd.Depth0m)) # these aren't identical

# create a list of those sites where it hasn't worked
na.sites <- sort(unique(c(which(is.na(sal.bfd$mean.Depth0m)), which(is.na(sal.bfd$sd.Depth0m)))))

# add a column to the dataframe showing which are exact
sal.bfd$exact <- "Y"
sal.bfd$exact[na.sites] <- "N"

# for those sites where it hasn't worked, calculate an average
for(i in na.sites)
{
  tryCatch(sal.bfd[i, grep("mean", colnames(sal.bfd))] <- WOA.na.temp(sal.bfd[i,], sal.mean.depth, sal.mean.depth[, 3:ncol(sal.mean.depth)]), error = function(e) print(i))
  
  tryCatch(sal.bfd[i, grep("sd", colnames(sal.bfd))] <- WOA.na.temp(sal.bfd[i,], sal.sd.depth, sal.sd.depth[, 3:ncol(sal.sd.depth)]), error = function(e) print(i))
}

rm(na.sites)

# check whether this has worked
which(is.na(sal.bfd$mean.Depth0m)) # there now aren't any
which(is.na(sal.bfd$sd.Depth0m)) # or of these

# plot the data to check it looks sensible when compared with the global plots
with(sal.bfd[sal.bfd$mean.Depth0m > 20, ], distrib.map(Long, Lat, mean.Depth0m))
with(sal.bfd[sal.bfd$mean.Depth0m > 20 & sal.bfd$exact == "N", ], distrib.map(Long, Lat, mean.Depth0m))

## 7iv. Add the surface salinity data to ldg.env --------------------------
ldg.env$meanSal.0m <- sal.bfd$mean.Depth0m
ldg.env$sdSal.0m <- sal.bfd$sd.Depth0m
ldg.env$sal.exact <- sal.bfd$exact


## 8. save the data --------------------------------------------------------
setwd("C:/Documents/Science/PhD/Work/1311 LDGPaper/")
# save complete workspace
save.image(file = "Output/140522_Environment_ws.Rdata")
# save ldg.env
save(ldg.env, file = "Output/140522_ldg_env.Rdata")
# save a file with both environmental and species data
comb <- merge(ldg.env, ldg.data)
write.csv(comb, "data_env.csv", row.names = FALSE)
