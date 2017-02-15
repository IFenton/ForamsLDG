## Created: 18 / 8 / 2016
## Last edited: 21 / 8 / 2016
## Isabel Fenton
## Reanalysis for LDG paper response to reviewers
##
## Previous file: MARGO/Code/Environmental_variables.R
## Next file: 1311 LDGPaper/Reanalysis/Prediction_environment.R

## Inputs ------------------------------------------------------------------
## Environmental_variables.Rdata - containing ldg.margo.data, ldg.margo.env and other files

## Outputs ----------------------------------------------------------------
## 150601 ldg_margo_dup.RData - the dataset with duplicates (on a 1degree scale) included
## ldg_margo_mod.RData - The modelling dataframe
## Images saved into Figures with prefix Ana
## Atlantic_simplification.RData - models enroute to simplification
## Indian_simplification.RData - models enroute to simplification
## Pacific_simplification.RData - models enroute to simplification
## Atlantic_simplified.RData - Simplified model
## Indian_simplified.RData - Simplified model
## Pacific_simplified.RData - simplified model
## mod_hres.RData - models for comparing data resolution
## Evenness_coding.RData - testing different coding styles for evenness
## Lineage_coding.RData - testing different coding styles for lineage age
## Metabolic_hypothesis.RData - testing the metabolic hypothesis
## Richness_model.RData - the final (complete) model for richness
## Richness_model_simplified.RData - the final simplified model for richness
## Evenness_model.RData - the final (complete) model for evenness
## Lineage_model.RData - the final (complete) model for lineage age
## op0_summary.csv - table of full model cofficients for richness
## eve0_summary.csv - table of full model cofficients for evenness
## lna0_summary.csv - table of full model cofficients for lineage age
## lr_out.csv - likelihood ratios for the richness model

## libraries ---------------------------------------------------------------
setwd("C:/Users/isabf/Dropbox/Documents/PhD/Work/1311 LDGPaper/Reanalysis2")
library(spdep) # for SAR models
library(HH) # for vif
library(ncf) # for Moran's I
library(lmtest) # for likelihood ratio tests
library(mgcv) # for gams
library(colorRamps) # for matlab.like palette
library(scales) # for alpha
source("../Code/140420SARerrOptimising_NC.R") # for the optimising function
source("../../../Code/plot_spline_correlog_n.R") # for plotting spline correlog (with axis titles etc.)
source("../../../Code/sar_plot.R") # for producing sar plots
source("../../../Code/sar_predict.R") # for predicting with sar
source("../../../Code/maps.R") # for maps
source("../../../Code/palettes.R") # ditto
source("../../../Code/lr_calculations.R") # code for calculating likelihood ratios
source("../../../Code/sar_predict.R") # for predicting with poly
source("../../../Code/compare.R") # for checking species names
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling


## 0i. Setting up the datasets -----------------------------------------------
rm(db.traits, db.traits.cons, WOA.depths)

# create a dataset for modelling with 
tmp <- c("Core", "Latitude", "Longitude", "Water.Depth", "Total_Planktics", "Ocean2", "sp.rich", "rarefy.sr", "simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
ldg.margo.mod <- merge(ldg.margo.env, ldg.margo.data[, tmp], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env, tmp)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]

# remove the extra factor level of the mediterranean
table(ldg.margo.mod$Ocean2)
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
table(ldg.margo.mod$Ocean2)

## which points are to be excluded because of delta_carb_ion?
with(ldg.margo.mod, plot(delta_carb_ion, rarefy.sr, pch = 16, col = Ocean2))
legend("topleft", levels(ldg.margo.mod$Ocean2), pch = 16, col = 1:3)
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -10.908), ]
with(ldg.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(ldg.margo.mod$Ocean2)
# based on using a 3500m / 4500m cut-off

dim(ldg.margo.mod)

# for dissolution, only want to use it to account for dissolution. As there is no dissolution at sites with delta_carb_ion > 0, set all these to zero. It now becomes a measure of delta_carb_ion below zero. 
ldg.margo.mod$delta_carb_ion[ldg.margo.mod$delta_carb_ion > 0] <- 0
with(ldg.margo.mod, plot(delta_carb_ion, rarefy.sr, pch = 16, col = Ocean2))

# expect salinity to impact both at the top and the bottom of the range. 
with(ldg.margo.mod, plot(meanSal.0m, rarefy.sr, pch = 16, col = Ocean2))
# n.b. much of the pattern in the Atlantic is driven by SST, so although we see a relationship it isn't what is seems. 
# work out ocean average
load("../../../Project/BFD/Environmental/Salinity/150414_salinity.RData")
summary(sal.mean.depth$Depth0m)
# ocean average is 34.28, but if we look at PF optima (see Be & Hutson, 1977), find the average optimal salinity of 35.1. Given we care about how the salinity is differing from optimum for species this is a more logical number to use.
rm(sal.margo, sal.mean.depth, sal.sd.depth)
ldg.margo.mod$absMnSal.0m <- abs(ldg.margo.mod$meanSal.0m - 35.1)
with(ldg.margo.mod, plot(absMnSal.0m, rarefy.sr, pch = 16, col = Ocean2))

## remove sites which have the wrong taxonomy
# check the ranges of the environmental variables with and without these sites
summary(ldg.margo.mod)
summary(ldg.margo.mod[ldg.margo.mod$cons_tax == "Y", ])
# they are the same, so hopefully won't make any difference to the modelling results

ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$cons_tax == "Y", ]

## save this dataset
save(ldg.margo.mod, file = "Outputs/ldg_margo_mod2.RData")

## create three different datasets for Rarefied, Evenness & LineageAge, & FRic
## for richness
cols <- c("simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
rsr.margo.mod <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
# remove other NAs
summary(rsr.margo.mod)
rsr.margo.mod <- na.omit(rsr.margo.mod)

## for evenness / lineage age
cols <- c("sp.rich", "rarefy.sr", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "cons_tax")
eve.margo.mod <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
# remove other NAs
summary(eve.margo.mod)
eve.margo.mod <- na.omit(eve.margo.mod)

## for FRic
cols <- c("sp.rich", "rarefy.sr", "simpson", "simpsonEve", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
fric.margo.mod <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
# remove other NAs
summary(fric.margo.mod)
fric.margo.mod <- na.omit(fric.margo.mod)

# Consider dimensions
dim(rsr.margo.mod)
dim(eve.margo.mod)
dim(fric.margo.mod)

## 0ii. How to handle multiple points in a grid cell --------------------------
# for richness
head(rsr.margo.mod)
# n.b. don't want to include the finer resolution SST, and productivity is at 1/6 degree, so exclude that as well
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK", "mean.bvf", "max.bvf", "depth.bvf")
dim(unique(rsr.margo.mod[, cols]))
tmp.1 <- which(duplicated(rsr.margo.mod[, cols]))

rsr.margo.mod$uni <- NA

# add a column for each unique set
rsr.margo.mod$uni[!duplicated(rsr.margo.mod[, cols])] <- 1:length(rsr.margo.mod$uni[!duplicated(rsr.margo.mod[, cols])])

## the sort through the duplicated rows to match with the unique sets
# for each row in the data
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(rsr.margo.mod[i, ], rsr.margo.mod, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  rsr.margo.mod$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows)
# make uni a factor as basically each value represents a unique grid cell
rsr.margo.mod$uni <- factor(rsr.margo.mod$uni)

rsr.margo.dup <- rsr.margo.mod

# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
# create a dataframe of the right length
rsr.margo.mod <- rsr.margo.dup[!duplicated(rsr.margo.dup$uni), ]

# for all the relevant columns calculate means, and replace these in the dataset
for (i in 1:ncol(rsr.margo.mod)) {
  if (!is.factor(rsr.margo.mod[,i]) & !is.character(rsr.margo.mod[, i])) {
    rsr.margo.mod[, i] <- as.numeric(tapply(rsr.margo.dup[, i], rsr.margo.dup$uni, mean, na.rm = TRUE))
  }
}
rm(i)
save(rsr.margo.dup, file = "Outputs/160821 rsr_margo_dup.RData")
rm(rsr.margo.dup, cols, tmp.1)


# for evenness / lineage age
head(eve.margo.mod)
# n.b. don't want to include the finer resolution SST, and productivity is at 1/6 degree, so exclude that as well
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK", "mean.bvf", "max.bvf", "depth.bvf")
dim(unique(eve.margo.mod[, cols]))
tmp.1 <- which(duplicated(eve.margo.mod[, cols]))

eve.margo.mod$uni <- NA

# add a column for each unique set
eve.margo.mod$uni[!duplicated(eve.margo.mod[, cols])] <- 1:length(eve.margo.mod$uni[!duplicated(eve.margo.mod[, cols])])

## the sort through the duplicated rows to match with the unique sets
# for each row in the data
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(eve.margo.mod[i, ], eve.margo.mod, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  eve.margo.mod$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows)
# make uni a factor as basically each value represents a unique grid cell
eve.margo.mod$uni <- factor(eve.margo.mod$uni)

eve.margo.dup <- eve.margo.mod

# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
# create a dataframe of the right length
eve.margo.mod <- eve.margo.dup[!duplicated(eve.margo.dup$uni), ]

# for all the relevant columns calculate means, and replace these in the dataset
for (i in 1:ncol(eve.margo.mod)) {
  if (!is.factor(eve.margo.mod[,i]) & !is.character(eve.margo.mod[, i])) {
    eve.margo.mod[, i] <- as.numeric(tapply(eve.margo.dup[, i], eve.margo.dup$uni, mean, na.rm = TRUE))
  }
}
rm(i)
save(eve.margo.dup, file = "Outputs/160821 eve_margo_dup.RData")
rm(eve.margo.dup, cols, tmp.1)


# for FRic
head(fric.margo.mod)
# n.b. don't want to include the finer resolution SST, and productivity is at 1/6 degree, so exclude that as well
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK", "mean.bvf", "max.bvf", "depth.bvf")
dim(unique(fric.margo.mod[, cols]))
tmp.1 <- which(duplicated(fric.margo.mod[, cols]))

fric.margo.mod$uni <- NA

# add a column for each unique set
fric.margo.mod$uni[!duplicated(fric.margo.mod[, cols])] <- 1:length(fric.margo.mod$uni[!duplicated(fric.margo.mod[, cols])])

## the sort through the duplicated rows to match with the unique sets
# for each row in the data
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(fric.margo.mod[i, ], fric.margo.mod, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  fric.margo.mod$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows)
# make uni a factor as basically each value represents a unique grid cell
fric.margo.mod$uni <- factor(fric.margo.mod$uni)

fric.margo.dup <- fric.margo.mod

# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
# create a dataframe of the right length
fric.margo.mod <- fric.margo.dup[!duplicated(fric.margo.dup$uni), ]

# for all the relevant columns calculate means, and replace these in the dataset
for (i in 1:ncol(fric.margo.mod)) {
  if (!is.factor(fric.margo.mod[,i]) & !is.character(fric.margo.mod[, i])) {
    fric.margo.mod[, i] <- as.numeric(tapply(fric.margo.dup[, i], fric.margo.dup$uni, mean, na.rm = TRUE))
  }
}
rm(i)
save(fric.margo.dup, file = "Outputs/160821 fric_margo_dup.RData")
rm(fric.margo.dup, cols, tmp.1)

## 0iii. Consider relationships ---------------------------------------------
# use most complete dataset i.e. ldg.margo.mod
par(ask = TRUE)
for(i in c(5:34, ncol(ldg.margo.mod))) {
  if (!is.character(ldg.margo.mod[, i]))
    with(ldg.margo.mod, plot(ldg.margo.mod[, i], rarefy.sr, pch = 16, col = Ocean2, main = names(ldg.margo.mod)[i]))
}
par(ask = FALSE)
rm(i)

# check whether anything else should be logged (i.e. the histogram of each environmental variable)
summary(ldg.margo.mod)
par(ask = TRUE)
for(i in c(5:34, ncol(ldg.margo.mod))) {
  if (!is.character(ldg.margo.mod[, i]))
    with(ldg.margo.mod, hist(ldg.margo.mod[, i], main = names(ldg.margo.mod)[i]))
}
par(ask = FALSE)
rm(i)
# seem reasonable

## 0iv. Create plots ------------------------------------------------------
png("Figures/Ana_0iii_map_rsr_2.png", 700, 500)
with(rsr.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", main = "Rarefied species richness", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_eve_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", main = "Simpson's Evenness", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_lna_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like", main = "Average Community Age", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_fric_2.png", 700, 500)
with(fric.margo.mod, distrib.map(Longitude, Latitude, FRic, palette = "matlab.like", main = "Functional richness", col.water = "white", col.land = "black"))
dev.off()

# tidy up
rm(ldg.margo.mod)


## 1. Check for correlation between explanatory variables ----------------------
names(rsr.margo.mod)

# variables are: mean/sd SST, mean/sd MLD, 10deg contour, mean/sd logProd, mean/sd Sal, Ocean2, carbonate ion 
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "sd.mld.t", "depth10deg", "logProd.mn.ann", "logProd.sd.ann", "absMnSal.0m", "sdSal.0m", "Ocean2", "delta_carb_ion", "meanOxy", "prop2.oxy", "mean.bvf", "max.bvf")

# pairs plot
png("Figures/Ana_1_pairs_2.png", 1200, 1200)
pairs(rsr.margo.mod[, names(rsr.margo.mod) %in% env.var])
dev.off()

# variance inflation factor
vif(rsr.margo.mod[, names(rsr.margo.mod) %in% env.var])

# currently the only ones that are potentially correlated (VIF > 5 and look correlated on the pairs plot) are mean and sd in Prod, and mean / sd mld and mean / max bvf, also meanOxy and meanSST.1deg seem to be correlated (as suggested). Look into this in a bit more detail
png("Figures/Ana_1_mnprodsdprod1deg_2.png")
with(rsr.margo.mod, plot(logProd.mn.ann, logProd.sd.ann, pch = 16)) # highly correlated.
dev.off()

# does the same hold for for the mld
png("Figures/Ana_1_mnMLDsdMLD_2.png")
with(rsr.margo.mod, plot(mean.mld.t, sd.mld.t, pch = 16)) # yes
dev.off()

png("Figures/Ana_1_mnBVFmxBVF_2.png")
with(rsr.margo.mod, plot(mean.bvf, max.bvf, pch = 16)) # yes
dev.off()

png("Figures/Ana_1_mnSST1degmnOxy_2.png")
with(rsr.margo.mod, plot(meanSST.1deg, meanOxy, pch = 16)) # highly correlated.
dev.off()

# therefore suggest exclusion of logProd.sd.ann, sd.mld.t, max.bvf and meanOxy from models, so
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "Ocean2", "delta_carb_ion", "prop2.oxy")

# check this has fixed the problem
pairs(rsr.margo.mod[, names(rsr.margo.mod) %in% env.var] )
vif(rsr.margo.mod[, names(rsr.margo.mod) %in% env.var] )
# it has, so now have new EVs


## 2. Model simplification -------------------------------------------------

## 2i. Create an OLS model and check for SAC --------------------------------
# an OLS model
mod.l0 <- lm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, data = rsr.margo.mod)

# check model plots
png("Figures/Ana_2i_modl0_2.png", 600, 600)
par(mfrow = c(2, 2))
plot(mod.l0)
par(mfrow = c(1, 1))
dev.off()

# look for spatial autocorrelation in the residuals
# using spline.correlog
# haven't run this properly - should be resamp = 1000
mod.l0.sac <- with(rsr.margo.mod, spline.correlog(Longitude, Latitude, mod.l0$residuals, latlon = TRUE, resamp = 1))
summary(mod.l0.sac)
png("Figures/Ana_2i_modl0SAC_2.png")
plot.spline.correlog.n(mod.l0.sac, xlab = "Distance / km")
dev.off()

# using correlog
mod.l0.SACcor <- with(rsr.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.l0), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2i_modl0SACcor_2.png")
plot(mod.l0.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()
rm(mod.l0.SACcor)

## look at the residuals plots
png("Figures/Ana_2i_modl0resid_2.png", 700, 500)
with(rsr.margo.mod, distrib.map(Longitude, Latitude, mod.l0$residuals, palette = "rwb"))
dev.off()

rm(mod.l0)

## 2ii. Create a GAM to check complexity / SAC -----------------------------
# n.b. GAMs can't do interactions (as additive models)
mod.g0 <- with(rsr.margo.mod, gam(rarefy.sr ~ s(Longitude, Latitude, k = 80, by = Ocean2) + s(meanSST.1deg, by = Ocean2) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2 + delta_carb_ion, gamma = 1.4))
summary(mod.g0)

png("Figures/Ana_2ii_modg0_2.png")
gam.check(mod.g0) # n.b. this gives an error for factors
dev.off()
par(mfrow = c(1,1))

## calculate SAC
# using spline.correlog
mod.g0.SAC <- with(rsr.margo.mod, spline.correlog(Longitude, Latitude, mod.g0$residuals, latlon = TRUE, resamp = 1))
summary(mod.g0.SAC)
png("Figures/Ana_2ii_modg0SAC_2.png")
plot.spline.correlog.n(mod.g0.SAC, xlab = "Distance / km")
dev.off()

# using correlog
mod.g0.SACcor <- with(rsr.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.g0), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2ii_modg0SACcor_2.png")
plot(mod.g0.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()

rm(mod.g0, mod.g0.SAC, mod.g0.SACcor)

## 2iii. Create an optimised SARerror model --------------------------------
# Make a matrix of coordinates (X and Y coordinates)
ldg.coords <- cbind(rsr.margo.mod$Longitude, rsr.margo.mod$Latitude)
ldg.coords <- as.matrix(ldg.coords)

# run model optimisation
# getting problems with Error in solve.default(asyvar, tol = tol.solve) : 
# system is computationally singular: reciprocal condition number = 8.20242e-19 
# The suggested answer is to rescale if the variables are on very different scales, so go from
# mod.sar.opW <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
# to
summary(rsr.margo.mod) # check that ranges are roughly equivalent after scaling
mod.sar.opW <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
# reckon that delta_carb_ion should not interact with things as its relationship with species richness (at least the bit we care about) shouldn't depend on anything else
summary(mod.sar.opW$obj, Nagelkerke = TRUE)

## check SAC has been removed
# using spline.correlog
mod.sar.opW.SAC <- with(rsr.margo.mod, spline.correlog(Longitude, Latitude, mod.sar.opW$obj$residuals, latlon = TRUE, resamp = 1))
summary(mod.sar.opW.SAC)
png("Figures/Ana_2iii_modSarOp0SAC_2.png")
plot.spline.correlog.n(mod.sar.opW.SAC, xlab = "Distance / km")
dev.off()
rm(mod.sar.opW.SAC)

# using correlog
mod.sar.opW.SACcor <- with(rsr.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.sar.opW$obj), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2iii_modSarOp0SACcor_2.png")
plot(mod.sar.opW.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()
rm(mod.sar.opW.SACcor)

# check whether different coding methods improve the AIC
mod.sar.opB <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 4253.397
AIC(mod.sar.opB$obj) # 4286.374

mod.sar.opS <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 4253.397
AIC(mod.sar.opS$obj) # 4259.442

mod.sar.opW
# So "W" is the best coding style and the best neighbourhood distance is 507.1066
rm(mod.sar.opS, mod.sar.opB)

## 2iv. Run model simplification -------------------------------------------
summary(mod.sar.opW$obj)

# re-run this optimised model through errorsarlm, so things like anova and lr.calc work
op.nb <- dnearneigh(ldg.coords, 0, mod.sar.opW$dist, longlat = TRUE)
op.w <- nb2listw(op.nb, glist = NULL, style = "W", zero.policy = TRUE)
mod.sar.op0 <- errorsarlm(mod.sar.opW$mod, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18)
rm(op.nb)

# start model simplification
summary(mod.sar.op0, Nagelkerke = TRUE) # 0.89216
summary(mod.sar.op0)$Coef[order(summary(mod.sar.op0)$Coef[, 4]),]
mod.sar.op1 <- update(mod.sar.op0, ~. -I(depth10deg/100):logProd.mn.ann)
anova(mod.sar.op0, mod.sar.op1)
summary(mod.sar.op1, Nagelkerke = TRUE) # 0.89216
AIC(mod.sar.op1) # 4251.402
summary(mod.sar.op1)$Coef[order(summary(mod.sar.op1)$Coef[, 4]),]

mod.sar.op2 <- update(mod.sar.op1, ~. -poly(meanSST.1deg, 3):sdSST.1deg + poly(meanSST.1deg, 2):sdSST.1deg)
anova(mod.sar.op2, mod.sar.op1)
summary(mod.sar.op2, Nagelkerke = TRUE) # 0.89216 
AIC(mod.sar.op2) # 4249.422
rm(mod.sar.op1)
summary(mod.sar.op2)$Coef[order(summary(mod.sar.op2)$Coef[, 4]),]

mod.sar.op3 <- update(mod.sar.op2, ~. -logProd.mn.ann:prop2.oxy )
anova(mod.sar.op3, mod.sar.op2)
summary(mod.sar.op3, Nagelkerke = TRUE) # 0.89215 
AIC(mod.sar.op3) # 4247.492
rm(mod.sar.op2)
summary(mod.sar.op3)$Coef[order(summary(mod.sar.op3)$Coef[, 4]),]

mod.sar.op4 <- update(mod.sar.op3, ~. -sdSST.1deg:poly(meanSST.1deg, 2) + sdSST.1deg:poly(meanSST.1deg, 1))
anova(mod.sar.op4, mod.sar.op3)
summary(mod.sar.op4, Nagelkerke = TRUE) # 0.89214
AIC(mod.sar.op4) # 4245.561
rm(mod.sar.op3)
summary(mod.sar.op4)$Coef[order(summary(mod.sar.op4)$Coef[, 4]),]

mod.sar.op5 <- update(mod.sar.op4, ~. -poly(meanSST.1deg, 3):I(mean.mld.t/10) + poly(meanSST.1deg, 2):I(mean.mld.t/10))
anova(mod.sar.op5, mod.sar.op4)
summary(mod.sar.op5, Nagelkerke = TRUE) # 0.89213
AIC(mod.sar.op5) # 4243.659
rm(mod.sar.op4)
summary(mod.sar.op5)$Coef[order(summary(mod.sar.op5)$Coef[, 4]),]

mod.sar.op6 <- update(mod.sar.op5, ~. -sdSST.1deg:absMnSal.0m)
anova(mod.sar.op6, mod.sar.op5)
summary(mod.sar.op6, Nagelkerke = TRUE) # 0.89212
AIC(mod.sar.op6) # 4241.756
rm(mod.sar.op5)
summary(mod.sar.op6)$Coef[order(summary(mod.sar.op6)$Coef[, 4]),]

mod.sar.op7 <- update(mod.sar.op6, ~. -poly(meanSST.1deg, 3):I(depth10deg/100) + poly(meanSST.1deg, 2):I(depth10deg/100))
anova(mod.sar.op7, mod.sar.op6)
summary(mod.sar.op7, Nagelkerke = TRUE) # 0.89211
AIC(mod.sar.op7) # 4239.864
rm(mod.sar.op6)
summary(mod.sar.op7)$Coef[order(summary(mod.sar.op7)$Coef[, 4]),]

mod.sar.op8 <- update(mod.sar.op7, ~. -I(mean.mld.t/10):prop2.oxy )
anova(mod.sar.op8, mod.sar.op7)
summary(mod.sar.op8, Nagelkerke = TRUE) # 0.89209
AIC(mod.sar.op8) # 4238.092
rm(mod.sar.op7)
summary(mod.sar.op8)$Coef[order(summary(mod.sar.op8)$Coef[, 4]),]

mod.sar.op9 <- update(mod.sar.op8, ~. -sdSST.1deg:I(depth10deg/100) )
anova(mod.sar.op9, mod.sar.op8)
summary(mod.sar.op9, Nagelkerke = TRUE) # 0.89206
AIC(mod.sar.op9) # 4236.386
rm(mod.sar.op8)
summary(mod.sar.op9)$Coef[order(summary(mod.sar.op9)$Coef[, 4]),]

mod.sar.op10 <- update(mod.sar.op9, ~. -sdSST.1deg:logProd.mn.ann )
anova(mod.sar.op10, mod.sar.op9)
summary(mod.sar.op10, Nagelkerke = TRUE) # 0.892
AIC(mod.sar.op10) # 4234.991
rm(mod.sar.op9)
summary(mod.sar.op10)$Coef[order(summary(mod.sar.op10)$Coef[, 4]),]

mod.sar.op11 <- update(mod.sar.op10, ~. -sdSal.0m:Ocean2)
anova(mod.sar.op11, mod.sar.op10)
summary(mod.sar.op11, Nagelkerke = TRUE) # 0.89194
AIC(mod.sar.op11) # 4231.56
rm(mod.sar.op10)
summary(mod.sar.op11)$Coef[order(summary(mod.sar.op11)$Coef[, 4]),]

mod.sar.op12 <- update(mod.sar.op11, ~. -absMnSal.0m:sdSal.0m )
anova(mod.sar.op12, mod.sar.op11)
summary(mod.sar.op12, Nagelkerke = TRUE) # 0.89189
AIC(mod.sar.op12) # 4230.036
rm(mod.sar.op11)
summary(mod.sar.op12)$Coef[order(summary(mod.sar.op12)$Coef[, 4]),]

mod.sar.op13 <- update(mod.sar.op12, ~. -I(mean.mld.t/10):sdSal.0m )
anova(mod.sar.op13, mod.sar.op12)
summary(mod.sar.op13, Nagelkerke = TRUE) # 0.89184
AIC(mod.sar.op13) #  4228.506
rm(mod.sar.op12)
summary(mod.sar.op13)$Coef[order(summary(mod.sar.op13)$Coef[, 4]),]

mod.sar.op14 <- update(mod.sar.op13, ~. -I(depth10deg/100):Ocean2)
anova(mod.sar.op14, mod.sar.op13)
summary(mod.sar.op14, Nagelkerke = TRUE) # 0.89176
AIC(mod.sar.op14) # 4225.307
rm(mod.sar.op13)
summary(mod.sar.op14)$Coef[order(summary(mod.sar.op14)$Coef[, 4]),]

mod.sar.op15 <- update(mod.sar.op14, ~. -I(mean.mld.t/10):I(depth10deg/100) )
anova(mod.sar.op15, mod.sar.op14)
summary(mod.sar.op15, Nagelkerke = TRUE) # 0.89168
AIC(mod.sar.op15) # 4224.159
rm(mod.sar.op14)
summary(mod.sar.op15)$Coef[order(summary(mod.sar.op15)$Coef[, 4]),]

mod.sar.op16 <- update(mod.sar.op15, ~. -poly(meanSST.1deg, 3):absMnSal.0m + poly(meanSST.1deg, 2):absMnSal.0m)
anova(mod.sar.op16, mod.sar.op15)
summary(mod.sar.op16, Nagelkerke = TRUE) # 0.89157
AIC(mod.sar.op16) # 4223.157
rm(mod.sar.op15)
summary(mod.sar.op16)$Coef[order(summary(mod.sar.op16)$Coef[, 4]),]

mod.sar.op17 <- update(mod.sar.op16, ~. -sdSST.1deg:prop2.oxy  )
anova(mod.sar.op17, mod.sar.op16)
summary(mod.sar.op17, Nagelkerke = TRUE) # 0.89145
AIC(mod.sar.op17) # 4222.382
rm(mod.sar.op16)
summary(mod.sar.op17)$Coef[order(summary(mod.sar.op17)$Coef[, 4]),]

mod.sar.op18 <- update(mod.sar.op17, ~. -I(depth10deg/100):poly(meanSST.1deg, 2)  + I(depth10deg/100):poly(meanSST.1deg, 1))
anova(mod.sar.op18, mod.sar.op17)
summary(mod.sar.op18, Nagelkerke = TRUE) # 0.89128
AIC(mod.sar.op18) # 4222.049
rm(mod.sar.op17)
summary(mod.sar.op18)$Coef[order(summary(mod.sar.op18)$Coef[, 4]),]

mod.sar.op19 <- update(mod.sar.op18, ~. -I(depth10deg/100):poly(meanSST.1deg, 1))
anova(mod.sar.op19, mod.sar.op18)
summary(mod.sar.op19, Nagelkerke = TRUE) # 0.89124
AIC(mod.sar.op19) # 4220.432
rm(mod.sar.op18)
summary(mod.sar.op19)$Coef[order(summary(mod.sar.op19)$Coef[, 4]),]

mod.sar.op20 <- update(mod.sar.op19, ~. -I(depth10deg/100):prop2.oxy)
anova(mod.sar.op20, mod.sar.op19)
summary(mod.sar.op20, Nagelkerke = TRUE) # 0.89111
AIC(mod.sar.op20) # 4219.697
rm(mod.sar.op19)
summary(mod.sar.op20)$Coef[order(summary(mod.sar.op20)$Coef[, 4]),]

mod.sar.op21 <- update(mod.sar.op20, ~. -poly(meanSST.1deg, 3):prop2.oxy + poly(meanSST.1deg, 2):prop2.oxy)
anova(mod.sar.op21, mod.sar.op20)
summary(mod.sar.op21, Nagelkerke = TRUE) # 0.89095
AIC(mod.sar.op21) # 4219.224
rm(mod.sar.op20)
summary(mod.sar.op21)$Coef[order(summary(mod.sar.op21)$Coef[, 4]),]

mod.sar.op22 <- update(mod.sar.op21, ~. -poly(meanSST.1deg, 3):logProd.mn.ann + poly(meanSST.1deg, 2):logProd.mn.ann)
anova(mod.sar.op22, mod.sar.op21)
summary(mod.sar.op22, Nagelkerke = TRUE) # 0.89075 
AIC(mod.sar.op22) # 4219.205
rm(mod.sar.op21)
summary(mod.sar.op22)$Coef[order(summary(mod.sar.op22)$Coef[, 4]),]

mod.sar.op23 <- update(mod.sar.op22, ~. -logProd.mn.ann:poly(meanSST.1deg, 2) + logProd.mn.ann:poly(meanSST.1deg, 1))
anova(mod.sar.op23, mod.sar.op22)
summary(mod.sar.op23, Nagelkerke = TRUE) # 0.89068
AIC(mod.sar.op23) # 4217.853  
rm(mod.sar.op22)
summary(mod.sar.op23)$Coef[order(summary(mod.sar.op23)$Coef[, 4]),]

mod.sar.op24 <- update(mod.sar.op23, ~. -I(mean.mld.t/10):poly(meanSST.1deg, 2) + I(mean.mld.t/10):poly(meanSST.1deg, 1))
anova(mod.sar.op24, mod.sar.op23)
summary(mod.sar.op24, Nagelkerke = TRUE) # 0.89047
AIC(mod.sar.op24) # 4217.919
rm(mod.sar.op23)
summary(mod.sar.op24)$Coef[order(summary(mod.sar.op24)$Coef[, 4]),]

mod.sar.op25 <- update(mod.sar.op24, ~. -I(mean.mld.t/10):poly(meanSST.1deg, 1))
anova(mod.sar.op25, mod.sar.op24)
summary(mod.sar.op25, Nagelkerke = TRUE) # 0.89047
AIC(mod.sar.op25) # 4215.922
rm(mod.sar.op24)
summary(mod.sar.op25)$Coef[order(summary(mod.sar.op25)$Coef[, 4]),]

mod.sar.op26 <- update(mod.sar.op25, ~. -sdSST.1deg:Ocean2)
anova(mod.sar.op26, mod.sar.op25)
summary(mod.sar.op26, Nagelkerke = TRUE) # 0.89024
AIC(mod.sar.op26) # 4214.076
rm(mod.sar.op25)
summary(mod.sar.op26)$Coef[order(summary(mod.sar.op26)$Coef[, 4]),]

# having got to here (which is that all values are <0.05 or can't be removed), then check whether any others should be removed (particularly those with ocean)
mod.sar.op27 <- update(mod.sar.op26, ~. -I(mean.mld.t/10):absMnSal.0m )
anova(mod.sar.op27, mod.sar.op26)
summary(mod.sar.op27, Nagelkerke = TRUE) # 0.88987
AIC(mod.sar.op27) # 4215.684
rm(mod.sar.op26)
summary(mod.sar.op27)$Coef[order(summary(mod.sar.op27)$Coef[, 4]),]

mod.sar.op28 <- update(mod.sar.op27, ~. -logProd.mn.ann:absMnSal.0m )
anova(mod.sar.op28, mod.sar.op27)
summary(mod.sar.op28, Nagelkerke = TRUE) # 0.88968
AIC(mod.sar.op28) # 4215.515
rm(mod.sar.op27)
summary(mod.sar.op28)$Coef[order(summary(mod.sar.op28)$Coef[, 4]),]

# having got to here (which is that all values are <0.05 or can't be removed), then check whether any others should be removed (particularly those with ocean)
mod.sar.op29 <- update(mod.sar.op28, ~. -sdSST.1deg:sdSal.0m)
anova(mod.sar.op29, mod.sar.op28)
# nothing else can be removed

mod.sar.opf <- mod.sar.op28
rm(mod.sar.op28, mod.sar.op29)

## 2v. Create a plot of model parameters --------------------------------------
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.88968

# generate a dataframe of coefficients
(ms.coef <- data.frame(names = names(mod.sar.opf$coefficients), coef.sar = mod.sar.opf$coefficients, row.names = 1:length(mod.sar.opf$coefficients), stars = NA))

# reorder the rows to something more sensible
order.coef.ms <- c(1:23, 40:45, 24:39)
(ms.coef <- ms.coef[order.coef.ms,])
rm(order.coef.ms)

# add a column of significance stars
stars <- c(0.001, 0.01, 0.05, 0.1)
names(stars) <- c("***", "**", "*", ".")
for (i in 1:length(stars)) {
  ms.coef$stars[which(summary(mod.sar.opf)$Coef[, 4] <= stars[i] & is.na(ms.coef$stars))] <- names(stars)[i]
}
rm(i)

# plot the absolute coefficients
png("Figures/Ana_2v_coef_modsaropf_2.png", width = 1000, height = 750)
plt.def <- par("plt")
par(plt = c(plt.def[1:2], 0.5, plt.def[4]))
tmp.x <- barplot(abs(ms.coef$coef.sar), names = ms.coef$names, las = 2, ylim = c(0, max(abs(ms.coef$coef.sar)) + 50))
text(tmp.x, abs(ms.coef$coef.sar) + 20, ms.coef$stars)
par(plt = plt.def)
dev.off()
rm(plt.def, tmp.x, ms.coef)

## 2vi. Calculate likelihood ratios for the SAR model ----------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.88968
AIC(mod.sar.opf) # 4215.515

# removing mean temp^3
mod.sar.lr.mnt3 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 2) - poly(meanSST.1deg, 3):sdSal.0m + poly(meanSST.1deg, 2):sdSal.0m - poly(meanSST.1deg, 3):Ocean2 + poly(meanSST.1deg, 2):Ocean2)
summary(mod.sar.lr.mnt3, Nagelkerke = T) # r2 = 0.87958
AIC(mod.sar.lr.mnt3) # 4300.394
lrtest(mod.sar.opf, mod.sar.lr.mnt3) # n.b. this is basically an anova
# LR = < 2.338e-16 ***

# removing mean temp^2
mod.sar.lr.mnt2 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 1) - poly(meanSST.1deg, 3):sdSal.0m + poly(meanSST.1deg, 1):sdSal.0m - Ocean2:poly(meanSST.1deg, 3) + Ocean2:poly(meanSST.1deg, 1) - prop2.oxy:poly(meanSST.1deg, 2) + prop2.oxy:poly(meanSST.1deg, 1) - absMnSal.0m:poly(meanSST.1deg, 2) + absMnSal.0m:poly(meanSST.1deg, 1) )
summary(mod.sar.lr.mnt2, Nagelkerke = T) # r2 = 0.86771
AIC(mod.sar.lr.mnt2) # 4388.045
lrtest(mod.sar.opf, mod.sar.lr.mnt2)
# LR = < 2.2e-16 ***

# removing mean temp
mod.sar.lr.mnt <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) - poly(meanSST.1deg, 3):sdSal.0m - Ocean2:poly(meanSST.1deg, 3) - prop2.oxy:poly(meanSST.1deg, 2) - absMnSal.0m:poly(meanSST.1deg, 2) - sdSST.1deg:poly(meanSST.1deg, 1) - logProd.mn.ann:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mnt, Nagelkerke = T) # r2 = 0.83188
AIC(mod.sar.lr.mnt) # 4626.084
lrtest(mod.sar.opf, mod.sar.lr.mnt)
# LR = < 2.2e-16 ***

# removing sd temp
mod.sar.lr.sdt <- update(mod.sar.opf, ~. -sdSST.1deg - sdSST.1deg:I(mean.mld.t/10) - sdSST.1deg:sdSal.0m - sdSST.1deg:poly(meanSST.1deg, 1))
summary(mod.sar.lr.sdt, Nagelkerke = T) # r2 = 0.88817
AIC(mod.sar.lr.sdt) # 4221.922
lrtest(mod.sar.opf, mod.sar.lr.sdt)
# LR = 0.006105 **

# removing mld temp
mod.sar.lr.mld <- update(mod.sar.opf, ~. -I(mean.mld.t/10) - sdSST.1deg:I(mean.mld.t/10) - I(mean.mld.t/10):logProd.mn.ann - I(mean.mld.t/10):Ocean2)
summary(mod.sar.lr.mld, Nagelkerke = T) # r2 = 0.88612
AIC(mod.sar.lr.mld) # 4239.213
lrtest(mod.sar.opf, mod.sar.lr.mld)
# LR = 2.735e-06 ***

# removing depth of 10 degree contour
mod.sar.lr.d10 <- update(mod.sar.opf, ~. -I(depth10deg/100) - I(depth10deg/100):absMnSal.0m - I(depth10deg/100):sdSal.0m)
summary(mod.sar.lr.d10, Nagelkerke = T) # r2 = 0.88675 
AIC(mod.sar.lr.d10) # 4237.339
lrtest(mod.sar.opf, mod.sar.lr.d10)
# LR =  3.955e-06 ***

# removing mean log Prod
mod.sar.lr.prod <- update(mod.sar.opf, ~. -logProd.mn.ann - I(mean.mld.t/10):logProd.mn.ann - logProd.mn.ann:sdSal.0m - logProd.mn.ann:Ocean2 - logProd.mn.ann:poly(meanSST.1deg, 1))
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.88361
AIC(mod.sar.lr.prod) # 4260.249
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = 2.068e-10 ***

# removing mean salinity
mod.sar.lr.msal <- update(mod.sar.opf, ~. -absMnSal.0m - I(depth10deg/100):absMnSal.0m - absMnSal.0m:prop2.oxy - absMnSal.0m:Ocean2 - absMnSal.0m:poly(meanSST.1deg, 2))
summary(mod.sar.lr.msal, Nagelkerke = T) # r2 = 0.88497
AIC(mod.sar.lr.msal) # 4245.839
lrtest(mod.sar.opf, mod.sar.lr.msal)
# LR = 1.85e-07 ***

# removing sd salinity
mod.sar.lr.sdsal <- update(mod.sar.opf, ~. -sdSal.0m - poly(meanSST.1deg, 3):sdSal.0m - sdSST.1deg:sdSal.0m - I(depth10deg/100):sdSal.0m - logProd.mn.ann:sdSal.0m - sdSal.0m:prop2.oxy)
summary(mod.sar.lr.sdsal, Nagelkerke = T) # r2 = 0.88583
AIC(mod.sar.lr.sdsal) # 4235.882
lrtest(mod.sar.opf, mod.sar.lr.sdsal)
# LR = 1.504e-05 ***

# removing prop2.oxy
mod.sar.lr.oxy <- update(mod.sar.opf, ~. -prop2.oxy - absMnSal.0m:prop2.oxy - sdSal.0m:prop2.oxy - prop2.oxy:Ocean2 - prop2.oxy:poly(meanSST.1deg, 2))
summary(mod.sar.lr.oxy, Nagelkerke = T) # r2 = 0.8874
AIC(mod.sar.lr.oxy) # 4223.214
lrtest(mod.sar.opf, mod.sar.lr.oxy)
# LR = 0.002863 **

# removing Ocean2
mod.sar.lr.oce <- update(mod.sar.opf, ~. -Ocean2 - I(mean.mld.t/10):Ocean2 - logProd.mn.ann:Ocean2 - absMnSal.0m:Ocean2 - prop2.oxy:Ocean2 - Ocean2:poly(meanSST.1deg, 3))
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.87856
AIC(mod.sar.lr.oce) # 4285.282
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = 1.615e-14 ***

# removing delta_carb_ion
mod.sar.lr.dis <- update(mod.sar.opf, ~. -delta_carb_ion)
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.88888
AIC(mod.sar.lr.dis) # 4221.173
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.005652 **

# create a vector with these likelihood values
ms.lr <- data.frame(names = c("mnt3", "mnt2", "mnt", "sdt", "mld", "d10", "prod", "msal", "sdsal", "oxy", "oce", "dis"), lr = NA, p = NA, stars = NA)

ms.lr$lr <- sapply(ms.lr$names, function (x) lrtest(mod.sar.opf, eval(parse(text = paste("mod.sar.lr.", x, sep = ""))))$Chisq[2])

ms.lr$p <- sapply(ms.lr$names, function (x) lrtest(mod.sar.opf, eval(parse(text = paste("mod.sar.lr.", x, sep = ""))))$Pr[2])

ms.lr

# add a column of significance stars
stars # a vector of signifance stars
for (i in 1:length(stars)) {
  ms.lr$stars[which(ms.lr$p <= stars[i] & is.na(ms.lr$stars))] <- names(stars)[i]
}
rm(i)

# plot these as a barplot
png("Figures/Ana_2vi_LRatio_sar_opf_2.png", width = 550)
tmp.x <- barplot(ms.lr$lr, names = ms.lr$names, las = 2, ylim = c(0, max(ms.lr$lr) + 30))
text(tmp.x, ms.lr$lr + 20, ms.lr$stars)
dev.off()
rm(tmp.x)

## 2vii. Calculate likelihood ratios for groups of EVs ---------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.88968
AIC(mod.sar.opf) # 4215.515

# removing temp
mod.sar.lr.temp <- mod.sar.lr.mnt
summary(mod.sar.lr.temp, Nagelkerke = T) # r2 = 0.83188
AIC(mod.sar.lr.temp) # 4626.084
lrtest(mod.sar.opf, mod.sar.lr.temp)
# LR = < 2.2e-16 ***

# removing structure
mod.sar.lr.str <- update(mod.sar.opf, ~. -I(mean.mld.t/10) - I(depth10deg/100) - sdSST.1deg:I(mean.mld.t/10) - I(mean.mld.t/10):logProd.mn.ann - I(mean.mld.t/10):Ocean2 - I(depth10deg/100):absMnSal.0m - I(depth10deg/100):sdSal.0m)
summary(mod.sar.lr.str, Nagelkerke = T) # r2 = 0.88397
AIC(mod.sar.lr.str) # 4252.977
lrtest(mod.sar.opf, mod.sar.lr.str)
# LR = 8.777e-09 ***

# removing stability
mod.sar.lr.stable <- update(mod.sar.opf, ~. -sdSST.1deg - sdSal.0m - sdSST.1deg:I(mean.mld.t/10) - sdSST.1deg:sdSal.0m - sdSST.1deg:poly(meanSST.1deg, 1) - poly(meanSST.1deg, 3):sdSal.0m - I(depth10deg/100):sdSal.0m - logProd.mn.ann:sdSal.0m - sdSal.0m:prop2.oxy)
summary(mod.sar.lr.stable, Nagelkerke = T) # r2 = 0.88512
AIC(mod.sar.lr.stable) # 4236.42
lrtest(mod.sar.opf, mod.sar.lr.stable)
# LR =  1.128e-05 ***

# removing productivity
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.88361
AIC(mod.sar.lr.prod) # 4260.249
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = 2.068e-10 ***

# removing stress
mod.sar.lr.sts <- update(mod.sar.opf, ~. -absMnSal.0m - prop2.oxy - I(depth10deg/100):absMnSal.0m - absMnSal.0m:prop2.oxy - absMnSal.0m:Ocean2 - sdSal.0m:prop2.oxy - prop2.oxy:Ocean2 - absMnSal.0m:poly(meanSST.1deg, 2) - prop2.oxy:poly(meanSST.1deg, 2))
summary(mod.sar.lr.sts, Nagelkerke = T) # r2 = 0.88383
AIC(mod.sar.lr.sts) # 4244.241
lrtest(mod.sar.opf, mod.sar.lr.sts)
# LR =  4.514e-07 ***

# removing Ocean2
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.87856
AIC(mod.sar.lr.oce) # 4285.282
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = 1.615e-14 ***

# removing delta_carb_ion
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.88888
AIC(mod.sar.lr.dis) # 4221.173
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.005652 **

# create a vector with these likelihood values
ms.lr.group <- data.frame(names = c("temp", "str", "stable", "prod", "sts", "oce", "dis"), lr = NA, p = NA, stars = NA)

ms.lr.group$lr <- sapply(ms.lr.group$names, function (x) lrtest(mod.sar.opf, eval(parse(text = paste("mod.sar.lr.", x, sep = ""))))$Chisq[2])

ms.lr.group$p <- sapply(ms.lr.group$names, function (x) lrtest(mod.sar.opf, eval(parse(text = paste("mod.sar.lr.", x, sep = ""))))$Pr[2])

ms.lr.group

# add a column of significance stars
stars # a vector of signifance stars
for (i in 1:length(stars)) {
  ms.lr.group$stars[which(ms.lr.group$p <= stars[i] & is.na(ms.lr.group$stars))] <- names(stars)[i]
}
rm(i)

# plot these as a barplot
png("Figures/Ana_2vii_LRatio_sar_opf_group_2.png", width = 550)
tmp.x <- barplot(ms.lr.group$lr, names = ms.lr.group$names, las = 2, ylim = c(0, max(ms.lr.group$lr) + 30))
text(tmp.x, ms.lr.group$lr + 20, ms.lr.group$stars)
dev.off()
rm(tmp.x)

rm(mod.sar.lr.d10, mod.sar.lr.dis, mod.sar.lr.mld, mod.sar.lr.mnt, mod.sar.lr.mnt2, mod.sar.lr.mnt3, mod.sar.lr.msal, mod.sar.lr.oce, mod.sar.lr.prod, mod.sar.lr.oxy, mod.sar.lr.sts, mod.sar.lr.sdsal, mod.sar.lr.sdt, mod.sar.lr.stable, mod.sar.lr.str, mod.sar.lr.temp)

## 2viii. Could we just use full models? -----------------------------------
# compare full model and simplified model for all coefficients
lr.sar.op0 <- lr.calc(mod.sar.op0)
par(mfrow = c(1,2))
sar.plot(mod.sar.op0)

lr.sar.opf <- lr.calc(mod.sar.opf) # n.b. this produces the same result as ms.lr
sar.plot(mod.sar.opf)
par(mfrow = c(1,1))

png("Figures/Ana_2viii_LRatio_opfop0_2.png", width = 800)
# get order from running without order first
lr.plot(lr.sar.op0, lr.sar.opf, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.op0g <- lr.calc(mod.sar.op0, tmp)
rm(tmp)

(tmp <- data.frame(names = model.evs(mod.sar.opf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.opfg <- lr.calc(mod.sar.opf, tmp)
rm(tmp)

png("Figures/Ana_2viii_LRatio_g_opfop0_2.png", width = 800)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(6:7, 4:3, 5, 2:1), leg.x = 17, leg.y = 400, leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

## 2ix. Does optimisation differ for simplified? ---------------------------
mod.fw.rop <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fw.rop # 506.9646
AIC(mod.fw.rop$obj)
# 4215.42

mod.fb.rop <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fb.rop # 506.3027
AIC(mod.fb.rop$obj)
# 4256.436

mod.fs.rop <- with(rsr.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fs.rop # 506.9893
AIC(mod.fs.rop$obj)
# 4222.384

# therefore justified in using W
rm(mod.fw.rop, mod.fb.rop, mod.fs.rop, ms.lr, ms.lr.group)


## 3. How do the different Oceans compare? ---------------------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.88968
AIC(mod.sar.opf) # 4215.515

# do I have sufficient points for each ocean?
table(rsr.margo.mod$Ocean2) # should do
# c.f. Atlantic: 372   Indian: 157  Pacific: 146, which were the values for bfd

## 3i. set up model for Atlantic -----------------------------
# run model based on best model for complete data with only Atlantic
atl.nb <- dnearneigh(ldg.coords[rsr.margo.mod$Ocean2 == "Atlantic", ], 0, mod.sar.opW$dist, longlat = TRUE)
atl.w <- nb2listw(atl.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(atl.nb)

## 3ii. set up model for Indian -----------------------------------
# run model based on best model for complete data with only Indian
ind.nb <- dnearneigh(ldg.coords[rsr.margo.mod$Ocean2 == "Indian", ], 0, mod.sar.opW$dist, longlat = TRUE)
ind.w <- nb2listw(ind.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(ind.nb)

## 3iii. set up model for Pacific  ---------------------------
pac.nb <- dnearneigh(ldg.coords[rsr.margo.mod$Ocean2 == "Pacific", ], 0, mod.sar.opW$dist, longlat = TRUE)
pac.w <- nb2listw(pac.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(pac.nb)

# n.b. removed the simplification in 3i-3iii as I'm only simplifying from the model with significant interactions (see 3vii). Therefore, also didn't need 3iv - vi

## 3vii. only significant interactions for Atlantic -----------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)
# get formula without the Ocean2
op.formula <- update(mod.sar.opf$call$formula, ~.-Ocean2 - poly(meanSST.1deg, 3):Ocean2 - I(mean.mld.t/10):Ocean2 - logProd.mn.ann:Ocean2 - absMnSal.0m:Ocean2 - prop2.oxy:Ocean2)

# create a model with only significant values
mod.sar.atlI <- errorsarlm(op.formula, listw = atl.w, zero.policy = TRUE, tol.solve = 1e-18, data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.atlI)
par(mfrow = c(1, 1))

# try adding in the other interactions
mod.sar.atlI2 <- update(mod.sar.atlI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI2)
anova(mod.sar.atlI, mod.sar.atlI2) # 0.20469
rm(mod.sar.atlI2)

mod.sar.atlI3 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.atlI3)
anova(mod.sar.atlI, mod.sar.atlI3) # 0.19816
rm(mod.sar.atlI3)

mod.sar.atlI4 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.atlI4)
anova(mod.sar.atlI, mod.sar.atlI4) # 0.66431
rm(mod.sar.atlI4)

mod.sar.atlI5 <- update(mod.sar.atlI, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI5)
anova(mod.sar.atlI, mod.sar.atlI5) # 0.78754
rm(mod.sar.atlI5)

# significant at the 0.05 level
mod.sar.atlI6 <- update(mod.sar.atlI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.atlI6)
anova(mod.sar.atlI, mod.sar.atlI6) # 0.040002
rm(mod.sar.atlI6)

mod.sar.atlI7 <- update(mod.sar.atlI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.atlI7)
anova(mod.sar.atlI, mod.sar.atlI7) # 0.58331
rm(mod.sar.atlI7)

mod.sar.atlI10 <- update(mod.sar.atlI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.atlI10)
anova(mod.sar.atlI, mod.sar.atlI10) # 0.43145
rm(mod.sar.atlI10)

mod.sar.atlI11 <- update(mod.sar.atlI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.atlI11)
anova(mod.sar.atlI, mod.sar.atlI11) # 0.56986
rm(mod.sar.atlI11)

mod.sar.atlI12 <- update(mod.sar.atlI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.atlI12)
anova(mod.sar.atlI, mod.sar.atlI12) # 0.45564
rm(mod.sar.atlI12)

mod.sar.atlI13 <- update(mod.sar.atlI, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.atlI13)
anova(mod.sar.atlI, mod.sar.atlI13) # 0.60732
rm(mod.sar.atlI13)

mod.sar.atlI14 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.atlI14)
anova(mod.sar.atlI, mod.sar.atlI14) # 0.32214
rm(mod.sar.atlI14)

mod.sar.atlI15 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.atlI15)
anova(mod.sar.atlI, mod.sar.atlI15) # 0.76089
rm(mod.sar.atlI15)

# significant at the 0.05 level
mod.sar.atlI16 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.atlI16)
anova(mod.sar.atlI, mod.sar.atlI16) # 0.54891
rm(mod.sar.atlI16)

mod.sar.atlI17 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.atlI17)
anova(mod.sar.atlI, mod.sar.atlI17) # 0.8422
rm(mod.sar.atlI17)

mod.sar.atlI18 <- update(mod.sar.atlI, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.atlI18)
anova(mod.sar.atlI, mod.sar.atlI18) # 0.19662
rm(mod.sar.atlI18)

mod.sar.atlI19 <- update(mod.sar.atlI, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.atlI19)
anova(mod.sar.atlI, mod.sar.atlI19) # 0.7125
rm(mod.sar.atlI19)

# significant at the 0.05 level
mod.sar.atlI22 <- update(mod.sar.atlI, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.atlI22)
anova(mod.sar.atlI, mod.sar.atlI22) # 0.028451

mod.sar.atlI23 <- update(mod.sar.atlI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.atlI23)
anova(mod.sar.atlI, mod.sar.atlI23) # 0.72343
rm(mod.sar.atlI23)

mod.sar.atlI24 <- update(mod.sar.atlI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.atlI24)
anova(mod.sar.atlI, mod.sar.atlI24) # 0.28158
rm(mod.sar.atlI24)


## adding back in finds logProd.mn.ann:absMnSal.0m should be included in the model, also potentially poly(meanSST.1deg, 3):absMnSal.0m
mod.sar.atl1.I <- mod.sar.atlI22
rm(mod.sar.atlI22)

# check whether any of the interactions are now significant
mod.sar.atl1.I2 <- update(mod.sar.atl1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atl1.I2)
anova(mod.sar.atl1.I, mod.sar.atl1.I2) #  0.14994
rm(mod.sar.atl1.I2)

mod.sar.atl1.I3 <- update(mod.sar.atl1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.atl1.I3)
anova(mod.sar.atl1.I, mod.sar.atl1.I3) # 0.16992
rm(mod.sar.atl1.I3)

mod.sar.atl1.I4 <- update(mod.sar.atl1.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.atl1.I4)
anova(mod.sar.atl1.I, mod.sar.atl1.I4) # 0.43498
rm(mod.sar.atl1.I4)

mod.sar.atl1.I5 <- update(mod.sar.atl1.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atl1.I5)
anova(mod.sar.atl1.I, mod.sar.atl1.I5) # 0.85363
rm(mod.sar.atl1.I5)

mod.sar.atl1.I6 <- update(mod.sar.atl1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.atl1.I6)
anova(mod.sar.atl1.I, mod.sar.atl1.I6) # 0.079059
rm(mod.sar.atl1.I6)

mod.sar.atl1.I7 <- update(mod.sar.atl1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.atl1.I7)
anova(mod.sar.atl1.I, mod.sar.atl1.I7) # 0.57468
rm(mod.sar.atl1.I7)

mod.sar.atl1.I10 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.atl1.I10)
anova(mod.sar.atl1.I, mod.sar.atl1.I10) # 0.34045
rm(mod.sar.atl1.I10)

mod.sar.atl1.I11 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.atl1.I11)
anova(mod.sar.atl1.I, mod.sar.atl1.I11) # 0.6716
rm(mod.sar.atl1.I11)

mod.sar.atl1.I12 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.atl1.I12)
anova(mod.sar.atl1.I, mod.sar.atl1.I12) # 0.19993
rm(mod.sar.atl1.I12)

mod.sar.atl1.I13 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.atl1.I13)
anova(mod.sar.atl1.I, mod.sar.atl1.I13) # 0.53017
rm(mod.sar.atl1.I13)

mod.sar.atl1.I14 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.atl1.I14)
anova(mod.sar.atl1.I, mod.sar.atl1.I14) # 0.23246
rm(mod.sar.atl1.I14)

mod.sar.atl1.I15 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.atl1.I15)
anova(mod.sar.atl1.I, mod.sar.atl1.I15) # 0.48123
rm(mod.sar.atl1.I15)

mod.sar.atl1.I16 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.atl1.I16)
anova(mod.sar.atl1.I, mod.sar.atl1.I16) # 0.63946
rm(mod.sar.atl1.I16)

mod.sar.atl1.I17 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.atl1.I17)
anova(mod.sar.atl1.I, mod.sar.atl1.I17) # 0.95827
rm(mod.sar.atl1.I17)

mod.sar.atl1.I18 <- update(mod.sar.atl1.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.atl1.I18)
anova(mod.sar.atl1.I, mod.sar.atl1.I18) # 0.32186
rm(mod.sar.atl1.I18)

mod.sar.atl1.I19 <- update(mod.sar.atl1.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.atl1.I19)
anova(mod.sar.atl1.I, mod.sar.atl1.I19) # 0.62866
rm(mod.sar.atl1.I19)

mod.sar.atl1.I23 <- update(mod.sar.atl1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.atl1.I23)
anova(mod.sar.atl1.I, mod.sar.atl1.I23) # 0.52192
rm(mod.sar.atl1.I23)

mod.sar.atl1.I24 <- update(mod.sar.atl1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.atl1.I24)
anova(mod.sar.atl1.I, mod.sar.atl1.I24) # 0.92306
rm(mod.sar.atl1.I24)

# nothing is significant, so 
mod.sar.atlIf <- mod.sar.atl1.I
summary(mod.sar.atlIf, Nagelkerke = TRUE) # 0.9217 
AIC(mod.sar.atlIf) # 2509.672

# Check for correlation
env.var.atl <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion", "prop2.oxy")

# pairs plot
png("Figures/Ana_3vii_atl_pairs_2.png", 1200, 1200)
pairs(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", names(rsr.margo.mod) %in% env.var.atl])
dev.off()

# variance inflation factor
vif(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", names(rsr.margo.mod) %in% env.var.atl])
cor(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", names(rsr.margo.mod) %in% env.var.atl])

lr.sar.atlIf <- lr.calc(mod.sar.atlIf)

png("Figures/Ana_3vii_LRatio_atlIf_2.png", width = 800)
lr.plot(lr.sar.atlIf, order = c(7:9, 10:11, 1, 4, 2:3, 5:6), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.atlIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.atlIfg <- lr.calc(mod.sar.atlIf, tmp)
rm(tmp)

png("Figures/Ana_3vii_LRatio_g_atlIf_2.png", width = 800)
lr.plot(lr.sar.atlIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.atlI, mod.sar.atl1.I, file = "Outputs/Atlantic_simplification2.RData")
rm(env.var.atl, mod.sar.atlI, mod.sar.atl1.I)

## 3viii. Only significant interactions for Indian -------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.indI <- errorsarlm(op.formula, listw = ind.w, zero.policy = TRUE, tol.solve = 1e-18, data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.indI)
par(mfrow = c(1, 1))
summary(mod.sar.indI, Nagelkerke = TRUE) # 0.90439 
AIC(mod.sar.indI) # 617.1263

# try adding in the other interactions
mod.sar.indI2 <- update(mod.sar.indI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.indI2)
anova(mod.sar.indI, mod.sar.indI2) # 0.19973
rm(mod.sar.indI2)

mod.sar.indI3 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.indI3)
anova(mod.sar.indI, mod.sar.indI3) # 0.53559
rm(mod.sar.indI3)

mod.sar.indI4 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.indI4)
anova(mod.sar.indI, mod.sar.indI4) # 0.57014
rm(mod.sar.indI4)

mod.sar.indI5 <- update(mod.sar.indI, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.indI5)
anova(mod.sar.indI, mod.sar.indI5) # 0.92648
rm(mod.sar.indI5)

mod.sar.indI6 <- update(mod.sar.indI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.indI6)
anova(mod.sar.indI, mod.sar.indI6) # 0.65188
rm(mod.sar.indI6)

mod.sar.indI7 <- update(mod.sar.indI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.indI7)
anova(mod.sar.indI, mod.sar.indI7) # 0.61373
rm(mod.sar.indI7)

mod.sar.indI10 <- update(mod.sar.indI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.indI10)
anova(mod.sar.indI, mod.sar.indI10) # 0.13818
rm(mod.sar.indI10)

mod.sar.indI11 <- update(mod.sar.indI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.indI11)
anova(mod.sar.indI, mod.sar.indI11) # 0.70847
rm(mod.sar.indI11)

mod.sar.indI12 <- update(mod.sar.indI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.indI12)
anova(mod.sar.indI, mod.sar.indI12) # 0.87202
rm(mod.sar.indI12)

mod.sar.indI13 <- update(mod.sar.indI, ~. + sdSST.1deg:prop2.oxy)
summary(mod.sar.indI13)
anova(mod.sar.indI, mod.sar.indI13) # 0.48376
rm(mod.sar.indI13)

mod.sar.indI14 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.indI14)
anova(mod.sar.indI, mod.sar.indI14) # 0.42464
rm(mod.sar.indI14)

mod.sar.indI15 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.indI15)
anova(mod.sar.indI, mod.sar.indI15) # 0.61519
rm(mod.sar.indI15)

mod.sar.indI16 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.indI16)
anova(mod.sar.indI, mod.sar.indI16) # 0.36236
rm(mod.sar.indI16)

# significant at the 0.05 level
mod.sar.indI17 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.indI17)
anova(mod.sar.indI, mod.sar.indI17) # 0.007008

mod.sar.indI18 <- update(mod.sar.indI, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.indI18)
anova(mod.sar.indI, mod.sar.indI18) # 0.77059
rm(mod.sar.indI18)

mod.sar.indI19 <- update(mod.sar.indI, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.indI19)
anova(mod.sar.indI, mod.sar.indI19) # 0.80413
rm(mod.sar.indI19)

mod.sar.indI22 <- update(mod.sar.indI, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.indI22)
anova(mod.sar.indI, mod.sar.indI22) # 0.72173
rm(mod.sar.indI22)

mod.sar.indI23 <- update(mod.sar.indI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.indI23)
anova(mod.sar.indI, mod.sar.indI23) # 0.27164
rm(mod.sar.indI23)

mod.sar.indI24 <- update(mod.sar.indI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.indI24)
anova(mod.sar.indI, mod.sar.indI24) # 0.46073
rm(mod.sar.indI24)


# I(mean.mld.t/10):prop2.oxy is significant at the 0.05 level, so add it in and recheck
mod.sar.ind1.I <- mod.sar.indI17
rm(mod.sar.indI17)

mod.sar.ind1.I2 <- update(mod.sar.ind1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind1.I2)
anova(mod.sar.ind1.I, mod.sar.ind1.I2) # 0.19606
rm(mod.sar.ind1.I2)

# significant at the 0.05 level
mod.sar.ind1.I3 <- update(mod.sar.ind1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.ind1.I3)
anova(mod.sar.ind1.I, mod.sar.ind1.I3) # 0.021843

mod.sar.ind1.I4 <- update(mod.sar.ind1.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind1.I4)
anova(mod.sar.ind1.I, mod.sar.ind1.I4) # 0.5202
rm(mod.sar.ind1.I4)

mod.sar.ind1.I5 <- update(mod.sar.ind1.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind1.I5)
anova(mod.sar.ind1.I, mod.sar.ind1.I5) # 0.76029
rm(mod.sar.ind1.I5)

mod.sar.ind1.I6 <- update(mod.sar.ind1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind1.I6)
anova(mod.sar.ind1.I, mod.sar.ind1.I6) # 0.43989
rm(mod.sar.ind1.I6)

mod.sar.ind1.I7 <- update(mod.sar.ind1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind1.I7)
anova(mod.sar.ind1.I, mod.sar.ind1.I7) # 0.29931
rm(mod.sar.ind1.I7)

mod.sar.ind1.I10 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind1.I10)
anova(mod.sar.ind1.I, mod.sar.ind1.I10) # 0.097313
rm(mod.sar.ind1.I10)

mod.sar.ind1.I11 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind1.I11)
anova(mod.sar.ind1.I, mod.sar.ind1.I11) # 0.99599
rm(mod.sar.ind1.I11)

mod.sar.ind1.I12 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind1.I12)
anova(mod.sar.ind1.I, mod.sar.ind1.I12) # 0.33137
rm(mod.sar.ind1.I12)

mod.sar.ind1.I13 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.ind1.I13)
anova(mod.sar.ind1.I, mod.sar.ind1.I13) # 0.93528
rm(mod.sar.ind1.I13)

mod.sar.ind1.I14 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind1.I14)
anova(mod.sar.ind1.I, mod.sar.ind1.I14) # 0.431 
rm(mod.sar.ind1.I14)

# significant at the 0.05 level
mod.sar.ind1.I15 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind1.I15)
anova(mod.sar.ind1.I, mod.sar.ind1.I15) # 0.049032
rm(mod.sar.ind1.I15)

mod.sar.ind1.I16 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind1.I16)
anova(mod.sar.ind1.I, mod.sar.ind1.I16) # 0.83814
rm(mod.sar.ind1.I16)

mod.sar.ind1.I18 <- update(mod.sar.ind1.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.ind1.I18)
anova(mod.sar.ind1.I, mod.sar.ind1.I18) # 0.87338 
rm(mod.sar.ind1.I18)

mod.sar.ind1.I19 <- update(mod.sar.ind1.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.ind1.I19)
anova(mod.sar.ind1.I, mod.sar.ind1.I19) # 0.75593
rm(mod.sar.ind1.I19)

mod.sar.ind1.I22 <- update(mod.sar.ind1.I, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.ind1.I22)
anova(mod.sar.ind1.I, mod.sar.ind1.I22) # 0.96835
rm(mod.sar.ind1.I22)

mod.sar.ind1.I23 <- update(mod.sar.ind1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind1.I23)
anova(mod.sar.ind1.I, mod.sar.ind1.I23) # 0.56976
rm(mod.sar.ind1.I23)

mod.sar.ind1.I24 <- update(mod.sar.ind1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind1.I24)
anova(mod.sar.ind1.I, mod.sar.ind1.I24) # 0.53796
rm(mod.sar.ind1.I24)


# poly(meanSST.1deg, 3):I(mean.mld.t/10) is significant at the 0.05 level, as is I(mean.mld.t/10):absMnSal.0m so add it in and recheck
mod.sar.ind2.I <- mod.sar.ind1.I3
rm(mod.sar.ind1.I3)

mod.sar.ind2.I2 <- update(mod.sar.ind2.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind2.I2)
anova(mod.sar.ind2.I, mod.sar.ind2.I2) # 0.68073
rm(mod.sar.ind2.I2)

mod.sar.ind2.I4 <- update(mod.sar.ind2.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind2.I4)
anova(mod.sar.ind2.I, mod.sar.ind2.I4) # 0.79285
rm(mod.sar.ind2.I4)

mod.sar.ind2.I5 <- update(mod.sar.ind2.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind2.I5)
anova(mod.sar.ind2.I, mod.sar.ind2.I5) # 0.88454
rm(mod.sar.ind2.I5)

mod.sar.ind2.I6 <- update(mod.sar.ind2.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind2.I6)
anova(mod.sar.ind2.I, mod.sar.ind2.I6) # 0.427
rm(mod.sar.ind2.I6)

mod.sar.ind2.I7 <- update(mod.sar.ind2.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind2.I7)
anova(mod.sar.ind2.I, mod.sar.ind2.I7) # 0.13061
rm(mod.sar.ind2.I7)

mod.sar.ind2.I10 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind2.I10)
anova(mod.sar.ind2.I, mod.sar.ind2.I10) # 0.15332
rm(mod.sar.ind2.I10)

mod.sar.ind2.I11 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind2.I11)
anova(mod.sar.ind2.I, mod.sar.ind2.I11) # 0.87863
rm(mod.sar.ind2.I11)

mod.sar.ind2.I12 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind2.I12)
anova(mod.sar.ind2.I, mod.sar.ind2.I12) # 0.51937
rm(mod.sar.ind2.I12)

mod.sar.ind2.I13 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.ind2.I13)
anova(mod.sar.ind2.I, mod.sar.ind2.I13) # 0.51173
rm(mod.sar.ind2.I13)

mod.sar.ind2.I14 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind2.I14)
anova(mod.sar.ind2.I, mod.sar.ind2.I14) # 0.16624
rm(mod.sar.ind2.I14)

mod.sar.ind2.I15 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind2.I15)
anova(mod.sar.ind2.I, mod.sar.ind2.I15) # 0.8972
rm(mod.sar.ind2.I15)

mod.sar.ind2.I16 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind2.I16)
anova(mod.sar.ind2.I, mod.sar.ind2.I16) # 0.57742
rm(mod.sar.ind2.I16)

mod.sar.ind2.I18 <- update(mod.sar.ind2.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.ind2.I18)
anova(mod.sar.ind2.I, mod.sar.ind2.I18) # 0.5929
rm(mod.sar.ind2.I18)

mod.sar.ind2.I19 <- update(mod.sar.ind2.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.ind2.I19)
anova(mod.sar.ind2.I, mod.sar.ind2.I19) # 0.30637
rm(mod.sar.ind2.I19)

mod.sar.ind2.I22 <- update(mod.sar.ind2.I, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.ind2.I22)
anova(mod.sar.ind2.I, mod.sar.ind2.I22) # 0.9621
rm(mod.sar.ind2.I22)

mod.sar.ind2.I23 <- update(mod.sar.ind2.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind2.I23)
anova(mod.sar.ind2.I, mod.sar.ind2.I23) # 0.72524
rm(mod.sar.ind2.I23)

mod.sar.ind2.I24 <- update(mod.sar.ind2.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind2.I24)
anova(mod.sar.ind2.I, mod.sar.ind2.I24) # 0.50764
rm(mod.sar.ind2.I24)



# Nothing more is significant
mod.sar.indIf <- mod.sar.ind2.I
summary(mod.sar.indIf, Nagelkerke = TRUE) # 0.91428
AIC(mod.sar.indIf) # 608.2109

# Check for correlation
env.var.ind <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion")

# pairs plot
png("Figures/Ana_3viii_ind_pairs_2.png", 1200, 1200)
pairs(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", names(rsr.margo.mod) %in% env.var.ind])
dev.off()

# variance inflation factor
vif(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", names(rsr.margo.mod) %in% env.var.ind])
cor(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", names(rsr.margo.mod) %in% env.var.ind])
# correlation between SST and MLD

lr.sar.indIf <- lr.calc(mod.sar.indIf)

png("Figures/Ana_3viii_LRatio_indIf_2.png", width = 800)
lr.plot(lr.sar.indIf, order = c(7:9, 10:11, 1, 4, 2:3, 5:6), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.indIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.indIfg <- lr.calc(mod.sar.indIf, tmp)
rm(tmp)

png("Figures/Ana_3viii_LRatio_g_indIf_2.png", width = 800)
lr.plot(lr.sar.indIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

save(mod.sar.indI, mod.sar.ind1.I, mod.sar.ind2.I, file = "Outputs/Indian_simplification2.RData")
rm(env.var.ind, mod.sar.indI, mod.sar.ind1.I, mod.sar.ind2.I)

## 3ix. Only significant interactions for Pacific --------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.pacI <- errorsarlm(op.formula, listw = pac.w, zero.policy = TRUE, tol.solve = 1e-18, data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.pacI)
par(mfrow = c(1, 1))

# try adding in the other interactions
mod.sar.pacI2 <- update(mod.sar.pacI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI2)
anova(mod.sar.pacI, mod.sar.pacI2) # 0.1552
rm(mod.sar.pacI2)

mod.sar.pacI3 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.pacI3)
anova(mod.sar.pacI, mod.sar.pacI3) # 0.054415
rm(mod.sar.pacI3)

# signficant at the 0.05 level
mod.sar.pacI4 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.pacI4)
anova(mod.sar.pacI, mod.sar.pacI4) # 0.0012343

# signficant at the 0.05 level
mod.sar.pacI5 <- update(mod.sar.pacI, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pacI5)
anova(mod.sar.pacI, mod.sar.pacI5) # 0.01917
rm(mod.sar.pacI5)

mod.sar.pacI6 <- update(mod.sar.pacI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pacI6)
anova(mod.sar.pacI, mod.sar.pacI6) # 0.12394
rm(mod.sar.pacI6)

mod.sar.pacI7 <- update(mod.sar.pacI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pacI7)
anova(mod.sar.pacI, mod.sar.pacI7) # 0.86147
rm(mod.sar.pacI7)

mod.sar.pacI10 <- update(mod.sar.pacI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pacI10)
anova(mod.sar.pacI, mod.sar.pacI10) # 0.62706
rm(mod.sar.pacI10)

mod.sar.pacI11 <- update(mod.sar.pacI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pacI11)
anova(mod.sar.pacI, mod.sar.pacI11) # 0.47695
rm(mod.sar.pacI11)

mod.sar.pacI12 <- update(mod.sar.pacI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pacI12)
anova(mod.sar.pacI, mod.sar.pacI12) # 0.051491
rm(mod.sar.pacI12)

mod.sar.pacI13 <- update(mod.sar.pacI, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.pacI13)
anova(mod.sar.pacI, mod.sar.pacI13) # 0.85563
rm(mod.sar.pacI13)

# significant at the 0.05 level
mod.sar.pacI14 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pacI14)
anova(mod.sar.pacI, mod.sar.pacI14) # 0.044108
rm(mod.sar.pacI14)

mod.sar.pacI15 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pacI15)
anova(mod.sar.pacI, mod.sar.pacI15) # 0.047649
rm(mod.sar.pacI15)

mod.sar.pacI16 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pacI16)
anova(mod.sar.pacI, mod.sar.pacI16) # 0.47901
rm(mod.sar.pacI16)

mod.sar.pacI17 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pacI17)
anova(mod.sar.pacI, mod.sar.pacI17) # 0.49358
rm(mod.sar.pacI17)

mod.sar.pacI18 <- update(mod.sar.pacI, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.pacI18)
anova(mod.sar.pacI, mod.sar.pacI18) # 0.33104
rm(mod.sar.pacI18)

mod.sar.pacI19 <- update(mod.sar.pacI, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.pacI19)
anova(mod.sar.pacI, mod.sar.pacI19) # 0.52564
rm(mod.sar.pacI19)

mod.sar.pacI22 <- update(mod.sar.pacI, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.pacI22)
anova(mod.sar.pacI, mod.sar.pacI22) # 0.42922
rm(mod.sar.pacI22)

mod.sar.pacI23 <- update(mod.sar.pacI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pacI23)
anova(mod.sar.pacI, mod.sar.pacI23) # 0.5351
rm(mod.sar.pacI23)

mod.sar.pacI24 <- update(mod.sar.pacI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pacI24)
anova(mod.sar.pacI, mod.sar.pacI24) # 0.47197
rm(mod.sar.pacI24)


# multiple interactions were significant. Add in the most significant (poly(meanSST.1deg, 3):I(depth10deg/100)) and then recalculate
mod.sar.pac1.I <- mod.sar.pacI4
rm(mod.sar.pacI4)

# try adding in the other interactions
mod.sar.pac1.I2 <- update(mod.sar.pac1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac1.I2)
anova(mod.sar.pac1.I, mod.sar.pac1.I2) # 0.75949
rm(mod.sar.pac1.I2)

# significant at the 0.05 level
mod.sar.pac1.I3 <- update(mod.sar.pac1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.pac1.I3)
anova(mod.sar.pac1.I, mod.sar.pac1.I3) # 0.015983

mod.sar.pac1.I5 <- update(mod.sar.pac1.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac1.I5)
anova(mod.sar.pac1.I, mod.sar.pac1.I5) # 0.21003
rm(mod.sar.pac1.I5)

mod.sar.pac1.I6 <- update(mod.sar.pac1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac1.I6)
anova(mod.sar.pac1.I, mod.sar.pac1.I6) # 0.92639
rm(mod.sar.pac1.I6)

mod.sar.pac1.I7 <- update(mod.sar.pac1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac1.I7)
anova(mod.sar.pac1.I, mod.sar.pac1.I7) # 0.75021
rm(mod.sar.pac1.I7)

mod.sar.pac1.I10 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac1.I10)
anova(mod.sar.pac1.I, mod.sar.pac1.I10) # 0.13952
rm(mod.sar.pac1.I10)

mod.sar.pac1.I11 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac1.I11)
anova(mod.sar.pac1.I, mod.sar.pac1.I11) # 0.20474
rm(mod.sar.pac1.I11)

mod.sar.pac1.I12 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac1.I12)
anova(mod.sar.pac1.I, mod.sar.pac1.I12) # 0.54202
rm(mod.sar.pac1.I12)

mod.sar.pac1.I13 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.pac1.I13)
anova(mod.sar.pac1.I, mod.sar.pac1.I13) # 0.54202
rm(mod.sar.pac1.I13)

mod.sar.pac1.I14 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac1.I14)
anova(mod.sar.pac1.I, mod.sar.pac1.I14) # 0.25033
rm(mod.sar.pac1.I14)

mod.sar.pac1.I15 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac1.I15)
anova(mod.sar.pac1.I, mod.sar.pac1.I15) # 0.18371
rm(mod.sar.pac1.I15)

mod.sar.pac1.I16 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac1.I16)
anova(mod.sar.pac1.I, mod.sar.pac1.I16) # 0.68586
rm(mod.sar.pac1.I16)

mod.sar.pac1.I17 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac1.I17)
anova(mod.sar.pac1.I, mod.sar.pac1.I17) # 0.53033
rm(mod.sar.pac1.I17)

mod.sar.pac1.I18 <- update(mod.sar.pac1.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.pac1.I18)
anova(mod.sar.pac1.I, mod.sar.pac1.I18) # 0.20484
rm(mod.sar.pac1.I18)

mod.sar.pac1.I19 <- update(mod.sar.pac1.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.pac1.I19)
anova(mod.sar.pac1.I, mod.sar.pac1.I19) # 0.34684
rm(mod.sar.pac1.I19)

# significant at the 0.05 level
mod.sar.pac1.I22 <- update(mod.sar.pac1.I, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.pac1.I22)
anova(mod.sar.pac1.I, mod.sar.pac1.I22) # 0.024018
rm(mod.sar.pac1.I22)

mod.sar.pac1.I23 <- update(mod.sar.pac1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac1.I23)
anova(mod.sar.pac1.I, mod.sar.pac1.I23) # 0.62774
rm(mod.sar.pac1.I23)

mod.sar.pac1.I24 <- update(mod.sar.pac1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac1.I24)
anova(mod.sar.pac1.I, mod.sar.pac1.I24) # 0.80476
rm(mod.sar.pac1.I24)


# add in poly(meanSST.1deg, 3):I(mean.mld.t/10) and recalculate
mod.sar.pac2.I <- mod.sar.pac1.I3
rm(mod.sar.pac1.I3)

# try adding in the other interactions
mod.sar.pac2.I2 <- update(mod.sar.pac2.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac2.I2)
anova(mod.sar.pac2.I, mod.sar.pac2.I2) # 0.93009
rm(mod.sar.pac2.I2)

mod.sar.pac2.I5 <- update(mod.sar.pac2.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac2.I5)
anova(mod.sar.pac2.I, mod.sar.pac2.I5) # 0.19948
rm(mod.sar.pac2.I5)

mod.sar.pac2.I6 <- update(mod.sar.pac2.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac2.I6)
anova(mod.sar.pac2.I, mod.sar.pac2.I6) # 0.17008
rm(mod.sar.pac2.I6)

mod.sar.pac2.I7 <- update(mod.sar.pac2.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac2.I7)
anova(mod.sar.pac2.I, mod.sar.pac2.I7) # 0.79691
rm(mod.sar.pac2.I7)

mod.sar.pac2.I10 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac2.I10)
anova(mod.sar.pac2.I, mod.sar.pac2.I10) # 0.2253
rm(mod.sar.pac2.I10)

mod.sar.pac2.I11 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac2.I11)
anova(mod.sar.pac2.I, mod.sar.pac2.I11) # 0.19855
rm(mod.sar.pac2.I11)

mod.sar.pac2.I12 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac2.I12)
anova(mod.sar.pac2.I, mod.sar.pac2.I12) # 0.74415
rm(mod.sar.pac2.I12)

mod.sar.pac2.I13 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.pac2.I13)
anova(mod.sar.pac2.I, mod.sar.pac2.I13) # 0.32927
rm(mod.sar.pac2.I13)

# significant at the 0.05 level
mod.sar.pac2.I14 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac2.I14)
anova(mod.sar.pac2.I, mod.sar.pac2.I14) # 0.034302
rm(mod.sar.pac2.I14)

mod.sar.pac2.I15 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac2.I15)
anova(mod.sar.pac2.I, mod.sar.pac2.I15) # 0.50493
rm(mod.sar.pac2.I15)

mod.sar.pac2.I16 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac2.I16)
anova(mod.sar.pac2.I, mod.sar.pac2.I16) # 0.50493
rm(mod.sar.pac2.I16)

mod.sar.pac2.I17 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac2.I17)
anova(mod.sar.pac2.I, mod.sar.pac2.I17) # 0.93129
rm(mod.sar.pac2.I17)

mod.sar.pac2.I18 <- update(mod.sar.pac2.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.pac2.I18)
anova(mod.sar.pac2.I, mod.sar.pac2.I18) # 0.16032
rm(mod.sar.pac2.I18)

mod.sar.pac2.I19 <- update(mod.sar.pac2.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.pac2.I19)
anova(mod.sar.pac2.I, mod.sar.pac2.I19) # 0.24975
rm(mod.sar.pac2.I19)

# significant at the 0.05 level
mod.sar.pac2.I22 <- update(mod.sar.pac2.I, ~. + logProd.mn.ann:absMnSal.0m)
summary(mod.sar.pac2.I22)
anova(mod.sar.pac2.I, mod.sar.pac2.I22) # 0.0054852

mod.sar.pac2.I23 <- update(mod.sar.pac2.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac2.I23)
anova(mod.sar.pac2.I, mod.sar.pac2.I23) # 0.96204
rm(mod.sar.pac2.I23)

mod.sar.pac2.I24 <- update(mod.sar.pac2.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac2.I24)
anova(mod.sar.pac2.I, mod.sar.pac2.I24) # 0.91605
rm(mod.sar.pac2.I24)


# add in logProd.mn.ann:absMnSal.0m and recalculate
mod.sar.pac3.I <- mod.sar.pac2.I22
rm(mod.sar.pac2.I22)

# try adding in the other interactions
mod.sar.pac3.I2 <- update(mod.sar.pac3.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac3.I2)
anova(mod.sar.pac3.I, mod.sar.pac3.I2) # 0.99835
rm(mod.sar.pac3.I2)

mod.sar.pac3.I5 <- update(mod.sar.pac3.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac3.I5)
anova(mod.sar.pac3.I, mod.sar.pac3.I5) # 0.29008
rm(mod.sar.pac3.I5)

mod.sar.pac3.I6 <- update(mod.sar.pac3.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac3.I6)
anova(mod.sar.pac3.I, mod.sar.pac3.I6) # 0.24558
rm(mod.sar.pac3.I6)

mod.sar.pac3.I7 <- update(mod.sar.pac3.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac3.I7)
anova(mod.sar.pac3.I, mod.sar.pac3.I7) # 0.75806
rm(mod.sar.pac3.I7)

mod.sar.pac3.I10 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac3.I10)
anova(mod.sar.pac3.I, mod.sar.pac3.I10) # 0.43313
rm(mod.sar.pac3.I10)

mod.sar.pac3.I11 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac3.I11)
anova(mod.sar.pac3.I, mod.sar.pac3.I11) # 0.11737
rm(mod.sar.pac3.I11)

mod.sar.pac3.I12 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac3.I12)
anova(mod.sar.pac3.I, mod.sar.pac3.I12) # 0.75137
rm(mod.sar.pac3.I12)

mod.sar.pac3.I13 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.pac3.I13)
anova(mod.sar.pac3.I, mod.sar.pac3.I13) # 0.12088
rm(mod.sar.pac3.I13)

# significant at the 0.05 level
mod.sar.pac3.I14 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac3.I14)
anova(mod.sar.pac3.I, mod.sar.pac3.I14) # 0.031459

mod.sar.pac3.I15 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac3.I15)
anova(mod.sar.pac3.I, mod.sar.pac3.I15) # 0.39928
rm(mod.sar.pac3.I15)

mod.sar.pac3.I16 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac3.I16)
anova(mod.sar.pac3.I, mod.sar.pac3.I16) # 0.49168
rm(mod.sar.pac3.I16)

mod.sar.pac3.I17 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac3.I17)
anova(mod.sar.pac3.I, mod.sar.pac3.I17) # 0.87656
rm(mod.sar.pac3.I17)

mod.sar.pac3.I18 <- update(mod.sar.pac3.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.pac3.I18)
anova(mod.sar.pac3.I, mod.sar.pac3.I18) # 0.48776
rm(mod.sar.pac3.I18)

mod.sar.pac3.I19 <- update(mod.sar.pac3.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.pac3.I19)
anova(mod.sar.pac3.I, mod.sar.pac3.I19) # 0.1172
rm(mod.sar.pac3.I19)

mod.sar.pac3.I23 <- update(mod.sar.pac3.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac3.I23)
anova(mod.sar.pac3.I, mod.sar.pac3.I23) # 0.66761
rm(mod.sar.pac3.I23)

mod.sar.pac3.I24 <- update(mod.sar.pac3.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac3.I24)
anova(mod.sar.pac3.I, mod.sar.pac3.I24) # 0.32953
rm(mod.sar.pac3.I24)


# add in logProd.mn.ann:absMnSal.0m and recalculate
mod.sar.pac4.I <- mod.sar.pac3.I14
rm(mod.sar.pac3.I14)

# try adding in the other interactions
mod.sar.pac4.I2 <- update(mod.sar.pac4.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac4.I2)
anova(mod.sar.pac4.I, mod.sar.pac4.I2) # 0.93151
rm(mod.sar.pac4.I2)

mod.sar.pac4.I5 <- update(mod.sar.pac4.I, ~. -logProd.mn.ann:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac4.I5)
anova(mod.sar.pac4.I, mod.sar.pac4.I5) # 0.18488
rm(mod.sar.pac4.I5)

mod.sar.pac4.I6 <- update(mod.sar.pac4.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac4.I6)
anova(mod.sar.pac4.I, mod.sar.pac4.I6) # 0.13364
rm(mod.sar.pac4.I6)

mod.sar.pac4.I7 <- update(mod.sar.pac4.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac4.I7)
anova(mod.sar.pac4.I, mod.sar.pac4.I7) # 0.75704
rm(mod.sar.pac4.I7)

mod.sar.pac4.I10 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac4.I10)
anova(mod.sar.pac4.I, mod.sar.pac4.I10) # 0.62662
rm(mod.sar.pac4.I10)

mod.sar.pac4.I11 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac4.I11)
anova(mod.sar.pac4.I, mod.sar.pac4.I11) # 0.18277
rm(mod.sar.pac4.I11)

mod.sar.pac4.I12 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac4.I12)
anova(mod.sar.pac4.I, mod.sar.pac4.I12) # 0.55785
rm(mod.sar.pac4.I12)

mod.sar.pac4.I13 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:prop2.oxy )
summary(mod.sar.pac4.I13)
anova(mod.sar.pac4.I, mod.sar.pac4.I13) # 0.13717
rm(mod.sar.pac4.I13)

mod.sar.pac4.I15 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac4.I15)
anova(mod.sar.pac4.I, mod.sar.pac4.I15) # 0.15108
rm(mod.sar.pac4.I15)

mod.sar.pac4.I16 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac4.I16)
anova(mod.sar.pac4.I, mod.sar.pac4.I16) # 0.66722
rm(mod.sar.pac4.I16)

mod.sar.pac4.I17 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac4.I17)
anova(mod.sar.pac4.I, mod.sar.pac4.I17) # 0.91997
rm(mod.sar.pac4.I17)

mod.sar.pac4.I18 <- update(mod.sar.pac4.I, ~. + I(depth10deg/100):logProd.mn.ann)
summary(mod.sar.pac4.I18)
anova(mod.sar.pac4.I, mod.sar.pac4.I18) # 0.88878
rm(mod.sar.pac4.I18)

mod.sar.pac4.I19 <- update(mod.sar.pac4.I, ~. + I(depth10deg/100):prop2.oxy)
summary(mod.sar.pac4.I19)
anova(mod.sar.pac4.I, mod.sar.pac4.I19) # 0.37104
rm(mod.sar.pac4.I19)

mod.sar.pac4.I23 <- update(mod.sar.pac4.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac4.I23)
anova(mod.sar.pac4.I, mod.sar.pac4.I23) # 0.7875
rm(mod.sar.pac4.I23)

mod.sar.pac4.I24 <- update(mod.sar.pac4.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac4.I24)
anova(mod.sar.pac4.I, mod.sar.pac4.I24) # 0.55957
rm(mod.sar.pac4.I24)


# don't need to add anything else
mod.sar.pacIf <- mod.sar.pac4.I

summary(mod.sar.pacIf, Nagelkerke = TRUE) # 0.77228 
AIC(mod.sar.pacIf) # 1024.041

# Check for correlation
env.var.pac <- c("meanSST.1deg", "sdSST.1deg", "depth10deg", "mean.mld.t", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion")

# pairs plot
png("Figures/Ana_3ix_pac_pairs_2.png", 1200, 1200)
pairs(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", names(rsr.margo.mod) %in% env.var.pac])
dev.off()

# variance inflation factor
vif(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", names(rsr.margo.mod) %in% env.var.pac])
cor(rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", names(rsr.margo.mod) %in% env.var.pac])

lr.sar.pacIf <- lr.calc(mod.sar.pacIf)

png("Figures/Ana_3ix_LRatio_pacIf_2.png", width = 800)
lr.plot(lr.sar.pacIf, order = c(7:11, 1, 4, 2:3, 5:6), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 5)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.pacIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.pacIfg <- lr.calc(mod.sar.pacIf, tmp)
rm(tmp)

png("Figures/Ana_3ix_LRatio_g_pacIf_2.png", width = 800)
lr.plot(lr.sar.pacIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.pacI, mod.sar.pac1.I, mod.sar.pac2.I, mod.sar.pac3.I, file = "Outputs/Pacific_simplification2.RData")
rm(env.var.pac, mod.sar.pacI, mod.sar.pac1.I, mod.sar.pac2.I, mod.sar.pac3.I)

## 3x. Comparison of the oceans --------------------------------------------
png("Figures/Ana_3x_LRatio_oceIf_2.png", width = 1200)
lr.plot(lr.sar.op0, lr.sar.atlIf, lr.sar.indIf, lr.sar.pacIf, order = c(9:7, 3:4, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
png("Figures/Ana_3x_LRatio_g_OceIf_2.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.atlIfg, lr.sar.indIfg, lr.sar.pacIfg, order = c(6:7, 4, 3, 5, 2:1), leg.txt = c("All", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 15, srt = 50)
dev.off()

## 3xi. Comparing randomly subsampled points ---------------------------------------
table(rsr.margo.mod$Ocean2)

## for the Atlantic
atl.sample <- lr.sar.atlIf
atl.sample.g <- lr.sar.atlIfg
names(atl.sample)[2:4] <- names(atl.sample.g)[2:4] <- paste(names(atl.sample)[2:4], "all", sep = "_")

atl.data <- rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", ]

(var.dat <- data.frame(names = model.evs(mod.sar.atlIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))

# repeat 100 times
for (i in 1:100) {
  # create a subsample of the dataset
  tmp <- atl.data[sample(1:nrow(atl.data), 168), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long,tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.opW$dist, longlat = TRUE)
  tmp.w <- nb2listw(tmp.nb, glist = NULL, style = "W", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.atlIf$call$formula, listw = tmp.w, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  atl.sample <- cbind(atl.sample, tmp2[2:4])
  atl.sample.g <- cbind(atl.sample.g, tmp3[2:4])
  names(atl.sample)[(3 * i + 2):(3 * i + 4)] <- names(atl.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.w, tmp2, tmp3, tmp.mod, i)

atl.final <- lr.sar.atlIf
atl.final$lr <- rowMeans(atl.sample[, seq(5, ncol(atl.sample), by = 3)])
atl.final$p <- rowMeans(atl.sample[, seq(6, ncol(atl.sample), by = 3)])
atl.final$stars <- NA

for (i in 1:length(stars)) {
  atl.final$stars[which(atl.final$p <= stars[i] & is.na(atl.final$stars))] <- names(stars)[i]
}
rm(i)

se.atl <- apply(atl.sample[, seq(5, ncol(atl.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.atl <- apply(atl.sample[, seq(5, ncol(atl.sample), by = 3)], 1, sd)

lr.plot(atl.final, order = c(9:7, 10:11, 1, 4, 2:3, 5:6), ylab = "Log Likelihood ratio", star.pos = 10, se.val = se.atl)

atl.final.g <- lr.sar.atlIfg
atl.final.g$lr <- rowMeans(atl.sample.g[, seq(5, ncol(atl.sample.g), by = 3)])
atl.final.g$p <- rowMeans(atl.sample.g[, seq(6, ncol(atl.sample.g), by = 3)])
atl.final.g$stars <- NA

for (i in 1:length(stars)) {
  atl.final.g$stars[which(atl.final.g$p <= stars[i] & is.na(atl.final.g$stars))] <- names(stars)[i]
}
rm(i)

se.atl.g <- apply(atl.sample.g[, seq(5, ncol(atl.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.atl.g <- apply(atl.sample.g[, seq(5, ncol(atl.sample.g), by = 3)], 1, sd)


lr.plot(atl.final.g, ylab = "Log Likelihood ratio", order = c(5:6, 1:4), star.pos = 10, se.val = se.atl.g)

rm(atl.data)

## for the Pacific
pac.sample <- lr.sar.pacIf
pac.sample.g <- lr.sar.pacIfg
names(pac.sample)[2:4] <- names(pac.sample.g)[2:4] <- paste(names(pac.sample)[2:4], "all", sep = "_")

pac.data <- rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", ]

# repeat 100 times
for (i in 1:100) {
  # create a subsample of the dataset
  tmp <- pac.data[sample(1:nrow(pac.data), 168), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long,tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.opW$dist, longlat = TRUE)
  tmp.w <- nb2listw(tmp.nb, glist = NULL, style = "W", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.pacIf$call$formula, listw = tmp.w, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  pac.sample <- cbind(pac.sample, tmp2[2:4])
  pac.sample.g <- cbind(pac.sample.g, tmp3[2:4])
  names(pac.sample)[(3 * i + 2):(3 * i + 4)] <- names(pac.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.w, tmp2, tmp3, tmp.mod, i)

pac.final <- lr.sar.pacIf
pac.final$lr <- rowMeans(pac.sample[, seq(5, ncol(pac.sample), by = 3)])
pac.final$p <- rowMeans(pac.sample[, seq(6, ncol(pac.sample), by = 3)])
pac.final$stars <- NA

for (i in 1:length(stars)) {
  pac.final$stars[which(pac.final$p <= stars[i] & is.na(pac.final$stars))] <- names(stars)[i]
}

se.pac <- apply(pac.sample[, seq(5, ncol(pac.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.pac <- apply(pac.sample[, seq(5, ncol(pac.sample), by = 3)], 1, sd)

lr.plot(pac.final, ylab = "Log Likelihood ratio", order = c(9:7, 10:11, 1, 4, 2:3, 5:6), star.pos = 5, se.val = se.pac)

pac.final.g <- lr.sar.pacIfg
pac.final.g$lr <- rowMeans(pac.sample.g[, seq(5, ncol(pac.sample.g), by = 3)])
pac.final.g$p <- rowMeans(pac.sample.g[, seq(6, ncol(pac.sample.g), by = 3)])
pac.final.g$stars <- NA

for (i in 1:length(stars)) {
  pac.final.g$stars[which(pac.final.g$p <= stars[i] & is.na(pac.final.g$stars))] <- names(stars)[i]
}

se.pac.g <- apply(pac.sample.g[, seq(5, ncol(pac.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.pac.g <- apply(pac.sample.g[, seq(5, ncol(pac.sample.g), by = 3)], 1, sd)

lr.plot(pac.final.g, ylab = "Log Likelihood ratio", order = c(5:6, 1:4), star.pos = 10, se.val = se.pac.g)

## the last bit of this section won't run until I have also run section 13

# plot the whole set up
sd.val <- rbind(sd.atl, rep(0, nrow(lr.sar.indIf)), sd.pac)
colnames(sd.val) <- atl.final$names
sd.val <- cbind(sd.val, Ocean2 = rep(0, 3))
sd.val <- sd.val[, order(colnames(sd.val))]
sd.val <- rbind(sd.rsr[order(rsr.sample$names)], sd.val)

png("Figures/Ana_3xi_LRatio_OceR_2.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0, atl.final, lr.sar.indIf, pac.final, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.val)
dev.off()

sd.val.g <- rbind(sd.atl.g, rep(0, nrow(lr.sar.indIfg)), sd.pac.g)
colnames(sd.val.g) <- atl.final.g$names
sd.val.g <- cbind(sd.val.g, Ocean = rep(0, 3))
sd.val.g <- sd.val.g[, order(colnames(sd.val.g))]
sd.val.g <- rbind(sd.rsr.g[order(rsr.sample.g$names)], sd.val.g)

png("Figures/Ana_3xi_LRatio_g_OceR_sd_2.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, atl.final.g, lr.sar.indIfg, pac.final.g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.val.g)
dev.off()

rm(pac.data)

## 3xii. Tidy up the results -----------------------------------------------
save(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, atl.w, file = "Outputs/Atlantic_simplified2.RData")
save(lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, ind.w, file = "Outputs/Indian_simplified2.RData")
save(lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, pac.w, file = "Outputs/Pacific_simplified2.RData")

rm(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, op.formula, ind.w, atl.w, pac.w)


## 4. Does resolution of variables matter? ---------------------------------
# Don't do this


## 5. Does delta_carb_ion cut-off matter for rarefied richness? ------------------------------

## 5i. Create a dataset for modelling with higher or lower cut-offs --------
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling
rm(db.traits, margo.traits, margo.traits.cons)

# create a dataset for modelling with 
tmp <- c("Core", "Latitude", "Longitude", "Water.Depth", "Total_Planktics", "Ocean2", "sp.rich", "rarefy.sr", "simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
ldg.margo.mod <- merge(ldg.margo.env, ldg.margo.data[, tmp], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env, tmp)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]
# remove the extra factor level of the mediterranean
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
## which points are to be excluded because of delta_carb_ion?
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -15), ]
ldg.margo.mod$delta_carb_ion[ldg.margo.mod$delta_carb_ion > 0] <- 0
# salinity
ldg.margo.mod$absMnSal.0m <- abs(ldg.margo.mod$meanSal.0m - 35.1)
# wrong taxonomy
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$cons_tax == "Y", ]

## for richness
cols <- c("simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun")
rsr.margo.mod15 <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
rsr.margo.mod15 <- na.omit(rsr.margo.mod15)

rm(ldg.margo.mod)

# remove unique data
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK")
dim(unique(rsr.margo.mod15[, cols]))
tmp.1 <- which(duplicated(rsr.margo.mod15[, cols]))
rsr.margo.mod15$uni <- NA
rsr.margo.mod15$uni[!duplicated(rsr.margo.mod15[, cols])] <- 1:length(rsr.margo.mod15$uni[!duplicated(rsr.margo.mod15[, cols])])
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(rsr.margo.mod15[i, ], rsr.margo.mod15, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  rsr.margo.mod15$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows, cols, tmp.1)
rsr.margo.mod15$uni <- factor(rsr.margo.mod15$uni)
rsr.margo.dup15 <- rsr.margo.mod15
# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
rsr.margo.mod15 <- rsr.margo.dup15[!duplicated(rsr.margo.dup15$uni), ]
for (i in 1:ncol(rsr.margo.mod15)) {
  if (!is.factor(rsr.margo.mod15[,i]) & !is.character(rsr.margo.mod15[, i])) {
    rsr.margo.mod15[, i] <- as.numeric(tapply(rsr.margo.dup15[, i], rsr.margo.dup15$uni, mean, na.rm = TRUE))
  }
}
rm(i)

## Also a cut-off of -5
rsr.margo.mod5 <- rsr.margo.mod[which(rsr.margo.mod$delta_carb_ion >= -5), ]

## compare these different cutoffs
# -5
with(rsr.margo.mod5, distrib.map(Longitude, Latitude, Ocean2))
table(rsr.margo.mod5$Ocean2)

# -10.9
with(rsr.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(rsr.margo.mod$Ocean2)

# -15
with(rsr.margo.mod15, distrib.map(Longitude, Latitude, Ocean2))
table(rsr.margo.mod15$Ocean2)

rm(rsr.margo.dup15)

## 5ii. Run the model with these different cutoffs -------------------------
# -5
ldg.coords.5 <- cbind(rsr.margo.mod5$Long,rsr.margo.mod5$Lat)
ldg.coords.5 <- as.matrix(ldg.coords.5)
op5.nb <- dnearneigh(ldg.coords.5, 0, mod.sar.opW$dist, longlat = TRUE)
op5.w <- nb2listw(op5.nb, glist = NULL, style = "W", zero.policy = TRUE)

mod.dis5.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op5.w, zero.policy = TRUE, tol.solve = 1e-18, data = rsr.margo.mod5)

# -15
ldg.coords.15 <- cbind(rsr.margo.mod15$Long,rsr.margo.mod15$Lat)
ldg.coords.15 <- as.matrix(ldg.coords.15)
op15.nb <- dnearneigh(ldg.coords.15, 0, mod.sar.opW$dist, longlat = TRUE)
op15.w <- nb2listw(op15.nb, glist = NULL, style = "W", zero.policy = TRUE)

mod.dis15.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op15.w, zero.policy = TRUE, tol.solve = 1e-18, data = rsr.margo.mod15)

## 5iii. Compare LRs -------------------------------------------------------
lr.dis5.op0 <- lr.calc(mod.dis5.op0)
lr.dis15.op0 <- lr.calc(mod.dis15.op0)

png("Figures/Ana_5iii_LRatio_opdis5op0_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.dis5.op0, lr.sar.op0, lr.dis15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.dis5.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.dis5.op0g <- lr.calc(mod.dis5.op0, tmp)
lr.dis15.op0g <- lr.calc(mod.dis15.op0, tmp)
rm(tmp)

png("Figures/Ana_5iii_LRatio_g_opdis5op0_2.png", width = 800)
lr.plot(lr.dis5.op0g, lr.sar.op0g, lr.dis15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

save(rsr.margo.mod5, mod.dis5.op0, lr.dis5.op0, lr.dis5.op0g, op5.w, ldg.coords.5, file = "Outputs/mod_dis5_2.RData")
save(rsr.margo.mod15, mod.dis15.op0, lr.dis15.op0, lr.dis15.op0g, op15.w, ldg.coords.15, file = "Outputs/mod_dis15_2.RData")
rm(rsr.margo.mod5, mod.dis5.op0, lr.dis5.op0, lr.dis5.op0g, rsr.margo.mod15, mod.dis15.op0, lr.dis15.op0, lr.dis15.op0g, ldg.coords.5, ldg.coords.15, op5.nb, op15.nb, op5.w, op15.w)


## 6. Does averaging explanatory variables matter? -------------------------
# 
#  
# ldg.coords.ran <- cbind(tmp7$Longitude, tmp7$Latitude)
# ldg.coords.ran <- as.matrix(ldg.coords.ran)
# opran.nb <- dnearneigh(ldg.coords.ran, 0, mod.sar.opW$dist, longlat = TRUE)
# opran.s <- nb2listw(opran.nb, glist = NULL, style = "W", zero.policy = TRUE)
# 
# mod.disran.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = opran.s, zero.policy = TRUE, tol.solve = 1e-18, data = tmp7)

## 6iii. Compare LRs -------------------------------------------------------
# lr.dis.ran.op0 <- lr.calc(mod.disran.op0)
# 
# png("Figures/Ana_6iii_LRatio_opranop0ii_2.png", width = 1000)
# # get order from running without order first
# lr.plot(lr.dis.ran.op0, lr.sar.op0, order = c(7:9, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Random", "Full"), ylab = "Log Likelihood ratio", star.pos = 20)
# dev.off()
# 
# # also for groups of variables
# (tmp <- data.frame(names = model.evs(mod.disran.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
# lr.dis.ran.op0g <- lr.calc(mod.disran.op0, tmp)
# rm(tmp)
# 
# png("Figures/Ana_6iii_LRatio_g_opranop0ii_2.png", width = 800)
# lr.plot(lr.dis.ran.op0g, lr.sar.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Random", "Full"), ylab = "Log Likelihood ratio", star.pos = 20)
# dev.off()
# 
# save(tmp7, mod.ran.op0, lr.ran.op0, lr.ran.op0g, file = "Outputs/mod_ran.RData")
# rm(rsr.margo.modran, mod.ran.op0, lr.ran.op0, lr.ran.op0g, ldg.coords.ran, opran.nb, opran.s)


## 7. Does the number of points in the dataset make a difference? ----------------------


## 8. Evenness -------------------------------------------------------------

## 8i. Create an OLS model ------------------------------------------------
mod.eve.l0 <- lm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, data = eve.margo.mod)

# check model plots
png("Figures/Ana_8i_modevel0_2.png", 600, 600)
par(mfrow = c(2, 2))
plot(mod.eve.l0)
par(mfrow = c(1, 1))
dev.off()

summary(mod.eve.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
mod.eve.l0.sac <- with(eve.margo.mod, spline.correlog(Longitude, Latitude, mod.eve.l0$residuals, latlon = TRUE, resamp = 1))
summary(mod.eve.l0.sac)
png("Figures/Ana_8i_modevel0SAC_2.png")
plot.spline.correlog.n(mod.eve.l0.sac, xlab = "Distance / km")
dev.off()

## 8ii. Run model optimisation ----------------------------------------------
ldg.coords.eve <- cbind(eve.margo.mod$Longitude, eve.margo.mod$Latitude)
ldg.coords.eve <- as.matrix(ldg.coords.eve)

mod.sar.eveW <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveW$obj, Nagelkerke = TRUE) # 0.44411
AIC(mod.sar.eveW$obj) # -2576.519

# check other coding styles
mod.sar.eveB <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveB$obj, Nagelkerke = TRUE) # 0.44539
AIC(mod.sar.eveB$obj) # 2579.608

mod.sar.eveS <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveS$obj, Nagelkerke = TRUE) # 0.45847
AIC(mod.sar.eveS$obj) # -2611.709

mod.sar.eveS
# So "S" is the best coding style and the best neighbourhood distance is 538.1047

eve.nb <- dnearneigh(ldg.coords.eve, 0, mod.sar.eveS$dist, longlat = TRUE)
eve.s <- nb2listw(eve.nb, glist = NULL, style = "S", zero.policy = TRUE)
mod.sar.eve0 <- errorsarlm(mod.sar.eveS$mod, listw = eve.s, zero.policy = TRUE, tol.solve = 1e-18)
rm(eve.nb)

save(mod.sar.eveB, mod.sar.eveS, mod.sar.eveW, file = "Outputs/Evenness_coding2.RData")
rm(mod.sar.eveB, mod.sar.eveW)

## 8iii. Create a plot of model parameters --------------------------------------
summary(mod.sar.eve0, Nagelkerke = T) # r2 = 0.45847

# generate a dataframe of coefficients
m.eve.coef <- data.frame(names = names(mod.sar.eve0$coefficients), coef.sar = mod.sar.eve0$coefficients, row.names = 1:length(mod.sar.eve0$coefficients), stars = NA)

# add a column of significance stars
for (i in 1:length(stars)) {
  m.eve.coef$stars[which(summary(mod.sar.eve0)$Coef[, 4] <= stars[i] & is.na(m.eve.coef$stars))] <- names(stars)[i]
}
rm(i)

# plot the absolute coefficients
png("Figures/Ana_8iii_coef_modsareve0_2.png", width = 1000, height = 750)
plt.def <- par("plt")
par(plt = c(plt.def[1:2], 0.5, plt.def[4]))
tmp.x <- barplot(abs(m.eve.coef$coef.sar), names = m.eve.coef$names, las = 2, ylim = c(0, max(m.eve.coef$coef.sar) + 50))
text(tmp.x, abs(m.eve.coef$coef.sar) + 10, m.eve.coef$stars)
par(plt = plt.def)
dev.off()
rm(plt.def, tmp.x, m.eve.coef)

## 8iv. Calculate likelihood ratios for the SAR model ----------------------
# full model is
summary(mod.sar.eve0, Nagelkerke = T) # r2 = 0.45847 
AIC(mod.sar.eve0) # -2611.709

lr.sar.eve0 <- lr.calc(mod.sar.eve0)

png("Figures/Ana_8iv_LRatio_eve0_2.png", width = 800)
lr.plot(lr.sar.eve0, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

## 8v. Calculate likelihood ratios for groups of EVs ---------------------

# for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.eve0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.eve0g <- lr.calc(mod.sar.eve0, tmp)
rm(tmp)

png("Figures/Ana_8v_LRatio_g_eve0_2.png", width = 800)
lr.plot(lr.sar.eve0g, order = c(6:7, 1:5), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

## 8vi. Effect of delta_carb_ion ----------------------------------------------
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling
rm(db.traits, margo.traits)

# create a dataset for modelling with 
tmp <- c("Core", "Latitude", "Longitude", "Water.Depth", "Total_Planktics", "Ocean2", "sp.rich", "rarefy.sr", "simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
ldg.margo.mod <- merge(ldg.margo.env, ldg.margo.data[, tmp], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env, tmp)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]
# remove the extra factor level of the mediterranean
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
## which points are to be excluded because of delta_carb_ion?
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -15), ]
ldg.margo.mod$delta_carb_ion[ldg.margo.mod$delta_carb_ion > 0] <- 0
# salinity
ldg.margo.mod$absMnSal.0m <- abs(ldg.margo.mod$meanSal.0m - 35.1)
# wrong taxonomy
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$cons_tax == "Y", ]

## for eveness
cols <- c("sp.rich", "rarefy.sr", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun")
eve.margo.mod15 <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
eve.margo.mod15 <- na.omit(eve.margo.mod15)

rm(ldg.margo.mod)

# remove unique data
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK")
dim(unique(eve.margo.mod15[, cols]))
tmp.1 <- which(duplicated(eve.margo.mod15[, cols]))
eve.margo.mod15$uni <- NA
eve.margo.mod15$uni[!duplicated(eve.margo.mod15[, cols])] <- 1:length(eve.margo.mod15$uni[!duplicated(eve.margo.mod15[, cols])])
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(eve.margo.mod15[i, ], eve.margo.mod15, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  eve.margo.mod15$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows, cols, tmp.1)
eve.margo.mod15$uni <- factor(eve.margo.mod15$uni)
eve.margo.dup15 <- eve.margo.mod15
# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
eve.margo.mod15 <- eve.margo.dup15[!duplicated(eve.margo.dup15$uni), ]
for (i in 1:ncol(eve.margo.mod15)) {
  if (!is.factor(eve.margo.mod15[,i]) & !is.character(eve.margo.mod15[, i])) {
    eve.margo.mod15[, i] <- as.numeric(tapply(eve.margo.dup15[, i], eve.margo.dup15$uni, mean, na.rm = TRUE))
  }
}
rm(i)

## Also a cut-off of -5
eve.margo.mod5 <- eve.margo.mod[which(eve.margo.mod$delta_carb_ion >= -5), ]

## compare these different cutoffs
# -5
with(eve.margo.mod5, distrib.map(Longitude, Latitude, Ocean2))
table(eve.margo.mod5$Ocean2)

# -10.9
with(eve.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(eve.margo.mod$Ocean2)

# -15
with(eve.margo.mod15, distrib.map(Longitude, Latitude, Ocean2))
table(eve.margo.mod15$Ocean2)

rm(eve.margo.dup15)

## Run the model with these different cutoffs
# -5
ldg.coords.eve.5 <- cbind(eve.margo.mod5$Long, eve.margo.mod5$Lat)
ldg.coords.eve.5 <- as.matrix(ldg.coords.eve.5)
op.eve5.nb <- dnearneigh(ldg.coords.eve.5, 0, mod.sar.eveS$dist, longlat = TRUE)
op.eve5.S <- nb2listw(op.eve5.nb, glist = NULL, style = "S", zero.policy = TRUE)

mod.eve5.op0 <- errorsarlm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.eve5.S, zero.policy = TRUE, tol.solve = 1e-18, data = eve.margo.mod5)

# -15
ldg.coords.eve.15 <- cbind(eve.margo.mod15$Long,eve.margo.mod15$Lat)
ldg.coords.eve.15 <- as.matrix(ldg.coords.eve.15)
op.eve15.nb <- dnearneigh(ldg.coords.eve.15, 0, mod.sar.eveS$dist, longlat = TRUE)
op.eve15.S <- nb2listw(op.eve15.nb, glist = NULL, style = "S", zero.policy = TRUE)

mod.eve15.op0 <- errorsarlm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.eve15.S, zero.policy = TRUE, tol.solve = 1e-18, data = eve.margo.mod15)

lr.eve5.op0 <- lr.calc(mod.eve5.op0)
lr.eve15.op0 <- lr.calc(mod.eve15.op0)

png("Figures/Ana_8vi_LRatio_op_eve5op0_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.eve5.op0, lr.sar.eve0, lr.eve15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.eve5.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.eve5.op0g <- lr.calc(mod.eve5.op0, tmp)
lr.eve15.op0g <- lr.calc(mod.eve15.op0, tmp)
rm(tmp)

png("Figures/Ana_8vi_LRatio_g_op_eve5op0_2.png", width = 800)
lr.plot(lr.eve5.op0g, lr.sar.eve0g, lr.eve15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

save(eve.margo.mod5, mod.eve5.op0, lr.eve5.op0, lr.eve5.op0g, op.eve5.S, file = "Outputs/mod_eve5_2.RData")
save(eve.margo.mod15, mod.eve15.op0, lr.eve15.op0, lr.eve15.op0g, op.eve15.S, file = "Outputs/mod_eve15_2.RData")
rm(mod.eve5.op0, lr.eve5.op0, lr.eve5.op0g, mod.eve15.op0, lr.eve15.op0, lr.eve15.op0g, ldg.coords.eve.5, ldg.coords.eve.15, op.eve5.nb, op.eve15.nb, op.eve5.S, op.eve15.S)


## 9. Lineage age -------------------------------------------------------------

## 9i. Lineage age models -------------------------------------------------
mod.sar.lnaW <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, LinAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaW$obj, Nagelkerke = TRUE) # 0.53665
AIC(mod.sar.lnaW$obj) # 6006.896

# check other coding styles
mod.sar.lnaB <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, LinAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaB$obj, Nagelkerke = TRUE) # 0.50455 
AIC(mod.sar.lnaB$obj) #  6096.969

mod.sar.lnaS <- with(eve.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, LinAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.eve, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaS$obj, Nagelkerke = TRUE) # 0.53137 
AIC(mod.sar.lnaS$obj) # 6022.125

mod.sar.lnaW
# So "W" is the best coding style and the best neighbourhood distance is 534.7033

lna.nb <- dnearneigh(ldg.coords.eve, 0, mod.sar.lnaW$dist, longlat = TRUE)
lna.w <- nb2listw(lna.nb, glist = NULL, style = "W", zero.policy = TRUE)
mod.sar.lna0 <- errorsarlm(mod.sar.lnaW$mod, listw = lna.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mod.sar.lna0, Nagelkerke = TRUE) # 0.53665 

save(mod.sar.lnaB, mod.sar.lnaS, mod.sar.lnaW, file = "Outputs/Lineage_coding_2.RData")
rm(lna.nb, mod.sar.lnaB, mod.sar.lnaS)

lr.sar.lna0 <- lr.calc(mod.sar.lna0)

png("Figures/Ana_9i_LRatio_lna0_2.png", width = 800)
lr.plot(lr.sar.lna0, ylab = "Log Likelihood ratio", order = c(8:12, 1, 4, 2:3, 5:7), star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.lna0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.lna0g <- lr.calc(mod.sar.lna0, tmp)
rm(tmp)

png("Figures/Ana_9i_LRatio_g_lna0_2.png", width = 800)
lr.plot(lr.sar.lna0g, order = c(6:7, 1:5), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

## 9ii. Dissolution cutoffs ------------------------------------------------
# datasets are the same as for evenness

## Run the model with these different cutoffs
# -5
ldg.coords.lna.5 <- cbind(eve.margo.mod5$Long, eve.margo.mod5$Lat)
ldg.coords.lna.5 <- as.matrix(ldg.coords.lna.5)
op.lna5.nb <- dnearneigh(ldg.coords.lna.5, 0, mod.sar.lnaW$dist, longlat = TRUE)
op.lna5.W <- nb2listw(op.lna5.nb, glist = NULL, style = "W", zero.policy = TRUE)

mod.lna5.op0 <- errorsarlm(LinAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.lna5.W, zero.policy = TRUE, tol.solve = 1e-18, data = eve.margo.mod5)

# -15
ldg.coords.lna.15 <- cbind(eve.margo.mod15$Long,eve.margo.mod15$Lat)
ldg.coords.lna.15 <- as.matrix(ldg.coords.lna.15)
op.lna15.nb <- dnearneigh(ldg.coords.lna.15, 0, mod.sar.lnaW$dist, longlat = TRUE)
op.lna15.W <- nb2listw(op.lna15.nb, glist = NULL, style = "W", zero.policy = TRUE)

mod.lna15.op0 <- errorsarlm(LinAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.lna15.W, zero.policy = TRUE, tol.solve = 1e-18, data = eve.margo.mod15)

lr.lna5.op0 <- lr.calc(mod.lna5.op0)
lr.lna15.op0 <- lr.calc(mod.lna15.op0)

png("Figures/Ana_9ii_LRatio_op_lna5op0_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.lna5.op0, lr.sar.lna0, lr.lna15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.lna5.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.lna5.op0g <- lr.calc(mod.lna5.op0, tmp)
lr.lna15.op0g <- lr.calc(mod.lna15.op0, tmp)
rm(tmp)

png("Figures/Ana_9ii_LRatio_g_op_lna5op0_2.png", width = 800)
lr.plot(lr.lna5.op0g, lr.sar.lna0g, lr.lna15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

save(eve.margo.mod5, mod.lna5.op0, lr.lna5.op0, lr.lna5.op0g, op.lna5.W, file = "Outputs/mod_lna5_2.RData")
save(eve.margo.mod15, mod.lna15.op0, lr.lna15.op0, lr.lna15.op0g, op.lna15.W, file = "Outputs/mod_lna15_2.RData")
rm(eve.margo.mod5, mod.lna5.op0, lr.lna5.op0, lr.lna5.op0g, eve.margo.mod15, mod.lna15.op0, lr.lna15.op0, lr.lna15.op0g, ldg.coords.lna.5, ldg.coords.lna.15, op.lna5.nb, op.lna15.nb, op.lna5.W, op.lna15.W)


## 10. Functional richness -----------------------------------------------

## 10i. Create an OLS model ------------------------------------------------
mod.fric.l0 <- lm(FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, data = fric.margo.mod)

# check model plots
png("Figures/Ana_10i_modfricl0_2.png", 600, 600)
par(mfrow = c(2, 2))
plot(mod.fric.l0)
par(mfrow = c(1, 1))
dev.off()

summary(mod.fric.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
mod.fric.l0.sac <- with(fric.margo.mod, spline.correlog(Longitude, Latitude, mod.fric.l0$residuals, latlon = TRUE, resamp = 1))
summary(mod.fric.l0.sac)
png("Figures/Ana_10i_modfricl0SAC_2.png")
plot.spline.correlog.n(mod.fric.l0.sac, xlab = "Distance / km")
dev.off()

## 10ii. Run model optimisation --------------------------------------------
ldg.coords.fric <- cbind(fric.margo.mod$Long, fric.margo.mod$Lat)
ldg.coords.fric <- as.matrix(ldg.coords.fric)

mod.sar.fricW <- with(fric.margo.mod, sar.optimised(mod.fric.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.fric, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricW$obj, Nagelkerke = TRUE) # 0.84902
AIC(mod.sar.fricW$obj) # -1584.785

# check other coding styles
mod.sar.fricB <- with(fric.margo.mod, sar.optimised(mod.fric.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.fric, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricB$obj, Nagelkerke = TRUE) #  0.84943 
AIC(mod.sar.fricB$obj) # -1588.43

mod.sar.fricS <- with(fric.margo.mod, sar.optimised(mod.fric.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, ldg.coords.fric, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricS$obj, Nagelkerke = TRUE) #  0.85458
AIC(mod.sar.fricS$obj) # -1634.781

mod.sar.fricS
# So "S" is the best coding style and the best neighbourhood distance is 576.8785

fric.nb <- dnearneigh(ldg.coords.fric, 0, mod.sar.fricS$dist, longlat = TRUE)
fric.s <- nb2listw(fric.nb, glist = NULL, style = "S", zero.policy = TRUE)
mod.sar.fric0 <- errorsarlm(mod.sar.fricS$mod, listw = fric.s, zero.policy = TRUE, tol.solve = 1e-18)
summary(mod.sar.fric0, Nagelkerke = TRUE) # 0.85458

save(mod.sar.fricB, mod.sar.fricS, mod.sar.fricW, file = "Outputs/Lineage_coding_2.RData")
rm(mod.sar.fricB, mod.sar.fricW)

## 10iii. Create model plots -----------------------------------------------
lr.sar.fric0 <- lr.calc(mod.sar.fric0)

png("Figures/Ana_10iii_LRatio_fric0_2.png", width = 800)
lr.plot(lr.sar.fric0, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.fric0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.fric0g <- lr.calc(mod.sar.fric0, tmp)
rm(tmp)

png("Figures/Ana_10iii_LRatio_g_fric0_2.png", width = 800)
lr.plot(lr.sar.fric0g, order = c(6:7, 1:5), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

## 10iv. Test for affect of dissolution ------------------------------------
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling
rm(db.traits, margo.traits)

# create a dataset for modelling with 
tmp <- c("Core", "Latitude", "Longitude", "Water.Depth", "Total_Planktics", "Ocean2", "sp.rich", "rarefy.sr", "simpson", "simpsonEve", "FRic", "symbionts_obl", "symbionts_obl_abun", "symbionts_all", "symbionts_all_abun", "surface", "surface_subsurface", "subsurface", "subsurface_deep", "deep", "surfaceAbun", "surface_subsurfaceAbun", "subsurfaceAbun", "subsurface_deepAbun", "deepAbun", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun", "cons_tax")
ldg.margo.mod <- merge(ldg.margo.env, ldg.margo.data[, tmp], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env, tmp)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]
# remove the extra factor level of the mediterranean
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
## which points are to be excluded because of delta_carb_ion?
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -15), ]
ldg.margo.mod$delta_carb_ion[ldg.margo.mod$delta_carb_ion > 0] <- 0
# salinity
ldg.margo.mod$absMnSal.0m <- abs(ldg.margo.mod$meanSal.0m - 35.1)
# wrong taxonomy
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$cons_tax == "Y", ]

## for FRic
cols <- c("sp.rich", "rarefy.sr", "simpson", "simpsonEve", "MorphoAge", "LinAge", "MorphoAgeAbun", "LinAgeAbun")
fric.margo.mod15 <- ldg.margo.mod[, !(names(ldg.margo.mod) %in% cols)]
rm(cols)
fric.margo.mod15 <- na.omit(fric.margo.mod15)

rm(ldg.margo.mod)

# remove unique data
cols <- c("meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.mld.t", "sd.mld.t", "mean.mld.d", "sd.mld.d", "mean.mld.v", "sd.mld.v", "mld.exact", "depth10deg", "meanSal.0m", "sdSal.0m", "sal.exact", "meanOxy", "sdOxy", "prop2.oxy", "oxy.exact", "delta_carb_ion", "delta_carb_ion.OK")
dim(unique(fric.margo.mod15[, cols]))
tmp.1 <- which(duplicated(fric.margo.mod15[, cols]))
fric.margo.mod15$uni <- NA
fric.margo.mod15$uni[!duplicated(fric.margo.mod15[, cols])] <- 1:length(fric.margo.mod15$uni[!duplicated(fric.margo.mod15[, cols])])
for (i in tmp.1) {
  # identify the matching rows (based on the relevant columns)
  match.rows <- merge(fric.margo.mod15[i, ], fric.margo.mod15, by.x = cols, by.y = cols)
  # extract the value for uni for these rows add that value to the duplicated row
  fric.margo.mod15$uni[i] <- unique(match.rows$uni.y)[!is.na(unique(match.rows$uni.y))]
}
rm(i, match.rows, cols, tmp.1)
fric.margo.mod15$uni <- factor(fric.margo.mod15$uni)
fric.margo.dup15 <- fric.margo.mod15
# having got a column that identifies each unique grid cell, now need to create a dataframe that contains the mean value for each of these
fric.margo.mod15 <- fric.margo.dup15[!duplicated(fric.margo.dup15$uni), ]
for (i in 1:ncol(fric.margo.mod15)) {
  if (!is.factor(fric.margo.mod15[,i]) & !is.character(fric.margo.mod15[, i])) {
    fric.margo.mod15[, i] <- as.numeric(tapply(fric.margo.dup15[, i], fric.margo.dup15$uni, mean, na.rm = TRUE))
  }
}
rm(i)

## Also a cut-off of -5
fric.margo.mod5 <- fric.margo.mod[which(fric.margo.mod$delta_carb_ion >= -5), ]

## compare these different cutoffs
# -5
with(fric.margo.mod5, distrib.map(Longitude, Latitude, Ocean2))
table(fric.margo.mod5$Ocean2)

# -10.9
with(fric.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(fric.margo.mod$Ocean2)

# -15
with(fric.margo.mod15, distrib.map(Longitude, Latitude, Ocean2))
table(fric.margo.mod15$Ocean2)

rm(fric.margo.dup15)

## Run the model with these different cutoffs
# -5
ldg.coords.fric.5 <- cbind(fric.margo.mod5$Long, fric.margo.mod5$Lat)
ldg.coords.fric.5 <- as.matrix(ldg.coords.fric.5)
op.fric5.nb <- dnearneigh(ldg.coords.fric.5, 0, mod.sar.fricS$dist, longlat = TRUE)
op.fric5.S <- nb2listw(op.fric5.nb, glist = NULL, style = "S", zero.policy = TRUE)

mod.fric5.op0 <- errorsarlm(FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.fric5.S, zero.policy = TRUE, tol.solve = 1e-18, data = fric.margo.mod5)

# -15
ldg.coords.fric.15 <- cbind(fric.margo.mod15$Long,fric.margo.mod15$Lat)
ldg.coords.fric.15 <- as.matrix(ldg.coords.fric.15)
op.fric15.nb <- dnearneigh(ldg.coords.fric.15, 0, mod.sar.fricS$dist, longlat = TRUE)
op.fric15.S <- nb2listw(op.fric15.nb, glist = NULL, style = "S", zero.policy = TRUE)

mod.fric15.op0 <- errorsarlm(FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2)^2 + delta_carb_ion, listw = op.fric15.S, zero.policy = TRUE, tol.solve = 1e-18, data = fric.margo.mod15)

lr.fric5.op0 <- lr.calc(mod.fric5.op0)
lr.fric15.op0 <- lr.calc(mod.fric15.op0)

png("Figures/Ana_10iv_LRatio_op_fric5op0_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.fric5.op0, lr.sar.fric0, lr.fric15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.fric5.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.fric5.op0g <- lr.calc(mod.fric5.op0, tmp)
lr.fric15.op0g <- lr.calc(mod.fric15.op0, tmp)
rm(tmp)

png("Figures/Ana_10iv_LRatio_g_op_fric5op0_2.png", width = 800)
lr.plot(lr.fric5.op0g, lr.sar.fric0g, lr.fric15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5", "Cutoff: -10.9", "Cutoff: -15"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

save(fric.margo.mod5, mod.fric5.op0, lr.fric5.op0, lr.fric5.op0g, op.fric5.S, file = "Outputs/mod_fric5_2.RData")
save(fric.margo.mod15, mod.fric15.op0, lr.fric15.op0, lr.fric15.op0g, op.fric15.S, file = "Outputs/mod_fric15_2.RData")
rm(mod.fric5.op0, lr.fric5.op0, lr.fric5.op0g, mod.fric15.op0, lr.fric15.op0, lr.fric15.op0g, ldg.coords.fric.5, ldg.coords.fric.15, op.fric5.nb, op.fric15.nb, op.fric5.S, op.fric15.S)


## 11. Paper outputs -------------------------------------------------------

## 11i. Table of full model cofficients for the models -----------
# for rarefied species richness
# create table
op0.summary <- as.data.frame(summary(mod.sar.op0, Nagelkerke = T)$Coef)

# create a function for rounding
paper.round <- function(x) {
  ifelse(abs(x) > 1, round(x, 2), signif(x, 3)) 
}

# apply it to all the columns
for (i in 1:ncol(op0.summary)) {
  op0.summary[, i] <- paper.round(op0.summary[, i])
}
rm(i)

# add significance stars
stars <- c(0.001, 0.01, 0.05, 0.1)
names(stars) <- c("***", "**", "*", ".")
op0.summary$Significance <- NA
for (i in 1:length(stars)) {
  op0.summary$Significance[which(op0.summary[,"Pr(>|z|)"] <= stars[i] & is.na(op0.summary$Significance))] <- names(stars)[i]
}
rm(i)

# write to csv
write.csv(op0.summary, "Outputs/op0_summary.csv")
rm(op0.summary)

# for evenness
eve0.summary <- as.data.frame(summary(mod.sar.eve0, Nagelkerke = T)$Coef)

# apply it to all the columns
for (i in 1:ncol(eve0.summary)) {
  eve0.summary[, i] <- paper.round(eve0.summary[, i])
}
rm(i)

# add significance stars
stars <- c(0.001, 0.01, 0.05, 0.1)
names(stars) <- c("***", "**", "*", ".")
eve0.summary$Significance <- NA
for (i in 1:length(stars)) {
  eve0.summary$Significance[which(eve0.summary[,"Pr(>|z|)"] <= stars[i] & is.na(eve0.summary$Significance))] <- names(stars)[i]
}
rm(i)

# write to csv
write.csv(eve0.summary, "Outputs/eve0_summary.csv")
rm(eve0.summary)

# for community age
lna0.summary <- as.data.frame(summary(mod.sar.lna0, Nagelkerke = T)$Coef)

# apply it to all the columns
for (i in 1:ncol(lna0.summary)) {
  lna0.summary[, i] <- paper.round(lna0.summary[, i])
}
rm(i)

# add significance stars
stars <- c(0.001, 0.01, 0.05, 0.1)
names(stars) <- c("***", "**", "*", ".")
lna0.summary$Significance <- NA
for (i in 1:length(stars)) {
  lna0.summary$Significance[which(lna0.summary[,"Pr(>|z|)"] <= stars[i] & is.na(lna0.summary$Significance))] <- names(stars)[i]
}
rm(i)

# write to csv
write.csv(lna0.summary, "Outputs/lna0_summary.csv")
rm(lna0.summary)

# for FRic
fric0.summary <- as.data.frame(summary(mod.sar.fric0, Nagelkerke = T)$Coef)

# apply it to all the columns
for (i in 1:ncol(fric0.summary)) {
  fric0.summary[, i] <- paper.round(fric0.summary[, i])
}
rm(i)

# add significance stars
stars <- c(0.001, 0.01, 0.05, 0.1)
names(stars) <- c("***", "**", "*", ".")
fric0.summary$Significance <- NA
for (i in 1:length(stars)) {
  fric0.summary$Significance[which(fric0.summary[,"Pr(>|z|)"] <= stars[i] & is.na(fric0.summary$Significance))] <- names(stars)[i]
}
rm(i)

# write to csv
write.csv(fric0.summary, "Outputs/fric0_summary.csv")
rm(fric0.summary)

## 11ii. likelihood ratios for the models ----------------------------------------
(lr.out <- cbind(lr.sar.op0, lr.sar.eve0[2:4], lr.sar.lna0[2:4], lr.sar.fric0[2:4]))
names(lr.out)[2:4] <- paste("sr", names(lr.out)[2:4], sep = "_")
names(lr.out)[5:7] <- paste("eve", names(lr.out)[5:7], sep = "_")
names(lr.out)[8:10] <- paste("lna", names(lr.out)[8:10], sep = "_")

write.csv(lr.out, "Outputs/lr_out.csv")
rm(lr.out)

## 11iii. log likelihood ratio plot for full vs simplified model rarified SR --------
png("Figures/Ana_11iii_LRatio_op0f_2.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(6:7, 4:3, 5, 2:1), ylab = "Log Likelihood ratio", leg.txt = c("Full", "Simplified"), cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

## 11iv. log likelihood ratio plot comparing different delta_carb_ion cut-off --------
load("Outputs/mod_dis5_2.RData")
load("Outputs/mod_dis15_2.RData")
load("Outputs/mod_eve5_2.RData")
load("Outputs/mod_eve15_2.RData")
load("Outputs/mod_lna5_2.RData")
load("Outputs/mod_lna15_2.RData")
load("Outputs/mod_fric5_2.RData")
load("Outputs/mod_fric15_2.RData")

# for rsr
png("Figures/Ana_11iv_LRatio_rsr_5op15_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.dis5.op0, lr.sar.op0, lr.dis15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

png("Figures/Ana_11iv_LRatio_g_rsr_5op15_2.png", width = 800)
lr.plot(lr.dis5.op0g, lr.sar.op0g, lr.dis15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

# for evenness
png("Figures/Ana_11iv_LRatio_eve_5op15_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.eve5.op0, lr.sar.eve0, lr.eve15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

png("Figures/Ana_11iv_LRatio_g_eve_5op15_2.png", width = 800)
lr.plot(lr.eve5.op0g, lr.sar.eve0g, lr.eve15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

# for lineage age
png("Figures/Ana_11iv_LRatio_lna_5op15_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.lna5.op0, lr.sar.lna0, lr.lna15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

png("Figures/Ana_11iv_LRatio_g_lna_5op15_2.png", width = 800)
lr.plot(lr.lna5.op0g, lr.sar.lna0g, lr.lna15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

# for FRic
png("Figures/Ana_11iv_LRatio_fric_5op15_2.png", width = 1000)
# get order from running without order first
lr.plot(lr.fric5.op0, lr.sar.fric0, lr.fric15.op0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

png("Figures/Ana_11iv_LRatio_g_fric_5op15_2.png", width = 800)
lr.plot(lr.fric5.op0g, lr.sar.fric0g, lr.fric15.op0g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Cutoff: -5.0", "Cutoff: -10.9", "Cutoff: -15.0"), ylab = "Log Likelihood ratio", star.pos = 20, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

rm(rsr.margo.mod5, mod.dis5.op0, lr.dis5.op0, lr.dis5.op0g, rsr.margo.mod15, mod.dis15.op0, lr.dis15.op0, lr.dis15.op0g, eve.margo.mod5, mod.eve5.op0, lr.eve5.op0, lr.eve5.op0g, eve.margo.mod15, mod.eve15.op0, lr.eve15.op0, lr.eve15.op0g, mod.lna5.op0, lr.lna5.op0, lr.lna5.op0g, mod.lna15.op0, lr.lna15.op0, lr.lna15.op0g, fric.margo.mod5, mod.fric5.op0, lr.fric5.op0, lr.fric5.op0g, fric.margo.mod15, mod.fric15.op0, lr.fric15.op0, lr.fric15.op0g)

## 11v. residuals map for full model --------------------------------------------
png("Figures/Ana_11v_rsp_res_2.png", 700, 500)
with(rsr.margo.mod, distrib.map(Longitude, Latitude, mod.sar.op0$residuals, palette = "rwbt", col.water = "white", col.land = "black", maintitle = "Residuals for rarefied richness", min.col = -8, max.col = 8))
dev.off()

png("Figures/Ana_11v_eve_res_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, mod.sar.eve0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -0.7, max.col = 0.7, maintitle = "Residuals for evenness"))
dev.off()

png("Figures/Ana_11v_lna_res_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, mod.sar.lna0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -11, max.col = 11, maintitle = "Residuals for average community age"))
dev.off()

png("Figures/Ana_11v_fric_res_2.png", 700, 500)
with(fric.margo.mod, distrib.map(Longitude, Latitude, mod.sar.fric0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -0.8, max.col = 0.8, maintitle = "Residuals for functional richness"))
dev.off()

## 11vi. coplot ------------------------------------------------------------------

# plot marginal effects (hold all else constant) of coefficients of T polynomial in oceans. Plot across temperature range observed in the oceans (3 colours / plotting symbols)

## 11vii. comparing all the models ------------------------------------------------
png("Figures/Ana_11vii_LRatio_RELF_2.png", width = 800)
lr.plot(lr.sar.op0, lr.sar.eve0, lr.sar.lna0, lr.sar.fric0, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age", "Functional richness"), cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

png("Figures/Ana_11vii_LRatio_g_RELF_2.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.eve0g, lr.sar.lna0g, lr.sar.fric0g, order = c(6:7, 4:3, 5, 2:1), ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age", "Functional richness"), cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

## 11viii. compare this to a random model ----------------------------------------
# rdm.vals <- data.frame(model = rep(NA, 1000 * 7), names = NA, lr = NA, p = NA, stars = NA)
# tmp <- data.frame(names = model.evs(mod.rdm), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature"))
# 
# for(i in 1:1000) {
#   print(i)
#   mod.rdm <- errorsarlm(rnorm(nrow(ldg.margo.mod)) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2)^2 + delta_carb_ion, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)
#   tmp.lr <- lr.calc(mod.rdm, tmp)
#   rdm.vals[1:7 + (i - 1) * 7, 1] <- rep(i, 7)
#   rdm.vals[1:7 + (i - 1) * 7, 2] <- as.character(tmp.lr[, 1])
#   rdm.vals[1:7 + (i - 1) * 7, 3:5] <- tmp.lr[, 2:4]
# }
# rm(i)
# 
# rdm.mn <- with(rdm.vals, tapply(lr, names, mean))
# rdm.sd <- with(rdm.vals, tapply(lr, names, sd))
# tmp.x <- barplot(rdm.mn)


## 12. Predicting models ---------------------------------------------------
# at this point, need to run Prediction_environment.R
# load("Outputs/ldg_p_margo.RData")

## 12i. Set up the dataset for prediction ----------------------------------
with(rsr.margo.mod, distrib.map(Longitude, Latitude, sar.predict(mod.sar.op0)))

# check NAs and remove them
summary(ldg.p.margo)
par(ask = TRUE)
for (i in 3:ncol(ldg.p.margo))
{
  with(ldg.p.margo[is.na(ldg.p.margo[,i]), ], distrib.map(Longitude, Latitude, Longitude, main = names(ldg.p.margo)[i]))
}
par(ask = FALSE)
rm(i)

ldg.p.margo <- na.omit(ldg.p.margo)

# don't predict outside the range of the margo data
# if I don't do this, get predictions of -500 species in the north of Russia
par(ask = TRUE)
for (i in 4:(ncol(ldg.p.margo) - 1)) {
  tmp.col <- which(names(eve.margo.mod) == names(ldg.p.margo)[i])
  ldg.p.margo[which(ldg.p.margo[, i] < min(eve.margo.mod[, tmp.col], na.rm = TRUE) | ldg.p.margo[, i] > max(eve.margo.mod[, tmp.col], na.rm = TRUE)), i] <- NA
  if (sum(is.na(ldg.p.margo[, i])) > 0)
    with(ldg.p.margo[is.na(ldg.p.margo[, i]), ], distrib.map(Longitude, Latitude, Ocean2, pch = 15, cex = 0.5, main = names(ldg.p.margo)[i]))
}
par(ask = FALSE)
# n.b. get warnings where there are no points outside the range
rm(i, tmp.col)

# simplify ldg.p.margo down to only contain those variables used in the analysis
names(ldg.p.margo)
ldg.p.margo <- ldg.p.margo[, names(ldg.p.margo) %in% c("Longitude", "Latitude", env.var)]
ldg.p.margo <- na.omit(ldg.p.margo)

# plot these up
for (i in env.var) {
  if (i != "delta_carb_ion" & i != "Ocean2")  # n.b. ignore these as neither plot well
  {
    png(paste("Figures/Ana_12i_map_", i, "_2.png", sep = ""), width = 800, height = 450)
    with(ldg.p.margo, distrib.map(Longitude, Latitude, ldg.p.margo[, i], palette = "matlab.like", pch = 15, cex = 0.5, main = i, col.land = "black", col.water = "white"))
    dev.off()
  }
}
rm(i)

## 12ii. predict species richness for this dataset -------------------------------
# necessary to use sar.predict not predict as the ordinary predict.sarlm function cannot cope with poly variables
ldg.p.margo$rarefy.sr <- sar.predict(mod.sar.op0, newdata = ldg.p.margo, olddata = rsr.margo.mod)
summary(ldg.p.margo$rarefy.sr)

with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with negative richness - where are these
with(ldg.p.margo[ldg.p.margo$rarefy.sr <= 0, ], distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, cex = 0.4)) # around coastlines

# compare with observed
png("Figures/Ana_12ii_rsr_pred_2.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$rarefy.sr > 0 & ldg.p.margo$rarefy.sr <= 27, ], distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", maintitle = "Rarefied species richness"))
dev.off()

png("Figures/Ana_12ii_rsr_obs_2.png", 700, 500)
with(rsr.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", min.col = 0, max.col = 27, col.water = "white", col.land = "black"))
dev.off()

## 12iii. predict evenness for this dataset ---------------------------------------
ldg.p.margo$simpsonEve <- sar.predict(mod.sar.eve0, newdata = ldg.p.margo, olddata = eve.margo.mod)
summary(ldg.p.margo$simpsonEve)

with(ldg.p.margo, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4)) 

# compare with observed
png("Figures/Ana_12iii_eve_pred_2.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$simpsonEve <= 1, ], distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", maintitle = "Simpson's evenness"))
dev.off()

png("Figures/Ana_12iii_eve_obs_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 12iv. predict average community age for this dataset --------------------
ldg.p.margo$LinAgeAbun <- sar.predict(mod.sar.lna0, newdata = ldg.p.margo, olddata = eve.margo.mod)
summary(ldg.p.margo$LinAgeAbun)

with(ldg.p.margo, distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with age outside sensible limits - where are these
plot(sort(ldg.p.margo$LinAgeAbun))
summary(eve.margo.mod$LinAgeAbun)
with(ldg.p.margo[ldg.p.margo$LinAgeAbun < 5, ], distrib.map(Longitude, Latitude, LinAgeAbun, pch = 15, cex = 0.4)) # around coastlines - very few points

# compare with observed
png("Figures/Ana_12iv_lna_pred_2.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$LinAgeAbun >= 5, ], distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", maintitle = "Average community age"))
dev.off()

png("Figures/Ana_12iv_lna_obs_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 12v. predict functional richness for this dataset --------------------
ldg.p.margo$FRic <- sar.predict(mod.sar.fric0, newdata = ldg.p.margo, olddata = fric.margo.mod)
summary(ldg.p.margo$FRic)
summary(fric.margo.mod$FRic)

with(ldg.p.margo, distrib.map(Longitude, Latitude, FRic, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with FRic outside sensible limits - where are these
with(ldg.p.margo[ldg.p.margo$FRic < 0, ], distrib.map(Longitude, Latitude, FRic, pch = 15, cex = 0.4)) # around coastlines

# compare with observed
png("Figures/Ana_12v_fric_pred_2.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$FRic >= 0, ], distrib.map(Longitude, Latitude, FRic, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", maintitle = "Functional richness"))
dev.off()

png("Figures/Ana_12v_fric_obs_2.png", 700, 500)
with(fric.margo.mod, distrib.map(Longitude, Latitude, FRic, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 12vi. Difference maps ----------------------------------------------------
# set up centred and scaled
rsr.margo.mod$rarefy.sr.cs <- (rsr.margo.mod$rarefy.sr - mean(rsr.margo.mod$rarefy.sr)) / sd(rsr.margo.mod$rarefy.sr)
eve.margo.mod$simpsonEve.cs <- (eve.margo.mod$simpsonEve - mean(eve.margo.mod$simpsonEve)) / sd(eve.margo.mod$simpsonEve)
eve.margo.mod$LinAgeAbun.cs <- (eve.margo.mod$LinAgeAbun - mean(eve.margo.mod$LinAgeAbun)) / sd(eve.margo.mod$LinAgeAbun)
fric.margo.mod$FRic.cs <- (fric.margo.mod$FRic - mean(fric.margo.mod$FRic)) / sd(fric.margo.mod$FRic)

ldg.p.margo$rarefy.sr.cs <- (ldg.p.margo$rarefy.sr - mean(ldg.p.margo$rarefy.sr)) / sd(ldg.p.margo$rarefy.sr)
ldg.p.margo$simpsonEve.cs <- (ldg.p.margo$simpsonEve - mean(ldg.p.margo$simpsonEve)) / sd(ldg.p.margo$simpsonEve)
ldg.p.margo$LinAgeAbun.cs <- (ldg.p.margo$LinAgeAbun - mean(ldg.p.margo$LinAgeAbun)) / sd(ldg.p.margo$LinAgeAbun)
ldg.p.margo$FRic.cs <- (ldg.p.margo$FRic - mean(ldg.p.margo$FRic)) / sd(ldg.p.margo$FRic)

# can't do most of the comparisons with the raw data, as the datasets are different sizes

## differences
# b/w eve and lna
# red - high eve & low lna
# blue - low eve & high lna
png("Figures/Ana_12vi_eve_lna_2.png", 700, 500)
with(eve.margo.mod, distrib.map(Longitude, Latitude, simpsonEve.cs - LinAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", min.col = -8, max.col = 8))
dev.off()

# predicted
# b/w sr and eve
# red - high sr & low eve
# blue - low sr & high eve
png("Figures/Ana_12vi_sr_eve_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr.cs - simpsonEve.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, min.col = -8, max.col = 8, cex = 0.4))
dev.off()

# b/w sr and lna
# red - high sr & low lna
# blue - low sr & high lna
png("Figures/Ana_12vi_sr_lna_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr.cs - LinAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -6, max.col = 6))
dev.off()

# b/w sr and fric
# red - high sr & low fric
# blue - low sr & high fric
png("Figures/Ana_12vi_sr_fric_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr.cs - FRic.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -5, max.col = 5))
dev.off()

# b/w eve and lna
# red - high eve & low lna
# blue - low eve & high lna
png("Figures/Ana_12vi_eve_lna_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, simpsonEve.cs - LinAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -10, max.col = 10))
dev.off()

# b/w eve and fric
# red - high eve & low fric
# blue - low eve & high fric
png("Figures/Ana_12vi_eve_fric_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, simpsonEve.cs - FRic.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -8, max.col = 8))
dev.off()

# b/w lna and fric
# red - high lna & low fric
# blue - low lna & high fric
png("Figures/Ana_12vi_lna_fric_pred_2.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, LinAgeAbun.cs - FRic.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -6, max.col = 6))
dev.off()

## 12vii. Comparison plots --------------------------------------------------
# Don't do these for grouped models, as there are multiple possible variables for the horizontal axis

# for full rarefied model
lr.calc(mod.sar.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = rsr.margo.mod, file.nm = "FullRSR/op0")

# for simplified rarefied model
lr.calc(mod.sar.opf, plots = TRUE, pred.data = ldg.p.margo, mod.data = rsr.margo.mod, file.nm = "SimpRSR/opf")


# for Oceans
load("Outputs/Atlantic_simplified2.RData")
load("Outputs/Indian_simplified2.RData")
load("Outputs/Pacific_simplified2.RData")

lr.sar.atlIf <- lr.calc(mod.sar.atlIf, plots = TRUE, pred.data = ldg.p.margo[ldg.p.margo$Ocean2 == "Atlantic",], mod.data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", ], file.nm = "Atlantic/atlIf")

lr.sar.indIf <- lr.calc(mod.sar.indIf, plots = TRUE, pred.data = ldg.p.margo[ldg.p.margo$Ocean2 == "Indian",], mod.data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", ], file.nm = "Indian/indIf")

lr.sar.pacIf <- lr.calc(mod.sar.pacIf, plots = TRUE, pred.data = ldg.p.margo[ldg.p.margo$Ocean2 == "Pacific",], mod.data = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", ], file.nm = "Pacific/pacIf")

rm(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, ind.w, atl.w, pac.w)


# for richness dissolution levels
load("Outputs/mod_dis5_2.RData")
load("Outputs/mod_dis15_2.RData")

lr.calc(mod.dis5.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = rsr.margo.mod5, file.nm = "FullRSR/op05")

lr.calc(mod.dis15.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = rsr.margo.mod15, file.nm = "FullRSR/op015")

rm(rsr.margo.mod5, mod.dis5.op0, lr.dis5.op0, lr.dis5.op0g, rsr.margo.mod15, mod.dis15.op0, lr.dis15.op0, lr.dis15.op0g, ldg.coords.5, ldg.coords.15, op5.w, op15.w)


# for evenness
lr.calc(mod.sar.eve0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod, file.nm = "Eve/eve")

# evenness dissolution
load("Outputs/mod_eve5_2.RData")
load("Outputs/mod_eve15_2.RData")

lr.calc(mod.eve5.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod5, file.nm = "Eve/op05")

lr.calc(mod.eve15.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod15, file.nm = "Eve/op015")

rm(mod.eve5.op0, lr.eve5.op0, lr.eve5.op0g, mod.eve15.op0, lr.eve15.op0, lr.eve15.op0g, op.eve5.S, op.eve15.S)


# for lineage age
lr.calc(mod.sar.lna0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod, file.nm = "Lna/lna")

# lineage age dissolution
load("Outputs/mod_lna5_2.RData")
load("Outputs/mod_lna15_2.RData")

lr.calc(mod.lna5.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod5, file.nm = "Lna/op05")

lr.calc(mod.lna15.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = eve.margo.mod15, file.nm = "Lna/op015")

rm(mod.lna5.op0, lr.lna5.op0, lr.lna5.op0g, mod.lna15.op0, lr.lna15.op0, lr.lna15.op0g, op.lna5.W, op.lna15.W)


# for FRic
lr.calc(mod.sar.fric0, plots = TRUE, pred.data = ldg.p.margo, mod.data = fric.margo.mod, file.nm = "FRic/fric")

# functional richness dissolution
load("Outputs/mod_fric5_2.RData")
load("Outputs/mod_fric15_2.RData")

lr.calc(mod.fric5.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = fric.margo.mod5, file.nm = "FRic/op05")

lr.calc(mod.fric15.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = fric.margo.mod15, file.nm = "FRic/op015")

rm(mod.fric5.op0, lr.fric5.op0, lr.fric5.op0g, mod.fric15.op0, lr.fric15.op0, lr.fric15.op0g, op.fric5.S, op.fric15.S)

## 12viii. RMSE ------------------------------------------------------------
rmse <- function(obs, pred)
{
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

## calculate predicted values 
# for richness
rsr.margo.mod$rarefy.sr.op0 <- sar.predict(mod.sar.op0, newdata = rsr.margo.mod, olddata = rsr.margo.mod)
rsr.margo.mod$rarefy.sr.opf <- sar.predict(mod.sar.opf, newdata = rsr.margo.mod, olddata = rsr.margo.mod)

load("Outputs/Atlantic_simplified2.RData")
load("Outputs/Indian_simplified2.RData")
load("Outputs/Pacific_simplified2.RData")

rsr.margo.mod$rarefy.sr.atl <- NA
rsr.margo.mod$rarefy.sr.atl[rsr.margo.mod$Ocean2 == "Atlantic"] <- sar.predict(mod.sar.atlIf, newdata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", ], olddata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Atlantic", ])

rsr.margo.mod$rarefy.sr.ind <- NA
rsr.margo.mod$rarefy.sr.ind[rsr.margo.mod$Ocean2 == "Indian"] <- sar.predict(mod.sar.indIf, newdata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", ], olddata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Indian", ])

rsr.margo.mod$rarefy.sr.pac <- NA
rsr.margo.mod$rarefy.sr.pac[rsr.margo.mod$Ocean2 == "Pacific"] <- sar.predict(mod.sar.pacIf, newdata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", ], olddata = rsr.margo.mod[rsr.margo.mod$Ocean2 == "Pacific", ])

rm(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, ind.w, atl.w, pac.w)

load("Outputs/mod_dis5_2.RData")
load("Outputs/mod_dis15_2.RData")

rsr.margo.mod5$rarefy.sr.dis5 <- sar.predict(mod.dis5.op0, newdata = rsr.margo.mod5, olddata = rsr.margo.mod5)
rsr.margo.mod15$rarefy.sr.dis15 <- sar.predict(mod.dis15.op0, newdata = rsr.margo.mod15, olddata = rsr.margo.mod15)

rm(mod.dis5.op0, lr.dis5.op0, lr.dis5.op0g, mod.dis15.op0, lr.dis15.op0, lr.dis15.op0g, ldg.coords.5, ldg.coords.15, op5.w, op15.w)

# for evenness
eve.margo.mod$simpsonEve.op0 <- sar.predict(mod.sar.eve0, newdata = eve.margo.mod, olddata = eve.margo.mod)

load("Outputs/mod_eve5_2.RData")
load("Outputs/mod_eve15_2.RData")
load("Outputs/mod_lna5_2.RData")
load("Outputs/mod_lna15_2.RData")
eve.margo.mod5$simpsonEve.dis5 <- sar.predict(mod.eve5.op0, newdata = eve.margo.mod5, olddata = eve.margo.mod5)
eve.margo.mod15$simpsonEve.dis15 <- sar.predict(mod.eve15.op0, newdata = eve.margo.mod15, olddata = eve.margo.mod15)
rm(mod.eve5.op0, lr.eve5.op0, lr.eve5.op0g, mod.eve15.op0, lr.eve15.op0, lr.eve15.op0g, op.eve5.S, op.eve15.S)

# for lineage age
eve.margo.mod$LinAgeAbun.op0 <- sar.predict(mod.sar.lna0, newdata = eve.margo.mod, olddata = eve.margo.mod)
eve.margo.mod5$LinAgeAbun.dis5 <- sar.predict(mod.lna5.op0, newdata = eve.margo.mod5, olddata = eve.margo.mod5)
eve.margo.mod15$LinAgeAbun.dis15 <- sar.predict(mod.lna15.op0, newdata = eve.margo.mod15, olddata = eve.margo.mod15)
rm(mod.lna5.op0, lr.lna5.op0, lr.lna5.op0g, mod.lna15.op0, lr.lna15.op0, lr.lna15.op0g, op.lna5.W, op.lna15.W)

# for functional richness
fric.margo.mod$FRic.op0 <- sar.predict(mod.sar.fric0, newdata = fric.margo.mod, olddata = fric.margo.mod)
load("Outputs/mod_fric5_2.RData")
load("Outputs/mod_fric15_2.RData")
fric.margo.mod5$FRic.dis5 <- sar.predict(mod.fric5.op0, newdata = fric.margo.mod5, olddata = fric.margo.mod5)
fric.margo.mod15$FRic.dis15 <- sar.predict(mod.fric15.op0, newdata = fric.margo.mod15, olddata = fric.margo.mod15)
rm(mod.fric5.op0, lr.fric5.op0, lr.fric5.op0g, mod.fric15.op0, lr.fric15.op0, lr.fric15.op0g, op.fric5.S, op.fric15.S)

save(rsr.margo.mod, rsr.margo.mod5, rsr.margo.mod15, eve.margo.mod, eve.margo.mod5, eve.margo.mod15, fric.margo.mod, fric.margo.mod5, fric.margo.mod15, file = "Outputs/margo_mod_predict_2.RData")

## calculate rmse for the different models
# rsr
with(rsr.margo.mod, rmse(rarefy.sr, rarefy.sr.op0)) # 1.7428
with(rsr.margo.mod, rmse(rarefy.sr, rarefy.sr.opf)) # 1.769167
with(rsr.margo.mod, rmse(rarefy.sr, rarefy.sr.atl)) # 1.591344
with(rsr.margo.mod, rmse(rarefy.sr, rarefy.sr.ind)) # 1.380638
with(rsr.margo.mod, rmse(rarefy.sr, rarefy.sr.pac)) # 1.82437
with(rsr.margo.mod5, rmse(rarefy.sr, rarefy.sr.dis5)) # 1.698643
with(rsr.margo.mod15, rmse(rarefy.sr, rarefy.sr.dis15)) # 1.764601

# evenness
with(eve.margo.mod, rmse(simpsonEve, simpsonEve.op0)) # 0.09586449
with(eve.margo.mod5, rmse(simpsonEve, simpsonEve.dis5)) # 0.09346534
with(eve.margo.mod15, rmse(simpsonEve, simpsonEve.dis15)) # 0.09761757

# lineage age
with(eve.margo.mod, rmse(LinAgeAbun, LinAgeAbun.op0)) #  2.349567
with(eve.margo.mod5, rmse(LinAgeAbun, LinAgeAbun.dis5)) # 2.368997
with(eve.margo.mod15, rmse(LinAgeAbun, LinAgeAbun.dis15)) # 2.368435

# functional richness
with(fric.margo.mod, rmse(FRic, FRic.op0)) #  0.1507324
with(fric.margo.mod5, rmse(FRic, FRic.dis5)) # 0.1470422
with(fric.margo.mod15, rmse(FRic, FRic.dis15)) # 0.1511181


## 13. Adding error bars ---------------------------------------------------
duplicated.random <- function(x, incomparables = FALSE, ...) 
{ 
  if ( is.vector(x) ) 
  { 
    permutation <- sample(length(x)) 
    x.perm      <- x[permutation] 
    result.perm <- duplicated(x.perm, incomparables, ...) 
    result      <- result.perm[order(permutation)] 
    return(result) 
  } 
  else if ( is.matrix(x) ) 
  { 
    permutation <- sample(nrow(x)) 
    x.perm      <- x[permutation,] 
    result.perm <- duplicated(x.perm, incomparables, ...) 
    result      <- result.perm[order(permutation)] 
    return(result) 
  } 
  else 
  { 
    stop(paste("duplicated.random() only supports vectors & matrices for now.")) 
  } 
} 

## 13i. Rarefied richness models -------------------------------------------
rsr.sample <- lr.sar.op0
rsr.sample.g <- lr.sar.op0g
names(rsr.sample)[2:4] <- names(rsr.sample.g)[2:4] <- paste(names(rsr.sample)[2:4], "all", sep = "_")

(var.dat <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))

load("Outputs/160821 rsr_margo_dup.RData")

# repeat 100 times
for (i in 1:100) {
  # create a subsample of the dataset
  tmp <- rsr.margo.dup[!duplicated.random(as.numeric(rsr.margo.dup$uni)), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long, tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.opW$dist, longlat = TRUE)
  tmp.w <- nb2listw(tmp.nb, glist = NULL, style = "W", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.opW$mod, listw = tmp.w, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  rsr.sample <- cbind(rsr.sample, tmp2[2:4])
  rsr.sample.g <- cbind(rsr.sample.g, tmp3[2:4])
  names(rsr.sample)[(3 * i + 2):(3 * i + 4)] <- names(rsr.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.w, tmp2, tmp3, tmp.mod, i)

rsr.final <- lr.sar.op0
rsr.final$lr <- rowMeans(rsr.sample[, seq(5, ncol(rsr.sample), by = 3)])
rsr.final$p <- rowMeans(rsr.sample[, seq(6, ncol(rsr.sample), by = 3)])
rsr.final$stars <- NA

for (i in 1:length(stars)) {
  rsr.final$stars[which(rsr.final$p <= stars[i] & is.na(rsr.final$stars))] <- names(stars)[i]
}

se.rsr <- apply(rsr.sample[, seq(5, ncol(rsr.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.rsr <- apply(rsr.sample[, seq(5, ncol(rsr.sample), by = 3)], 1, sd)

lr.plot(rsr.final, order = c(10:8, 11:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.rsr)

rsr.final.g <- lr.sar.op0g
rsr.final.g$lr <- rowMeans(rsr.sample.g[, seq(5, ncol(rsr.sample.g), by = 3)])
rsr.final.g$p <- rowMeans(rsr.sample.g[, seq(6, ncol(rsr.sample.g), by = 3)])
rsr.final.g$stars <- NA

for (i in 1:length(stars)) {
  rsr.final.g$stars[which(rsr.final.g$p <= stars[i] & is.na(rsr.final.g$stars))] <- names(stars)[i]
}

se.rsr.g <- apply(rsr.sample.g[, seq(5, ncol(rsr.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.rsr.g <- apply(rsr.sample.g[, seq(5, ncol(rsr.sample.g), by = 3)], 1, sd)

lr.plot(rsr.final.g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 25, se.val = se.rsr.g)

lr.plot(lr.sar.op0g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 25, se.val = se.rsr.g)


lr.plot(lr.sar.op0g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 25, se.val = sd.rsr.g)

load("Outputs/Atlantic_simplified2.RData")
load("Outputs/Indian_simplified2.RData")
load("Outputs/Pacific_simplified2.RData")

# need to re-run section from 3xi. now I have rsr.sd
# plot the whole set up
sd.val <- rbind(sd.atl, rep(0, nrow(lr.sar.indIf)), sd.pac)
colnames(sd.val) <- atl.final$names
sd.val <- cbind(sd.val, Ocean2 = rep(0, 3))
sd.val <- sd.val[, order(colnames(sd.val))]
sd.val <- rbind(sd.rsr[order(rsr.sample$names)], sd.val)

png("Figures/Ana_3xi_LRatio_OceR_2.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0, atl.final, lr.sar.indIf, pac.final, order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.val)
dev.off()

sd.val.g <- rbind(sd.atl.g, rep(0, nrow(lr.sar.indIfg)), sd.pac.g)
colnames(sd.val.g) <- atl.final.g$names
sd.val.g <- cbind(sd.val.g, Ocean = rep(0, 3))
sd.val.g <- sd.val.g[, order(colnames(sd.val.g))]
sd.val.g <- rbind(sd.rsr.g[order(rsr.sample.g$names)], sd.val.g)


png("Figures/Ana_3xi_LRatio_g_OceR_sd_2.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, atl.final.g, lr.sar.indIfg, pac.final.g, order = c(6:7, 4:3, 5, 2:1), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.val.g)
dev.off()

save(rsr.sample, rsr.sample.g, rsr.final, rsr.final.g, file = "Outputs/RSR_resampling_2.RData")

## 13ii. Evenness models ---------------------------------------------------
eve.sample <- lr.sar.eve0
eve.sample.g <- lr.sar.eve0g
names(eve.sample)[2:4] <- names(eve.sample.g)[2:4] <- paste(names(eve.sample)[2:4], "all", sep = "_")

(var.dat <- data.frame(names = model.evs(mod.sar.eve0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))

load("Outputs/160821 eve_margo_dup.RData")

# repeat 100 times
for (i in 27:100) {
  # create a subsample of the dataset
  print(i)
  tmp <- eve.margo.dup[!duplicated.random(as.numeric(eve.margo.dup$uni)), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long, tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.eveS$dist, longlat = TRUE)
  tmp.S <- nb2listw(tmp.nb, glist = NULL, style = "S", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.eveS$mod, listw = tmp.S, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  eve.sample <- cbind(eve.sample, tmp2[2:4])
  eve.sample.g <- cbind(eve.sample.g, tmp3[2:4])
  names(eve.sample)[(3 * i + 2):(3 * i + 4)] <- names(eve.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.S, tmp2, tmp3, tmp.mod, i)

eve.final <- lr.sar.eve0
eve.final$lr <- rowMeans(eve.sample[, seq(5, ncol(eve.sample), by = 3)])
eve.final$p <- rowMeans(eve.sample[, seq(6, ncol(eve.sample), by = 3)])
eve.final$stars <- NA

for (i in 1:length(stars)) {
  eve.final$stars[which(eve.final$p <= stars[i] & is.na(eve.final$stars))] <- names(stars)[i]
}

se.eve <- apply(eve.sample[, seq(5, ncol(eve.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.eve <- apply(eve.sample[, seq(5, ncol(eve.sample), by = 3)], 1, sd)

lr.plot(eve.final, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.eve)

lr.plot(lr.sar.eve0, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, se.val = sd.eve)

eve.final.g <- lr.sar.eve0g
eve.final.g$lr <- rowMeans(eve.sample.g[, seq(5, ncol(eve.sample.g), by = 3)])
eve.final.g$p <- rowMeans(eve.sample.g[, seq(6, ncol(eve.sample.g), by = 3)])
eve.final.g$stars <- NA

for (i in 1:length(stars)) {
  eve.final.g$stars[which(eve.final.g$p <= stars[i] & is.na(eve.final.g$stars))] <- names(stars)[i]
}

se.eve.g <- apply(eve.sample.g[, seq(5, ncol(eve.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.eve.g <- apply(eve.sample.g[, seq(5, ncol(eve.sample.g), by = 3)], 1, sd)

lr.plot(eve.final.g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 10, se.val = sd.eve.g)

lr.plot(lr.sar.eve0g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 10, se.val = sd.eve.g)

save(eve.sample, eve.sample.g, eve.final, eve.final.g, file = "Outputs/Eve_resampling_2.RData")

## 13iii. Lineage age models -----------------------------------------------
lna.sample <- lr.sar.lna0
lna.sample.g <- lr.sar.lna0g
names(lna.sample)[2:4] <- names(lna.sample.g)[2:4] <- paste(names(lna.sample)[2:4], "all", sep = "_")

(var.dat <- data.frame(names = model.evs(mod.sar.lna0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))

# repeat 100 times
for (i in 1:100) {
  # create a subsample of the dataset
  tmp <- eve.margo.dup[!duplicated.random(as.numeric(eve.margo.dup$uni)), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long, tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.lnaW$dist, longlat = TRUE)
  tmp.W <- nb2listw(tmp.nb, glist = NULL, style = "W", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.lnaW$mod, listw = tmp.W, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  lna.sample <- cbind(lna.sample, tmp2[2:4])
  lna.sample.g <- cbind(lna.sample.g, tmp3[2:4])
  names(lna.sample)[(3 * i + 2):(3 * i + 4)] <- names(lna.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.S, tmp2, tmp3, tmp.mod, i)

lna.final <- lr.sar.lna0
lna.final$lr <- rowMeans(lna.sample[, seq(5, ncol(lna.sample), by = 3)])
lna.final$p <- rowMeans(lna.sample[, seq(6, ncol(lna.sample), by = 3)])
lna.final$stars <- NA

for (i in 1:length(stars)) {
  lna.final$stars[which(lna.final$p <= stars[i] & is.na(lna.final$stars))] <- names(stars)[i]
}

se.lna <- apply(lna.sample[, seq(5, ncol(lna.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.lna <- apply(lna.sample[, seq(5, ncol(lna.sample), by = 3)], 1, sd)

lr.plot(lna.final, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, se.val = se.lna)

lna.final.g <- lr.sar.lna0g
lna.final.g$lr <- rowMeans(lna.sample.g[, seq(5, ncol(lna.sample.g), by = 3)])
lna.final.g$p <- rowMeans(lna.sample.g[, seq(6, ncol(lna.sample.g), by = 3)])
lna.final.g$stars <- NA

for (i in 1:length(stars)) {
  lna.final.g$stars[which(lna.final.g$p <= stars[i] & is.na(lna.final.g$stars))] <- names(stars)[i]
}

se.lna.g <- apply(lna.sample.g[, seq(5, ncol(lna.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.lna.g <- apply(lna.sample.g[, seq(5, ncol(lna.sample.g), by = 3)], 1, sd)

lr.plot(lna.final.g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 10, se.val = se.lna.g)

save(lna.sample, lna.sample.g, lna.final, lna.final.g, file = "Outputs/lna_resampling_2.RData")

## 13iv. Functional richness models ----------------------------------------
fric.sample <- lr.sar.fric0
fric.sample.g <- lr.sar.fric0g
names(fric.sample)[2:4] <- names(fric.sample.g)[2:4] <- paste(names(fric.sample)[2:4], "all", sep = "_")

(var.dat <- data.frame(names = model.evs(mod.sar.fric0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))

load("Outputs/160821 fric_margo_dup.RData")

# repeat 100 times
for (i in 1:100) {
  # create a subsample of the dataset
  tmp <- fric.margo.dup[!duplicated.random(as.numeric(fric.margo.dup$uni)), ]
  # run the model & calculate the lr values
  ldg.coords.tmp <- as.matrix(cbind(tmp$Long, tmp$Lat))
  tmp.nb <- dnearneigh(ldg.coords.tmp, 0, mod.sar.fricS$dist, longlat = TRUE)
  tmp.S <- nb2listw(tmp.nb, glist = NULL, style = "S", zero.policy = TRUE)
  tmp.mod <- errorsarlm(mod.sar.fricS$mod, listw = tmp.S, zero.policy = TRUE, tol.solve = 1e-18, data = tmp)
  tmp2 <- lr.calc(tmp.mod)
  tmp3 <- lr.calc(tmp.mod, var.dat)
  # add to the table
  fric.sample <- cbind(fric.sample, tmp2[2:4])
  fric.sample.g <- cbind(fric.sample.g, tmp3[2:4])
  names(fric.sample)[(3 * i + 2):(3 * i + 4)] <- names(fric.sample.g)[(3 * i + 2):(3 * i + 4)] <- paste(names(tmp2)[2:4], i, sep = "_")
}
rm(tmp, ldg.coords.tmp, tmp.nb, tmp.S, tmp2, tmp3, tmp.mod, i)

fric.final <- lr.sar.fric0
fric.final$lr <- rowMeans(fric.sample[, seq(5, ncol(fric.sample), by = 3)])
fric.final$p <- rowMeans(fric.sample[, seq(6, ncol(fric.sample), by = 3)])
fric.final$stars <- NA

for (i in 1:length(stars)) {
  fric.final$stars[which(fric.final$p <= stars[i] & is.na(fric.final$stars))] <- names(stars)[i]
}

se.fric <- apply(fric.sample[, seq(5, ncol(fric.sample), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.fric <- apply(fric.sample[, seq(5, ncol(fric.sample), by = 3)], 1, sd)

lr.plot(fric.final, order = c(8:12, 1, 4, 2:3, 5:7), ylab = "Log Likelihood ratio", star.pos = 10, se.val = se.fric)

fric.final.g <- lr.sar.fric0g
fric.final.g$lr <- rowMeans(fric.sample.g[, seq(5, ncol(fric.sample.g), by = 3)])
fric.final.g$p <- rowMeans(fric.sample.g[, seq(6, ncol(fric.sample.g), by = 3)])
fric.final.g$stars <- NA

for (i in 1:length(stars)) {
  fric.final.g$stars[which(fric.final.g$p <= stars[i] & is.na(fric.final.g$stars))] <- names(stars)[i]
}

se.fric.g <- apply(fric.sample.g[, seq(5, ncol(fric.sample.g), by = 3)], 1, function (x) sqrt(var(x) / length(x)))

sd.fric.g <- apply(fric.sample.g[, seq(5, ncol(fric.sample.g), by = 3)], 1, sd)


lr.plot(fric.final.g, ylab = "Log Likelihood ratio", order = c(6:7, 1:5), star.pos = 10, se.val = se.fric.g)

save(fric.sample, fric.sample.g, fric.final, fric.final.g, file = "Outputs/FRic_resampling_2.RData")

sd.val <- rbind(sd.rsr, sd.eve, sd.lna, sd.fric)
colnames(sd.val) <- rsr.final$names

png("Figures/Ana_13iv_LRatio_all_sd_2.png", 1200, 500)
lr.plot(rsr.final, eve.final, lna.final, fric.final, ylab = "Log Likelihood ratio", order = c(9:7, 4:3, 12:11, 5, 1, 10, 6, 2), star.pos = 20, se.val = sd.val, leg.txt = c("Species richness", "Evenness", "Lineage Age", "Functional richness"))
dev.off()


sd.val.g <- rbind(sd.rsr.g, sd.eve.g, sd.lna.g, sd.fric.g)
colnames(sd.val.g) <- rsr.final.g$names

png("Figures/Ana_13iv_LRatio_g_all_sd_2.png", 800, 500)
lr.plot(rsr.final.g, eve.final.g, lna.final.g, fric.final.g, ylab = "Log Likelihood ratio", order = c(6:7, 4:3, 5, 2:1), star.pos = 20, se.val = sd.val.g, leg.txt = c("Species richness", "Evenness", "Lineage Age", "Functional richness"), srt = 45)
dev.off()

png("Figures/Ana_13iv_LRatio_g_lna_sd_2.png", width = 6, height = 6, units = 'in', res = 300)
lr.plot(lna.final.g, ylab = "Log Likelihood ratio", star.pos = 20, order = c(6, 1, 3, 2, 7, 4, 5),  se.val = sd.val.g[3, ], cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, legend = FALSE, col = "navyblue", las = 3)
dev.off()



## 14. Metabolic theory of ecology -----------------------------------------
# log transformed SR
Ln_SR <- log(rsr.margo.mod$rarefy.sr)

# 1 / kT (boltzman constant in eV K-1) * absolute T
MTE_SST <- 1 / (8.6173324E−5 * (rsr.margo.mod$meanSST.1deg + 273.15))

# calculate the relationship
mte.mod <- errorsarlm(Ln_SR ~ MTE_SST, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mte.mod)

# get predictions
p.MTE_SST <- seq(min(MTE_SST), max(MTE_SST), length.out = 100)
p.Ln_SR <- predict(mte.mod, data.frame(MTE_SST = p.MTE_SST))[, 1]

# what would MTE predict?
mn.MTE_SST <- mean(MTE_SST)
mn.Ln_SR <- mean(Ln_SR)
# mte.MTE_SST == p.MTE_SST
mte.Ln_SR <- -0.65 * (p.MTE_SST - mn.MTE_SST) + mn.Ln_SR

# what if we model it separately for the three oceans
MTE_oce <- rsr.margo.mod$Ocean2
mte.oce.mod <- errorsarlm(Ln_SR ~ MTE_SST*MTE_oce, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mte.oce.mod)
anova(mte.oce.mod, mte.mod) # significant difference, implying that dropping ocean produces a significantly worse model
p.ocean <- data.frame(MTE_SST = rep(p.MTE_SST, 3), MTE_oce = rep(levels(MTE_oce), each = 100))
p.ocean$p.Ln_SR <- predict(mte.oce.mod, p.ocean)[, 1]

# plot them
png("Figures/Ana_14_mte_plot_2.png")
plot(Ln_SR ~ MTE_SST, pch = 16, col = rsr.margo.mod$Ocean2, xlab = "Temperature (1 / kT)", ylab = "ln (rarefied species richness)", bty = "l", las = 1, cex.lab = 1.2, cex.axis = 1.2)
points(p.MTE_SST, p.Ln_SR, type = "l", lwd = 3) # observed
with(p.ocean[p.ocean$MTE_oce == "Atlantic",], points(MTE_SST, p.Ln_SR, type = "l", col = 1)) # Atlantic
with(p.ocean[p.ocean$MTE_oce == "Indian",], points(MTE_SST, p.Ln_SR, type = "l", col = 2)) # Indian
with(p.ocean[p.ocean$MTE_oce == "Pacific",], points(MTE_SST, p.Ln_SR, type = "l", col = 3)) # Pacific
points(p.MTE_SST, mte.Ln_SR, type = "l", lty = 2, lwd = 3) # predicted
legend("topright", levels(rsr.margo.mod$Ocean2), pch = 16, col = 1:3)
dev.off()

# plot this relationship on the SST ~ SR relationship
tiff("Figures/Ana_14_sst_mte_plot_2.tif", width = 4, height = 4, units = "in")
with(rsr.margo.mod, plot(rarefy.sr ~ meanSST.1deg, pch = 16, col = c("gray50", "springgreen3", "orchid3")[rsr.margo.mod$Ocean2], xlab = expression(paste("SST / ", degree, "C")), ylab = "Rarefied species richness", bty = "l", las = 1, cex.lab = 1.2, cex.axis = 1.2))
points((1 / (8.6173324E−5 * p.MTE_SST) - 273.15), exp(mte.Ln_SR), type = "l", lwd = 2)
legend(0, 23, levels(rsr.margo.mod$Ocean2)[1:3], pch = 16, col = c("gray50", "springgreen3", "orchid3"))
dev.off()

save(Ln_SR, MTE_SST, mte.mod, mte.oce.mod, p.Ln_SR, p.MTE_SST, file = "Outputs/Metabolic_hypothesis_2.RData")
rm(Ln_SR, MTE_SST, MTE_oce, mn.Ln_SR, mn.MTE_SST, mte.Ln_SR, mte.mod, mte.oce.mod, p.Ln_SR, p.MTE_SST, p.ocean)


## 15. Speciation rates ----------------------------------------------------
# some morphospecies are in the same lineage, so need a trimmed down trait dataset that only has unique lineages

sp.traits <- unique(margo.traits[, c(21, 32:33, 35)])
sp.all.traits <- unique(pf.traits[, c(21, 25:26, 32:33)])

# 15i. Add up the total length of ALL the extant species (#1) ------------------
total.extant <- sum(sp.traits$aL.age)

# 15ii. Number of species they’ve given rise to  ---------------------------
# Tot up the number of species they’ve given rise to (#2) 
# and the number of those that are still extant (#3).
speciations <- sum(sp.traits$spec)
ext.speciations <- sum(sp.traits$extant.spec)

# 15iii. Get an average speciation rate for the extant species ---------------
# (#2) / (#1)
speciations / total.extant
# 0.073 speciations / Ma

# how does this compare with the total?
load("../../../Project/Foraminifera/Outputs/150318_pf_traits.RData")
load("../../../Project/Foraminifera/Data/2011-04-11 aL.Rdata")

# total speciation rate is 2.97 sp / Ma (across the entire clade)
length(aL$nm) / max(sp.all.traits$aL.start)

# this differs from a lineage speciation rate as if there are multiple lineages present then their individual speciation rate at any one time must sum to produce the total speciation rate
# lineage speciation rate is 0.126 sp / Ma (for an individual lineage)
sum(sp.all.traits$spec) / sum(sp.all.traits$aL.start - sp.all.traits$aL.end)

# so Recent is lower than average 

# 15iv. Get an average net diversification rate for the extant spec --------
# (#3) / (#1)
ext.speciations / total.extant
# diversification rate of 0.033 sp / Ma

# how does this compare with the total?
# total diversification rate is 0.453 sp / Ma (across the entire clade)
sum(aL$en == 0) / max(sp.all.traits$aL.start)

# this differs from a lineage diversification rate as if there are multiple lineages present then their individual diversification rate at any one time must sum to produce the total diversification rate
# lineage diversification rate is 0.018 sp / Ma (for an individual lineage)
sum(sp.all.traits$extant.spec) / sum(sp.all.traits$aL.start - sp.all.traits$aL.end)

# Recent is higher but I would expect that, as extant species are more likely to have extant offspring.

# 15v. Are subpolar assemblages particularly different? -------------------------
# Now do the same set of steps for just the species that are in, or perhaps just those that are dominant within, the subpolar assemblages.

# n.b. I got these lists from BFD/Code/Zones.R
polar <- "Neogloboquadrina pachyderma"
sapply(polar, compare)

subpolar <- c("Neogloboquadrina incompta","Globigerina quinqueloba","Globigerina bulloides","Globigerinita bradyi","Globorotalia scitula")
sapply(subpolar, compare)
subpolar <- as.character(sapply(subpolar, compare))
subpolar[1] <- "Neogloboquadrina incompta"
subpolar <- subpolar[subpolar != "Micro"]

temperate <- "Globorotalia inflata"
sapply(temperate, compare)
temperate <- as.character(sapply(temperate, compare))

subtropical <- c("Globigerinoides ruber","Globigerinoides conglobatus","Hastigerina pelagica","Globigerinita glutinata","Globorotalia truncatulinoides", "Globorotalia hirsuta","Globigerina rubescens","Globigerinella aequilateralis","Orbulina universa","Globoquadrina dutertrei","Globigerina falconensis","Globorotalia crassaformis")
sapply(subtropical, compare)
subtropical <- as.character(sapply(subtropical, compare))
subtropical <- subtropical[subtropical != "Micro"]

tropical <- c("Globigerinoides sacculifer","Sphaeroidinella dehiscens","Globorotalia menardii","Globorotalia tumida","Pulleniatina obliquiloculata","Candeina nitida", "Hastigerinella digitata","Globoquadrina conglomerata","Globigerinella adamsi","Globoquadrina hexagona")
sapply(tropical, compare)
tropical <- as.character(sapply(tropical, compare))
tropical <- tropical[tropical != "Micro"]

# create a set of traits for each zone
polar.traits <- margo.traits[rownames(margo.traits) %in% polar,]
subpolar.traits <- margo.traits[rownames(margo.traits) %in% subpolar,]
temperate.traits <- margo.traits[rownames(margo.traits) %in% temperate,]
subtropical.traits <- margo.traits[rownames(margo.traits) %in% subtropical,]
tropical.traits <- margo.traits[rownames(margo.traits) %in% tropical,]

# create a dataframe to hold this information
zone.spec <- data.frame(zone = c("polar", "subpolar", "temperate", "subtropical", "tropical", "extant", "macroperforates"), no.spec = c(length(polar), length(subpolar), length(temperate), length(subtropical), length(tropical), nrow(sp.traits), length(aL$nm)))

# column of total extant lineage length
zone.spec$total.age <- c(sum(polar.traits$aL.age), sum(subpolar.traits$aL.age), sum(temperate.traits$aL.age), sum(subtropical.traits$aL.age), sum(tropical.traits$aL.age), total.extant, sum(sp.all.traits$aL.start - sp.all.traits$aL.end))

# column of number of speciations
zone.spec$spec <- c(sum(polar.traits$spec), sum(subpolar.traits$spec), sum(temperate.traits$spec), sum(subtropical.traits$spec), sum(tropical.traits$spec), speciations, 211)

# cheated for the all species values, as there is too much duplication, which i was finding it hard to track down, so I used an excel spreadsheet to removed duplicated (I think) rows. 

# column of number of extant speciations
zone.spec$ext.spec <- c(sum(polar.traits$extant.spec), sum(subpolar.traits$extant.spec), sum(temperate.traits$extant.spec), sum(subtropical.traits$extant.spec), sum(tropical.traits$extant.spec), ext.speciations, 29)

# calculate column of speciation rates
zone.spec$spec.rate <- zone.spec$spec / zone.spec$total.age

# calculate column of extant specation rate
zone.spec$ext.spec.rate <- zone.spec$ext.spec / zone.spec$total.age

# how do speciation rates compare
zone.spec
# 0 for polar; 0.061 for subpolar; 0 for temperate; 0.081 for subtropical; 0.085 for tropical
# for extant species 0.072 speciations / Ma
# total speciation rate is 2.97 sp / Ma (across the entire clade)
# lineage speciation rate is 0.109 sp / Ma (for an individual lineage)
# rates are lower today (?are species more likely to speciate later in life)
with(zone.spec, plot(zone, spec.rate))

# what about diversification rates
zone.spec
# 0; for polar; 0.030 for subpolar; 0 for temperate; 0.37 for subtropical; 0.26 for tropical
# for extant species diversification rate of 0.032 sp / Ma
# total diversification rate is 0.453 sp / Ma (across the entire clade)
# lineage diversification rate is 0.017 sp / Ma (for an individual lineage)
with(zone.spec, plot(zone, ext.spec.rate))

# tropical might have higher speciation rate than subpolar, but it has fewer extant offspring

## 15. Tidy up -------------------------------------------------------------
save(rsr.margo.mod, eve.margo.mod, fric.margo.mod, file = "Outputs/margo_mod_2.RData")
save(lr.sar.op0, lr.sar.op0g, mod.l0.sac, mod.sar.op0, mod.sar.opW, op.w, rsr.margo.mod, file = "Outputs/Richness_model_2.RData")
save(lr.sar.opf, lr.sar.opfg, mod.sar.opf, rsr.margo.mod, op.w, file = "Outputs/Richness_model_simplified_2.RData")
save(lr.sar.eve0, lr.sar.eve0g, mod.eve.l0, mod.eve.l0.sac, mod.sar.eve0, eve.margo.mod, eve.s, file = "Outputs/Evenness_model_2.RData")
save(lr.sar.lna0, lr.sar.lna0g, mod.sar.lna0, eve.margo.mod, lna.w, file = "Outputs/Lineage_model_2.RData")
save(lr.sar.fric0, lr.sar.fric0g, mod.sar.fric0, fric.margo.mod, fric.s, file = "Outputs/FRic_model_2.RData")
save(ldg.p.margo, file = "Outputs/ldg_p_margo_pred_2.RData")

write.csv(ldg.margo.mod, file = "Outputs/ldg_margo_mod.csv")


rm(ldg.coords, stars, env.var, op.w, eve.s, lna.w, fric.s, mod.sar.eveS, mod.sar.lnaW)

