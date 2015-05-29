## Created: 15 / 4 / 2015
## Last edited: 29 / 5 / 2015
## Isabel Fenton
## Reanalysis for LDG paper
##
## Previous file: MARGO/Code/Environmental_variables.R
## Next file: 1311 LDGPaper/Reanalysis/Prediction_environment.R

## Inputs ------------------------------------------------------------------
## Environmental_variables.Rdata - containing ldg.margo.data, ldg.margo.env

## Outputs ----------------------------------------------------------------
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
setwd("C:/Documents/Science/PhD/Work/1311 LDGPaper/Reanalysis")
library(spdep) # for SAR models
library(HH) # for vif
library(ncf) # for Moran's I
library(lmtest) # for likelihood ratio tests
library(mgcv) # for gams
library(colorRamps) # for matlab.like palette
source("../Code/140420SARerrOptimising_NC.R") # for the optimising function
source("../../../Code/plot_spline_correlog_n.R") # for plotting spline correlog (with axis titles etc.)
source("../../../Code/sar_plot.R") # for producing sar plots
source("../../../Code/sar_predict.R") # for predicting with sar
source("../../../Code/maps.R") # for maps
source("../../../Code/palettes.R") # ditto
source("../../../Code/lr_calculations.R") # code for calculating likelihood ratios
source("../../../Code/sar_predict.R") # for predicting with poly
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling


## 0i. Setting up the dataset -----------------------------------------------
rm(db.traits, margo.traits)

# create a dataset for modelling with 
ldg.margo.mod <- merge(ldg.margo.env, ldg.margo.data[, c(1, 3:5, 48, 54:58, 60, 65:82)], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]
# set the FRic NAs to an arbitrary amount, as don't want to remove those rows
ldg.margo.mod$FRic[is.na(ldg.margo.mod$FRic)] <- 999
# remove other NAs
summary(ldg.margo.mod)
# most of the NAs are in rarefy.SR as it can't be calculated for data without total planktics.
ldg.margo.mod <- na.omit(ldg.margo.mod)
# reset these back to NA
ldg.margo.mod$FRic[ldg.margo.mod$FRic == 999] <- NA
  
# remove the extra factor level of the mediterranean
table(ldg.margo.mod$Ocean2)
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
table(ldg.margo.mod$Ocean2)

## which points are to be excluded because of delta_carb_ion?
with(ldg.margo.mod, plot(delta_carb_ion, rarefy.sr, pch = 16, col = Ocean2))
legend("topleft", levels(ldg.margo.mod$Ocean2), pch = 16, col = 1:3)
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -14), ]
with(ldg.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(ldg.margo.mod$Ocean2)
# based on using a 3500m / 4500m cut-off

dim(ldg.margo.mod)

## 0ii. Consider relationships ---------------------------------------------
par(ask = TRUE)
for(i in 5:30) {
  if (!is.character(ldg.margo.mod[, i]))
    with(ldg.margo.mod, plot(ldg.margo.mod[, i], rarefy.sr, pch = 16, col = Ocean2, main = names(ldg.margo.mod)[i]))
}
par(ask = FALSE)
rm(i)

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

# check whether anything else should be logged (i.e. the histogram of each environmental variable)
summary(ldg.margo.mod)
par(ask = TRUE)
for(i in c(5:30, ncol(ldg.margo.mod))) {
  if (!is.character(ldg.margo.mod[, i]))
    with(ldg.margo.mod, hist(ldg.margo.mod[, i], main = names(ldg.margo.mod)[i]))
}
par(ask = FALSE)
rm(i)
# seem reasonable

## 0iii. Create plots ------------------------------------------------------
png("Figures/Ana_0iii_map_rsr.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", main = "Rarefied species richness", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_eve.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", main = "Simpson's Evenness", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_lna.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", main = "Average Community Age", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_0iii_map_fric.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, FRic, palette = "matlab.like", main = "Functional richness", col.water = "white", col.land = "black"))
dev.off()


## 1. Check for correlation between explanatory variables ----------------------
names(ldg.margo.mod)

# variables are: mean/sd SST, mean/sd MLD, 10deg contour, mean/sd logProd, mean/sd Sal, Ocean2, carbonate ion 
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "sd.mld.t", "depth10deg", "logProd.mn.ann", "logProd.sd.ann", "absMnSal.0m", "sdSal.0m", "Ocean2", "delta_carb_ion", "meanOxy", "prop2.oxy")

# pairs plot
png("Figures/Ana_1_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var])

# currently the only ones that are potentially correlated (VIF > 5 and look correlated on the pairs plot) are mean and sd in Prod, and mean / sd mld, also meanOxy and meanSST.1deg seem to be correlated (as suggested). Look into this in a bit more detail
png("Figures/Ana_1_mnprodsdprod1deg.png")
with(ldg.margo.mod, plot(logProd.mn.ann, logProd.sd.ann, pch = 16)) # highly correlated.
dev.off()

# does the same hold for for the mld
png("Figures/Ana_1_mnMLDsdMLD.png")
with(ldg.margo.mod, plot(mean.mld.t, sd.mld.t, pch = 16)) # yes
dev.off()

png("Figures/Ana_1_mnSST1degmnOxy.png")
with(ldg.margo.mod, plot(meanSST.1deg, meanOxy, pch = 16)) # highly correlated.
dev.off()

# therefore suggest exclusion of logProd.sd.ann, sd.mld.t and meanOxy from models, so
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "Ocean2", "delta_carb_ion", "prop2.oxy")

# check this has fixed the problem
pairs(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var], )
vif(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var], )
# it has, so now have new EVs


## 2. Model simplification -------------------------------------------------

## 2i. Create an OLS model and check for SAC --------------------------------
# an OLS model
mod.l0 <- lm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, data = ldg.margo.mod)

# check model plots
png("Figures/Ana_2i_modl0.png", 600, 600)
par(mfrow = c(2, 2))
plot(mod.l0)
par(mfrow = c(1, 1))
dev.off()

# look for spatial autocorrelation in the residuals
# using spline.correlog
# haven't run this properly - should be resamp = 1000
mod.l0.sac <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.l0$residuals, latlon = TRUE, resamp = 1))
summary(mod.l0.sac)
png("Figures/Ana_2i_modl0SAC.png")
plot.spline.correlog.n(mod.l0.sac, xlab = "Distance / km")
dev.off()

# using correlog
mod.l0.SACcor <- with(ldg.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.l0), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2i_modl0SACcor.png")
plot(mod.l0.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()
rm(mod.l0.SACcor)

## look at the residuals plots
png("Figures/Ana_2i_modl0resid.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.l0$residuals, palette = "rwb"))
dev.off()

rm(mod.l0)

## 2ii. Create a GAM to check complexity / SAC -----------------------------
# n.b. GAMs can't do interactions (as additive models)
mod.g0 <- with(ldg.margo.mod, gam(rarefy.sr ~ s(Longitude, Latitude, k = 80, by = Ocean2) + s(meanSST.1deg, by = Ocean2) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2 + delta_carb_ion, gamma = 1.4))
summary(mod.g0)

png("Figures/Ana_2ii_modg0.png")
gam.check(mod.g0) # n.b. this gives an error for factors
dev.off()
par(mfrow = c(1,1))

# ## calculate SAC
# using spline.correlog
mod.g0.SAC <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.g0$residuals, latlon = TRUE, resamp = 1))
summary(mod.g0.SAC)
png("Figures/Ana_2ii_modg0SAC.png")
plot.spline.correlog.n(mod.g0.SAC, xlab = "Distance / km")
dev.off()

# using correlog
mod.g0.SACcor <- with(ldg.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.g0), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2ii_modg0SACcor.png")
plot(mod.g0.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()

rm(mod.g0, mod.g0.SAC, mod.g0.SACcor)

## 2iii. Create an optimised SARerror model --------------------------------
# Make a matrix of coordinates (X and Y coordinates)
ldg.coords <- cbind(ldg.margo.mod$Longitude, ldg.margo.mod$Latitude)
ldg.coords <- as.matrix(ldg.coords)

# run model optimisation
# getting problems with Error in solve.default(asyvar, tol = tol.solve) : 
# system is computationally singular: reciprocal condition number = 8.20242e-19 
# The suggested answer is to rescale if the variables are on very different scales, so go from
# mod.sar.opW <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
# to
summary(ldg.margo.mod) # check that ranges are roughly equivalent after scaling
mod.sar.opW <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.opW$obj, Nagelkerke = TRUE)

## check SAC has been removed
# using spline.correlog
mod.sar.opW.SAC <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.sar.opW$obj$residuals, latlon = TRUE, resamp = 1))
summary(mod.sar.opW.SAC)
png("Figures/Ana_2iii_modSarOp0SAC.png")
plot.spline.correlog.n(mod.sar.opW.SAC, xlab = "Distance / km")
dev.off()
rm(mod.sar.opW.SAC)

# using correlog
mod.sar.opW.SACcor <- with(ldg.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.sar.opW$obj), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2iii_modSarOp0SACcor.png")
plot(mod.sar.opW.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()
rm(mod.sar.opW.SACcor)

# check whether different coding methods improve the AIC
mod.sar.opB <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 7629.261
AIC(mod.sar.opB$obj) # 7678.733

mod.sar.opS <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 7629.261
AIC(mod.sar.opS$obj) # 7624.853

mod.sar.opS
# So "S" is the best coding style and the best neighbourhood distance is 501.5814
rm(mod.sar.opW, mod.sar.opB)

## 2iv. Run model simplification -------------------------------------------
summary(mod.sar.opS$obj)

# re-run this optimised model through errorsarlm, so things like anova and lr.calc work
op.nb <- dnearneigh(ldg.coords, 0, mod.sar.opS$dist, longlat = TRUE)
op.s <- nb2listw(op.nb, glist = NULL, style = "S", zero.policy = TRUE)
mod.sar.op0 <- errorsarlm(mod.sar.opS$mod, listw = op.s, zero.policy = TRUE, tol.solve = 1e-18)
rm(op.nb)

# start model simplification
summary(mod.sar.op0, Nagelkerke = TRUE) # 0.86714
summary(mod.sar.op0)$Coef[order(summary(mod.sar.op0)$Coef[, 4]),]
mod.sar.op1 <- update(mod.sar.op0, ~. -poly(meanSST.1deg, 3):sdSST.1deg + poly(meanSST.1deg, 2):sdSST.1deg)
anova(mod.sar.op0, mod.sar.op1)
summary(mod.sar.op1, Nagelkerke = TRUE) # 0.86714
AIC(mod.sar.op1) # 7622.862
summary(mod.sar.op1)$Coef[order(summary(mod.sar.op1)$Coef[, 4]),]

mod.sar.op2 <- update(mod.sar.op1, ~. -I(depth10deg/100):delta_carb_ion)
anova(mod.sar.op2, mod.sar.op1)
summary(mod.sar.op2, Nagelkerke = TRUE) #  0.86714
AIC(mod.sar.op2) # 7620.911
rm(mod.sar.op1)
summary(mod.sar.op2)$Coef[order(summary(mod.sar.op2)$Coef[, 4]),]

mod.sar.op3 <- update(mod.sar.op2, ~. -absMnSal.0m:delta_carb_ion)
anova(mod.sar.op3, mod.sar.op2)
summary(mod.sar.op3, Nagelkerke = TRUE) # 0.86714
AIC(mod.sar.op3) # 7618.946
rm(mod.sar.op2)
summary(mod.sar.op3)$Coef[order(summary(mod.sar.op3)$Coef[, 4]),]

mod.sar.op4 <- update(mod.sar.op3, ~. -poly(meanSST.1deg, 3):I(depth10deg/100) + poly(meanSST.1deg, 2):I(depth10deg/100))
anova(mod.sar.op4, mod.sar.op3)
summary(mod.sar.op4, Nagelkerke = TRUE) # 0.86713
AIC(mod.sar.op4) # 7617.004
rm(mod.sar.op3)
summary(mod.sar.op4)$Coef[order(summary(mod.sar.op4)$Coef[, 4]),]

mod.sar.op5 <- update(mod.sar.op4, ~. -I(mean.mld.t/10):I(depth10deg/100))
anova(mod.sar.op5, mod.sar.op4)
summary(mod.sar.op5, Nagelkerke = TRUE) # 0.86713
AIC(mod.sar.op5) # 7615.052
rm(mod.sar.op4)
summary(mod.sar.op5)$Coef[order(summary(mod.sar.op5)$Coef[, 4]),]

mod.sar.op6 <- update(mod.sar.op5, ~. -absMnSal.0m:sdSal.0m)
anova(mod.sar.op6, mod.sar.op5)
summary(mod.sar.op6, Nagelkerke = TRUE) # 0.86712
AIC(mod.sar.op6) # 7613.137
rm(mod.sar.op5)
summary(mod.sar.op6)$Coef[order(summary(mod.sar.op6)$Coef[, 4]),]

mod.sar.op7 <- update(mod.sar.op6, ~. -logProd.mn.ann:delta_carb_ion )
anova(mod.sar.op7, mod.sar.op6)
summary(mod.sar.op7, Nagelkerke = TRUE) # 0.86711
AIC(mod.sar.op7) # 7611.282
rm(mod.sar.op6)
summary(mod.sar.op7)$Coef[order(summary(mod.sar.op7)$Coef[, 4]),]

mod.sar.op8 <- update(mod.sar.op7, ~. -I(depth10deg/100):poly(meanSST.1deg, 2) + I(depth10deg/100):poly(meanSST.1deg, 1))
anova(mod.sar.op8, mod.sar.op7)
summary(mod.sar.op8, Nagelkerke = TRUE) # 0.8671
AIC(mod.sar.op8) # 7609.506
rm(mod.sar.op7)
summary(mod.sar.op8)$Coef[order(summary(mod.sar.op8)$Coef[, 4]),]

mod.sar.op9 <- update(mod.sar.op8, ~. -I(depth10deg/100):poly(meanSST.1deg, 1))
anova(mod.sar.op9, mod.sar.op8)
summary(mod.sar.op9, Nagelkerke = TRUE) # 0.8671
AIC(mod.sar.op9) # 7607.521
rm(mod.sar.op8)
summary(mod.sar.op9)$Coef[order(summary(mod.sar.op9)$Coef[, 4]),]

mod.sar.op10 <- update(mod.sar.op9, ~. -sdSST.1deg:I(depth10deg/100))
anova(mod.sar.op10, mod.sar.op9)
summary(mod.sar.op10, Nagelkerke = TRUE) # 0.86708
AIC(mod.sar.op10) # 7605.722
rm(mod.sar.op9)
summary(mod.sar.op10)$Coef[order(summary(mod.sar.op10)$Coef[, 4]),]

mod.sar.op11 <- update(mod.sar.op10, ~. -sdSST.1deg:poly(meanSST.1deg, 2) + sdSST.1deg:poly(meanSST.1deg, 1))
anova(mod.sar.op11, mod.sar.op10)
summary(mod.sar.op11, Nagelkerke = TRUE) # 0.86707
AIC(mod.sar.op11) # 7603.855
rm(mod.sar.op10)
summary(mod.sar.op11)$Coef[order(summary(mod.sar.op11)$Coef[, 4]),]

mod.sar.op12 <- update(mod.sar.op11, ~. -poly(meanSST.1deg, 3):I(mean.mld.t/10) + poly(meanSST.1deg, 2):I(mean.mld.t/10))
anova(mod.sar.op12, mod.sar.op11)
summary(mod.sar.op12, Nagelkerke = TRUE) # 0.86704
AIC(mod.sar.op12) # 7602.283
rm(mod.sar.op11)
summary(mod.sar.op12)$Coef[order(summary(mod.sar.op12)$Coef[, 4]),]

mod.sar.op13 <- update(mod.sar.op12, ~. -I(depth10deg/100):Ocean2)
anova(mod.sar.op13, mod.sar.op12)
summary(mod.sar.op13, Nagelkerke = TRUE) # 0.867
AIC(mod.sar.op13) # 7598.817 
rm(mod.sar.op12)
summary(mod.sar.op13)$Coef[order(summary(mod.sar.op13)$Coef[, 4]),]

mod.sar.op14 <- update(mod.sar.op13, ~. -poly(meanSST.1deg, 3):prop2.oxy + poly(meanSST.1deg, 2):prop2.oxy)
anova(mod.sar.op14, mod.sar.op13)
summary(mod.sar.op14, Nagelkerke = TRUE) # 0.86695
AIC(mod.sar.op14) # 7597.631
rm(mod.sar.op13)
summary(mod.sar.op14)$Coef[order(summary(mod.sar.op14)$Coef[, 4]),]

mod.sar.op15 <- update(mod.sar.op14, ~. -sdSal.0m:Ocean2)
anova(mod.sar.op15, mod.sar.op14)
summary(mod.sar.op15, Nagelkerke = TRUE) # 0.86681
AIC(mod.sar.op15) # 7595.539
rm(mod.sar.op14)
summary(mod.sar.op15)$Coef[order(summary(mod.sar.op15)$Coef[, 4]),]

mod.sar.op16 <- update(mod.sar.op15, ~. -logProd.mn.ann:prop2.oxy)
anova(mod.sar.op16, mod.sar.op15)
summary(mod.sar.op16, Nagelkerke = TRUE) # 0.86677
AIC(mod.sar.op16) # 7594.115
rm(mod.sar.op15)
summary(mod.sar.op16)$Coef[order(summary(mod.sar.op16)$Coef[, 4]),]

mod.sar.op17 <- update(mod.sar.op16, ~. -poly(meanSST.1deg, 3):delta_carb_ion + poly(meanSST.1deg, 2):delta_carb_ion)
anova(mod.sar.op17, mod.sar.op16)
summary(mod.sar.op17, Nagelkerke = TRUE) # 0.86672
AIC(mod.sar.op17) # 7592.81
rm(mod.sar.op16)
summary(mod.sar.op17)$Coef[order(summary(mod.sar.op17)$Coef[, 4]),]

mod.sar.op18 <- update(mod.sar.op17, ~. -poly(meanSST.1deg, 3):absMnSal.0m + poly(meanSST.1deg, 2):absMnSal.0m)
anova(mod.sar.op18, mod.sar.op17)
summary(mod.sar.op18, Nagelkerke = TRUE) # 0.86667
AIC(mod.sar.op18) # 7591.527
rm(mod.sar.op17)
summary(mod.sar.op18)$Coef[order(summary(mod.sar.op18)$Coef[, 4]),]

mod.sar.op19 <- update(mod.sar.op18, ~. -poly(meanSST.1deg, 3):logProd.mn.ann + poly(meanSST.1deg, 2):logProd.mn.ann)
anova(mod.sar.op19, mod.sar.op18)
summary(mod.sar.op19, Nagelkerke = TRUE) # 0.86662
AIC(mod.sar.op19) # 7590.173
rm(mod.sar.op18)
summary(mod.sar.op19)$Coef[order(summary(mod.sar.op19)$Coef[, 4]),]

mod.sar.op20 <- update(mod.sar.op19, ~. -logProd.mn.ann:poly(meanSST.1deg, 2) + logProd.mn.ann:poly(meanSST.1deg, 1))
anova(mod.sar.op20, mod.sar.op19)
summary(mod.sar.op20, Nagelkerke = TRUE) # 0.86662
AIC(mod.sar.op20) # 7588.197
rm(mod.sar.op19)
summary(mod.sar.op20)$Coef[order(summary(mod.sar.op20)$Coef[, 4]),]

mod.sar.op21 <- update(mod.sar.op20, ~. -I(mean.mld.t/10):prop2.oxy )
anova(mod.sar.op21, mod.sar.op20)
summary(mod.sar.op21, Nagelkerke = TRUE) # 0.86657
AIC(mod.sar.op21) # 7586.887
rm(mod.sar.op20)
summary(mod.sar.op21)$Coef[order(summary(mod.sar.op21)$Coef[, 4]),]

mod.sar.op22 <- update(mod.sar.op21, ~. -sdSST.1deg:delta_carb_ion)
anova(mod.sar.op22, mod.sar.op21)
summary(mod.sar.op22, Nagelkerke = TRUE) # 0.86651
AIC(mod.sar.op22) # 7585.703
rm(mod.sar.op21)
summary(mod.sar.op22)$Coef[order(summary(mod.sar.op22)$Coef[, 4]),]

mod.sar.op23 <- update(mod.sar.op22, ~. -poly(meanSST.1deg, 3):Ocean2 + poly(meanSST.1deg, 2):Ocean2)
anova(mod.sar.op23, mod.sar.op22)
summary(mod.sar.op23, Nagelkerke = TRUE) # 0.86642
AIC(mod.sar.op23) # 7583.039
rm(mod.sar.op22)
summary(mod.sar.op23)$Coef[order(summary(mod.sar.op23)$Coef[, 4]),]

mod.sar.op24 <- update(mod.sar.op23, ~. -prop2.oxy:delta_carb_ion)
anova(mod.sar.op24, mod.sar.op23)
summary(mod.sar.op24, Nagelkerke = TRUE) # 0.86634
AIC(mod.sar.op24) # 7582.212
rm(mod.sar.op23)
summary(mod.sar.op24)$Coef[order(summary(mod.sar.op24)$Coef[, 4]),]

mod.sar.op25 <- update(mod.sar.op24, ~. -I(mean.mld.t/10):poly(meanSST.1deg, 2) + I(mean.mld.t/10):poly(meanSST.1deg, 1))
anova(mod.sar.op25, mod.sar.op24)
summary(mod.sar.op25, Nagelkerke = TRUE) # 0.86626
AIC(mod.sar.op25) #  7581.221
rm(mod.sar.op24)
summary(mod.sar.op25)$Coef[order(summary(mod.sar.op25)$Coef[, 4]),]

mod.sar.op26 <- update(mod.sar.op25, ~. -I(mean.mld.t/10):poly(meanSST.1deg, 1))
anova(mod.sar.op26, mod.sar.op25)
summary(mod.sar.op26, Nagelkerke = TRUE) # 0.86626 
AIC(mod.sar.op26) # 7579.293
rm(mod.sar.op25)
summary(mod.sar.op26)$Coef[order(summary(mod.sar.op26)$Coef[, 4]),]

mod.sar.op27 <- update(mod.sar.op26, ~. -I(mean.mld.t/10):sdSal.0m)
anova(mod.sar.op27, mod.sar.op26)
summary(mod.sar.op27, Nagelkerke = TRUE) # 0.86621 
AIC(mod.sar.op27) # 7577.955
rm(mod.sar.op26)
summary(mod.sar.op27)$Coef[order(summary(mod.sar.op27)$Coef[, 4]),]

mod.sar.op28 <- update(mod.sar.op27, ~. -sdSST.1deg:logProd.mn.ann)
anova(mod.sar.op28, mod.sar.op27)
summary(mod.sar.op28, Nagelkerke = TRUE) # 0.86613
AIC(mod.sar.op28) # 7577.055
rm(mod.sar.op27)
summary(mod.sar.op28)$Coef[order(summary(mod.sar.op28)$Coef[, 4]),]

mod.sar.op29 <- update(mod.sar.op28, ~. -sdSST.1deg:absMnSal.0m)
anova(mod.sar.op29, mod.sar.op28)
summary(mod.sar.op29, Nagelkerke = TRUE) # 0.86605
AIC(mod.sar.op29) # 7576.18
rm(mod.sar.op28)
summary(mod.sar.op29)$Coef[order(summary(mod.sar.op29)$Coef[, 4]),]

mod.sar.op30 <- update(mod.sar.op29, ~. -logProd.mn.ann:poly(meanSST.1deg, 1))
anova(mod.sar.op30, mod.sar.op29)
summary(mod.sar.op30, Nagelkerke = TRUE) # 0.8659
AIC(mod.sar.op30) # 7576.258
rm(mod.sar.op29)
summary(mod.sar.op30)$Coef[order(summary(mod.sar.op30)$Coef[, 4]),]

mod.sar.op31 <- update(mod.sar.op30, ~. -sdSal.0m:delta_carb_ion)
anova(mod.sar.op31, mod.sar.op30)
summary(mod.sar.op31, Nagelkerke = TRUE) # 0.86567
AIC(mod.sar.op31) # 7577.488
rm(mod.sar.op30)
summary(mod.sar.op31)$Coef[order(summary(mod.sar.op31)$Coef[, 4]),]

mod.sar.op32 <- update(mod.sar.op31, ~. -I(mean.mld.t/10):delta_carb_ion)
anova(mod.sar.op32, mod.sar.op31)
summary(mod.sar.op32, Nagelkerke = TRUE) # 0.86546
AIC(mod.sar.op32) # 7578.527
rm(mod.sar.op31)
summary(mod.sar.op32)$Coef[order(summary(mod.sar.op32)$Coef[, 4]),]

mod.sar.op33 <- update(mod.sar.op32, ~. -I(depth10deg/100):absMnSal.0m)
anova(mod.sar.op33, mod.sar.op32)
summary(mod.sar.op33, Nagelkerke = TRUE) # 0.8652
AIC(mod.sar.op33) # 7580.091
rm(mod.sar.op32)
summary(mod.sar.op33)$Coef[order(summary(mod.sar.op33)$Coef[, 4]),]

# having got to here (which is that all values are <0.05 or can't be removed), then check whether any others should be removed (particularly those with ocean)
mod.sar.op34 <- update(mod.sar.op33, ~. -absMnSal.0m:Ocean2)
anova(mod.sar.op34, mod.sar.op33) # this is 0.06, so the variable should be removed
summary(mod.sar.op34, Nagelkerke = TRUE) # 0.8648
AIC(mod.sar.op34) # 7581.642
rm(mod.sar.op33)
summary(mod.sar.op34)$Coef[order(summary(mod.sar.op34)$Coef[, 4]),]

mod.sar.op35 <- update(mod.sar.op34, ~. -I(mean.mld.t/10):absMnSal.0m)
anova(mod.sar.op35, mod.sar.op34)
summary(mod.sar.op35, Nagelkerke = TRUE) # 0.86457
AIC(mod.sar.op35) # 7582.792
rm(mod.sar.op34)
summary(mod.sar.op35)$Coef[order(summary(mod.sar.op35)$Coef[, 4]),]

mod.sar.op36 <- update(mod.sar.op35, ~. -absMnSal.0m:prop2.oxy)
anova(mod.sar.op36, mod.sar.op35)
summary(mod.sar.op36, Nagelkerke = TRUE) # 0.86434
AIC(mod.sar.op36) # 7583.963
rm(mod.sar.op35)
summary(mod.sar.op36)$Coef[order(summary(mod.sar.op36)$Coef[, 4]),]

mod.sar.op37 <- update(mod.sar.op36, ~. -sdSal.0m:prop2.oxy)
anova(mod.sar.op37, mod.sar.op36) # just not significant (0.050159)
summary(mod.sar.op37, Nagelkerke = TRUE) # 0.86407
AIC(mod.sar.op37) # 7585.799
rm(mod.sar.op36)
summary(mod.sar.op37)$Coef[order(summary(mod.sar.op37)$Coef[, 4]),]

# removing anything else doesn't help, so this is as simplified as possible
mod.sar.op38 <- update(mod.sar.op37, ~. -sdSST.1deg:prop2.oxy)
anova(mod.sar.op38, mod.sar.op37)

mod.sar.opf <- mod.sar.op37
rm(mod.sar.op37, mod.sar.op38)

## 2v. Create a plot of model parameters --------------------------------------
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.86407

# generate a dataframe of coefficients
(ms.coef <- data.frame(names = names(mod.sar.opf$coefficients), coef.sar = mod.sar.opf$coefficients, row.names = 1:length(mod.sar.opf$coefficients), stars = NA))

# reorder the rows to something more sensible
order.coef.ms <- c(1:17, 37:47, 18:36)
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
png("Figures/Ana_2v_coef_modsaropf.png", width = 1000, height = 750)
plt.def <- par("plt")
par(plt = c(plt.def[1:2], 0.5, plt.def[4]))
tmp.x <- barplot(abs(ms.coef$coef.sar), names = ms.coef$names, las = 2, ylim = c(0, max(abs(ms.coef$coef.sar)) + 50))
text(tmp.x, abs(ms.coef$coef.sar) + 20, ms.coef$stars)
par(plt = plt.def)
dev.off()
rm(plt.def, tmp.x, ms.coef)

## 2vi. Calculate likelihood ratios for the SAR model ----------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.86407
AIC(mod.sar.opf) # 7585.799

# removing mean temp^3
mod.sar.lr.mnt3 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 2) - poly(meanSST.1deg, 3):sdSal.0m + poly(meanSST.1deg, 2):sdSal.0m)
summary(mod.sar.lr.mnt3, Nagelkerke = T) # r2 = 0.85219
AIC(mod.sar.lr.mnt3) # 7738.831
lrtest(mod.sar.opf, mod.sar.lr.mnt3) # n.b. this is basically an anova
# LR = < 2.2e-16 ***

# removing mean temp^2
mod.sar.lr.mnt2 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 1) - poly(meanSST.1deg, 3):sdSal.0m + poly(meanSST.1deg, 1):sdSal.0m - prop2.oxy:poly(meanSST.1deg, 2) + prop2.oxy:poly(meanSST.1deg, 1) - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 1) - absMnSal.0m:poly(meanSST.1deg, 2) + absMnSal.0m:poly(meanSST.1deg, 1) - Ocean2:poly(meanSST.1deg, 2) + Ocean2:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mnt2, Nagelkerke = T) # r2 = 0.83775
AIC(mod.sar.lr.mnt2) # 7899.629
lrtest(mod.sar.opf, mod.sar.lr.mnt2)
# LR = < 2.2e-16 ***

# removing mean temp
mod.sar.lr.mnt <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) - poly(meanSST.1deg, 3):sdSal.0m - prop2.oxy:poly(meanSST.1deg, 2) - delta_carb_ion:poly(meanSST.1deg, 2) - absMnSal.0m:poly(meanSST.1deg, 2) - Ocean2:poly(meanSST.1deg, 2) - sdSST.1deg:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mnt, Nagelkerke = T) # r2 = 0.78477
AIC(mod.sar.lr.mnt) # 8413.431
lrtest(mod.sar.opf, mod.sar.lr.mnt)
# LR = < 2.2e-16 ***

# removing sd temp
mod.sar.lr.sdt <- update(mod.sar.opf, ~. -sdSST.1deg - sdSST.1deg:I(mean.mld.t/10) - sdSST.1deg:sdSal.0m - sdSST.1deg:prop2.oxy - sdSST.1deg:Ocean2 - sdSST.1deg:poly(meanSST.1deg, 1))
summary(mod.sar.lr.sdt, Nagelkerke = T) # r2 = 0.86156
AIC(mod.sar.lr.sdt) # 7606.092
lrtest(mod.sar.opf, mod.sar.lr.sdt)
# LR = 1.518e-05 ***

# removing mld temp
mod.sar.lr.mld <- update(mod.sar.opf, ~. -I(mean.mld.t/10) - sdSST.1deg:I(mean.mld.t/10) - I(mean.mld.t/10):logProd.mn.ann - I(mean.mld.t/10):Ocean2)
summary(mod.sar.lr.mld, Nagelkerke = T) # r2 = 0.85942
AIC(mod.sar.lr.mld) # 7638.864
lrtest(mod.sar.opf, mod.sar.lr.mld)
# LR = 2.822e-12 ***

# removing depth of 10 degree contour
mod.sar.lr.d10 <- update(mod.sar.opf, ~. -I(depth10deg/100) - I(depth10deg/100):logProd.mn.ann - I(depth10deg/100):sdSal.0m - I(depth10deg/100):prop2.oxy)
summary(mod.sar.lr.d10, Nagelkerke = T) # r2 = 0.86106
AIC(mod.sar.lr.d10) # 7618.825
lrtest(mod.sar.opf, mod.sar.lr.d10)
# LR = 2.654e-08 ***

# removing mean log Prod
mod.sar.lr.prod <- update(mod.sar.opf, ~. -logProd.mn.ann - I(mean.mld.t/10):logProd.mn.ann - I(depth10deg/100):logProd.mn.ann - logProd.mn.ann:absMnSal.0m - logProd.mn.ann:sdSal.0m - logProd.mn.ann:Ocean2)
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.8544
AIC(mod.sar.lr.prod) # 7700.535
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = < 2.2e-16 ***

# removing mean salinity
mod.sar.lr.msal <- update(mod.sar.opf, ~. -absMnSal.0m - logProd.mn.ann:absMnSal.0m - absMnSal.0m:poly(meanSST.1deg, 2))
summary(mod.sar.lr.msal, Nagelkerke = T) # r2 = 0.85924
AIC(mod.sar.lr.msal) # 7643.248
lrtest(mod.sar.opf, mod.sar.lr.msal)
# LR = 2.07e-13 ***

# removing sd salinity
mod.sar.lr.sdsal <- update(mod.sar.opf, ~. -sdSal.0m - poly(meanSST.1deg, 3):sdSal.0m - sdSST.1deg:sdSal.0m - I(depth10deg/100):sdSal.0m - logProd.mn.ann:sdSal.0m)
summary(mod.sar.lr.sdsal, Nagelkerke = T) # r2 = 0.8616
AIC(mod.sar.lr.sdsal) # 7605.455
lrtest(mod.sar.opf, mod.sar.lr.sdsal)
# LR = 1.997e-05 ***

# removing prop2.oxy
mod.sar.lr.oxy <- update(mod.sar.opf, ~. -prop2.oxy - sdSST.1deg:prop2.oxy - I(depth10deg/100):prop2.oxy - prop2.oxy:Ocean2 - prop2.oxy:poly(meanSST.1deg, 2))
summary(mod.sar.lr.oxy, Nagelkerke = T) # r2 = 0.86101
AIC(mod.sar.lr.oxy) # 7613.419
lrtest(mod.sar.opf, mod.sar.lr.oxy)
# LR = 6.155e-07 ***

# removing Ocean2
mod.sar.lr.oce <- update(mod.sar.opf, ~. -Ocean2 - sdSST.1deg:Ocean2 - I(mean.mld.t/10):Ocean2 - logProd.mn.ann:Ocean2 - prop2.oxy:Ocean2 - Ocean2:delta_carb_ion - Ocean2:poly(meanSST.1deg, 2))
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.85089
AIC(mod.sar.lr.oce) # 7727.268
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = < 2.2e-16 ***

# removing delta_carb_ion
mod.sar.lr.dis <- update(mod.sar.opf, ~. -delta_carb_ion - Ocean2:delta_carb_ion - delta_carb_ion:poly(meanSST.1deg, 2))
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.86315
AIC(mod.sar.lr.dis) # 7588.323
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.02828 *

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
png("Figures/Ana_2vi_LRatio_sar_opf.png", width = 550)
tmp.x <- barplot(ms.lr$lr, names = ms.lr$names, las = 2, ylim = c(0, max(ms.lr$lr) + 30))
text(tmp.x, ms.lr$lr + 20, ms.lr$stars)
dev.off()
rm(tmp.x)

## 2vii. Calculate likelihood ratios for groups of EVs ---------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.86407
AIC(mod.sar.opf) # 7585.799

# removing temp
mod.sar.lr.temp <- mod.sar.lr.mnt
summary(mod.sar.lr.temp, Nagelkerke = T) # r2 = 0.78477
AIC(mod.sar.lr.temp) # 8413.431
lrtest(mod.sar.opf, mod.sar.lr.temp)
# LR = < 2.2e-16 ***

# removing structure
mod.sar.lr.str <- update(mod.sar.opf, ~. -I(mean.mld.t/10) - I(depth10deg/100) - sdSST.1deg:I(mean.mld.t/10) - I(mean.mld.t/10):logProd.mn.ann - I(mean.mld.t/10):Ocean2 - I(depth10deg/100):logProd.mn.ann - I(depth10deg/100):sdSal.0m - I(depth10deg/100):prop2.oxy)
summary(mod.sar.lr.str, Nagelkerke = T) # r2 = 0.85729
AIC(mod.sar.lr.str) # 7659
lrtest(mod.sar.opf, mod.sar.lr.str)
# LR = 9.344e-16 ***

# removing stability
mod.sar.lr.stable <- update(mod.sar.opf, ~. -sdSST.1deg - sdSal.0m - sdSST.1deg:I(mean.mld.t/10) - sdSST.1deg:sdSal.0m - sdSST.1deg:prop2.oxy - sdSST.1deg:Ocean2 - sdSST.1deg:poly(meanSST.1deg, 1) - poly(meanSST.1deg, 3):sdSal.0m - sdSST.1deg:sdSal.0m - I(depth10deg/100):sdSal.0m - logProd.mn.ann:sdSal.0m)
summary(mod.sar.lr.stable, Nagelkerke = T) # r2 = 0.85953
AIC(mod.sar.lr.stable) # 7621.302
lrtest(mod.sar.opf, mod.sar.lr.stable)
# LR = 2.825e-08 ***

# removing productivity
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.8544
AIC(mod.sar.lr.prod) # 7700.535
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = < 2.2e-16 ***

# removing stress
mod.sar.lr.sts <- update(mod.sar.opf, ~. -absMnSal.0m - prop2.oxy - logProd.mn.ann:absMnSal.0m - absMnSal.0m:poly(meanSST.1deg, 2) - sdSST.1deg:prop2.oxy - I(depth10deg/100):prop2.oxy - prop2.oxy:Ocean2 - prop2.oxy:poly(meanSST.1deg, 2))
summary(mod.sar.lr.sts, Nagelkerke = T) # r2 = 0.85654
AIC(mod.sar.lr.sts) # 7664.827
lrtest(mod.sar.opf, mod.sar.lr.sts)
# LR =  < 2.2e-16 ***

# removing Ocean2
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.85089
AIC(mod.sar.lr.oce) # 7727.268
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = < 2.2e-16 ***

# removing delta_carb_ion
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.86315
AIC(mod.sar.lr.dis) # 7588.323
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.02828 *

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
png("Figures/Ana_2vii_LRatio_sar_opf_group.png", width = 550)
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

png("Figures/Ana_2viii_LRatio_opfop0.png", width = 800)
# get order from running without order first
lr.plot(lr.sar.op0, lr.sar.opf, order = c(7:9, 4:3, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.op0g <- lr.calc(mod.sar.op0, tmp)
rm(tmp)

(tmp <- data.frame(names = model.evs(mod.sar.opf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.opfg <- lr.calc(mod.sar.opf, tmp)
rm(tmp)

png("Figures/Ana_2viii_LRatio_g_opfop0.png", width = 800)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(6:7, 4:3, 5, 2:1), leg.x = 17, leg.y = 800, leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

## 2ix. Does optimisation differ for simplified? ---------------------------
mod.fw.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fw.rop # 565.2149
AIC(mod.fw.rop$obj)
# 7590.583

mod.fb.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fb.rop # 501.6069
AIC(mod.fb.rop$obj)
# 7647.389

mod.fs.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fs.rop # 505.3868
AIC(mod.fs.rop$obj)
# 7586.166

# therefore justified in using S
rm(mod.fw.rop, mod.fb.rop, mod.fs.rop, ms.lr, ms.lr.group)


## 3. How do the different Oceans compare? ---------------------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.86407
AIC(mod.sar.opf) # 7585.799

# do I have suffient points for each ocean?
table(ldg.margo.mod$Ocean2) # should do
# c.f. Atlantic: 372   Indian: 157  Pacific: 146, which were the values for bfd

## 3i. set up model for Atlantic -----------------------------
# run model based on best model for complete data with only Atlantic
atl.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Atlantic", ], 0, mod.sar.opS$dist, longlat = TRUE)
atl.s <- nb2listw(atl.nb, glist = NULL, style = "S", zero.policy = TRUE)
rm(atl.nb)

## 3ii. set up model for Indian -----------------------------------
# run model based on best model for complete data with only Indian
ind.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Indian", ], 0, mod.sar.opS$dist, longlat = TRUE)
ind.s <- nb2listw(ind.nb, glist = NULL, style = "S", zero.policy = TRUE)
rm(ind.nb)

## 3iii. set up model for Pacific  ---------------------------
pac.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Pacific", ], 0, mod.sar.opS$dist, longlat = TRUE)
pac.s <- nb2listw(pac.nb, glist = NULL, style = "S", zero.policy = TRUE)
rm(pac.nb)

# n.b. removed the simplification in 3i-3iii as I'm only simplifying from the model with significant interactions (see 3vii). Therefore, also didn't need 3iv - vi

## 3vii. only significant interactions for Atlantic -----------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)
# get formula without the Ocean2
op.formula <- update(mod.sar.opf$call$formula, ~.-Ocean2 - sdSST.1deg:Ocean2 - I(mean.mld.t/10):Ocean2 - logProd.mn.ann:Ocean2 - prop2.oxy:Ocean2 - Ocean2:delta_carb_ion - Ocean2:poly(meanSST.1deg, 2))

# create a model with only significant values
mod.sar.atlI <- errorsarlm(op.formula, listw = atl.s, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.atlI)
par(mfrow = c(1, 1))

# try adding in the other interactions
mod.sar.atlI2 <- update(mod.sar.atlI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI2)
anova(mod.sar.atlI, mod.sar.atlI2) # 0.29206
rm(mod.sar.atlI2)

mod.sar.atlI3 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.atlI3)
anova(mod.sar.atlI, mod.sar.atlI3) # 0.34832
rm(mod.sar.atlI3)

mod.sar.atlI4 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.atlI4)
anova(mod.sar.atlI, mod.sar.atlI4) # 0.38968
rm(mod.sar.atlI4)

mod.sar.atlI5 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI5)
anova(mod.sar.atlI, mod.sar.atlI5) # 0.57059
rm(mod.sar.atlI5)

mod.sar.atlI6 <- update(mod.sar.atlI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.atlI6)
anova(mod.sar.atlI, mod.sar.atlI6) # 0.12442
rm(mod.sar.atlI6)

mod.sar.atlI7 <- update(mod.sar.atlI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.atlI7)
anova(mod.sar.atlI, mod.sar.atlI7) # 0.66236
rm(mod.sar.atlI7)

mod.sar.atlI9 <- update(mod.sar.atlI, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atlI9)
anova(mod.sar.atlI, mod.sar.atlI9) # 0.084804
rm(mod.sar.atlI9)

mod.sar.atlI10 <- update(mod.sar.atlI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.atlI10)
anova(mod.sar.atlI, mod.sar.atlI10) # 0.4846
rm(mod.sar.atlI10)

mod.sar.atlI11 <- update(mod.sar.atlI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.atlI11)
anova(mod.sar.atlI, mod.sar.atlI11) # 0.95632
rm(mod.sar.atlI11)

mod.sar.atlI12 <- update(mod.sar.atlI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.atlI12)
anova(mod.sar.atlI, mod.sar.atlI12) # 0.060058
rm(mod.sar.atlI12)

mod.sar.atlI13 <- update(mod.sar.atlI, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.atlI13)
anova(mod.sar.atlI, mod.sar.atlI13) # 0.60343
rm(mod.sar.atlI13)

mod.sar.atlI14 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.atlI14)
anova(mod.sar.atlI, mod.sar.atlI14) # 0.79249
rm(mod.sar.atlI14)

mod.sar.atlI15 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.atlI15)
anova(mod.sar.atlI, mod.sar.atlI15) # 0.91004
rm(mod.sar.atlI15)

mod.sar.atlI16 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.atlI16)
anova(mod.sar.atlI, mod.sar.atlI16) # 0.38174
rm(mod.sar.atlI16)

mod.sar.atlI17 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.atlI17)
anova(mod.sar.atlI, mod.sar.atlI17) # 0.54658
rm(mod.sar.atlI17)

mod.sar.atlI18 <- update(mod.sar.atlI, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.atlI18)
anova(mod.sar.atlI, mod.sar.atlI18) # 0.82099
rm(mod.sar.atlI18)

mod.sar.atlI19 <- update(mod.sar.atlI, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.atlI19)
anova(mod.sar.atlI, mod.sar.atlI19) # 0.59692
rm(mod.sar.atlI19)

mod.sar.atlI21 <- update(mod.sar.atlI, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.atlI21)
anova(mod.sar.atlI, mod.sar.atlI21) # 0.39067
rm(mod.sar.atlI21)

# significant at the 0.05 level
mod.sar.atlI22 <- update(mod.sar.atlI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.atlI22)
anova(mod.sar.atlI, mod.sar.atlI22) # 0.03498
rm(mod.sar.atlI22)

mod.sar.atlI23 <- update(mod.sar.atlI, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.atlI23)
anova(mod.sar.atlI, mod.sar.atlI23) # 0.74448
rm(mod.sar.atlI23)

mod.sar.atlI24 <- update(mod.sar.atlI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.atlI24)
anova(mod.sar.atlI, mod.sar.atlI24) # 0.77829
rm(mod.sar.atlI24)

mod.sar.atlI25 <- update(mod.sar.atlI, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.atlI25)
anova(mod.sar.atlI, mod.sar.atlI25) # 0.85475
rm(mod.sar.atlI25)

mod.sar.atlI27 <- update(mod.sar.atlI, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.atlI27)
anova(mod.sar.atlI, mod.sar.atlI27) # 0.89812
rm(mod.sar.atlI27)

mod.sar.atlI28 <- update(mod.sar.atlI, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.atlI28)
anova(mod.sar.atlI, mod.sar.atlI28) # 0.84054
rm(mod.sar.atlI28)

# significant at the 0.05 level
mod.sar.atlI30 <- update(mod.sar.atlI, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.atlI30)
anova(mod.sar.atlI, mod.sar.atlI30) # 0.019602

mod.sar.atlI31 <- update(mod.sar.atlI, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.atlI31)
anova(mod.sar.atlI, mod.sar.atlI31) # 0.67247
rm(mod.sar.atlI31)

## adding back in finds sdSal.0m:delta_carb_ion should be included in the model, also potentially logProd.mn.ann:prop2.oxy
mod.sar.atl1.I <- mod.sar.atlI30
rm(mod.sar.atlI30)

# check whether any of the interactions are now significant
mod.sar.atl1.I2 <- update(mod.sar.atl1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atl1.I2)
anova(mod.sar.atl1.I, mod.sar.atl1.I2) # 0.35061
rm(mod.sar.atl1.I2)

mod.sar.atl1.I3 <- update(mod.sar.atl1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.atl1.I3)
anova(mod.sar.atl1.I, mod.sar.atl1.I3) # 0.36778
rm(mod.sar.atl1.I3)

mod.sar.atl1.I4 <- update(mod.sar.atl1.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.atl1.I4)
anova(mod.sar.atl1.I, mod.sar.atl1.I4) # 0.4248
rm(mod.sar.atl1.I4)

mod.sar.atl1.I5 <- update(mod.sar.atl1.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atl1.I5)
anova(mod.sar.atl1.I, mod.sar.atl1.I5) # 0.53555
rm(mod.sar.atl1.I5)

mod.sar.atl1.I6 <- update(mod.sar.atl1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.atl1.I6)
anova(mod.sar.atl1.I, mod.sar.atl1.I6) # 0.12247
rm(mod.sar.atl1.I6)

mod.sar.atl1.I7 <- update(mod.sar.atl1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.atl1.I7)
anova(mod.sar.atl1.I, mod.sar.atl1.I7) # 0.59396
rm(mod.sar.atl1.I7)

mod.sar.atl1.I9 <- update(mod.sar.atl1.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atl1.I9)
anova(mod.sar.atl1.I, mod.sar.atl1.I9) # 0.23402
rm(mod.sar.atl1.I9)

mod.sar.atl1.I10 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.atl1.I10)
anova(mod.sar.atl1.I, mod.sar.atl1.I10) # 0.45743
rm(mod.sar.atl1.I10)

mod.sar.atl1.I11 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.atl1.I11)
anova(mod.sar.atl1.I, mod.sar.atl1.I11) # 0.91279
rm(mod.sar.atl1.I11)

mod.sar.atl1.I12 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.atl1.I12)
anova(mod.sar.atl1.I, mod.sar.atl1.I12) # 0.055326
rm(mod.sar.atl1.I12)

mod.sar.atl1.I13 <- update(mod.sar.atl1.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.atl1.I13)
anova(mod.sar.atl1.I, mod.sar.atl1.I13) # 0.90115
rm(mod.sar.atl1.I13)

mod.sar.atl1.I14 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.atl1.I14)
anova(mod.sar.atl1.I, mod.sar.atl1.I14) # 0.78888
rm(mod.sar.atl1.I14)

mod.sar.atl1.I15 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.atl1.I15)
anova(mod.sar.atl1.I, mod.sar.atl1.I15) # 0.96784
rm(mod.sar.atl1.I15)

mod.sar.atl1.I16 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.atl1.I16)
anova(mod.sar.atl1.I, mod.sar.atl1.I16) # 0.38872
rm(mod.sar.atl1.I16)

mod.sar.atl1.I17 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.atl1.I17)
anova(mod.sar.atl1.I, mod.sar.atl1.I17) # 0.51591
rm(mod.sar.atl1.I17)

mod.sar.atl1.I18 <- update(mod.sar.atl1.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.atl1.I18)
anova(mod.sar.atl1.I, mod.sar.atl1.I18) # 0.61322
rm(mod.sar.atl1.I18)

mod.sar.atl1.I19 <- update(mod.sar.atl1.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.atl1.I19)
anova(mod.sar.atl1.I, mod.sar.atl1.I19) # 0.60771
rm(mod.sar.atl1.I19)

mod.sar.atl1.I21 <- update(mod.sar.atl1.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.atl1.I21)
anova(mod.sar.atl1.I, mod.sar.atl1.I21) # 0.43443
rm(mod.sar.atl1.I21)

# significant at the 0.05 level
mod.sar.atl1.I22 <- update(mod.sar.atl1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.atl1.I22)
anova(mod.sar.atl1.I, mod.sar.atl1.I22) # 0.019167

mod.sar.atl1.I23 <- update(mod.sar.atl1.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.atl1.I23)
anova(mod.sar.atl1.I, mod.sar.atl1.I23) # 0.53117
rm(mod.sar.atl1.I23)

mod.sar.atl1.I24 <- update(mod.sar.atl1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.atl1.I24)
anova(mod.sar.atl1.I, mod.sar.atl1.I24) # 0.66135
rm(mod.sar.atl1.I24)

mod.sar.atl1.I25 <- update(mod.sar.atl1.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.atl1.I25)
anova(mod.sar.atl1.I, mod.sar.atl1.I25) # 0.79875
rm(mod.sar.atl1.I25)

mod.sar.atl1.I27 <- update(mod.sar.atl1.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.atl1.I27)
anova(mod.sar.atl1.I, mod.sar.atl1.I27) # 0.43224
rm(mod.sar.atl1.I27)

mod.sar.atl1.I28 <- update(mod.sar.atl1.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.atl1.I28)
anova(mod.sar.atl1.I, mod.sar.atl1.I28) # 0.85621
rm(mod.sar.atl1.I28)

mod.sar.atl1.I31 <- update(mod.sar.atl1.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.atl1.I31)
anova(mod.sar.atl1.I, mod.sar.atl1.I31) # 0.46566
rm(mod.sar.atl1.I31)

## adding back in finds logProd.mn.ann:prop2.oxy
mod.sar.atl2.I <- mod.sar.atl1.I22
rm(mod.sar.atl1.I22)

# check whether any of the interactions are now significant
mod.sar.atl2.I2 <- update(mod.sar.atl2.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atl2.I2)
anova(mod.sar.atl2.I, mod.sar.atl2.I2) # 0.25463
rm(mod.sar.atl2.I2)

mod.sar.atl2.I3 <- update(mod.sar.atl2.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.atl2.I3)
anova(mod.sar.atl2.I, mod.sar.atl2.I3) # 0.40342
rm(mod.sar.atl2.I3)

mod.sar.atl2.I4 <- update(mod.sar.atl2.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.atl2.I4)
anova(mod.sar.atl2.I, mod.sar.atl2.I4) # 0.36932
rm(mod.sar.atl2.I4)

mod.sar.atl2.I5 <- update(mod.sar.atl2.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atl2.I5)
anova(mod.sar.atl2.I, mod.sar.atl2.I5) # 0.40448
rm(mod.sar.atl2.I5)

mod.sar.atl2.I6 <- update(mod.sar.atl2.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.atl2.I6)
anova(mod.sar.atl2.I, mod.sar.atl2.I6) # 0.12673
rm(mod.sar.atl2.I6)

mod.sar.atl2.I7 <- update(mod.sar.atl2.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.atl2.I7)
anova(mod.sar.atl2.I, mod.sar.atl2.I7) # 0.32729
rm(mod.sar.atl2.I7)

mod.sar.atl2.I9 <- update(mod.sar.atl2.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atl2.I9)
anova(mod.sar.atl2.I, mod.sar.atl2.I9) # 0.32372
rm(mod.sar.atl2.I9)

mod.sar.atl2.I10 <- update(mod.sar.atl2.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.atl2.I10)
anova(mod.sar.atl2.I, mod.sar.atl2.I10) # 0.38486
rm(mod.sar.atl2.I10)

mod.sar.atl2.I11 <- update(mod.sar.atl2.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.atl2.I11)
anova(mod.sar.atl2.I, mod.sar.atl2.I11) # 0.55001
rm(mod.sar.atl2.I11)

mod.sar.atl2.I12 <- update(mod.sar.atl2.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.atl2.I12)
anova(mod.sar.atl2.I, mod.sar.atl2.I12) # 0.050214
rm(mod.sar.atl2.I12)

mod.sar.atl2.I13 <- update(mod.sar.atl2.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.atl2.I13)
anova(mod.sar.atl2.I, mod.sar.atl2.I13) # 0.78259
rm(mod.sar.atl2.I13)

mod.sar.atl2.I14 <- update(mod.sar.atl2.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.atl2.I14)
anova(mod.sar.atl2.I, mod.sar.atl2.I14) # 0.93954
rm(mod.sar.atl2.I14)

mod.sar.atl2.I15 <- update(mod.sar.atl2.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.atl2.I15)
anova(mod.sar.atl2.I, mod.sar.atl2.I15) # 0.96913
rm(mod.sar.atl2.I15)

mod.sar.atl2.I16 <- update(mod.sar.atl2.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.atl2.I16)
anova(mod.sar.atl2.I, mod.sar.atl2.I16) # 0.31357
rm(mod.sar.atl2.I16)

mod.sar.atl2.I17 <- update(mod.sar.atl2.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.atl2.I17)
anova(mod.sar.atl2.I, mod.sar.atl2.I17) # 0.38476
rm(mod.sar.atl2.I17)

mod.sar.atl2.I18 <- update(mod.sar.atl2.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.atl2.I18)
anova(mod.sar.atl2.I, mod.sar.atl2.I18) # 0.90752
rm(mod.sar.atl2.I18)

mod.sar.atl2.I19 <- update(mod.sar.atl2.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.atl2.I19)
anova(mod.sar.atl2.I, mod.sar.atl2.I19) # 0.60374
rm(mod.sar.atl2.I19)

mod.sar.atl2.I21 <- update(mod.sar.atl2.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.atl2.I21)
anova(mod.sar.atl2.I, mod.sar.atl2.I21) # 0.41062
rm(mod.sar.atl2.I21)

mod.sar.atl2.I23 <- update(mod.sar.atl2.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.atl2.I23)
anova(mod.sar.atl2.I, mod.sar.atl2.I23) # 0.93479
rm(mod.sar.atl2.I23)

mod.sar.atl2.I24 <- update(mod.sar.atl2.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.atl2.I24)
anova(mod.sar.atl2.I, mod.sar.atl2.I24) # 0.45112
rm(mod.sar.atl2.I24)

mod.sar.atl2.I25 <- update(mod.sar.atl2.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.atl2.I25)
anova(mod.sar.atl2.I, mod.sar.atl2.I25) # 0.64701
rm(mod.sar.atl2.I25)

mod.sar.atl2.I27 <- update(mod.sar.atl2.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.atl2.I27)
anova(mod.sar.atl2.I, mod.sar.atl2.I27) # 0.51563
rm(mod.sar.atl2.I27)

mod.sar.atl2.I28 <- update(mod.sar.atl2.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.atl2.I28)
anova(mod.sar.atl2.I, mod.sar.atl2.I28) # 0.26055
rm(mod.sar.atl2.I28)

mod.sar.atl2.I31 <- update(mod.sar.atl2.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.atl2.I31)
anova(mod.sar.atl2.I, mod.sar.atl2.I31) # 0.74811
rm(mod.sar.atl2.I31)

# nothing is significant, so 
mod.sar.atlIf <- mod.sar.atl2.I
summary(mod.sar.atlIf, Nagelkerke = TRUE) # 0.90757 
AIC(mod.sar.atlIf) # 4413.503

# Check for correlation
env.var.atl <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion", "prop2.oxy")

# pairs plot
png("Figures/Ana_3vii_atl_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])

lr.sar.atlIf <- lr.calc(mod.sar.atlIf)

png("Figures/Ana_3vii_LRatio_atlIf.png", width = 800)
lr.plot(lr.sar.atlIf, order = c(7:9, 10:11, 1, 4, 2:3, 5:6), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.atlIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.atlIfg <- lr.calc(mod.sar.atlIf, tmp)
rm(tmp)

png("Figures/Ana_3vii_LRatio_g_atlIf.png", width = 800)
lr.plot(lr.sar.atlIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.atlI, mod.sar.atl1.I, mod.sar.atl2.I, file = "Outputs/Atlantic_simplification.RData")
rm(env.var.atl, mod.sar.atlI, mod.sar.atl1.I, mod.sar.atl2.I)

## 3viii. Only significant interactions for Indian -------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.indI <- errorsarlm(op.formula, listw = ind.s, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.indI)
par(mfrow = c(1, 1))
summary(mod.sar.indI, Nagelkerke = TRUE) # 0.8496 
AIC(mod.sar.indI) # 1050.043

# try adding in the other interactions
mod.sar.indI2 <- update(mod.sar.indI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.indI2)
anova(mod.sar.indI, mod.sar.indI2) # 0.23316
rm(mod.sar.indI2)

mod.sar.indI3 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.indI3)
anova(mod.sar.indI, mod.sar.indI3) # 0.16103
rm(mod.sar.indI3)

mod.sar.indI4 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.indI4)
anova(mod.sar.indI, mod.sar.indI4) # 0.71998
rm(mod.sar.indI4)

mod.sar.indI5 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.indI5)
anova(mod.sar.indI, mod.sar.indI5) # 0.82565
rm(mod.sar.indI5)

mod.sar.indI6 <- update(mod.sar.indI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.indI6)
anova(mod.sar.indI, mod.sar.indI6) # 0.89703
rm(mod.sar.indI6)

mod.sar.indI7 <- update(mod.sar.indI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.indI7)
anova(mod.sar.indI, mod.sar.indI7) # 0.78843
rm(mod.sar.indI7)

# significant at the 0.05 level
mod.sar.indI9 <- update(mod.sar.indI, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.indI9)
anova(mod.sar.indI, mod.sar.indI9) # 0.017565
rm(mod.sar.indI9)

mod.sar.indI10 <- update(mod.sar.indI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.indI10)
anova(mod.sar.indI, mod.sar.indI10) # 0.17926
rm(mod.sar.indI10)

mod.sar.indI11 <- update(mod.sar.indI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.indI11)
anova(mod.sar.indI, mod.sar.indI11) # 0.95716
rm(mod.sar.indI11)

mod.sar.indI12 <- update(mod.sar.indI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.indI12)
anova(mod.sar.indI, mod.sar.indI12) # 0.75331 
rm(mod.sar.indI12)

mod.sar.indI13 <- update(mod.sar.indI, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.indI13)
anova(mod.sar.indI, mod.sar.indI13) # 0.11561
rm(mod.sar.indI13)

mod.sar.indI14 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.indI14)
anova(mod.sar.indI, mod.sar.indI14) # 0.3944
rm(mod.sar.indI14)

mod.sar.indI15 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.indI15)
anova(mod.sar.indI, mod.sar.indI15) # 0.55653
rm(mod.sar.indI15)

mod.sar.indI16 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.indI16)
anova(mod.sar.indI, mod.sar.indI16) # 0.2993
rm(mod.sar.indI16)

# significant at the 0.05 level
mod.sar.indI17 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.indI17)
anova(mod.sar.indI, mod.sar.indI17) # 0.0052249

mod.sar.indI18 <- update(mod.sar.indI, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.indI18)
anova(mod.sar.indI, mod.sar.indI18) # 0.20616
rm(mod.sar.indI18)

mod.sar.indI19 <- update(mod.sar.indI, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.indI19)
anova(mod.sar.indI, mod.sar.indI19) # 0.17104
rm(mod.sar.indI19)

mod.sar.indI21 <- update(mod.sar.indI, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.indI21)
anova(mod.sar.indI, mod.sar.indI21) # 0.51582
rm(mod.sar.indI21)

mod.sar.indI22 <- update(mod.sar.indI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.indI22)
anova(mod.sar.indI, mod.sar.indI22) # 0.83483
rm(mod.sar.indI22)

mod.sar.indI23 <- update(mod.sar.indI, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.indI23)
anova(mod.sar.indI, mod.sar.indI23) # 0.12629
rm(mod.sar.indI23)

mod.sar.indI24 <- update(mod.sar.indI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.indI24)
anova(mod.sar.indI, mod.sar.indI24) # 0.73389
rm(mod.sar.indI24)

mod.sar.indI25 <- update(mod.sar.indI, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.indI25)
anova(mod.sar.indI, mod.sar.indI25) # 0.71142
rm(mod.sar.indI25)

mod.sar.indI27 <- update(mod.sar.indI, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.indI27)
anova(mod.sar.indI, mod.sar.indI27) # 0.30243
rm(mod.sar.indI27)

# significant at the 0.05 level
mod.sar.indI28 <- update(mod.sar.indI, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.indI28)
anova(mod.sar.indI, mod.sar.indI28) # 0.048386
rm(mod.sar.indI28)

mod.sar.indI30 <- update(mod.sar.indI, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.indI30)
anova(mod.sar.indI, mod.sar.indI30) # 0.32894
rm(mod.sar.indI30)

mod.sar.indI31 <- update(mod.sar.indI, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.indI31)
anova(mod.sar.indI, mod.sar.indI31) # 0.053794
rm(mod.sar.indI31)


# I(mean.mld.t/10):prop2.oxy is significant at the 0.05 level (as is sdSal.0m:prop2.oxy and delta_carb_ion:poly(meanSST.1deg, 3) ), so add it in and recheck
mod.sar.ind1.I <- mod.sar.indI17
rm(mod.sar.indI17)

mod.sar.ind1.I2 <- update(mod.sar.ind1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind1.I2)
anova(mod.sar.ind1.I, mod.sar.ind1.I2) # 0.18636
rm(mod.sar.ind1.I2)

# significant at the 0.05 level
mod.sar.ind1.I3 <- update(mod.sar.ind1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.ind1.I3)
anova(mod.sar.ind1.I, mod.sar.ind1.I3) # 0.009128

mod.sar.ind1.I4 <- update(mod.sar.ind1.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind1.I4)
anova(mod.sar.ind1.I, mod.sar.ind1.I4) # 0.69709
rm(mod.sar.ind1.I4)

mod.sar.ind1.I5 <- update(mod.sar.ind1.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind1.I5)
anova(mod.sar.ind1.I, mod.sar.ind1.I5) # 0.67025
rm(mod.sar.ind1.I5)

mod.sar.ind1.I6 <- update(mod.sar.ind1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind1.I6)
anova(mod.sar.ind1.I, mod.sar.ind1.I6) # 0.7565
rm(mod.sar.ind1.I6)

mod.sar.ind1.I7 <- update(mod.sar.ind1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind1.I7)
anova(mod.sar.ind1.I, mod.sar.ind1.I7) # 0.65098
rm(mod.sar.ind1.I7)

# significant at the 0.05 level
mod.sar.ind1.I9 <- update(mod.sar.ind1.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.ind1.I9)
anova(mod.sar.ind1.I, mod.sar.ind1.I9) # 0.012566
rm(mod.sar.ind1.I9)

mod.sar.ind1.I10 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind1.I10)
anova(mod.sar.ind1.I, mod.sar.ind1.I10) # 0.18451
rm(mod.sar.ind1.I10)

mod.sar.ind1.I11 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind1.I11)
anova(mod.sar.ind1.I, mod.sar.ind1.I11) # 0.82692
rm(mod.sar.ind1.I11)

mod.sar.ind1.I12 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind1.I12)
anova(mod.sar.ind1.I, mod.sar.ind1.I12) #  0.63618
rm(mod.sar.ind1.I12)

mod.sar.ind1.I13 <- update(mod.sar.ind1.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.ind1.I13)
anova(mod.sar.ind1.I, mod.sar.ind1.I13) # 0.062204
rm(mod.sar.ind1.I13)

mod.sar.ind1.I14 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind1.I14)
anova(mod.sar.ind1.I, mod.sar.ind1.I14) # 0.6673
rm(mod.sar.ind1.I14)

mod.sar.ind1.I15 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind1.I15)
anova(mod.sar.ind1.I, mod.sar.ind1.I15) # 0.059822
rm(mod.sar.ind1.I15)

mod.sar.ind1.I16 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind1.I16)
anova(mod.sar.ind1.I, mod.sar.ind1.I16) # 0.73619
rm(mod.sar.ind1.I16)

mod.sar.ind1.I18 <- update(mod.sar.ind1.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.ind1.I18)
anova(mod.sar.ind1.I, mod.sar.ind1.I18) # 0.13316
rm(mod.sar.ind1.I18)

mod.sar.ind1.I19 <- update(mod.sar.ind1.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.ind1.I19)
anova(mod.sar.ind1.I, mod.sar.ind1.I19) # 0.082024
rm(mod.sar.ind1.I19)

mod.sar.ind1.I21 <- update(mod.sar.ind1.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.ind1.I21)
anova(mod.sar.ind1.I, mod.sar.ind1.I21) # 0.34191
rm(mod.sar.ind1.I21)

mod.sar.ind1.I22 <- update(mod.sar.ind1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind1.I22)
anova(mod.sar.ind1.I, mod.sar.ind1.I22) # 0.91231
rm(mod.sar.ind1.I22)

mod.sar.ind1.I23 <- update(mod.sar.ind1.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.ind1.I23)
anova(mod.sar.ind1.I, mod.sar.ind1.I23) # 0.15113
rm(mod.sar.ind1.I23)

mod.sar.ind1.I24 <- update(mod.sar.ind1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind1.I24)
anova(mod.sar.ind1.I, mod.sar.ind1.I24) # 0.80927
rm(mod.sar.ind1.I24)

mod.sar.ind1.I25 <- update(mod.sar.ind1.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.ind1.I25)
anova(mod.sar.ind1.I, mod.sar.ind1.I25) # 0.4924
rm(mod.sar.ind1.I25)

mod.sar.ind1.I27 <- update(mod.sar.ind1.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.ind1.I27)
anova(mod.sar.ind1.I, mod.sar.ind1.I27) # 0.67617
rm(mod.sar.ind1.I27)

mod.sar.ind1.I28 <- update(mod.sar.ind1.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.ind1.I28)
anova(mod.sar.ind1.I, mod.sar.ind1.I28) # 0.1948
rm(mod.sar.ind1.I28)

mod.sar.ind1.I30 <- update(mod.sar.ind1.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.ind1.I30)
anova(mod.sar.ind1.I, mod.sar.ind1.I30) # 0.22831
rm(mod.sar.ind1.I30)

mod.sar.ind1.I31 <- update(mod.sar.ind1.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.ind1.I31)
anova(mod.sar.ind1.I, mod.sar.ind1.I31) # 0.10392
rm(mod.sar.ind1.I31)


# poly(meanSST.1deg, 3):I(mean.mld.t/10) is significant at the 0.05 level (as is delta_carb_ion:poly(meanSST.1deg, 3), so add it in and recheck
# 
mod.sar.ind2.I <- mod.sar.ind1.I3
rm(mod.sar.ind1.I3)

mod.sar.ind2.I2 <- update(mod.sar.ind2.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind2.I2)
anova(mod.sar.ind2.I, mod.sar.ind2.I2) # 0.39787
rm(mod.sar.ind2.I2)

mod.sar.ind2.I4 <- update(mod.sar.ind2.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind2.I4)
anova(mod.sar.ind2.I, mod.sar.ind2.I4) # 0.76191
rm(mod.sar.ind2.I4)

mod.sar.ind2.I5 <- update(mod.sar.ind2.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind2.I5)
anova(mod.sar.ind2.I, mod.sar.ind2.I5) # 0.42822
rm(mod.sar.ind2.I5)

mod.sar.ind2.I6 <- update(mod.sar.ind2.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind2.I6)
anova(mod.sar.ind2.I, mod.sar.ind2.I6) # 0.97361
rm(mod.sar.ind2.I6)

mod.sar.ind2.I7 <- update(mod.sar.ind2.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind2.I7)
anova(mod.sar.ind2.I, mod.sar.ind2.I7) # 0.54549
rm(mod.sar.ind2.I7)

# significant at the 0.05 level
mod.sar.ind2.I9 <- update(mod.sar.ind2.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.ind2.I9)
anova(mod.sar.ind2.I, mod.sar.ind2.I9) # 0.023986

mod.sar.ind2.I10 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind2.I10)
anova(mod.sar.ind2.I, mod.sar.ind2.I10) # 0.14827
rm(mod.sar.ind2.I10)

mod.sar.ind2.I11 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind2.I11)
anova(mod.sar.ind2.I, mod.sar.ind2.I11) # 0.99448
rm(mod.sar.ind2.I11)

mod.sar.ind2.I12 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind2.I12)
anova(mod.sar.ind2.I, mod.sar.ind2.I12) # 0.9036 
rm(mod.sar.ind2.I12)

mod.sar.ind2.I13 <- update(mod.sar.ind2.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.ind2.I13)
anova(mod.sar.ind2.I, mod.sar.ind2.I13) # 0.1692
rm(mod.sar.ind2.I13)

mod.sar.ind2.I14 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind2.I14)
anova(mod.sar.ind2.I, mod.sar.ind2.I14) # 0.76513
rm(mod.sar.ind2.I14)

mod.sar.ind2.I15 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind2.I15)
anova(mod.sar.ind2.I, mod.sar.ind2.I15) # 0.79312
rm(mod.sar.ind2.I15)

mod.sar.ind2.I16 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind2.I16)
anova(mod.sar.ind2.I, mod.sar.ind2.I16) # 0.8453
rm(mod.sar.ind2.I16)

mod.sar.ind2.I18 <- update(mod.sar.ind2.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.ind2.I18)
anova(mod.sar.ind2.I, mod.sar.ind2.I18) # 0.3734
rm(mod.sar.ind2.I18)

# significant at the 0.05 level
mod.sar.ind2.I19 <- update(mod.sar.ind2.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.ind2.I19)
anova(mod.sar.ind2.I, mod.sar.ind2.I19) # 0.042127
rm(mod.sar.ind2.I19)

mod.sar.ind2.I21 <- update(mod.sar.ind2.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.ind2.I21)
anova(mod.sar.ind2.I, mod.sar.ind2.I21) # 0.43901
rm(mod.sar.ind2.I21)

mod.sar.ind2.I22 <- update(mod.sar.ind2.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind2.I22)
anova(mod.sar.ind2.I, mod.sar.ind2.I22) # 0.87582
rm(mod.sar.ind2.I22)

mod.sar.ind2.I23 <- update(mod.sar.ind2.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.ind2.I23)
anova(mod.sar.ind2.I, mod.sar.ind2.I23) # 0.14845
rm(mod.sar.ind2.I23)

mod.sar.ind2.I24 <- update(mod.sar.ind2.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind2.I24)
anova(mod.sar.ind2.I, mod.sar.ind2.I24) # 0.40109
rm(mod.sar.ind2.I24)

mod.sar.ind2.I25 <- update(mod.sar.ind2.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.ind2.I25)
anova(mod.sar.ind2.I, mod.sar.ind2.I25) # 0.17482
rm(mod.sar.ind2.I25)

mod.sar.ind2.I27 <- update(mod.sar.ind2.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.ind2.I27)
anova(mod.sar.ind2.I, mod.sar.ind2.I27) # 0.7013
rm(mod.sar.ind2.I27)

mod.sar.ind2.I28 <- update(mod.sar.ind2.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.ind2.I28)
anova(mod.sar.ind2.I, mod.sar.ind2.I28) # 0.96966
rm(mod.sar.ind2.I28)

mod.sar.ind2.I30 <- update(mod.sar.ind2.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.ind2.I30)
anova(mod.sar.ind2.I, mod.sar.ind2.I30) # 0.52635
rm(mod.sar.ind2.I30)

mod.sar.ind2.I31 <- update(mod.sar.ind2.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.ind2.I31)
anova(mod.sar.ind2.I, mod.sar.ind2.I31) # 0.28305
rm(mod.sar.ind2.I31)


# delta_carb_ion:poly(meanSST.1deg, 3) is significant at the 0.05 level (as is I(depth10deg/100):absMnSal.0m, so add it in and recheck
# 
mod.sar.ind3.I <- mod.sar.ind2.I9
rm(mod.sar.ind2.I9)

mod.sar.ind3.I2 <- update(mod.sar.ind3.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind3.I2)
anova(mod.sar.ind3.I, mod.sar.ind3.I2) # 0.31159
rm(mod.sar.ind3.I2)

mod.sar.ind3.I4 <- update(mod.sar.ind3.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind3.I4)
anova(mod.sar.ind3.I, mod.sar.ind3.I4) # 0.55417
rm(mod.sar.ind3.I4)

mod.sar.ind3.I5 <- update(mod.sar.ind3.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind3.I5)
anova(mod.sar.ind3.I, mod.sar.ind3.I5) # 0.22841
rm(mod.sar.ind3.I5)

mod.sar.ind3.I6 <- update(mod.sar.ind3.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind3.I6)
anova(mod.sar.ind3.I, mod.sar.ind3.I6) # 0.75827
rm(mod.sar.ind3.I6)

mod.sar.ind3.I7 <- update(mod.sar.ind3.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind3.I7)
anova(mod.sar.ind3.I, mod.sar.ind3.I7) # 0.46962
rm(mod.sar.ind3.I7)

mod.sar.ind3.I10 <- update(mod.sar.ind3.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind3.I10)
anova(mod.sar.ind3.I, mod.sar.ind3.I10) # 0.078003
rm(mod.sar.ind3.I10)

mod.sar.ind3.I11 <- update(mod.sar.ind3.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind3.I11)
anova(mod.sar.ind3.I, mod.sar.ind3.I11) # 0.78386
rm(mod.sar.ind3.I11)

mod.sar.ind3.I12 <- update(mod.sar.ind3.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind3.I12)
anova(mod.sar.ind3.I, mod.sar.ind3.I12) #  0.96801
rm(mod.sar.ind3.I12)

mod.sar.ind3.I13 <- update(mod.sar.ind3.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.ind3.I13)
anova(mod.sar.ind3.I, mod.sar.ind3.I13) # 0.9163
rm(mod.sar.ind3.I13)

mod.sar.ind3.I14 <- update(mod.sar.ind3.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind3.I14)
anova(mod.sar.ind3.I, mod.sar.ind3.I14) # 0.63716
rm(mod.sar.ind3.I14)

mod.sar.ind3.I15 <- update(mod.sar.ind3.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind3.I15)
anova(mod.sar.ind3.I, mod.sar.ind3.I15) # 0.8968
rm(mod.sar.ind3.I15)

mod.sar.ind3.I16 <- update(mod.sar.ind3.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind3.I16)
anova(mod.sar.ind3.I, mod.sar.ind3.I16) # 0.99215
rm(mod.sar.ind3.I16)

mod.sar.ind3.I18 <- update(mod.sar.ind3.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.ind3.I18)
anova(mod.sar.ind3.I, mod.sar.ind3.I18) # 0.53561
rm(mod.sar.ind3.I18)

# significant at the 0.05 level
mod.sar.ind3.I19 <- update(mod.sar.ind3.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.ind3.I19)
anova(mod.sar.ind3.I, mod.sar.ind3.I19) # 0.048813

mod.sar.ind3.I21 <- update(mod.sar.ind3.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.ind3.I21)
anova(mod.sar.ind3.I, mod.sar.ind3.I21) # 0.8154
rm(mod.sar.ind3.I21)

mod.sar.ind3.I22 <- update(mod.sar.ind3.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind3.I22)
anova(mod.sar.ind3.I, mod.sar.ind3.I22) # 0.75453
rm(mod.sar.ind3.I22)

mod.sar.ind3.I23 <- update(mod.sar.ind3.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.ind3.I23)
anova(mod.sar.ind3.I, mod.sar.ind3.I23) # 0.26439
rm(mod.sar.ind3.I23)

mod.sar.ind3.I24 <- update(mod.sar.ind3.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind3.I24)
anova(mod.sar.ind3.I, mod.sar.ind3.I24) # 0.40043
rm(mod.sar.ind3.I24)

mod.sar.ind3.I25 <- update(mod.sar.ind3.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.ind3.I25)
anova(mod.sar.ind3.I, mod.sar.ind3.I25) # 0.27059
rm(mod.sar.ind3.I25)

mod.sar.ind3.I27 <- update(mod.sar.ind3.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.ind3.I27)
anova(mod.sar.ind3.I, mod.sar.ind3.I27) # 0.8925
rm(mod.sar.ind3.I27)

mod.sar.ind3.I28 <- update(mod.sar.ind3.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.ind3.I28)
anova(mod.sar.ind3.I, mod.sar.ind3.I28) # 0.883
rm(mod.sar.ind3.I28)

mod.sar.ind3.I30 <- update(mod.sar.ind3.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.ind3.I30)
anova(mod.sar.ind3.I, mod.sar.ind3.I30) # 0.93081
rm(mod.sar.ind3.I30)

mod.sar.ind3.I31 <- update(mod.sar.ind3.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.ind3.I31)
anova(mod.sar.ind3.I, mod.sar.ind3.I31) # 0.58186
rm(mod.sar.ind3.I31)


# I(depth10deg/100):absMnSal.0m is significant at the 0.05 level so add it in and recheck
# 
mod.sar.ind4.I <- mod.sar.ind3.I19
rm(mod.sar.ind3.I19)

mod.sar.ind4.I2 <- update(mod.sar.ind4.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.ind4.I2)
anova(mod.sar.ind4.I, mod.sar.ind4.I2) # 0.46367
rm(mod.sar.ind4.I2)

mod.sar.ind4.I4 <- update(mod.sar.ind4.I, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.ind4.I4)
anova(mod.sar.ind4.I, mod.sar.ind4.I4) # 0.31166
rm(mod.sar.ind4.I4)

mod.sar.ind4.I5 <- update(mod.sar.ind4.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.ind4.I5)
anova(mod.sar.ind4.I, mod.sar.ind4.I5) # 0.18473
rm(mod.sar.ind4.I5)

mod.sar.ind4.I6 <- update(mod.sar.ind4.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.ind4.I6)
anova(mod.sar.ind4.I, mod.sar.ind4.I6) # 0.8312
rm(mod.sar.ind4.I6)

mod.sar.ind4.I7 <- update(mod.sar.ind4.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.ind4.I7)
anova(mod.sar.ind4.I, mod.sar.ind4.I7) # 0.3646
rm(mod.sar.ind4.I7)

mod.sar.ind4.I10 <- update(mod.sar.ind4.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.ind4.I10)
anova(mod.sar.ind4.I, mod.sar.ind4.I10) # 0.17899
rm(mod.sar.ind4.I10)

mod.sar.ind4.I11 <- update(mod.sar.ind4.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.ind4.I11)
anova(mod.sar.ind4.I, mod.sar.ind4.I11) # 0.92685
rm(mod.sar.ind4.I11)

mod.sar.ind4.I12 <- update(mod.sar.ind4.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.ind4.I12)
anova(mod.sar.ind4.I, mod.sar.ind4.I12) #  0.97611
rm(mod.sar.ind4.I12)

mod.sar.ind4.I13 <- update(mod.sar.ind4.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.ind4.I13)
anova(mod.sar.ind4.I, mod.sar.ind4.I13) # 0.91465
rm(mod.sar.ind4.I13)

mod.sar.ind4.I14 <- update(mod.sar.ind4.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.ind4.I14)
anova(mod.sar.ind4.I, mod.sar.ind4.I14) # 0.38342
rm(mod.sar.ind4.I14)

mod.sar.ind4.I15 <- update(mod.sar.ind4.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.ind4.I15)
anova(mod.sar.ind4.I, mod.sar.ind4.I15) # 0.77137
rm(mod.sar.ind4.I15)

mod.sar.ind4.I16 <- update(mod.sar.ind4.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.ind4.I16)
anova(mod.sar.ind4.I, mod.sar.ind4.I16) # 0.98177
rm(mod.sar.ind4.I16)

mod.sar.ind4.I18 <- update(mod.sar.ind4.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.ind4.I18)
anova(mod.sar.ind4.I, mod.sar.ind4.I18) # 0.53644
rm(mod.sar.ind4.I18)

mod.sar.ind4.I21 <- update(mod.sar.ind4.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.ind4.I21)
anova(mod.sar.ind4.I, mod.sar.ind4.I21) # 0.8847
rm(mod.sar.ind4.I21)

mod.sar.ind4.I22 <- update(mod.sar.ind4.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.ind4.I22)
anova(mod.sar.ind4.I, mod.sar.ind4.I22) # 0.84796
rm(mod.sar.ind4.I22)

mod.sar.ind4.I23 <- update(mod.sar.ind4.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.ind4.I23)
anova(mod.sar.ind4.I, mod.sar.ind4.I23) # 0.45701
rm(mod.sar.ind4.I23)

mod.sar.ind4.I24 <- update(mod.sar.ind4.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.ind4.I24)
anova(mod.sar.ind4.I, mod.sar.ind4.I24) # 0.47642
rm(mod.sar.ind4.I24)

mod.sar.ind4.I25 <- update(mod.sar.ind4.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.ind4.I25)
anova(mod.sar.ind4.I, mod.sar.ind4.I25) # 0.62473
rm(mod.sar.ind4.I25)

mod.sar.ind4.I27 <- update(mod.sar.ind4.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.ind4.I27)
anova(mod.sar.ind4.I, mod.sar.ind4.I27) # 0.93019
rm(mod.sar.ind4.I27)

mod.sar.ind4.I28 <- update(mod.sar.ind4.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.ind4.I28)
anova(mod.sar.ind4.I, mod.sar.ind4.I28) # 0.94016
rm(mod.sar.ind4.I28)

mod.sar.ind4.I30 <- update(mod.sar.ind4.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.ind4.I30)
anova(mod.sar.ind4.I, mod.sar.ind4.I30) # 0.96117
rm(mod.sar.ind4.I30)

mod.sar.ind4.I31 <- update(mod.sar.ind4.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.ind4.I31)
anova(mod.sar.ind4.I, mod.sar.ind4.I31) # 0.85767
rm(mod.sar.ind4.I31)

# Nothing more is significant

mod.sar.indIf <- mod.sar.ind4.I
summary(mod.sar.indIf, Nagelkerke = TRUE) # 0.86535
AIC(mod.sar.indIf) # 1033.724

# Check for correlation
env.var.ind <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion")

# pairs plot
png("Figures/Ana_3viii_ind_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])
# correlation between SST and MLD

lr.sar.indIf <- lr.calc(mod.sar.indIf)

png("Figures/Ana_3viii_LRatio_indIf.png", width = 800)
lr.plot(lr.sar.indIf, order = c(7:9, 10:11, 1, 4, 2:3, 5:6), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.indIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.indIfg <- lr.calc(mod.sar.indIf, tmp)
rm(tmp)

png("Figures/Ana_3viii_LRatio_g_indIf.png", width = 800)
lr.plot(lr.sar.indIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

save(mod.sar.indI, mod.sar.ind1.I, mod.sar.ind2.I, mod.sar.ind3.I, mod.sar.ind4.I, file = "Outputs/Indian_simplification.RData")
rm(ind.s, env.var.ind, mod.sar.indI, mod.sar.ind1.I, mod.sar.ind2.I, mod.sar.ind3.I, mod.sar.ind4.I)

## 3ix. Only significant interactions for Pacific --------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.pacI <- errorsarlm(op.formula, listw = pac.s, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", ])
par(mfrow = c(1, 2))
sar.plot(mod.sar.pacI)
par(mfrow = c(1, 1))

# try adding in the other interactions
mod.sar.pacI2 <- update(mod.sar.pacI, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI2)
anova(mod.sar.pacI, mod.sar.pacI2) # 0.49757
rm(mod.sar.pacI2)

# significant at the 0.05 level
mod.sar.pacI3 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.pacI3)
anova(mod.sar.pacI, mod.sar.pacI3) # 0.0054992
rm(mod.sar.pacI3)

# significant at the 0.05 level
mod.sar.pacI4 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):I(depth10deg/100))
summary(mod.sar.pacI4)
anova(mod.sar.pacI, mod.sar.pacI4) # 0.0048417

mod.sar.pacI5 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pacI5)
anova(mod.sar.pacI, mod.sar.pacI5) # 0.20129
rm(mod.sar.pacI5)

mod.sar.pacI6 <- update(mod.sar.pacI, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pacI6)
anova(mod.sar.pacI, mod.sar.pacI6) # 0.4505
rm(mod.sar.pacI6)

mod.sar.pacI7 <- update(mod.sar.pacI, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pacI7)
anova(mod.sar.pacI, mod.sar.pacI7) # 0.97162
rm(mod.sar.pacI7)

mod.sar.pacI9 <- update(mod.sar.pacI, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI9)
anova(mod.sar.pacI, mod.sar.pacI9) # 0.72598
rm(mod.sar.pacI9)

mod.sar.pacI10 <- update(mod.sar.pacI, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pacI10)
anova(mod.sar.pacI, mod.sar.pacI10) # 0.19835
rm(mod.sar.pacI10)

mod.sar.pacI11 <- update(mod.sar.pacI, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pacI11)
anova(mod.sar.pacI, mod.sar.pacI11) # 0.12939
rm(mod.sar.pacI11)

# significant at the 0.05 level
mod.sar.pacI12 <- update(mod.sar.pacI, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pacI12)
anova(mod.sar.pacI, mod.sar.pacI12) # 0.031738
rm(mod.sar.pacI12)

mod.sar.pacI13 <- update(mod.sar.pacI, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.pacI13)
anova(mod.sar.pacI, mod.sar.pacI13) # 0.064652
rm(mod.sar.pacI13)

mod.sar.pacI14 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pacI14)
anova(mod.sar.pacI, mod.sar.pacI14) # 0.74371
rm(mod.sar.pacI14)

mod.sar.pacI15 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pacI15)
anova(mod.sar.pacI, mod.sar.pacI15) # 0.21602
rm(mod.sar.pacI15)

mod.sar.pacI16 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pacI16)
anova(mod.sar.pacI, mod.sar.pacI16) # 0.36686
rm(mod.sar.pacI16)

mod.sar.pacI17 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pacI17)
anova(mod.sar.pacI, mod.sar.pacI17) # 0.34832
rm(mod.sar.pacI17)

# significant at the 0.05 level
mod.sar.pacI18 <- update(mod.sar.pacI, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.pacI18)
anova(mod.sar.pacI, mod.sar.pacI18) # 0.015415
rm(mod.sar.pacI18)

mod.sar.pacI19 <- update(mod.sar.pacI, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.pacI19)
anova(mod.sar.pacI, mod.sar.pacI19) # 0.092581
rm(mod.sar.pacI19)

mod.sar.pacI21 <- update(mod.sar.pacI, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.pacI21)
anova(mod.sar.pacI, mod.sar.pacI21) # 0.13341
rm(mod.sar.pacI21)

mod.sar.pacI22 <- update(mod.sar.pacI, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pacI22)
anova(mod.sar.pacI, mod.sar.pacI22) # 0.28196
rm(mod.sar.pacI22)

mod.sar.pacI23 <- update(mod.sar.pacI, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.pacI23)
anova(mod.sar.pacI, mod.sar.pacI23) # 0.37956
rm(mod.sar.pacI23)

mod.sar.pacI24 <- update(mod.sar.pacI, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pacI24)
anova(mod.sar.pacI, mod.sar.pacI24) # 0.9209
rm(mod.sar.pacI24)

mod.sar.pacI25 <- update(mod.sar.pacI, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.pacI25)
anova(mod.sar.pacI, mod.sar.pacI25) # 0.40176
rm(mod.sar.pacI25)

mod.sar.pacI27 <- update(mod.sar.pacI, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.pacI27)
anova(mod.sar.pacI, mod.sar.pacI27) # 0.81533
rm(mod.sar.pacI27)

mod.sar.pacI28 <- update(mod.sar.pacI, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.pacI28)
anova(mod.sar.pacI, mod.sar.pacI28) # 0.69828
rm(mod.sar.pacI28)

mod.sar.pacI30 <- update(mod.sar.pacI, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.pacI30)
anova(mod.sar.pacI, mod.sar.pacI30) # 0.63095
rm(mod.sar.pacI30)

mod.sar.pacI31 <- update(mod.sar.pacI, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.pacI31)
anova(mod.sar.pacI, mod.sar.pacI31) # 0.46578
rm(mod.sar.pacI31)


# multiple interactions were significant. Add in the most significant and then recalculate
mod.sar.pac1.I <- mod.sar.pacI4
rm(mod.sar.pacI4)

# try adding in the other interactions
mod.sar.pac1.I2 <- update(mod.sar.pac1.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac1.I2)
anova(mod.sar.pac1.I, mod.sar.pac1.I2) # 0.92108
rm(mod.sar.pac1.I2)

# significant at the 0.05 level
mod.sar.pac1.I3 <- update(mod.sar.pac1.I, ~. + poly(meanSST.1deg, 3):I(mean.mld.t/10))
summary(mod.sar.pac1.I3)
anova(mod.sar.pac1.I, mod.sar.pac1.I3) # 0.0078896

mod.sar.pac1.I5 <- update(mod.sar.pac1.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac1.I5)
anova(mod.sar.pac1.I, mod.sar.pac1.I5) # 0.62158
rm(mod.sar.pac1.I5)

mod.sar.pac1.I6 <- update(mod.sar.pac1.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac1.I6)
anova(mod.sar.pac1.I, mod.sar.pac1.I6) # 0.93069
rm(mod.sar.pac1.I6)

mod.sar.pac1.I7 <- update(mod.sar.pac1.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac1.I7)
anova(mod.sar.pac1.I, mod.sar.pac1.I7) # 0.87748
rm(mod.sar.pac1.I7)

mod.sar.pac1.I9 <- update(mod.sar.pac1.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pac1.I9)
anova(mod.sar.pac1.I, mod.sar.pac1.I9) # 0.51157 
rm(mod.sar.pac1.I9)

mod.sar.pac1.I10 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac1.I10)
anova(mod.sar.pac1.I, mod.sar.pac1.I10) # 0.52793
rm(mod.sar.pac1.I10)

mod.sar.pac1.I11 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac1.I11)
anova(mod.sar.pac1.I, mod.sar.pac1.I11) # 0.11369
rm(mod.sar.pac1.I11)

# significant at the 0.05 level
mod.sar.pac1.I12 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac1.I12)
anova(mod.sar.pac1.I, mod.sar.pac1.I12) # 0.18555
rm(mod.sar.pac1.I12)

# significant at the 0.05 level
mod.sar.pac1.I13 <- update(mod.sar.pac1.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.pac1.I13)
anova(mod.sar.pac1.I, mod.sar.pac1.I13) # 0.040274
rm(mod.sar.pac1.I13)

mod.sar.pac1.I14 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac1.I14)
anova(mod.sar.pac1.I, mod.sar.pac1.I14) # 0.52655
rm(mod.sar.pac1.I14)

mod.sar.pac1.I15 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac1.I15)
anova(mod.sar.pac1.I, mod.sar.pac1.I15) # 0.87029
rm(mod.sar.pac1.I15)

mod.sar.pac1.I16 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac1.I16)
anova(mod.sar.pac1.I, mod.sar.pac1.I16) # 0.33075
rm(mod.sar.pac1.I16)

mod.sar.pac1.I17 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac1.I17)
anova(mod.sar.pac1.I, mod.sar.pac1.I17) # 0.43237
rm(mod.sar.pac1.I17)

# significant at the 0.05 level
mod.sar.pac1.I18 <- update(mod.sar.pac1.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.pac1.I18)
anova(mod.sar.pac1.I, mod.sar.pac1.I18) # 0.01561
rm(mod.sar.pac1.I18)

mod.sar.pac1.I19 <- update(mod.sar.pac1.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.pac1.I19)
anova(mod.sar.pac1.I, mod.sar.pac1.I19) # 0.067521
rm(mod.sar.pac1.I19)

mod.sar.pac1.I21 <- update(mod.sar.pac1.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.pac1.I21)
anova(mod.sar.pac1.I, mod.sar.pac1.I21) # 0.2749
rm(mod.sar.pac1.I21)

mod.sar.pac1.I22 <- update(mod.sar.pac1.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac1.I22)
anova(mod.sar.pac1.I, mod.sar.pac1.I22) # 0.20087
rm(mod.sar.pac1.I22)

mod.sar.pac1.I23 <- update(mod.sar.pac1.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.pac1.I23)
anova(mod.sar.pac1.I, mod.sar.pac1.I23) # 0.20485
rm(mod.sar.pac1.I23)

mod.sar.pac1.I24 <- update(mod.sar.pac1.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac1.I24)
anova(mod.sar.pac1.I, mod.sar.pac1.I24) # 0.74585
rm(mod.sar.pac1.I24)

mod.sar.pac1.I25 <- update(mod.sar.pac1.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.pac1.I25)
anova(mod.sar.pac1.I, mod.sar.pac1.I25) # 0.57501
rm(mod.sar.pac1.I25)

mod.sar.pac1.I27 <- update(mod.sar.pac1.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.pac1.I27)
anova(mod.sar.pac1.I, mod.sar.pac1.I27) # 0.98376
rm(mod.sar.pac1.I27)

mod.sar.pac1.I28 <- update(mod.sar.pac1.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.pac1.I28)
anova(mod.sar.pac1.I, mod.sar.pac1.I28) # 0.55538
rm(mod.sar.pac1.I28)

mod.sar.pac1.I30 <- update(mod.sar.pac1.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.pac1.I30)
anova(mod.sar.pac1.I, mod.sar.pac1.I30) # 0.86085
rm(mod.sar.pac1.I30)

mod.sar.pac1.I31 <- update(mod.sar.pac1.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.pac1.I31)
anova(mod.sar.pac1.I, mod.sar.pac1.I31) # 0.45622
rm(mod.sar.pac1.I31)


# add in poly(meanSST.1deg, 3):mean.mld.t and recalculate
mod.sar.pac2.I <- mod.sar.pac1.I3
rm(mod.sar.pac1.I3)

# try adding in the other interactions
mod.sar.pac2.I2 <- update(mod.sar.pac2.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac2.I2)
anova(mod.sar.pac2.I, mod.sar.pac2.I2) # 0.47797
rm(mod.sar.pac2.I2)

mod.sar.pac2.I5 <- update(mod.sar.pac2.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac2.I5)
anova(mod.sar.pac2.I, mod.sar.pac2.I5) # 0.91649
rm(mod.sar.pac2.I5)

mod.sar.pac2.I6 <- update(mod.sar.pac2.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac2.I6)
anova(mod.sar.pac2.I, mod.sar.pac2.I6) # 0.47812
rm(mod.sar.pac2.I6)

mod.sar.pac2.I7 <- update(mod.sar.pac2.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac2.I7)
anova(mod.sar.pac2.I, mod.sar.pac2.I7) # 0.57602
rm(mod.sar.pac2.I7)

mod.sar.pac2.I9 <- update(mod.sar.pac2.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pac2.I9)
anova(mod.sar.pac2.I, mod.sar.pac2.I9) #  0.88527
rm(mod.sar.pac2.I9)

mod.sar.pac2.I10 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac2.I10)
anova(mod.sar.pac2.I, mod.sar.pac2.I10) # 0.55538
rm(mod.sar.pac2.I10)

# significant at the 0.05 level
mod.sar.pac2.I11 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:logProd.mn.ann)
summary(mod.sar.pac2.I11)
anova(mod.sar.pac2.I, mod.sar.pac2.I11) # 0.01735

mod.sar.pac2.I12 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac2.I12)
anova(mod.sar.pac2.I, mod.sar.pac2.I12) # 0.31569
rm(mod.sar.pac2.I12)

mod.sar.pac2.I13 <- update(mod.sar.pac2.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.pac2.I13)
anova(mod.sar.pac2.I, mod.sar.pac2.I13) # 0.054336
rm(mod.sar.pac2.I13)

mod.sar.pac2.I14 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac2.I14)
anova(mod.sar.pac2.I, mod.sar.pac2.I14) # 0.69515
rm(mod.sar.pac2.I14)

mod.sar.pac2.I15 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac2.I15)
anova(mod.sar.pac2.I, mod.sar.pac2.I15) # 0.60551
rm(mod.sar.pac2.I15)

mod.sar.pac2.I16 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac2.I16)
anova(mod.sar.pac2.I, mod.sar.pac2.I16) # 0.63653
rm(mod.sar.pac2.I16)

mod.sar.pac2.I17 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac2.I17)
anova(mod.sar.pac2.I, mod.sar.pac2.I17) # 0.35766
rm(mod.sar.pac2.I17)

mod.sar.pac2.I18 <- update(mod.sar.pac2.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.pac2.I18)
anova(mod.sar.pac2.I, mod.sar.pac2.I18) # 0.1286
rm(mod.sar.pac2.I18)

mod.sar.pac2.I19 <- update(mod.sar.pac2.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.pac2.I19)
anova(mod.sar.pac2.I, mod.sar.pac2.I19) # 0.13658
rm(mod.sar.pac2.I19)

mod.sar.pac2.I21 <- update(mod.sar.pac2.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.pac2.I21)
anova(mod.sar.pac2.I, mod.sar.pac2.I21) # 0.40787
rm(mod.sar.pac2.I21)

mod.sar.pac2.I22 <- update(mod.sar.pac2.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac2.I22)
anova(mod.sar.pac2.I, mod.sar.pac2.I22) # 0.086263
rm(mod.sar.pac2.I22)

mod.sar.pac2.I23 <- update(mod.sar.pac2.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.pac2.I23)
anova(mod.sar.pac2.I, mod.sar.pac2.I23) # 0.15135
rm(mod.sar.pac2.I23)

mod.sar.pac2.I24 <- update(mod.sar.pac2.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac2.I24)
anova(mod.sar.pac2.I, mod.sar.pac2.I24) # 0.34378
rm(mod.sar.pac2.I24)

mod.sar.pac2.I25 <- update(mod.sar.pac2.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.pac2.I25)
anova(mod.sar.pac2.I, mod.sar.pac2.I25) # 0.91648
rm(mod.sar.pac2.I25)

mod.sar.pac2.I27 <- update(mod.sar.pac2.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.pac2.I27)
anova(mod.sar.pac2.I, mod.sar.pac2.I27) # 0.98268
rm(mod.sar.pac2.I27)

mod.sar.pac2.I28 <- update(mod.sar.pac2.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.pac2.I28)
anova(mod.sar.pac2.I, mod.sar.pac2.I28) # 0.27949
rm(mod.sar.pac2.I28)

mod.sar.pac2.I30 <- update(mod.sar.pac2.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.pac2.I30)
anova(mod.sar.pac2.I, mod.sar.pac2.I30) # 0.76627
rm(mod.sar.pac2.I30)

mod.sar.pac2.I31 <- update(mod.sar.pac2.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.pac2.I31)
anova(mod.sar.pac2.I, mod.sar.pac2.I31) # 0.32665
rm(mod.sar.pac2.I31)

# add in sdSST.1deg:logProd.mn.ann and recalculate
mod.sar.pac3.I <- mod.sar.pac2.I11
rm(mod.sar.pac2.I11)

# try adding in the other interactions
mod.sar.pac3.I2 <- update(mod.sar.pac3.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac3.I2)
anova(mod.sar.pac3.I, mod.sar.pac3.I2) # 0.91224
rm(mod.sar.pac3.I2)

mod.sar.pac3.I5 <- update(mod.sar.pac3.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac3.I5)
anova(mod.sar.pac3.I, mod.sar.pac3.I5) # 0.25957
rm(mod.sar.pac3.I5)

mod.sar.pac3.I6 <- update(mod.sar.pac3.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac3.I6)
anova(mod.sar.pac3.I, mod.sar.pac3.I6) # 0.5501
rm(mod.sar.pac3.I6)

mod.sar.pac3.I7 <- update(mod.sar.pac3.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac3.I7)
anova(mod.sar.pac3.I, mod.sar.pac3.I7) # 0.75835
rm(mod.sar.pac3.I7)

mod.sar.pac3.I9 <- update(mod.sar.pac3.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pac3.I9)
anova(mod.sar.pac3.I, mod.sar.pac3.I9) # 0.8027 
rm(mod.sar.pac3.I9)

mod.sar.pac3.I10 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac3.I10)
anova(mod.sar.pac3.I, mod.sar.pac3.I10) # 0.89013
rm(mod.sar.pac3.I10)

mod.sar.pac3.I12 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac3.I12)
anova(mod.sar.pac3.I, mod.sar.pac3.I12) # 0.19926
rm(mod.sar.pac3.I12)

# significant at the 0.05 level
mod.sar.pac3.I13 <- update(mod.sar.pac3.I, ~. + sdSST.1deg:delta_carb_ion)
summary(mod.sar.pac3.I13)
anova(mod.sar.pac3.I, mod.sar.pac3.I13) # 0.037206

mod.sar.pac3.I14 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac3.I14)
anova(mod.sar.pac3.I, mod.sar.pac3.I14) # 0.81682
rm(mod.sar.pac3.I14)

mod.sar.pac3.I15 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac3.I15)
anova(mod.sar.pac3.I, mod.sar.pac3.I15) # 0.65358
rm(mod.sar.pac3.I15)

mod.sar.pac3.I16 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac3.I16)
anova(mod.sar.pac3.I, mod.sar.pac3.I16) # 0.91734
rm(mod.sar.pac3.I16)

mod.sar.pac3.I17 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac3.I17)
anova(mod.sar.pac3.I, mod.sar.pac3.I17) # 0.46314
rm(mod.sar.pac3.I17)

mod.sar.pac3.I18 <- update(mod.sar.pac3.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.pac3.I18)
anova(mod.sar.pac3.I, mod.sar.pac3.I18) # 0.15594
rm(mod.sar.pac3.I18)

mod.sar.pac3.I19 <- update(mod.sar.pac3.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.pac3.I19)
anova(mod.sar.pac3.I, mod.sar.pac3.I19) # 0.20299
rm(mod.sar.pac3.I19)

mod.sar.pac3.I21 <- update(mod.sar.pac3.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.pac3.I21)
anova(mod.sar.pac3.I, mod.sar.pac3.I21) # 0.34057
rm(mod.sar.pac3.I21)

mod.sar.pac3.I22 <- update(mod.sar.pac3.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac3.I22)
anova(mod.sar.pac3.I, mod.sar.pac3.I22) # 0.20734
rm(mod.sar.pac3.I22)

mod.sar.pac3.I23 <- update(mod.sar.pac3.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.pac3.I23)
anova(mod.sar.pac3.I, mod.sar.pac3.I23) # 0.28197
rm(mod.sar.pac3.I23)

mod.sar.pac3.I24 <- update(mod.sar.pac3.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac3.I24)
anova(mod.sar.pac3.I, mod.sar.pac3.I24) # 0.46697
rm(mod.sar.pac3.I24)

mod.sar.pac3.I25 <- update(mod.sar.pac3.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.pac3.I25)
anova(mod.sar.pac3.I, mod.sar.pac3.I25) # 0.688
rm(mod.sar.pac3.I25)

mod.sar.pac3.I27 <- update(mod.sar.pac3.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.pac3.I27)
anova(mod.sar.pac3.I, mod.sar.pac3.I27) # 0.99263
rm(mod.sar.pac3.I27)

mod.sar.pac3.I28 <- update(mod.sar.pac3.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.pac3.I28)
anova(mod.sar.pac3.I, mod.sar.pac3.I28) # 0.42745
rm(mod.sar.pac3.I28)

mod.sar.pac3.I30 <- update(mod.sar.pac3.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.pac3.I30)
anova(mod.sar.pac3.I, mod.sar.pac3.I30) # 0.95914
rm(mod.sar.pac3.I30)

mod.sar.pac3.I31 <- update(mod.sar.pac3.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.pac3.I31)
anova(mod.sar.pac3.I, mod.sar.pac3.I31) # 0.44131
rm(mod.sar.pac3.I31)


# add in sdSST.1deg:delta_carb_ion and recalculate
mod.sar.pac4.I <- mod.sar.pac3.I13
rm(mod.sar.pac3.I13)

# try adding in the other interactions
mod.sar.pac4.I2 <- update(mod.sar.pac4.I, ~. - poly(meanSST.1deg, 1):sdSST.1deg + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pac4.I2)
anova(mod.sar.pac4.I, mod.sar.pac4.I2) # 0.73057
rm(mod.sar.pac4.I2)

mod.sar.pac4.I5 <- update(mod.sar.pac4.I, ~. + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pac4.I5)
anova(mod.sar.pac4.I, mod.sar.pac4.I5) # 0.19794
rm(mod.sar.pac4.I5)

mod.sar.pac4.I6 <- update(mod.sar.pac4.I, ~. - absMnSal.0m:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):absMnSal.0m)
summary(mod.sar.pac4.I6)
anova(mod.sar.pac4.I, mod.sar.pac4.I6) # 0.60861
rm(mod.sar.pac4.I6)

mod.sar.pac4.I7 <- update(mod.sar.pac4.I, ~. - prop2.oxy:poly(meanSST.1deg, 2) + poly(meanSST.1deg, 3):prop2.oxy)
summary(mod.sar.pac4.I7)
anova(mod.sar.pac4.I, mod.sar.pac4.I7) # 0.76467
rm(mod.sar.pac4.I7)

mod.sar.pac4.I9 <- update(mod.sar.pac4.I, ~. - delta_carb_ion:poly(meanSST.1deg, 2) + delta_carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pac4.I9)
anova(mod.sar.pac4.I, mod.sar.pac4.I9) # 0.45163 
rm(mod.sar.pac4.I9)

mod.sar.pac4.I10 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:I(depth10deg/100))
summary(mod.sar.pac4.I10)
anova(mod.sar.pac4.I, mod.sar.pac4.I10) # 0.85379
rm(mod.sar.pac4.I10)

mod.sar.pac4.I12 <- update(mod.sar.pac4.I, ~. + sdSST.1deg:absMnSal.0m)
summary(mod.sar.pac4.I12)
anova(mod.sar.pac4.I, mod.sar.pac4.I12) # 0.18787
rm(mod.sar.pac4.I12)

mod.sar.pac4.I14 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):I(depth10deg/100))
summary(mod.sar.pac4.I14)
anova(mod.sar.pac4.I, mod.sar.pac4.I14) # 0.71541
rm(mod.sar.pac4.I14)

mod.sar.pac4.I15 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):absMnSal.0m)
summary(mod.sar.pac4.I15)
anova(mod.sar.pac4.I, mod.sar.pac4.I15) # 0.44578
rm(mod.sar.pac4.I15)

mod.sar.pac4.I16 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):sdSal.0m)
summary(mod.sar.pac4.I16)
anova(mod.sar.pac4.I, mod.sar.pac4.I16) # 0.94996
rm(mod.sar.pac4.I16)

mod.sar.pac4.I17 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):prop2.oxy)
summary(mod.sar.pac4.I17)
anova(mod.sar.pac4.I, mod.sar.pac4.I17) # 0.37635
rm(mod.sar.pac4.I17)

mod.sar.pac4.I18 <- update(mod.sar.pac4.I, ~. + I(mean.mld.t/10):delta_carb_ion)
summary(mod.sar.pac4.I18)
anova(mod.sar.pac4.I, mod.sar.pac4.I18) # 0.33318
rm(mod.sar.pac4.I18)

mod.sar.pac4.I19 <- update(mod.sar.pac4.I, ~. + I(depth10deg/100):absMnSal.0m)
summary(mod.sar.pac4.I19)
anova(mod.sar.pac4.I, mod.sar.pac4.I19) # 0.14677
rm(mod.sar.pac4.I19)

mod.sar.pac4.I21 <- update(mod.sar.pac4.I, ~. + I(depth10deg/100):delta_carb_ion)
summary(mod.sar.pac4.I21)
anova(mod.sar.pac4.I, mod.sar.pac4.I21) # 0.71721
rm(mod.sar.pac4.I21)

mod.sar.pac4.I22 <- update(mod.sar.pac4.I, ~. + logProd.mn.ann:prop2.oxy)
summary(mod.sar.pac4.I22)
anova(mod.sar.pac4.I, mod.sar.pac4.I22) # 0.21723
rm(mod.sar.pac4.I22)

mod.sar.pac4.I23 <- update(mod.sar.pac4.I, ~. + logProd.mn.ann:delta_carb_ion)
summary(mod.sar.pac4.I23)
anova(mod.sar.pac4.I, mod.sar.pac4.I23) # 0.27182
rm(mod.sar.pac4.I23)

mod.sar.pac4.I24 <- update(mod.sar.pac4.I, ~. + absMnSal.0m:sdSal.0m)
summary(mod.sar.pac4.I24)
anova(mod.sar.pac4.I, mod.sar.pac4.I24) # 0.45903
rm(mod.sar.pac4.I24)

mod.sar.pac4.I25 <- update(mod.sar.pac4.I, ~. + absMnSal.0m:prop2.oxy )
summary(mod.sar.pac4.I25)
anova(mod.sar.pac4.I, mod.sar.pac4.I25) # 0.5496
rm(mod.sar.pac4.I25)

mod.sar.pac4.I27 <- update(mod.sar.pac4.I, ~. + absMnSal.0m:delta_carb_ion)
summary(mod.sar.pac4.I27)
anova(mod.sar.pac4.I, mod.sar.pac4.I27) # 0.95035
rm(mod.sar.pac4.I27)

mod.sar.pac4.I28 <- update(mod.sar.pac4.I, ~. + sdSal.0m:prop2.oxy)
summary(mod.sar.pac4.I28)
anova(mod.sar.pac4.I, mod.sar.pac4.I28) # 0.45047
rm(mod.sar.pac4.I28)

mod.sar.pac4.I30 <- update(mod.sar.pac4.I, ~. + sdSal.0m:delta_carb_ion)
summary(mod.sar.pac4.I30)
anova(mod.sar.pac4.I, mod.sar.pac4.I30) # 0.9957
rm(mod.sar.pac4.I30)

mod.sar.pac4.I31 <- update(mod.sar.pac4.I, ~. + prop2.oxy:delta_carb_ion)
summary(mod.sar.pac4.I31)
anova(mod.sar.pac4.I, mod.sar.pac4.I31) # 0.2942
rm(mod.sar.pac4.I31)

# don't need to add anything else
mod.sar.pacIf <- mod.sar.pac4.I

summary(mod.sar.pacIf, Nagelkerke = TRUE) # 0.71599 
AIC(mod.sar.pacIf) # 2012.546

# Check for correlation
env.var.pac <- c("meanSST.1deg", "sdSST.1deg", "depth10deg", "mean.mld.t", "logProd.mn.ann", "absMnSal.0m", "sdSal.0m", "delta_carb_ion")

# pairs plot
png("Figures/Ana_3ix_pac_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])

lr.sar.pacIf <- lr.calc(mod.sar.pacIf)

png("Figures/Ana_3ix_LRatio_pacIf.png", width = 800)
lr.plot(lr.sar.pacIf, order = c(7:11, 1, 4, 2:3, 5:6), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 5)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.pacIf), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.pacIfg <- lr.calc(mod.sar.pacIf, tmp)
rm(tmp)

png("Figures/Ana_3ix_LRatio_g_pacIf.png", width = 800)
lr.plot(lr.sar.pacIfg, order = c(5:6, 1:4), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.pacI, mod.sar.pac1.I, mod.sar.pac2.I, mod.sar.pac3.I, mod.sar.pac4.I, file = "Outputs/Pacific_simplification.RData")
rm(env.var.pac, mod.sar.pacI, mod.sar.pac1.I, mod.sar.pac2.I, mod.sar.pac3.I, mod.sar.pac4.I)

## 3x. Comparison of the oceans --------------------------------------------
png("Figures/Ana_3x_LRatio_oceIf.png", width = 1200)
lr.plot(lr.sar.op0, lr.sar.atlIf, lr.sar.indIf, lr.sar.pacIf, order = c(7:9, 3:4, 12:11, 5, 1, 10, 6, 2), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 20)
dev.off()

# also for groups of variables
png("Figures/Ana_3x_LRatio_g_OceIf.png", width = 8, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.atlIfg, lr.sar.indIfg, lr.sar.pacIfg, order = c(6:7, 4, 3, 5, 2:1), leg.txt = c("All", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 15, srt = 50)
dev.off()

save(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, file = "Outputs/Atlantic_simplified.RData")
save(lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, file = "Outputs/Indian_simplified.RData")
save(lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, file = "Outputs/Pacific_simplified.RData")

rm(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, op.formula, ind.s, atl.s, pac.s)


## 4. Does resolution of variables matter? ---------------------------------

## 4i. Run full model for higher resolution ----------------------
mod.hres.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.4km, 3) + sdSST.4km + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2 + delta_carb_ion)^2, listw = op.s, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)

## 4iii. Compare LRs -------------------------------------------------------
lr.hres.op0 <- lr.calc(mod.hres.op0)

png("Figures/Ana_4iii_LRatio_ophresop0.png", width = 800)
# get order from running without order first
lr.plot(lr.sar.op0, lr.hres.op0, order = c(7, 13, 8, 14, 9, 15, 4:3, 12, 16, 11, 5, 1, 10, 6, 2), leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.hres.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.hres.op0g <- lr.calc(mod.hres.op0, tmp)
rm(tmp)

png("Figures/Ana_4iii_LRatio_g_ophresop0.png", width = 800)
lr.plot(lr.sar.op0g, lr.hres.op0g, order = c(6:7, 4, 3, 5, 2:1), leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

save(mod.hres.op0, lr.hres.op0, lr.hres.op0g, file = "Outputs/mod_hres.RData")
rm(mod.hres.op0, lr.hres.op0, lr.hres.op0g)


## 5. Does delta_carb_ion cut-off matter? -------------------------------------

## 5i. Create a dataset for modelling with higher or lower cut-offs --------
load("../../../Project/MARGO/Outputs/Environmental_variables.Rdata") # the datasets for the modelling
rm(db.traits, margo.traits)

# create a dataset for modelling with 
ldg.margo.mod.20 <- merge(ldg.margo.env, ldg.margo.data[, c(1, 3:5, 48, 54:58, 60, 65:82)], by.x = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"), by.y = c("Core", "Latitude", "Longitude", "Water.Depth", "Ocean2"))
rm(ldg.margo.data, ldg.margo.env)

ldg.margo.data <- ldg.mar

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0
# remove Water.Depth again (as it has NAs)
ldg.margo.mod <- ldg.margo.mod[, -which(names(ldg.margo.mod) == "Water.Depth")]
# set the FRic NAs to an arbitrary amount, as don't want to remove those rows
ldg.margo.mod$FRic[is.na(ldg.margo.mod$FRic)] <- 999
# remove other NAs
summary(ldg.margo.mod)
# most of the NAs are in rarefy.SR as it can't be calculated for data without total planktics.
ldg.margo.mod <- na.omit(ldg.margo.mod)
# reset these back to NA
ldg.margo.mod$FRic[ldg.margo.mod$FRic == 999] <- NA

# remove the extra factor level of the mediterranean
table(ldg.margo.mod$Ocean2)
ldg.margo.mod <- ldg.margo.mod[ldg.margo.mod$Ocean2 != "Mediterranean", ]
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)
table(ldg.margo.mod$Ocean2)

## which points are to be excluded because of delta_carb_ion?
with(ldg.margo.mod, plot(delta_carb_ion, rarefy.sr, pch = 16, col = Ocean2))
legend("topleft", levels(ldg.margo.mod$Ocean2), pch = 16, col = 1:3)
ldg.margo.mod <- ldg.margo.mod[which(ldg.margo.mod$delta_carb_ion >= -14), ]
with(ldg.margo.mod, distrib.map(Longitude, Latitude, Ocean2))
table(ldg.margo.mod$Ocean2)
# based on using a 3500m / 4500m cut-off

## 5ii. Run the model with these different cutoffs -------------------------
mod.hres.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.4km, 3) + sdSST.4km + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy  + Ocean2 + delta_carb_ion)^2, listw = op.s, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)

## 5iii. Compare LRs -------------------------------------------------------
lr.hres.op0 <- lr.calc(mod.hres.op0)

png("Figures/Ana_4iii_LRatio_ophresop0.png", width = 800)
# get order from running without order first
lr.plot(lr.sar.op0, lr.hres.op0, order = c(7, 13, 8, 14, 9, 15, 4:3, 12, 16, 11, 5, 1, 10, 6, 2), leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.hres.op0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.hres.op0g <- lr.calc(mod.hres.op0, tmp)
rm(tmp)

png("Figures/Ana_4iii_LRatio_g_ophresop0.png", width = 800)
lr.plot(lr.sar.op0g, lr.hres.op0g, order = c(6:7, 5, 3, 4, 2:1), leg.x = 17, leg.y = 275, leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

save(mod.hres.op0, lr.hres.op0, lr.hres.op0g, file = "Outputs/mod_hres.RData")
rm(mod.hres.op0, lr.hres.op0, lr.hres.op0g)


## 6. Does averaging explanatory variables matter? -------------------------

## 7. Does the number of points in the dataset make a difference? ----------------------


## 8. Evenness -------------------------------------------------------------

## 8i. Create an OLS model ------------------------------------------------
mod.eve.l0 <- lm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, data = ldg.margo.mod)

# check model plots
png("Figures/Ana_8i_modevel0.png", 600, 600)
par(mfrow = c(2, 2))
plot(mod.eve.l0)
par(mfrow = c(1, 1))
dev.off()

summary(mod.eve.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
mod.eve.l0.sac <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.eve.l0$residuals, latlon = TRUE, resamp = 1))
summary(mod.eve.l0.sac)
png("Figures/Ana_8i_modevel0SAC.png")
plot.spline.correlog.n(mod.eve.l0.sac, xlab = "Distance / km")
dev.off()

## 8ii. Run model optimisation ----------------------------------------------
mod.sar.eveW <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveW$obj, Nagelkerke = TRUE) # 0.47911
AIC(mod.sar.eveW$obj) # -3860.016

# check other coding styles
mod.sar.eveB <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveB$obj, Nagelkerke = TRUE) # 0.4892
AIC(mod.sar.eveB$obj) # -3896.698

mod.sar.eveS <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveS$obj, Nagelkerke = TRUE) # 0.50095
AIC(mod.sar.eveS$obj) # -3940.305

# So "S" is the best coding style and the best neighbourhood distance is 544.0801

eve.nb <- dnearneigh(ldg.coords, 0, mod.sar.eveS$dist, longlat = TRUE)
eve.s <- nb2listw(eve.nb, glist = NULL, style = "S", zero.policy = TRUE)
mod.sar.eve0 <- errorsarlm(mod.sar.eveS$mod, listw = eve.s, zero.policy = TRUE, tol.solve = 1e-18)
rm(eve.nb)

save(mod.sar.eveB, mod.sar.eveS, mod.sar.eveW, file = "Outputs/Evenness_coding.RData")
rm(mod.sar.eveB, mod.sar.eveW)

## 8ii. Run model simplification -------------------------------------------

## 8iii. Create a plot of model parameters --------------------------------------
summary(mod.sar.eve0, Nagelkerke = T) # r2 = 0.45442

# generate a dataframe of coefficients
m.eve.coef <- data.frame(names = names(mod.sar.eve0$coefficients), coef.sar = mod.sar.eve0$coefficients, row.names = 1:length(mod.sar.eve0$coefficients), stars = NA)

# add a column of significance stars
for (i in 1:length(stars)) {
  m.eve.coef$stars[which(summary(mod.sar.eve0)$Coef[, 4] <= stars[i] & is.na(m.eve.coef$stars))] <- names(stars)[i]
}

# plot the absolute coefficients
png("Figures/Ana_8iii_coef_modsareve0.png", width = 1000, height = 750)
plt.def <- par("plt")
par(plt = c(plt.def[1:2], 0.5, plt.def[4]))
tmp.x <- barplot(abs(m.eve.coef$coef.sar), names = m.eve.coef$names, las = 2, ylim = c(0, max(m.eve.coef$coef.sar) + 50))
text(tmp.x, abs(m.eve.coef$coef.sar) + 10, m.eve.coef$stars)
par(plt = plt.def)
dev.off()
rm(plt.def, tmp.x, m.eve.coef)

## 8iv. Calculate likelihood ratios for the SAR model ----------------------
# full model is
summary(mod.sar.eve0, Nagelkerke = T) # r2 = 0.50095 
AIC(mod.sar.eve0) # -3940.305

lr.sar.eve0 <- lr.calc(mod.sar.eve0)

png("Figures/Ana_8iv_LRatio_eve0.png", width = 800)
lr.plot(lr.sar.eve0, order = c(7:9, 3:2, 12:11, 4:5, 10, 6, 1), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

## 8v. Calculate likelihood ratios for groups of EVs ---------------------

# for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.eve0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.eve0g <- lr.calc(mod.sar.eve0, tmp)
rm(tmp)

png("Figures/Ana_8v_LRatio_g_eve0.png", width = 800)
lr.plot(lr.sar.eve0g, order = c(6:7, 5, 3, 4, 1:2), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

# mnt, d10, prod (can't tell if seasonality or productivity), seasonality (sdt, sdsal), salinity, delta_carb_ion, ocean

# full model lr plots for each ocean? map differences of RV between simplified model with and without ocean (by points not by layer) - Supplementary info

## 8vi. Effect of delta_carb_ion ----------------------------------------------


## 9. Lineage age / FRic -------------------------------------------------------------

## 9i. Lineage age models -------------------------------------------------
mod.sar.lnaW <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaW$obj, Nagelkerke = TRUE) # 0.63489
AIC(mod.sar.lnaW$obj) # 5565.102

# check other coding styles
mod.sar.lnaB <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaB$obj, Nagelkerke = TRUE) # 0.61397 
AIC(mod.sar.lnaB$obj) #  5669.584

mod.sar.lnaS <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaS$obj, Nagelkerke = TRUE) # 0.63364 
AIC(mod.sar.lnaS$obj) # 5571.538

mod.sar.lnaW
# So "W" is the best coding style and the best neighbourhood distance is 537.6619

lna.nb <- dnearneigh(ldg.coords, 0, mod.sar.lnaW$dist, longlat = TRUE)
lna.w <- nb2listw(lna.nb, glist = NULL, style = "W", zero.policy = TRUE)
mod.sar.lna0 <- errorsarlm(mod.sar.lnaW$mod, listw = lna.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mod.sar.lna0, Nagelkerke = TRUE) # 0.63489

save(mod.sar.lnaB, mod.sar.lnaS, mod.sar.lnaW, file = "Outputs/Lineage_coding.RData")
rm(lna.nb, mod.sar.lnaB, mod.sar.lnaS)

lr.sar.lna0 <- lr.calc(mod.sar.lna0)

png("Figures/Ana_9i_LRatio_lna0.png", width = 800)
lr.plot(lr.sar.lna0, ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.lna0), group = c("Stability", "Productivity", "Stress", "Stability", "Stress", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature", "Vertical niche structure", "Vertical niche structure")))
lr.sar.lna0g <- lr.calc(mod.sar.lna0, tmp)
rm(tmp)

png("Figures/Ana_9i_LRatio_g_lna0.png", width = 800)
lr.plot(lr.sar.lna0g, order = c(6:7, 5, 3, 4, 1:2), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

## 9ii. Dissolution cutoffs ------------------------------------------------

## 9iii. Functional richness -----------------------------------------------
mod.sar.fricW <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricW$obj, Nagelkerke = TRUE) # 0.87776
AIC(mod.sar.fricW$obj) # -2791.866

# check other coding styles
mod.sar.fricB <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricB$obj, Nagelkerke = TRUE) # 0.87663 
AIC(mod.sar.fricB$obj) #  -2774.851

mod.sar.fricS <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, FRic ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.fricS$obj, Nagelkerke = TRUE) # 0.88194 
AIC(mod.sar.fricS$obj) # -2856.622

mod.sar.fricS
# So "S" is the best coding style and the best neighbourhood distance is 537.6619

fric.nb <- dnearneigh(ldg.coords, 0, mod.sar.fricS$dist, longlat = TRUE)
fric.s <- nb2listw(fric.nb, glist = NULL, style = "S", zero.policy = TRUE)
mod.sar.fric0 <- errorsarlm(mod.sar.fricS$mod, listw = fric.s, zero.policy = TRUE, tol.solve = 1e-18)
summary(mod.sar.fric0, Nagelkerke = TRUE) # 0.58086
rm(fric.nb, fric.s)

save(mod.sar.fricB, mod.sar.fricS, mod.sar.fricW, file = "Outputs/Lineage_coding.RData")
rm(fric.nb, mod.sar.fricB, mod.sar.fricS)

lr.sar.fric0 <- lr.calc(mod.sar.fric0)

png("Figures/Ana_9iii_LRatio_fric0.png", width = 800)
lr.plot(lr.sar.fric0, order = c(7:9, 3:2, 12:11, 4:5, 10, 6, 1), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.fric0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.fric0g <- lr.calc(mod.sar.fric0, tmp)
rm(tmp)

png("Figures/Ana_9iii_LRatio_g_fric0.png", width = 800)
lr.plot(lr.sar.fric0g, order = c(6:7, 5, 3, 4, 1:2), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()



## 10. Paper outputs -------------------------------------------------------

## 10i. Table of full model cofficients for the models -----------
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

## 10ii. likelihood ratios for the models ----------------------------------------
(lr.out <- cbind(lr.sar.op0, lr.sar.eve0[2:4], lr.sar.lna0[2:4]))
names(lr.out)[2:4] <- paste("sr", names(lr.out)[2:4], sep = "_")
names(lr.out)[5:7] <- paste("eve", names(lr.out)[5:7], sep = "_")
names(lr.out)[8:10] <- paste("lna", names(lr.out)[8:10], sep = "_")

write.csv(lr.out, "Outputs/lr_out.csv")
rm(lr.out)

# (lr.oce.out <- merge(lr.sar.op0, lr.sar.atlf, by = "names", all = T))
# (lr.oce.out <- merge(lr.oce.out, lr.sar.indf, by = "names", all = T))
# (lr.oce.out <- merge(lr.oce.out, lr.sar.pacf, by = "names", all = T))
# names(lr.oce.out)[2:4] <- paste("All", names(lr.sar.op0)[2:4], sep = "_")
# names(lr.oce.out)[5:7] <- paste("atl", names(lr.sar.op0)[2:4], sep = "_")
# names(lr.oce.out)[8:10] <- paste("ind", names(lr.sar.op0)[2:4], sep = "_")
# names(lr.oce.out)[11:13] <- paste("pac", names(lr.sar.op0)[2:4], sep = "_")

# write.csv(lr.oce.out, "Outputs/lr_oce_out.csv")

## 10iii. log likelihood ratio plot for full vs simplified model rarified SR --------

png("Figures/Ana_10iii_LRatio_op0f.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(5, 7, 6, 3:4, 2:1), ylab = "Log Likelihood ratio", leg.txt = c("Full", "Simplified"), cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

## 10iv. log likelihood ratio plot comparing different delta_carb_ion cut-of --------
# (tmp <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
# lr.sar.op0g <- lr.calc(mod.sar.op0, tmp)
# rm(tmp)
# 
# (tmp <- data.frame(names = model.evs(mod.sar8.op0),group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
# lr.sar8.op0g <- lr.calc(mod.sar8.op0, tmp)
# rm(tmp)
# 
# (tmp <- data.frame(names = model.evs(mod.sar4.op0),group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
# lr.sar4.op0g <- lr.calc(mod.sar4.op0, tmp)
# rm(tmp)
# 
# png("Figures/Ana_10iv_LRatio_g_op468.png", width = 800)
# lr.plot(lr.sar8.op0g, lr.sar.op0g, lr.sar4.op0g, order = c(6, 7, 5, 3, 4, 2, 1), leg.txt = c("8%", "6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2)
# dev.off()
# 
# png("Figures/Ana_10iv_LRatio_g_op46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.op0g, lr.sar4.op0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, srt = 45, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2)
# dev.off()
# 
# 
# png("Figures/Ana_10iv_LRatio_g_eve46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.eve0g, lr.sar4.eve0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 10, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
# dev.off()
# 
# png("Figures/Ana_10iv_LRatio_g_lna46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.lna0g, lr.sar4.lna0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
# dev.off()

## 10v. residuals map for full model --------------------------------------------
png("Figures/Ana_10v_rsp_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.op0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -8, max.col = 8, maintitle = "Residuals for rarefied richness"))
dev.off()

png("Figures/Ana_10v_eve_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.eve0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -0.5, max.col = 0.5, maintitle = "Residuals for evenness"))
dev.off()

png("Figures/Ana_10v_lna_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.lna0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -6, max.col = 6, maintitle = "Residuals for average community age"))
dev.off()

## 10vi. coplot ------------------------------------------------------------------

# plot marginal effects (hold all else constant) of coefficients of T polynomial in oceans. Plot across temperature range observed in the oceans (3 colours / plotting symbols)

## 10vii. comparing all the models ------------------------------------------------
png("Figures/Ana_10vii_LRatio_REM.png", width = 800)
lr.plot(lr.sar.op0, lr.sar.eve0, lr.sar.lna0, order = c(7:9, 4, 1, 11:10, 3, 5, 2, 6), ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age"))
dev.off()

png("Figures/Ana_10vii_LRatio_g_REM.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.eve0g, lr.sar.lna0g, order = c(6:7, 5, 3:4, 2:1),  ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age"), srt = 45)
dev.off()

## 10viii. compare this to a random model ----------------------------------------
# rdm.vals <- data.frame(model = rep(NA, 1000 * 7), names = NA, lr = NA, p = NA, stars = NA)
# tmp <- data.frame(names = model.evs(mod.rdm), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature"))
# 
# for(i in 1:1000) {
#   print(i)
#   mod.rdm <- errorsarlm(rnorm(nrow(ldg.margo.mod)) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + I(mean.mld.t/10) + I(depth10deg/100) + logProd.mn.ann + absMnSal.0m + sdSal.0m + prop2.oxy + Ocean2 + delta_carb_ion)^2, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)
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


## 11. Predicting models ---------------------------------------------------

## 11i. Set up the dataset for prediction ----------------------------------
with(ldg.margo.mod, distrib.map(Longitude, Latitude, predict(mod.sar.op0)))

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
  tmp.col <- which(names(ldg.margo.mod) == names(ldg.p.margo)[i])
  ldg.p.margo[which(ldg.p.margo[, i] < min(ldg.margo.mod[, tmp.col], na.rm = TRUE) | ldg.p.margo[, i] > max(ldg.margo.mod[, tmp.col], na.rm = TRUE)), i] <- NA
  if (sum(is.na(ldg.p.margo[, i])) > 0)
    with(ldg.p.margo[is.na(ldg.p.margo[, i]), ], distrib.map(Longitude, Latitude, Ocean2, pch = 15, cex = 0.5, main = names(ldg.p.margo)[i]))
}
par(ask = FALSE)
rm(i, tmp.col)
ldg.p.margo <- na.omit(ldg.p.margo)

# plot these up
for (i in env.var) {
  png(paste("Figures/Ana_11i_map_", i, ".png", sep = ""), width = 800, height = 450)
  with(ldg.p.margo, distrib.map(Longitude, Latitude, ldg.p.margo[, i], palette = "matlab.like", pch = 15, cex = 0.5, main = i, col.land = "black", col.water = "white"))
  dev.off()
}
rm(i)

## 11ii. predict species richness for this dataset -------------------------------
# necessary to use sar.predict not predict as the ordinary predict.sarlm function cannot cope with poly variables
ldg.p.margo$rarefy.sr <- sar.predict(mod.sar.op0, newdata = ldg.p.margo, olddata = ldg.margo.mod)
summary(ldg.p.margo$rarefy.sr)

with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with negative richness - where are these
with(ldg.p.margo[ldg.p.margo$rarefy.sr <= 0, ], distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, cex = 0.4)) # around coastlines

# compare with observed
png("Figures/Ana_11ii_rsr_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$rarefy.sr > 0, ], distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", min.col = 0, max.col = 27))
dev.off()

png("Figures/Ana_11ii_rsr_obs.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", min.col = 0, max.col = 27, col.water = "white", col.land = "black"))
dev.off()

## 11iii. predict evenness for this dataset ---------------------------------------
ldg.p.margo$simpsonEve <- sar.predict(mod.sar.eve0, newdata = ldg.p.margo, olddata = ldg.margo.mod)
summary(ldg.p.margo$simpsonEve)

with(ldg.p.margo, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with evenness outside sensible limits - where are these
with(ldg.p.margo[ldg.p.margo$simpsonEve > 1, ], distrib.map(Longitude, Latitude, simpsonEve, pch = 15, cex = 0.4)) # around arctic / antarctica coastlines - not a problem

# compare with observed
png("Figures/Ana_11iii_eve_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$simpsonEve <= 1, ], distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_11iii_eve_obs.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 11iv. predict average community age for this dataset --------------------
ldg.p.margo$MorphoAgeAbun <- sar.predict(mod.sar.lna0, newdata = ldg.p.margo, olddata = ldg.margo.mod)
summary(ldg.p.margo$MorphoAgeAbun)

with(ldg.p.margo, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with evenness outside sensible limits - where are these
plot(sort(ldg.p.margo$MorphoAgeAbun))
summary(ldg.margo.mod$MorphoAgeAbun)
with(ldg.p.margo[ldg.p.margo$MorphoAgeAbun < 5, ], distrib.map(Longitude, Latitude, MorphoAgeAbun, pch = 15, cex = 0.4)) # around coastlines - very few points

# compare with observed
png("Figures/Ana_11iv_lna_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$MorphoAgeAbun >= 5 & ldg.p.margo$MorphoAgeAbun <= 17, ], distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_11iv_lna_obs.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 11v. Difference maps ----------------------------------------------------
# set up centred and scaled
ldg.margo.mod$rarefy.sr.cs <- (ldg.margo.mod$rarefy.sr - mean(ldg.margo.mod$rarefy.sr)) / sd(ldg.margo.mod$rarefy.sr)
ldg.margo.mod$simpsonEve.cs <- (ldg.margo.mod$simpsonEve - mean(ldg.margo.mod$simpsonEve)) / sd(ldg.margo.mod$simpsonEve)
ldg.margo.mod$MorphoAgeAbun.cs <- (ldg.margo.mod$MorphoAgeAbun - mean(ldg.margo.mod$MorphoAgeAbun)) / sd(ldg.margo.mod$MorphoAgeAbun)

ldg.p.margo$rarefy.sr.cs <- (ldg.p.margo$rarefy.sr - mean(ldg.p.margo$rarefy.sr)) / sd(ldg.p.margo$rarefy.sr)
ldg.p.margo$simpsonEve.cs <- (ldg.p.margo$simpsonEve - mean(ldg.p.margo$simpsonEve)) / sd(ldg.p.margo$simpsonEve)
ldg.p.margo$MorphoAgeAbun.cs <- (ldg.p.margo$MorphoAgeAbun - mean(ldg.p.margo$MorphoAgeAbun)) / sd(ldg.p.margo$MorphoAgeAbun)

## differences
# b/w sr and eve
# red - high sr & low eve
# blue - low sr & high eve
png("Figures/Ana_11v_sr_eve.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr.cs - simpsonEve.cs, palette = "rwb", col.water = "white", col.land = "black", min.col = -10, max.col = 10))
dev.off()

# b/w sr and lna
# red - high sr & low lna
# blue - low sr & high lna
png("Figures/Ana_11v_sr_lna.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr.cs - MorphoAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black"))
dev.off()

# b/w eve and lna
# red - high eve & low lna
# blue - low eve & high lna
png("Figures/Ana_11v_eve_lna.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, simpsonEve.cs - MorphoAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", min.col = -8, max.col = 8))
dev.off()

# predicted
# b/w sr and eve
# red - high sr & low eve
# blue - low sr & high eve
png("Figures/Ana_11v_sr_eve_pred.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr.cs - simpsonEve.cs, palette = "rwb", col.water = "white", col.land = "black", min.col = -15, max.col = 15, pch = 15, cex = 0.4))
dev.off()

# b/w sr and lna
# red - high sr & low lna
# blue - low sr & high lna
png("Figures/Ana_11v_sr_lna_pred.png", 700, 500)
with(ldg.p.margo, distrib.map(Longitude, Latitude, rarefy.sr.cs - MorphoAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4, min.col = -15, max.col = 15))
dev.off()

# b/w eve and lna
# red - high eve & low lna
# blue - low eve & high lna
png("Figures/Ana_11v_eve_lna_pred.png", 700, 500)
with(ldg.p.margo[abs(ldg.p.margo$simpsonEve.cs - ldg.p.margo$MorphoAgeAbun.cs) <= 15, ], distrib.map(Longitude, Latitude, simpsonEve.cs - MorphoAgeAbun.cs, palette = "rwb", col.water = "white", col.land = "black", pch = 15, cex = 0.4))
dev.off()

## 11vi. Comparison plots --------------------------------------------------
load("Outputs/Atlantic_simplified.RData")
load("Outputs/Indian_simplified.RData")
load("Outputs/Pacific_simplified.RData")
load("Outputs/mod_hres.RData")

# for full rarefied model
lr.calc(mod.sar.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "op0")
(tmp <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.op0g <- lr.calc(mod.sar.op0, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "op0g")
rm(tmp)

# for simplified rarefied model
lr.calc(mod.sar.opf, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "opf")
(tmp <- data.frame(names = model.evs(mod.sar.opf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.opfg <- lr.calc(mod.sar.opf, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "opfg")
rm(tmp)

# for Atlantic
lr.sar.atlIf <- lr.calc(mod.sar.atlIf, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "atlIf")
# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.atlIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.atlIfg <- lr.calc(mod.sar.atlIf, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "atlIfg")
rm(tmp)

# for Indian
lr.sar.indIf <- lr.calc(mod.sar.indIf, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "indIf")
(tmp <- data.frame(names = model.evs(mod.sar.indIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.indIfg <- lr.calc(mod.sar.indIf, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "indIfg")
rm(tmp)

# for Pacific
lr.sar.pacIf <- lr.calc(mod.sar.pacIf, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "pacIf")
(tmp <- data.frame(names = model.evs(mod.sar.pacIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.pacIfg <- lr.calc(mod.sar.pacIf, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "pacIfg")
rm(tmp)

# for lineage age
lr.calc(mod.sar.lna0, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "lna")
(tmp <- data.frame(names = model.evs(mod.sar.lna0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.lna0g <- lr.calc(mod.sar.lna0, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "lnag")
rm(tmp)


# for evenness
lr.calc(mod.sar.eve0, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "eve")
(tmp <- data.frame(names = model.evs(mod.sar.eve0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.eve0g <- lr.calc(mod.sar.eve0, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "eveg")
rm(tmp)

# for high resolution
lr.hres.op0 <- lr.calc(mod.hres.op0, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "hres")
(tmp <- data.frame(names = model.evs(mod.hres.op0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.hres.op0g <- lr.calc(mod.hres.op0, tmp, plots = TRUE, pred.data = ldg.p.margo, mod.data = ldg.margo.mod, file.nm = "hresg")
rm(tmp)


# 12. Metabolic theory of ecology -----------------------------------------
# log transformed SR
Ln_SR <- log(ldg.margo.mod$rarefy.sr)

# 1 / kT (boltzman constant in eV K-1) * absolute T
MTE_SST <- 1 / (8.6173324E−5 * (ldg.margo.mod$meanSST.1deg + 273.15))

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
MTE_oce <- ldg.margo.mod$Ocean2
mte.oce.mod <- errorsarlm(Ln_SR ~ MTE_SST*MTE_oce, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mte.oce.mod)
anova(mte.oce.mod, mte.mod) # significant difference, implying that dropping ocean produces a significantly worse model
p.ocean <- data.frame(MTE_SST = rep(p.MTE_SST, 3), MTE_oce = rep(levels(MTE_oce), each = 100))
p.ocean$p.Ln_SR <- predict(mte.oce.mod, p.ocean)[, 1]

# plot them
png("Figures/Ana_12_mte_plot.png")
plot(Ln_SR ~ MTE_SST, pch = 16, col = ldg.margo.mod$Ocean2, xlab = "Temperature (1 / kT)", ylab = "ln (rarefied species richness)", bty = "l", las = 1, cex.lab = 1.2, cex.axis = 1.2)
points(p.MTE_SST, p.Ln_SR, type = "l", lwd = 3) # observed
with(p.ocean[p.ocean$MTE_oce == "Atlantic",], points(MTE_SST, p.Ln_SR, type = "l", col = 1)) # Atlantic
with(p.ocean[p.ocean$MTE_oce == "Indian",], points(MTE_SST, p.Ln_SR, type = "l", col = 2)) # Indian
with(p.ocean[p.ocean$MTE_oce == "Pacific",], points(MTE_SST, p.Ln_SR, type = "l", col = 3)) # Pacific
points(p.MTE_SST, mte.Ln_SR, type = "l", lty = 2, lwd = 3) # predicted
legend("topright", levels(ldg.margo.mod$Ocean2), pch = 16, col = 1:3)
dev.off()

# plot this relationship on the SST ~ SR relationship
png("Figures/Ana_12_sst_mte_plot.png")
with(ldg.margo.mod, plot(rarefy.sr ~ meanSST.1deg, pch = 16, col = ldg.margo.mod$Ocean2, xlab = expression(paste("SST / ", degree, "C")), ylab = "Rarefied species richness", bty = "l", las = 1, cex.lab = 1.2, cex.axis = 1.2))
points((1 / (8.6173324E−5 * p.MTE_SST) - 273.15), exp(mte.Ln_SR), type = "l", lwd = 2)
legend(0, 24, levels(ldg.margo.mod$Ocean2)[1:3], pch = 16, col = 1:3)
dev.off()

save(Ln_SR, MTE_SST, mte.mod, mte.oce.mod, p.Ln_SR, p.MTE_SST, file = "Outputs/Metabolic_hypothesis.RData")
rm(Ln_SR, MTE_SST, MTE_oce, mn.Ln_SR, mn.MTE_SST, mte.Ln_SR, mte.mod, mte.oce.mod, p.Ln_SR, p.MTE_SST, p.ocean)


# 13. Tidy up -------------------------------------------------------------
save(ldg.margo.mod, file = "Outputs/ldg_margo_mod.RData")
save(lr.sar.op0, lr.sar.op0g, mod.l0.sac, mod.sar.op0, mod.sar.opW, file = "Outputs/Richness_model.RData")
save(lr.sar.opf, lr.sar.opfg, ms.lr, ms.lr.group, mod.sar.opf, file = "Outputs/Richness_model_simplified.RData")
save(lr.sar.eve0, lr.sar.eve0g, mod.eve.l0, mod.eve.l0.sac, mod.sar.eve0, file = "Outputs/Evenness_model.RData")
save(lr.sar.lna0, lr.sar.lna0g, mod.sar.lna0, file = "Outputs/Lineage_model.RData")

rm(ldg.coords, stars, env.var, op.w, eve.s, lna.w, mod.sar.eveS, mod.sar.lnaW)

