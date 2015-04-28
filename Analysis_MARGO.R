## Created: 15 / 4 / 2015
## Last edited: 27 / 4 / 2015
## Isabel Fenton
## Reanalysis for LDG paper

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


## 0. Setting up the dataset -----------------------------------------------
rm(db.traits, ldg.margo, margo.traits, WOA.depths)

# create a dataset for modelling with 
ldg.margo.mod <- ldg.margo.env[, c(1:3, 5:26, 31:33)]
ldg.margo.mod <- cbind(ldg.margo.mod, ldg.margo.data[, c(48, 55:58, 65:82)])
rm(ldg.margo.data, ldg.margo.env)

# set NAs in 10 deg depth to 0, so it can be modelled
ldg.margo.mod$depth10deg[is.na(ldg.margo.mod$depth10deg)] <- 0

# remove the extra factor level of the mediterranean
ldg.margo.mod$Ocean2 <- droplevels(ldg.margo.mod$Ocean2)

# remove other NAs
summary(ldg.margo.mod)
ldg.margo.mod <- na.exclude(ldg.margo.mod)
  
# # which points are to be excluded because of carb_ion?
# ldg.margo.mod <- ldg.margo.mod[-which(ldg.margo.mod$carb_ion > 6), ] 
# currently not doing this

dim(ldg.margo.mod)

# create plots
png("Figures/map_rsr.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", main = "Rarefied species richness", col.water = "white", col.land = "black"))
dev.off()

png("Figures/map_eve.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", main = "Simpson's Evenness", col.water = "white", col.land = "black", min.col = 0.1))
dev.off()

png("Figures/map_lna.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", main = "Average Community Age", col.water = "white", col.land = "black", max.col = 16.3))
dev.off()


## 1. Check for correlation between explanatory variables ----------------------
names(ldg.margo.mod)

# variables are: mean/sd SST, mean/sd MLD, 10deg contour, mean/sd logProd, mean/sd Sal, Ocean2, carbonate ion 
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "sd.mld.t", "depth10deg", "logProd.mn.ann", "logProd.sd.ann", "meanSal.0m", "sdSal.0m", "Ocean2", "carb_ion")

# pairs plot
png("Figures/Ana_1_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var])

# currently the only ones that are potentially correlated (VIF > 5 and look correlated on the pairs plot) are mean and sd in Prod, and mean / sd mld. Look into this in a bit more detail
png("Figures/Ana_1_mnprodsdprod1deg.png")
with(ldg.margo.mod, plot(logProd.mn.ann, logProd.sd.ann, pch = 16)) # highly correlated.
dev.off()

# does the same hold for for the mld
png("Figures/Ana_1_mnMLDsdMLD.png")
with(ldg.margo.mod, plot(mean.mld.t, sd.mld.t, pch = 16)) # yes
dev.off()

# therefore suggest exclusion of logProd.sd.ann and sd.mld.t from models, so
env.var <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "meanSal.0m", "sdSal.0m", "Ocean2", "carb_ion")

# check this has fixed the problem
pairs(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var], )
vif(ldg.margo.mod[, names(ldg.margo.mod) %in% env.var], )
# it has, so now have new EVs


## 2. Model simplification -------------------------------------------------

## 2i. Create an OLS model and check for SAC --------------------------------
# an OLS model
mod.l0 <- lm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, data = ldg.margo.mod)

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
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.l0$residuals - min(mod.l0$residuals), key = FALSE))
dev.off()

rm(mod.l0)

## 2ii. Create a GAM to check complexity / SAC -----------------------------
# n.b. GAMs can't do interactions (as additive models)
mod.g0 <- with(ldg.margo.mod, gam(rarefy.sr ~ s(Longitude, Latitude, k = 80, by = Ocean2) + s(meanSST.1deg, by = Ocean2) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion, gamma = 1.4))
summary(mod.g0)

png("Figures/Ana_2ii_modg0.png")
gam.check(mod.g0) # n.b. this gives an error for factors
dev.off()

# ## calculate SAC
# # haven't run this
# # using spline.correlog
# mod.g0.SAC <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.g0$residuals, latlon = TRUE))
# summary(mod.g0.SAC)
# png("Figures/Ana_2ii_modg0SAC.png")
# plot.spline.correlog.n(mod.g0.SAC, xlab = "Distance / km")
# dev.off()

# using correlog
mod.g0.SACcor <- with(ldg.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.g0), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2ii_modg0SACcor.png")
plot(mod.g0.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()

rm(mod.g0, mod.g0.SACcor)

## 2iii. Create an optimised SARerror model --------------------------------
# Make a matrix of coordinates (X and Y coordinates)
ldg.coords <- cbind(ldg.margo.mod$Longitude,ldg.margo.mod$Latitude)
ldg.coords <- as.matrix(ldg.coords)

# run model optimisation
mod.sar.opW <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.opW$obj)

## check SAC has been removed
# using spline.correlog
# haven't run this
# mod.sar.opW.SAC <- with(ldg.margo.mod, spline.correlog(Longitude, Latitude, mod.sar.opW$obj$residuals, latlon = TRUE))
# summary(mod.sar.opW.SAC)
# png("Figures/Ana_2iii_modSarOp0SAC.png")
# plot.spline.correlog.n(mod.sar.opW.SAC, xlab = "Distance / km")
# dev.off()

# using correlog
mod.sar.opW.SACcor <- with(ldg.margo.mod, correlog(Longitude, Latitude, z = residuals(mod.sar.opW$obj), na.rm = T, increment = 100, resamp = 1, latlon = T))
png("Figures/Ana_2iii_modSarOp0SACcor.png")
plot(mod.sar.opW.SACcor$correlation, type = "b", pch = 1, cex = 1.2, lwd = 1.5, ylim = c(-0.5, 1), xlab = "distance", ylab = "Moran's I", cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0)
dev.off()
rm(mod.sar.opW.SACcor)

# check whether different coding methods improve the AIC
mod.sar.opB <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 9652.784
AIC(mod.sar.opB$obj) # 9760.968
# the AIC is the first number printed in best.dist

mod.sar.opS <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
AIC(mod.sar.opW$obj) # 9652.784
AIC(mod.sar.opS$obj) # 9676.573

# So "W" is the best coding style and the best neighbourhood distance is 619.3581
rm(mod.sar.opS, mod.sar.opB)

## 2iii.i. Collinearity for optimised models -------------------------------
#vif.sar(errorsarlm(, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18))

## 2iv. Run model simplification -------------------------------------------
summary(mod.sar.opW$obj)

# re-run this optimised model through errorsarlm, so things like anova and lr.calc work
op.nb <- dnearneigh(ldg.coords, 0, mod.sar.opW$dist, longlat = TRUE)
op.w <- nb2listw(op.nb, glist = NULL, style = "W", zero.policy = TRUE)
mod.sar.op0 <- errorsarlm(mod.sar.opW$mod, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18)
rm(op.nb)

# start model simplification
summary(mod.sar.op0, Nagelkerke = TRUE) # 
mod.sar.op1 <- update(mod.sar.op0, ~. -depth10deg:carb_ion)
anova(mod.sar.op0, mod.sar.op1)
summary(mod.sar.op1, Nagelkerke = TRUE) # 0.81801 
AIC(mod.sar.op1) # 9650.784

mod.sar.op2 <- update(mod.sar.op1, ~. -sdSST.1deg:sdSal.0m)
anova(mod.sar.op2, mod.sar.op1)
summary(mod.sar.op2, Nagelkerke = TRUE) # 0.818 
AIC(mod.sar.op2) # 9648.786
rm(mod.sar.op1)

mod.sar.op3 <- update(mod.sar.op2, ~. -poly(meanSST.1deg, 3):sdSST.1deg + poly(meanSST.1deg, 2):sdSST.1deg)
anova(mod.sar.op3, mod.sar.op2)
summary(mod.sar.op3, Nagelkerke = TRUE) # 0.818 
AIC(mod.sar.op3) # 9646.9
rm(mod.sar.op2)

mod.sar.op4 <- update(mod.sar.op3, ~. -sdSST.1deg:carb_ion)
anova(mod.sar.op4, mod.sar.op3)
summary(mod.sar.op4, Nagelkerke = TRUE) # 0.81798
AIC(mod.sar.op4) # 9645.047
rm(mod.sar.op3)

mod.sar.op5 <- update(mod.sar.op4, ~. -mean.mld.t:meanSal.0m)
anova(mod.sar.op5, mod.sar.op4)
summary(mod.sar.op5, Nagelkerke = TRUE) #  0.81797
AIC(mod.sar.op5) # 9643.232
rm(mod.sar.op4)

mod.sar.op6 <- update(mod.sar.op5, ~. -mean.mld.t:carb_ion)
anova(mod.sar.op6, mod.sar.op5)
summary(mod.sar.op6, Nagelkerke = TRUE) # 0.81794
AIC(mod.sar.op6) # 9641.616
rm(mod.sar.op5)

mod.sar.op7 <- update(mod.sar.op6, ~. -poly(meanSST.1deg, 3):carb_ion + poly(meanSST.1deg, 2):carb_ion)
anova(mod.sar.op7, mod.sar.op6)
summary(mod.sar.op7, Nagelkerke = TRUE) # 0.81791
AIC(mod.sar.op7) # 9639.988
rm(mod.sar.op6)

mod.sar.op8 <- update(mod.sar.op7, ~. -meanSal.0m:sdSal.0m)
anova(mod.sar.op8, mod.sar.op7)
summary(mod.sar.op8, Nagelkerke = TRUE) # 0.81786
AIC(mod.sar.op8) # 9638.584
rm(mod.sar.op7)

mod.sar.op9 <- update(mod.sar.op8, ~. -mean.mld.t:sdSal.0m)
anova(mod.sar.op9, mod.sar.op8)
summary(mod.sar.op9, Nagelkerke = TRUE) # 0.81783
AIC(mod.sar.op9) # 9636.994
rm(mod.sar.op8)

mod.sar.op10 <- update(mod.sar.op9, ~. -poly(meanSST.1deg, 3):sdSal.0m + poly(meanSST.1deg, 2):sdSal.0m)
anova(mod.sar.op10, mod.sar.op9)
summary(mod.sar.op10, Nagelkerke = TRUE) # 0.81777
AIC(mod.sar.op10) # 9635.708
rm(mod.sar.op9)

mod.sar.op11 <- update(mod.sar.op10, ~. -sdSal.0m:poly(meanSST.1deg, 2) + sdSal.0m:poly(meanSST.1deg, 1))
anova(mod.sar.op11, mod.sar.op10)
summary(mod.sar.op11, Nagelkerke = TRUE) # 0.81777
AIC(mod.sar.op11) # 9633.732
rm(mod.sar.op10)

mod.sar.op12 <- update(mod.sar.op11, ~. -sdSal.0m:Ocean2)
anova(mod.sar.op12, mod.sar.op11)
summary(mod.sar.op12, Nagelkerke = TRUE) # 0.81773
AIC(mod.sar.op12) # 9630.254
rm(mod.sar.op11)

mod.sar.op13 <- update(mod.sar.op12, ~. -depth10deg:logProd.mn.ann )
anova(mod.sar.op13, mod.sar.op12)
summary(mod.sar.op13, Nagelkerke = TRUE) # 0.81766
AIC(mod.sar.op13) # 9629.058
rm(mod.sar.op12)

mod.sar.op14 <- update(mod.sar.op13, ~. -mean.mld.t:depth10deg)
anova(mod.sar.op14, mod.sar.op13)
summary(mod.sar.op14, Nagelkerke = TRUE) # 0.81761
AIC(mod.sar.op14) # 9627.75
rm(mod.sar.op13)

mod.sar.op15 <- update(mod.sar.op14, ~. -poly(meanSST.1deg, 3):depth10deg + poly(meanSST.1deg, 2):depth10deg)
anova(mod.sar.op15, mod.sar.op14)
summary(mod.sar.op15, Nagelkerke = TRUE) # 0.81755
AIC(mod.sar.op15) # 9626.433
rm(mod.sar.op14)

mod.sar.op16 <- update(mod.sar.op15, ~. -poly(meanSST.1deg, 3):mean.mld.t + poly(meanSST.1deg, 2):mean.mld.t)
anova(mod.sar.op16, mod.sar.op15)
summary(mod.sar.op16, Nagelkerke = TRUE) # 0.8175
AIC(mod.sar.op16) # 9625.098
rm(mod.sar.op15)

mod.sar.op17 <- update(mod.sar.op16, ~. -poly(meanSST.1deg, 3):logProd.mn.ann + poly(meanSST.1deg, 2):logProd.mn.ann)
anova(mod.sar.op17, mod.sar.op16)
summary(mod.sar.op17, Nagelkerke = TRUE) # 0.81743
AIC(mod.sar.op17) # 9623.913
rm(mod.sar.op16)

mod.sar.op18 <- update(mod.sar.op17, ~. -sdSST.1deg:logProd.mn.ann)
anova(mod.sar.op18, mod.sar.op17)
summary(mod.sar.op18, Nagelkerke = TRUE) # 0.81739
AIC(mod.sar.op18) # 9622.435
rm(mod.sar.op17)

mod.sar.op19 <- update(mod.sar.op18, ~. -sdSST.1deg:poly(meanSST.1deg, 2) + sdSST.1deg:poly(meanSST.1deg, 1))
anova(mod.sar.op19, mod.sar.op18)
summary(mod.sar.op19, Nagelkerke = TRUE) # 0.81733
AIC(mod.sar.op19) # 9621.179
rm(mod.sar.op18)

mod.sar.op20 <- update(mod.sar.op19, ~. -sdSST.1deg:poly(meanSST.1deg, 1))
anova(mod.sar.op20, mod.sar.op19)
summary(mod.sar.op20, Nagelkerke = TRUE) # 0.81725
AIC(mod.sar.op20) # 9620.159
rm(mod.sar.op19)

mod.sar.op21 <- update(mod.sar.op20, ~. -depth10deg:Ocean2)
anova(mod.sar.op21, mod.sar.op20)
summary(mod.sar.op21, Nagelkerke = TRUE) # 0.81711
AIC(mod.sar.op21) # 9617.83
rm(mod.sar.op20)

mod.sar.op22 <- update(mod.sar.op21, ~. -sdSal.0m:poly(meanSST.1deg, 1))
anova(mod.sar.op22, mod.sar.op21)
summary(mod.sar.op22, Nagelkerke = TRUE) # 0.81701
AIC(mod.sar.op22) # 9617.086
rm(mod.sar.op21)

mod.sar.op23 <- update(mod.sar.op22, ~. -depth10deg:meanSal.0m)
anova(mod.sar.op23, mod.sar.op22)
summary(mod.sar.op23, Nagelkerke = TRUE) # 0.81692
AIC(mod.sar.op23) # 9616.269
rm(mod.sar.op22)

mod.sar.op24 <- update(mod.sar.op23, ~. -depth10deg:poly(meanSST.1deg, 2) + depth10deg:poly(meanSST.1deg, 1))
anova(mod.sar.op24, mod.sar.op23)
summary(mod.sar.op24, Nagelkerke = TRUE) # 0.81685
AIC(mod.sar.op24) # 9615.119
rm(mod.sar.op23)

mod.sar.op25 <- update(mod.sar.op24, ~. -depth10deg:poly(meanSST.1deg, 1))
anova(mod.sar.op25, mod.sar.op24)
summary(mod.sar.op25, Nagelkerke = TRUE) # 0.81675
AIC(mod.sar.op25) # 9614.297
rm(mod.sar.op24)

mod.sar.op26 <- update(mod.sar.op25, ~. -mean.mld.t:poly(meanSST.1deg, 2) + mean.mld.t:poly(meanSST.1deg, 1))
anova(mod.sar.op26, mod.sar.op25)
summary(mod.sar.op26, Nagelkerke = TRUE) # 0.81654
AIC(mod.sar.op26) # 9614.961
rm(mod.sar.op25)

mod.sar.op27 <- update(mod.sar.op26, ~. -mean.mld.t:Ocean2)
anova(mod.sar.op27, mod.sar.op26)
summary(mod.sar.op27, Nagelkerke = TRUE) # 0.8163
AIC(mod.sar.op27) # 9613.85
rm(mod.sar.op26)

mod.sar.op28 <- update(mod.sar.op27, ~. -mean.mld.t:logProd.mn.ann )
anova(mod.sar.op28, mod.sar.op27)
summary(mod.sar.op28, Nagelkerke = TRUE) # 0.81616
AIC(mod.sar.op28) # 9613.658
rm(mod.sar.op27)

mod.sar.op29 <- update(mod.sar.op28, ~. -poly(meanSST.1deg, 3):meanSal.0m + poly(meanSST.1deg, 2):meanSal.0m)
anova(mod.sar.op29, mod.sar.op28)
summary(mod.sar.op29, Nagelkerke = TRUE) # 0.81591
AIC(mod.sar.op29) # 9614.715
rm(mod.sar.op28)

mod.sar.op30 <- update(mod.sar.op29, ~. -sdSal.0m:carb_ion)
anova(mod.sar.op30, mod.sar.op29)
summary(mod.sar.op30, Nagelkerke = TRUE) # 0.81562
AIC(mod.sar.op30) # 9616.216
rm(mod.sar.op29)

mod.sar.opf <- mod.sar.op30
rm(mod.sar.op30)

## 2v. Create a plot of model parameters --------------------------------------
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.81562

# generate a dataframe of coefficients
(ms.coef <- data.frame(names = names(mod.sar.opf$coefficients), coef.sar = mod.sar.opf$coefficients, row.names = 1:length(mod.sar.opf$coefficients), stars = NA))

# reorder the rows to something more sensible
order.coef.ms <- c(1:19, 36:42, 20:35)
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
tmp.x <- barplot(abs(ms.coef$coef.sar), names = ms.coef$names, las = 2, ylim = c(0, max(ms.coef$coef.sar) + 50))
text(tmp.x, abs(ms.coef$coef.sar) + 50, ms.coef$stars)
par(plt = plt.def)
dev.off()
rm(plt.def, tmp.x, ms.coef)

## 2vi. Calculate likelihood ratios for the SAR model ----------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.81562
AIC(mod.sar.opf) # 9616.216

# removing mean temp^3
mod.sar.lr.mnt3 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 2) - poly(meanSST.1deg, 3):Ocean2 + poly(meanSST.1deg, 2):Ocean2)
summary(mod.sar.lr.mnt3, Nagelkerke = T) # r2 = 0.81123
AIC(mod.sar.lr.mnt3) # 9663.48
lrtest(mod.sar.opf, mod.sar.lr.mnt3) # n.b. this is basically an anova
# LR = 1.61e-11 ***

# removing mean temp^2
mod.sar.lr.mnt2 <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) + poly(meanSST.1deg, 1) - poly(meanSST.1deg, 3):Ocean2 + poly(meanSST.1deg, 1):Ocean2 - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 1) - logProd.mn.ann:poly(meanSST.1deg, 2) + logProd.mn.ann:poly(meanSST.1deg, 1) - meanSal.0m:poly(meanSST.1deg, 2) + meanSal.0m:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mnt2, Nagelkerke = T) # r2 = 0.79471
AIC(mod.sar.lr.mnt2) # 9841.466
lrtest(mod.sar.opf, mod.sar.lr.mnt2)
# LR = < 2.2e-16 ***

# removing mean temp
mod.sar.lr.mnt <- update(mod.sar.opf, ~. -poly(meanSST.1deg, 3) - poly(meanSST.1deg, 3):Ocean2 - carb_ion:poly(meanSST.1deg, 2) - logProd.mn.ann:poly(meanSST.1deg, 2) - meanSal.0m:poly(meanSST.1deg, 2) - mean.mld.t:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mnt, Nagelkerke = T) # r2 = 0.77267
AIC(mod.sar.lr.mnt) # 10058.4
lrtest(mod.sar.opf, mod.sar.lr.mnt)
# LR = < 2.2e-16 ***

# removing sd temp
mod.sar.lr.sdt <- update(mod.sar.opf, ~. -sdSST.1deg - sdSST.1deg:mean.mld.t - sdSST.1deg:depth10deg - sdSST.1deg:meanSal.0m - sdSST.1deg:Ocean2)
summary(mod.sar.lr.sdt, Nagelkerke = T) # r2 = 0.81298
AIC(mod.sar.lr.sdt) # 9636.458
lrtest(mod.sar.opf, mod.sar.lr.sdt)
# LR = 1.467e-05 ***

# removing mld temp
mod.sar.lr.mld <- update(mod.sar.opf, ~. -mean.mld.t - sdSST.1deg:mean.mld.t - mean.mld.t:poly(meanSST.1deg, 1))
summary(mod.sar.lr.mld, Nagelkerke = T) # r2 = 0.81447
AIC(mod.sar.lr.mld) # 9624.395
lrtest(mod.sar.opf, mod.sar.lr.mld)
# LR = 0.002672 **

# removing depth of 10 degree contour
mod.sar.lr.d10 <- update(mod.sar.opf, ~. -depth10deg - sdSST.1deg:depth10deg - depth10deg:sdSal.0m)
summary(mod.sar.lr.d10, Nagelkerke = T) # r2 = 0.81418
AIC(mod.sar.lr.d10) # 9627.893
lrtest(mod.sar.opf, mod.sar.lr.d10)
# LR = 0.0005128 ***

# removing mean log Prod
mod.sar.lr.prod <- update(mod.sar.opf, ~. -logProd.mn.ann - logProd.mn.ann:meanSal.0m - logProd.mn.ann:sdSal.0m - logProd.mn.ann:Ocean2 - logProd.mn.ann:carb_ion - logProd.mn.ann:poly(meanSST.1deg, 2))
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.80651
AIC(mod.sar.lr.prod) # 9709.413
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = < 2.2e-16 ***

# removing mean salinity
mod.sar.lr.msal <- update(mod.sar.opf, ~. -meanSal.0m - sdSST.1deg:meanSal.0m - logProd.mn.ann:meanSal.0m - meanSal.0m:Ocean2 - meanSal.0m:carb_ion - meanSal.0m:poly(meanSST.1deg, 2))
summary(mod.sar.lr.msal, Nagelkerke = T) # r2 = 0.81018
AIC(mod.sar.lr.msal) # 9666.071
lrtest(mod.sar.opf, mod.sar.lr.msal)
# LR = 3.27e-11 ***

# removing sd salinity
mod.sar.lr.sdsal <- update(mod.sar.opf, ~. -sdSal.0m - depth10deg:sdSal.0m - logProd.mn.ann:sdSal.0m)
summary(mod.sar.lr.sdsal, Nagelkerke = T) # r2 =  0.81354
AIC(mod.sar.lr.sdsal) # 9635.657
lrtest(mod.sar.opf, mod.sar.lr.sdsal)
# LR = 1.249e-05 ***

# removing Ocean2
mod.sar.lr.oce <- update(mod.sar.opf, ~. -Ocean2 - poly(meanSST.1deg, 3):Ocean2 - sdSST.1deg:Ocean2 - logProd.mn.ann:Ocean2 - meanSal.0m:Ocean2 - Ocean2:carb_ion )
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.80581
AIC(mod.sar.lr.oce) # 9701.616
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = < 2.2e-16 ***

# removing carb_ion
mod.sar.lr.dis <- update(mod.sar.opf, ~. -carb_ion - logProd.mn.ann:carb_ion - meanSal.0m:carb_ion - Ocean2:carb_ion - carb_ion:poly(meanSST.1deg, 2))
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.81381
AIC(mod.sar.lr.dis) # 9624.4
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.002362 **

# create a vector with these likelihood values
ms.lr <- data.frame(names = c("mnt3", "mnt2", "mnt", "sdt", "mld", "d10", "prod", "msal", "sdsal", "oce", "dis"), lr = NA, p = NA, stars = NA)

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
text(tmp.x, ms.lr$lr + 10, ms.lr$stars)
dev.off()
rm(tmp.x)

## 2vii. Calculate likelihood ratios for groups of EVs ---------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 0.81562
AIC(mod.sar.opf) # 9616.216

# removing temp
mod.sar.lr.temp <- mod.sar.lr.mnt
summary(mod.sar.lr.temp, Nagelkerke = T) # r2 = 0.77267
AIC(mod.sar.lr.temp) # 10058.4
lrtest(mod.sar.opf, mod.sar.lr.temp)
# LR = < 2.2e-16 ***

# removing structure
mod.sar.lr.str <- update(mod.sar.opf, ~. -mean.mld.t - depth10deg - sdSST.1deg:mean.mld.t - sdSST.1deg:depth10deg - depth10deg:sdSal.0m - mean.mld.t:poly(meanSST.1deg, 1))
summary(mod.sar.lr.str, Nagelkerke = T) # r2 = 0.81338
AIC(mod.sar.lr.str) # 9631.646
lrtest(mod.sar.opf, mod.sar.lr.str)
# LR = 0.0001203 ***

# removing stability
mod.sar.lr.stable <- update(mod.sar.opf, ~. -sdSST.1deg - sdSal.0m - sdSST.1deg:mean.mld.t - sdSST.1deg:depth10deg - sdSST.1deg:meanSal.0m - sdSST.1deg:Ocean2 - depth10deg:sdSal.0m - logProd.mn.ann:sdSal.0m)
summary(mod.sar.lr.stable, Nagelkerke = T) # r2 = 0.81119
AIC(mod.sar.lr.stable) # 9652.073
lrtest(mod.sar.opf, mod.sar.lr.stable)
# LR = 2.01e-08 ***

# removing productivity
summary(mod.sar.lr.prod, Nagelkerke = T) # r2 = 0.80651
AIC(mod.sar.lr.prod) # 9709.413
lrtest(mod.sar.opf, mod.sar.lr.prod)
# LR = < 2.2e-16 ***

# removing salinity
mod.sar.lr.sal <- mod.sar.lr.msal
summary(mod.sar.lr.sal, Nagelkerke = T) # r2 = 0.81018
AIC(mod.sar.lr.sal) # 9666.071
lrtest(mod.sar.opf, mod.sar.lr.sal)
# LR = 3.27e-11 ***

# removing Ocean2
summary(mod.sar.lr.oce, Nagelkerke = T) # r2 = 0.80581
AIC(mod.sar.lr.oce) # 9701.616
lrtest(mod.sar.opf, mod.sar.lr.oce)
# LR = < 2.2e-16 ***

# removing carb_ion
summary(mod.sar.lr.dis, Nagelkerke = T) # r2 = 0.81381
AIC(mod.sar.lr.dis) # 9624.4
lrtest(mod.sar.opf, mod.sar.lr.dis)
# LR = 0.002362 **

# create a vector with these likelihood values
ms.lr.group <- data.frame(names = c("temp", "str", "stable", "prod", "sal", "oce", "dis"), lr = NA, p = NA, stars = NA)

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
text(tmp.x, ms.lr.group$lr + 10, ms.lr.group$stars)
dev.off()
rm(tmp.x)

rm(mod.sar.lr.d10, mod.sar.lr.dis, mod.sar.lr.mld, mod.sar.lr.mnt, mod.sar.lr.mnt2, mod.sar.lr.mnt3, mod.sar.lr.msal, mod.sar.lr.oce, mod.sar.lr.prod, mod.sar.lr.sal, mod.sar.lr.sdsal, mod.sar.lr.sdt, mod.sar.lr.stable, mod.sar.lr.str, mod.sar.lr.temp)

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
lr.plot(lr.sar.op0, lr.sar.opf, order = c(7:9, 4, 2, 11:10, 3, 5:6, 1), leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio")
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.op0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.op0g <- lr.calc(mod.sar.op0, tmp)
rm(tmp)

(tmp <- data.frame(names = model.evs(mod.sar.opf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.opfg <- lr.calc(mod.sar.opf, tmp)
rm(tmp)

png("Figures/Ana_2viii_LRatio_g_opfop0.png", width = 800)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(6:7, 5, 3, 4, 1:2), leg.x = 17, leg.y = 275, leg.txt = c("Full", "Simplified"), ylab = "Log Likelihood ratio")
dev.off()

# test for heteroscedasticity
bptest.sarlm(mod.sar.op0) 

## 2ix. Does optimisation differ for simplified? ---------------------------
mod.fw.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fw.rop # 620.1094
AIC(mod.fw.rop$obj)
# 9616.738

mod.fb.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fb.rop # 620.6914
AIC(mod.fb.rop$obj)
# 9751.51

mod.fs.rop <- with(ldg.margo.mod, sar.optimised(mod.l0.sac$real$x.intercept, mod.sar.opf$call$formula, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
mod.fs.rop # 620.6914
AIC(mod.fs.rop$obj)
# 9649.923

# therefore justified in using w
rm(mod.fw.rop, mod.fb.rop, mod.fs.rop)


## 3. How do the different Oceans compare? ---------------------------------
# best model is
summary(mod.sar.opf, Nagelkerke = T) # r2 = 
AIC(mod.sar.opf) # 

# do I have suffient points for each ocean?
table(ldg.margo.mod$Ocean2) # should do

## 3i. set up model for Atlantic -----------------------------
# run model based on best model for complete data with only Atlantic
atl.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Atlantic", ], 0, mod.sar.opW$dist, longlat = TRUE)
atl.w <- nb2listw(atl.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(atl.nb)

## 3ii. set up model for Indian -----------------------------------
# run model based on best model for complete data with only Indian
ind.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Indian", ], 0, mod.sar.opW$dist, longlat = TRUE)
ind.w <- nb2listw(ind.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(ind.nb)

## 3iii. set up model for Pacific  ---------------------------
pac.nb <- dnearneigh(ldg.coords[ldg.margo.mod$Ocean2 == "Pacific", ], 0, mod.sar.opW$dist, longlat = TRUE)
pac.w <- nb2listw(pac.nb, glist = NULL, style = "W", zero.policy = TRUE)
rm(pac.nb)

# n.b. removed the simplification in 3i-3iii as I'm only simplifying from the model with significant interactions (see 3vii). Therefore, also didn't need 3iv - vi

## 3vii. only significant interactions for Atlantic -----------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)
# get formula without the Ocean2
op.formula <- update(mod.sar.opf$call$formula, ~.-Ocean2 - poly(meanSST.1deg, 3):Ocean2 - sdSST.1deg:Ocean2 - logProd.mn.ann:Ocean2 - meanSal.0m:Ocean2 - Ocean2:carb_ion)

# create a model with only significant values
mod.sar.atlI <- errorsarlm(op.formula, listw = atl.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", ])
sar.plot(mod.sar.atlI)

# try adding in the other interactions

mod.sar.atlI2 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI2)
anova(mod.sar.atlI, mod.sar.atlI2) # 0.83789
rm(mod.sar.atlI2)

mod.sar.atlI3 <- update(mod.sar.atlI, ~. -poly(meanSST.1deg, 1):mean.mld.t + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.atlI3)
anova(mod.sar.atlI, mod.sar.atlI3) # 0.8938
rm(mod.sar.atlI3)

mod.sar.atlI4 <- update(mod.sar.atlI, ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.atlI4)
anova(mod.sar.atlI, mod.sar.atlI4) # 0.82195
rm(mod.sar.atlI4)

mod.sar.atlI5 <- update(mod.sar.atlI, ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI5)
anova(mod.sar.atlI, mod.sar.atlI5) # 0.32516
rm(mod.sar.atlI5)

# significant at the 0.05 level
mod.sar.atlI6 <- update(mod.sar.atlI, ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.atlI6)
anova(mod.sar.atlI, mod.sar.atlI6) # 0.0020231

mod.sar.atlI7 <- update(mod.sar.atlI, ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.atlI7)
anova(mod.sar.atlI, mod.sar.atlI7) # 0.50761
rm(mod.sar.atlI7)

mod.sar.atlI8 <- update(mod.sar.atlI, ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atlI8)
anova(mod.sar.atlI, mod.sar.atlI8) # 0.53505
rm(mod.sar.atlI8)

mod.sar.atlI9 <- update(mod.sar.atlI, ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.atlI9)
anova(mod.sar.atlI, mod.sar.atlI9) # 0.23541
rm(mod.sar.atlI9)

# significant at the 0.05 level
mod.sar.atlI10 <- update(mod.sar.atlI, ~. + sdSST.1deg:carb_ion)
summary(mod.sar.atlI10)
anova(mod.sar.atlI, mod.sar.atlI10) # 0.041153
rm(mod.sar.atlI10)

mod.sar.atlI11 <- update(mod.sar.atlI, ~. + mean.mld.t:depth10deg)
summary(mod.sar.atlI11)
anova(mod.sar.atlI, mod.sar.atlI11) # 0.94255
rm(mod.sar.atlI11)

mod.sar.atlI12 <- update(mod.sar.atlI, ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.atlI12)
anova(mod.sar.atlI, mod.sar.atlI12) # 0.85666
rm(mod.sar.atlI12)

mod.sar.atlI13 <- update(mod.sar.atlI, ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.atlI13)
anova(mod.sar.atlI, mod.sar.atlI13) # 0.16798
rm(mod.sar.atlI13)

mod.sar.atlI14 <- update(mod.sar.atlI, ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.atlI14)
anova(mod.sar.atlI, mod.sar.atlI14) # 0.86732
rm(mod.sar.atlI14)

mod.sar.atlI15 <- update(mod.sar.atlI, ~. + mean.mld.t:carb_ion )
summary(mod.sar.atlI15)
anova(mod.sar.atlI, mod.sar.atlI15) # 0.28498
rm(mod.sar.atlI15)

# significant at the 0.05 level
mod.sar.atlI16 <- update(mod.sar.atlI, ~. + depth10deg:logProd.mn.ann )
summary(mod.sar.atlI16)
anova(mod.sar.atlI, mod.sar.atlI16) # 0.011923
rm(mod.sar.atlI16)

mod.sar.atlI17 <- update(mod.sar.atlI, ~. + depth10deg:meanSal.0m)
summary(mod.sar.atlI17)
anova(mod.sar.atlI, mod.sar.atlI17) # 0.58715
rm(mod.sar.atlI17)

mod.sar.atlI18 <- update(mod.sar.atlI, ~. + depth10deg:carb_ion )
summary(mod.sar.atlI18)
anova(mod.sar.atlI, mod.sar.atlI18) # 0.8071
rm(mod.sar.atlI18)

mod.sar.atlI19 <- update(mod.sar.atlI, ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.atlI19)
anova(mod.sar.atlI, mod.sar.atlI19) # 0.1013
rm(mod.sar.atlI19)

mod.sar.atlI20 <- update(mod.sar.atlI, ~. + sdSal.0m:carb_ion)
summary(mod.sar.atlI20)
anova(mod.sar.atlI, mod.sar.atlI20) # 0.20287
rm(mod.sar.atlI20)

## adding back in finds poly(meanSST.1deg, 2):meanSal.0m should be included in the model, also potentially sdSST.1deg:carb_ion and depth10deg:logProd.mn.ann
mod.sar.atlI1. <- mod.sar.atlI6
rm(mod.sar.atlI6)

# check whether any of the interactions are now significant
mod.sar.atlI1.1 <- update(mod.sar.atlI1., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI1.1)
anova(mod.sar.atlI1., mod.sar.atlI1.1) # 0.57295
rm(mod.sar.atlI1.1)

mod.sar.atlI1.2 <- update(mod.sar.atlI1., ~. -mean.mld.t:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.atlI1.2)
anova(mod.sar.atlI1., mod.sar.atlI1.2) # 0.64088
rm(mod.sar.atlI1.2)

mod.sar.atlI1.3 <- update(mod.sar.atlI1., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.atlI1.3)
anova(mod.sar.atlI1., mod.sar.atlI1.3) # 0.10308
rm(mod.sar.atlI1.3)

mod.sar.atlI1.5 <- update(mod.sar.atlI1., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI1.5)
anova(mod.sar.atlI1., mod.sar.atlI1.5) # 0.64746
rm(mod.sar.atlI1.5)

mod.sar.atlI1.7 <- update(mod.sar.atlI1., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.atlI1.7)
anova(mod.sar.atlI1., mod.sar.atlI1.7) # 0.55092
rm(mod.sar.atlI1.7)

mod.sar.atlI1.8 <- update(mod.sar.atlI1., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atlI1.8)
anova(mod.sar.atlI1., mod.sar.atlI1.8) # 0.59264
rm(mod.sar.atlI1.8)

mod.sar.atlI1.9 <- update(mod.sar.atlI1., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.atlI1.9)
anova(mod.sar.atlI1., mod.sar.atlI1.9) # 0.26767
rm(mod.sar.atlI1.9)

# significant at the 0.05 level
mod.sar.atlI1.10 <- update(mod.sar.atlI1., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.atlI1.10)
anova(mod.sar.atlI1., mod.sar.atlI1.10) # 0.017228
rm(mod.sar.atlI1.10)

mod.sar.atlI1.11 <- update(mod.sar.atlI1., ~. + mean.mld.t:depth10deg)
summary(mod.sar.atlI1.11)
anova(mod.sar.atlI1., mod.sar.atlI1.11) # 0.62318
rm(mod.sar.atlI1.11)

mod.sar.atlI1.12 <- update(mod.sar.atlI1., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.atlI1.12)
anova(mod.sar.atlI1., mod.sar.atlI1.12) # 0.24256
rm(mod.sar.atlI1.12)

mod.sar.atlI1.13 <- update(mod.sar.atlI1., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.atlI1.13)
anova(mod.sar.atlI1., mod.sar.atlI1.13) # 0.46506
rm(mod.sar.atlI1.13)

mod.sar.atlI1.14 <- update(mod.sar.atlI1., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.atlI1.14)
anova(mod.sar.atlI1., mod.sar.atlI1.14) # 0.69291
rm(mod.sar.atlI1.14)

mod.sar.atlI1.15 <- update(mod.sar.atlI1., ~. + mean.mld.t:carb_ion )
summary(mod.sar.atlI1.15)
anova(mod.sar.atlI1., mod.sar.atlI1.15) # 0.35889
rm(mod.sar.atlI1.15)

# significant at the 0.05 level
mod.sar.atlI1.16 <- update(mod.sar.atlI1., ~. + depth10deg:logProd.mn.ann )
summary(mod.sar.atlI1.16)
anova(mod.sar.atlI1., mod.sar.atlI1.16) # 0.005155

mod.sar.atlI1.17 <- update(mod.sar.atlI1., ~. + depth10deg:meanSal.0m)
summary(mod.sar.atlI1.17)
anova(mod.sar.atlI1., mod.sar.atlI1.17) # 0.34126
rm(mod.sar.atlI1.17)

mod.sar.atlI1.18 <- update(mod.sar.atlI1., ~. + depth10deg:carb_ion )
summary(mod.sar.atlI1.18)
anova(mod.sar.atlI1., mod.sar.atlI1.18) # 0.97048
rm(mod.sar.atlI1.18)

mod.sar.atlI1.19 <- update(mod.sar.atlI1., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.atlI1.19)
anova(mod.sar.atlI1., mod.sar.atlI1.19) # 0.47582
rm(mod.sar.atlI1.19)

mod.sar.atlI1.20 <- update(mod.sar.atlI1., ~. + sdSal.0m:carb_ion)
summary(mod.sar.atlI1.20)
anova(mod.sar.atlI1., mod.sar.atlI1.20) # 0.142
rm(mod.sar.atlI1.20)

## adding back in finds depth10deg:logProd.mn.ann should be included in the model, also potentially sdSST.1deg:carb_ion
mod.sar.atlI2. <- mod.sar.atlI1.16
rm(mod.sar.atlI1.16)

# check whether any of the interactions are now significant
mod.sar.atlI2.1 <- update(mod.sar.atlI2., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI2.1)
anova(mod.sar.atlI2., mod.sar.atlI2.1) # 0.49526
rm(mod.sar.atlI2.1)

mod.sar.atlI2.2 <- update(mod.sar.atlI2., ~. -mean.mld.t:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.atlI2.2)
anova(mod.sar.atlI2., mod.sar.atlI2.2) # 0.49451
rm(mod.sar.atlI2.2)

mod.sar.atlI2.3 <- update(mod.sar.atlI2., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.atlI2.3)
anova(mod.sar.atlI2., mod.sar.atlI2.3) # 0.18022
rm(mod.sar.atlI2.3)

mod.sar.atlI2.5 <- update(mod.sar.atlI2., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI2.5)
anova(mod.sar.atlI2., mod.sar.atlI2.5) # 0.97819
rm(mod.sar.atlI2.5)

mod.sar.atlI2.7 <- update(mod.sar.atlI2., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.atlI2.7)
anova(mod.sar.atlI2., mod.sar.atlI2.7) # 0.53435
rm(mod.sar.atlI2.7)

mod.sar.atlI2.8 <- update(mod.sar.atlI2., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atlI2.8)
anova(mod.sar.atlI2., mod.sar.atlI2.8) # 0.69915
rm(mod.sar.atlI2.8)

mod.sar.atlI2.9 <- update(mod.sar.atlI2., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.atlI2.9)
anova(mod.sar.atlI2., mod.sar.atlI2.9) # 0.37692
rm(mod.sar.atlI2.9)

# significant at the 0.05 level
mod.sar.atlI2.10 <- update(mod.sar.atlI2., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.atlI2.10)
anova(mod.sar.atlI2., mod.sar.atlI2.10) # 0.025392

mod.sar.atlI2.11 <- update(mod.sar.atlI2., ~. + mean.mld.t:depth10deg)
summary(mod.sar.atlI2.11)
anova(mod.sar.atlI2., mod.sar.atlI2.11) # 0.30343
rm(mod.sar.atlI2.11)

mod.sar.atlI2.12 <- update(mod.sar.atlI2., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.atlI2.12)
anova(mod.sar.atlI2., mod.sar.atlI2.12) # 0.25065
rm(mod.sar.atlI2.12)

mod.sar.atlI2.13 <- update(mod.sar.atlI2., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.atlI2.13)
anova(mod.sar.atlI2., mod.sar.atlI2.13) # 0.44838
rm(mod.sar.atlI2.13)

mod.sar.atlI2.14 <- update(mod.sar.atlI2., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.atlI2.14)
anova(mod.sar.atlI2., mod.sar.atlI2.14) # 0.73246
rm(mod.sar.atlI2.14)

mod.sar.atlI2.15 <- update(mod.sar.atlI2., ~. + mean.mld.t:carb_ion )
summary(mod.sar.atlI2.15)
anova(mod.sar.atlI2., mod.sar.atlI2.15) # 0.31739
rm(mod.sar.atlI2.15)

mod.sar.atlI2.17 <- update(mod.sar.atlI2., ~. + depth10deg:meanSal.0m)
summary(mod.sar.atlI2.17)
anova(mod.sar.atlI2., mod.sar.atlI2.17) # 0.69122
rm(mod.sar.atlI2.17)

mod.sar.atlI2.18 <- update(mod.sar.atlI2., ~. + depth10deg:carb_ion )
summary(mod.sar.atlI2.18)
anova(mod.sar.atlI2., mod.sar.atlI2.18) # 0.78031
rm(mod.sar.atlI2.18)

mod.sar.atlI2.19 <- update(mod.sar.atlI2., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.atlI2.19)
anova(mod.sar.atlI2., mod.sar.atlI2.19) # 0.65969
rm(mod.sar.atlI2.19)

mod.sar.atlI2.20 <- update(mod.sar.atlI2., ~. + sdSal.0m:carb_ion)
summary(mod.sar.atlI2.20)
anova(mod.sar.atlI2., mod.sar.atlI2.20) # 0.12301
rm(mod.sar.atlI2.20)

## adding back in finds sdSST.1deg:carb_ion should be included in the model
mod.sar.atlI3. <- mod.sar.atlI2.10
rm(mod.sar.atlI2.10)

# check whether any of the interactions are now significant
mod.sar.atlI3.1 <- update(mod.sar.atlI3., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.atlI3.1)
anova(mod.sar.atlI3., mod.sar.atlI3.1) # 0.39792
rm(mod.sar.atlI3.1)

mod.sar.atlI3.2 <- update(mod.sar.atlI3., ~. -mean.mld.t:poly(meanSST.1deg, 1) + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.atlI3.2)
anova(mod.sar.atlI3., mod.sar.atlI3.2) # 0.64641
rm(mod.sar.atlI3.2)

mod.sar.atlI3.3 <- update(mod.sar.atlI3., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.atlI3.3)
anova(mod.sar.atlI3., mod.sar.atlI3.3) # 0.15599
rm(mod.sar.atlI3.3)

mod.sar.atlI3.5 <- update(mod.sar.atlI3., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.atlI3.5)
anova(mod.sar.atlI3., mod.sar.atlI3.5) # 0.79261
rm(mod.sar.atlI3.5)

mod.sar.atlI3.7 <- update(mod.sar.atlI3., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.atlI3.7)
anova(mod.sar.atlI3., mod.sar.atlI3.7) # 0.62064
rm(mod.sar.atlI3.7)

mod.sar.atlI3.8 <- update(mod.sar.atlI3., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.atlI3.8)
anova(mod.sar.atlI3., mod.sar.atlI3.8) # 0.2012
rm(mod.sar.atlI3.8)

mod.sar.atlI3.9 <- update(mod.sar.atlI3., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.atlI3.9)
anova(mod.sar.atlI3., mod.sar.atlI3.9) # 0.18594
rm(mod.sar.atlI3.9)

mod.sar.atlI3.11 <- update(mod.sar.atlI3., ~. + mean.mld.t:depth10deg)
summary(mod.sar.atlI3.11)
anova(mod.sar.atlI3., mod.sar.atlI3.11) # 0.29551
rm(mod.sar.atlI3.11)

mod.sar.atlI3.12 <- update(mod.sar.atlI3., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.atlI3.12)
anova(mod.sar.atlI3., mod.sar.atlI3.12) # 0.34265
rm(mod.sar.atlI3.12)

mod.sar.atlI3.13 <- update(mod.sar.atlI3., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.atlI3.13)
anova(mod.sar.atlI3., mod.sar.atlI3.13) # 0.66359
rm(mod.sar.atlI3.13)

mod.sar.atlI3.14 <- update(mod.sar.atlI3., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.atlI3.14)
anova(mod.sar.atlI3., mod.sar.atlI3.14) # 0.56702
rm(mod.sar.atlI3.14)

mod.sar.atlI3.15 <- update(mod.sar.atlI3., ~. + mean.mld.t:carb_ion )
summary(mod.sar.atlI3.15)
anova(mod.sar.atlI3., mod.sar.atlI3.15) # 0.91687
rm(mod.sar.atlI3.15)

mod.sar.atlI3.17 <- update(mod.sar.atlI3., ~. + depth10deg:meanSal.0m)
summary(mod.sar.atlI3.17)
anova(mod.sar.atlI3., mod.sar.atlI3.17) # 0.60124
rm(mod.sar.atlI3.17)

mod.sar.atlI3.18 <- update(mod.sar.atlI3., ~. + depth10deg:carb_ion )
summary(mod.sar.atlI3.18)
anova(mod.sar.atlI3., mod.sar.atlI3.18) # 0.95038
rm(mod.sar.atlI3.18)

mod.sar.atlI3.19 <- update(mod.sar.atlI3., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.atlI3.19)
anova(mod.sar.atlI3., mod.sar.atlI3.19) # 0.72213
rm(mod.sar.atlI3.19)

mod.sar.atlI3.20 <- update(mod.sar.atlI3., ~. + sdSal.0m:carb_ion)
summary(mod.sar.atlI3.20)
anova(mod.sar.atlI3., mod.sar.atlI3.20) # 0.34207
rm(mod.sar.atlI3.20)

rm(mod.sar.atl)

# nothing is significant, so 
mod.sar.atlIf <- mod.sar.atlI3.
summary(mod.sar.atlIf, Nagelkerke = TRUE) # 0.89508 
AIC(mod.sar.atlIf) # 4681.261

# Check for correlation
env.var.atl <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "meanSal.0m", "sdSal.0m", "carb_ion")

# pairs plot
png("Figures/Ana_3vii_atl_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", names(ldg.margo.mod) %in% env.var.atl])

lr.sar.atlIf <- lr.calc(mod.sar.atlIf)

png("Figures/Ana_3vii_LRatio_atlIf.png", width = 800)
lr.plot(lr.sar.atlIf, order = c(8:10, 2:3, 1, 6, 4:5, 7), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.atlIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.atlIfg <- lr.calc(mod.sar.atlIf, tmp)
rm(tmp)

png("Figures/Ana_3vii_LRatio_g_atlIf.png", width = 800)
lr.plot(lr.sar.atlIfg, order = c(6, 2, 1, 3:5), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.atlI, mod.sar.atlI1., mod.sar.atlI2., mod.sar.atlI3., file = "Outputs/Atlantic_simplification.RData")
rm(env.var.atl, mod.sar.atlI, mod.sar.atlI1., mod.sar.atlI2., mod.sar.atlI3.)

## 3viii. Only significant interactions for Indian -------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.indI <- errorsarlm(op.formula, listw = ind.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", ])
sar.plot(mod.sar.indI)
summary(mod.sar.indI, Nagelkerke = TRUE) # 0.789 
AIC(mod.sar.indI) # 1516.811

# try adding in the other interactions
mod.sar.indI2 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.indI2)
anova(mod.sar.indI, mod.sar.indI2) # 0.70595
rm(mod.sar.indI2)

mod.sar.indI3 <- update(mod.sar.indI, ~. -poly(meanSST.1deg, 1):mean.mld.t + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.indI3)
anova(mod.sar.indI, mod.sar.indI3) # 0.097567
rm(mod.sar.indI3)

mod.sar.indI4 <- update(mod.sar.indI, ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.indI4)
anova(mod.sar.indI, mod.sar.indI4) # 0.5006
rm(mod.sar.indI4)

mod.sar.indI5 <- update(mod.sar.indI, ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.indI5)
anova(mod.sar.indI, mod.sar.indI5) # 0.79888
rm(mod.sar.indI5)

mod.sar.indI6 <- update(mod.sar.indI, ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.indI6)
anova(mod.sar.indI, mod.sar.indI6) # 0.42223
rm(mod.sar.indI6)

mod.sar.indI7 <- update(mod.sar.indI, ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.indI7)
anova(mod.sar.indI, mod.sar.indI7) # 0.31537
rm(mod.sar.indI7)

# significant at the 0.05 level
mod.sar.indI8 <- update(mod.sar.indI, ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.indI8)
anova(mod.sar.indI, mod.sar.indI8) # 0.0016597

mod.sar.indI9 <- update(mod.sar.indI, ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.indI9)
anova(mod.sar.indI, mod.sar.indI9) # 0.17245
rm(mod.sar.indI9)

mod.sar.indI10 <- update(mod.sar.indI, ~. + sdSST.1deg:carb_ion)
summary(mod.sar.indI10)
anova(mod.sar.indI, mod.sar.indI10) # 0.55745
rm(mod.sar.indI10)

mod.sar.indI11 <- update(mod.sar.indI, ~. + mean.mld.t:depth10deg)
summary(mod.sar.indI11)
anova(mod.sar.indI, mod.sar.indI11) # 0.63037
rm(mod.sar.indI11)

mod.sar.indI12 <- update(mod.sar.indI, ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.indI12)
anova(mod.sar.indI, mod.sar.indI12) # 0.3051
rm(mod.sar.indI12)

mod.sar.indI13 <- update(mod.sar.indI, ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.indI13)
anova(mod.sar.indI, mod.sar.indI13) # 0.14565
rm(mod.sar.indI13)

mod.sar.indI14 <- update(mod.sar.indI, ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.indI14)
anova(mod.sar.indI, mod.sar.indI14) # 0.76213
rm(mod.sar.indI14)

mod.sar.indI15 <- update(mod.sar.indI, ~. + mean.mld.t:carb_ion )
summary(mod.sar.indI15)
anova(mod.sar.indI, mod.sar.indI15) # 0.75948
rm(mod.sar.indI15)

mod.sar.indI16 <- update(mod.sar.indI, ~. + depth10deg:logProd.mn.ann )
summary(mod.sar.indI16)
anova(mod.sar.indI, mod.sar.indI16) # 0.45732
rm(mod.sar.indI16)

mod.sar.indI17 <- update(mod.sar.indI, ~. + depth10deg:meanSal.0m)
summary(mod.sar.indI17)
anova(mod.sar.indI, mod.sar.indI17) # 0.62925
rm(mod.sar.indI17)

mod.sar.indI18 <- update(mod.sar.indI, ~. + depth10deg:carb_ion )
summary(mod.sar.indI18)
anova(mod.sar.indI, mod.sar.indI18) # 0.39993
rm(mod.sar.indI18)

mod.sar.indI19 <- update(mod.sar.indI, ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.indI19)
anova(mod.sar.indI, mod.sar.indI19) # 0.20655
rm(mod.sar.indI19)

mod.sar.indI20 <- update(mod.sar.indI, ~. + sdSal.0m:carb_ion)
summary(mod.sar.indI20)
anova(mod.sar.indI, mod.sar.indI20) # 0.1801
rm(mod.sar.indI20)

# carb_ion:poly(meanSST.1deg, 3) is significant at the 0.05 level, so add it in and recheck
mod.sar.indI1. <- mod.sar.indI8
rm(mod.sar.indI8)

mod.sar.indI1.2 <- update(mod.sar.indI1., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.indI1.2)
anova(mod.sar.indI1., mod.sar.indI1.2) # 0.77728
rm(mod.sar.indI1.2)

mod.sar.indI1.3 <- update(mod.sar.indI1., ~. -poly(meanSST.1deg, 1):mean.mld.t + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.indI1.3)
anova(mod.sar.indI1., mod.sar.indI1.3) # 0.59566
rm(mod.sar.indI1.3)

mod.sar.indI1.4 <- update(mod.sar.indI1., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.indI1.4)
anova(mod.sar.indI1., mod.sar.indI1.4) # 0.4446
rm(mod.sar.indI1.4)

mod.sar.indI1.5 <- update(mod.sar.indI1., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.indI1.5)
anova(mod.sar.indI1., mod.sar.indI1.5) # 0.97587
rm(mod.sar.indI1.5)

mod.sar.indI1.6 <- update(mod.sar.indI1., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.indI1.6)
anova(mod.sar.indI1., mod.sar.indI1.6) # 0.85363
rm(mod.sar.indI1.6)

mod.sar.indI1.7 <- update(mod.sar.indI1., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.indI1.7)
anova(mod.sar.indI1., mod.sar.indI1.7) # 0.74069
rm(mod.sar.indI1.7)

mod.sar.indI1.9 <- update(mod.sar.indI1., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.indI1.9)
anova(mod.sar.indI1., mod.sar.indI1.9) # 0.13813
rm(mod.sar.indI1.9)

mod.sar.indI1.10 <- update(mod.sar.indI1., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.indI1.10)
anova(mod.sar.indI1., mod.sar.indI1.10) # 0.53888
rm(mod.sar.indI1.10)

mod.sar.indI1.11 <- update(mod.sar.indI1., ~. + mean.mld.t:depth10deg)
summary(mod.sar.indI1.11)
anova(mod.sar.indI1., mod.sar.indI1.11) # 0.56359
rm(mod.sar.indI1.11)

mod.sar.indI1.12 <- update(mod.sar.indI1., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.indI1.12)
anova(mod.sar.indI1., mod.sar.indI1.12) # 0.67234
rm(mod.sar.indI1.12)

mod.sar.indI1.13 <- update(mod.sar.indI1., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.indI1.13)
anova(mod.sar.indI1., mod.sar.indI1.13) # 0.32334
rm(mod.sar.indI1.13)

mod.sar.indI1.14 <- update(mod.sar.indI1., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.indI1.14)
anova(mod.sar.indI1., mod.sar.indI1.14) # 0.7707
rm(mod.sar.indI1.14)

mod.sar.indI1.15 <- update(mod.sar.indI1., ~. + mean.mld.t:carb_ion )
summary(mod.sar.indI1.15)
anova(mod.sar.indI1., mod.sar.indI1.15) # 0.093634
rm(mod.sar.indI1.15)

mod.sar.indI1.16 <- update(mod.sar.indI1., ~. + depth10deg:logProd.mn.ann )
summary(mod.sar.indI1.16)
anova(mod.sar.indI1., mod.sar.indI1.16) # 0.52883
rm(mod.sar.indI1.16)

mod.sar.indI1.17 <- update(mod.sar.indI1., ~. + depth10deg:meanSal.0m)
summary(mod.sar.indI1.17)
anova(mod.sar.indI1., mod.sar.indI1.17) # 0.59609
rm(mod.sar.indI1.17)

mod.sar.indI1.18 <- update(mod.sar.indI1., ~. + depth10deg:carb_ion )
summary(mod.sar.indI1.18)
anova(mod.sar.indI1., mod.sar.indI1.18) # 0.48651
rm(mod.sar.indI1.18)

mod.sar.indI1.19 <- update(mod.sar.indI1., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.indI1.19)
anova(mod.sar.indI1., mod.sar.indI1.19) # 0.29317
rm(mod.sar.indI1.19)

mod.sar.indI1.20 <- update(mod.sar.indI1., ~. + sdSal.0m:carb_ion)
summary(mod.sar.indI1.20)
anova(mod.sar.indI1., mod.sar.indI1.20) # 0.41853
rm(mod.sar.indI1.20)

## nothing else needs adding in 
mod.sar.indIf <- mod.sar.indI1.
summary(mod.sar.indIf, Nagelkerke = TRUE) # 0.79458
AIC(mod.sar.indIf) # 1508.919

# Check for correlation
env.var.ind <- c("meanSST.1deg", "sdSST.1deg", "mean.mld.t", "depth10deg", "logProd.mn.ann", "meanSal.0m", "sdSal.0m", "carb_ion")

# pairs plot
png("Figures/Ana_3vii_ind_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", names(ldg.margo.mod) %in% env.var.ind])

lr.sar.indIf <- lr.calc(mod.sar.indIf)

png("Figures/Ana_3viii_LRatio_indIf.png", width = 800)
lr.plot(lr.sar.indIf, order = c(8:10, 2:3, 1, 6, 4:5, 7), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.indIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.indIfg <- lr.calc(mod.sar.indIf, tmp)
rm(tmp)

png("Figures/Ana_3viii_LRatio_g_indIf.png", width = 800)
lr.plot(lr.sar.indIfg, order = c(6, 2, 1, 3:5), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 3)
dev.off()

save(mod.sar.indI, mod.sar.indI1., file = "Outputs/Indian_simplification.RData")
rm(ind.w, env.var.ind, mod.sar.indI, mod.sar.indI1.)

## 3ix. Only significant interactions for Pacific --------------------------
summary(mod.sar.op0)
summary(mod.sar.opf)

# create a model with only significant values
mod.sar.pacI <- errorsarlm(op.formula, listw = pac.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", ])
sar.plot(mod.sar.pacI)

# try adding in the other interactions
mod.sar.pacI2 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI2)
anova(mod.sar.pacI, mod.sar.pacI2) # 0.28549
rm(mod.sar.pacI2)

# significant at the 0.05 level
mod.sar.pacI3 <- update(mod.sar.pacI, ~. -poly(meanSST.1deg, 1):mean.mld.t + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.pacI3)
anova(mod.sar.pacI, mod.sar.pacI3) # 0.002242
rm(mod.sar.pacI3)

mod.sar.pacI4 <- update(mod.sar.pacI, ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI4)
anova(mod.sar.pacI, mod.sar.pacI4) # 0.16161
rm(mod.sar.pacI4)

mod.sar.pacI5 <- update(mod.sar.pacI, ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pacI5)
anova(mod.sar.pacI, mod.sar.pacI5) # 0.23137
rm(mod.sar.pacI5)

mod.sar.pacI6 <- update(mod.sar.pacI, ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI6)
anova(mod.sar.pacI, mod.sar.pacI6) # 0.79634
rm(mod.sar.pacI6)

mod.sar.pacI7 <- update(mod.sar.pacI, ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI7)
anova(mod.sar.pacI, mod.sar.pacI7) # 0.60898
rm(mod.sar.pacI7)

mod.sar.pacI8 <- update(mod.sar.pacI, ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI8)
anova(mod.sar.pacI, mod.sar.pacI8) # 0.33694
rm(mod.sar.pacI8)

# significant at the 0.05 level
mod.sar.pacI9 <- update(mod.sar.pacI, ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.pacI9)
anova(mod.sar.pacI, mod.sar.pacI9) # 0.0051662
rm(mod.sar.pacI9)

mod.sar.pacI10 <- update(mod.sar.pacI, ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI10)
anova(mod.sar.pacI, mod.sar.pacI10) # 0.22407
rm(mod.sar.pacI10)

mod.sar.pacI11 <- update(mod.sar.pacI, ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI11)
anova(mod.sar.pacI, mod.sar.pacI11) # 0.12321
rm(mod.sar.pacI11)

mod.sar.pacI12 <- update(mod.sar.pacI, ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI12)
anova(mod.sar.pacI, mod.sar.pacI12) # 0.32431
rm(mod.sar.pacI12)

mod.sar.pacI13 <- update(mod.sar.pacI, ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI13)
anova(mod.sar.pacI, mod.sar.pacI13) # 0.64123
rm(mod.sar.pacI13)

mod.sar.pacI14 <- update(mod.sar.pacI, ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI14)
anova(mod.sar.pacI, mod.sar.pacI14) # 0.65683
rm(mod.sar.pacI14)

mod.sar.pacI15 <- update(mod.sar.pacI, ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI15)
anova(mod.sar.pacI, mod.sar.pacI15) # 0.29536
rm(mod.sar.pacI15)

# significant at the 0.05 level
mod.sar.pacI16 <- update(mod.sar.pacI, ~. + depth10deg:logProd.mn.ann )
summary(mod.sar.pacI16)
anova(mod.sar.pacI, mod.sar.pacI16) # 0.00045056

# significant at the 0.05 level
mod.sar.pacI17 <- update(mod.sar.pacI, ~. + depth10deg:meanSal.0m)
summary(mod.sar.pacI17)
anova(mod.sar.pacI, mod.sar.pacI17) # 0.03202
rm(mod.sar.pacI17)

mod.sar.pacI18 <- update(mod.sar.pacI, ~. + depth10deg:carb_ion )
summary(mod.sar.pacI18)
anova(mod.sar.pacI, mod.sar.pacI18) # 0.53884
rm(mod.sar.pacI18)

mod.sar.pacI19 <- update(mod.sar.pacI, ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI19)
anova(mod.sar.pacI, mod.sar.pacI19) # 0.81047
rm(mod.sar.pacI19)

mod.sar.pacI20 <- update(mod.sar.pacI, ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI20)
anova(mod.sar.pacI, mod.sar.pacI20) # 0.33898
rm(mod.sar.pacI20)

# multiple interactions were significant. Add in the most significant and then recalculate
mod.sar.pacI1. <- mod.sar.pacI16
rm(mod.sar.pacI16)

# try adding in the other interactions
mod.sar.pacI1.2 <- update(mod.sar.pacI1., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI1.2)
anova(mod.sar.pacI1., mod.sar.pacI1.2) # 0.35635
rm(mod.sar.pacI1.2)

# significant at the 0.05 level
mod.sar.pacI1.3 <- update(mod.sar.pacI1., ~. -poly(meanSST.1deg, 1):mean.mld.t + poly(meanSST.1deg, 3):mean.mld.t)
summary(mod.sar.pacI1.3)
anova(mod.sar.pacI1., mod.sar.pacI1.3) # 0.0062143

# significant at the 0.05 level
mod.sar.pacI1.4 <- update(mod.sar.pacI1., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI1.4)
anova(mod.sar.pacI1., mod.sar.pacI1.4) # 0.031591
rm(mod.sar.pacI1.4)

# significant at the 0.05 level
mod.sar.pacI1.5 <- update(mod.sar.pacI1., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pacI1.5)
anova(mod.sar.pacI1., mod.sar.pacI1.5) # 0.027259
rm(mod.sar.pacI1.5)

mod.sar.pacI1.6 <- update(mod.sar.pacI1., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI1.6)
anova(mod.sar.pacI1., mod.sar.pacI1.6) # 0.86779
rm(mod.sar.pacI1.6)

mod.sar.pacI1.7 <- update(mod.sar.pacI1., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI1.7)
anova(mod.sar.pacI1., mod.sar.pacI1.7) # 0.65464
rm(mod.sar.pacI1.7)

mod.sar.pacI1.8 <- update(mod.sar.pacI1., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI1.8)
anova(mod.sar.pacI1., mod.sar.pacI1.8) # 0.63417
rm(mod.sar.pacI1.8)

mod.sar.pacI1.9 <- update(mod.sar.pacI1., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.pacI1.9)
anova(mod.sar.pacI1., mod.sar.pacI1.9) # 0.052901
rm(mod.sar.pacI1.9)

mod.sar.pacI1.10 <- update(mod.sar.pacI1., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI1.10)
anova(mod.sar.pacI1., mod.sar.pacI1.10) # 0.3855
rm(mod.sar.pacI1.10)

mod.sar.pacI1.11 <- update(mod.sar.pacI1., ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI1.11)
anova(mod.sar.pacI1., mod.sar.pacI1.11) # 0.2786
rm(mod.sar.pacI1.11)

mod.sar.pacI1.12 <- update(mod.sar.pacI1., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI1.12)
anova(mod.sar.pacI1., mod.sar.pacI1.12) # 0.36162
rm(mod.sar.pacI1.12)

mod.sar.pacI1.13 <- update(mod.sar.pacI1., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI1.13)
anova(mod.sar.pacI1., mod.sar.pacI1.13) # 0.93635
rm(mod.sar.pacI1.13)

mod.sar.pacI1.14 <- update(mod.sar.pacI1., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI1.14)
anova(mod.sar.pacI1., mod.sar.pacI1.14) # 0.89079
rm(mod.sar.pacI1.14)

mod.sar.pacI1.15 <- update(mod.sar.pacI1., ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI1.15)
anova(mod.sar.pacI1., mod.sar.pacI1.15) # 0.22944
rm(mod.sar.pacI1.15)

# significant at the 0.05 level
mod.sar.pacI1.17 <- update(mod.sar.pacI1., ~. + depth10deg:meanSal.0m)
summary(mod.sar.pacI1.17)
anova(mod.sar.pacI1., mod.sar.pacI1.17) # 0.044703
rm(mod.sar.pacI1.17)

mod.sar.pacI1.18 <- update(mod.sar.pacI1., ~. + depth10deg:carb_ion )
summary(mod.sar.pacI1.18)
anova(mod.sar.pacI1., mod.sar.pacI1.18) # 0.14779
rm(mod.sar.pacI1.18)

mod.sar.pacI1.19 <- update(mod.sar.pacI1., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI1.19)
anova(mod.sar.pacI1., mod.sar.pacI1.19) # 0.90735
rm(mod.sar.pacI1.19)

mod.sar.pacI1.20 <- update(mod.sar.pacI1., ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI1.20)
anova(mod.sar.pacI1., mod.sar.pacI1.20) # 0.26247
rm(mod.sar.pacI1.20)

# add in poly(meanSST.1deg, 3):mean.mld.t and recalculate
mod.sar.pacI2. <- mod.sar.pacI1.3
rm(mod.sar.pacI1.3)

# try adding in the other interactions
mod.sar.pacI2.2 <- update(mod.sar.pacI2., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI2.2)
anova(mod.sar.pacI2., mod.sar.pacI2.2) # 0.59391
rm(mod.sar.pacI2.2)

# significant at the 0.05 level
mod.sar.pacI2.4 <- update(mod.sar.pacI2., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI2.4)
anova(mod.sar.pacI2., mod.sar.pacI2.4) # 0.027863
rm(mod.sar.pacI2.4)

# significant at the 0.05 level
mod.sar.pacI2.5 <- update(mod.sar.pacI2., ~. - poly(meanSST.1deg, 2):logProd.mn.ann + poly(meanSST.1deg, 3):logProd.mn.ann)
summary(mod.sar.pacI2.5)
anova(mod.sar.pacI2., mod.sar.pacI2.5) # 0.0049051

mod.sar.pacI2.6 <- update(mod.sar.pacI2., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI2.6)
anova(mod.sar.pacI2., mod.sar.pacI2.6) # 0.7395
rm(mod.sar.pacI2.6)

mod.sar.pacI2.7 <- update(mod.sar.pacI2., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI2.7)
anova(mod.sar.pacI2., mod.sar.pacI2.7) # 0.67182
rm(mod.sar.pacI2.7)

mod.sar.pacI2.8 <- update(mod.sar.pacI2., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI2.8)
anova(mod.sar.pacI2., mod.sar.pacI2.8) # 0.53985
rm(mod.sar.pacI2.8)

# significant at the 0.05 level
mod.sar.pacI2.9 <- update(mod.sar.pacI2., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.pacI2.9)
anova(mod.sar.pacI2., mod.sar.pacI2.9) # 0.032021
rm(mod.sar.pacI2.9)

mod.sar.pacI2.10 <- update(mod.sar.pacI2., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI2.10)
anova(mod.sar.pacI2., mod.sar.pacI2.10) # 0.45616
rm(mod.sar.pacI2.10)

mod.sar.pacI2.11 <- update(mod.sar.pacI2., ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI2.11)
anova(mod.sar.pacI2., mod.sar.pacI2.11) # 0.98345
rm(mod.sar.pacI2.11)

mod.sar.pacI2.12 <- update(mod.sar.pacI2., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI2.12)
anova(mod.sar.pacI2., mod.sar.pacI2.12) # 0.32097
rm(mod.sar.pacI2.12)

mod.sar.pacI2.13 <- update(mod.sar.pacI2., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI2.13)
anova(mod.sar.pacI2., mod.sar.pacI2.13) # 0.49106
rm(mod.sar.pacI2.13)

mod.sar.pacI2.14 <- update(mod.sar.pacI2., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI2.14)
anova(mod.sar.pacI2., mod.sar.pacI2.14) # 0.9384
rm(mod.sar.pacI2.14)

mod.sar.pacI2.15 <- update(mod.sar.pacI2., ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI2.15)
anova(mod.sar.pacI2., mod.sar.pacI2.15) # 0.41799
rm(mod.sar.pacI2.15)

# significant at the 0.05 level
mod.sar.pacI2.17 <- update(mod.sar.pacI2., ~. + depth10deg:meanSal.0m)
summary(mod.sar.pacI2.17)
anova(mod.sar.pacI2., mod.sar.pacI2.17) # 0.044389
rm(mod.sar.pacI2.17)

mod.sar.pacI2.18 <- update(mod.sar.pacI2., ~. + depth10deg:carb_ion )
summary(mod.sar.pacI2.18)
anova(mod.sar.pacI2., mod.sar.pacI2.18) # 0.18737
rm(mod.sar.pacI2.18)

mod.sar.pacI2.19 <- update(mod.sar.pacI2., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI2.19)
anova(mod.sar.pacI2., mod.sar.pacI2.19) # 0.84719
rm(mod.sar.pacI2.19)

mod.sar.pacI2.20 <- update(mod.sar.pacI2., ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI2.20)
anova(mod.sar.pacI2., mod.sar.pacI2.20) # 0.28964
rm(mod.sar.pacI2.20)


# add in poly(meanSST.1deg, 3):logProd.mn.ann and recalculate
mod.sar.pacI3. <- mod.sar.pacI2.5
rm(mod.sar.pacI2.5)

# try adding in the other interactions
mod.sar.pacI3.2 <- update(mod.sar.pacI3., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI3.2)
anova(mod.sar.pacI3., mod.sar.pacI3.2) # 0.90999
rm(mod.sar.pacI3.2)

# significant at the 0.05 level
mod.sar.pacI3.4 <- update(mod.sar.pacI3., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI3.4)
anova(mod.sar.pacI3., mod.sar.pacI3.4) # 0.046959
rm(mod.sar.pacI3.4)

mod.sar.pacI3.6 <- update(mod.sar.pacI3., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI3.6)
anova(mod.sar.pacI3., mod.sar.pacI3.6) # 0.79265
rm(mod.sar.pacI3.6)

mod.sar.pacI3.7 <- update(mod.sar.pacI3., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI3.7)
anova(mod.sar.pacI3., mod.sar.pacI3.7) # 0.89132
rm(mod.sar.pacI3.7)

mod.sar.pacI3.8 <- update(mod.sar.pacI3., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI3.8)
anova(mod.sar.pacI3., mod.sar.pacI3.8) # 0.54568
rm(mod.sar.pacI3.8)

# significant at the 0.05 level
mod.sar.pacI3.9 <- update(mod.sar.pacI3., ~. + sdSST.1deg:logProd.mn.ann )
summary(mod.sar.pacI3.9)
anova(mod.sar.pacI3., mod.sar.pacI3.9) # 0.00076959

mod.sar.pacI3.10 <- update(mod.sar.pacI3., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI3.10)
anova(mod.sar.pacI3., mod.sar.pacI3.10) # 0.70217
rm(mod.sar.pacI3.10)

mod.sar.pacI3.11 <- update(mod.sar.pacI3., ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI3.11)
anova(mod.sar.pacI3., mod.sar.pacI3.11) # 0.65807
rm(mod.sar.pacI3.11)

mod.sar.pacI3.12 <- update(mod.sar.pacI3., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI3.12)
anova(mod.sar.pacI3., mod.sar.pacI3.12) # 0.17522
rm(mod.sar.pacI3.12)

mod.sar.pacI3.13 <- update(mod.sar.pacI3., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI3.13)
anova(mod.sar.pacI3., mod.sar.pacI3.13) # 0.55169
rm(mod.sar.pacI3.13)

mod.sar.pacI3.14 <- update(mod.sar.pacI3., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI3.14)
anova(mod.sar.pacI3., mod.sar.pacI3.14) # 0.82764
rm(mod.sar.pacI3.14)

mod.sar.pacI3.15 <- update(mod.sar.pacI3., ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI3.15)
anova(mod.sar.pacI3., mod.sar.pacI3.15) # 0.66249
rm(mod.sar.pacI3.15)

# significant at the 0.05 level
mod.sar.pacI3.17 <- update(mod.sar.pacI3., ~. + depth10deg:meanSal.0m)
summary(mod.sar.pacI3.17)
anova(mod.sar.pacI3., mod.sar.pacI3.17) # 0.013773
rm(mod.sar.pacI3.17)

mod.sar.pacI3.18 <- update(mod.sar.pacI3., ~. + depth10deg:carb_ion )
summary(mod.sar.pacI3.18)
anova(mod.sar.pacI3., mod.sar.pacI3.18) # 0.22181
rm(mod.sar.pacI3.18)

mod.sar.pacI3.19 <- update(mod.sar.pacI3., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI3.19)
anova(mod.sar.pacI3., mod.sar.pacI3.19) # 0.68451
rm(mod.sar.pacI3.19)

mod.sar.pacI3.20 <- update(mod.sar.pacI3., ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI3.20)
anova(mod.sar.pacI3., mod.sar.pacI3.20) # 0.30835
rm(mod.sar.pacI3.20)


# add in sdSST.1deg:logProd.mn.ann and recalculate
mod.sar.pacI4. <- mod.sar.pacI3.9
rm(mod.sar.pacI3.9)

# try adding in the other interactions
mod.sar.pacI4.2 <- update(mod.sar.pacI4., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI4.2)
anova(mod.sar.pacI4., mod.sar.pacI4.2) # 0.57088
rm(mod.sar.pacI4.2)

mod.sar.pacI4.4 <- update(mod.sar.pacI4., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI4.4)
anova(mod.sar.pacI4., mod.sar.pacI4.4) # 0.071393
rm(mod.sar.pacI4.4)

mod.sar.pacI4.6 <- update(mod.sar.pacI4., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI4.6)
anova(mod.sar.pacI4., mod.sar.pacI4.6) # 0.85021
rm(mod.sar.pacI4.6)

mod.sar.pacI4.7 <- update(mod.sar.pacI4., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI4.7)
anova(mod.sar.pacI4., mod.sar.pacI4.7) # 0.55951
rm(mod.sar.pacI4.7)

mod.sar.pacI4.8 <- update(mod.sar.pacI4., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI4.8)
anova(mod.sar.pacI4., mod.sar.pacI4.8) # 0.32664
rm(mod.sar.pacI4.8)

mod.sar.pacI4.10 <- update(mod.sar.pacI4., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI4.10)
anova(mod.sar.pacI4., mod.sar.pacI4.10) # 0.94065
rm(mod.sar.pacI4.10)

mod.sar.pacI4.11 <- update(mod.sar.pacI4., ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI4.11)
anova(mod.sar.pacI4., mod.sar.pacI4.11) # 0.71158
rm(mod.sar.pacI4.11)

mod.sar.pacI4.12 <- update(mod.sar.pacI4., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI4.12)
anova(mod.sar.pacI4., mod.sar.pacI4.12) # 0.42457
rm(mod.sar.pacI4.12)

mod.sar.pacI4.13 <- update(mod.sar.pacI4., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI4.13)
anova(mod.sar.pacI4., mod.sar.pacI4.13) # 0.76326
rm(mod.sar.pacI4.13)

mod.sar.pacI4.14 <- update(mod.sar.pacI4., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI4.14)
anova(mod.sar.pacI4., mod.sar.pacI4.14) # 0.52006
rm(mod.sar.pacI4.14)

mod.sar.pacI4.15 <- update(mod.sar.pacI4., ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI4.15)
anova(mod.sar.pacI4., mod.sar.pacI4.15) # 0.80309
rm(mod.sar.pacI4.15)

# significant at the 0.05 level
mod.sar.pacI4.17 <- update(mod.sar.pacI4., ~. + depth10deg:meanSal.0m)
summary(mod.sar.pacI4.17)
anova(mod.sar.pacI4., mod.sar.pacI4.17) # 0.024867

mod.sar.pacI4.18 <- update(mod.sar.pacI4., ~. + depth10deg:carb_ion )
summary(mod.sar.pacI4.18)
anova(mod.sar.pacI4., mod.sar.pacI4.18) # 0.33344
rm(mod.sar.pacI4.18)

mod.sar.pacI4.19 <- update(mod.sar.pacI4., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI4.19)
anova(mod.sar.pacI4., mod.sar.pacI4.19) # 0.56724
rm(mod.sar.pacI4.19)

mod.sar.pacI4.20 <- update(mod.sar.pacI4., ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI4.20)
anova(mod.sar.pacI4., mod.sar.pacI4.20) # 0.70814
rm(mod.sar.pacI4.20)


# add in depth10deg:meanSal.0m and recalculate
mod.sar.pacI5. <- mod.sar.pacI4.17
rm(mod.sar.pacI4.17)

# try adding in the other interactions
mod.sar.pacI5.2 <- update(mod.sar.pacI5., ~. + poly(meanSST.1deg, 3):sdSST.1deg)
summary(mod.sar.pacI5.2)
anova(mod.sar.pacI5., mod.sar.pacI5.2) # 0.33918
rm(mod.sar.pacI5.2)

mod.sar.pacI5.4 <- update(mod.sar.pacI5., ~. + poly(meanSST.1deg, 3):depth10deg)
summary(mod.sar.pacI5.4)
anova(mod.sar.pacI5., mod.sar.pacI5.4) # 0.39369
rm(mod.sar.pacI5.4)

mod.sar.pacI5.6 <- update(mod.sar.pacI5., ~. -poly(meanSST.1deg, 2):meanSal.0m + poly(meanSST.1deg, 3):meanSal.0m)
summary(mod.sar.pacI5.6)
anova(mod.sar.pacI5., mod.sar.pacI5.6) # 0.97758
rm(mod.sar.pacI5.6)

mod.sar.pacI5.7 <- update(mod.sar.pacI5., ~. +poly(meanSST.1deg, 3):sdSal.0m )
summary(mod.sar.pacI5.7)
anova(mod.sar.pacI5., mod.sar.pacI5.7) # 0.62082
rm(mod.sar.pacI5.7)

mod.sar.pacI5.8 <- update(mod.sar.pacI5., ~. - carb_ion:poly(meanSST.1deg, 2) + carb_ion:poly(meanSST.1deg, 3) )
summary(mod.sar.pacI5.8)
anova(mod.sar.pacI5., mod.sar.pacI5.8) # 0.32956
rm(mod.sar.pacI5.8)

mod.sar.pacI5.10 <- update(mod.sar.pacI5., ~. + sdSST.1deg:carb_ion)
summary(mod.sar.pacI5.10)
anova(mod.sar.pacI5., mod.sar.pacI5.10) # 0.83809
rm(mod.sar.pacI5.10)

mod.sar.pacI5.11 <- update(mod.sar.pacI5., ~. + mean.mld.t:depth10deg)
summary(mod.sar.pacI5.11)
anova(mod.sar.pacI5., mod.sar.pacI5.11) # 0.38446
rm(mod.sar.pacI5.11)

mod.sar.pacI5.12 <- update(mod.sar.pacI5., ~. + mean.mld.t:logProd.mn.ann)
summary(mod.sar.pacI5.12)
anova(mod.sar.pacI5., mod.sar.pacI5.12) # 0.57214
rm(mod.sar.pacI5.12)

mod.sar.pacI5.13 <- update(mod.sar.pacI5., ~. + mean.mld.t:meanSal.0m )
summary(mod.sar.pacI5.13)
anova(mod.sar.pacI5., mod.sar.pacI5.13) # 0.71183
rm(mod.sar.pacI5.13)

mod.sar.pacI5.14 <- update(mod.sar.pacI5., ~. + mean.mld.t:sdSal.0m )
summary(mod.sar.pacI5.14)
anova(mod.sar.pacI5., mod.sar.pacI5.14) # 0.52065
rm(mod.sar.pacI5.14)

mod.sar.pacI5.15 <- update(mod.sar.pacI5., ~. + mean.mld.t:carb_ion )
summary(mod.sar.pacI5.15)
anova(mod.sar.pacI5., mod.sar.pacI5.15) # 0.80717
rm(mod.sar.pacI5.15)

mod.sar.pacI5.18 <- update(mod.sar.pacI5., ~. + depth10deg:carb_ion )
summary(mod.sar.pacI5.18)
anova(mod.sar.pacI5., mod.sar.pacI5.18) # 0.53719
rm(mod.sar.pacI5.18)

mod.sar.pacI5.19 <- update(mod.sar.pacI5., ~. + meanSal.0m:sdSal.0m)
summary(mod.sar.pacI5.19)
anova(mod.sar.pacI5., mod.sar.pacI5.19) # 0.61537
rm(mod.sar.pacI5.19)

mod.sar.pacI5.20 <- update(mod.sar.pacI5., ~. + sdSal.0m:carb_ion)
summary(mod.sar.pacI5.20)
anova(mod.sar.pacI5., mod.sar.pacI5.20) # 0.76832
rm(mod.sar.pacI5.20)

# don't need to add anything else
mod.sar.pacIf <- mod.sar.pacI5.

summary(mod.sar.pacIf, Nagelkerke = TRUE) # 0.60717 
AIC(mod.sar.pacIf) # 3204.398

# Check for correlation
env.var.pac <- c("meanSST.1deg", "sdSST.1deg", "depth10deg", "mean.mld.t", "logProd.mn.ann", "meanSal.0m", "sdSal.0m", "carb_ion")

# pairs plot
png("Figures/Ana_3vii_pac_pairs.png", 1200, 1200)
pairs(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])
dev.off()

# variance inflation factor
vif(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])
cor(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", names(ldg.margo.mod) %in% env.var.pac])

lr.sar.pacIf <- lr.calc(mod.sar.pacIf)

png("Figures/Ana_3vii_LRatio_pacIf.png", width = 800)
lr.plot(lr.sar.pacIf, order = c(8:10, 2:3, 1, 6, 4:5, 7), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 5)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.pacIf), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.pacIfg <- lr.calc(mod.sar.pacIf, tmp)
rm(tmp)

png("Figures/Ana_3ix_LRatio_g_pacIf.png", width = 800)
lr.plot(lr.sar.pacIfg, order = c(6, 2, 1, 3:5), legend = FALSE, ylab = "Log Likelihood ratio", star.pos = 7)
dev.off()

save(mod.sar.pacI, mod.sar.pacI1., mod.sar.pacI2., mod.sar.pacI3., mod.sar.pacI4., mod.sar.pacI5., file = "Outputs/Pacific_simplification.RData")
rm(env.var.pac, mod.sar.pacI, mod.sar.pacI1., mod.sar.pacI2., mod.sar.pacI3., mod.sar.pacI4., mod.sar.pacI5.)

## 3x. Comparison of the oceans --------------------------------------------
png("Figures/Ana_3x_LRatio_oceIf.png", width = 800)
lr.plot(lr.sar.op0, lr.sar.atlIf, lr.sar.indIf, lr.sar.pacIf, order = c(7:9, 4, 2, 11:10, 3, 5, 1, 6), leg.txt = c("Full", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 5)
dev.off()

# also for groups of variables
png("Figures/Ana_3x_LRatio_g_OceIf.png", width = 8, height = 6, units = 'in', res = 300)
tmp <- lr.plot(lr.sar.op0g, lr.sar.atlIfg, lr.sar.indIfg, lr.sar.pacIfg, order = c(6:7, 5, 3, 4, 1:2), leg.txt = c("All", "Atlantic", "Indian", "Pacific"), ylab = "Log Likelihood ratio", star.pos = 8, srt = 50)
dev.off()

save(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, file = "Outputs/Atlantic_simplified.RData")
save(lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, file = "Outputs/Indian_simplified.RData")
save(lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, file = "Outputs/Pacific_simplified.RData")

rm(lr.sar.atlIf, lr.sar.atlIfg, mod.sar.atlIf, lr.sar.indIf, lr.sar.indIfg, mod.sar.indIf, lr.sar.pacIf, lr.sar.pacIfg, mod.sar.pacIf, op.formula, atl.nb, atl.w, pac.w)


## 4. Does resolution of variables matter? ---------------------------------

## 4i. Run full model for higher resolution ----------------------
mod.hres.op0 <- errorsarlm(rarefy.sr ~ (poly(meanSST.4km, 3) + sdSST.4km + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)

## 4iii. Compare LRs -------------------------------------------------------
lr.hres.op0 <- lr.calc(mod.hres.op0)

png("Figures/Ana_4iii_LRatio_ophresop0.png", width = 800)
# get order from running without order first
lr.plot(lr.sar.op0, lr.hres.op0, order = c(7, 12, 8, 13, 9, 14, 4, 2, 10:11, 15, 3, 5:6, 1), leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.hres.op0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.hres.op0g <- lr.calc(mod.hres.op0, tmp)
rm(tmp)

png("Figures/Ana_4iii_LRatio_g_ophresop0.png", width = 800)
lr.plot(lr.sar.op0g, lr.hres.op0g, order = c(6:7, 5, 3, 4, 2:1), leg.x = 17, leg.y = 275, leg.txt = c("Coarser resolution", "Finer resolution"), ylab = "Log Likelihood ratio")
dev.off()

save(mod.hres.op0, lr.hres.op0, lr.hres.op0g, file = "Outputs/mod_hres.RData")
rm(mod.hres.op0, lr.hres.op0, lr.hres.op0g)


## 5. Does carb_ion cut-off matter? -------------------------------------

## 6. Does averaging explanatory variables matter? -------------------------

## 7. Does using rarefied richness make a difference? ----------------------


## 8. Evenness -------------------------------------------------------------

## 8i. Create an OLS model ------------------------------------------------
mod.eve.l0 <- lm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, data = ldg.margo.mod)

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
mod.sar.eveW <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveW$obj, Nagelkerke = TRUE) # 0.43969
AIC(mod.sar.eveW$obj) # 4595.396

# check other coding styles
mod.sar.eveB <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveB$obj, Nagelkerke = TRUE) # 0.4369
AIC(mod.sar.eveB$obj) # -4584.184

mod.sar.eveS <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.eveS$obj, Nagelkerke = TRUE) # 0.45442
AIC(mod.sar.eveS$obj) # -4655.743

# So "S" is the best coding style and the best neighbourhood distance is 543.5247

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
summary(mod.sar.eve0, Nagelkerke = T) # r2 = 0.45442 
AIC(mod.sar.eve0) # -4655.743

lr.sar.eve0 <- lr.calc(mod.sar.eve0)

png("Figures/Ana_8iv_LRatio_eve0.png", width = 800)
lr.plot(lr.sar.eve0, order = c(9:11, 2:3, 1, 6, 4:5, 7:8), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

## 8v. Calculate likelihood ratios for groups of EVs ---------------------

# for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.eve0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.eve0g <- lr.calc(mod.sar.eve0, tmp)
rm(tmp)

png("Figures/Ana_8iv_LRatio_g_eve0.png", width = 800)
lr.plot(lr.sar.eve0g, order = c(7, 2:1, 3:6), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

# mnt, d10, prod (can't tell if seasonality or productivity), seasonality (sdt, sdsal), salinity, carb_ion, ocean

# full model lr plots for each ocean? map differences of RV between simplified model with and without ocean (by points not by layer) - Supplementary info

## 8vi. Effect of carb_ion ----------------------------------------------


## 9. Lineage age -------------------------------------------------------------

## 9i. Lineage age models -------------------------------------------------
mod.sar.lnaW <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaW$obj, Nagelkerke = TRUE) # 0.58086
AIC(mod.sar.lnaW$obj) # 7506.566

# check other coding styles
mod.sar.lnaB <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "B", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaB$obj, Nagelkerke = TRUE) # 0.52895 
AIC(mod.sar.lnaB$obj) #  7770.955

mod.sar.lnaS <- with(ldg.margo.mod, sar.optimised(mod.eve.l0.sac$real$x.intercept, MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, ldg.coords, style = "S", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(mod.sar.lnaS$obj, Nagelkerke = TRUE) # 0.56298 
AIC(mod.sar.lnaS$obj) # 7601.177

mod.sar.lnaW
# So "W" is the best coding style and the best neighbourhood distance is 537.6619

lna.nb <- dnearneigh(ldg.coords, 0, mod.sar.lnaW$dist, longlat = TRUE)
lna.w <- nb2listw(lna.nb, glist = NULL, style = "W", zero.policy = TRUE)
mod.sar.lna0 <- errorsarlm(mod.sar.lnaS$mod, listw = lna.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(mod.sar.lna0, Nagelkerke = TRUE) # 0.58086
rm(lna.nb, lna.w)

save(mod.sar.lnaB, mod.sar.lnaS, mod.sar.lnaW, file = "Outputs/Lineage_coding.RData")
rm(lna.nb, mod.sar.lnaB, mod.sar.lnaW)

lr.sar.lna0 <- lr.calc(mod.sar.lna0)

png("Figures/Ana_9_LRatio_lna0.png", width = 800)
lr.plot(lr.sar.lna0, order = c(9:11, 2:3, 1, 6, 4:5, 7:8), ylab = "Log Likelihood ratio", star.pos = 10, legend = FALSE)
dev.off()

# also for groups of variables
(tmp <- data.frame(names = model.evs(mod.sar.lna0), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature")))
lr.sar.lna0g <- lr.calc(mod.sar.lna0, tmp)
rm(tmp)

png("Figures/Ana_9_LRatio_g_lna0.png", width = 800)
lr.plot(lr.sar.lna0g, order = c(7, 2:1, 3:6), ylab = "Log Likelihood ratio", star.pos = 7, legend = FALSE)
dev.off()

## 9ii. Dissolution cutoffs ------------------------------------------------


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

png("Figures/Ana_10_LRatio_op0f.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.opfg, order = c(5, 7, 6, 3:4, 2:1), ylab = "Log Likelihood ratio", leg.txt = c("Full", "Simplified"), cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
dev.off()

## 10iv. log likelihood ratio plot comparing different carb_ion cut-of --------
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
# png("Figures/Ana_10_LRatio_g_op468.png", width = 800)
# lr.plot(lr.sar8.op0g, lr.sar.op0g, lr.sar4.op0g, order = c(6, 7, 5, 3, 4, 2, 1), leg.txt = c("8%", "6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2)
# dev.off()
# 
# png("Figures/Ana_10_LRatio_g_op46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.op0g, lr.sar4.op0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, srt = 45, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2)
# dev.off()
# 
# 
# png("Figures/Ana_10_LRatio_g_eve46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.eve0g, lr.sar4.eve0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 10, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
# dev.off()
# 
# png("Figures/Ana_10_LRatio_g_lna46.png", width = 10, height = 6, units = 'in', res = 300)
# lr.plot(lr.sar.lna0g, lr.sar4.lna0g, order = c(5, 7, 6, 3:4, 2:1), leg.txt = c("6%", "4%"), ylab = "Log Likelihood ratio", star.pos = 7, cex.names = 1.2, cex.lab = 1.2, cex.axis = 1.2, leg.cex = 1.2, srt = 45)
# dev.off()

## 10v. residuals map for full model --------------------------------------------
png("Figures/Ana_10_rsp_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.op0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -8, max.col = 8, maintitle = "Residuals for rarefied richness"))
dev.off()

png("Figures/Ana_10_eve_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.eve0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -0.5, max.col = 0.5, maintitle = "Residuals for evenness"))
dev.off()

png("Figures/Ana_10_lna_res.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, mod.sar.lna0$residuals, palette = "rwbt", col.water = "white", col.land = "black", min.col = -6, max.col = 6, maintitle = "Residuals for average community age"))
dev.off()

## 10vi. coplot ------------------------------------------------------------------

# plot marginal effects (hold all else constant) of coefficients of T polynomial in oceans. Plot across temperature range observed in the oceans (3 colours / plotting symbols)

## 10vii. comparing all the models ------------------------------------------------
png("Figures/Ana_10_LRatio_REM.png", width = 800)
lr.plot(lr.sar.op0, lr.sar.eve0, lr.sar.lna0, order = c(7:9, 4, 1, 11:10, 3, 5, 2, 6), ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age"))
dev.off()

png("Figures/Ana_10_LRatio_g_REM.png", width = 10, height = 6, units = 'in', res = 300)
lr.plot(lr.sar.op0g, lr.sar.eve0g, lr.sar.lna0g, order = c(6:7, 5, 3:4, 2:1),  ylab = "Log Likelihood ratio", star.pos = 10, leg.txt = c("Rarefied species richness", "Evenness", "Average community age"), srt = 45)
dev.off()

## 10viii. compare this to a random model ----------------------------------------
# rdm.vals <- data.frame(model = rep(NA, 1000 * 7), names = NA, lr = NA, p = NA, stars = NA)
# tmp <- data.frame(names = model.evs(mod.rdm), group = c("Stability", "Vertical niche structure", "Vertical niche structure", "Productivity", "Salinity", "Stability", "Ocean", "Dissolution", "Temperature", "Temperature", "Temperature"))
# 
# for(i in 1:1000) {
#   print(i)
#   mod.rdm <- errorsarlm(rnorm(nrow(ldg.margo.mod)) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.mld.t + depth10deg + logProd.mn.ann + meanSal.0m + sdSal.0m + Ocean2 + carb_ion)^2, listw = op.w, zero.policy = TRUE, tol.solve = 1e-18, data = ldg.margo.mod)
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
  png(paste("Figures/map_", i, ".png", sep = ""), width = 800, height = 450)
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
png("Figures/Ana_11_rsr_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$rarefy.sr > 0, ], distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black", min.col = 0, max.col = 27))
dev.off()

png("Figures/Ana_11_rsr_obs.png", 700, 500)
with(ldg.margo.mod, distrib.map(Longitude, Latitude, rarefy.sr, palette = "matlab.like", min.col = 0, max.col = 27, col.water = "white", col.land = "black"))
dev.off()

## 11iii. predict evenness for this dataset ---------------------------------------
ldg.p.margo$simpsonEve <- sar.predict(mod.sar.eve0, newdata = ldg.p.margo, olddata = ldg.margo.mod)
summary(ldg.p.margo$simpsonEve)

with(ldg.p.margo, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4)) # get some sites with evenness outside sensible limits - where are these
with(ldg.p.margo[ldg.p.margo$simpsonEve > 1, ], distrib.map(Longitude, Latitude, simpsonEve, pch = 15, cex = 0.4)) # around arctic / antarctica coastlines - not a problem

# compare with observed
png("Figures/Ana_11_eve_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$simpsonEve <= 1, ], distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_11_eve_obs.png", 700, 500)
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
png("Figures/Ana_11_lna_pred.png", 700, 500)
with(ldg.p.margo[ldg.p.margo$MorphoAgeAbun >= 5 & ldg.p.margo$MorphoAgeAbun <= 17, ], distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", pch = 15, cex = 0.4, col.water = "white", col.land = "black"))
dev.off()

png("Figures/Ana_11_lna_obs.png", 700, 500)
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
# generate the mean values for the data
data.names <- names(ldg.p.margo)
data.names <- data.names[-c(1:3, 22:27)]
data.means <- colMeans(ldg.margo.mod[, data.names], na.rm = TRUE)

# create a dataframe from this
pred.data <- data.frame(Ocean2 = rep(c("Atlantic", "Indian", "Pacific"), each = 100))
for (i in 1:length(data.names)) {
  pred.data <- cbind(pred.data, rep(data.means[i], 300))
  names(pred.data)[i + 1] <- data.names[i]
}
rm(i)
pred.data$carb_ion <- 0
head(pred.data)

data.xlab <- c(expression(paste("Mean SST / ", degree, "C")), expression(paste("sd SST / ", degree, "C")), "Mixed layer depth / m", expression(paste("10", degree, "C Contour Depth / m")), expression(paste("Mean log (productivity) / log(mg /", m^3, ")")), expression(paste("sd log (Productivity) / log(mg / ", m^3, ")")), "Mean Surface Salinity", "sd Surface Salinity", "Dissolution")
names(data.xlab) <- names(pred.data)[-1]

# create a function for plots of each ocean
oce.plt <- function(p.data, file.name, model, xlab, ylab, data.reset, main) {
  # create plots for each ocean
  for (i in names(p.data)[-1]) {
    p.data[p.data$Ocean2 == "Atlantic", i] <- seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), length.out = 100)
    p.data[p.data$Ocean2 == "Indian", i] <- seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), length.out = 100)
    p.data[p.data$Ocean2 == "Pacific", i] <- seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), length.out = 100)
    pred.atl <- sar.predict(model, newdata = p.data, olddata = ldg.margo.mod)[1:100]
    pred.ind <- sar.predict(model, newdata = p.data, olddata = ldg.margo.mod)[101:200]
    pred.pac <- sar.predict(model, newdata = p.data, olddata = ldg.margo.mod)[201:300]
    png(paste(file.name, i, ".png"), width = 300, height = 300)
    plot(p.data[, i], c(pred.atl, pred.ind, pred.pac), type = "n", ylab = ylab, xlab = xlab[i], bty = "l", main = main)
    points(p.data[p.data$Ocean2 == "Atlantic", i], pred.atl, type = "l", lty = 1)
    points(p.data[p.data$Ocean2 == "Indian", i], pred.ind, type = "l", lty = 2)
    points(p.data[p.data$Ocean2 == "Pacific", i], pred.pac, type = "l", lty = 3)
    #legend("topright", c("Atlantic", "Indian", "Pacific"), lty = 1:3)
    dev.off()
    ifelse (i == "carb_ion", p.data[, i] <- 0, p.data[, i] <- data.reset[i])
  }
}

# run plots for the means of the three models
oce.plt(pred.data, "Figures/Ana_11vi_rsr_", mod.sar.op0, data.xlab, "Rarefied species richness", data.means, "Mean")
oce.plt(pred.data, "Figures/Ana_11vi_eve_", mod.sar.eve0, data.xlab, "Simpsons evenness", data.means, "Mean")
oce.plt(pred.data, "Figures/Ana_11vi_lna_", mod.sar.lna0, data.xlab, "Average community age", data.means, "Mean")

# check the 3d plots to see if the patterns are the same across
# 3d plots
pred.data3d <- data.frame(Ocean2 = rep(rep(c("Atlantic", "Indian", "Pacific"), each = 100), 100))
for (i in 1:length(data.names)) {
  pred.data3d <- cbind(pred.data3d, rep(data.means[i], 30000))
  names(pred.data3d)[i + 1] <- data.names[i]
}
rm(i)
pred.data3d$carb_ion <- 0

i <- names(pred.data3d)[-1][1]
pred.data3d[pred.data3d$Ocean2 == "Atlantic",  i]<- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), length.out = 100), 100)
pred.data3d[pred.data3d$Ocean2 == "Indian", i] <- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), length.out = 100), 100)
pred.data3d[pred.data3d$Ocean2 == "Pacific", i] <- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), length.out = 100), 100)

i <- names(pred.data3d)[-1][2]
pred.data3d[pred.data3d$Ocean2 == "Atlantic",  i] <- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Atlantic", i], na.rm = TRUE), length.out = 100), each = 100)
pred.data3d[pred.data3d$Ocean2 == "Indian", i] <- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), length.out = 100), each = 100)
pred.data3d[pred.data3d$Ocean2 == "Pacific", i] <- rep(seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Pacific", i], na.rm = TRUE), length.out = 100), each = 100)

pred.data3d.atl <- pred.data3d[pred.data3d$Ocean2 == "Atlantic", ]

pred.val3d <- sar.predict(mod.sar.op0, newdata = pred.data3d, olddata = ldg.margo.mod)[which(pred.data3d$Ocean2 == "Atlantic")]
pred.val3d <- matrix(pred.val3d, nrow = 100)

persp(pred.data3d.atl$meanSST.1deg[1:100], seq(min(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), max(ldg.margo.mod[ldg.margo.mod$Ocean2 == "Indian", i], na.rm = TRUE), length.out = 100), pred.val3d)
# suggests they aren't so also run plots for max and min values

# values for min
data.min <- apply(ldg.margo.mod[, data.names], 2, min, na.rm = TRUE)
pred.data.min <- data.frame(Ocean2 = rep(c("Atlantic", "Indian", "Pacific"), each = 100))
for (i in 1:length(data.names)) {
  pred.data.min <- cbind(pred.data.min, rep(data.min[i], 300))
  names(pred.data.min)[i + 1] <- data.names[i]
}
rm(i)
pred.data.min$carb_ion <- 0
head(pred.data.min)

oce.plt(pred.data.min, "Figures/Ana_11vi_rsr_min_", mod.sar.op0, data.xlab, "Rarefied species richness", data.min, "Minimum")
oce.plt(pred.data.min, "Figures/Ana_11vi_eve_min_", mod.sar.eve0, data.xlab, "Simpsons evenness", data.min, "Minimum")
oce.plt(pred.data.min, "Figures/Ana_11vi_lna_min_", mod.sar.lna0, data.xlab, "Average community age", data.min, "Minimum")

# values for max
data.max <- apply(ldg.margo.mod[, data.names], 2, max, na.rm = TRUE)
pred.data.max <- data.frame(Ocean2 = rep(c("Atlantic", "Indian", "Pacific"), each = 100))
for (i in 1:length(data.names)) {
  pred.data.max <- cbind(pred.data.max, rep(data.max[i], 300))
  names(pred.data.max)[i + 1] <- data.names[i]
}
rm(i)
pred.data.max$carb_ion <- 0
head(pred.data.max)

oce.plt(pred.data.max, "Figures/Ana_11vi_rsr_max_", mod.sar.op0, data.xlab, "Rarefied species richness", data.max, "Maximum")
oce.plt(pred.data.max, "Figures/Ana_11vi_eve_max_", mod.sar.eve0, data.xlab, "Simpsons evenness", data.max, "Maximum")
oce.plt(pred.data.max, "Figures/Ana_11vi_lna_max_", mod.sar.lna0, data.xlab, "Average community age", data.max, "Maximum")


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
png("Figures/mte_plot.png")
plot(Ln_SR ~ MTE_SST, pch = 16, col = ldg.margo.mod$Ocean2, xlab = "Temperature (1 / kT)", ylab = "ln (rarefied species richness)", bty = "l", las = 1, cex.lab = 1.2, cex.axis = 1.2)
points(p.MTE_SST, p.Ln_SR, type = "l", lwd = 3) # observed
with(p.ocean[p.ocean$MTE_oce == "Atlantic",], points(MTE_SST, p.Ln_SR, type = "l", col = 1)) # Atlantic
with(p.ocean[p.ocean$MTE_oce == "Indian",], points(MTE_SST, p.Ln_SR, type = "l", col = 2)) # Indian
with(p.ocean[p.ocean$MTE_oce == "Pacific",], points(MTE_SST, p.Ln_SR, type = "l", col = 3)) # Pacific
points(p.MTE_SST, mte.Ln_SR, type = "l", lty = 2, lwd = 3) # predicted
legend("topright", levels(ldg.margo.mod$Ocean2), pch = 16, col = 1:3)
dev.off()

# plot this relationship on the SST ~ SR relationship
png("Figures/sst_mte_plot.png")
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

rm(ldg.coords, stars, env.var, op.w, eve.s, lna.w, mod.sar.eveS, mod.sar.lnaS)

