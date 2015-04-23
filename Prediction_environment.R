## Creating a dataset for prediction with values every degree
## Created: 16 / 4 / 15
## Last edited: 23 / 4 / 15
## Isabel Fenton

## Based on the code from Analysis_MARGO.R

## Inputs ------------------------------------------------------------------
# requires ldg.margo.mod

## Outputs -----------------------------------------------------------------
# 150421_OceanBoundaries.RData - points for defining the ocean boundaries

## Libraries ---------------------------------------------------------------
library(fields) # map
data(world.dat) # map
library(sp) # point.in.polygon
source("../../../Code/maps.R") # for maps


## 1. Generate dataframe for predicting ------------------------------------
# meanSST.1deg, sdSST.1deg, mean.mld.t, depth10deg, logProd.mn.ann, meanSal.0m, sdSal.0m, Ocean2, carb_ion
ldg.p.margo <- data.frame(Longitude = rep(-179.5:179.5, 180), Latitude = rep(-89.5:89.5, each = 360))
names(ldg.margo.mod)


## 2. Ocean2 ---------------------------------------------------------------
with(ldg.margo.mod, distrib.map(Longitude, Latitude, Ocean2, key = "FALSE"))

## 2i. Create polygons of the ocean boundaries -----------------------------
# Each ocean is defined using the locator() function. Save these results
# save(pacific, atlantic, atlantic.2, indian, pacific.2, file = "Outputs/150421_OceanBoundaries.RData")
# with these boundaries, then define polygons of the oceans. Exclude the Mediterranean and the Gulf of Mexico, and high northern latitudes where I have no data
pacific.1 <- pacific
pacific.1$x <- c(pacific$x[length(pacific$x)], -180, -180, pacific$x[2:length(pacific$x)])
pacific.1$y <- c(pacific$y[length(pacific$y)], -90, pacific$y[2], pacific$y[2:length(pacific$y)])
points(pacific.1$x, pacific.1$y, type = "l", col = "blue")

atlantic.1 <- atlantic
atlantic.1$x <- c(atlantic$x, pacific$x[9:length(pacific$x)], rev(atlantic.2$x), atlantic$x[1])
atlantic.1$y <- c(atlantic$y, pacific$y[9:length(pacific$x)], rev(atlantic.2$y), atlantic$y[1])
points(atlantic.1$x, atlantic.1$y, type = "l", col = "yellow")

indian.1 <- indian
indian.1$x <- c(indian$x[length(indian$x)], rev(atlantic.2$x)[1:3], indian$x)
indian.1$y <- c(indian$y[length(indian$y)],rev(atlantic.2$y)[1:3], indian$y)
points(indian.1$x, indian.1$y, type = "l", col = "blue")

pacific.3 <- pacific.2
pacific.3$x <- c(180, rev(indian$x)[1:(length(indian$x) - 7)], rev(pacific.2$x), 180)
pacific.3$y <- c(-90, rev(indian$y)[1:(length(indian$y) - 7)], rev(pacific.2$y), -90)
points(pacific.3$x, pacific.3$y, type = "l", col = "green")

rm(pacific, atlantic, atlantic.2, indian, pacific.2)

## 2ii. Add column to ldg.p.margo ------------------------------------------
# calculate the ocean for each point
ldg.p.margo$Ocean2 <- NA
ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, pacific.1$x, pacific.1$y) == 1)] <- "Pacific"
ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, pacific.3$x, pacific.3$y) == 1)] <- "Pacific"
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Pacific", ], points(Longitude, Latitude, col = "yellow"))

ldg.p.margo$Ocean2[which(point.in.polygon(ldg.p.margo$Longitude, ldg.p.margo$Latitude, atlantic.1$x, atlantic.1$y) == 1)] <- "Atlantic"
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
with(ldg.margo.mod, distrib.map(Longitude, Latitude, Ocean2, key = "FALSE"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Atlantic", ], points(Longitude, Latitude, col = "red"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Indian", ], points(Longitude, Latitude, col = "orange"))
with(ldg.p.margo[ldg.p.margo$Ocean2 == "Pacific", ], points(Longitude, Latitude, col = "yellow"))

rm(tmp, land, atlantic.1, indian.1, pacific.1, pacific.3)


## 3. meanSST.4km -------------------------------------------------------------


## 4. sdSST.4km ------------------------------------------------------------


# 5. meanSST.1deg ---------------------------------------------------------
# temperature
load("")
sum(mean.t.depth$Longitude[order(mean.t.depth$Latitude, mean.t.depth$Longitude)] != ldg.p.margo$Longitude)
sum(mean.t.depth$Latitude[order(mean.t.depth$Latitude, mean.t.depth$Longitude)] != ldg.p.margo$Latitude)
ldg.p.margo$meanSST.1deg <- mean.t.depth$depth0m[order(mean.t.depth$Latitude, mean.t.depth$Longitude)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, meanSST.1deg))


# 6. sdSST.1deg --------------------------------------------------------------
sum(sd.t.depth$Longitude[order(sd.t.depth$Latitude, sd.t.depth$Longitude)] != ldg.p.margo$Longitude)
sum(sd.t.depth$Latitude[order(sd.t.depth$Latitude, sd.t.depth$Longitude)] != ldg.p.margo$Latitude)
ldg.p.margo$sdSST.1deg <- sd.t.depth$depth0m[order(sd.t.depth$Latitude, sd.t.depth$Longitude)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, sdSST.1deg))


# 7. mean.mld.t ------------------------------------------------------------


# 8. sd.mld.t -------------------------------------------------------------


# 9. mean.mld.d -----------------------------------------------------------


# 10. sd.mld.d ------------------------------------------------------------


# 11. mean.mld.v ----------------------------------------------------------


# 12. sd.mld.v ------------------------------------------------------------


# 13. depth10deg ----------------------------------------------------------
# depth10deg
head(strat)
sum(strat$Longitude[order(strat$Latitude, strat$Longitude)] != ldg.p.margo$Longitude)
sum(strat$Latitude[order(strat$Latitude, strat$Longitude)] != ldg.p.margo$Latitude)
ldg.p.margo$depth10deg <- strat$depth10deg[order(strat$Latitude, strat$Longitude)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, depth10deg))


# 14. prod.mn.ann ---------------------------------------------------------


# 15. prod.sd.ann ---------------------------------------------------------


# 16. logProd.mn.ann ------------------------------------------------------


# 17. ldgProd.sd.ann ------------------------------------------------------


# 18. meanSal.0m ----------------------------------------------------------
head(sal.mean.depth)
ldg.p.margo$meanSal.0m <- NA
for (i in 1:nrow(sal.mean.depth)) {
  ldg.p.margo$meanSal.0m[sal.mean.depth$Longitude[i] == ldg.p.margo$Longitude & sal.mean.depth$Latitude[i] == ldg.p.margo$Latitude] <- sal.mean.depth$Depth0m[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, meanSal.0m))


# 19. sdSal.0m ------------------------------------------------------------
# salinity

ldg.p.margo$sdSal.0m <- NA
for (i in 1:nrow(sal.sd.depth)) {
  ldg.p.margo$sdSal.0m[sal.sd.depth$Longitude[i] == ldg.p.margo$Longitude & sal.sd.depth$Latitude[i] == ldg.p.margo$Latitude] <- sal.sd.depth$Depth0m[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, sdSal.0m))


# 20. carb_ion ------------------------------------------------------------
# carb_ion
ldg.p.margo$carb_ion <- 0


# 21. Save the data -------------------------------------------------------
save(ldg.p.margo, file = "Outputs/ldg_p_margo.RData")






# mld
ldg.p.margo$mean.mld.t <- NA
for (i in 1:nrow(mld.ann)) {
  ldg.p.margo$mean.mld.t[mld.ann$Longitude[i] == ldg.p.margo$Longitude & mld.ann$Latitude[i] == ldg.p.margo$Latitude] <- mld.ann$mean.mld.t[i]
}
rm(i)
with(ldg.p.margo, distrib.map(Longitude, Latitude, mean.mld.t))


# productivity
head(chla.1deg)
sum(chla.1deg$Longitude[order(chla.1deg$Latitude, chla.1deg$Longitude)] != ldg.p.margo$Longitude)
sum(chla.1deg$Latitude[order(chla.1deg$Latitude, chla.1deg$Longitude)] != ldg.p.margo$Latitude)
ldg.p.margo$logProd.mn.ann <- chla.1deg$mean.logChl[order(chla.1deg$Latitude, chla.1deg$Longitude)]
ldg.p.margo$logProd.sd.ann <- chla.1deg$sd.logChl[order(chla.1deg$Latitude, chla.1deg$Longitude)]
with(ldg.p.margo, distrib.map(Longitude, Latitude, logProd.mn.ann))
with(ldg.p.margo, distrib.map(Longitude, Latitude, logProd.sd.ann))

