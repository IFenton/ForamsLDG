## Created: 11 / 3 / 2015
## Last edited: 2 / 4 / 2015
## Isabel Fenton
## creation of a dataset for the LDG paper from MARGO
## This was split into MARGO/Code/Dataset_setup & MARGO/Code/Diversity_measures, so that the creation and initial analysis of the MARGO dataset was separated from the reanalysis.
##
## Requires running of Environment.RData
##
## this will contain all the columns required for the LDG analysis (and will be complete)
## EVs - SST, mean log chla, sd log chla, mid, width, Ocean, sdSST, 
## RVs - SR, PSV, FRic, FEve, FDiv, simpsonEve
##
## N.B. based on Dataset.R (for LDG analysis from BFD)
## 
## Inputs:
## MARGO.csv - the combined MARGO dataset
## margo_names.csv - the species names found in the MARGO dataset
## 140522_ldg_data.RData - the bfd dataset
## 140523_bfd_tree.Rdata - phylogeny of bfd species for calculating PSV
## 150318_pf_traits.RData - traits of all pfs
##
## Outputs:
## 
## Figures:
## Dat_1v_totalPFs.png - figure showing which sites have a value for the total PF count
## Dat_1vi_BFD.png - figure showing which sites are in BFD (at least approximately)
## Dat_1vi_BFD_wo_percent.png - figure showing which sites are in BFD and have a value of the total PF count
## Dat_1vi_MARGO_cf_BFD.png - how does the BFD data compare to the MARGO data
## Dat_2i_sprich.png - species richness maps with different colour schemes
## Dat_2ii_simpson.png - simpsons diversity maps with different colour schemes
## Dat_2ii_simpsonEve_all.png - simpsons eveness maps with different colour schemes
## Dat_2ii_simpsonEve_0.7.png - simpsons eveness maps with different colour schemes, with only values < 0.7
## Dat_2iii_helmus_psv.png - helmus PSV map
## Dat_2iv_FRic.png - functional richness maps with different colour schemes
## Dat_2iv_FEve.png - functional evenness maps with different colour schemes
## Dat_2iv_FDiv.png - functional divergence maps with different colour schemes
## Dat_2v_MorphoAge.png - average community morphospecies age
## Dat_2v_MorphoAgeAbun.png - average community morphospecies age weighted by abundance
## Dat_2v_LinAge.png - average community lineage age
## Dat_2v_LinAgeAbun.png - average community lineage age weighted by abundance


source("C:/Documents/Science/PhD/Code/compare.R") # the compare function
source("C:/Documents/Science/PhD/Code/maps.R") # for maps
source("C:/Documents/Science/PhD/Code/palettes.R") # ditto
library(fields) # world map 
data(world.dat) # load in world data identifying points on land
library(sp) # point.in.polygon
library(vegan) # simpsons index
library(picante) # PSV
library(FD) # functional diversity measures
library(colorRamps) # matlab colours
setwd("C:/Documents/Science/PhD/Work/1311 LDGPaper/Reanalysis/")

## 1. Load in MARGO and correct names ----------------------------------------

## 1i. Load in the data ----------------------------------------------------
# load in the bfd and the list of names
ldg.margo <- read.csv("C:/Documents/Science/PhD/Project/MARGO/IF_data/MARGO.csv")
ldg.forams.margo <- read.csv("C:/Documents/Science/PhD/Project/MARGO/IF_data/margo_names.csv", header = F)

## 1ii. Correct the names --------------------------------------------------
# run the compare function and check the output
margo.names <- sapply(as.character(ldg.forams.margo[, 1]), compare, USE.NAMES = FALSE)
corr.margo.names <- cbind(margo.names, ldg.forams.margo, 14:58) # third column is position in original margo dataframe
colnames(corr.margo.names) <- c("Corrected.margo", "Original.margo", "position")
corr.margo.names

# create a new data frame that is to be edited
ldg.margo.data <- ldg.margo

# replace the relevant names in the dataset
colnames(ldg.margo.data)[14:58] <- as.vector(corr.margo.names$Corrected.margo)

# Globigerinoides ruber. 16 and 17 are summed in 18. Therefore keep only column 18
# Globigerinoides sacculifer. 20 and 21 are summed in 22. Therefore keep only column 22
# Neogloboquadrina incompta. 34 and 36 are basically the same. But where the don't match the difference should be added to N. dutertrei.
# Truncorotalia truncatulinoides. 41 and 42 are summed in 43. Therefore keep only column 43
# 48, 53, 54, 56 and 58 are micro. Therefore sum
# Globorotalia tumida. 52 is a sum of 49, 50 and 51. Don't want to merge these so keep them separate
# T. humulis. 57 is basically in 31. Therefore remove this column

# deal with the pachyderma-dutertrei intergrade (see MARGO/metadata.docx & Chen et al, 2005)
# identify which are have extra
tmp <- which(ldg.margo.data[, 34] != ldg.margo.data[, 36])
# add these to dutertrei
ldg.margo.data[tmp, c(34, 36)]
ldg.margo.data[tmp, 35] <- ldg.margo.data[tmp, 36] - ldg.margo.data[tmp, 34] + ldg.margo.data[tmp, 35]

rm(tmp)

ldg.margo.data[, 48] <- ldg.margo.data[, 48] + ldg.margo.data[, 53] + ldg.margo.data[, 54] + ldg.margo.data[, 56] + ldg.margo.data[, 58] # merge micro

# remove columns that are duplicated elsewhere
ldg.margo.data <- ldg.margo.data[, -c(16, 17, 20, 21, 36, 41, 42, 52, 53, 54, 56, 57, 58)]

# move Micro to the end
colnames(ldg.margo.data)
ldg.margo.data <- ldg.margo.data[, c(1:40, 42:45, 41, 46:ncol(ldg.margo.data))]

# check they have replaced correctly
head(ldg.margo.data)

rm(ldg.forams.margo, margo.names, corr.margo.names, ldg.margo)

## 1iii. Convert latitude and longitude ------------------------------------
ldg.margo.data$Longitude[which(ldg.margo.data$Longitude > 180)] <- ldg.margo.data$Longitude[which(ldg.margo.data$Longitude > 180)] - 360

## 1iv. Check and remove points on land --------------
ldg.margo.data <- ldg.margo.data[-which(point.in.polygon(ldg.margo.data$Longitude, ldg.margo.data$Latitude, world.dat$x, world.dat$y) != 0), ]

## 1v. Calculate total planktics from percent -----------------------------
# for those points which are percent, but have total PFs, calculate the absolute values

# show which datapoints will be lost
png("Figures/Dataset/Dat_1v_totalPFs.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, Total))
dev.off()

# identify rows
table(ldg.margo.data[, c(49, 51)])
tmp <- which(ldg.margo.data$percent_or_Raw == "P" & ldg.margo.data$Total == "Y")

# identify the columns with species names
margo.macro <- colnames(ldg.margo.data)[14:44]
margo.all.species <- colnames(ldg.margo.data)[14:47]

# divide by current sum & multiply by total
for(i in tmp) 
{
  ldg.margo.data[i, margo.all.species] <- round(ldg.margo.data[i, margo.all.species] / ldg.margo.data$Current_sum[i] * ldg.margo.data$Total_Planktics[i], 0)
}
rm(i, tmp)

## 1vi. Compare with the BFD -----------------------------------------------
# add a column for the BFD
ldg.margo.data$BFD <- "N"
ldg.margo.data$BFD[which(ldg.margo.data$Publication == 'Prell_et_al.,_1999_("Brown_Database")')] <- "Y"
ldg.margo.data$BFD <- factor(ldg.margo.data$BFD)

# plot this up
png("Figures/Dataset/Dat_1vi_BFD.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, BFD))
dev.off()

png("Figures/Dataset/Dat_1vi_BFD_wo_percent.png", 800, 500)
with(ldg.margo.data[ldg.margo.data$Total == "Y",], distrib.map(Longitude, Latitude, BFD))
dev.off()

# check that points from the BFD match those from MARGO
load("C:/Documents/Science/PhD/Work/1311 LDGPaper/Output/140522_ldg_data.RData")

tmp <- which(as.character(ldg.margo.data$Core) == ldg.data$Core.ID[2]) # first row was from a different source
# check source
ldg.margo.data$Publication[tmp]
# compare
merge(ldg.margo.data[tmp,], ldg.data[2, ], all = TRUE) # only difference is whether incompta is separated out

# not all bfd points are in the MARGO dataset. To see which are missing
png("Figures/Dataset/Dat_1vi_MARGO_cf_BFD.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, Total, key = FALSE))
with(ldg.data, points(Long, Lat))
text(0, -97, "red - total PFs, yellow - percentages, black - BFD")
dev.off()

rm(tmp, ldg.data)

## 1vii. Check for duplicates ----------------------------------------------
dim(ldg.margo.data)

# identify exact duplicates
tmp <- which(duplicated(ldg.margo.data))
# check for several that these are duplicates
ldg.margo.data[which(ldg.margo.data$Core == ldg.margo.data$Core[tmp[30]]), ]
# remove duplicated rows
ldg.margo.data <- ldg.margo.data[-tmp, ]
rm(tmp)

# identify duplicates in all the raw data
tmp <- which(duplicated(ldg.margo.data[, c(3:5, 14:47)]))
# check these are duplicates
ldg.margo.data[ldg.margo.data$Core == ldg.margo.data$Core[tmp[5]], ]
ldg.margo.data <- ldg.margo.data[-tmp, ]
rm(tmp)

# same core and lat / long / depth but different data
tmp <- which(duplicated(ldg.margo.data[, c(1, 3:5)]))
write.csv(ldg.margo.data[ldg.margo.data$Core %in% ldg.margo.data$Core[tmp], c(1, 3:5, 11:13)], file = "duplicates.csv")
duplicates <- which(ldg.margo.data$Core %in% ldg.margo.data$Core[tmp])
# to see reasons for choice see "duplicates_edited.csv".
ldg.margo.data <- ldg.margo.data[-duplicates[c(3:4, 6, 9:106, 108:122, 129:134, 137:143, 145:157, 159:161, 209, 275:280, 308, 310:311)], ]

rm(tmp, duplicates)

## 1vii. Remove unnecessary columns ------------------------------------
colnames(ldg.margo.data)

# don't need Coring_device, Sample_depth_upper, Sample_depth_lower, chronozone_level, sedimentation_rate, Publication, added_by, date_of_addition, 
ldg.margo.data <- ldg.margo.data[, -c(2, 7:13)]


## 2. Calculate diversity metrics ------------------------------------------

## 2i. Species Richness ----------------------------------------------------

# calculate the species richness
ldg.margo.data$sp.rich <- rowSums(ldg.margo.data[, which(colnames(ldg.margo.data) %in% margo.macro)] > 0)

png("Figures/Diversity/Dat_2i_sprich.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, sp.rich))
dev.off()

png("Figures/Diversity/Dat_2i_sprich2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, sp.rich, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2i_sprich3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, sp.rich, palette = "matlab.like", col.land = "black", col.water = "white"))
dev.off()

## 2ii. Simpsons index -----------------------------------------------------
# simpson index / simpsons evenness to measure evenness 
# This calculates the sum of the proportions of the square of each species.
# In this context, we assume the infinite version is valid (i.e. removing a specimen doesn't affect the chances of finding another one of the same species)
# 1 is equal abundances of all species, small is very unequal abundances
# 32 species so if equal abundances max would be 0.96875: Actual max = 0.0.9121; min = 0

# calculate simpsons index
ldg.margo.data$simpson <- sapply(1:nrow(ldg.margo.data), function (i) diversity(ldg.margo.data[i, which(colnames(ldg.margo.data) %in% margo.macro)], "simpson"))
# check
1 - sum(ldg.margo.data[2, margo.macro]^2/sum(ldg.margo.data[2, margo.macro])^2)

png("Figures/Diversity/Dat_2ii_simpson.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpson))
dev.off()

png("Figures/Diversity/Dat_2ii_simpson2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpson, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2ii_simpson3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpson, palette = "matlab.like", col.land = "black", col.water = "white"))
dev.off()

# calculate simpsons evenness
ldg.margo.data$simpsonEve <- sapply(1:nrow(ldg.margo.data), function (i) diversity(ldg.margo.data[i, which(colnames(ldg.margo.data) %in% margo.macro)], "invsimpson") / ldg.margo.data$sp.rich[i])
# check
(1 / sum(ldg.margo.data[2, margo.macro]^2/sum(ldg.margo.data[2, margo.macro])^2)) / ldg.margo.data$sp.rich[2]

png("Figures/Diversity/Dat_2ii_simpsonEve_all.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpsonEve))
dev.off()

png("Figures/Diversity/Dat_2ii_simpsonEve_0.7.png", 800, 500)
with(ldg.margo.data[ldg.margo.data$simpsonEve < 0.7, ], distrib.map(Longitude, Latitude, simpsonEve))
dev.off()

png("Figures/Diversity/Dat_2ii_simpsonEve_all2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2ii_simpsonEve_0.7_2.png", 800, 500)
with(ldg.margo.data[ldg.margo.data$simpsonEve < 0.7, ], distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2ii_simpsonEve_all3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", col.land = "black", col.water = "white"))
dev.off()

png("Figures/Diversity/Dat_2ii_simpsonEve_0.7_3.png", 800, 500)
with(ldg.margo.data[ldg.margo.data$simpsonEve < 0.7, ], distrib.map(Longitude, Latitude, simpsonEve, palette = "matlab.like", col.land = "black", col.water = "white"))
dev.off()

## 2iii. Phylogenetic diversity --------------------------------------------
# load in the bfd phylogeny
load("C:/Documents/Science/PhD/Project/Foraminifera/Outputs/140523_bfd_tree.Rdata")

# calculate helmus psv
ldg.margo.data$helmus.psv <- psv(as.matrix(ldg.margo.data[, which(colnames(ldg.margo.data) %in% margo.macro)]), bfd.tree, compute.var=TRUE)$PSV

png("Figures/Diversity/Dat_2iii_helmus_psv.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, helmus.psv))
dev.off()

png("Figures/Diversity/Dat_2iii_helmus_psv2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, helmus.psv, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iii_helmus_psv3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, helmus.psv, palette = "matlab.like", col.land = "black", col.water = "white"))
dev.off()

rm(bfd.tree)

## 2iv. Functional diversity -----------------------------------------------
# generate margo.traits
load("C:/Documents/Science/PhD/Project/Foraminifera/Outputs/150318_pf_traits.RData")
margo.traits <- pf.traits[match(margo.macro, rownames(pf.traits)), ]
# check names
cbind(margo.macro, rownames(margo.traits))

# traits for dbFD - only those that are useful for trait comparisons
db.traits <- margo.traits[, c("spinose", "structure", "mph", "eco", "sqrt.area", "symbionts_R", "depth_R", "dissolution_R")]

# convert characters and integers (i.e. mph / eco) to factors (need to do this to use dbFD)
str(db.traits)
for (i in colnames(db.traits)) 
{
  if (is.character(db.traits[, i]) | is.integer(db.traits[, i])) {
    db.traits[, i] <- factor(db.traits[, i])
  }
}
rm(i)

str(db.traits)

# calculate functional diversity
ldg.FD <- dbFD(db.traits, ldg.margo.data[, margo.macro], corr = "cailliez")

str(ldg.FD)

ldg.margo.data$FRic <- ldg.FD$FRic
ldg.margo.data$FEve <- ldg.FD$FEve
ldg.margo.data$FDiv <- ldg.FD$FDiv

png("Figures/Diversity/Dat_2iv_FRic.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FRic))
dev.off()

png("Figures/Diversity/Dat_2iv_FRic2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FRic, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iv_FRic3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FRic, col.water = "white", col.land = "black", palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iv_FEve.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FEve))
dev.off()

png("Figures/Diversity/Dat_2iv_FEve2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FEve, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iv_FEve3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FEve, col.water = "white", col.land = "black", palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iv_FDiv.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FDiv))
dev.off()

png("Figures/Diversity/Dat_2iv_FDiv2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FDiv, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2iv_FDiv3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, FDiv, col.water = "white", col.land = "black", palette = "matlab.like"))
dev.off()

rm(pf.traits, ldg.FD)

## 2v. Average clade age ---------------------------------------------------
# as MARGO data is from the present day, then we can calculate morphospecies and lineage ages as values relative to recent

# add columns to margo.traits
colnames(margo.traits)

## morphospecies age
# check that all species are extant 
margo.traits$mEn 
rownames(margo.traits)[which(margo.traits$mEn != 0)] # get 2 false's, Truncorotalia crassula and Globorotalia flexuosa. Don't do anything about this

margo.traits$aM.age <- margo.traits$mSt - margo.traits$mEn

## lineage age
# check that all species are extant 
margo.traits$aL.end 
rownames(margo.traits)[which(margo.traits$aL.end != 0)] # Truncorotalia crassula - don't do anything about this

margo.traits$aL.age <- margo.traits$aL.start -  margo.traits$aL.end

# generate a function to calculate average age
ave.age <- function(names, ages, data, abun = FALSE) {
  # names - a vector of the species names
  # ages - a vector of the ages (in the same order as the species names)
  # data - a row of site specific data, with name() the species names
  # abun - if true, calculate abundance weighted age, otherwise just use p/a
  # output - a value of the average age
  
  # which species are greater than 0
  sp <- names(data)[data > 0]
  # what are their ages
  sp.ages <-ages[names %in% sp]
  if(abun) {
    # weight by abundances
    return(sum(sp.ages * data[data > 0] / sum(data)))  
  }
  return(mean(sp.ages))  
}

# example
ave.age(rownames(margo.traits), margo.traits$aM.age, ldg.margo.data[1, colnames(ldg.margo.data) %in% rownames(margo.traits)])

# calculate the average age of each site (both lineage and morphospecies)
ldg.margo.data$MorphoAge <- sapply(1:nrow(ldg.margo.data), function (i) ave.age(rownames(margo.traits), margo.traits$aM.age, ldg.margo.data[i, colnames(ldg.margo.data) %in% rownames(margo.traits)]))

ldg.margo.data$LinAge <- sapply(1:nrow(ldg.margo.data), function (i) ave.age(rownames(margo.traits), margo.traits$aL.age, ldg.margo.data[i, colnames(ldg.margo.data) %in% rownames(margo.traits)]))

# calculate this weighted by evenness
ldg.margo.data$MorphoAgeAbun <- sapply(1:nrow(ldg.margo.data), function (i) ave.age(rownames(margo.traits), margo.traits$aM.age, ldg.margo.data[i, colnames(ldg.margo.data) %in% rownames(margo.traits)], abun = T))

ldg.margo.data$LinAgeAbun <- sapply(1:nrow(ldg.margo.data), function (i) ave.age(rownames(margo.traits), margo.traits$aL.age, ldg.margo.data[i, colnames(ldg.margo.data) %in% rownames(margo.traits)], abun = T))

# try plotting these up
summary(ldg.margo.data$MorphoAge)
summary(ldg.margo.data$LinAge)
summary(ldg.margo.data$MorphoAgeAbun)
summary(ldg.margo.data$LinAgeAbun)

png("Figures/Diversity/Dat_2v_MorphoAge.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAge))
dev.off()

png("Figures/Diversity/Dat_2v_MorphoAge2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAge, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2v_MorphoAge3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAge, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Diversity/Dat_2v_MorphoAgeAbun.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAgeAbun))
dev.off()

png("Figures/Diversity/Dat_2v_MorphoAgeAbun2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2v_MorphoAgeAbun3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, MorphoAgeAbun, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Diversity/Dat_2v_LinAge.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAge))
dev.off()

png("Figures/Diversity/Dat_2v_LinAge2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAge, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2v_LinAge3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAge, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

png("Figures/Diversity/Dat_2v_LinAgeAbun.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAgeAbun))
dev.off()

png("Figures/Diversity/Dat_2v_LinAgeAbun2.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like"))
dev.off()

png("Figures/Diversity/Dat_2v_LinAgeAbun3.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, LinAgeAbun, palette = "matlab.like", col.water = "white", col.land = "black"))
dev.off()

## 2vi. Percent symbionts ---------------------------------------------------
colnames(margo.traits)

## how do I define symbionts?
summary(margo.traits$symbionts_R)
# probably want two percentages, one for faculative and one for obligatory
# assume that ?None is equivalent to None

# generate a function for calculating traits
# worth checking how I did the other stuff for BFD

## calculate the percent of symbiotic species by site
# all symbionts


ldg.margo.data$all_symb <- 
  
  # for each row
  # which species are symbionts (nchar > 5)
  # sum that, or multiply by abundance

# calculate the abundance weighted percent symbionts



## 3. Add environmental variables --------------------------------------------
# Environment.R should be run at this point
load("Output/140522_ldg_env.Rdata")

# check they match (should be FALSE)
which(ldg.margo.data$Latitude != ldg.env$Latitude)


## 4. Produce a dataset without dissolved sites ----------------------------

## 4i. Use the Kucera equation to identify the level of dissolution --------
## Fls (p249, Hillaire-Marcel & de Vernal (2007)) 
head(ldg.margo.data)

# create the function
frag.LS <- function(num.frag, num.comp) {
  (num.frag / 8) / ((num.frag / 8) + num.comp) * 100
}
# n.b. the factor of 8 was introduced as 1 foram doesn't provide one fragment

# play round to discover behaviour
frag.LS(0, 300) # no fragments - 0
frag.LS(300,0) # all fragments - 100
frag.LS(300, 300) # as many fragments as complete - 11.1 or 1/9

# then calculate this for the bfd data
ldg.margo.data$dissolution <- frag.LS(ldg.margo.data$Total.Fragments, ldg.margo.data$Total.Planktics)
# also add to ldg.env
ldg.env$dissolution <- ldg.margo.data$dissolution

## 4ii. Look at the nature of the dissolution values -----------------------
summary(ldg.margo.data$dissolution)

png("Figures/Diversity/Dat_4ii_sortDissol.png")
plot(sort(ldg.margo.data$dissolution))
dev.off()

png("Figures/Diversity/Dat_4ii_dissolution.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, dissolution))
dev.off()

png("Figures/Diversity/Dat_4ii_u15dissolution.png", 800, 500)
with(ldg.margo.data[ldg.margo.data$dissolution < 15, ], distrib.map(Longitude, Latitude, dissolution)) # based on De Villiers (2003) A 425 kyr record of foraminiferal shell weight variability in the western equatorial Pacific
dev.off()

png("Figures/Diversity/Dat_4ii_u6dissolution.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, as.factor(dissolution < 6)))
dev.off()
dim(ldg.margo.data[ldg.margo.data$dissolution < 6, ]) # gives us 965 poinst with a spread across all Oceans

# how does dissolution correlate with depth?
png("Figures/Diversity/Dat_4ii_dissolDepth.png")
with(ldg.margo.data, plot(dissolution, Water.Depth, pch = 16, col = Ocean))
dev.off()

with(ldg.margo.data[ldg.margo.data$Ocean == "Atlantic",], plot(dissolution, Water.Depth, pch = 16, col = Ocean))
# not that clear cut, get high values at greater depths, but a lot of spread


## 5. Rarefaction species richness---------------------------------------------
# calculate rarefied species richness
ldg.margo.data$rarefy.sr <- NA
rare.min <- sort(ldg.margo.data$Total.Planktics)[2]
ldg.margo.data$rarefy.sr[-order(ldg.margo.data$Total.Planktics)[1]] <- rarefy(ldg.margo.data[-order(ldg.margo.data$Total.Planktics)[1], 7:36], rare.min)

summary(ldg.margo.data$rarefy.sr)
# check the one set to NA is the one with the lowest number of planktics (92)
ldg.margo.data[which(is.na(ldg.margo.data$rarefy.sr)), ]

# plot this up and compare with species richness
png("Figures/Diversity/Dat_5_rarefySR.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, rarefy.sr))
dev.off()

png("Figures/Diversity/Dat_5_F_rarefySR.png", 800, 500)
with(ldg.margo.data[!is.na(ldg.margo.data$rarefy.sr), ], distrib.filled(Longitude, Latitude, rarefy.sr))
dev.off()

png("Figures/Diversity/Dat_5_F_SR.png", 800, 500)
with(ldg.margo.data[!is.na(ldg.margo.data$rarefy.sr), ], distrib.filled(Longitude, Latitude, sp.rich))
dev.off()

png("Figures/Diversity/Dat_5_rarefySRmSR.png", 800, 500)
with(ldg.margo.data, distrib.map(Longitude, Latitude, sp.rich - rarefy.sr))
dev.off()

png("Figures/Diversity/Dat_5_F_rarefySRmSR.png", 800, 500)
with(ldg.margo.data[!is.na(ldg.margo.data$rarefy.sr), ], distrib.filled(Longitude, Latitude, sp.rich - rarefy.sr, nlevels = 100))
dev.off()

png("Figures/Diversity/Dat_5_SRLat.png")
with(ldg.margo.data, plot(Latitude, sp.rich, pch = 16, col = Ocean))
dev.off()

png("Figures/Diversity/Dat_5_rarefySRLat.png")
with(ldg.margo.data, plot(Latitude, rarefy.sr, pch = 16, col = Ocean))
dev.off()


## 6. Create a dataframe for modelling -------------------------------------
tmp <- c("Core.ID", "Latitude", "Longitude", "Ocean", "Water.Depth", "dissolution", "sp.rich", "rarefy.sr", "simpson", "simpsonEve", "MorphoAge", "MorphoAgeAbun")
ldg.m.data <- ldg.margo.data[, tmp]

tmp <- c("meanSST.4km", "sdSST.4km", "meanSST.1deg", "sdSST.1deg", "SST.1deg.exact", "mean.pt", "sd.pt", "mld.exact", "depth10deg", "meanChla.4km", "sdChla.4km", "meanChla.1deg", "sdChla.1deg",  "mean.logChla.1deg", "sd.logChla.1deg", "chl.exact", "meanSal.0m", "sdSal.0m", "sal.exact")

ldg.m.data <- cbind(ldg.m.data, ldg.env[, tmp])

rm(tmp)


## 7. Save out the data ----------------------------------------------------
# save complete workspace
save.image(file = "Output/140522_Dataset_ws.Rdata")
# save ldg.margo.data & ldg.m.data
save(ldg.margo.data, file = "Output/140522_ldg_data.Rdata")
save(ldg.m.data, file = "Output/140522_ldg_m_data.Rdata")
