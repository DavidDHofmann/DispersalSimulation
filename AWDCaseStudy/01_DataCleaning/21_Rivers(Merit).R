################################################################################
#### Preparation of the River Widths By Merit
################################################################################
# Description: Preparation and cleaning of the Merit river layers

# Clear R's brain
rm(list = ls())

# Specify the working directories
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(raster) # To handle raster data
library(terra)  # To handle raster data

################################################################################
#### Stitch Tiles
################################################################################
# Load the data
files <- dir(
    path        = "03_Data/01_RawData/MERIT"
  , pattern     = ".tif$"
  , full.names  = T
)
dat <- lapply(files, raster)

# Merge the tiles
dat <- mosaic(dat[[1]], dat[[2]], dat[[3]], dat[[4]], dat[[5]]
  , dat[[6]], dat[[7]], dat[[8]], dat[[9]], fun = max)

# Convert to terra
dat <- rast(dat)

# Reclassify pixel values and keep only those rivers with a desired width
rcl <- data.frame(from = c(-Inf, 10), to = c(10, Inf), new = c(0, 1))
dat <- classify(dat, rcl)

# Load the reference raster in order to reproject the river layer
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# First, we aggregate the layer to match the resolution of the reference raster
# The reference raster is resolve around 250 meters, the river layer at around
# 90 meters
fact <- res(r)[1] / res(dat)[1]
dat <- aggregate(dat, fact = round(fact), fun = max)

# Secondly, we ca resample the river layer to match the origin and extent of the
# reference raster
dat <- resample(dat, r, "near")

# Store the final result to file
writeRaster(
    raster(dat)
  , "03_Data/02_CleanData/03_LandscapeFeatures_Rivers_MERIT.tif"
  , overwrite = TRUE
)
