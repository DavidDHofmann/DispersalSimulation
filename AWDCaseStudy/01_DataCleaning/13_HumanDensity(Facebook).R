################################################################################
#### Cleaning Facebook Human Density Data
################################################################################
# Description: Preparation of the tiles that were downloaded from facebook
# https://data.humdata.org/dataset/highresolutionpopulationdensitymaps
# Includes stitching, aggregating and reprojecting the tiles

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)     # To handle raster data
library(rgdal)      # To handle spatial data

# Make use of multicore abilities
beginCluster()

################################################################################
#### Stitch Tiles
################################################################################
# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/FACEBOOK"
  , pattern    = "population.*tif$"
  , full.names = T
)
dat <- lapply(dat, raster)

# Merge all tiles together
merged <- do.call(merge, dat)

# Load the reference shapefile and raster
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile.shp")
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Crop the merged tiles to our reference shapefile
merged <- crop(merged, s)

# Store the raster to file
writeRaster(
    x         = merged
  , filename  = "03_Data/01_RawData/FACEBOOK/HumanDensity.tif"
  , overwrite = TRUE
)

################################################################################
#### Aggregate and Resample
################################################################################
# Aggregate the layer to 250m
fact <- res(r)[1] / res(merged)[1]
coarse <- aggregate(merged, fact = round(fact), fun = sum)

# Resample the population density layer to the reference raster
coarse <- resample(coarse, raster(r), method = "bilinear")

# Replace NAs and values below 0 with 0s
coarse <- reclassify(coarse, rcl = c(NA, NA, 0))
coarse <- reclassify(coarse, rcl = c(-Inf, 0, 0))

# Store the result
writeRaster(
    x         = coarse
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_FACEBOOK.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
