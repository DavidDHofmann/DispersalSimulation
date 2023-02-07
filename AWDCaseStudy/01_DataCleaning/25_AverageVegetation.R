################################################################################
#### Averaged Vegetation Layer
################################################################################
# Description: Prepare layer that depicts the average vegetation layer. Note
# that this is basically the counterpart of the average water layer

# clear r's brain
rm(list = ls())

# Load required packages
library(tidyverse)  # For data wrangling
library(raster)     # For manipulating spatial data
library(terra)      # For manipulating spatial data

# Load the vegetation layers again
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B/Stitched"
  , pattern     = ".*MODIS.tif$"
  , full.names  = T
)
names <- substr(
    x     = basename(files)
  , start = 14
  , stop  = nchar(basename(files)) - 10
)
modis <- rast(files)
names(modis) <- names

# Make sure values range from 0 to 1
modis <- modis / 100

# Extract the separate layers
modis_shrub <- modis[[1]]
modis_noveg <- modis[[2]]
modis_trees <- modis[[3]]

# Values above 1 are water. Let's reclassify those to 0% Vegetation, i.e. 100%
# NonVegetated
values(modis_shrub)[values(modis_shrub) > 1] <- 0
values(modis_noveg)[values(modis_noveg) > 1] <- 1
values(modis_trees)[values(modis_trees) > 1] <- 0

# Visualize again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Load averaged water layer
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")

# Replace values below water to 0
modis_shrub <- mask(modis_shrub
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 0
)
modis_noveg <- mask(modis_noveg
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 1
)
modis_trees <- mask(modis_trees
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 0
)

# Visualize layers again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Put the stacks into a list
modis <- list(modis_shrub, modis_noveg, modis_trees)

# Prepare filenames
names <- c(
    "03_Data/02_CleanData/01_LandCover_NonTreeVegetationAveraged_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_NonVegetatedAveraged_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_TreeCoverAveraged_MODIS.tif"
)

# Store the rasterstacks
for (i in 1:length(names)){
  writeRaster(raster(modis[[i]]), names[i], overwrite = TRUE)
}
