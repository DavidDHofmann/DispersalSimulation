################################################################################
#### Averaged Water
################################################################################
# Description: Prepare averaged watermap. This will be a static watermap that we
# can use for our simulations, for plots etc. In contrast to the dynamic
# watermap, we'll prepare it for the entire extent of the KAZA. We'll also use
# this script to prepare a corresponding "DistanceToWater" layer.

# clear r's brain
rm(list = ls())

# Load required packages
library(spatstat)     # To calculate distances quickly
library(maptools)     # To calculate distances quickly
library(tidyverse)    # For data wrangling
library(raster)       # To handle spatial data
library(terra)        # To handle spatial data
library(davidoff)     # Custom functions

# Load and merge static layers
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCover_GLOBELAND.tif")
merit <- rast("03_Data/02_CleanData/03_LandscapeFeatures_Rivers_MERIT.tif")

# Load dyanmic layers
flood <- rast("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# We want to create an "average" floodmap. For this we calculate how often
# each cell was covered by water. If this is more than a desired threshold, we
# will use the cell for our "averaged" map. Let's calculate the threshold
nlayers <- nlyr(flood)
threshold <- nlayers * 0.1

# Sum up all layers to get the number of times each cell was flooded
flood <- sum(flood)

# Every cell that was flooded more than x% of the times will be kept for our
# averaged map. Let's prepare the reclassification table for this operation
rcl <- data.frame(
    oldfrom = c(-Inf, threshold)
  , oldto   = c(threshold, Inf)
  , new     = c(0, 1)
)

# Apply the reclassification
flood <- classify(flood, rcl)

# Combine dynamic with static water layer
flood <- expand(flood, water, value = NA)
water <- cover(flood, water)

# Add merit rivers
water <- max(water, merit)

# Visualize
plot(raster(water), col = c("white", "blue"))

# Store the raster to file
writeRaster(raster(water)
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif"
  , overwrite = TRUE
)

# Use the averaged layer to calculate a distance to water layer
distance <- distanceTo(raster(water), value = 1)

# Plot the resulting layer
plot(distance)

# Store the layer to file
writeRaster(distance
  , filename  = "03_Data/02_CleanData/01_LandCover_DistanceToWaterAveraged_MERGED.tif"
  , overwrite = TRUE
)
