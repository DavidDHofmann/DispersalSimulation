################################################################################
#### Combining all Sources for Water Layers
################################################################################
# Description: In this script I combine the water layers from all different
# sources (Globeland, ORI, OSM, Dominik, MERIT) to create dynamic watermaps.
# Note that the dynamic representation will be limited to the smaller extent of
# the Okavango Delta, whereas the remainder will be based on static watermaps.

# Clean environment
rm(list = ls())

# load packages
library(tidyverse)  # For data wrangling
library(raster)     # To handle spatial data
library(terra)      # To handle spatial data
library(rgdal)      # To handle spatial data
library(pbmcapply)  # To use multiple cores

# Load the layers we want to merge
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(path = ., pattern = ".tif$", full.names = T) %>%
  rast()
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCover_GLOBELAND.tif")
merit <- rast("03_Data/02_CleanData/03_LandscapeFeatures_Rivers_MERIT.tif")

# We only need dynamic watermaps for the extent on which we have GPS data. So
# let's crop with a slight buffer
tracks <- vect("03_Data/02_CleanData/00_General_Dispersers_POPECOL.shp")
extent <- ext(tracks)
extent <- ext(c(xmin(extent)-1, xmax(extent)+1, ymin(extent)-1, ymax(extent)+1))

# Crop watermaps
water <- crop(water, extent)
flood <- crop(flood, extent)
merit <- crop(merit, extent)

# We need to remove the waterbodies for the extent of the dynamic floodmaps.
# Let's get a polygon for the extent for which we have dynamic floodmaps first
p <- ext(trim(flood[[1]]))

# Replace the values below the polygon to 0
water[cells(water, p)] <- 0

# Before we add the dynamic floodmaps, let's merge the globeland and merit data
water <- max(water, merit)

# Coerce to raster
water <- raster(water)
flood <- stack(flood)

# "Expand" floodmaps and fill values with values from static floodmaps
water <- suppressMessages(
  pbmclapply(
      X                  = 1:nlayers(flood)
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x){
      extended <- extend(flood[[x]], water, value = NA)
      filled <- mask(water, extended, maskvalue = 1, updatevalue = 1)
      filled <- writeRaster(filled, tempfile())
      return(filled)
  }) %>% stack()
)

# Let's also transfer the layernames
names(water) <- names(flood)

# Visualize some maps
plot(water[[1:4]], col = c("white", "blue"))

# Save the result to file. We'll store them uncompressed which allows faster
# reading times
writeRaster(
    x         = water
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd"
  , overwrite = TRUE
  , options   = c("COMPRESSION=NONE")
)
