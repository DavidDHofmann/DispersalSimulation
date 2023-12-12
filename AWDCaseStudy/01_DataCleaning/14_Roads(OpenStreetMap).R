################################################################################
#### Cleaning and Preparation of the Geofabrik Roads Data
################################################################################
# Description: In this script I clean the road files that I downloaded from
# Geofabrik (www.geofabrik.de), which provides ready to download shapefiles
# originating from the open stree maps project. I also cut down the number of
# roads to include only large tar roads.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)    # To handle raster data
library(terra)     # To handle raster data
library(tidyverse) # For data wrangling
library(davidoff)  # Access to custom functions

# Load road data
roads <- vect("03_Data/01_RawData/GEOFABRIK/Roads.shp")

# Check out the description of the classes as derived from the OSM webpage
# (https://wiki.openstreetmap.org/wiki/Key:highway).
legend <- read_csv2("03_Data/01_RawData/GEOFABRIK/RoadsDescription.csv")

# Keep only the largest roads (1-4) and their links (9-12)
roads <- subset(roads, roads$fclass %in% legend$Value[c(1:4, 9:12)])

# Crop the shapefile to our study extent
s <- vect("03_Data/02_CleanData/00_General_Shapefile.shp")
roads <- crop(roads, ext(s))

# Save the cropped shapefile
writeVector(
    x         = roads
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp"
  , overwrite = T
)

################################################################################
#### Rasterize Roads
################################################################################
# Let's load the reference raster so that we can rasterize our roads
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Rasterize the roads
roads$value <- 1
roads_r <- rasterize(roads, r, field = "value", background = 0)

# Store the layer
writeRaster(
    x         = raster(roads_r)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif"
  , overwrite = T
)
