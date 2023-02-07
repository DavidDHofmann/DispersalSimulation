################################################################################
#### Preparation of the River Shapefiles (from Dominik)
################################################################################
# Description: Here I cut the rivers that I received from Dominik

# Clear R's brain
rm(list = ls())

# Load packages
library(raster) # For handling spatial data
library(rgdal)  # Fro handling spatial data

# Import shapefile
riv <- readOGR("03_Data/01_RawData/DOMINIK/Rivers.shp")

# Save the file
writeOGR(riv
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_Rivers_DOMINIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
