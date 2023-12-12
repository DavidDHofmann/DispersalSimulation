################################################################################
#### Preparation of the River Shapefiles (from Geofrabrik)
################################################################################
# Description: Here I cut the rivers and water areas shapefiles downloaded from
# Geofabrik.

# Clear R's brain
rm(list = ls())

# Load packages
library(raster) # For handling spatial data
library(rgdal)  # Fro handling spatial data

# Import shapefile
riv <- readOGR("03_Data/01_RawData/GEOFABRIK/Rivers.shp")

# We only want to keep bodies that are classified as rivers
riv <- subset(riv, fclass == "river")

# Crop the data according to the reference shapefile
s <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")
riv <- crop(riv, s)

# Save the file
writeOGR(riv
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_Rivers_GEOFABRIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
