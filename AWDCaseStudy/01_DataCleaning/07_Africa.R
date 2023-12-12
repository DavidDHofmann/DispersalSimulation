################################################################################
#### Preparation of the Africa Shapefile
################################################################################
# Description: Preparation of the Africa shapefile that was downloaded from
# http://www.maplibrary.org/library/stacks/Africa/index.htm

# Clean environment
rm(list = ls())

# load packages
library(rgdal)    # To handle spatial data
library(raster)   # To handle spatial data
library(rgeos)    # To manipulate spatial data

# Load the Africa shapefile
africa <- readOGR("03_Data/01_RawData/ESRI/Africa.shp")

# The file currently misses a correct projection. Let's assign it.
crs(africa) <- CRS("+init=epsg:4326")

# Check if the geometry is valid
gIsValid(africa)

# Make it valid
africa <- gBuffer(africa, width = 0, byid = TRUE)

# Check again
gIsValid(africa)

# Visualize
plot(africa)

# Now save the projected file to our data folder
writeOGR(africa
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Africa_ESRI"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
