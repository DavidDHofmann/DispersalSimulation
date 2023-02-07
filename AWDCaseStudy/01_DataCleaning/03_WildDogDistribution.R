################################################################################
#### Preparation of the Wild Dog Distribution Data
################################################################################
# Description: Preparation of a shapefile depicting the current spatial
# distribution of wild dogs.

# Clean environment
rm(list = ls())

# Load required packages
library(rgdal)    # To load spatial data

# Load the shapefile
dogs <- readOGR("03_Data/01_RawData/IUCN/data_0.shp")

# Visualize
plot(dogs, col = "purple", border = "black", lwd = 0.5)

# Store the file
writeOGR(
    dogs
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_WildDogs_IUCN"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
