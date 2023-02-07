################################################################################
#### Rasterization of Farms by Gabriele
################################################################################
# Description: Here I rasterize all of the farms from the shapefile that
# Gabriele provided

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)     # To handle raster data
library(rgdal)      # To handle vector data
library(gdalUtils)  # Some helpful tools

# Import Gabriele's farms
crops <- readOGR("03_Data/01_RawData/GABRIELE/Farmland.shp")

# Write the layer to the cleaned data
writeOGR(
    crops
  , dsn       = "03_Data/02_CleanData"
  , layer     = "04_AnthropogenicFeatures_Farms_GABRIELE"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)

# Load reference raster for 250 meters
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
values(r) <- 0

# Use gdal to rasterize the farms
writeRaster(
    r
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GABRIELE.tif"
  , overwrite = TRUE
)

# Rasterize the farms
gdal_rasterize(
    src_datasource = "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_GABRIELE.shp"
  , dst_filename   = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GABRIELE.tif"
  , burn           = 1
  , at             = TRUE
)
