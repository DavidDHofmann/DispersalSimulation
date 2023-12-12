################################################################################
#### Preparation of a Reference Shapefile
################################################################################
# Description: Here I create a reference shapefile according to which I will
# crop all the other shapefiles. This is comparable to the reference raster, yet
# for vector data. I will thus use the reference raster to define the spatial
# extent.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)   # To handle raster data
library(rgdal)    # To handle vector data

# I will use the raster created in the previous script for this
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Extract the extent from the raster file and create a polygon with the same
# extent. Note that we also need to reassign the correct crs.
extent <- extent(r)
s <- as(extent, "SpatialPolygons")
crs(s) <- crs(r)

# Plot the shapefile to make sure that it properly frames the reference raster
plot(r)
plot(s, border = "red", add = T, lwd = 2)
area(s) / 1000000

# Assign an ID
s <- SpatialPolygonsDataFrame(s, data = data.frame(Name = "ReferenceShapefile"))

# Save the final object
writeOGR(
    obj       = s
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Shapefile"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
