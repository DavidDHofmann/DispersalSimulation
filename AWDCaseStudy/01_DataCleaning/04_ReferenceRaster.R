################################################################################
#### Preparation of the Reference Raster (250m)
################################################################################
# Description: In order to facilitate the manipulation of raster files, I create
# a reference raster according to which I can crop and resample all other
# rasters. This will ensure that all rasters ultimately have the same
# resolution, origin and extent.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)   # To handle raster data
library(rgdal)    # To handle vector data
library(davidoff) # To access custom functions

# Load in the KAZA shapefile and use it as a reference for the spatial extent
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# We want an extent that is slightly larger than the KAZA extent
new_ext <- extent(kaza) + c(-0.5, +0.25, -0.63, +0.4)

# Now we create a raster that is defined by this new extent
r <- raster(new_ext, crs = CRS("+init=epsg:4326"))

# Fill the raster with some random values (nicer to plot)
values(r) <- runif(ncell(r))

# Plot the raster and the kaza together.
plot(r)
plot(kaza, border = "red", add = TRUE)

# Create a 250m raster with the same projection and extent
r250  <- raster(r)
res(r250) <- metersToDegrees(250)

# We also fill it with random data so we can eventually plot the layer. For the
# 30 meters layer this does not work due to the massive amount of pixels.
values(r250) <- runif(ncell(r250))

# Visualize
plot(r250)
plot(kaza, add = T)

# Save the rasters to file
writeRaster(
    x         = r250
  , filename  = "03_Data/02_CleanData/00_General_Raster.tif"
  , overwrite = TRUE
)
