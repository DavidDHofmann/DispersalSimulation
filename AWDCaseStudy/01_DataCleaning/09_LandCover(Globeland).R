################################################################################
#### Preparation of Globeland Data
################################################################################
# Description: The script stitches together the seperated tiles downloaded from
# Globeland. After stitching the resulting raster is cropped, simplified
# (reclassified) and resampled to 250m.

# Clear R's brain
rm(list = ls())

# Load packages
library(raster)     # To handle raster data
library(rgdal)      # To handle spatial data
library(gdalUtils)  # To stitch raster tiles
library(tidyverse)  # For data wrangling

# Make use of multiple cores
beginCluster()

################################################################################
#### Stitching the Tiles
################################################################################
# Identify all the tifs that need to be stitched
files <- dir(
    path        = "03_Data/01_RawData/GLOBELAND"
  , pattern     = "030.tif$"
  , full.names  = T
)

# As of now, the rasters come with different projections. We need to adjust this
# and make them equal. Note that by default the algorithm uses "nearest
# neighbor" for resampling
new_files <- c()
for (i in 1:length(files)){
  new_files[i] <- tempfile(fileext = ".tif")
  gdalwarp(
      srcfile = files[i]
    , dstfile = new_files[i]
    , t_srs   = CRS("+init=epsg:4326")
  )
}

# Put the reprojected files into a virtual raster
new_vrt <- tempfile(fileext = ".vrt")
gdalbuildvrt(
    gdalfile    = new_files
  , srcnodata   = 0
  , output.vrt  = new_vrt
)

# Store the virtual raster to file
gdal_translate(
    src_dataset   = new_vrt
  , dst_dataset   = "03_Data/01_RawData/GLOBELAND/Globeland.tif"
  , output_Raster = TRUE
  , options       = c("BIGTIFFS=YES")
)

################################################################################
#### Cropping, Aggregating, and Simplifying the Stitched Raster
################################################################################
# Load the merged file
merged <- raster("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Load the reference shapefile
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile.shp")
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Crop the merged tiles to our reference shapefile
merged <- crop(merged, s)

# Store the merged object to file
writeRaster(
    x         = merged
  , filename  = "03_Data/01_RawData/GLOBELAND/Globeland.tif"
  , overwrite = T
)

# Aggregate to coarser resolution
fact <- res(r)[1] / res(merged)[1]
coarse <- aggregate(merged, fact = round(fact), fun = modal)

# Check out the distribution of values
freq(coarse, useNA = "ifany")

# Prepare classes. Note that I'm going to use the same class description and
# codes as for the copernicus dataset. Note that I'll replace NAs with grassland
info <- data.frame(
    Code = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 255)
  , Class = c(
      "CultivatedLand"
    , "Forest"
    , "Grassland"
    , "Shrubland"
    , "Wetland"
    , "WaterBodies"
    , "Tundra"
    , "ArtificialSurfaces"
    , "Bareland"
    , "PermanentIceSnow"
    , "NA"
  )
  , CodeNew = c(3, 4, 6, 5, 1, 1, 6, 2, 7, 7, 6)
  , ClassNew = c(
      "Cropland"
    , "Forest"
    , "Grassland"
    , "Shrubs"
    , "Water"
    , "Water"
    , "Grassland"
    , "Urban"
    , "Bare"
    , "Bare"
    , "Grassland"
  )
)

# Arrange
info <- arrange(info, CodeNew)

# Assign a color to the new classes
info$Color[info$CodeNew == 1] <- "blue"
info$Color[info$CodeNew == 2] <- "red"
info$Color[info$CodeNew == 3] <- "pink"
info$Color[info$CodeNew == 4] <- "darkgreen"
info$Color[info$CodeNew == 5] <- "orange"
info$Color[info$CodeNew == 6] <- "beige"
info$Color[info$CodeNew == 7] <- "grey"

# Reclassify raster
rcl <- dplyr::select(info, c(Code, CodeNew))
new <- reclassify(coarse, rcl)

# Resample to reference raster
new <- resample(new, r, method = "ngb")

# Check out the frequency of different values (especially NAs)
freq(new, useNA = "ifany")
sum(is.na(values(new)))

# Visualize it
plot(new, col = unique(info$Color), breaks = 0:6)

################################################################################
#### Store Final Raster
################################################################################
# Store the raster
writeRaster(
    x         = new
  , filename  = "03_Data/02_CleanData/01_LandCover_LandCover_GLOBELAND.tif"
  , overwrite = TRUE
)

# Also store the water-cover layer seperately
water <- new == 1
writeRaster(
    x         = water
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_GLOBELAND.tif"
  , overwrite = TRUE
)

# Store the information table
info %>%
  dplyr::select(Class = ClassNew, Code = CodeNew, Color) %>%
  distinct() %>%
  write_csv("03_Data/02_CleanData/01_LandCover_LandCover_GLOBELAND.csv")

# End cluster
endCluster()
