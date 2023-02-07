################################################################################
#### Resampling and Cropping of the Globeland Agriculture Data
################################################################################
# Description: We can use globeland's agricultural fields to augment our
# croplands layer. In this script I therefore extract agricultural fields so we
# can merge them.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)   # To handle raster data

# Make use of multiple cores
beginCluster()

# Import the Globelands dataset
crops <- raster("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Keep only crops
crops <- crops == 10

# Aggregate the globelands dataset to match the resolution of the reference
# raster
fact <- res(r)[1] / res(crops)[1]
crops <- aggregate(crops, fact = round(fact), fun = max)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r, "ngb")

# Check NAs
sum(is.na(values(crops_res)))

# Save the result to file
writeRaster(
    x         = crops_res
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GLOBELAND.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
