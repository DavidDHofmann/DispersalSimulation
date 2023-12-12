################################################################################
#### Step Selection Function - Extraction of Covariates
################################################################################
# Description: In this script we will extract all covariates underlying the
# generated steps

# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(lubridate)    # To handle dates
library(amt)          # To coerce gps fixes to steps
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(pbmcapply)    # For multicore abilities
library(spatstat)     # For point patter analysis (calculate distance)
library(maptools)     # For point patter analysis (calculate distance)
library(rgeos)        # For spatial manipulation

# Load the generated steps
lines <- readOGR("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF).shp")

################################################################################
#### LandCover - DistanceToWater
################################################################################
# Transform the lines to utm
lines_utm <- spTransform(lines, CRS("+init=epsg:32734"))

# To extract distances we need to coerce the lines to a psp object
linesppp <- lapply(lines_utm@lines, function(x){lapply(x@Lines, as.psp)})
linesppp <- do.call("c", linesppp)

# We also need to create points on the lines to extract average distances. We
# can do so by setting a regular distance (100 meters in this case)
linesppp <- lapply(linesppp, pointsOnLines, eps = 100)

# Load the merged water cover dataset
dat <- stack("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Prepare a list that stores a ppp layer for water (Code 1) for each floodmap
datppp <- suppressMessages(
  pbmclapply(
      X                   = 1:nlayers(dat)
    , mc.cores            = detectCores() / 2
    , ignore.interactive  = T
    , FUN                 = function(x){
      points <- rasterToPoints(dat[[x]], fun = function(z){z == 1}, spatial = T)
      points <- spTransform(points, CRS("+init=epsg:32734"))
      points <- as(points, "ppp")
      return(points)
      gc()
  })
)

# Calculate the average distance to water on the ppp object that is closest in
# date to the actual step
lines$DistanceToWater <- suppressMessages(
  pbmclapply(
      X                  = 1:nrow(lines)
    , mc.cores           = detectCores() / 2
    , ignore.interactive = T
    , FUN                = function(x){
      index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
      distance <- nncross(linesppp[[x]], datppp[[index]])
      distance <- mean(distance$dist)
      return(distance)
      gc()
  }) %>% do.call(rbind, .)
)

# Remove objects we don't need (release some memory)
rm(lines_utm, linesppp, datppp)
gc()

################################################################################
#### LandCover - Water
################################################################################
# Load the merged water cover dataset
dat <- stack("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Now we can extract the percentage cover of Water along each step
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# We only want to keep the values from the dates that are closest in time to the
# steps
lines$Water <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() / 2
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### LandCover - Trees
################################################################################
# Load the treecover map
dat <- stack("03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Now we can extract the percentage cover of water along each step
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Keep only values closest in date
lines$Trees <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() / 2
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### LandCover - Shrubs/Grassland
################################################################################
# Load the shrubcover map
dat <- stack("03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Extract the average shrub cover along each line
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Shrubs <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() / 2
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### Anthropogenic - Human Influence
################################################################################
# Load human influence data
dat <- stack("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")

# Extract values
extracted <- extrCov(dat, lines)

# Assign names
names(extracted) <- paste0("HumanInfluence", names(dat))

# Put the extracted values into the dataframe
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Storing
################################################################################
# Remove undesired columns
lines$drctn_p <- NULL

# Reorder the columns
lines@data <- dplyr::select(lines@data, c(
  , id
  , step_id_ = step_d_
  , State
  , case_
  , x1_
  , x2_
  , y1_
  , y2_
  , t1_
  , t2_
  , dt_
  , sl_
  , ta_
  , absta_
  , inactive = inactiv
  , everything()
  )
)

# To store the files we need to coerce the duration column to a numeric
lines$DistanceToWater <- as.vector(lines$DistanceToWater)

# Make "case_" and "inactive" logical
lines$case_ <- as.logical(lines$case_)
lines$inactive <- as.logical(lines$inactive)

# Prepare filenames
filename <- "00_General_Dispersers_POPECOL(SSF_Extracted)"

# Save the lines to a spatial lines dataframe
writeOGR(lines
  , "03_Data/02_CleanData"
  , filename
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

# Let's also store the data to a regular csv. We can use this file to restore
# the original column names since the ESRI shapefiles will store abbreviated
# names
write.csv(lines@data, paste0("03_Data/02_CleanData/", filename, ".csv"))
