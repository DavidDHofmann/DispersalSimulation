################################################################################
#### Resampling MODIS Maps and Preparing for Merge
################################################################################
# Description: We now have two different sources for floodmaps: ORI and MODIS.
# We want to merge them afterwards, so let's prepare everything.

# clear r's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates nicely
library(raster)       # To handle raster data
library(terra)        # To handle reaster data quickly

################################################################################
#### Identify Floodmaps to Resample
################################################################################
# Load all classified maps
files <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = "*.tif$"
  , full.names  = T
)
flood <- rast(files)

# Extract the dates from the layer descriptions
flood_dates <- names(flood) %>%
  substr(start = 1, stop = 10) %>%
  as.Date(., format = "%Y-%m-%d")

# Load all dispersal tracks
disp <- "03_Data/02_CleanData/00_General_Dispersers_Popecol.csv" %>%
  read.csv() %>%
  subset(., State == "Disperser")

# Extract the unique dates from the dispersers' dataframe
dis_dates <- unique(disp$Timestamp) %>% as.Date(.) %>% sort(.)

# For all these dispersal dates we want to know the closest flood dates
closest1 <- as.Date(NA)
closest2 <- as.Date(NA)
for (i in 1:length(dis_dates)){
  closest1[i] <- flood_dates[which(abs(dis_dates[i] - flood_dates) ==
    min(abs(dis_dates[i] - flood_dates)))][1]
  closest2[i] <- flood_dates[which(abs(dis_dates[i] - flood_dates) ==
    min(abs(dis_dates[i] - flood_dates)))][2]
}

# Put the dates together
dates <- data.frame(
    Dispersal   = dis_dates
  , Closest1    = closest1
  , Closest2    = closest2
  , Difference  = abs(dis_dates - closest1)
)
arrange(dates, -Difference)

# Look at the result
summary(as.numeric(dates$Difference))

# Let's now subset to the floodmaps which are closest to some dispersal event.
# Let's select only these maps that are closest to some dispersal event
indices <- which(flood_dates %in% dates$Closest1)
length(indices)

# I also want to conserve the last floodmap
indices <- c(indices, length(flood_dates))
flood <- flood[[indices]]
files <- files[indices]

# Reclassify the floodmaps to values 0 and 1
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
flood <- classify(flood, rcl)

################################################################################
#### Resample Floodmaps
################################################################################
# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
r <- crop(r, flood)

# Create directory for resampled floodmaps
dir.create("03_Data/02_CleanData/00_Floodmaps/02_Resampled")

# Loop through the floodmaps, resample and store them
for (i in 1:length(files)){
  map <- terra::resample(flood[[i]], r250, "near")
  newname <- substr(files[i], start = 47, stop = 56) %>% paste0(., ".tif")
  newname <- paste0("03_Data/02_CleanData/00_Floodmaps/02_Resampled/", newname)
  map <- raster(map)
  writeRaster(map, newname, overwrite = TRUE)
  cat(i, "of", length(files), "done...\n")
}
