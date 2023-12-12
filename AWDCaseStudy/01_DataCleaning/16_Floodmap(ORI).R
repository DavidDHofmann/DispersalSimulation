################################################################################
#### Preparation of the Ori Floodmaps
################################################################################
# Description: Resampling and preparation of the ORI floodmaps. We will also
# check for which dates we need to download additional floodmaps.

# Clean environment
rm(list = ls())

# load packages
library(tidyverse)    # For data wrangling
library(raster)       # To handle raster data
library(lubridate)    # To handle dates

# Allow for multicore use
beginCluster()

################################################################################
#### Compare Dispersal and Ori Dates
################################################################################
# Identify the ori files in the working directory
files <- dir(
    path        = "03_Data/01_RawData/ORI/01_Maps"
  , pattern     = "*.tif$"
  , full.names  = T
)

# Create a vector that identifies the dates for which we have ORI data
ori_dates <- files %>%
  basename() %>%
  substr(start = 1, stop = 10) %>%
  as.Date()

# Put filenames of their respective date into a dataframe
files <- data.frame(Filename = files, Date = ori_dates)

# Load the disperser's tracks
dispersal <- "03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv" %>%
  read_csv() %>%
  subset(State == "Disperser")

# Extract unique dates
dis_dates <- unique(dispersal$Timestamp) %>% as.Date()

# For all these dispersal dates we want to know the closest ORI dates
closest1 <- as.Date(NA)
closest2 <- as.Date(NA)
for (i in 1:length(dis_dates)){
  closest1[i] <- ori_dates[which(abs(dis_dates[i] - ori_dates) ==
    min(abs(dis_dates[i] - ori_dates)))][1]
  closest2[i] <- ori_dates[which(abs(dis_dates[i] - ori_dates) ==
    min(abs(dis_dates[i] - ori_dates)))][2]
}

# Put the dates together
dates <- data.frame(
    Dispersal   = dis_dates
  , ClosestORI1 = closest1
  , ClosestORI2 = closest2
  , Difference  = abs(dis_dates - closest1)
)

# Look at the result
dates

# Let's see how close we get with the ori data
hist(as.numeric(dates$Difference))

# For some fixes we are really far away
summary(as.numeric(dates$Difference))

# Let's see for which periods we are most off
off <- subset(dates, Difference > 50)
arrange(off, -Difference)

# Let's identify the dates for which we would prefer more precise data. Let's
# say that 3 weeks should be the largest gap
missing <- subset(dates, Difference > 21) %>% arrange(-Difference)

# We can store the dates for which we would like to have additional floodmaps.
# This will be helpful in the script in which we download additional MODIS
# images
write.csv(missing, "03_Data/01_RawData/ORI/MissingDates.csv")

################################################################################
#### Resample ORI layers to match MODIS layers
################################################################################
# We need to make sure that ORI and MODIS images have the same extent and origin
# To achieve this I created a reference raster.
files <- dir(
    path        = "03_Data/01_RawData/ORI/01_Maps"
  , pattern     = "*.tif$"
  , full.names  = T
)
ori <- stack(files, bands = 1)

# Also load the reference map
ref <- dir(
    path        = "03_Data/01_RawData/ORI"
  , pattern     = "Floodmap"
  , full.names  = T
)[1] %>% raster()

# Resample and store each ORI layer
dir.create("03_Data/02_CleanData/00_Floodmaps")
dir.create("03_Data/02_CleanData/00_Floodmaps/01_Original")
newnames <- substr(names(ori), start = 2, stop = 11) %>% as.Date(format = "%Y.%m.%d")
newnames <- paste0("03_Data/02_CleanData/00_Floodmaps/01_Original/", newnames, ".tif")
for (i in 1:nlayers(ori)){
  ori_res <- resample(ori[[i]], ref, "ngb")
  writeRaster(ori_res, newnames[i], overwrite = TRUE)
  gc()
  cat(i, "out of", nlayers(ori), "done...\n")
}

# Terminate the cluster
endCluster()
