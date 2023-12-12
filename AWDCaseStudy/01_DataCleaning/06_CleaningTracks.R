################################################################################
#### Cleaning and Preparing Tracks from GPS Fixes
################################################################################
# Description: In this script I put together all gps fixes, regardless of
# whether an individual is a disperser or not. I then clean the data and add
# additional information, such as dispersal dates and pack-affiliation.

# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates easily
library(pbmcapply)  # To make use of multiple cores
library(davidoff)   # Custom functions
library(raster)     # To handle spatial data
library(rgdal)      # To handle spatial data
library(plotKML)    # To store the cleaned tracks to a kml
library(spacetime)  # To store the cleaned tracks to a kml (required for slider)
library(viridis)    # For nice colors

################################################################################
#### Data Source 1: POPECOL
################################################################################
# Identify all files containing GPS data
files <- dir(
    path        = "03_Data/01_RawData/POPECOL/01_GPS"
  , pattern     = "GPS_Collar"
  , full.names = T
)

# There are some inconsistencies we need to take care of. Firstly, sometimes "."
# is used as decimal, sometimes ",". Moreover, there are some files containing
# "°" symbols. We need to get rid of this stuff. This only needs to be run once
# as it overwrites the default csv files. This also makes sure that every file
# is stored exactly the same way
lapply(files, function(x){

  # Identify separator
  sep <- readLines(x, n = 1)
  sep <- if_else(grepl(x = sep, pattern = ";", useBytes = T), ";", ",")

  # Load the file as plain text and remove funny characters
  dat <- read_file(x, local = locale(encoding = "latin1"))
  if (sep == ","){
      dat <- gsub(dat, pattern = ",", replacement = ";")
    } else {
      dat <- gsub(dat, pattern = ",", replacement = ".")
  }
  dat <- gsub(dat, pattern = "°", replacement = "")
  dat <- gsub(dat, pattern = "/", replacement = ".")
  write(dat, x)
})

# Load all of them and do some cleaning
dat1 <- lapply(files, function(x){

  # Load data...
  dat <- x %>%

    # ...using readr
    read_delim(local = locale(encoding = "latin1"), delim = ";") %>%

    # Retrieve DogName from filename. For this we use regex to identify a
    # pattern preceeded by five digits (?<=\\d{5}) and a "_", but is then
    # followed by a "_" and whatever
    mutate(DogName = str_extract(x, pattern = "(?<=\\d{5}_).*(?=_)")) %>%

    # Retrieve timestamp
    mutate(Timestamp = as.POSIXct(
      paste(UTC_Date, UTC_Time), tz = "UTC", format = "%d.%m.%Y %H:%M:%S")
    ) %>%

    # Remove special characters like [°]
    setNames(gsub(names(.), pattern = " \\[*", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "°", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "\\]", replacement = "")) %>%

    # Keep only desired columns
    dplyr::select(.
      , DogName   = DogName
      , DOP       = DOP
      , CollarID  = CollarID
      , x         = Longitude
      , y         = Latitude
      , Timestamp = Timestamp
    ) %>%

    # Add a column indicating the data source
    mutate(Source = "Popecol")
}) %>% do.call(rbind, .)

# Check validity of columns
summary(dat1$x)
summary(dat1$y)
sum(is.na(dat1$Timestamp))
sum(is.na(dat1$CollarID))
range(dat1$Timestamp)

# The GPS devices typcially record GPS locations even before and after collaring
# an animal. These GPS fixes need to be removed since they provide wrong
# location data. To do so, we load the dataframe that Dominik prepared. It shows
# for each collar and dog the first and last fix retrieved, i.e. when an
# individual was collared or uncollared. The table also contains valuable
# information about pack-affiliations
collars <- read_csv2("03_Data/01_RawData/POPECOL/CollarSettings.csv") %>%

  # Remove rows where the Dog Names are NA
  subset(., !is.na(`Dog Name`)) %>%

  # Keep only the desired columns and rename them nicely. Note that we want to
  # keep the Dog Code because it provides information about the birth pack of
  # each dog.
  dplyr::select(.
    , CollarID  = `Collar Nr.`
    , DogName   = `Dog Name`
    , DogCode   = `Dog Code`
    , Sex       = `Sex`
    , FirstDate = `Collaring Date`
    , LastDate1 = `Stop recording date`
    , LastDate2 = `Last fix date`
  )

# We are missing exact times for some of the timestamps. This will cause errors
# when we convert those to posixct. We Will therefore assign very conservative
# times (i.e. 24:00 for first dates, 00:01 for last dates)
for (i in 1:nrow(collars)){
  if (!is.na(collars$FirstDate[i]) & nchar(collars$FirstDate[i]) == 10){
    collars$FirstDate[i] <- paste0(collars$FirstDate[i], " 24:00")
  }
  if (!is.na(collars$LastDate1[i]) & nchar(collars$LastDate1[i]) == 10){
    collars$LastDate1[i] <- paste0(collars$LastDate1[i], " 00:01")
  }
  if (!is.na(collars$LastDate2[i]) & nchar(collars$LastDate2[i]) == 10){
    collars$LastDate2[i] <- paste0(collars$LastDate2[i], " 00:01")
  }
}

# Taryn's collar number is wrong
collars$CollarID[collars$DogName == "Taryn" & collars$CollarID == 22028] <- 20228

# Now we can coerce the date columns to true dates. Note that we subtract two
# hours because the original format was in local time (which is UTC + 2 hours).
# Also note that there are two possible "last dates". The first one refers to
# the last fix, the other to the date when the dog was uncollared. We need to
# combine both because depending on whether the dog is still collared we will
# use on or the other. Again we subtract two hours to make sure that the times
# are in UTC
collars <- collars %>%
  mutate(.
    , FirstDate = as.POSIXct(FirstDate
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
    , LastDate1 = as.POSIXct(LastDate1
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
    , LastDate2 = as.POSIXct(LastDate2
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
  )

# Merge the columns for the "LastDates". We only want to keep the latest of the
# two.
collars$LastDate <- pmax(collars$LastDate1, collars$LastDate2, na.rm = T)
collars$LastDate <- collars$LastDate + minutes(5)
collars$LastDate1 <- NULL
collars$LastDate2 <- NULL

# Logical check that LastDate > FirstDate and make sure there are no duplicates
table(collars$LastDate > collars$FirstDate)
table(paste(collars$CollarID, collars$DogName))
table(table(paste(collars$CollarID, collars$DogName)))

# Create a column that indicates if the respective fix lies within the first and
# last date
dat1 <- left_join(dat1, collars, by = c("CollarID", "DogName"))
dat1$Keep <- ifelse(
    test  = dat1$Timestamp >= dat1$FirstDate & dat1$Timestamp <= dat1$LastDate
  , yes   = T
  , no    = F
)

# Let's check how many datapoints we loose (per dog and collar)
table(dat1$Keep)
table(dat1$DogName, dat1$Keep)
table(dat1$CollarID, dat1$Keep)

# Subset data and remove unnecessary columns
dat1 <- subset(dat1, Keep)
dat1 <- dplyr::select(dat1, -c("FirstDate", "LastDate", "Keep", "DogCode"))

# Prepare data for visualization
vis <- dat1 %>% group_by(DogName, CollarID) %>% summarize(
    First = range(Timestamp)[1]
  , Last  = range(Timestamp)[2]
)

# Verify that collars are not overlapping
ggplot(vis, aes(color = factor(CollarID))) +
  geom_segment(
      aes(x = First, xend = Last, y = DogName, yend = DogName)
    , size = 3
    , alpha = 0.6
  )

################################################################################
#### Data Source 2: Abrahms
################################################################################
# We also want to import the GPS fixes provided by Abrahms. Note that I will not
# clean this data anymore, as this is the cleaned data that Abrahms used for her
# publication. Her timestamps are in UTC already.
dat2 <- read_csv("03_Data/01_RawData/ABRAHMS/DispersalPaths.csv") %>%

  # Remove rows with missing fixes
  filter(., !is.na(`Longitude..deg.`)) %>%

  # Select only the columns of interest
  dplyr::select(.,
      DogName   = `id`
    , x         = `Longitude..deg.`
    , y         = `Latitude..deg.`
    , Timestamp = `Timestamp`
  ) %>%

  # Add columns that are also in dat1 (so we can bind the data afterwards). Note
  # that all individuals from Abrahms are males.
  mutate(., CollarID = NA, DOP = NA, Sex = "M", Source = "Abrahms")

################################################################################
#### Data Source 3: Botswana Camp
################################################################################
# Lastly, we retrieved fixes from the camp in Botswana that were collected prior
# to own our project in Botswana. This is data from residents only! I only
# include it for completeness. This data requires some intensive cleaning. Also
# note that the data is in UTC already.
files <- dir(
    path        = "03_Data/01_RawData/DOMINIK/GPS"
  , pattern     = ".txt$"
  , full.names  = T
)

# Load all the files
dat3 <- lapply(files, function(x){

  # Identify the DogName
  name <- substr(basename(x), start = 6, stop = nchar(basename(x)) - 7)

  # Load the data (skip the first two rows as they contain stuff we dont care
  # about)
  read_delim(x, delim = "\t", skip = 2) %>%

    # Keep only the columns of interest
    dplyr::select(
        x                   = `Longitude (deg)`
      , y                   = `Latitude (deg)`
      , Timestamp           = `UTC time (yyyy-mm-dd HH:MM:SS)`
      , Height              = `Height above MSL (m)`
      , HorizontalAccuracy  = `Horizontal accuracy (m)`
      , VerticalAccuracy    = `Vertical accuracy (m)`
    ) %>%

    # Add dog name
    mutate(DogName = name) %>%

    # Prepare columns present in dat1
    mutate(., CollarID = NA, DOP = NA, Sex = "M", Source = "Botswana") %>%

    # Remove rows with missing fixes
    filter(., !is.na(x)) %>%

    # Order the data by timestamp
    arrange(., Timestamp) %>%

    # Create Column that indicates the timelag and distance between two fixes
    # This allows us to also calculate the speed
    mutate(x,
        dt = as.numeric(Timestamp - lag(Timestamp), units = "hours")
      , dl = sqrt((x - lag(x))**2 + (y - lag(y)) ** 2) * 111
      , speed = dl / dt
    )
}) %>% do.call(rbind, .)

# There are some pretty beefy outliers that we need to remove. We can use the
# Accuracy indicators, the height and speed in order to remove any unreasonable
# gps fix (we are going to be pretty rough here)
dat3 <- dat3 %>%
  filter(., Height > quantile(Height, 0.1)) %>%
  filter(., Height < quantile(Height, 0.9)) %>%
  filter(., VerticalAccuracy < quantile(VerticalAccuracy, 0.9)) %>%
  filter(., HorizontalAccuracy < quantile(HorizontalAccuracy, 0.9)) %>%
  filter(., year(Timestamp) != 1970) %>%
  filter(., speed < quantile(speed, 0.9, na.rm = TRUE)) %>%
  dplyr::select(.
    , c("DogName", "CollarID", "x", "y", "Timestamp", "DOP", "Sex", "Source")
  )

# Put data of all sources together
data <- rbind(dat1, dat2, dat3)

# Check sources
table(data$Source)

# Check if there are any duplicates
dups_complete <- duplicated(data)
sum(dups_complete)
dups_incomplete <- duplicated(data[, c("DogName", "Timestamp")])
sum(dups_incomplete)

# Remove them
data <- subset(data, !dups_incomplete)

################################################################################
#### Adding Cutoff Dates
################################################################################
# For each GPS fix we want to know whether it was collected during dispersal or
# residency. We can use the following table for this.
cut <- read_csv("03_Data/01_RawData/POPECOL/CutoffDates.csv") %>%

  # Remove undesired columns
  dplyr::select(c(DogName, StartDate = StartDate_UTC, EndDate = EndDate_UTC)) %>%

  # Sort
  arrange(DogName, StartDate)

# As you can see for some individuals there are multiple dispersal phases. Let's
# create a counter for each individual
cut <- cut %>%
  group_by(DogName) %>%
  mutate(DispersalNo = row_number()) %>%
  ungroup()

# We can now check hof often each individual "dispersed"
cut %>%
  group_by(DogName) %>%
  summarize(max(DispersalNo))

# Store the merged cutoff dates
write.csv(cut, "03_Data/02_CleanData/00_General_CutoffDates_POPECOL.csv")

# We can use the table to create a new column that indicates the state of the
# animal. By default we say that our individuals are resident.
data$State <- "Resident"

# Now loop through all dogs and use the cutoff table to specify whether a fix
# was taken during dispersal or not
names <- unique(cut$DogName)
for (i in seq_along(names)){
  cutoff <- subset(cut, DogName == names[i])
  index <- which(data$DogName == names[i])
  for (h in 1:nrow(cutoff)){
    data$State[index][data$Timestamp[index] >= cutoff$StartDate[h] &
    data$Timestamp[index] <= cutoff$EndDate[h]] <- "Disperser"
  }
}

# Look at the final table
head(data)
tail(data)

# Check the number of data per individual
table(data$DogName, data$State)

# Order the data nicely
data <- data %>% arrange(DogName, Timestamp)

# Note: For Abrahms individuals there is an overlap between data collected by
# Abrahms and data collected by the Staff in Botswana (its actually the same
# data). Because of a tiny mismatch in the timestamps they are not recognized as
# duplicates. However, the coordinates align perfectly and the temporal mismatch
# is minor. Thus, the issue will be resolved once we resample the data to a
# coarser resolution.

################################################################################
#### Visualize Dispersal Phases
################################################################################
# We now want to create a plot illustrating the dispersal phases. For this, we
# need to create a table indicating the start and end of each "phase"
phases <- data %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    phase <- x$State != lag(x$State)
    phase <- replace_na(phase, F)
    x$Phase <- cumsum(phase)
    return(x)
  })) %>%
  unnest(cols = data) %>%
  group_by(DogName, State, Phase) %>%
  summarize(FirstDate = min(Timestamp), LastDate = max(Timestamp)) %>%
  mutate(DogName = as.factor(DogName))

# Reorder the factors
# phases$DogName <- fct_reorder(phases$DogName, phases$DogName, .desc = T)

# Visualize them
ggplot(phases, aes(x = rbind(FirstDate, LastDate), y = DogName)) +

  # Add segments for the resident phase
  geom_segment(data = subset(phases, State == "Resident"), aes(
      x     = FirstDate
    , xend  = LastDate
    , y     = DogName
    , yend  = DogName
  ), size = 2.5, color = "cornflowerblue") +

  # Add segments for the resident phase
  geom_segment(data = subset(phases, State == "Disperser"), aes(
      x     = FirstDate
    , xend  = LastDate
    , y     = DogName
    , yend  = DogName
  ), size = 1.0, color = "black") +

  # Reverse the y scale
  scale_y_discrete(limits = rev(levels(phases$DogName))) +

  # Put a useful title and axis labels
  ggtitle("GPS Observations") +
  xlab("Date") +
  ylab("Name")

################################################################################
#### Resampling Data
################################################################################
# Only keep individuals that eventually dispersed
dispersers <- unique(subset(data, State == "Disperser")$DogName)
data <- subset(data, DogName %in% dispersers)

# Remove NA fixes
table(is.na(data$x))
table(is.na(data$y))
data <- subset(data, !is.na(x) & !is.na(y))

# Sort data
data <- arrange(data, DogName, Timestamp)

# Some of the data has a very high resolution that we don't need. We will
# therefore subsample to a resolution we can work with.
data <- data %>% group_by(DogName) %>% nest()

# Resmaple
data$data <- suppressMessages(
  pbmclapply(1:nrow(data)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x){
      resFix(data$data[[x]], hours = 1, start = 1, tol = 0.5)
    }
  )
)
data <- unnest(data)

# Check out nrow of remaining data
nrow(data)

################################################################################
#### Store the Output (as csv and shapefile)
################################################################################
# Write the data to file
write_csv(data, "03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

# Create SpatialLinesDataFrame for visualization
data <- data %>% group_by(DogName) %>% nest()
data$Tracks <- suppressMessages(
  pbmclapply(1:nrow(data)
    , ignore.interactive = T
    , mc.cores = detectCores() - 1
    , function(x){
      coords <- data$data[[x]]
      coords$DogName <- data$DogName[x]
      x <- coords
      coordinates(x) <- c("x", "y")
      crs(x) <- CRS("+init=epsg:4326")
      lines <- spLines(x)
      lines <- createSegments(lines)
      lines <- as(lines, "SpatialLinesDataFrame")
      lines@data <- x@data[1:(nrow(x) - 1), ]
      crs(lines) <- CRS("+init=epsg:4326")
      return(lines)
    }
  )
)

# Create SpatialPointsDataFrame for visualization
data$Points <- lapply(data$data, function(x){
  coordinates(x) <- c("x", "y")
  crs(x) <- CRS("+init=epsg:4326")
  return(x)
})

# Store shapefile
tracks <- do.call(rbind, data$Tracks)
writeOGR(tracks
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Dispersers_POPECOL"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Create KML Files
################################################################################
# We may also want to create kml files for visualization in google earth. Let's
# first point to the shapes visualized in google earth.
shape1 <- "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"
shape2 <- "http://maps.google.com/mapfiles/kml/shapes/square.png"

# Prepare a folder into which we will store the kmls
dir.create("03_Data/02_CleanData/00_KML")

# Loop through the individuals and prepare kml files for them
lapply(1:nrow(data), function(x){

  # Access required data
  name    <- data$DogName[[x]]
  points  <- data$Points[[x]]
  track   <- data$Tracks[[x]]

  # Identify start and endpoints
  first <- as(points[1, ], "SpatialPoints")
  last  <- as(points[nrow(points), ], "SpatialPoints")

  # Make row-names valid
  row.names(points) <- as.character(1:nrow(points))
  row.names(track) <- as.character(1:nrow(track))

  # Create spacetime object from points
  points_st <- STIDF(
      sp   = as(points, "SpatialPoints")
    , time = points$Timestamp
    , data = points@data
  )

  # Create spacetime object from track
  track_st <- STIDF(
      sp   = as(track, "SpatialLines")
    , time = track$Timestamp
    , data = track@data
  )

  # Generate a name for the kml file
  filename <- paste0("03_Data/02_CleanData/00_KML/", name, ".kml")

  # Generate kml file (note that for some individuals we can't produce an
  # info-table because there are too many GPS fixes)
  kml_open(filename)
  kml_layer(first
    , colour     = "green"
    , size       = 1
    , shape      = shape2
    , LabelScale = 0
  )
  kml_layer(last
    , colour     = "red"
    , size       = 1
    , shape      = shape2
    , LabelScale = 0
  )
  kml_layer(track_st
    , colour       = State
    , colour_scale = rev(viridis(2, begin = 0.6))
    , width        = 2.5
  )
  tryCatch(kml_layer(points_st
    , colour       = State
    , colour_scale = rev(viridis(2, begin = 0.6))
    , size         = 0.75
    , balloon      = T
    , shape        = shape1
    , LabelScale   = 0
  ), error = function(e){
    kml_layer(points_st
      , colour       = State
      , colour_scale = rev(viridis(2, begin = 0.6))
      , size         = 0.75
      , balloon      = F
      , shape        = shape1
      , LabelScale   = 0
    )
  })
  kml_close(filename)
})
