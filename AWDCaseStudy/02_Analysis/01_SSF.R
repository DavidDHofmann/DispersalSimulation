################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we coerce our gps data to steps and generate
# random steps for the (time varying) integrated step selection analysis.

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

# Set seed for reproducability (need to change random sampler for parallel)
set.seed(1234)

################################################################################
#### Data Cleaning
################################################################################
# Load the gps data
data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

# We only need dispersers' data
data <- subset(data, State == "Disperser")

# Let's create a timestamp that is rounded to the nearest hour
data$TimestampRounded <- round_date(data$Timestamp, "1 hour")

# Make sure there are no duplicates!
table(duplicated(data[, c("DogName", "TimestampRounded")]))

# Count number of dogs
length(unique(data$DogName))

# Nest data by dog
data <- data %>% group_by(DogName) %>% nest()

# Resample to two (minimally) 2 hours
data$data <- suppressMessages(
  pbmclapply(1:nrow(data)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x){
      resFix(data$data[[x]], hours = 2, start = 1, tol = 0.5)
    }
  )
)
data <- data %>% unnest(data)

# Keep only fixes on our regular scheme
data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))

# We're going to convert the fixes to steps, yet we want to indicate a new
# burst if a step takes longer than 4.25 hours (or longer than 8.25 hours
# between 07:00 and 15:00). Let's thus calculate the timelags.
data <- data %>%
  group_by(DogName) %>%
  mutate(dt_ = Timestamp - lag(Timestamp)) %>%
  mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
  ungroup()

# Indicate when a new burst starts (a new burst also starts whenever dt_ = NA)
# In case we're using iSSF
data <- data %>%
  mutate(NewBurst = ifelse(
    dt_ > 8.25 |
    (dt_ > 4.25 & hour(TimestampRounded) != 15) |
    is.na(dt_), yes = 1, no = 0)
  ) %>%
  mutate(BurstID = cumsum(NewBurst)) %>%
  dplyr::select(-c(NewBurst))

# We can only work with bursts that contain at least three fixes
data <- data %>%
  group_by(DogName, BurstID) %>%
  nest() %>%
  mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
  subset(Nrow >= 3) %>%
  unnest(data) %>%
  ungroup() %>%
  dplyr::select(-c(Nrow, dt_))

# Create burst IDs that are unique
data$BurstID <- data %>% group_indices(DogName, BurstID)

# Now we can coerce the data to proper steps. We'll trick amt::make_track by
# using the burstID as animalID
tracks <- data %>%
  make_track(.
    , .x      = x
    , .y      = y
    , .t      = Timestamp
    , id      = BurstID
    , crs     = CRS("+init=epsg:4326")
    , State   = State
    , DogName = DogName
  ) %>%

  # Transform the tracks to utm
  transform_coords(CRS("+init=epsg:32734")) %>%

  # Nest the tracks
  nest(data = -"id") %>%

  # Turn to a step representation (look up the amt vignette for details)
  mutate(data = map(data, function(x){
    x %>%

      # The function creates steps from the fixes in each tibble row. The
      # option "keep_cols" allows to keep the time of day (tod) column that we
      # added
      steps(keep_cols = "start") %>%

      # Transform the difftime column to a numeric column to avoid that there
      # are heterogeneous units
      mutate(., dt_ = as.numeric(dt_, units = "hours"))
  })) %>%

  # Unnest the tibble
  unnest(data) %>%

  # Multiply turning angles with negative one (for some reason the package
  # calculates the turning angles conuterclockwise but we want them clockwise)
  mutate(ta_ = ta_ * (-1)) %>%

  # Add a column indicating the absolute turning angle (important for the
  # ellipses). We can use the function we created above.
  mutate(absta_ = absAngle(.)) %>%

  # We can only work with steps for which we have a turning angle. Let's get
  # rid of any steps where the turning angle is NA
  filter(!is.na(ta_))

# Generate unique step id and indicate that steps are observed (case_) steps
tracks <- tracks %>% mutate(
    step_id_ = 1:nrow(.)
  , case_    = T
)

# Remove burst ids
tracks$id <- tracks$DogName
tracks$DogName <- NULL

# Indicate if a fix was taken during time of activity or inactivity (we define
# inactivity as anything between 07:00 and 15:00)
tracks <- tracks %>% mutate(
  inactive =
    hour(round_date(t1_, "1 hour")) >= 7 &
    hour(round_date(t2_, "1 hour")) <= 15 &
    hour(round_date(t1_, "1 hour")) <
    hour(round_date(t2_, "1 hour"))
)

################################################################################
#### Generation of Random Steps
################################################################################
# We can't work with 0 step lengths or step speeds
tracks$sl_[tracks$sl_ == 0] <- 1

# Fit a Gamma distribution to the step speed (again, fitting a gamma to step
# speeds or fitting a gamma to (normalized) step lengths yields the same under
# iSSF) -> almost at least
sl <- fit_distr(tracks$sl_, "gamma")

# Let's visualize the fit
x <- seq(0, max(tracks$sl_), 1)
y <- dgamma(x, scale = sl$params$scale, shape = sl$params$shape)
hist(tracks$sl_, freq = F, breaks = 100)
lines(y ~ x, col = "red")

# Write the distribution to file
write_rds(sl, "03_Data/03_Results/99_GammaDistribution.rds")

# Define the number of random steps
nsteps <- 24

# Prepare an empty list in which we can store the random steps of each individual
randomSteps <- list()

# Now we use the fitted gamma distribution to sample step lengths. For the
# turning angles we follow Avgar et al. 2016 and use a uniform distribution
for (i in 1:nrow(tracks)){

  # Draw random turning angles
  ta_new <- runif(nsteps, min = -pi, max = pi)

  # Draw random step lengths with the fitted parameters
  sl_new <- rgamma(nsteps
    , shape = sl$params$shape
    , scale = sl$params$scale
  )

  # Put the step lengths and turning angles into a new dataframe and calculate the
  # new endpoints of each random step. We also indicate that the steps are control
  # steps (i.e. 0)
  randomSteps[[i]] <- tracks[rep(i, nsteps), ] %>%
    mutate(.
      , absta_  = absta_ + (ta_new - ta_)
      , ta_     = ta_new
      , sl_     = sl_new
      , case_   = 0
      , x2_     = x1_ + sin(absta_) * sl_
      , y2_     = y1_ + cos(absta_) * sl_
    )
}

# Collapse the list of dataframes into a single dataframe
randomSteps <- do.call(rbind, randomSteps)

# We need to make sure that the absolute turning angle ranges from 0 to 2 * pi
randomSteps$absta_[randomSteps$absta_ > 2 * pi] <-
  randomSteps$absta_[randomSteps$absta_ > 2 * pi] - 2 * pi
randomSteps$absta_[randomSteps$absta_ < 0] <-
  randomSteps$absta_[randomSteps$absta_ < 0] + 2 * pi

# Merge the dataframes of the observed and random steps
ssf <- rbind(tracks, randomSteps) %>%
  arrange(step_id_, desc(case_)) %>%
  transform(case_ = as.logical(case_))

# Prepare a filename
filename <- "00_General_Dispersers_POPECOL(SSF)"

# Store to file
write_csv(ssf, paste0("03_Data/02_CleanData/", filename, ".csv"))

# Coerce the steps to spatial lines. We can use the function we defined
# earlier for this
lines <- lineTrack(ssf, CRS("+init=epsg:32734"))

# Transform the lines to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

# Coerce the duration column to a format that can be stored (numeric rather
# than difftime)
lines$dt_ <- as.numeric(lines$dt_)

# Store the lines
writeOGR(lines
  , dsn       = "03_Data/02_CleanData"
  , layer     = filename
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
