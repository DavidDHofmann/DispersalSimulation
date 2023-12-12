################################################################################
#### Connectivity Validation
################################################################################
# Description: Here, we use independent dispersal data to validate our
# predictions of landscape connectivity. Specifically, we will use
# step-selection functions to determine whether dispersing individuals have a
# tendency to stick to areas of high connectivity.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(lubridate)      # To handle timestamps
library(wilddogr)       # To download dispersal data
library(amt)            # Some movement metrics
library(davidoff)       # Access to custom functions
library(terra)          # To handle spatial data
library(raster)         # To handle spatial data
library(pbmcapply)      # For multicore abilities


# Define number of random paths / random steps and the maximally allowed
# distance at which paths may be shifted
n_rand   <- 50
shiftmax <- 50

################################################################################
#### Download Latest Dispersal Data
################################################################################
# # Let's access data from the most recent dispersal events from dropbox
# # rdrop2::drop_auth(cache = F)
# files <- dog_files(rvc = F)
#
# # We are only interested in some of the dispersers
# keep <- c("Aspen", "Carson", "Chiounard", "Dell", "Earth", "Encinitas", "Ripley"
#   , "Saturday", "Sishen")
# # keep <- c("Carson", "Dell", "Earth", "Encinitas", "Ripley", "Saturday")
# dispersers <- subset(files, DogName %in% keep)
#
# # Download their data to a temporary file
# downloaded <- dog_download(dispersers
#   , outdir    = tempdir()
#   , printpath = T
#   , overwrite = T
#   , clean     = T
# )
#
# # Note that we are only interested in dispersal phases
# dat <- downloaded %>%
#   read_csv() %>%
#   subset(State == "Disperser" & !is.na(x) & !is.na(y))
#
# # Let's only keep individuals for which there are at least 10 datapoints
# keep <- dat %>%
#   count(DogName) %>%
#   subset(n >= 10) %>%
#   pull(DogName)
# dat <- subset(dat, DogName %in% keep)
#
# # Store those to file
# write_csv(dat, "03_Data/02_CleanData/99_ValidationDispersalData.csv")
dat <- read_csv("03_Data/02_CleanData/99_ValidationDispersalData.csv")

# Summarize number of fixes by coalition
dat %>%
  group_by(DogName) %>%
  summarize(n = n()) %>%
  summarize(mean = mean(n), sd = sd(n))

################################################################################
#### Preparation for Path Selection Analysis
################################################################################
# Create observed paths for each individual
psf <- dat %>%
  nest(Data = -DogName) %>%
  mutate(
      FirstDate   = map(Data, function(x) {min(x$Timestamp)}) %>% do.call(c, .)
    , LastDate    = map(Data, function(x) {max(x$Timestamp)}) %>% do.call(c, .)
    , Duration    = round(difftime(LastDate, FirstDate, units = "days"))
    , NumberFixes = map(Data, nrow) %>% do.call(c, .)
  ) %>%
  mutate(Path = map(Data, function(x) {
    path <- spLines(cbind(x$x, x$y), crs = "+init=epsg:4326")
    return(path)
  })) %>%
  mutate(case_ = 1)

# Function to randomize a path
randomizePath <- function(line, shiftmax) {
  coords <- coordinates(line)[[1]][[1]]
  mean_x <- mean(coords[, 1])
  mean_y <- mean(coords[, 2])
  shiftdist  <- runif(1, 0, shiftmax ** 2)
  angle      <- runif(1, 0, 2*pi)
  shift_x    <- sqrt(shiftdist) * cos(angle)
  shift_y    <- sqrt(shiftdist) * sin(angle)
  rotate     <- runif(1, 0, 360)
  newpath <- maptools::elide(line
    , rotate = rotate
    , center = c(mean_x, mean_y)
    , shift  = c(shift_x, shift_y)
  )
  crs(newpath) <- crs(line)
  return(newpath)
}

# Generate 100 random paths for each observed path
psf_random <- lapply(1:nrow(psf), function(x) {
  path  <- psf$Path[[x]]
  paths <- lapply(1:n_rand, function(y) {
    randomizePath(path, shiftmax = shiftmax / 111)
  })
  result       <- psf[rep(x, n_rand), ]
  result$Path  <- paths
  result$case_ <- F
  return(result)
}) %>% do.call(rbind, .)

# Combine observed and random paths
psf <- rbind(psf, psf_random) %>% arrange(DogName, desc(case_))

# Combine them into a single object
psf_data  <- dplyr::select(psf, -c(Data, Path))
psf       <- as(do.call(rbind, psf$Path), "SpatialLinesDataFrame")
psf@data  <- psf_data

# Let's write the observed paths to file
writeVector(vect(psf), "03_Data/02_CleanData/99_ValidationPaths.shp", overwrite = T)

################################################################################
#### Prepare Data for Step Selection Analysis
################################################################################
# Resample the data to 4 hours and only keep those fixes that were collected at
# certain times. Finally, we compute the duration of each step and determine
# bursts.
ssf <- dat %>%
  resampleFixes(hours = 4, start = 3, tol = 0.5) %>%
  subset(!is.na(Timestamp)) %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  subset(hour(TimestampRounded) %in% c(3, 7, 15, 19, 23)) %>%
  group_by(DogName) %>%
  mutate(dt_ = difftime(Timestamp, lag(Timestamp), units = "hours")) %>%
  ungroup() %>%
  mutate(NewBurst = ifelse(
    dt_ > 8.25 |
    (dt_ > 4.25 & hour(TimestampRounded) != 15) |
    is.na(dt_), yes = 1, no = 0)
  ) %>%
  mutate(BurstID = cumsum(NewBurst)) %>%
  dplyr::select(-c(NewBurst)) %>%
  subset(!is.na(x) & !is.na(y))

# How many datapoints per dog and how many fixes per burst?
count(ssf, DogName)
count(ssf, BurstID)

# We can only work with bursts that contain at least three fixes
ssf <- ssf %>%
  group_by(DogName, BurstID) %>%
  nest() %>%
  mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
  subset(Nrow >= 3) %>%
  unnest(data) %>%
  ungroup() %>%
  dplyr::select(-c(Nrow, dt_)) %>%
  group_by(DogName, BurstID) %>%
  mutate(BurstID = cur_group_id()) %>%
  ungroup()

# How many datapoints per dog?
ssf_count <- count(ssf, DogName)

# Use amt to generate random steps
ssf <- ssf %>%
  make_track(.
    , .x      = x
    , .y      = y
    , .t      = Timestamp
    , id      = BurstID
    , crs     = CRS("+init=epsg:4326")
    , State   = State
    , DogName = DogName
  ) %>%
  transform_coords(CRS("+init=epsg:32734")) %>%
  nest(data = -"id") %>%
  mutate(data = map(data, function(x){
    x %>%
      steps(keep_cols = "start") %>%
      mutate(., dt_ = as.numeric(dt_, units = "hours"))
  })) %>%
  unnest(data) %>%
  subset(!is.na(ta_)) %>%
  random_steps(n_control = n_rand) %>%
  as.data.frame()

# Convert coordinates again
ssf[, c("x1_", "y1_")] <- reprojCoords(xy = as.matrix(ssf[, c("x1_", "y1_")])
  , from = "+init=epsg:32734"
  , to = "+init=epsg:4326"
)
ssf[, c("x2_", "y2_")] <- reprojCoords(xy = as.matrix(ssf[, c("x2_", "y2_")])
  , from = "+init=epsg:32734"
  , to = "+init=epsg:4326"
)

# Generate spatial lines that represent the different steps
ssf_lines <- pbmclapply(1:nrow(ssf), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  line <- spLines(
    rbind(
        c(ssf$x1_[x], ssf$y1_[x])
      , c(ssf$x2_[x], ssf$y2_[x])
    )
  )
  return(line)
}) %>% do.call(rbind, .)
crs(ssf_lines) <- CRS("+init=epsg:4326")
ssf_lines <- as(ssf_lines, "SpatialLinesDataFrame")
ssf_lines@data <- ssf
ssf <- ssf_lines

# Let's store those to file as well
writeVector(vect(ssf), "03_Data/02_CleanData/99_ValidationSteps.shp", overwrite = T)

# Combine data for path and step selection into a single tibble
dat <- tibble(
    ModelType     = c("PSF", "SSF")
  , DispersalData = list(psf, ssf)
)

################################################################################
#### Covariate Extraction
################################################################################
# Load the heatmaps
heat <- read_rds("03_Data/03_Results/99_Heatmaps.rds")

# Put together the buffer and the main areas
heat <- lapply(unique(heat$steps), function(x) {
  merged <- sum(stack(subset(heat, steps == x)$heatmap))
  merged <- tibble(Steps = x, Heatmap = list(merged))
  return(merged)
}) %>% do.call(rbind, .)

# Combine the maps with the dispersal data
dat <- expand_grid(dat, heat)

# Go through the different combinations and extract underlaying covariates
dat$DispersalData <- pbmclapply(1:nrow(dat), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  extracted <- extrCov(dat$Heatmap[[x]], dat$DispersalData[[x]])[, 1]
  extracted <- (extracted - mean(extracted)) / sd(extracted)
  newdat <- cbind(dat$DispersalData[[x]]@data, Heatmap = extracted)
  return(newdat)
})

# Let's run conditional logistic regression using the extracted data
dat$Model <- lapply(1:nrow(dat), function(x) {
  if (dat$ModelType[[x]] == "PSF") {
      formula <- case_ ~ Heatmap + strata(DogName)
      mod <- fit_clogit(dat$DispersalData[[x]], formula = formula)
      mod <- summary(mod)$coefficients
      mod <- as.data.frame(mod)
      mod <- mod[, c(1, 3, 5)]
      names(mod) <- c("coef", "coef_se", "pr>z")
    } else {
      formula <- case_ ~ Heatmap + strata(step_id_)
      dat_nested <- dat$DispersalData[[x]] %>% nest(Data = -DogName)
      dat_nested <- mutate(dat_nested, Model = map(Data, function(y) {
        mod <- fit_clogit(y, formula = formula)
        mod <- summary(mod)$coefficients
        mod <- as.data.frame(mod)
        mod <- mod[, c(1, 3, 5)]
        names(mod) <- c("coef", "coef_se", "pr>z")
        return(mod)
      }))
      mod <- dat_nested %>% dplyr::select(DogName, Model) %>% unnest(Model)
  }
  return(mod)
})

# Store the model results to file
dat <- dplyr::select(dat, -c(DispersalData, Heatmap))
write_rds(dat, "03_Data/03_Results/99_ConnectivityValidation.rds")

# Visualize results from psf
dat %>%
  subset(ModelType == "PSF") %>%
  dplyr::select(Steps, Model) %>%
  unnest(Model) %>%
  ggplot(aes(x = as.factor(Steps), y = coef, ymin = coef - 1.96 * coef_se, ymax = coef + 1.96 * coef_se)) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    geom_point() +
    geom_errorbar(width = 0.1) +
    theme_minimal() +
    ylim(c(-1, 5))

# Visualize results from ssf
dat %>%
  subset(ModelType == "SSF") %>%
  dplyr::select(Steps, Model) %>%
  unnest(Model) %>%
  ggplot(aes(x = as.factor(Steps), y = coef, ymin = coef - 1.96 * coef_se, ymax = coef + 1.96 * coef_se, col = DogName)) +
    geom_hline(yintercept = 0, lty = 2, col = "gray") +
    geom_point(position = position_dodge(width = 0.75)) +
    geom_errorbar(position = position_dodge(width = 0.75)) +
    theme_minimal() +
    scale_color_viridis_d()
