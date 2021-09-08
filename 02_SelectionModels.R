################################################################################
#### Step Selection Functions
################################################################################
# Description: Here, we combine the simulated movement data with the simulated
# covariates and apply step selection functions to estimate habitat and movement
# preferences.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(amt)           # For quick resource selection
library(broom)         # To clean a model
library(rgdal)         # To load shapefiles
library(sf)            # For nice plots
library(parallel)      # For multicore use
library(pbmcapply)     # For progress bar multicore use

# # Set working directory
# setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("00_Functions.R")

# Load observed movement data and covariate layers as well as the true
# preferences
dat   <- read_rds("99_ObservedMovements.rds")
cov   <- read_rds("99_CovariateLayers.rds")
nps   <- read_rds("99_NationalParks.rds")
truth <- read_rds("99_TruePreferences.rds")

# Make sure covariate layers are loaded into memory for quicker compuations on
# the layers.
cov <- readAll(cov)
inMemory(cov)

# How many GPS observations are there?
nrow(dat)

# How many observations per individual?
table(dat$ID)

# Visualize covariates and the collected GPS data
cov %>%
  as.data.frame(xy = T) %>%
  gather(key = covariate, value = value, 3:ncol(.)) %>%
  ggplot() +
    geom_raster(aes(x = x, y = y, fill = value)) +
    geom_point(data = dat, inherit.aes = F, aes(x = x, y = y, col = ID), size = 0.1) +
    geom_path(data = dat, inherit.aes = F, aes(x = x, y = y, col = ID), size = 0.1) +
    geom_sf(data = st_as_sf(nps), col = "white", fill = NA) +
    geom_sf_text(data = st_as_sf(nps), aes(label = ParkName), col = "white"
      , nudge_y = 10, size = 3) +
    facet_wrap("covariate") +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf() +
    theme(legend.position = "none")

################################################################################
#### Step Selection Function
################################################################################
# Nest the data by individual and coerce data to amt-tracks. We generate one
# track per individual.
dat_track <- dat %>%
  nest(data = -ID) %>%
  mutate(track = map(data, function(x){
    suppressMessages(
      make_track(x, .x = x, .y = y, .t = timestamp
      , step_number = step_number, step_id = step_id)
    )
  }))

# Let's see how this looks like.
print(dat_track)

# What's the sampling rate (in hours) of the GPS data?
sapply(dat_track$track, summarize_sampling_rate, time_unit = "hour")

# It looks like the data was collected at a resolution of 1 fix per hour. Hence,
# let's identify bursts where the fixes are not separated by more than one hour
# (+- 15 minutes of tolerance). We then also remove bursts with fewer than two
# consequcitve relocations
dat_track <- mutate(dat_track, track = map(track, function(x){
  temp <- track_resample(x, rate = hours(1), tolerance = minutes(15))
  temp <- filter_min_n_burst(temp, 2)
  return(temp)
}))

# We can now move from a point to a step representation. That is, we convert
# consecutive relocations into a line.
dat_track <- mutate(dat_track, track = map(track, function(x){
  steps_by_burst(x)
}))

# Unnest the steps so that we can fit distributions to step lengths and turning
# angles (we will use the same distributions for all individuals)
unnested <- dat_track %>%
  select(ID, track) %>%
  unnest(track) %>%
  select(ta_, sl_)

# Visualize distributions for turning angles and step lengths
hist(unnested$sl_, breaks = 30, main = "Step Length Distribution", xlab = "Step Length")
hist(unnested$ta_, breaks = 30, main = "Turning Angle Distribution", xlab = "Turning Angle")

# We can actually fit the scale and shape parameters of the gamma distribution
# and the concentration parameter of the vonmises distribution
sl_dist <- fit_distr(unnested$sl_, dist_name = "gamma")
ta_dist <- fit_distr(unnested$ta_, dist_name = "vonmises")

# Check the parameters of the fitted distributions. Note that the value of kappa
# is close to beta_cos_ta (= 1). This is exactly what Avgar et al. 2016 showed!
sl_dist
ta_dist

# However, Avgar et al. 2016 suggest to use a uniform distribution (i.e. a von
# Mises distribution with kappa = 0) for the turning angles when proposing
# random steps. Hence, we'll force kappa to 0
ta_dist$params$kappa <- 0

# We don't need the unnested data anymore and remove it
rm(unnested)

# We can now use the distributions to generate random steps. That is, for each
# "observed" or "realized" step, we now generate a set of 25 alternative steps.
# To generate said steps, we sample step lengths and turning angles from the
# specified distributions.
dat_track <- mutate(dat_track, ssf = map(track, function(x){
  random_steps(x
    , n_control = 25
    , sl_distr  = sl_dist
    , ta_distr  = ta_dist
  )
}))

# Unnest the steps (observed + random)
ssf <- dat_track %>%
  select(ID, ssf) %>%
  unnest(ssf)

# Compute the cosine of each relative turning angle -> cos_ta
ssf$cos_ta <- cos(ssf$ta_)

# Note that step_id_ is currently not unique across individuals as it was.
# Hence, we replace it with a unique identifier so that each stratum (a realized
# step + its 25 random steps) gets a unique identiier.
ssf <- ssf %>%
  group_by(ID, step_id_) %>%
  mutate(step_id = cur_group_id()) %>%
  ungroup() %>%
  select(-step_id_)

# Finally, we extract covariates along each step and calculate average values
# along the steps and bind the extracted data to the respective step
extracted <- extract_covariates_along_interpolated(ssf, cov, by = 0.1)
ssf <- cbind(ssf, extracted)

################################################################################
#### Estimate Beta: USING CLOGIT
################################################################################
# Now we can contrast realized and random steps using conditional logistic
# regression analysis. For simplicity, we do not consider random effects (in our
# case individuals have the same preferences anyways). However, if you're
# interested on how to fit a mixed effects conditional logistic regression,
# check Fieberg et al 2020
mod <- fit_clogit(case_ ~
  + elev
  + dist
  + cos_ta
  + strata(step_id)
  , data = ssf
)

# Check the summary and extract the coefficients
summary(mod)
beta_cl <- coef(mod)

# Visualize the model and add the true preferences for comparison
ggplot(tidy(mod$model), aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), height = 0.1) +
  geom_point(data = truth, inherit.aes = F, aes(x = Coefficient, y = Covariate), col = "red", pch = 8, size = 3) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 2)

# It appears that the model picks up the preferences correctly. Increasing the
# number of simulated individuals would of course improve estimates.
cbind(truth, beta_cl)

# Store the estimated preferences to a file.
results <- data.frame(
    Coefficient = names(beta_cl)
  , Estimate    = beta_cl
)
rownames(results) <- NULL
write_rds(results, "99_Estimates.rds")

# We also want to keep track of the step length and turn angle distirbutions
write_rds(sl_dist, "99_StepLengthDistribution.rds")
write_rds(ta_dist, "99_TurningAngleDistribution.rds")
