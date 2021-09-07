################################################################################
#### Selection Functions
################################################################################
# Description: Use SSF analysis to estimate selection parameters using simulated
# data

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(amt)           # For quick resource selection
library(broom)         # To clean a model
library(lubridate)     # To handle dates
library(rgdal)         # To load shapefiles
library(sf)            # For nice plots
library(parallel)      # For multicore use
library(pbmcapply)     # For progress bar multicore use
library(jagsUI)        # To fit bayesian model

# Set working directory
setwd("/home/david/ownCloud/DispersalSimulation")

# Load custom functions
source("00_Functions.R")

# Load observed movement data and covariate layers as well as the true
# preferences
dat <- read_csv("ObservedMovements.csv")
cov <- stack("CovariateLayers.grd")
nps <- readOGR("NationalParks.shp")
truth <- read_csv("TruePreferences.csv")

# Convert IDs to factors
dat$ID <- as.factor(dat$ID)

# Make sure covariate layers are loaded into memory for maximal efficiency
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
    facet_wrap("covariate") +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    coord_sf() +
    theme(legend.position = "none")

################################################################################
#### Step Selection Function
################################################################################
# Nest the data by individual and coerce all data to amt-tracks
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

# What is the sampling rate at which GPS data was collected?
sapply(dat_track$track, summarize_sampling_rate, time_unit = "min")

# It looks like the data was collected at a resolution of 1 fix every hour.
# Hence, let's identify bursts where the fixes are not separated by more than
# one hour (+- tolerance). We then also remove bursts with fewer than two
# consequcitve relocations
dat_track <- mutate(dat_track, track = map(track, function(x){
  temp <- track_resample(x, rate = hours(1), tolerance = minutes(15))
  temp <- filter_min_n_burst(temp, 2)
  return(temp)
}))

# We can now move from a point to a step representation
dat_track <- mutate(dat_track, track = map(track, function(x){
  steps_by_burst(x)
}))

# Unnest the steps so that we can fit distributions to step lengths and turning
# angles (we will use the same distribution for all individuals)
unnested <- dat_track %>%
  select(ID, track) %>%
  unnest(track) %>%
  select(ta_, sl_)

# Visualize distributions for turning angles and step lengths
hist(unnested$sl_, breaks = 30)
hist(unnested$ta_, breaks = 30)

# We can actually fit the scale and shape parameters of the gamma distribution
# and the concentration parameter of the vonmises distribution
sl_dist <- fit_distr(unnested$sl_, dist_name = "gamma")
ta_dist <- fit_distr(unnested$ta_, dist_name = "vonmises")

# Check the parameters of the fitted distributions. Note that the value of kappa
# is equal to beta_cos_ta
sl_dist
ta_dist

# However, Avgar et al. 2016 suggest to use a uniform distribution (i.e. a von
# Mises distribution with kappa = 0) for the turning angles when proposing
# random steps. Hence, we'll force kappa to 0
ta_dist$params$kappa <- 0

# We don't need the unnested data anymore and remove it
rm(unnested)

# We can now use the distributions to generate random steps
dat_track <- mutate(dat_track, ssf = map(track, function(x){
  random_steps(x
    , n_control = 25
    , sl_distr  = sl_dist
    , ta_distr  = ta_dist
  )
}))

# Unnest the steps
ssf <- dat_track %>%
  select(ID, ssf) %>%
  unnest(ssf)

# Compute cos_ta
ssf$cos_ta <- cos(ssf$ta_)

# Note that step_id_ is not unique across individuals. Hence, we replace it with
# a unique identifier in each stratum
ssf <- ssf %>%
  group_by(ID, step_id_) %>%
  mutate(step_id = cur_group_id()) %>%
  ungroup() %>%
  select(-step_id_)

# Finally, we extract covariates along each step and calculate average values
# along the steps
extracted <- extract_covariates_along_interpolated(ssf, cov, by = 0.1)
ssf <- cbind(ssf, extracted)

################################################################################
#### Estimate Beta: USING CLOGIT
################################################################################
# Run conditional logistic regression model to estimate betas. For details on
# how to fit a mixed effects conditional logistic regression, check Fieberg et
# al 2020
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

# Visualize the model
ggplot(tidy(mod$model), aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), height = 0.1) +
  geom_point(data = truth, inherit.aes = F, aes(x = Coefficient, y = Covariate), col = "red", pch = 8, size = 3) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 2)

# ################################################################################
# #### Estimate Beta: MAXIMUM LIKELIHOOD BY HAND
# ################################################################################
# # Function to calculate the log-likelihoood
# loglik <- function(beta) {
#
#   # Identify unique strata
#   strata <- unique(ssf$step_id)
#
#   # Go through each strata and calculate log-likelihood for proposed beta (here
#   # we use parallel computing)
#   loglik <- mclapply(strata, mc.cores = detectCores() - 1, function(x){
#     sub <- ssf[ssf$step_id == x, ]
#     mat <- model.matrix(case_ ~ elev + dist + cos_ta, sub)[, -1]
#     score <- exp(mat %*% beta)
#     lik <- score[1] / sum(score)
#     log(lik)
#   })
#   loglik <- do.call(rbind, loglik)
#
#   # Return the summed log-likelihood
#   sum <- sum(loglik)
#   return(sum)
#
# }
#
# # Calculate log-likelihood for different betas
# betas <- expand.grid(beta1 = seq(0, 1, by = 0.1), beta2 = seq(-1, 0, by = 0.1), beta3 = seq(0, 1, by = 0.1))
# betas <- as.matrix(betas)
# logliks <- sapply(1:nrow(betas), function(x){loglik(beta = betas[x, ])})
#
# # Find maximum likelihood estimate
# beta_ml <- betas[which.max(logliks), ]
#
# # Compare to model
# cbind(ML_clogit = beta_cl, ML_hand = beta_ml)
#
# # Plot
# plot(logliks ~ betas[, 1], type = "o", xlab = "Beta", ylab = "Log-Likelihood", axes = F, pch = 16)
# abline(v = beta_ml, lty = 2, col = "red")
# abline(v = coef(mod), lty = 2, col = "blue")
# legend("bottomright", col = c("red", "blue"), legend = c("Brute ML", "Model ML"), lty = c(2, 2))
# axis(1)
# axis(2)
#
# ################################################################################
# #### Estimate Beta: MAXIMUM LIKELIHOOD USING MLE
# ################################################################################
# # To use the mle function, we need to be able to calculate the negative
# # log-likelihood instead of the log-likelihood
# negloglik <- function(beta){
#   loglik <- loglik(beta)
#   negloglik <- -loglik
#   return(negloglik)
# }
# beta_mle <- coef(mle(negloglik, start = list(beta = 0), method = "BFGS"))
#
# # Compare to model
# cbind(ML_clogit = beta_cl, ML_hand = beta_ml, ML_mle = beta_mle)

################################################################################
#### Estimate Beta: BAYESIAN
################################################################################
ssf
dat <- ssf[, c("step_id", "case_", "elev", "dist", "cos_ta")]
y <- ssf %>%
  dplyr::select(step_id, case_) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = case_) %>%
  dplyr::select(-step_id) %>%
  as.matrix()
elev <- ssf %>%
  dplyr::select(step_id, elev) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = elev) %>%
  dplyr::select(-step_id) %>%
  as.matrix()
dist <- ssf %>%
  dplyr::select(step_id, dist) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = dist) %>%
  dplyr::select(-step_id) %>%
  as.matrix()
cos_ta <- ssf %>%
  dplyr::select(step_id, cos_ta) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = cos_ta) %>%
  dplyr::select(-step_id) %>%
  as.matrix()

# Bundle data
dat_jags <- list(
    y           = ifelse(y, 1, 0)
  , elev        = elev
  , dist        = dist
  , cos_ta      = cos_ta
  , strata      = nrow(y)
  , strata_size = ncol(y)
)

# Write JAGS model file
cat(file = "model.txt", "model {

  # Prior
  beta_1 ~ dnorm(0, 0.001)
  beta_2 ~ dnorm(0, 0.001)
  beta_3 ~ dnorm(0, 0.001)

  # Likelihood
  for (i in 1:strata) {
    for (j in 1:strata_size) {
      phi[i, j] <- exp(beta_1 * elev[i, j] + beta_2 * dist[i, j] + beta_3 * dist[i, j])
    }
    for (j in 1:strata_size){
      lambda[i, j] <- phi[i, j] / sum(phi[i, 1:strata_size])
    }
    y[i, 1:strata_size] ~ dmulti(lambda[i, 1:strata_size], 1)
  }

}")

# Function to sample initial values
inits <- function(){
  list(
      beta_1 = rnorm(1)
    , beta_2 = rnorm(1)
    , beta_3 = rnorm(1)
  )
}

# Define parameters to be monitored (i.e. estimated)
params <- c("beta_1", "beta_2", "beta_3")

# MCMC settings (usually defined using trial and error)
na <- 2000        # Number of iterations in the adaptive phase
ni <- 2500        # Number of draws from the posterior (in each chain)
nb <- 1000        # Number of draws to discard as burn-in
nc <- 5           # Number of chains
nt <- 5           # Thinning rate (nt = 1 means we do not thin)

# Run the model
mod_jags <- jags(
    data               = dat_jags
  , inits              = inits
  , parameters.to.save = params
  , model.file         = "model.txt"
  , n.iter             = ni
  , n.burnin           = nb
  , n.chains           = nc
  , n.thin             = nt
  , n.adapt            = na
  , parallel           = T
)

# Show traceplots
par(mfrow = c(2, 2))
traceplot(mod_jags)

# Visualize distribution of beta_1
hist(mod_jags$sims.list$beta_1
  , breaks = 20
  , col    = "cornflowerblue"
  , main   = "Distribution of Beta 1"
  , xlab   = expression(beta)
  , xlim   = c(-1.5, 1.5)
)
abline(v = 0, lty = 2, col = "red")

# Visualize distribution of beta_2
hist(mod_jags$sims.list$beta_2
  , breaks = 20
  , col    = "cornflowerblue"
  , main   = "Distribution of Beta 2"
  , xlab   = expression(beta)
  , xlim   = c(-3, 3)
)
abline(v = 0, lty = 2, col = "red")

# Visualize distribution of beta_3
hist(mod_jags$sims.list$beta_3
  , breaks = 20
  , col    = "cornflowerblue"
  , main   = "Distribution of Beta 3"
  , xlab   = expression(beta)
  , xlim   = c(-3, 3)
)
abline(v = 0, lty = 2, col = "red")

# Summary of output, rounded to 3 digits
print(mod_jags, 3)

# Let's look at the estimates
beta_jags <- c(mod_jags$mean$beta_1, mod_jags$mean$beta_2)

################################################################################
#### Compare Results to Truth
################################################################################
# Compare estimates to true parameters
cbind(
    ML_clogit = beta_cl
  , ML_hand   = beta_ml
  , ML_mle    = beta_mle
  , JAGS      = beta_jags
)

# From now on, we'll work with the betas as estimated using conditional logistic
# regression. Thus, store them to a file.
results <- data.frame(
    Coefficient = names(beta_cl)
  , Estimate    = beta_cl
)
rownames(results) <- NULL
write_csv(results, "Estimates.csv")
