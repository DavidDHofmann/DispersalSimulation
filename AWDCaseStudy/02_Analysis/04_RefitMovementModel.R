################################################################################
#### Refit Movement Model
################################################################################
# Description: Here, I'll refit the movement model, yet without scaling movement
# covariates. I'll need to convert step lengths to kilometers for the model to
# converge though.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(glmmTMB)      # For modelling

################################################################################
#### Loading Data
################################################################################
# Load step selection data
dat <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")

# Keep only desired columns
dat <- dat %>% select(c(
      id
    , step_id_
    , case_
    , sl_
    , ta_
    , inactive
    , Water
    , DistanceToWater
    , Trees
    , Shrubs
    , HumansBuff5000 = HumanInfluenceBuffer_5000
))

# This is a crucial step. In the original model, step lengths were included in
# meters. This resulted in convergence issues, which I overcame by standardizing
# also the movement metrics. Alternatively, I could have also converted step
# lengths to kilometers, yet I realized that the model fit was better using
# standardized step lengths. Here, I will convert the step lengths to kilometers
# which will allow me to refit the full model and see if my
# "back-transformation" of the full model is correct.
dat$sl_ <- dat$sl_ / 1000

# We want to add the log of the step speed and the cosine of the turning angle
# and calculate the sqrt of the DistanceToWater
dat <- dat %>% mutate(
    log_sl_             = log(sl_)
  , cos_ta_             = cos(ta_)
  , SqrtDistanceToWater = sqrt(DistanceToWater)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id_, case_, sl_, log_sl_, ta_, cos_ta_, inactive, everything()
))

# Create scaled habitat covariates
dat <- transform(dat
  , DistanceToWater     = scale(DistanceToWater)
  , SqrtDistanceToWater = scale(SqrtDistanceToWater)
  , Water               = scale(Water)
  , Trees               = scale(Trees)
  , Shrubs              = scale(Shrubs)
  , HumansBuff5000      = scale(HumansBuff5000)
)

################################################################################
#### Refit Model
################################################################################
# Let's load our movement model
mov <- read_rds("03_Data/03_Results/99_MovementModel.rds")
mov <- mov$Model[[1]]

# Extract the model function
form <- mov$call$formula

# Refit the movement model using the unscaled data
mod <- glmm_clogit(form, data = dat)

# Show model results
summary(mod)

# Let's store the model
write_rds(mod, "03_Data/03_Results/99_MovementModelUnscaled.rds")
