################################################################################
#### Movement Model
################################################################################
# Interpreting the movement model. I have fitted two different models. One were
# all covariates (habitat AND movement covariates) were scaled and one were only
# habitat covariates were scaled. In order for the second model to converge I
# needed to convert the step lengths to  kilometers. We'll use the first model
# to interpret the habitat kernel, yet the second model to investigate the
# movement kernel.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(glmmTMB)      # To handle glmm models
library(davidoff)     # Custom functions
library(amt)          # For updating distributions
library(lemon)        # For nice capped coords
library(viridis)      # For nice colors
library(ggpubr)       # To arrange plots
library(metR)         # To be able to get nice colorscales for factorial data

# Load original data used to fit the model
orig <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")
orig <- subset(orig, case_)

# Load our movement models (scaled and partly scaled with km steps)
model1 <- read_rds("03_Data/03_Results/99_MovementModel.rds")$Model[[1]]
model2 <- read_rds("03_Data/03_Results/99_MovementModelUnscaled.rds")

# Load the scaling parameters so that we can rescale stuff for visuals
scaling <- read_rds("03_Data/03_Results/99_Scaling.rds")
sc_cent <- scaling$center
sc_scal <- scaling$scale

# Load the tentative gamma distribution and adjust it to kilometers
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")
sl_dist$params$scale <- sl_dist$params$scale / 1000

# It may be worth pointing out that this is the same as the following
fit_distr(orig$sl_ / 1000, dist_name = "gamma")

# Prepare a turning angle distribution. We used a uniform distribution, which is
# similar to a vonmises distribution with kappa = 0
ta_dist <- list(
    name   = "vonmises"
  , params = list(kappa = 0)
)

# Extract data from the movement model and "unscale" the data. I'll refer to
# unscaled data as data_un
modeldat <- model1 %>%
  model.frame() %>%
  subset(case_) %>%
  dplyr::select(-c(case_, step_id_, id, inactive)) %>%
  mutate(
      sl_un                  = sl_ * sc_scal["sl_"] + sc_cent["sl_"]
    , log_sl_un              = log_sl_ * sc_scal["log_sl_"] + sc_cent["log_sl_"]
    , cos_ta_un              = cos_ta_ * sc_scal["cos_ta_"] + sc_cent["cos_ta_"]
    , Water_un               = Water * sc_scal["Water"] + sc_cent["Water"]
    , SqrtDistanceToWater_un = SqrtDistanceToWater * sc_scal["SqrtDistanceToWater"] + sc_cent["SqrtDistanceToWater"]
    , Shrubs_un              = Shrubs * sc_scal["Shrubs"] + sc_cent["Shrubs"]
    , Trees_un               = Trees * sc_scal["Trees"] + sc_cent["Trees"]
    , HumansBuff5000_un      = HumansBuff5000 * sc_scal["HumansBuff5000"] + sc_cent["HumansBuff5000"]
  ) %>%
  mutate_all(as.numeric) %>%
  dplyr::select(sort(names(.)))

# Identify the range on which we observed each covariate on the scaled and
# unscaled scale
ranges <- modeldat %>%
  gather(key = Covariate, value = Value) %>%
  group_by(Covariate) %>%
  summarize(
      Min    = min(Value)
    , Max    = max(Value)
    , Center = (Min + Max) / 2
    , Mean   = mean(Value)
    , Median = median(Value)
    , Q25    = quantile(Value, 0.25)
    , Q75    = quantile(Value, 0.75)
  )

################################################################################
#### Turning Angle vs Step Length
################################################################################
# Check out the movement model to see on which interactions the effect of
# turning angle depends
summary(model2)

# Extract model Coefficients from the partly scaled model in km
coefs <- fixef(model2)$cond

# Sequence for different step lengths (unscaled)
seq_sl_ <- seq(
    ranges$Min[ranges$Covariate == "sl_un"]
  , ranges$Max[ranges$Covariate == "sl_un"]
  , length.out = 100
) / 1000

# Show turning angle for different values of sl_
dat <- lapply(seq_sl_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        x +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"] +
      coefs["cos_ta_:HumansBuff5000"] *
        ranges$Center[ranges$Covariate == "HumansBuff5000"]
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the step length
  plot_ta$sl_ <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Backtransform the step lengths to km
dat$sl_ <- dat$sl_ * 1000

# Visualize using contour
p1a <- ggplot(dat, aes(x = ta_, y = sl_, z = prob)) +
  geom_contour_filled() +
  geom_rug(data = orig, aes(x = ta_)
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  geom_rug(data = orig, aes(y = sl_)
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  geom_contour(col = "black") +
  scale_fill_viridis(
      super  = metR::ScaleDiscretised
    , option = "magma"
    , name   = "Probability Density"
    , guide  = guide_colorsteps(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_y_continuous(
      breaks = seq(0, 35000, by = 5000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab("Step Length (m)") +
  theme(legend.position = "bottom")

# Visualize using x-y plot
p1b <- ggplot(dat, aes(x = ta_, y = prob, group = sl_, color = sl_)) +
  geom_line(size = 1) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = "Step Length (m)"
    , labels = function(x){format(x, big.mark = "'")}
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_y_continuous(
      breaks = seq(0, 0.3, by = 0.1)
    , labels = sprintf("%0.1f", seq(0, 0.3, by = 0.1))
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plots
p1a
p1b

################################################################################
#### Turning Angle vs SqrtDistanceToWater
################################################################################
# Sequence for different distances to water
seq_wat_ <- seq(
    ranges$Min[ranges$Covariate == "SqrtDistanceToWater"]
  , ranges$Max[ranges$Covariate == "SqrtDistanceToWater"]
  , length.out = 100
)

# Show turning angle for different values of sl_
dat <- lapply(seq_wat_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        ranges$Center[ranges$Covariate == "sl_un"] / 1000 +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        x +
      coefs["cos_ta_:HumansBuff5000"] *
        ranges$Center[ranges$Covariate == "HumansBuff5000"]
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the water cover
  plot_ta$SqrtDistanceToWater <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Backtransform the covariates
dat$SqrtDistanceToWater <- dat$SqrtDistanceToWater *
  sc_scal["SqrtDistanceToWater"] + sc_cent["SqrtDistanceToWater"]

# Visualize using contour
p2a <- ggplot(dat, aes(x = ta_, y = SqrtDistanceToWater, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  geom_rug(data = orig, aes(x = ta_)
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  geom_rug(data = orig, aes(y = sqrt(DistanceToWater))
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  scale_fill_viridis(
      super  = metR::ScaleDiscretised
    , option = "magma"
    , name   = "Probability Density"
    , guide  = guide_colorsteps(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_y_continuous(
      labels = function(x){format(x, big.mark = "'")}
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab(bquote(DistanceToWater^0.5 * .(" (m)"))) +
  theme(legend.position = "bottom")

# Visualize using x-y plot
p2b <- ggplot(dat, aes(x = ta_, y = prob, group = SqrtDistanceToWater, color = SqrtDistanceToWater)) +
  geom_line(size = 1) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = bquote(DistanceToWater^0.5 * .(" (m)"))
    , labels = function(x){format(x, big.mark = "'")}
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plots
p2a
p2b

################################################################################
#### Turning Angle vs Humans
################################################################################
# Sequence for different human influences
seq_hum_ <- seq(
    ranges$Min[ranges$Covariate == "HumansBuff5000"]
  , ranges$Max[ranges$Covariate == "HumansBuff5000"]
  , length.out = 100
)

# Show turning angle for different values of sl_
dat <- lapply(seq_hum_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        ranges$Center[ranges$Covariate == "sl_un"] / 1000 +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"] +
      coefs["cos_ta_:HumansBuff5000"] *
        x
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the water cover
  plot_ta$HumansBuff5000 <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Backtransform the covariates
dat$HumansBuff5000 <- dat$HumansBuff5000 *
  sc_scal["HumansBuff5000"] + sc_cent["HumansBuff5000"]

# Visualize using contour
p3a <- ggplot(dat, aes(x = ta_, y = HumansBuff5000, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  geom_rug(data = orig, aes(x = ta_)
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  geom_rug(data = orig, aes(y = HumanInfluenceBuffer_5000)
    , inherit.aes = F
    , size        = 0.3
    , alpha       = 0.3
  ) +
  scale_fill_viridis(
      super  = metR::ScaleDiscretised
    , option = "magma"
    , name   = "Probability Density"
    , guide  = guide_colorsteps(
        show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 10)
  ) +
  scale_y_continuous(
      breaks = seq(0, 10, by = 2)
    , labels = seq(0, 10, by = 2)
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab("Human Influence") +
  theme(legend.position = "bottom")

# Visualize using x-y plot
p3b <- ggplot(dat, aes(x = ta_, y = prob, group = HumansBuff5000, color = HumansBuff5000)) +
  geom_line(size = 1) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = "Human Influence"
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plots
p3a
p3b

################################################################################
#### Step Length vs Water
################################################################################
# Check out model to see on which variables the effect of step lengths depends
summary(model2)

# Sequence for different distances to water
seq_wat_ <- seq(
    ranges$Min[ranges$Covariate == "Water"]
  , ranges$Max[ranges$Covariate == "Water"]
  , length.out = 100
)

# Show sl_ for different values of water
dat <- lapply(seq_wat_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        0 +
      coefs["sl_:Water"] *
        x +
      coefs["sl_:Trees"] *
        ranges$Center[ranges$Covariate == "Trees"] +
      coefs["sl_:Shrubs"] *
        ranges$Center[ranges$Covariate == "Shrubs"] +
      coefs["sl_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"]
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_sl <- data.frame(sl_ = seq(from = 0, to = 35, length.out = 1000))

  # Insert the distance to water
  plot_sl$Water <- x

  # Get probabilities from updated distribution
  plot_sl$prob <- dgamma(plot_sl$sl_
    , scale = updated$params$scale
    , shape = updated$params$shape
  )

  # Return the final data
  return(plot_sl)

}) %>% do.call(rbind, .)

# Backtransform covariates
dat$sl_ <- dat$sl_ * 1000
dat$Water <- dat$Water *
  sc_scal["Water"] + sc_cent["Water"]

# Visualize using x-y plot
p4 <- ggplot(dat, aes(x = sl_, y = prob, group = Water, color = Water)) +
  geom_line(size = 0.2) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
    , xlim   = c(0, 35000)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = "Water Cover (%)"
    , limits = c(0, 1)
    , breaks = seq(0, 1, by = 0.2)
    , labels = seq(0, 1, by = 0.2)
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = seq(0, 35000, 10000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plot
p4

################################################################################
#### Step Length vs SqrtDistanceToWater
################################################################################
# Sequence for different distances to water
seq_wat_ <- seq(
    ranges$Min[ranges$Covariate == "SqrtDistanceToWater"]
  , ranges$Max[ranges$Covariate == "SqrtDistanceToWater"]
  , length.out = 100
)

# Show sl_ for different values of water
dat <- lapply(seq_wat_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        0 +
      coefs["sl_:Water"] *
        ranges$Center[ranges$Covariate == "Water"] +
      coefs["sl_:Trees"] *
        ranges$Center[ranges$Covariate == "Trees"] +
      coefs["sl_:Shrubs"] *
        ranges$Center[ranges$Covariate == "Shrubs"] +
      coefs["sl_:SqrtDistanceToWater"] *
        x
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_sl <- data.frame(sl_ = seq(from = 0, to = 35, length.out = 1000))

  # Insert the distance to water
  plot_sl$SqrtDistanceToWater <- x

  # Get probabilities from updated distribution
  plot_sl$prob <- dgamma(plot_sl$sl_
    , scale = updated$params$scale
    , shape = updated$params$shape
  )

  # Return the final data
  return(plot_sl)

}) %>% do.call(rbind, .)

# Backtransform covariates
dat$sl_ <- dat$sl_ * 1000
dat$SqrtDistanceToWater <- dat$SqrtDistanceToWater *
  sc_scal["SqrtDistanceToWater"] + sc_cent["SqrtDistanceToWater"]

# Visualize using x-y plot
p5 <- ggplot(dat, aes(x = sl_, y = prob, group = SqrtDistanceToWater, color = SqrtDistanceToWater)) +
  geom_line(size = 0.2) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
    , xlim   = c(0, 35000)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = bquote(DistanceToWater^0.5 * .(" (m)"))
    , limits = c(0, 200)
    , breaks = seq(0, 200, by = 50)
    , labels = seq(0, 200, by = 50)
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = seq(0, 35000, 10000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plot
p5

################################################################################
#### Step Length vs Shrubs
################################################################################
# Sequence for different shrub cover
seq_shrub_ <- seq(
    ranges$Min[ranges$Covariate == "Shrubs"]
  , ranges$Max[ranges$Covariate == "Shrubs"]
  , length.out = 100
)

# Show sl_ for different values of water
dat <- lapply(seq_shrub_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        0 +
      coefs["sl_:Water"] *
        ranges$Center[ranges$Covariate == "Water"] +
      coefs["sl_:Trees"] *
        ranges$Center[ranges$Covariate == "Trees"] +
      coefs["sl_:Shrubs"] *
        x +
      coefs["sl_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"]
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_sl <- data.frame(sl_ = seq(from = 1, to = 35, length.out = 1000))

  # Insert the distance to water
  plot_sl$Shrubs <- x

  # Get probabilities from updated distribution
  plot_sl$prob <- dgamma(plot_sl$sl_
    , scale = updated$params$scale
    , shape = updated$params$shape
  )

  # Return the final data
  return(plot_sl)

}) %>% do.call(rbind, .)

# Backtransform covariates
dat$sl_ <- dat$sl_ * 1000
dat$Shrubs <- dat$Shrubs *
  sc_scal["Shrubs"] + sc_cent["Shrubs"]

# Visualize using x-y plot
p6 <- ggplot(dat, aes(x = sl_, y = prob, group = Shrubs, color = Shrubs)) +
  geom_line(size = 0.2) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
    , xlim   = c(0, 35000)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = "Shrub/Grassland Cover (%)"
    , limits = c(0, 1)
    , breaks = seq(0, 1, by = 0.2)
    , labels = seq(0, 1, by = 0.2)
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = seq(0, 35000, 10000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plot
p6

################################################################################
#### Step Length vs Trees
################################################################################
# Sequence for different shrub cover
seq_tree_ <- seq(
    ranges$Min[ranges$Covariate == "Trees"]
  , ranges$Max[ranges$Covariate == "Trees"]
  , length.out = 100
)

# Show sl_ for different values of water
dat <- lapply(seq_tree_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        0 +
      coefs["sl_:Water"] *
        ranges$Center[ranges$Covariate == "Water"] +
      coefs["sl_:Trees"] *
        x +
      coefs["sl_:Shrubs"] *
        ranges$Center[ranges$Covariate == "Shrubs"] +
      coefs["sl_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"]
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_sl <- data.frame(sl_ = seq(from = 1, to = 35, length.out = 1000))

  # Insert the distance to water
  plot_sl$Trees <- x

  # Get probabilities from updated distribution
  plot_sl$prob <- dgamma(plot_sl$sl_
    , scale = updated$params$scale
    , shape = updated$params$shape
  )

  # Return the final data
  return(plot_sl)

}) %>% do.call(rbind, .)

# Backtransform covariates
dat$sl_ <- dat$sl_ * 1000
dat$Trees <- dat$Trees *
  sc_scal["Trees"] + sc_cent["Trees"]

# Visualize using x-y plot
p7 <- ggplot(dat, aes(x = sl_, y = prob, group = Trees, color = Trees)) +
  geom_line(size = 0.2) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.3)
    , xlim   = c(0, 35000)
  ) +
  scale_color_viridis_c(
      begin  = 0.2
    , end    = 0.9
    , option = "magma"
    , name   = "Woodland Cover (%)"
    , limits = c(0, 0.3)
    , breaks = seq(0, 0.3, by = 0.05)
    , labels = seq(0, 0.3, by = 0.05)
    , guide  = guide_colorbar(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_x_continuous(
      breaks = seq(0, 35000, 10000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme(legend.position = "bottom")

# Show plot
p7

################################################################################
#### Step Length: Inactive vs Active
################################################################################
# Sequence for activity
seq_active_ <- c(0, 1)

# Show sl_ for different values of water
dat <- lapply(seq_active_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        x +
      coefs["sl_:Water"] *
        ranges$Center[ranges$Covariate == "Water"] +
      coefs["sl_:Trees"] *
        ranges$Center[ranges$Covariate == "Trees"] +
      coefs["sl_:Shrubs"] *
        ranges$Center[ranges$Covariate == "Shrubs"] +
      coefs["sl_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"]
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_sl <- data.frame(sl_ = seq(from = 1, to = 35, length.out = 1000))

  # Insert the distance to water
  plot_sl$Inactive <- x

  # Get probabilities from updated distribution
  plot_sl$prob <- dgamma(plot_sl$sl_
    , scale = updated$params$scale
    , shape = updated$params$shape
  )

  # Return the final data
  return(plot_sl)

}) %>% do.call(rbind, .)

# Backtransform covariates
dat$sl_ <- dat$sl_ * 1000

# Make "Active" vs. "Inactive" factor
dat$Activity <- ifelse(dat$Inactive, "Inactive", "Active")

# Visualize using x-y plot
p8 <- ggplot(dat, aes(x = sl_, y = prob, color = factor(Activity))) +
  geom_line(size = 1) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , ylim   = c(0, 0.2)
    , xlim   = c(0, 35000)
  ) +
  scale_color_viridis_d(
      begin  = 0.2
    , end    = 0.8
    , option = "magma"
    , name   = "Activity Period"
    , guide  = guide_legend(
        title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.3, "cm")
      , barwidth       = unit(5.0, "cm")
    )
  ) +
  scale_y_continuous(
      breaks = seq(0, 0.3, by = 0.1)
    , labels = sprintf("%0.1f", seq(0, 0.3, by = 0.1))
  ) +
  scale_x_continuous(
      breaks = seq(0, 35000, 10000)
    , labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme(
      legend.position  = "bottom"
    , legend.key.width = unit(1, "cm")
  )

# Show plot
p8

################################################################################
#### Function to Predict the RSS
################################################################################
# Let's write a function that allows us to predict the relative selection score
predictRSS <- function(model, df1, df2, ci = NULL, return_data = F){

  # Use custom function to prepare model
  mod <- prepareModel(model)

  # Span model matrices
  df1_mm <- as.data.frame(model.matrix(mod$formula, df1))
  df2_mm <- as.data.frame(model.matrix(mod$formula, df2))

  # Remove intercept
  df1_mm <- df1_mm[, names(df1_mm) != "(Intercept)"]
  df2_mm <- df2_mm[, names(df2_mm) != "(Intercept)"]

  # Predict scores
  s1 <- predictScore(
      coefficients = mod$coefficients
    , formula      = mod$formula
    , data         = df1
  )
  s2 <- predictScore(
      coefficients = mod$coefficients
    , formula      = mod$formula
    , data         = df2
  )

  # Calculate relative selection scores
  rss <- s1 / s2

  # Calculate confidence interval if desired
  if (!is.null(ci)){

    # Get model variance matrix
    vars <- vcov(model)[[1]]

    # Remove intercept
    vars <- vars[rownames(vars) != "(Intercept)", colnames(vars) != "(Intercept)"]

    # Get difference matrix between df1 and df2
    diff <- sweep(as.matrix(df1_mm), 2, as.matrix(df2_mm))

    # Get variance of log-rss prediction
    varpred <- diag(diff %*% vars %*% t(diff))

    # Get the standard error
    logrss_se <- unname(sqrt(varpred))

    # Prepare table to store confidence intervals
    intervals <- lapply(1:length(ci), function(x){

      # Get critical value
      p <- 1 - ((1 - ci[x]) / 2)
      zstar <- qnorm(p)

      # Compute intervals
      lwr <- exp(log(rss) - zstar * logrss_se)
      upr <- exp(log(rss) + zstar * logrss_se)

      # Prepare dataframe
      cis <- as.data.frame(cbind(lwr, upr))
      names(cis) <- paste0(c("Lower_", "Upper_"), 100 * ci[x])

      # Return them
      return(cis)

    }) %>% do.call(cbind, .)

    rss <- data.frame(RSS = rss, intervals)

  }

  # Return data if desired
  if (return_data){
    rss <- cbind(df1, rss)
  }

  # Return the final data
  return(rss)

}

# df1 <- expand_grid(
#     sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
#   , cos_ta_             = seq(-2, 2, length.out = 100)
#   , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
#   , Shrubs              = 0
#   , Water               = 0
#   , SqrtDistanceToWater = 0
#   , Trees               = 0
#   , HumansBuff5000      = seq(-2, 2, length.out = 100)
#   , inactive            = F
# )
# df2 <- expand_grid(
#     sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
#   , cos_ta_             = 0
#   , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
#   , Shrubs              = 0
#   , Water               = 0
#   , SqrtDistanceToWater = 0
#   , Trees               = 0
#   , HumansBuff5000      = 0
#   , inactive            = F
# )
# library(raster)
# library(viridis)
# library(rasterVis)
# test <- predictRSS(model, df1, df2, ci = NULL, return_data = T)
# test <- test[, c("HumansBuff5000", "cos_ta_", "rss")]
# test$HumansBuff5000 <- test$HumansBuff5000 * scaling$scale["HumansBuff5000"] + scaling$center["HumansBuff5000"]
# test$cos_ta_ <- test$cos_ta_ * scaling$scale["cos_ta_"] + scaling$center["cos_ta_"]
# test <- rasterFromXYZ(test)
# plot(test, col = viridis(50), xlab = "HumanBuff5000", ylab = "cos_ta_")
# levelplot(test, xlab = "HumanBuff5000", ylab = "cos_ta_")
# contour(test, add = T)
# as.matrix(spread(test, cos_ta_, rss))

# Create a reference dataframe (all covariates need to be scaled)
ref <- data.frame(
    sl_                 = scale(2000, center = sc_cent["sl_"], scale = sc_scal["sl_"])
  , cos_ta_             = 0
  , log_sl_             = scale(log(2000), center = sc_cent["log_sl_"], scale = sc_scal["log_sl_"])
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Function to plot RSS against a covariate
showRSS <- function(
    model     = NULL
  , refdat    = NULL
  , covariate = NULL
  , ci        = c(0.99, 0.95, 0.9)
  , values    = NULL){

  # Create two dataframes. One as reference
  df1 <- refdat
  df2 <- refdat

  # Replace the values of the covariate for which we check the RSS
  df1 <- df1[rep(1, length(values)), ]
  df1[, c(covariate)] <- values

  # Predict scores
  pred <- predictRSS(
      model       = model
    , df1         = df1
    , df2         = df2
    , ci          = ci
    , return_data = F
  )

  # Add the covariate
  pred[, covariate] <- values

  # Backtransform the covariate
  pred[, covariate] <- pred[, covariate] * sc_scal[covariate] + sc_cent[covariate]
  pred$Covariate <- pred[, covariate]

  # Make tidy
  pred <- pred %>%
    gather(key = Interval, value = Boundary, c(contains("Lower"), contains("Upper"))) %>%
    separate(Interval, into = c("Type", "Level"), sep = "_") %>%
    spread(key = Type, value = Boundary)

  # Visualize
  ggplot(pred, aes(x = Covariate, y = RSS)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
      , linetype = "solid"
      , alpha    = 0.33
      , color    = "orange"
      , fill     = "orange"
      , lwd      = 0.1
    ) +
    geom_line(size = 1) +
    xlab(covariate) +
    ylab(paste0("RSS vs. ", covariate)) +
    theme_classic()

}

# Span vectors for all variables that we want to check
vars <- c("Water", "SqrtDistanceToWater", "Shrubs", "Trees", "HumansBuff5000")
vars <- lapply(vars, function(x){
  df <- seq(
      ranges$Min[ranges$Covariate == x]
    , ranges$Max[ranges$Covariate == x]
    , length.out = 1000
  ) %>% as.data.frame() %>% setNames(x)
  return(df)
}) %>% do.call(cbind, .)

# Visualize the RSS
p9 <-   showRSS(model = model1, refdat = ref
  , covariate = "Water", values = vars$Water)
p10 <-  showRSS(model = model1, refdat = ref
  , covariate = "SqrtDistanceToWater", values = vars$SqrtDistanceToWater)
p11 <-  showRSS(model = model1, refdat = ref
  , covariate = "Shrubs", values = vars$Shrubs)
p12 <-  showRSS(model = model1, refdat = ref
  , covariate = "Trees", values = vars$Trees)
p13 <-  showRSS(model = model1, refdat = ref
  , covariate = "HumansBuff5000", values = vars$HumansBuff5000)

# Adjust names of y and x axis
p9 <- p9 + xlab("Water Cover (%)") + ylab("RSS vs. Water Cover")
p10 <- p10 + xlab(expression(DistanceToWater^0.5)) + ylab(bquote(.("RSS vs.") ~ DistanceToWater^0.5))
p11 <- p11 + xlab("Shrub/Grassland Cover (%)") + ylab("RSS vs. Shrub/Grassland Cover")
p12 <- p12 + xlab("Woodland Cover (%)") + ylab("RSS vs. Woodland Cover")
p13 <- p13 + xlab("Human Influence") + ylab("RSS vs. Human Influence")

# Prepare a separate legend
legend <- data.frame(x = 0:100, y = 0:100)
legend$Level <- sample(c("99%", "95%", "90%"), size = nrow(legend), replace = T)
legend <- ggplot(legend, aes(x = x, y = y)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = y - 20, ymax = y + 20, fill = Level)
    , linetype = "solid"
    , color    = "orange"
    , lwd      = 0.1
  ) +
  geom_line(size = 1) +
  scale_fill_manual(name = "Confidence Level", values = c(
      adjustcolor("orange", alpha.f = 0.33)
    , adjustcolor("orange", alpha.f = 0.5)
    , adjustcolor("orange", alpha.f = 0.7)
    )
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(
      legend.position   = "bottom"
    , legend.text       = element_text(face = 3)
    , legend.title      = element_text(face = 3)
  )
legend <- get_legend(legend)

################################################################################
#### Putting Plots Together
################################################################################
# Put plots for movement and habitat kernel together
p_mov <- ggarrange(p1b, p2b, p3b, p4, p5, p6, p7, p8
  , ncol   = 4
  , nrow   = 2
  , labels = paste0("a", 1:8)
  , vjust  = -0.4
)
p_hab <- ggarrange(p9, p10, p11, p12, p13, legend
  , ncol   = 4
  , nrow   = 2
  , labels = paste0("b", 1:5)
  , vjust  = -0.4
)

# Specify plot titles
title_mov <- expression(atop(bold("Movement Kernel + Interactions")))
title_hab <- expression(atop(bold("Habitat Kernel")))

# Add to plots
p_mov <- annotate_figure(p_mov, top = text_grob(title_mov, size = 15))
p_hab <- annotate_figure(p_hab, top = text_grob(title_hab, size = 15))

# Put plots together
p <- ggarrange(p_mov, p_hab, nrow = 2)

# Store them
ggsave("04_Manuscript/99_MovementModelInterpretation.png", plot = p, scale = 1.7, bg = "white")
