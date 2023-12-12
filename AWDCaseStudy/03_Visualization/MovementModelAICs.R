################################################################################
#### Plot of Model AICs
################################################################################
# Description: Plot of model AICs

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates nicely
library(ggpubr)       # For nice plots
library(davidoff)     # Custom functions
library(glmmTMB)      # To handle model results
library(lemon)        # For nice axis labels
library(viridis)      # For nice colors
library(ggdark)       # Dark ggplot theme
library(latex2exp)    # For latex expressions
library(msir)         # For prediction interval around loess
library(RColorBrewer) # For nice colors
library(xtable)       # To write model results to a .tex file

# Suppress the scientific notion
options(scipen = 999)

################################################################################
#### Prepare Data
################################################################################
# Load the model aics
aics <- readRDS("03_Data/03_Results/99_MovementModelAICs.rds")

# Also remove the ModelID column
aics$ModelID <- NULL

# Remove models with 0 weight
aics <- subset(aics, Weight > 0)

# Replace covariates with short code and add cos(ta) and log(sl) as covariatess
aics$Covariates <- aics$Covariates %>%
  gsub(pattern = "\\bShrubs, SqrtDistanceToWater, Trees, Water, HumansBuff5000\\b", replacement = "Base Model") %>%
  gsub(pattern = "\\bWater\\b", replacement = "WA") %>%
  gsub(pattern = "\\bShrubs\\b", replacement = "SH") %>%
  gsub(pattern = "\\bTrees\\b", replacement = "WO") %>%
  gsub(pattern = "\\bHumansBuff5000\\b", replacement = "HI") %>%
  gsub(pattern = "\\bSqrtDistanceToWater\\b", replacement = "DTW") %>%
  gsub(pattern = "\\inactive\\b", replacement = "LA") %>%
  gsub(pattern = "sl\\_", replacement = "sl") %>%
  gsub(pattern = "log\\_sl", replacement = "log(sl)") %>%
  gsub(pattern = "cos\\_ta\\_", replacement = "cos(ta)") %>%
  gsub(pattern = ",", replacement = " +")

# # Add new model IDs
# aics$ID <- 1:nrow(aics)

# Write the table to a .tex table
print(xtable(aics, digits = 2)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_MovementModelAICs.tex"
)

# We also want to visualize the models
models <- readRDS("03_Data/03_Results/99_MovementModel.rds")

# Assign new model IDs
models$ModelID <- 1:nrow(models)

# Remove unneccessary columns
models$Coeffs <- lapply(models$Model, function(x){
  getCoeffs(x, pvalue = T)[-1, ]
})

# Unnest
coeffs <- models %>%
  select(ModelID, Coeffs) %>%
  unnest(Coeffs)

# Calculate confidence intervals
coeffs <- coeffs %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05){
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1){
    return("*")
  } else {
    return(" ")
  }
})

# Rename covariates
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "cos_ta_"
  , replacement = "cos(ta)"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "log_sl_"
  , replacement = "log(sl)"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "sl_"
  , replacement = "sl"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "Shrubs"
  , replacement = "Shrubs/Grassland"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "Trees"
  , replacement = "Woodland"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "HumansBuff5000"
  , replacement = "HumanInfluence"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "inactiveTRUE"
  , replacement = "LowActivity"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "SqrtDistanceToWater"
  , replacement = "DistanceToWater '*m^0.5'"
)
coeffs$Preference <- ifelse(coeffs$Coefficient > 0, "Preferred", "Avoided")
coeffs$Preference <- factor(coeffs$Preference, levels = c("Preferred", "Avoided"))

# Create a dataframe showing which covariates are missing
missing <- coeffs %>%
  dplyr::select(ModelID, Covariate, Coefficient) %>%
  spread(key = Covariate, value = Coefficient) %>%
  mutate(across(!ModelID, function(x){!is.na(x)})) %>%
  gather(key = Covariate, value = Present, 2:ncol(.)) %>%
  subset(!Present) %>%
  mutate(Coefficient = 0)

# Specify the order in which the coefficients should be plotted
order <- c(
      "Water"
    , "DistanceToWater '*m^0.5'"
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "sl"
    , "cos(ta)"
    , "log(sl)"
    , "cos(ta):sl"
    , "cos(ta):log(sl)"
    , "sl:LowActivity"
    , "log(sl):LowActivity"
    , "sl:Water"
    , "sl:DistanceToWater '*m^0.5'"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , "sl:HumanInfluence"
    , "cos(ta):Water"
    , "cos(ta):DistanceToWater '*m^0.5'"
    , "cos(ta):Woodland"
    , "cos(ta):Shrubs/Grassland"
    , "cos(ta):HumanInfluence"
    , "log(sl):Water"
    , "log(sl):DistanceToWater '*m^0.5'"
    , "log(sl):Woodland"
    , "log(sl):Shrubs/Grassland"
    , "log(sl):HumanInfluence"
)

# Prepare groups with "expressions"
groups <- c(
      "Water"
    , expression(DistanceToWater^0.5)
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "sl"
    , "cos(ta)"
    , "log(sl)"
    , "cos(ta):sl"
    , "cos(ta):log(sl)"
    , "sl:LowActivity"
    , "log(sl):LowActivity"
    , "sl:Water"
    , expression(sl:DistanceToWater^0.5)
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , "sl:HumanInfluence"
    , "cos(ta):Water"
    , expression(cos(ta):DistanceToWater^0.5)
    , "cos(ta):Woodland"
    , "cos(ta):Shrubs/Grassland"
    , "cos(ta):HumanInfluence"
    , "log(sl):Water"
    , expression(log(sl):DistanceToWater^0.5)
    , "log(sl):Woodland"
    , "log(sl):Shrubs/Grassland"
    , "log(sl):HumanInfluence"
)

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_vline(
      xintercept = 0
    , color      = "darkgrey"
    , lty        = 2
    , lwd        = 0.3
  ) +
  geom_point(data = missing
    , col   = "gray30"
    , shape = 4
    , size  = 0.8
  ) +
  geom_point(
      shape = 3
    , size  = 2.5
  ) +
  geom_errorbarh(aes(
      xmin = LCI_90
    , xmax = UCI_90)
    , height = 0, size = 2, alpha = 0.5
  ) +
  geom_errorbarh(aes(
      xmin = LCI_95
    , xmax = UCI_95)
    , height = 0, size = 1, alpha = 0.75
  ) +
  geom_errorbarh(aes(
      xmin = LCI_99
    , xmax = UCI_99)
    , height = 0, size = 0.3, alpha = 1
  ) +
  geom_text(
      aes(label = Significance, hjust = 0.5, vjust = 0)
    , show.legend = F
    , size        = 2.5
  ) +
  scale_y_discrete(
      labels = rev(groups)
    , limits = rev(order)
  ) +
  theme_classic() +
  xlim(c(-1.5, 0.6)) +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(values = c("#5B9BD5", "orange")) +
  theme(
      legend.title      = element_blank()
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 0, 0, 0)
    , legend.box.margin = margin(-10, -10, -10, -10)
  ) +
  facet_wrap("ModelID", nrow = 3, scales = "free_x")

# Store the result to file
ggsave("test.png", plot = p1, width = 10, height = 12)

# Add annotations
cols <- brewer.pal(n = 9, name = "BuPu")
p2 <- p1 + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 22.75
    , ymax  = 27.25
    , alpha = 0.8
    , fill  = cols[6]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 15.75
    , ymax  = 22.25
    , alpha = 0.8
    , fill  = cols[5]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 1
    , ymax  = 15.25
    , alpha = 0.8
    , fill  = cols[4]
  )
  # ) + annotate(geom = "text"
  #   , x        = -1.25
  #   , y        = 23
  #   , label    = "HABITAT KERNEL"
  #   , color    = darken(cols[6])
  #   , angle    = 90
  #   , fontface = 3
  #   , size     = 2.5
  # ) + annotate(geom = "text"
  #   , x        = -1.25
  #   , y        = 18
  #   , label    = "MOVEMENT KERNEL"
  #   , color    = darken(cols[5])
  #   , angle    = 90
  #   , fontface = 3
  #   , size     = 2.5
  # ) + annotate(geom = "text"
  #   , x        = -1.25
  #   , y        = 7
  #   , label    = "INTERACTIONS"
  #   , color    = darken(cols[4])
  #   , angle    = 90
  #   , fontface = 3
  #   , size     = 2.5
  # )

# Add lines to separate kernels
p3 <- p2 + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 22.5
    , yend   = 22.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  ) + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 15.5
    , yend   = 15.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  )

# Store the plot
ggsave("04_Manuscript/99_MovementModelAICs.png", plot = p3, width = 10, height = 12)
