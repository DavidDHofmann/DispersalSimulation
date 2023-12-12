################################################################################
#### Plot Movement Model Results
################################################################################
# Description: Plot of most parsimonious movement model and results from the
# model validation procedure

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
# Load the movement model
best <- readRDS("03_Data/03_Results/99_MovementModel.rds")
best <- best$Model[[1]]

# Extract model results (note that newer versions of glmmTMB and TMB are not
# compatible with older models and we need to manually extract the values of
# interest)
best <- tibble(
    Covariate   = names(fixef(best)$cond)
  , Coefficient = fixef(best)$cond
  , SE          = sqrt(diag(vcov(best, full = T)))[1:length(Coefficient)]
)

# Calculate confidence intervals
coeffs <- best[-1, ] %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Compute p-values
coeffs <- mutate(coeffs, pvalue = 2 * pnorm(-abs(Coefficient) / SE))

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05) {
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1) {
    return("*")
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
coeffs$Preference <- factor(coeffs$Preference, levels = c("Avoided", "Preferred"))

# Create a tidy dataframe for the confidence intervals
confs <- coeffs %>%
  dplyr::select(Covariate, Coefficient, LCI_90, UCI_90, LCI_95, UCI_95, LCI_99, UCI_99) %>%
  pivot_longer(!Covariate & !Coefficient) %>%
  separate(name, into = c("Limit", "Level"), sep = "_", convert = T) %>%
  spread(key = Limit, value = value) %>%
  mutate(alpha = case_when(
      Level == 90 ~ 0.50
    , Level == 95 ~ 0.75
    , Level == 99 ~ 1.00
  )) %>%
  mutate(size = case_when(
      Level == 90 ~ 2.0
    , Level == 95 ~ 1.0
    , Level == 99 ~ 0.3
  )) %>%
  mutate(Level = as.factor(Level))

# Join the information on preference
confs <- coeffs %>%
  dplyr::select(Covariate, Preference) %>%
  left_join(confs, ., by = "Covariate")

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
    , "sl:Water"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , "sl:DistanceToWater '*m^0.5'"
    , "cos(ta):HumanInfluence"
    , "cos(ta):DistanceToWater '*m^0.5'"
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
    , "sl:Water"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , expression(sl:DistanceToWater^0.5)
    , "cos(ta):HumanInfluence"
    , expression(cos(ta):DistanceToWater^0.5)
)

# Prepare dataset for plotting confidence intervals
coeffs2 <- coeffs %>%
  dplyr::select(Covariate, Coefficient, Preference, LCI_90:UCI_99) %>%
  gather(key = confidence_level, value = value, LCI_90:UCI_99) %>%
  separate(col = confidence_level, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = value) %>%
  mutate(Level = paste0(Level, "%"))

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 3, size = 2.5) +
  geom_errorbarh(
      aes(
        xmin = LCI
      , xmax = UCI
      , size = factor(Level)
    )
    , data = coeffs2
    , height = 0
    , alpha  = 0.5
  ) +
  geom_text(
      aes(label = Significance, hjust = 0.5, vjust = 0)
    , show.legend = F
  ) +
  geom_vline(
      xintercept = 0
    , color      = "darkgrey"
    , lty        = 2
    , lwd        = 0.3
  ) +
  scale_y_discrete(
      labels = rev(groups)
    , limits = rev(order)
  ) +
  theme_classic() +
  xlim(c(-1.3, 0.6)) +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(
      name   = "Preference"
    , values = c("#5B9BD5", "orange")
  ) +
  scale_size_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.3)
  ) +
  theme(
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 50, 0, -20)
    , legend.box.margin = margin(-5, -10, -5, -10)
    , legend.text       = element_text(face = 3)
    , legend.title      = element_text(face = 3)
  ) +
  guides(
      colour = guide_legend(title.position = "top", title.hjust = 0.5)
    , size   = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# # Legacy Plot
# p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
#   geom_point(shape = 3, size = 2.5) +
#   geom_errorbarh(aes(
#       xmin = LCI_90
#     , xmax = UCI_90)
#     , height = 0, size = 2, alpha = 0.5
#   ) +
#   geom_errorbarh(aes(
#       xmin = LCI_95
#     , xmax = UCI_95)
#     , height = 0, size = 1, alpha = 0.75
#   ) +
#   geom_errorbarh(aes(
#       xmin = LCI_99
#     , xmax = UCI_99)
#     , height = 0, size = 0.3, alpha = 1
#   ) +
#   geom_text(
#       aes(label = Significance, hjust = 0.5, vjust = 0)
#     , show.legend = F
#   ) +
#   geom_vline(
#       xintercept = 0
#     , color      = "darkgrey"
#     , lty        = 2
#     , lwd        = 0.3
#   ) +
#   scale_y_discrete(
#       labels = rev(groups)
#     , limits = rev(order)
#   ) +
#   theme_classic() +
#   xlim(c(-1.3, 0.6)) +
#   coord_capped_cart(
#       left   = "both"
#     , bottom = "both"
#   ) +
#   labs(x = expression(beta*"-Coefficient")) +
#   scale_color_manual(values = c("#5B9BD5", "orange")) +
#   theme(
#       legend.title      = element_blank()
#     , legend.position   = "bottom"
#     , legend.margin     = margin(0, 0, 0, 0)
#     , legend.box.margin = margin(-10, -10, -10, -10)
#   )

# Add annotations
cols <- brewer.pal(n = 9, name = "BuPu")
p2 <- p1 + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 12.75
    , ymax  = 17
    , alpha = 0.8
    , fill  = cols[6]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 6.75
    , ymax  = 12.25
    , alpha = 0.8
    , fill  = cols[5]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 1
    , ymax  = 6.25
    , alpha = 0.8
    , fill  = cols[4]
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 15
    , label    = "HABITAT KERNEL"
    , color    = darken(cols[6])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 9.5
    , label    = "MOVEMENT KERNEL"
    , color    = darken(cols[5])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 3.5
    , label    = "INTERACTIONS"
    , color    = darken(cols[4])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  )

# Add lines to separate kernels
p3 <- p2 + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 6.5
    , yend   = 6.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  ) + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 12.5
    , yend   = 12.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  )

# Show the plot
p3

# # Store the plot
# ggsave("04_Manuscript/99_MovementModel.png", device = "png", width = 7, height = 7, scale = 0.75)
# ggsave("04_Manuscript/99_MovementModel.pdf", device = "pdf", width = 7, height = 7, scale = 0.75)

# Let's also prepare a table to include in the manuscript
tex <- left_join(data.frame(Covariate = order), coeffs, by = "Covariate") %>%
  mutate(Group = rep(c("Habitat Kernel", "Movement Kernel", "Interaction"), c(5, 6, 6))) %>%
  select(Group, Covariate, Coefficient, SE, pvalue, Significance) %>%
  mutate(across(c(Coefficient, SE, pvalue), round, 3))

# Write it to a .tex file
print(xtable(tex, digits = 3)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "04_Manuscript/99_MovementModel.tex"
)

################################################################################
#### Validation
################################################################################
# Load required data
validation     <- read_rds("03_Data/03_Results/99_ModelValidation.rds")
dat_pref       <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[1]]
dat_rand       <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[2]]
validation_con <- read_rds("03_Data/03_Results/99_ConnectivityValidation.rds")

# Put the k-fold cross validation data from observed and random preferences
# together
dat_pref$Group <- "Realized"
dat_rand$Group <- "Random"
dat <- rbind(dat_pref, dat_rand)

# We want to plot this data and add the information from the validation table
# too. Let's prepare a column that indicates the text that we want to plot on
# top of the data. Let's first round the values from the validation table
validation[, c(2:4)] <- round(validation[, c(2:4)], 2)

# Match the Grouping names of the validation table to the groups above
validation$Group <- as.factor(c("Realized", "Random"))

# Create a dataframe which we use to annotate the two facets
text <- data.frame(
    Group = c("Realized", "Random")
  , Text = c(
        paste0(
            "$\\bar{r}_s = "
          , validation$Mean[1]
          , "$, $95\\%-CI = "
          , validation$LCL[1]
          , "$, $"
          , validation$UCL[1], "$"
        )
      , paste0(
            "$\\bar{r}_s = "
          , validation$Mean[2]
          , "$, $95\\%-CI = "
          , validation$LCL[2]
          , "$, $"
          , validation$UCL[2], "$"
        )
    #   "$\\bar{r}_s = -0.45$, $95%-CI = -0.48$, $-0.42$"
    # , "$\\bar{r}_s = 0.04$, $95%-CI = 0.01$, $0.08$"
  )
)

# Reorder the factors in the Group variable
dat$Group <- factor(dat$Group, levels = c("Realized", "Random"))

# Get loess
loess <- dat %>%
  group_by(Group) %>%
  nest() %>%
  mutate(data = lapply(data, function(x){
    l <- loess.sd(x$Frequency ~ x$Rank)
    df <- data.frame(
        Loess = l$y
      , Rank  = l$x
      , Upper = l$upper
      , Lower = l$lower
    )
    return(df)
  })) %>%
  unnest(data)

# Plot the data, once for realized, once for random preferences
p4 <- lapply(unique(loess$Group), function(x) {
  sub_loess <- subset(loess, Group == x)
  sub_text <- subset(text, Group == x)
  sub_data <- subset(dat, Group == x)
  ggplot(sub_loess, aes(x = Rank, y = Loess)) +
    geom_jitter(aes(x = Rank, y = Frequency)
      , data  = sub_data,
      , alpha = 0.1
      , size  = 1
    ) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper)
      , alpha = 0.4
      , fill  = "steelblue"
      , col   = "steelblue"
    ) +
    geom_line(
        size = 1
      , col  = "steelblue"
    ) +
    theme_classic() +
    coord_capped_cart(
        left   = "both"
      , bottom = "both"
      , ylim   = c(0, 50)
    ) +
    geom_text(data = sub_text
      , mapping = aes(
          x = -Inf
        , y = -Inf
        , label = TeX(Text, output = "character")
      )
      , hjust   = -0.05
      , vjust   = -0.5
      , parse   = TRUE
      , size    = 3
    ) +
    ylab("Frequency")
})

# Let's add some more information to the validation of our connectivity map
validation_con <- validation_con %>%
  subset(ModelType == "PSF") %>%
  dplyr::select(Steps, Model) %>%
  unnest(Model) %>%
  setNames(c("Steps", "Coefficient", "SE", "pvalue")) %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Add stars indicating the significance
validation_con$Significance <- sapply(1:nrow(validation_con), function(x){
  if (validation_con$pvalue[x] <= 0.01){
    return("***")
  } else if (validation_con$pvalue[x] <= 0.05){
    return("**")
  } else if (validation_con$pvalue[x] <= 0.1){
    return("*")
  }
})

# Prepare dataset for plotting confidence intervals
validation_confs <- validation_con %>%
  dplyr::select(Steps, Coefficient, LCI_90:UCI_99) %>%
  gather(key = confidence_level, value = value, LCI_90:UCI_99) %>%
  separate(col = confidence_level, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = value) %>%
  mutate(Level = paste0(Level, "%"))

# Plot of path selection model
p4[[3]] <- ggplot(data = validation_con, aes(y = Coefficient, x = as.factor(Steps))) +
  geom_point(shape = 3, size = 2.5, col = "steelblue") +
  geom_errorbar(
      aes(
        ymin = LCI
      , ymax = UCI
      , size = factor(Level)
    )
    , data   = validation_confs
    , width  = 0
    , col    = "steelblue"
    , alpha  = 0.5
  ) +
  geom_text(
      aes(label = Significance, hjust = 0.5, vjust = -0.5, angle = 90)
    , show.legend = F
  ) +
  geom_hline(
      yintercept = 0
    , color      = "darkgrey"
    , lty        = 2
    , lwd        = 0.3
  ) +
  theme_classic() +
  ylim(c(-1, 3)) +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  labs(
      x = "Number of Simulated Steps"
    # , y = expression(beta*"-Coefficient_{Connectivity}")
    , y = parse(text = TeX("$\\beta-Coefficient_{Connectivity}$"))
  ) +
  scale_size_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.3)
  ) +
  theme(
    , legend.position   = "none"
    , legend.text       = element_text(face = 3)
    , legend.title      = element_text(face = 3)
  ) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5, legend.hjust = 4))

# Put plots together
p5 <- ggarrange(p4[[1]], p4[[2]], p4[[3]], ncol = 1, labels = c("b1", "b2", "c"), label.x = -0.10)

################################################################################
#### Combine Plots
################################################################################
# Arrange the two plots nicely
p6 <- ggarrange(p3, p5, widths = c(1.5, 1), labels = "a")

# Show it
p6

# Save all
ggsave("04_Manuscript/99_MovementModel.pdf"
  , plot   = p6
  , width  = 12
  , height = 7
  , device = "pdf"
  , scale  = 0.7
)

################################################################################
#### Random Effects
################################################################################
# Rename model
mod <- best

# Check out coefficients per individual
coeffs <- coef(mod)
ranefs <- ranef(mod, condVar = T)
coeffs$cond$id
ranefs$cond$id

# Note: ranef yields the difference between the individual specific effect and
# the mean level effect. coef, on the other hand, yields the individual specific
# effect. Thus, the followin two lines yield (approximately) the same
mean(coeffs$cond$id$cos_ta_) + ranefs$cond$id$cos_ta_[3]

# We now want to visualize the individual variation. There are two possibilities
# for this: lme4::dotplot() or a ggplot. The dotplot is easier, yet not
# customizable. Let's first do the dotplot, then recreate it in ggplot.
lme4:::dotplot.ranef.mer(ranef(mod)$cond)

# Maybe scalefree?
lme4:::dotplot.ranef.mer(ranef(mod)$cond, scales = list(x = list(relation = "free")))

# Prepare dataframe that we need to plot the same in ggplot
rfs <- ranefs$cond$id %>%
  rownames_to_column() %>%
  gather(key = Covariate, value = Mean, 2:8)
names(rfs)[1] <- "id"

# We need to add the conditional variance
condVar <- attributes(ranefs$cond$id)$condVar
names(condVar) <- attributes(ranefs$cond$id)$names
condVar <- as.data.frame(do.call(rbind, condVar))
names(condVar) <- attributes(ranefs$cond$id)$row.names
condVar <- rownames_to_column(condVar)
names(condVar)[1] <- "Covariate"
condVar <- gather(condVar, key = id, value = Variance, 2:17)

# Join data to rfs dataframe
rfs <- left_join(rfs, condVar)

# Rename stuff nicely
rfs$Covariate <- gsub(rfs$Covariate, pattern = "cos_ta_", replacement = "cos(ta)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "log_sl_", replacement = "log(sl)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "sl_", replacement = "sl")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "HumansBuff5000", replacement = "HumanInfluence")

# Make covariates a factor
rfs$Covariate <- factor(rfs$Covariate, levels = c(
  "cos(ta)", "sl", "log(sl)", "Water", "DistanceToWater", "Shrubs"
  , "Woodland", "HumanInfluence"
))

# Visualize. Note that I am transforming the variance using mean - 2 *
# sqrt(Variance). This was taken from here: https://stackoverflow.com/questions
# /13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-
# how-to-mak
p <- ggplot(rfs, aes(x = Mean, y = id)) +
  geom_point() +
  facet_wrap("Covariate", nrow = 2) +
  geom_errorbarh(aes(
      xmin = Mean - 2 * sqrt(Variance)
    , xmax = Mean + 2 * sqrt(Variance)
  ), colour = "black", height = 0)
