################################################################################
#### Map Comparisons
################################################################################
# Description: In this script, we compare connectivity metrics inside and
# outside the KAZA-TFCA and produce a nice plot to visualize the results

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # For general raster manipulation
library(rgdal)        # To read spatial data
library(tidyverse)    # For data wrangling
library(viridis)      # For nicer colors
library(ggpubr)       # To arrange multiple plots

# Reload connectivity maps
heatmaps <- read_rds("03_Data/03_Results/99_Heatmaps.rds")
betweenness <- stack("03_Data/03_Results/99_Betweenness.grd")

# Apply focal filter to buffer/smooth maps
betweenness <- lapply(1:nlayers(betweenness), function(x){
  focal(betweenness[[x]], w = matrix(1, 3, 3), fun = mean)
}) %>% stack()

# Also, we want to square root
betweenness <- sqrt(betweenness)

# Subset to steps of interest
heatmaps <- subset(heatmaps, steps %in% c(125, 500, 2000))

# We also need to add the main area and the buffer
heatmaps <- lapply(unique(heatmaps$steps), function(x) {
  sub <- subset(heatmaps, steps == x)
  ras <- stack(sub$heatmap)
  ras <- sum(ras)
  result <- tibble(Steps = x, Heatmap = list(ras))
  return(result)
}) %>% do.call(rbind, .)

# Convert the betweenness maps into a tibble
betweenness <- tibble(
    Steps       = c(125, 500, 2000)
  , Betweenness = list(betweenness[[1]], betweenness[[2]], betweenness[[3]])
)

# Put the maps together
maps <- bind_cols(heatmaps, betweenness[, 2]) %>%
  pivot_longer(cols = 2:3, names_to = "Type", values_to = "Map")

################################################################################
#### Median Values Inside and Outside KAZA
################################################################################
# Load required shapefiles
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
stud <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")

# Loop through the different heatmaps and compare heat inside and outside kaza
maps <- mutate(maps, Comparison = map(Map, function(x) {
  sub <- crop(x, stud, snap = "in")
  inside <- na.omit(values(mask(sub, kaza)))
  outside <- na.omit(values(mask(sub, kaza, inverse = T)))
  dat <- data.frame(
      Where = c(rep("Inside KAZA-TFCA", length(inside)), rep("Outside KAZA-TFCA", length(outside)))
    , Value = c(inside, outside)
  )
  return(dat)
}))

# Generate a tidy dataframe for plotting
dat <- arrange(maps, Type, Steps) %>%
  dplyr::select(Steps, Type, Comparison) %>%
  unnest(Comparison) %>%
  mutate(Steps = as.factor(Steps), Type = as.factor(Type))

# Summarize values
summary <- dat %>%
  mutate(Value = ifelse(Type == "Betweenness", Value ** 2, Value)) %>%
  group_by(Steps, Type, Where) %>%
  summarize(
      Median  = median(Value)
    , IQR     = IQR(Value)
    , .groups = "drop"
  )
print(summary)

# Plot the comparison of the heatmap
p1 <- ggplot(subset(dat, Type == "Heatmap"), aes(x = Value, y = as.factor(Steps), color = Where, fill = Where)) +
  geom_boxplot(lwd = 0.5, outlier.size = 0.1, alpha = 0.2) +
  # geom_violin(scale = "area", lwd = 0.5, alpha = 0.2, adjust = 100, trim = F) +
  scale_color_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  scale_fill_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  theme_classic() +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.2)
    , legend.position  = c(0.9, 0.25)
  ) +
  xlab("# Traversing Trajectories") +
  ylab("# Simulated Steps") +
  ggtitle("Heatmap")

# Plot the comparison of the betweenness map
p2 <- ggplot(subset(dat, Type == "Betweenness"), aes(x = Value, y = as.factor(Steps), color = Where, fill = Where)) +
  geom_boxplot(lwd = 0.5, outlier.size = 0.1, alpha = 0.2) +
  # geom_violin(scale = "area", lwd = 0.5, alpha = 0.2, adjust = 100, trim = F) +
  scale_color_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  scale_fill_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  theme_classic() +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.2)
    , legend.position  = c(0.9, 0.26)
  ) +
  xlab(expression("Betweenness" ^ "0.5")) +
  ylab("# Simulated Steps") +
  ggtitle("Betweenness Map")

# Arrange the plots
p <- ggarrange(p1, p2, nrow = 2, labels = c("a", "b"))
# ggsave("04_Manuscript/99_MapComparison.png", plot = p, width = 8, height = 7)
ggsave("04_Manuscript/99_MapComparison.png", plot = p, width = 7, height = 5)

################################################################################
#### Median Values Inside and Outside KAZA by Country
################################################################################
# Load required shapefiles
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
stud <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")
afri <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")

# Keep only the countries of interest
afri <- subset(afri, COUNTRY %in% c("Angola", "Namibia", "Botswana", "Zimbabwe", "Zambia"))

# Aggregate by country
afri <- aggregate(afri, by = "COUNTRY")
afri <- crop(afri, stud)

# Loop through the different heatmaps and compare heat inside and outside kaza
# and inside or outside the different countries
maps <- mutate(maps, ComparisonCountry = map(Map, function(x) {
  sub <- crop(x, stud, snap = "in")
  values <- lapply(unique(afri$COUNTRY), function(z) {
    sub_country <- mask(x, subset(afri, COUNTRY == z))
    inside <- na.omit(values(mask(sub_country, kaza)))
    outside <- na.omit(values(mask(sub_country, kaza, inverse = T)))
    dat <- data.frame(
        Where   = c(rep("Inside KAZA-TFCA", length(inside)), rep("Outside KAZA-TFCA", length(outside)))
      , Value   = c(inside, outside)
      , Country = rep(z, (length(inside) + length(outside)))
    )
    return(dat)
  })
  values <- do.call(rbind, values)
}))

# Generate a tidy dataframe for plotting
dat <- arrange(maps, Type, Steps) %>%
  dplyr::select(Steps, Type, ComparisonCountry) %>%
  unnest(ComparisonCountry) %>%
  mutate(Steps = as.factor(Steps), Type = as.factor(Type))

# Plot the comparison of the heatmap
p1 <- ggplot(subset(dat, Type == "Heatmap"), aes(x = Value, y = as.factor(Steps), color = Where, fill = Where)) +
  geom_boxplot(lwd = 0.5, outlier.size = 0.1, alpha = 0.2) +
  # geom_violin(scale = "area", lwd = 0.5, alpha = 0.2, adjust = 100, trim = F) +
  scale_color_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  scale_fill_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  theme_classic() +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.2)
    , legend.position  = c(0.8, 0.05)
  ) +
  xlab("# Traversing Trajectories") +
  ylab("# Simulated Steps") +
  ggtitle("Heatmap") +
  facet_wrap(~ Country, ncol = 1)

# Plot the comparison of the betweenness map
p2 <- ggplot(subset(dat, Type == "Betweenness"), aes(x = Value, y = as.factor(Steps), color = Where, fill = Where)) +
  geom_boxplot(lwd = 0.5, outlier.size = 0.1, alpha = 0.2) +
  # geom_violin(scale = "area", lwd = 0.5, alpha = 0.2, adjust = 100, trim = F) +
  scale_color_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  scale_fill_manual(values = c("#5B9BD5", "orange"), name = "Location") +
  theme_classic() +
  theme(
      panel.grid.major = element_line(colour = "gray90", size = 0.2)
    , legend.position  = c(0.8, 0.05)
  ) +
  xlab(expression("Betweenness" ^ "0.5")) +
  ylab("# Simulated Steps") +
  ggtitle("Betweenness Map") +
  facet_wrap(~ Country, ncol = 1)

# Arrange the plots
p <- ggarrange(p1, p2, ncol = 2, labels = c("a", "b"))
ggsave("04_Manuscript/99_MapComparisonCountries.png", plot = p, width = 7, height = 9, scale = 1.2)

################################################################################
#### Mean Traversal Frequency South of Linyanti
################################################################################
# Create polygon around dispersal hotspot south of Linyanti and extract heat
ext <- as(extent(c(23.3241, 24.4825, -18.7901, -18.2496)), "SpatialPolygons")
ext <- trim(mask(maps$Map[[5]], ext))
plot(ext)
median(values(ext))
IQR(values(ext))

# Create polygon around dispersal hotspot south of Linyanti and extract
# betweenness
ext <- as(extent(c(23.3241, 24.4825, -18.7901, -18.2496)), "SpatialPolygons")
ext <- trim(mask(maps$Map[[6]], ext))
plot(ext)
median(values(ext))
IQR(values(ext))
