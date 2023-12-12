################################################################################
#### Plot of Study Area
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(rgdal)        # To load and store spatial data
library(raster)       # To manipulate raster data
library(rgeos)        # To manipulate spatial data
library(terra)        # To handle raster data
library(tmap)         # For beautiful maps
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(grid)         # To arrange multiple plots
library(Cairo)        # To store plots
library(elevatr)      # To create nice hillshade map

################################################################################
#### Plot of the Study Area
################################################################################
# Load required data
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
dogs    <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
water   <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")

# Clean up the water file for nicer illustration
water <- aggregate(water, fun = max, fact = 2)
water <- raster(water)
clump <- clump(water)
clump_freq <- data.frame(freq(clump))
clump_freq <- subset(clump_freq, count <= 40)
rcl <- data.frame(old = clump_freq$value, new = NA)
water <- reclassify(clump, rcl = rcl)
water <- !is.na(water)

# Remove 0 from water layer
water[water == 0] <- NA
plot(water)

# Simplify Protection zones
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"
prot$Desig <- as.factor(as.character(prot$Desig))

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza) + c(-1, +1, -1, +1) * 50 / 111, "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 44, 35)
  , y = c(-8, -30, -27, -14, -24)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(37, -14), c(29, -15)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Visualize
plot(labels_countries)
plot(lines_countries, add = T)
plot(lines_countries[4], add = T, col = "red")
plot(africa, add = T)

# Prepare a plot of the KAZA. Create labels for countries first.
labels_countries2 <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries2) <- c("x", "y")
crs(labels_countries2) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.6, 22.9, 27.8, 27.25)
  , y     = c(-19, -18.2, -17, -20.7, -15.0, -14.4, -15.58)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
    , "Barotse\nFloodplain", "Lukanga\nSwamp", "Kafue\nFlats"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Visualize
plot(water)
plot(labels_waters, add = T, col = "red")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78, 22.63, 27.92, 28.54)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
    , "Chizarira", "Matobo"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Visualize
plot(subset(prot, Desig == "National Park"))
plot(labels_nationalparks, add = T, col = "red")

# Create elevation map for africa
elev_africa <- get_elev_raster(africa, z = 4)
elev_africa <- crop(elev_africa, africa)
elev_africa <- mask(elev_africa, africa, updatevalue = NA)
terrain_africa <- terrain(elev_africa ** 1.6, opt = c("slope", "aspect"))
hill_africa <- hillShade(aspect = terrain_africa$aspect, slope = terrain_africa$slope)

# Visualize it
plot(hill_africa, col = gray(80:100/100))
plot(africa, add = T, border = "gray80")

# Create elevation map for kaza
ext <- as(extent(kaza) + c(-1, 1, -1, 1), "SpatialPolygons")
crs(ext) <- CRS("+init=epsg:4326")
elev_kaza <- get_elev_raster(kaza, z = 6)
elev_kaza <- crop(elev_kaza, ext)
elev_kaza <- mask(elev_kaza, ext, updatevalue = NA)
terrain_kaza <- terrain(elev_kaza ** 1.6, opt = c("slope", "aspect"))
hill_kaza <- hillShade(aspect = terrain_kaza$aspect, slope = terrain_kaza$slope)

# Visualize it
plot(hill_kaza, col = gray(80:100/100))
plot(africa, add = T, border = "gray80")
plot(kaza, add = T, border = "gray80")

# Prepare a map of Africa
p1 <- tm_shape(hill_africa) +
  tm_raster(
      style       = "cont"
    , palette     = gray(40:100/100)
    , alpha       = 0.2
    , legend.show = F
  ) +
    tm_shape(dogs) +
  tm_polygons(
      col        = "orange"
    , border.col = "orange"
  ) +
  tm_shape(africa, is.master = T) +
    tm_borders(
        col        = "gray40"
      , lwd        = 0.2
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 1.5
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "blue"
      , lty = 1
      , lwd = 0.5
    ) +
  # tm_shape(lines_countries) +
  #   tm_lines(
  #       col = "gray30"
  #     , lty = 1
  #     , lwd = 1
  #   ) +
  # tm_shape(labels_countries) +
  #   tm_text(
  #       text = "Label"
  #     , size = 0.6
  #     , col  = "black"
  #   ) +
  tm_layout(
      asp         = 0.8
    , frame       = F
    , legend.show = F
    , bg.color    = "transparent"
  ) +
  tm_credits("a"
    , position  = c("right", "top")
    , size      = 1.5
    , col       = "black"
    , fontface  = "bold"
)

# Plot of the KAZA
p2 <- tm_shape(hill_kaza) +
  tm_raster(
      style       = "cont"
    , palette     = gray(60:100/100)
    , alpha       = 0.2
    , legend.show = F
  ) +
  tm_shape(prot) +
    tm_polygons(
        col         = "Desig"
      , palette     = c("#70ab70", "#d9f0d3")
      , lwd         = 0.1
      , border.col  = "#6ba36b"
      , legend.show = F
    ) +
  tm_shape(water) +
    tm_raster(
        palette     = "#96d0ff"
      , legend.show = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray40"
      , lwd = 0.5
    ) +
  tm_shape(labels_countries2) +
    tm_text("Label"
      , col       = "black"
      , fontface  = 2
      , size      = 1.5
    ) +
  tm_shape(labels_waters) +
    tm_text("Label"
      , fontface  = 3
      , size      = 0.5
      , col       = "#064886"
    ) +
  tm_shape(labels_nationalparks) +
    tm_text("Label"
      , col       = "#004800"
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_graticules(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
    # , frame                   = "gray20"
    , frame                   = "blue"
    , frame.lwd               = 1
    , asp                     = 1.2
    , legend.outside          = TRUE
    , legend.outside.position = "left"
    , legend.stack            = "vertical"
    , legend.text.size        = 0.8
  ) +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "black"
      , width     = 0.125
  ) +
  tm_credits("b"
    , position = c("left", "top")
    , size     = 1.5
    , col      = "black"
    , fontface = "bold"
  ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
    , position    = c("left", "bottom")
)

################################################################################
#### Prepare a common Legend
################################################################################
# We now prepare a legend that the two plots share. However, we have to assign
# to one of the two plots. Let's use p2 for this
p3 <- p2 + tm_add_legend(
    type    = "fill"
  , labels  = c(
      "Wild Dog Populations"
    , "Water Bodies"
    , "Protected Areas"
    , "National Parks (NP)"
  )
  , col = c(
      "orange"
    , "#96d0ff"
    , "#d9f0d3"
    , "#70ab70"
  )) + tm_add_legend(
    type    = "line"
  , labels  = c("KAZA-TFCA Borders", "Country Borders")
  , col     = c("black", "gray40")
  , lty     = c(1, 1)
  , lwd     = c(2, 1)
) + tm_layout(legend.frame.lwd = 2, legend.text.size = 1.05)

# Show it
p3

# Store the plot
tmap_save(
    tm        = p3
  , filename  = "04_Manuscript/99_StudyArea.png"
  , insets_tm = p1
  , insets_vp = viewport(0.164, 0.32, width = 0.66, height = 0.66)
  , width     = 9
  , height    = 5.25
  , scale     = 1.1
)
