################################################################################
#### Preparation of the Protected Areas Layer
################################################################################
# Description: Preparation of the shapefile of protected areas as downloaded
# from the Peace Parks Foundation website
# (http://new-ppfmaps.opendata.arcgis.com/datasets/ppf-protected-areas-
# detailed?geometry=-13.87%2C-25.558%2C69.846%2C-11.001)

# Clear R's brain
rm(list = ls())

# Load packages
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(gdalUtils)    # To manipulate spatial data
library(cleangeo)     # To clean spatial data
library(tidyverse)    # For data wrangling
library(RColorBrewer) # To access color palettes
library(terra)        # For quicker rasterization

################################################################################
#### Simplify Categories
################################################################################
# Import file
prot <- readOGR("03_Data/01_RawData/PEACEPARKS/PPF_Protected_Areas_Detailed.shp")

# Use the reference shapefile to crop the areas
r <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")
prot <- crop(prot, r)

# Keep only the attributes of interest and rename them nicely
prot@data <- dplyr::select(prot@data
  , Name    = Name
  , Desig   = Designatio
  , IUCN    = IUCN
  , Country = Country
)

# We now want to simplify the protection categories. We created a
# reclassification table for this, so let's use it
desigs <- "03_Data/01_RawData/PEACEPARKS/Reclassification.csv" %>%
  read_csv() %>%
  set_names(., c("Nr", "Old", "New", "Comment"))

# Join the dataframes
prot@data <- left_join(prot@data, desigs, by = c("Desig" = "Old"))

# Check out the distribution of the new categories
table(prot$New)

# Game reserves in Botswana serve the same purpose as national parks. Let's thus
# reclassify them accordingly
prot$New[prot$Desig == "Game Reserve" & prot$Country == "Botswana"] <-
  "National Park"

# Remove columns that we dont need anymore
prot@data <- prot@data %>% dplyr::select(-c("Desig", "Nr", "Comment"))

# Rename the remaining columns
prot@data <- prot@data %>% rename(Desig = New)

# Delete the objects and attributes that are not needed
prot <- subset(prot, prot$Desig != "Delete")

# Check validity of shapefile
gIsValid(prot, reason = TRUE)

# Plot the stuff with each designation in a different colour. Note that the
# country borders will look pretty nasty due to the low resolution.
map <- shapefile("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
u <- unique(prot$Desig)
m <- match(prot$Desig, u)
n <- length(unique(prot$Desig))
pal <- brewer.pal(n, "Greens")
plot(prot, col = pal[m])
plot(map, add = TRUE)
text(subset(prot, prot$Desig == "National Park"), 'Name', cex = 0.5, halo = TRUE)
legend("topleft", legend = u, col = pal, pch = 19)

# Let's assign a value to each class
values <- data.frame(
    Desig  = c("Forest Reserve", "Protected", "National Park")
  , Values = 1:3
)

# Join the values to the shapefile
prot@data <- left_join(prot@data, values, by = "Desig")

# Store the layer
writeOGR(prot
  , "03_Data/02_CleanData"
  , "02_LandUse_Protected_PEACEPARKS"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Rasterize Protected Areas
################################################################################
# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Coerce protected areas to vect
prot <- vect(prot)

# Rasterize protected areas
prot_r <- rasterize(
    x           = prot
  , y           = r
  , field       = "Values"
  , background  = 0
)

# Create also a binary map
prot_rb <- prot_r > 0

# Create dataframe with class information
info <- data.frame(
    Class = c("Unprotected", "ForestReserve", "Protected", "NationalPark")
  , Code  = 0:3
)

# Visualize
plot(raster(prot_r))
plot(raster(prot_rb))

# Convert to raster and store
writeRaster(
    x         = raster(prot_rb)
  , filename  = "03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS(1Class).tif"
  , overwrite = T
)

# Store the categorical raster
writeRaster(
    x         = raster(prot_r)
  , filename  = "03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS(3Classes).tif"
  , overwrite = T
)

# Store the information table
write_csv(info, "03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS(3Classes).csv")
