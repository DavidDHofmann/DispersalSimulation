################################################################################
#### Some Descriptive Metrics
################################################################################
# Description: A simple plot of the Dispersal Durations

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # For spatial stuff
library(rgdal)        # For spatial stuff
library(rgeos)        # For spatial stuff
library(davidoff)     # Custom functions
library(pbmcapply)    # For multicoring

################################################################################
#### Dispersal Simulation
################################################################################
# Get an idea of the dispersal simulation duration
info <- read_csv("03_Data/03_Results/99_Simulations/Report.csv")

# Average and sd of simulation duration
hist(info$duration, breaks = 20)
mean(info$duration)
sd(info$duration)

# Load simulated and observed trajectories
sims  <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
obs   <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF).csv")

# Check how many simulated trajectories eventually hit a map boundary
sims %>%
  select(TrackID, BoundaryHit, Area) %>%
  group_by(TrackID) %>%
  summarize(BoundaryHit = max(BoundaryHit, na.rm = T), Area = unique(Area)) %>%
  group_by(Area) %>%
  summarize(BoundaryHit = mean(BoundaryHit))

# Get an idea of the step lengths and turning angles of simulated and observed dispersers
obs   <- subset(obs, case_)
sims <- select(sims, ta_, sl_)
obs <- select(obs, ta_, sl_)

# Put them together
sims$Type <- "Simulated"
obs$Type <- "Observed"
dat <- rbind(sims, obs)

# Average and sd of simulated step lengths
dat %>%
  group_by(Type) %>%
  summarize(
      mean_cos_ta = mean(cos(ta_), na.rm = T)
    , sd_cos_ta   = sd(cos(ta_), na.rm   = T)
    , mean_sl     = mean(sl_, na.rm      = T)
    , sd_sl       = sd(sl_, na.rm        = T)
  )

################################################################################
#### Areas Reached (Western Delta)
################################################################################
# Load simulated trajectories
sims  <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Define a polygon for the area representing the "south-western" okavango delta
x_coord <- c(23.12376278, 22.75065198, 22.56930762, 22.30250212, 22.21287215 , 22.32543072, 22.40672302, 22.54950565, 22.68603502, 22.90177227, 23.12376278)
y_coord <- c(-20.22675467, -20.49772899, -20.56859920, -20.27052744, -19.97870893 , -19.48678631, -19.29293544, -19.38048099, -19.53264350, -19.80778667, -20.22675467)
pol <- cbind(x_coord, y_coord)
pol <- Polygon(pol)
pol <- Polygons(list(pol), 1)
pol1 <- SpatialPolygons(list(pol))

# Define a polygon for the area representing the "south-western" okavango delta
x_coord <- c(23.11542511, 23.22589926, 23.55315287, 24.21599777, 24.43069281, 24.30979657, 23.92834810, 23.39056827, 23.11542511)
y_coord <- c(-18.51127872, -18.21946022, -17.95265472, -17.76714153, -18.02144052, -18.39246690, -18.70512959, -18.75932446, -18.51127872)
pol <- cbind(x_coord, y_coord)
pol <- Polygon(pol)
pol <- Polygons(list(pol), 2)
pol2 <- SpatialPolygons(list(pol))

# Bind them
pol <- rbind(pol1, pol2)

# Determine all simulations leaving from Moremi
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory leaves
first$From <- as.character(over(first, prot)$Name)
first <- first@data[, c("TrackID", "From")]

# Join information to simulated tracks
sims <- left_join(sims, first, by = "TrackID")

# Keep only simulations leaving from moremi
sims <- subset(sims, From == "Moremi")

# Coerce those to proper trajectories
tracks <- sims2tracks(sims)

# Visualize them
plot(tracks, lwd = 0.2)
plot(subset(prot, Name == "Moremi"), add = T, col = "red")

# Let's check the number of tracks that reach the western section of the delta
success <- gIntersects(pol, tracks, byid = T)
colSums(success)

################################################################################
#### Areas Reached (Duration)
################################################################################
# Reload data on areas reached
visits <- read_rds("03_Data/03_Results/99_AreasReached.rds")

# Load protected areas
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Keep only national parks
prot <- subset(prot, Desig == "National Park")

# Assign a unique ID to each of the areas
prot$ID <- 1:nrow(prot)

# Extract names
names <- prot@data[, c("Name", "ID")]

# Replace from to values with proper names
visits <- left_join(visits, names, by = c("From" = "ID"))
visits <- left_join(visits, names, by = c("To" = "ID"))
visits$From <- visits$Name.x
visits$To   <- visits$Name.y
visits$Name.x <- NULL
visits$Name.y <- NULL

# Find average duration from Moremi to Chobe and Hwange
subset(visits, From == "Moremi" & To == "Chobe")
subset(visits, From == "Moremi" & To == "Hwange")
