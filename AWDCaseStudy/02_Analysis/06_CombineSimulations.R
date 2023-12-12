################################################################################
#### Combine Simulations into Single Object
################################################################################
# Description: Take all simulations and put them into a single r file

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # For handling spatial data
library(rgdal)        # For handling spatial data
library(pbmcapply)    # For multicore processes with progress bar
library(davidoff)     # Access to custom functions
library(rgeos)        # For spatial data manipulation
library(lemon)        # For capped coordinates
library(adehabitatLT) # For some simple checks
library(Matrix)       # To store matrices more efficiently

################################################################################
#### Putting Data Together
################################################################################
# Identify all simulation files
files_main <- dir(
    "03_Data/03_Results/99_Simulations/Main"
  , pattern     = ".rds$"
  , full.names  = T
)
files_buffer <- dir(
    "03_Data/03_Results/99_Simulations/Buffer"
  , pattern     = ".rds$"
  , full.names  = T
)

# Put them together
files <- data.frame(
    file = c(files_main, files_buffer)
  , area = c(rep("Main", length(files_main)), rep("Buffer", length(files_buffer)))
  , stringsAsFactors = F
)

# How many are there?
nrow(files)

# Load the files and bind their rows together
sims <- lapply(1:nrow(files), function(x){
  dat <- read_rds(files$file[x])
  dat$Area <- files$area[x]
  dat$SimID <- x
  return(dat)
}) %>% do.call(rbind, .)

# Take a look at the data
head(sims)

# Count the number of simulated steps (in Mio)
nrow(sims) / 1e6

# Because in each simulation we start off with new IDs for trajectories, they
# are not unique across simulations. We thus combine ID and SimID to create an
# ID that is unqiue to each simulated path, across all simulations
sims <- sims %>% mutate(TrackID = group_indices(., SimID, TrackID))

# Make sure it worked
table(table(sims$TrackID))

# Collect garbage
gc()

# Let's also create a step counter, indicating the number of the step in its
# respective trajectory
sims <- sims %>% group_by(TrackID) %>% mutate(StepNumber = (row_number()))

# Check object size
format(object.size(sims), units = "Gb")

# Ungroup
sims <- ungroup(sims)

# Write to an rds
write_rds(sims, "03_Data/03_Results/99_DispersalSimulation.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Remove raw data of simulations to free some space
file.remove(files_main)
file.remove(files_buffer)

################################################################################
#### Interpolated Simulation
################################################################################
# To compute the betweenness, we will want to use simulations where the steps
# were interpolated. We will store those interpolated simulations to a special
# subdirectory.
dir.create("03_Data/03_Results/99_Simulations/Interpolated", showWarnings = F)

# Now let's loop through all simulated paths and create interpolated versions
sims <- nest(sims, Data = -TrackID)
sims$Interpolated <- pbmclapply(
    X                  = 1:nrow(sims)
  , ignore.interactive = T
  , mc.cores           = detectCores() / 2
  , FUN                = function(x) {
    filename <- paste0(
        "03_Data/03_Results/99_Simulations/Interpolated/Track_"
      , sprintf("%05d", sims$TrackID[x])
      , ".rds"
    )
    if (!file.exists(filename)) {
      path  <- sims$Data[[x]]
      inter <- lapply(1:(nrow(path) - 1), function(i) {
        xy_new <- interpolatePointsC(
            x1 = path$x[i]
          , x2 = path$x[i + 1]
          , y1 = path$y[i]
          , y2 = path$y[i + 1]
          , by = 1000 / 111000
        )
        int <- path[rep(i, nrow(xy_new)), ]
        int[, c("x", "y")] <- xy_new
        return(int)
      }) %>% do.call(rbind, .)
      inter$Segment <- 1:nrow(inter)
      write_rds(inter, filename)
    }
  return(filename)
})
sims$Interpolated <- NULL
sims <- unnest(sims, Data)

# ################################################################################
# #### Verify Turning Angles and Step Lengths
# ################################################################################
# # Subset to a few trajectories
# sub <- subsims(sims, nid = 200, area = "Main")
#
# # Coerce tracks to ltraj
# coordinates(sub) <- c("x", "y")
# lt <- as.ltraj(coordinates(sub), date = sub$Timestamp, id = sub$TrackID)
#
# # Extract absolute turning angles, relative turning angles and step lengths
# abstas <- sapply(lt, function(x){x$abs.angle})
# reltas <- sapply(lt, function(x){x$rel.angle})
# sls <- sapply(lt, function(x){x$dist})
#
# # Draw histograms
# hist(abstas, breaks = 100)
# hist(reltas, breaks = 100)
# hist(sls, breaks = 100)

################################################################################
#### Seek Convergence
################################################################################
# For this, we'll only look at the trajectories from the main area
sims <- sims[sims$Area == "Main", ]
gc()

# Coerce simulations to spatial lines
tracks <- sims2tracks(sims, keep.data = F, steps = 2000, cores = detectCores() / 2)
gc()

# Check data
show(tracks)

# Load reference shapefile
s <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")

# Distribute 1000 "chek points" in our study area. For this, generate 5km2
# squares in our study area
check <- raster(s, res = metersToDegrees(5000), vals = 0)
set.seed(12345)
values(check)[sample(1:ncell(check), 1000)] <- 1
check <- rasterToPolygons(check, fun = function(x){x == 1})

# Assign a checkpoint ID to each polygon
check$CheckID <- 1:length(check)

# Extract all trajectory ids
ids <- unique(tracks$TrackID)
length(ids)

# Write the checkpoints to file
check$layer <- NULL
writeOGR(check, "03_Data/03_Results", "99_Checkpoints", driver = "ESRI Shapefile")

# Load source areas
source_areas <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
buffer_areas <- readOGR("03_Data/03_Results/99_BufferArea.shp")

# Visualize them
plot(buffer_areas
  , axes   = T
  , las    = 1
  , col    = colTrans("red")
  , border = NA
  , main   = "Checkpoints in Study Area"
)
plot(source_areas
  , col    = colTrans("darkgreen")
  , border = NA
  , main   = "Checkpoints in Study Area"
  , add    = T
)
plot(check
  , cex    = 0.3
  , pch    = 20
  , col    = "blue"
  , border = NA
  , add    = T
)
plot(tracks[sample(nrow(tracks), 10), ]
  , add = T
  , col = "gold"
  , lwd = 0.6
)
legend("bottomright"
  , legend = c("Checkpoints", "Source Areas", "Buffer Zone", "Simulated Dispersers")
  , pch    = c(15, 15, 15, NA)
  , lty    = c(NA, NA, NA, 1)
  , col    = c("blue", colTrans("darkgreen"), colTrans("red"), "gold")
)

# Free space
rm(sims, buffer_areas, source_areas, files_main, files_buffer, ints, s, abstas, reltas, sls, sub)
gc()

# Check all intersections with checkpoints
ints <- gIntersects(check, tracks, byid = T, prepared = T) * 1

# Convert to sparse matrix
ints <- Matrix(ints)

# Store to file
write_rds(ints, "03_Data/03_Results/99_Intersections.rds")

# Free further space
rm(tracks)
gc()

# Prepare a design matrix through which we will loop
design <- expand_grid(
    NTracks    = seq(0, nrow(ints), by = 500)
  , Replicate  = c(1:100)
)

# Let's shuffle the matrix (makes prediction of calculation time in pbmclapply
# more accurate)
design <- design[sample(nrow(design), replace = F), ]

# Calculate traversal frequencies
design <- mutate(design, Frequency = pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() / 4
  , FUN                = function(x){

  # Sample rows
  index <- sample(nrow(ints), size = design$NTracks[x], replace = T)

  # Subset to corresponding rows (allows for duplicates)
  ints_sub <- ints[index, ]

  # Calculate rowsums
  ints_sub <- colSums(ints_sub)
  ints_sub <- enframe(ints_sub, name = "CheckID", value = "Traversals")

  # Return it
  return(ints_sub)
  gc()

}))

# Clean results and calculate relative traversal frequency
convergence <- design %>%
  unnest(cols = Frequency) %>%
  mutate(RelativeTraversals = Traversals / NTracks) %>%
  mutate(RelativeTraversals = ifelse(is.nan(RelativeTraversals)
  , 0
  , RelativeTraversals)
)

# Store results
write_rds(convergence, "03_Data/03_Results/99_Convergence.rds")
