################################################################################
#### Function to Interpolate Between Two Coordinates
################################################################################
# Extracting covariate values along spatial lines can be really really slow. To
# improve efficiency, you can therefore extract covariates along points that are
# distributed on the line instead. Note that you may want to modify this
# function if you work with non-projected data!

# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, by = 1){

  # Calculate length of line between points
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

  # Calculate how many segments we need
  nsegs <- max(ceiling(length / by), 1)

  # Interpolate between points
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)
  return(cbind(x, y))
}

# # Try it
# interpolatePoints(x1 = 0, x2 = 10, y1 = 0, y2 = 10, by = 2)
#
# # Let's generate a random line
# start <- runif(n = 2, min = 0, max = 100)
# end <- runif(n = 2, min = 0, max = 100)
# line <- spLines(SpatialPoints(rbind(start, end)))
#
# # Now distribute points on that line using our custom function
# # The smaller the value for "by", the more points are put on the line
# line_inter <- interpolatePoints(
#     x1 = start[1]
#   , x2 = end[1]
#   , y1 = start[2]
#   , y2 = end[2]
#   , by = 0.1
# )
#
# # Plot the line and interpolated coordinates
# plot(line)
# points(line_inter, col = "red", pch = 20, cex = 0.5)
#
# # Are results the same?
# a <- raster::extract(cov, line, fun = "mean")
# b <- colMeans(raster::extract(cov, line_inter))
# rbind(along_line = as.vector(a), at_points = b)
#
# # Let's see how extraction speeds compare
# benchmark <- microbenchmark(
#     AlongLine = raster::extract(cov, line, fun = "mean")
#   , AtPoints  = colMeans(raster::extract(cov, line_inter, fun = "mean"))
#   , times     = 20
# )
# summary(benchmark)


################################################################################
#### Function to Interpolate and entire Path
################################################################################
# Function that takes a path (sequence of steps) and interpolates to a certain
# maximum distance
interpolatePath <- function(x, y, eps = 0.1) {
  inter <- lapply(1:(length(x) - 1), function(i) {
    xy_new <- interpolatePoints(
        x1 = x[i]
      , x2 = x[i + 1]
      , y1 = y[i]
      , y2 = y[i + 1]
      , by = eps
    )
    xy_new <- as.data.frame(xy_new)
    xy_new$segment_id <- i
    return(xy_new)
  })
  inter <- do.call(rbind, inter)
  return(inter)
}

# Try it
# path <- data.frame(x = c(1, 5, 3), y = c(0, 2, 10))
# path_inter <- interpolatePath(x = path$x, y = path$y, eps = 0.25)
# plot(y ~ x, data = path, type = "o", pch = 20, xlim = c(0, 10), ylim = c(0, 10))
# points(y ~ x, data = path_inter, pch = 20, col = "red")

################################################################################
#### Von Mises Distribution Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu){
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01){
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

# # Let's ensure that the function works as expected
# x <- seq(-pi, +pi, by = 0.01)
# y1 <- dvonmises(x, mu = 0, kappa = 0)
# y2 <- dvonmises(x, mu = 0, kappa = 1)
# y3 <- dvonmises(x, mu = 0, kappa = 2)
#
# # Visualize the density for different kappas
# plot(NA, xlim = c(-pi, +pi), ylim = c(0, 0.6), xlab = "Turning Angle", ylab = "Prob. Density")
# abline(v = 0, lty = 2, col = "gray80")
# lines(y1 ~ x, col = "blue")
# lines(y2 ~ x, col = "purple")
# lines(y3 ~ x, col = "red")
# legend("topright"
#   , lty    = 1
#   , col    = c("blue", "purple", "red")
#   , legend = c("kappa = 0", "kappa = 1", "kappa = 2")
# )

################################################################################
#### Function to Simulate Movement from a fitted SSF
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL    # Source point (in matrix form -> n * 2)
    , covars   = NULL    # Stack of covariate layers (names need to match formula!)
    , formula  = NULL    # Model formula used to predict selection score
    , prefs    = NULL    # Preferences used to predict selection score
    , sl_dist  = NULL    # Parameters describing the step length distribution
    , ta_dist  = NULL    # Parameters describing the turning angle distribution
    , ext      = NULL    # Extent on which animals are allowed to move
    , n_steps  = 10      # Number of steps simulated
    , n_rsteps = 25      # Number of random steps proposed at each step
    , stop     = TRUE    # Should the simulation stop at boundaries?
  ){

  # For testing only
  # xy        <- coordinates(spsample(nps, type = "random", n = 1))
  # covars    <- cov
  # formula   <- ~ elev + dist
  # prefs     <- c(0.5, -0.2) # They need to match the formula!!!
  # sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 0.5))
  # ta_dist   <- list(name = "vonmises", params = list(kappa = 0.2, mu = 0))
  # n_rsteps  <- 25
  # n_steps   <- 200
  # stop      <- F

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = NA
    , sl    = NA
  )

  # Simulate random steps
  for (i in 2:n_steps){

    # # For testing only
    # i <- 2

    # Draw random turning angles
    ta_new <- rvonmises(n_rsteps
      , mu    = ta_dist$params$mu
      , kappa = ta_dist$params$kappa
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$params$shape
      , scale = sl_dist$params$scale
    )

    # Make sure that the steps cover at least a minimal distance (this is
    # relevant if we need to compute the log of it)
    sl_new[sl_new < 0.0001] <- 0.0001

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        absta  = track$absta[i - 1] + ta_new
      , ta     = ta_new
      , sl     = sl_new
    )

    # We need to make sure that the absolute turning angle ranges from 0 to 2 *
    # pi
    rand$absta[rand$absta > 2 * pi] <-
      rand$absta[rand$absta > 2 * pi] - 2 * pi
    rand$absta[rand$absta < 0] <-
      rand$absta[rand$absta < 0] + 2 * pi

    # Calculate new endpoints
    rand$x <- track$x[i] + sin(rand$absta) * rand$sl
    rand$y <- track$y[i] + cos(rand$absta) * rand$sl

    # Create spatial points from endpoints
    coordinates(rand) <- c("x", "y")

    # Depending on the answer in the beginning, the loop breaks if one of the
    # new coordinates is outside the map boundaries
    if (stop){
      if (nrow(rand[ext, ]) != n_rsteps){
        break
      }
    } else {
      rand <- rand[ext, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)
    rand$xy <- NULL

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]

    # Interpolate coordinates and extract covariates
    extracted <- sapply(1:nrow(endcoords), function(x){
      line <- interpolatePoints(
          x1 = begincoords[1, 1]
        , x2 = endcoords[x, 1]
        , y1 = begincoords[1, 2]
        , y2 = endcoords[x, 2]
        , by = 0.1
      )
      extr <- raster::extract(covars, line)
      extr <- colMeans(extr)
      return(extr)
    })

    # Bind with data on random steps
    rand <- cbind(rand, t(extracted))

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, rand)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Sample one of the steps based on predicted probabilities
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y
  }

  # Assign step numbers
  track$step_number <- 0:(nrow(track) - 1)

  # Return track, yet remove initial pseudo-fix
  return(track[-1, ])
}

################################################################################
#### Function to Extract Covariates Along Step
################################################################################
# Function to extract covariates at interpolated locations and compute their
# average
extract_covariates_along_interpolated <- function(x, covariates, by = 0.1){
  extracted <- pbmclapply(1:nrow(x), ignore.interactive = T, function(i){
    int <- interpolatePoints(
        x1 = x$x1_[i]
      , x2 = x$x2_[i]
      , y1 = x$y1_[i]
      , y2 = x$y2_[i]
      , by = 0.1
    )
    extr <- raster::extract(covariates, int)
    extr <- colMeans(extr)
    return(extr)
  })
  extracted <- do.call(rbind, extracted)
  extracted <- as.data.frame(extracted)
  return(extracted)
}

################################################################################
#### Function to Extend a Covariate Layer
################################################################################
# Function that extends a covariate layer and fills the added border with values
# sampled from the layer
extendRaster <- function(x, y){
  extended <- lapply(1:nlayers(x), function(z) {
    layer <- x[[z]]
    na_mask <- is.na(layer)
    vals <- values(layer)
    vals <- na.omit(vals)
    r <- extend(layer, y)
    na_mask <- extend(na_mask, y, value = 0)
    indices <- which(is.na(values(r)))
    values(r)[indices] <- vals[runif(length(indices), min = 1, max = length(vals))]
    r <- mask(r, na_mask, maskvalue = 1, updatevalue = NA)
    return(r)
  })
  return(stack(extended))
}

################################################################################
#### Function to Quickly Rasterize Spatial Lines
################################################################################
# To rasterize and count the lines, you could use rasterize(tracks, heatmap, fun
# = "count"). However, this takes ages. Let's write a better function that uses
# the "spatstat" library
rasterizeSpatstat <- function(l, r){
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im
  for (y in 1:length(l)){
    line    <- as.psp(l[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
  }
  return(raster(summed))
}
