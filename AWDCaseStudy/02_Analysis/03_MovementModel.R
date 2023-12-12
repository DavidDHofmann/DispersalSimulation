################################################################################
#### Step Selection Function - Model Selection
################################################################################
# Description: Here, we'll refit the base model described in Hofmann (2021) to
# verify that its still the "best" model for the new data. Afterwards, we'll
# expand it to a proper movement model by allowing interactions between
# environmental covariates and movement metrics.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(parallel)     # To run stuff in parallel
library(corrplot)     # To plot correlations
library(cowplot)      # For nice plots
library(ggpubr)       # For nice plots
library(glmmTMB)      # For modelling
library(pbmcapply)    # For progress bar parallel

# Identify number of cores for parallel computing
cores <- detectCores()

################################################################################
#### Loading Data
################################################################################
# Load data
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

################################################################################
#### Scaling Data
################################################################################
# Scale the covariates
dat <- transform(dat
  , log_sl_             = scale(log_sl_)
  , cos_ta_             = scale(cos_ta_)
  , sl_                 = scale(sl_)
  , DistanceToWater     = scale(DistanceToWater)
  , SqrtDistanceToWater = scale(SqrtDistanceToWater)
  , Water               = scale(Water)
  , Trees               = scale(Trees)
  , Shrubs              = scale(Shrubs)
  , HumansBuff5000      = scale(HumansBuff5000)
)

# Prepare dataframe into which we store the values
scaling <- data.frame(
    ColumnName  = names(dat)
  , Center      = NA
  , Scale       = NA
)

# Extract the values from the data
for (i in 1:ncol(dat)){
  center <- attr(dat[, i], "scaled:center")
  scale <- attr(dat[, i], "scaled:scale")
  scaling$Center[i] <- ifelse(is.null(center), NA, center)
  scaling$Scale[i] <- ifelse(is.null(scale), NA, scale)
}

# Prepare a list of named vectors
scaling <- list(
      center = setNames(scaling[, 2], scaling$ColumnName)
    , scale  = setNames(scaling[, 3], scaling$ColumnName)
)

# Store the object to an RDS
write_rds(scaling, "03_Data/03_Results/99_Scaling.rds")

################################################################################
#### Base Model: Forward Model Selection for Main Effects
################################################################################
# Here, I want to reassure that the model used in Hofmann et al. 2021 is still
# the most parsimonious model. Hence, we'll propose the same covariates and see
# if they all remain in the full model
covars <- c(
    "Water"
  , "SqrtDistanceToWater"
  , "Trees"
  , "Shrubs"
  , "HumansBuff5000"
)

# Initate vector of selected covars
selected <- c()

# Initiate vector into which we store the corresponding models
model_sel <- list()

# Prepare loop counter
i <- 1

# Run loop until no more covariates are left for selection
while(length(covars) > 0){

  # Prepare model formula
  forformulas <- lapply(covars, function(x){
    c(selected, x)
  })

  # Prepare a vector that contains all covariates. This will allow us to easily
  # spot the covariates of a specific model
  covariates <- forformulas %>%
    do.call(rbind, .) %>%
    as.data.frame(.) %>%
    unite(., col = "Covariates", sep = ", ")

  # Prepare a tibble that we want to fill
  models <- tibble(
      ModelID     = 1:length(covars)
    , Covariates  = covariates
    , Formula     = lapply(forformulas, writeForm)
    , Model       = NA # Result of the model
    , AIC         = NA # AIC of the model
  )

  # Run all models in parallel
  models$Model <- mclapply(models$Formula, mc.cores = cores - 1, function(x){
    glmm_clogit(x, data = dat)
  })

  # Extract the AIC value of each model
  models$AIC <- sapply(models$Model, AIC)

  # Identify the model with the lowest AIC and put the respective covariate
  # into the "selected" vector
  best_covar <- na.omit(covars[models$AIC == min(models$AIC, na.rm = TRUE)])
  selected <- c(selected, best_covar)

  # Update which covariates remain for selection
  covars <- covars[!(covars %in% selected)]

  # Store the model results into the list
  model_sel[[i]] <- models

  # Update the index
  i <- i + 1

  # Print how many covariates are left
  cat(length(covars), "covariates left...\n")

}

################################################################################
#### Base Model: Identify the "Best" Model
################################################################################
# Put all relevant information of the model selection into a single tibble
models <- tibble(
    Covariates  = model_sel %>%
    lapply(., function(x){x$Covariates}) %>% do.call(rbind, .)
  , Formula     = model_sel %>%
    lapply(., function(x){x$Formula}) %>% do.call(c, .)
  , Model       = model_sel %>%
    lapply(., function(x){x$Model}) %>% do.call(c, .)
  , AIC         = model_sel %>%
    lapply(., function(x){x$AIC}) %>% do.call(c, .)
  , DeltaAIC    = NA # Difference in AIC to the "best" model
  , Weight      = NA # Weight assigned according to AIC
  , LogLik      = NA # Loglikelihood of the model
)

# Extract the LogLiks
models$LogLik <- lapply(models$Model, logLik) %>%
  do.call(rbind, .) %>%
  as.vector()

# Add a unique model id
models$ModelID <- 1:nrow(models)

# Calculate the DeltaAIC to the most parsimonious model
models$DeltaAIC <- models$AIC - min(models$AIC, na.rm = TRUE)

# Sort the tibble by the DeltaAIC
models <- arrange(models, DeltaAIC)

# Calculate the AIC weights for the models with deltaAIC <= 2. Lets first write
# a function to get the AIC weight of the best models
AICweight <- function(x){
  (exp(-0.5 * x)) / sum(exp(-0.5 * x))
}

# Split our tibble into models with DeltaAIC <= 2 as the threshold
threshold <- 2
models1 <- subset(models, DeltaAIC <= threshold)
models2 <- subset(models, !DeltaAIC <= threshold | is.na(DeltaAIC))

# Calculate the model weights. Models of the second group will all get weight 0
models1$Weight <- AICweight(models1$DeltaAIC)
models2$Weight <- 0

# Put the models together in a tibble again
rownames(models1$Covariates) <- 1:nrow(models1)
rownames(models2$Covariates) <- (1 + nrow(models1)):(nrow(models2) + nrow(models1))
models <- rbind(models1, models2)

# Look at the result
print(models, n = 100)

# Look at the results of the best model
summary(models$Model[[1]])

# Plot the model coefficients of the best model
showCoeffs(getCoeffs(models$Model[[1]])[-1, ], xlim = c(-1.5, 1.5))

# Let's prepare a nice table to report in our results
table <- cbind(
    models$Covariates
  , select(models, -c(Covariates, Formula, Model))
)

# Store the table
write_rds(table, "03_Data/03_Results/99_BaseModelAICs.rds")

# To save space when saving the model results we subset to only those models
# that get positive weight according to AIC
models <- subset(models, Weight > 0)

# Store the result for later
write_rds(models, "03_Data/03_Results/99_BaseModel.rds")

################################################################################
#### Movement Model: Forward Model Selection for Interaction Effects
################################################################################
# Let's retrieve all covariates from the best model
base <- models$Covariates[1, ] %>% strsplit(., ", ") %>% .[[1]]

# Now combine all of these covariates with the cos(ta_) and the log(sl_)
covars <- c(
    paste0(base, ":", "sl_")
  , paste0(base, ":", "log_sl_")
  , paste0(base, ":", "cos_ta_")
)

# We also want to test for an interaction between the step length and the
# activity phase. In addition, we may want to check for interactions between
# step lengths and turning angles
covars <- c(covars
  , "log_sl_:inactive"
  , "sl_:inactive"
  , "sl_:cos_ta_"
  , "log_sl_:cos_ta_"
)

# Look at the covariates that we are going to test for
covars

# Prepare a vector into which we put the covariates that are selected during our
# model selection loop
selected <- c()

# Also prepare a list into which we store all model results
model_sel <- list()

# Initiate a counter
i <- 1

# Run a loop that iteratively adds interaction terms until none are left for
# selection
while (length(covars) > 0){

  # Write down all of the covariates for the different models
  forformulas <- lapply(covars, function(x){
    c(base, selected, x)
  })

  # Prepare a vector that makes the covariates of a model slightly better
  # visible
  covariates <- forformulas %>%

    # Collapse the list
    do.call(rbind, .) %>%

    # Convert the matrix to a dataframe
    as.data.frame(.) %>%

    # Combine the columns into a vector that shows all covariates of a specific
    # model
    unite(., col = "Covariates", sep = ", ")

  # Prepare an empty tibble to fill
  models <- tibble(
      ModelID     = 1:length(covars)
    , Covariates  = covariates
    , Formula     = lapply(covars, function(x){

      # Write the model for the baseline covariates
      formula <- writeForm(base)

      # Check if we selected some interaction terms in the previous iteration
      if (length(selected) > 0){

        # If so, add them to the models as seperate terms (without random slope)
        toadd <- selected %>%

          # Transpose the matrix
          t(.) %>%

          # Convert the matrix to a dataframe
          as.data.frame(.) %>%

          # Combine the interactions using + as a seperator
          unite(., col, sep = " + ") %>%

          # Get rid of the dataframe structure
          .[["col"]]

        # Update the baseline model with the selected interaction terms
        formula <- update(formula, paste("~ . +", toadd))
      }

      # Add the interaction terms we still need to check
      formula <- update(formula, paste("~ . +", x))

      # Return the model formula that we just created
      return(formula)
    })
    , Model       = NA # Result of the model
    , AIC         = NA # AIC of the model
  )

  # Run all models
  models$Model <- mclapply(models$Formula, mc.cores = cores - 1, function(x){
    glmm_clogit(x, data = dat)
  })

  # Extract the AIC value of each model
  models$AIC <- lapply(models$Model, AIC) %>%
    do.call(rbind, .) %>%
    as.vector()

  # Identify the model with the lowest AIC and put the respective covariate
  # into the "selected" vector
  best_covar <- na.omit(covars[models$AIC == min(models$AIC, na.rm = TRUE)])
  selected <- c(selected, best_covar)

  # Update which covariates remain for selection
  covars <- covars[!(covars %in% selected)]

  # Store the model results into the list
  model_sel[[i]] <- models

  # Update the index
  i <- i + 1

  # Print how many covariates are left
  cat(length(covars), "covariates left...\n")
}

################################################################################
#### Movement Model: Identify the "Best" Model
################################################################################
# Put all relevant information into a single tibble
models <- tibble(
    Covariates  = model_sel %>%
    lapply(., function(x){x$Covariates}) %>% do.call(rbind, .)
  , Formula     = model_sel %>%
    lapply(., function(x){x$Formula}) %>% do.call(c, .)
  , Model       = model_sel %>%
    lapply(., function(x){x$Model}) %>% do.call(c, .)
  , AIC         = model_sel %>%
    lapply(., function(x){x$AIC}) %>% do.call(c, .)
  , DeltaAIC    = NA # Difference in AIC to the "best" model
  , Weight      = NA # Weight assigned according to AIC
  , LogLik      = NA # Loglikelihood of the model
)

# Extract the LogLiks
models$LogLik <- lapply(models$Model, logLik) %>%
  do.call(rbind, .) %>%
  as.vector()

# Add a unique model id
models$ModelID <- 1:nrow(models)

# Calculate the DeltaAIC to the most parsimonious model
models$DeltaAIC <- models$AIC - min(models$AIC, na.rm = TRUE)

# Sort the tibble by the DeltaAIC
models <- arrange(models, DeltaAIC)

# Calculate the AIC weights for the models with deltaAIC <= 2. Lets first write
# a function to get the AIC weight of the best models
AICweight <- function(x){
  (exp(-0.5 * x)) / sum(exp(-0.5 * x))
}

# Split our tibble into models with DeltaAIC <= 2 as the threshold
threshold <- 2
models1 <- subset(models, DeltaAIC <= threshold)
models2 <- subset(models, !DeltaAIC <= threshold | is.na(DeltaAIC))

# Calculate the model weights. Models of the second group will all get weight 0
models1$Weight <- AICweight(models1$DeltaAIC)
models2$Weight <- 0

# Put the models together in a tibble again
rownames(models1$Covariates) <- 1:nrow(models1)
rownames(models2$Covariates) <- (1 + nrow(models1)):(nrow(models2) + nrow(models1))
models <- rbind(models1, models2)

# Look at the result
print(models, n = 100)

# Look at the results of the best model
summary(models$Model[[1]])

# Plot the model coefficients of the best model
showCoeffs(getCoeffs(models$Model[[1]])[-1, ], xlim = c(-1.5, 1.5))

# Let's prepare a nice table to report in our results
table <- cbind(
    models$Covariates
  , select(models, -c(Covariates, Formula, Model))
)

# Store the table
write_rds(table, "03_Data/03_Results/99_MovementModelAICs.rds")

# To save space when saving the model results we subset to only those models
# that get positive weight according to AIC
models <- subset(models, Weight > 0)

# Store the result for later
write_rds(models, "03_Data/03_Results/99_MovementModel.rds")

# ################################################################################
# #### Model Averaging
# ################################################################################
# # Check the paper by Symonds and Moussalli (2011) for details on the formulas
# # In our case model averaging isn't really necessary since the first model gets
# # a weight of 1 anyways. However, this might still be a useful exercise. First
# # we subset to the models with a postive weight (only one in our case)
# predis <- subset(models, Weight > 0)
#
# # From each model we need to extract the coefficients
# predis <- predis %>% mutate(Coeffs = map(Model, function(x){
#   getCoeffs(x)
# }))
#
# # Add the weights to the coefficient tables. Let's also use this step to add the
# # model ID into each dataframe
# for (i in 1:nrow(predis)){
#   predis$Coeffs[[i]]$Weight <- predis$Weight[i]
#   predis$Coeffs[[i]]$ModelID <- predis$ModelID[i]
# }
#
# # Bind all dataframes of coefficients together
# predis <- do.call(rbind, predis$Coeffs)
#
# # Make sure that the weights of the covariates from the base model add up to 1
# predis %>%
#   group_by(Covariate) %>%
#   summarize(TotalWeight = sum(Weight))
#
# # Calculate the weighted coefficients
# predis <- predis %>% mutate(WeightedCoefficient = Coefficient * Weight)
#
# # Sum up weighted coefficients
# avg_coeffs <- predis %>%
#   group_by(Covariate) %>%
#   summarize(AveragedCoefficient = sum(WeightedCoefficient))
#
# # We still need to calculate the averaged variance. To do so we firsts put the
# # averaged coefficients back into the predis table
# predis <- left_join(predis, avg_coeffs, by = "Covariate")
#
# # We can now calculate the weighted SEs
# predis <- predis %>%
#   mutate(WeightedSE = Weight * sqrt(
#     SE ** 2 + (Coefficient - AveragedCoefficient) ** 2
#   ))
#
# # Sum up the weighted SEs
# avg_ses <- predis %>%
#   group_by(Covariate) %>%
#   summarize(AveragedSEs = sum(WeightedSE))
#
# # Put everything into one table
# coeffs <- full_join(avg_coeffs, avg_ses) %>%
#   rename(Coefficient = AveragedCoefficient, SE = AveragedSEs)
#
# # Look at the results
# showCoeffs(coeffs[-1, ])

################################################################################
#### Model Validation: Necessary Functions
################################################################################
# Finally I validate the best model using the procedure described in Fortin et
# al. 2009. Let's reload the model results
models <- read_rds("03_Data/03_Results/99_MovementModel.rds")

# Set a seed for reproducibility
set.seed(1234)

# Write a function that allows us to run the validation process
crossVal <- function(model, data, ratio, random = FALSE){

  # Identify all unique step IDs
  step_ids <- unique(data$step_id_)

  # Split IDs according to ratio
  sampled <- sample(step_ids, size = length(step_ids) * ratio)

  # Split the data into training and validation datasets
  ssf_train <- subset(dat, step_id_ %in% sampled)
  ssf_valid <- subset(dat, !step_id_ %in% sampled)

  # Train the model using the training data
  formula <- formula(model$call)
  mod <- glmm_clogit(formula, data = ssf_train)

  # Use our "prepareModel" function to get the model coefficients and a cleaned
  # model formula
  mod <- prepareModel(mod)

  # Prepare model matrix for validation data
  modeldat <- model.matrix(mod$formula, data = ssf_valid)

  # Multiply the matrix with estimated selection coefficients to calculate
  # predicted selection scores
  ssf_valid$Scores <- exp(modeldat %*% mod$coefficients)

  # Depending on whether we want to randomize preferences, we include or exclude
  # realized steps
  if (random){

    # In case we want to randomize, get rid of case_ steps
    validation <- subset(ssf_valid, !case_) %>%

      # Group by cluster
      group_by(., step_id_) %>%

      # Identify the rank of each step
      mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%

      # Randomly select one of the steps per cluster
      sample_n(1)

  } else {

    # In case we dont want to randomize, we keep all steps
    validation <- ssf_valid %>%

      # Group by cluster
      group_by(., step_id_) %>%

      # Identify the rank of each step
      mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%

      # Subset to realized steps only
      subset(., case_)

  }

  # Now do some transformations
  validation <- validation %>%
    ungroup(.) %>%
    select(., Rank) %>%
    table(.) %>%
    as.data.frame(.) %>%
    set_names(., c("Rank", "Frequency")) %>%
    mutate(., Rank = as.numeric(Rank))

  # Run a pearson rank corrleation test
  spearman <- cor.test(
      x       = validation$Rank
    , y       = validation$Frequency
    , method  = "spearman"
    , exact   = FALSE
  )

  # Return the spearman test result as well as the data
  return(
    list(
        "Data"                = validation
      , "SpearmanCorrelation" = as.vector(spearman$estimate)
    )
  )
}

################################################################################
#### Validation
################################################################################
# Set seed for multicore processes
set.seed(1234)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

# Run validation repeadedly with observed preferences
valid_pref <- pbmclapply(
    X                  = 1:100
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){
    success <- F
    while (!success){
      output <- tryCatch(crossVal(
          model   = models$Model[[1]]
        , data    = dat
        , ratio   = 0.8
        , random  = FALSE
      ), warning = function(w){return("skip")})
      success <- ifelse (class(output) == "list", T, F)
    }
    return(output)
})

# Extract all the rhos
rho_pref <- lapply(valid_pref, function(x){
  x[[2]]
}) %>% unlist()

# Extract all the data
dat_pref <- lapply(valid_pref, function(x){
  x[[1]]
}) %>% do.call(rbind, .)

# Run validation repeadedly with random preferences
valid_rand <- pbmclapply(
    X                  = 1:100
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){
    success <- F
    while (!success){
      output <- tryCatch(crossVal(
          model   = models$Model[[1]]
        , data    = dat
        , ratio   = 0.8
        , random  = TRUE
      ), warning = function(w){return("skip")})
      success <- ifelse (class(output) == "list", T, F)
    }
    return(output)
})

# Extract all the rhos
rho_rand <- lapply(valid_rand, function(x){
  x[[2]]
}) %>% unlist()

# Extract all the data
dat_rand <- lapply(valid_rand, function(x){
  x[[1]]
}) %>% do.call(rbind, .)

################################################################################
#### Validation: Testing
################################################################################
# Find the means and confidence intervals of the rho's from the two datasets
test_pref <- t.test(rho_pref)
test_rand <- t.test(rho_rand)

# Put everything together into a single dataframe
validation <- data.frame(
    Group = c("rho_pref", "rho_rand")
  , Mean  = c(test_pref$estimate, test_rand$estimate)
  , LCL   = c(test_pref$conf.int[1], test_rand$conf.int[1])
  , UCL   = c(test_pref$conf.int[2], test_rand$conf.int[2])
)

# Put the results to file
write_rds(
    validation
  , "03_Data/03_Results/99_ModelValidation.rds"
)
write_rds(
    list(dat_pref, dat_rand)
  , "03_Data/03_Results/99_ModelValidation(Data).rds"
)
