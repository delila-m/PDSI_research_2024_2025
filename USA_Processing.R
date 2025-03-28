# working directory for Geology computer
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024")

# load in necessary libraries
source("drought_functions.R")


# load in PDSI data for whole country
pdsi <- rast("agg_met_pdsi_1979_CurrentYear_CONUS.nc")

# load in FIPS codes and filter for only continuous US, 
usa_fips <- read.csv("All Counties and FIPS Codes.csv")
contUS <- filter(usa_fips, ! State %in% c("AK","HI","PR"))
# make sure each fips code is 5 digits
contUS <- contUS %>%
  mutate(`AOI.Value` = sprintf("%05d", `AOI.Value`))
# fix naming so each only includes the actual name and not suffixes like 'county'
contUS <- contUS %>%
  mutate(
    AreaClean = case_when(
      # Handle special case
      Area == "District of Columbia" ~ "Columbia",
      # Handle common suffixes with regex to ensure exact matching at the end
      str_detect(Area, " County$") ~ str_replace(Area, " County$", ""),
      str_detect(Area, " Parish$") ~ str_replace(Area, " Parish$", ""),
      str_detect(Area, " City$") ~ str_replace(Area, " City$", ""),
      # Default case - keep original
      TRUE ~ Area), 
    AreaClean = str_trim(AreaClean)
  )

# # initialize data frame
# usa.data <- data.frame()
# 
# # loop through each county
# for (index in 1:nrow(contUS)){
#   county.name <- contUS$AreaClean[index]
#   state.name <- contUS$State[index]
#   fips <- contUS$AOI.Value[index]
# 
#   # get the file name
#   file.name <- paste0("countyData/USDM-", fips, ".csv")
# 
#   # read in the data
#   drought.data <- read.csv(file.name)
# 
#   # processing
#   clean.data <- clean.county.data(state.name, county.name, pdsi, drought.data, TRUE)
# 
#   # add state and county column for identification
#   clean.data <- clean.data %>% mutate(State = state.name,
#                                       County = county.name)
# 
#   # add all of the data together
#   usa.data <- rbind(usa.data, clean.data)
# }
# 
# 
# # got through 4 counties in Montana, then quit bc the dataset was too big
# filtered.usa <- usa.data %>% filter(State != "MT")
# binned.usa <- bin.lat.long(filtered.usa, 0.25)
# 


#####################################
# cleaning the US county data

states <- unique(contUS$State)
state.fips <- vector(mode = "list",length = length(states))
allStatesDf25 <- data.frame()
# loop through continental US states
for(index in seq_along(states)){
  
  # filter for just the state's data
  #state.fips[[index]] <- contUS %>% filter(State == states[index])
  state <- readRDS(paste0(states[index], " .RDS"))
  allStatesDf25 <- rbind(allStatesDf25, state)
}

saveRDS(allStatesDf25, file = "CleanedUS_0.25.RDS")


# Run in furrr ------------------------------------------------------------
library(future)
library(furrr)
future::plan(strategy = multisession, workers = 6)

allStates <- furrr::future_map(state.fips,bin.state.data,TRUE, 1, .progress = TRUE)

allStatesDf <- list_rbind(allStates)

saveRDS(allStatesDf, file = "CleanedUS_1.RDS")

######
# loading in data
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")

allStatesDf <- readRDS("CleanedUS_1.RDS")

# # filter for west-coast data only 
# westUS <- allStatesDf %>% filter(bin.x <= -100)
# # take sample fraction bc it's still too big 
# subsetWest <- westUS %>% sample_frac(0.50)

# rf modeling 
# split into training and testing 
train <- allStatesDf %>% sample_frac(0.80)
test <- anti_join(allStatesDf , train)


# using ranger
library(ranger)

rf.ranger.fit <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
  )

# perhaps three iterations of tests of the training set?
importance <- rf.ranger.fit$variable.importance.local
### trying to make sense of importance values 
# Check the dimensions and structure
dim(rf.ranger.fit$variable.importance.local)
str(rf.ranger.fit$variable.importance.local)

# Check the number of rows in your training data
nrow(train)

# Compare with the length of the importance object
length(rf.ranger.fit$variable.importance.local)


# predict 
preds.ranger <- predict(rf.ranger.fit, test)


# Add predictions to your test dataset
test$predicted <- preds.ranger$predictions


# general RMSE, and rmse for each observation 
RMSE(preds.ranger$predictions, test$USDM_Avg) # 0.5687996!!!

# calculate absolute error and squared error for each observation
test.binned <- test %>% group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            .groups = 'drop')


# save fit object, traing set, testing set, and predictions
save(rf.ranger.fit, train, test, preds.ranger, file = "RFAnalysis1.Rdata")

################ # plot some stuff

# load in the rf objects
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")
load("RFAnalysis1.Rdata")

library(viridis)  
library(ggplot2)
library(metR)


ggplot(test.binned, aes(x = bin.x, y = bin.y, z = MSE)) +
  geom_contour_fill() +  # Interpolated filled contours
  scale_fill_viridis(name = "Absolute Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude", 
       y = "Latitude")


# Option 1: Using geom_tile() for a heatmap
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# Option 2: Using geom_raster() which can be faster for large datasets
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = MSE)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")


###### Trying the model again with elevation as a predictor


####### trying model again without 2024..





####################################
# trying again at 0.25
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_0.25")


allStatesDf25 <- readRDS("CleanedUS_0.25.RDS")

train25 <- allStatesDf25 %>% sample_frac(0.80)
test25 <- anti_join(allStatesDf25, train25)


rf.ranger.fit.25 <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train25,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
)

# perhaps three iterations of tests of the training set?
importance25 <- rf.ranger.fit.25$variable.importance.local

# predict and find RMSE
preds.ranger.25 <- predict(rf.ranger.fit.25, test25)



RMSE(preds.ranger.25$predictions, test25$USDM_Avg)
