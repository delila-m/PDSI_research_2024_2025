# Libraries 
library(ranger)
library(viridis)  
library(ggplot2)
library(metR)
library(tigris)
library(patchwork)
library(tidyverse)
library(dplyr)
library(xgboost)
library(caret)
library(sf)
##
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024")
# load in necessary libraries
source("drought_functions.R")



# Loading in Data, cleaning, and binning to nearest degree of lat/long
#####

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
  state.fips[[index]] <- contUS %>% filter(State == states[index])
  
}

saveRDS(allStatesDf25, file = "CleanedUS_0.25.RDS")


# Run in furrr ------------------------------------------------------------
library(future)
library(furrr)
future::plan(strategy = multisession, workers = 6)

allStates <- furrr::future_map(state.fips,bin.state.data,TRUE, 0.5, .progress = TRUE)

allStatesDf_0.5 <- list_rbind(allStates)
# save cleaned data binned to nearest degree
saveRDS(allStatesDf_0.5, file = "CleanedUS_0.5.RDS")
#####


### RF Modeling

# Running RF model on data binned to one degree
######
# loading in data
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")

allStatesDf <- readRDS("CleanedUS_1.RDS")

# split into training and testing 
train <- allStatesDf %>% sample_frac(0.80)
test <- anti_join(allStatesDf , train)



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
importance <- data.frame(rf.ranger.fit$variable.importance.local)
train_ungrouped <- train %>% ungroup()

# Then add the importance columns
train.importance <- train_ungrouped %>% 
  mutate(
    x.bin.importance = importance$bin.x, 
    y.bin.importance = importance$bin.y, 
    PDSI.importance = importance$PDSI_Avg
  )

# predict 
preds.ranger <- predict(rf.ranger.fit, test)


# Add predictions to your test dataset
test$predicted <- preds.ranger$predictions


# general RMSE, and rmse for each observation 
RMSE(preds.ranger$predictions, test$USDM_Avg) # 0.5687996!!!


## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test.binned <- test %>% group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            .groups = 'drop')
# bin training set by lat/long to look at importance
train.binned <- train.importance %>% group_by(bin.x, bin.y) %>% 
  summarise(Mean.binx.importance = mean(x.bin.importance), 
            Mean.biny.importance = mean(y.bin.importance),
            Mean.PDSI.importance = mean(PDSI.importance), 
            .groups = 'drop')




# save fit object, traing set, testing set, and predictions
save(rf.ranger.fit, train, test, preds.ranger, file = "RFAnalysis1.Rdata")

### plot some stuff

# load in the rf objects
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")
load("RFAnalysis1.Rdata")



#### Plotting MSE
ggplot(test.binned, aes(x = bin.x, y = bin.y, z = MSE)) +
  geom_contour_fill() +  # Interpolated filled contours
  scale_fill_viridis(name = "Absolute Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude", 
       y = "Latitude")


# Option 1: Using geom_tile() for a heatmap
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = MSE)) +
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


##### Plotting variable importance 
ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean PDSI Importance") +
  theme_minimal() +
  labs(title = "PDSI Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.binx.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Longitude Importance") +
  theme_minimal() +
  labs(title = "Longitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.biny.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Latitude Importance") +
  theme_minimal() +
  labs(title = "Latitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")
#####


###### Trying the model again with elevation as a predictor
######
# loading in data
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")
allStatesDf <- readRDS("CleanedUS_1.RDS")
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/")


# Get continental US boundary using tigris
states <- states(cb = TRUE)
continental_us <- states[!states$STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "AS", "MP"), ]
continental_us_vect <- vect(continental_us)

# Get elevation data (adjust z for resolution)
elevation <- get_elev_raster(continental_us, z = 8, src = "aws")
# Note: Lower z value for initial test, then increase for more detail

# Convert to terra raster
elevation_terra <- rast(elevation)

plot(elevation_terra)

# Save as GeoTIFF
writeRaster(elevation_terra, "us_elevation.tif", overwrite=TRUE)

elevation_terra <- rast("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/us_elevation.tif")

# extract raster into a regular dataframe
elev.df <- extract(x = elevation_terra, continental_us_vect, xy = TRUE)
# filter out noise- low points in lakes and two random high points which likely got extracted wrong 
# filtering by the highest and lowest points in the continental US
elev.cleaned <- elev.df %>% 
  filter(file5b885a507540 >= -86 & file5b885a507540 <= 4421) %>% 
  mutate(Elevation_m = file5b885a507540) %>%
  select(ID, x, y, Elevation_m)
# save those guys
saveRDS(elev.df, file = "rawElevDataContUS.RDS")
saveRDS(elev.cleaned, file = "cleanedElevDataContUS.RDS")

# aggregate to .25 degree
binned.elev <- bin.elevation.data(elev.cleaned, 1)

# combine elevation data with usdm/pdsi data
allStatesElevDf_1 <- inner_join(allStatesDf, binned.elev, 
                                by = c("bin.x" = "bin.x", 
                                       "bin.y" = "bin.y")) %>% 
  select(c("bin.x", "bin.y", "Date", "PDSI_Avg", 
           "USDM_Avg", "Elev_Avg"))


# Now we can train the model with elevation 
# split into training and testing 
train <- allStatesElevDf_1 %>% sample_frac(0.80)
test <- anti_join(allStatesElevDf_1 , train)

rf.ranger.fit <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + Elev_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
)
# predict 
preds.ranger <- predict(rf.ranger.fit, test)

# Add predictions to your test dataset
test$predicted <- preds.ranger$predictions

# save it all for later use 
save(rf.ranger.fit, train, test, preds.ranger, file = "RFAnalysisElev1.Rdata")


# general RMSE, and rmse for each observation 
RMSE(preds.ranger$predictions, test$USDM_Avg) # 0.5613799!!!


## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test.binned <- test %>% 
  mutate(Abs_Error = abs(USDM_Avg - predicted), 
         Sq_Error = (USDM_Avg - predicted)^2) %>% 
  group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            .groups = 'drop')

# plotting MSE
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# perhaps three iterations of tests of the training set?
importance <- data.frame(rf.ranger.fit$variable.importance.local)
train_ungrouped <- train %>% ungroup()

# Then add the importance columns
train.importance <- train_ungrouped %>% 
  mutate(
    x.bin.importance = importance$bin.x, 
    y.bin.importance = importance$bin.y, 
    PDSI.importance = importance$PDSI_Avg,
    elev.importance = importance$Elev_Avg
  )

# bin training set by lat/long to look at importance
train.binned <- train.importance %>% group_by(bin.x, bin.y) %>% 
  summarise(Mean.binx.importance = mean(x.bin.importance), 
            Mean.biny.importance = mean(y.bin.importance),
            Mean.PDSI.importance = mean(PDSI.importance), 
            Mean.elev.importance = mean(elev.importance),
            .groups = 'drop')

### Plotting variable importance 
ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean PDSI Importance") +
  theme_minimal() +
  labs(title = "PDSI Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.binx.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Longitude Importance") +
  theme_minimal() +
  labs(title = "Longitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.biny.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Latitude Importance") +
  theme_minimal() +
  labs(title = "Latitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.elev.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Elevation Importance") +
  theme_minimal() +
  labs(title = "Elevation Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")
#####


# trying again with elevation but leaving out 20% of years, from 2020-2025
#####
# Extract the year from Date column
allStatesElevDf_1$Year <- as.numeric(format(allStatesElevDf_1$Date, "%Y"))

# Get unique years and sort them
years <- sort(unique(allStatesElevDf_1$Year))

# Calculate how many years make up approximately 20%
num_test_years <- ceiling(length(years) * 0.2)

# Use the most recent consecutive years as test set
test_years <- tail(years, num_test_years)
train_years <- head(years, length(years) - num_test_years)

# Create train and test datasets
train <- allStatesElevDf_1[allStatesElevDf_1$Year %in% train_years, ]
test <- allStatesElevDf_1[allStatesElevDf_1$Year %in% test_years, ]

rf.ranger.fit.yearsplit <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + Elev_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
)
# predict 
preds.ranger <- predict(rf.ranger.fit.yearsplit, test)

# save it all for later use 
save(rf.ranger.fit.yearsplit, train, test, preds.ranger, file = "RFAnalysisElev1YearSplit.Rdata")


# general RMSE, and rmse for each observation 
RMSE(preds.ranger$predictions, test$USDM_Avg) # 0.6284244!!!

## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test.binned <- test %>% 
  mutate(Abs_Error = abs(USDM_Avg - predicted), 
         Sq_Error = (USDM_Avg - predicted)^2) %>% 
  group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            .groups = 'drop')

# plotting MSE
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# perhaps three iterations of tests of the training set?
importance <- data.frame(rf.ranger.fit.yearsplit$variable.importance.local)
train_ungrouped <- train %>% ungroup()

# Then add the importance columns
train.importance <- train_ungrouped %>% 
  mutate(
    x.bin.importance = importance$bin.x, 
    y.bin.importance = importance$bin.y, 
    PDSI.importance = importance$PDSI_Avg,
    elev.importance = importance$Elev_Avg
  )


# bin training set by lat/long to look at importance
train.binned <- train.importance %>% group_by(bin.x, bin.y) %>% 
  summarise(Mean.binx.importance = mean(x.bin.importance), 
            Mean.biny.importance = mean(y.bin.importance),
            Mean.PDSI.importance = mean(PDSI.importance), 
            Mean.elev.importance = mean(elev.importance),
            .groups = 'drop')

### Plotting variable importance 
ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean PDSI Importance") +
  theme_minimal() +
  labs(title = "PDSI Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.binx.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Longitude Importance") +
  theme_minimal() +
  labs(title = "Longitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.biny.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Latitude Importance") +
  theme_minimal() +
  labs(title = "Latitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.elev.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Elevation Importance") +
  theme_minimal() +
  labs(title = "Elevation Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# plot 2020 predictions vs actual

# predict only 2020
test.2020 <- test %>% filter(Year == 2024)
preds.2020 <- predict(rf.ranger.fit.yearsplit, test.2020)
# add to test df for easy plotting 
test.2020$preds <- preds.2020$predictions

# create predicted and actual plots, add together and plot
pred <- ggplot(test.2020, aes(x = bin.x, y = bin.y, fill = preds)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Predicted Average USDM Across Continental US for 2024",
       x = "Longitude Bin", 
       y = "Latitude Bin")

actual <-  ggplot(test.2020, aes(x = bin.x, y = bin.y, fill = USDM_Avg)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Actual Average USDM Across Continental US for 2024",
       x = "Longitude Bin", 
       y = "Latitude Bin")

plot <- pred + actual
plot

#####


# trying to bin data by 0.5 degree 
#####
# loading in data
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_0.5")

allStatesDf_0.5 <- readRDS("CleanedUS_0.5.RDS")

# split into training and testing 
train_0.5 <- allStatesDf_0.5 %>% sample_frac(0.80)
test_0.5 <- anti_join(allStatesDf_0.5 , train)

# fit da model 
rf.ranger.fit.0.5 <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train_0.5,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
)

# grab variable importance  
importance_0.5 <- data.frame(rf.ranger.fit.0.5$variable.importance.local)
train_ungrouped <- train_0.5 %>% ungroup()

# Then add the importance columns
train.importance.0.5 <- train_ungrouped %>% 
  mutate(
    x.bin.importance = importance_0.5$bin.x, 
    y.bin.importance = importance_0.5$bin.y, 
    PDSI.importance = importance_0.5$PDSI_Avg
  )

# predict 
preds.ranger.0.5 <- predict(rf.ranger.fit.0.5, test_0.5, type = "response")

# Add predictions to your test dataset
test_0.5$predicted <- preds.ranger.0.5$predictions


# general RMSE, and rmse for each observation 
RMSE(preds.ranger.0.5$predictions, test_0.5$USDM_Avg) # 0.4755648!!!


## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test_0.5.binned <- test_0.5 %>% 
  mutate(Abs_Error = abs(USDM_Avg - predicted), 
         Sq_Error = (USDM_Avg - predicted)^2) %>%
  group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            Mean_USDM = mean(USDM_Avg),
            Mean_predicted = mean(predicted),
            .groups = 'drop')
# bin training set by lat/long to look at importance
train.binned.0.5 <- train.importance.0.5 %>% group_by(bin.x, bin.y) %>% 
  summarise(Mean.binx.importance = mean(x.bin.importance), 
            Mean.biny.importance = mean(y.bin.importance),
            Mean.PDSI.importance = mean(PDSI.importance), 
            .groups = 'drop')

# save fit object, traing set, testing set, and predictions
save(rf.ranger.fit.0.5, train_0.5, test_0.5, preds.ranger.0.5, file = "RFAnalysis0.5.Rdata")

# load object in again
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_0.5")

load("RFAnalysis0.5.Rdata")

# plotting MSE
ggplot(test_0.5.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

### Plotting variable importance 
ggplot(train.binned.0.5, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean PDSI Importance") +
  theme_minimal() +
  labs(title = "PDSI Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned.0.5, aes(x = bin.x, y = bin.y, fill = Mean.binx.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Longitude Importance") +
  theme_minimal() +
  labs(title = "Longitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned.0.5, aes(x = bin.x, y = bin.y, fill = Mean.biny.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Latitude Importance") +
  theme_minimal() +
  labs(title = "Latitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# creating plot of predicted vs actual USDM rating 
# create predicted and actual plots, add together and plot
pred <- ggplot(test_0.5.binned, aes(x = bin.x, y = bin.y, fill = Mean_predicted)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Predicted Average USDM Across Continental US for Testing set",
       x = "Longitude Bin", 
       y = "Latitude Bin")

actual <-  ggplot(test_0.5.binned, aes(x = bin.x, y = bin.y, fill = Mean_USDM)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Actual Average USDM Across Continental US for Testing Set",
       x = "Longitude Bin", 
       y = "Latitude Bin")

plot <- pred + actual
plot

#####


# leave out years with data binned to nearest 0.5 degree
#####

setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_0.5")

allStatesDf_0.5 <- readRDS("CleanedUS_0.5.RDS")

# Extract the year from Date column
allStatesDf_0.5$Year <- as.numeric(format(allStatesDf_0.5$Date, "%Y"))

# Get unique years and sort them
years <- sort(unique(allStatesDf_0.5$Year))

# Calculate how many years make up approximately 20%
num_test_years <- ceiling(length(years) * 0.2)

# Use the most recent consecutive years as test set
test_years <- tail(years, num_test_years)
train_years <- head(years, length(years) - num_test_years)

# Create train and test datasets
train.yearsplit.0.5 <- allStatesDf_0.5[allStatesDf_0.5$Year %in% train_years, ]
test.yearsplit.0.5 <- allStatesDf_0.5[allStatesDf_0.5$Year %in% test_years, ]

rf.ranger.fit.yearsplit.0.5 <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train.yearsplit.0.5,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  importance = 'permutation',  # To measure feature importance
  verbose = TRUE, 
  local.importance = TRUE,
  quantreg = TRUE
)
# predict and add predictions to testing set
preds.ranger.yearsplit.0.5 <- predict(rf.ranger.fit.yearsplit.0.5, test.yearsplit.0.5)

test.yearsplit.0.5$predictions <- preds.ranger.yearsplit.0.5$predictions

# save it all for later use 
save(rf.ranger.fit.yearsplit.0.5, train.yearsplit.0.5, test.yearsplit.0.5, preds.ranger.yearsplit.0.5, file = "RFAnalysisElev1YearSplit.Rdata")


# general RMSE, and rmse for each observation 
RMSE(preds.ranger.yearsplit.0.5$predictions, test.yearsplit.0.5$USDM_Avg) # 0.6701104!!!

## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test.binned <- test.yearsplit.0.5 %>% 
  mutate(Abs_Error = abs(USDM_Avg - predictions), 
         Sq_Error = (USDM_Avg - predictions)^2) %>% 
  group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            .groups = 'drop')

# plotting MSE
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# perhaps three iterations of tests of the training set?
importance <- data.frame(rf.ranger.fit.yearsplit.0.5$variable.importance.local)
train_ungrouped <- train.yearsplit.0.5 %>% ungroup()

# Then add the importance columns
train.importance <- train_ungrouped %>% 
  mutate(
    x.bin.importance = importance$bin.x, 
    y.bin.importance = importance$bin.y, 
    PDSI.importance = importance$PDSI_Avg,
  )


# bin training set by lat/long to look at importance
train.binned <- train.importance %>% group_by(bin.x, bin.y) %>% 
  summarise(Mean.binx.importance = mean(x.bin.importance), 
            Mean.biny.importance = mean(y.bin.importance),
            Mean.PDSI.importance = mean(PDSI.importance), 
            .groups = 'drop')

### Plotting variable importance 
ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean PDSI Importance") +
  theme_minimal() +
  labs(title = "PDSI Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.binx.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Longitude Importance") +
  theme_minimal() +
  labs(title = "Longitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.biny.importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Latitude Importance") +
  theme_minimal() +
  labs(title = "Latitude Importance Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")

# plot 2020 predictions vs actual

# predict only 2020
test.2020 <- test.yearsplit.0.5 %>% filter(Year == 2020)
preds.2020 <- predict(rf.ranger.fit.yearsplit.0.5, test.2020)
# add to test df for easy plotting 
test.2020$preds <- preds.2020$predictions

# create predicted and actual plots, add together and plot
pred <- ggplot(test.2020, aes(x = bin.x, y = bin.y, fill = preds)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Predicted Average USDM Across Continental US for 2020",
       x = "Longitude Bin", 
       y = "Latitude Bin")

actual <-  ggplot(test.2020, aes(x = bin.x, y = bin.y, fill = USDM_Avg)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Actual Average USDM Across Continental US for 2020",
       x = "Longitude Bin", 
       y = "Latitude Bin")

plot <- pred + actual
plot

#####


### XGBoost Modeling
# trying xgboost at 1 degree
#####
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")

allStatesDf_1 <- readRDS("CleanedUS_1.RDS")

# split into training and testing 
train_1 <- allStatesDf_1 %>% sample_frac(0.80)
test_1 <- anti_join(allStatesDf_1, train_1)

train_x = data.matrix(train_1 %>% select(-c(USDM_Avg, x_Avg, y_Avg, Date)))
train_y = train_1 %>% select(USDM_Avg)
train_y <- train_y %>% select(-c(bin.x, bin.y))

#define predictor and response variables in testing set
test_x = data.matrix(test_1 %>% select(-c(USDM_Avg, x_Avg, y_Avg, Date)))
test_y = test_1 %>% select(USDM_Avg)

#define final training and testing sets
xgb_train = xgb.DMatrix(data = train_x, label = train_y$USDM_Avg)
xgb_test = xgb.DMatrix(data = test_x, label = test_y$USDM_Avg)

#defining a watchlist
watchlist = list(train=xgb_train, test=xgb_test)

#fit XGBoost model and display training and testing data at each iteartion
model = xgb.train(data = xgb_train, max.depth = 3, 
                  watchlist=watchlist, nrounds = 500, 
                  verbose = TRUE)



# different approach using caret
grid_tune <- expand_grid(nrounds = c(500, 1000, 1500), 
                         max_depth = c(2,4,6), 
                         eta = 0.3, 
                         gamma = 0, 
                         colsample_bytree = 1, 
                         min_child_weight = 1, 
                         subsample = 1)

train_ctrl <- trainControl(method = "cv", 
                           number = 3, 
                           verboseIter = TRUE, 
                           allowParallel = TRUE)

xgb_tune <-  train(x = train_x, 
                   y = train_y$USDM_Avg, 
                   trControl = train_ctrl, 
                   tuneGrid = grid_tune, 
                   method = "xgbTree", 
                   verbose = TRUE)

# look at the best model 
xgb_tune$best
# predict
xgb.preds <- predict(xgb_tune, test_1)

test_1$predicted <- xgb.preds
save(xgb_tune, train_1, test_1, xgb.preds, file = "XGBAnalysis_1.RData")

# general RMSE, and rmse for each observation 
RMSE(xgb.preds, test_1$USDM_Avg) # 0.5602002!!!

## bin training and testing sets based off lat/long values 

# calculate absolute error and squared error for each observation
test.binned <- test_1 %>% 
  mutate(Abs_Error = abs(USDM_Avg - predicted), 
         Sq_Error = (USDM_Avg - predicted)^2) %>% 
  group_by(bin.x, bin.y) %>% 
  summarise(MAE = mean(Abs_Error), 
            MSE = mean(Sq_Error),
            RMSE = sqrt(MSE), 
            Mean_predicted = mean(predicted), 
            Mean_USDM = mean(USDM_Avg),
            .groups = 'drop')

# plotting MSE
ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = RMSE)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Squared Error") +
  theme_minimal() +
  labs(title = "Prediction Error Across Continental US",
       x = "Longitude Bin", 
       y = "Latitude Bin")


# creating plot of predicted vs actual USDM rating 
# create predicted and actual plots, add together and plot
pred <- ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = Mean_predicted)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Predicted Average USDM Across Continental US for Testing set",
       x = "Longitude Bin", 
       y = "Latitude Bin")

actual <-  ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = Mean_USDM)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM") +
  theme_minimal() +
  labs(title = "Actual Average USDM Across Continental US for Testing Set",
       x = "Longitude Bin", 
       y = "Latitude Bin")

plot <- pred + actual
plot
# rmse and variance maps
#####


# trying an initial Xgboost model with data binned to 0.5 degree
#####
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_0.5")

allStatesDf_0.5 <- readRDS("CleanedUS_0.5.RDS")

# split into training and testing 
train_0.5 <- allStatesDf_0.5 %>% sample_frac(0.80)
test_0.5 <- anti_join(allStatesDf_0.5 , train_0.5)

train_x = data.matrix(train_0.5 %>% select(-c(USDM_Avg, x_Avg, y_Avg, Date)))
train_y = train_0.5 %>% select(USDM_Avg)
train_y <- train_y %>% select(-c(bin.x, bin.y))

#define predictor and response variables in testing set
test_x = data.matrix(test_0.5 %>% select(-c(USDM_Avg, x_Avg, y_Avg, Date)))
test_y = test_0.5 %>% select(USDM_Avg)

#define final training and testing sets
xgb_train = xgb.DMatrix(data = train_x, label = train_y$USDM_Avg)
xgb_test = xgb.DMatrix(data = test_x, label = test_y$USDM_Avg)

#defining a watchlist
watchlist = list(train=xgb_train, test=xgb_test)

#fit XGBoost model and display training and testing data at each iteartion
model = xgb.train(data = xgb_train, max.depth = 3, 
                  watchlist=watchlist, nrounds = 500, 
                  verbose = TRUE)



# different approach using caret
grid_tune <- expand_grid(nrounds = c(500, 1000, 1500), 
                         max_depth = c(2,4,6), 
                         eta = 0.3, 
                         gamma = 0, 
                         colsample_bytree = 1, 
                         min_child_weight = 1, 
                         subsample = 1)

train_ctrl <- trainControl(method = "cv", 
                           number = 3, 
                           verboseIter = TRUE, 
                           allowParallel = TRUE)

xgb_tune <-  train(x = train_x, 
                 y = train_y$USDM_Avg, 
                 trControl = train_ctrl, 
                 tuneGrid = grid_tune, 
                 method = "xgbTree", 
                 verbose = TRUE)


xbg_tune$best




#####

# create a plot of sd of USDM to see where the most variability in the US is
#####
library(ggplot2)
library(dplyr)

# Calculate the standard deviation of USDM_Avg for each location for all years
usdm_sd_by_location <- allStatesDf_0.5 %>%
  group_by(bin.x, bin.y) %>%
  summarize(
    USDM_SD = sd(USDM_Avg, na.rm = TRUE),
    x_Avg = first(x_Avg),
    y_Avg = first(y_Avg),
    .groups = "drop"
  )

# Create the spatial plot
ggplot(usdm_sd_by_location, aes(x = bin.x, y = bin.y, fill = USDM_SD)) +
  geom_tile() +
  scale_fill_viridis_c(name = "USDM SD") +
  theme_minimal() +
  labs(
    title = "Spatial Variability of US Drought Monitor Values",
    subtitle = "Standard Deviation of USDM by Location",
    x = "Longitude",
    y = "Latitude"
  ) +
  coord_fixed() +
  theme(legend.position = "right")
#####


# create ts plots of USDM and PDSI to show their differences 
#####

pdsi <- rast("agg_met_pdsi_1979_CurrentYear_CONUS.nc")
# tigris is not working for whatever reason
pima.pdsi <- crop.county.pdsi("Arizona", "Pima", pdsi, TRUE) 
  
# Get the direct download URL
# url <- tigris:::tigris_url("COUNTY", "2022")
# Or use this specific URL for 2022 county data
url <- "https://www2.census.gov/geo/tiger/TIGER2022/COUNTY/tl_2022_us_county.zip"

# Set local filenames
zip_file <- file.path(tempdir(), "counties.zip")
unzip_dir <- file.path(tempdir(), "counties")

# Download the file
download.file(url, zip_file, mode = "wb", method = "auto")

# Create directory and unzip manually
dir.create(unzip_dir, showWarnings = FALSE, recursive = TRUE)
unzip(zip_file, exdir = unzip_dir)

# Read the shapefile
counties_sf <- st_read(unzip_dir)

# Filter for Arizona counties
arizona_counties <- counties_sf[counties_sf$STATEFP == "04", ]

pima <- arizona_counties[arizona_counties$NAME == "Pima", ]
pdsi.data <- terra::extract(pdsi, pima,
                                 cells = TRUE, xy = TRUE)
# assign grid cell IDs
pdsi.data$ID <- seq_len(nrow(pdsi.data))

goodData <- select(pdsi.data, ID, cell, x, y,
                   starts_with("daily_mean_palmer_drought_severity_index_day="))

# get rid of of the long name
allNames <- names(goodData)
newNames <- str_remove(allNames,
                       "daily_mean_palmer_drought_severity_index_day=")
names(goodData) <- newNames

goodLong <- pivot_longer(goodData, cols = -c(ID, x, y, cell),
                         names_to = "Day", values_to = "PDSI") %>%
  # deal with the date
  mutate(
    day = as.integer(Day),
    Date = as.Date(day, origin = "1900-01-01"), # starting date
    midPointDay = day - 2.5)





#####


# trying again at 0.25- not useable, 0.25 dataset is too large
#####
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
#####