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

################ # plot some stuff

# load in the rf objects
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")
load("RFAnalysis1.Rdata")

library(viridis)  
library(ggplot2)
library(metR)

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

###### Trying the model again with elevation as a predictor
# loading in data
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/CleanedUS_1")
allStatesDf <- readRDS("CleanedUS_1.RDS")

# split into training and testing 
train <- allStatesDf %>% sample_frac(0.80)
test <- anti_join(allStatesDf , train)

# Install and load the necessary packages
install.packages(c("elevatr", "sf", "raster"))
library(elevatr)
library(sf)
library(raster)

# Define the continental US bounding box
# Format: c(xmin, ymin, xmax, ymax)
continental_us_bbox <- data.frame(
  lon = c(-124.848974, -124.848974, -66.885444, -66.885444, -124.848974),
  lat = c(24.396308, 49.384358, 49.384358, 24.396308, 24.396308))

# Get elevation data using get_elev_raster
# zoom parameter controls resolution (higher = more detailed but larger file)
elevation_data <- get_elev_raster(locations = continental_us_bbox, 
                                  prj = "EPSG:4326", 
                                  z = 7,  # Adjust zoom level as needed
                                  src = "aws")  # Source: AWS Terrain Tiles

# Plot the elevation data
plot(elevation_data)

############################# wrong approach...
######
# 1. Define the geographic extent based on your bin ranges
x_range <- range(train$bin.x)
y_range <- range(train$bin.y)

# 2. Create a grid covering your entire study area
# Expand slightly to ensure coverage
extent <- c(
  x_range[1] - 0.5,
  x_range[2] + 0.5,
  y_range[1] - 0.5,
  y_range[2] + 0.5
)

# Create a polygon representing the study area
study_area <- st_polygon(list(rbind(
  c(extent[1], extent[3]),
  c(extent[2], extent[3]),
  c(extent[2], extent[4]),
  c(extent[1], extent[4]),
  c(extent[1], extent[3])
)))
study_area_sf <- st_sf(geometry = st_sfc(study_area, crs = 4326))

# 3. Get elevation data for the entire area
  # z = 8: ~ 600m per pixel 
  # z = 9: ~300 meters per pixel (good for continental analysis)
  # z = 10: ~150 meters per pixel
  # z = 11: ~75 meters per pixel
  # z = 12: ~38 meters per pixel
  # z = 13: ~19 meters per pixel
  # z = 14: ~9.5 meters per pixel (very detailed)
elev_raster <- get_elev_raster(study_area_sf, z = 9, src = "aws",
                               verbose = TRUE)
# crop bc degree- buffer is problematic 
# Define the original extent without the 0.5-degree buffer
original_extent <- c(
  x_range[1],       # Original min longitude
  x_range[2],       # Original max longitude
  y_range[1],       # Original min latitude
  y_range[2]        # Original max latitude
)

# Create a polygon representing just the original study area
original_area <- st_polygon(list(rbind(
  c(original_extent[1], original_extent[3]),
  c(original_extent[2], original_extent[3]),
  c(original_extent[2], original_extent[4]),
  c(original_extent[1], original_extent[4]),
  c(original_extent[1], original_extent[3])
)))
original_area_sf <- st_sf(geometry = st_sfc(original_area, crs = 4326))

# Convert sf object to Spatial object for use with raster
original_area_sp <- as(original_area_sf, "Spatial")

# Crop the elevation raster to the original extent
elev_raster_cropped <- crop(elev_raster, original_area_sp)

# 4. Aggregate the raster to 1-degree resolution to match your bins
# Create a target raster with 1-degree cells
target_raster <- raster(
  xmn = floor(original_extent[1]), 
  xmx = ceiling(original_extent[2]),
  ymn = floor(original_extent[3]), 
  ymx = ceiling(original_extent[4]),
  resolution = 1
)

# Aggregate the high-resolution elevation to the target
agg_elev <- aggregate(elev_raster_cropped, fact = c(
  res(target_raster)[1] / res(elev_raster_cropped)[1],  
  res(target_raster)[2] / res(elev_raster_cropped)[2]
), fun = mean)

# 5. Extract values for each bin
unique_bins <- train %>%
  select(bin.x, bin.y) %>%
  distinct()

# Create points for extraction
bin_points <- unique_bins %>%
  st_as_sf(coords = c("bin.x", "bin.y"), crs = 4326)

# Extract values
bin_elevations <- raster::extract(agg_elev, bin_points)
bin_elevations <- bin_elevations %>% ungroup()

unique_bins <- unique_bins %>% ungroup()

# Create final dataframe
elevation_df <- unique_bins %>%
  mutate(elevation = bin_elevations)

# 6. Join to your training data
train_with_elev <- train %>%
  left_join(bin_elevations, by = c("bin.x", "bin.y"))




# Check the range of elevation values
summary(bin_elevations$elevation)

# Look at extremely low values
low_elevations <- train_with_elev %>% 
  filter(elevation < -100) %>%
  arrange(elevation)

# View the first few extremely low points
head(low_elevations)

# Plot the points on a map to see where they are
library(ggplot2)
ggplot(train_with_elev, aes(x = bin.x, y = bin.y, color = elevation)) +
  geom_point() +
  scale_color_gradientn(
    colors = c("red", "orange", "yellow", "green", "blue"),
    values = scales::rescale(c(min(train_with_elev$elevation), 
                               -100, 0, 1000, max(train_with_elev$elevation)))
  ) +
  theme_minimal() +
  labs(title = "Elevation by Location", 
       subtitle = "Note: Extremely negative values shown in red")


# Check the resolutions of both rasters
res_elev <- res(elev_raster)
res_target <- res(target_raster)

# Print them to see what's happening
print(paste("Elevation raster resolution:", res_elev[1], "x", res_elev[2]))
print(paste("Target raster resolution:", res_target[1], "x", res_target[2]))

# Calculate aggregation factors
x_fact <- res_target[1] / res_elev[1]  # Note: DIVISION REVERSED from original code
y_fact <- res_target[2] / res_elev[2]

print(paste("Aggregation factors:", x_fact, "x", y_fact))
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
