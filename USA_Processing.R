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
# loop through continental US states
for(index in seq_along(states)){
  
  # filter for just the state's data
  state.fips[[index]] <- contUS %>% filter(State == states[index])
  
}

# saveRDS(state.fips, file = "CleanStateFips.RDS")


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

# filter for west-coast data only 
westUS <- allStatesDf %>% filter(bin.x <= -100)

subsetWest <- westUS %>% sample_frac(0.50)

# rf modeling 
# split into training and testing 
train <- subsetWest %>% sample_frac(0.80)
test <- anti_join(subsetWest , train)

# using caret
library(caret)

ctrl <- trainControl(method = "cv", number = 5, verboseIter = TRUE)

rf.caret.fit <- train(USDM_Avg ~ PDSI_Avg + bin.x + bin.y, 
                data = train,
                method = "rf", 
                ntree = 50,
                trControl = ctrl)
# using ranger 
library(ranger)

rf.ranger.fit <- ranger(
  formula = USDM_Avg ~ PDSI_Avg + bin.x + bin.y,  # Formula specifying the target and predictors
  data = train,  # Training dataset
  num.trees = 100,  # Number of trees in the forest
  mtry = 2,  # Number of features to consider at each split
  classification = TRUE,  # Specify that it's a classification problem
  importance = 'impurity',  # To measure feature importance
  probability = TRUE,  # To get probabilities for each class
  verbose = TRUE
  )

# predict and find RMSE
preds <- predict(rf.fit, test)


RMSE(preds, test$USDM_Avg)


ggplot(allStatesDf, aes(x=Date, y=PDSI_Avg))+
  geom_point()

