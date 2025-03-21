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



# future function 
bin.state.data <- function(state.fips, LocationFactor = TRUE){
  pdsi <- rast("agg_met_pdsi_1979_CurrentYear_CONUS.nc")
  
  # initialize data frame to hold cleaned state data
  state.data <- data.frame()
  
  # loop through all of the counties in the state
  for(index in 1:nrow(state.fips)){
    
    # grab the state, county, and fips code 
    county.name <- state.fips$AreaClean[index]
    state.name <- state.fips$State[index]
    fips <- state.fips$AOI.Value[index]
    
    # get the file name
    file.name <- paste0("countyData/USDM-", fips, ".csv")
    
    # read in the data
    drought.data <- read.csv(file.name)
    
    # processing
    clean.data <- clean.county.data(state.name, county.name, pdsi, drought.data, LocationFactor)
    
    # add state and county column for identification 
    clean.data <- clean.data %>% mutate(State = state.name,
                                        County = county.name)
    
    # bin data to the nearest .25 degree of lat/long
    binned.state <- bin.lat.long(clean.data, 0.25)
    
    
    # add all of the data together
    state.data <- rbind(binned.state, state.data)
  }
  saveRDS(state.data,file = paste(state.name,".RDS"))
  return(state.data)
}



#####################################
# new filtering loop 

states <- unique(contUS$State)
state.fips <- vector(mode = "list",length = length(states))
# loop through continental US states
for(index in seq_along(states)){
  
  # filter for just the state's data
  state.fips[[index]] <- contUS %>% filter(State == states[index])
  
  # # bin each county in the state and combine to df using our function 
  # binned.state <- bin.state.data(state.fips, pdsi, true)
  # 
  # # save state ?
  # savefilename <- paste0(state.name, "_Binned.csv")
  # write_csv(binned.state, file = savefilename)
  
  # add to US data frame?
  
}



# Run in furrr ------------------------------------------------------------
library(future)
library(furrr)
future::plan(strategy = multisession, workers = 6)

allStates <- furrr::future_map(state.fips,bin.state.data,pdsi,.progress = TRUE)









