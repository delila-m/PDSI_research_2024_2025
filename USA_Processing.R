# load in necessary libraries
source("C:/Users/delil/Desktop/Fall 2024/Research 2024/PDSI_research_2024/drought_functions.R")

# load in PDSI data for whole country
pdsi <- rast("agg_met_pdsi_1979_CurrentYear_CONUS.nc")

# load in FIPS codes and filter for only continuous US
usa_fips <- read.csv("All Counties and FIPS Codes.csv")
contUs <- filter(countryfips, ! State %in% c("AK","HI","PR"))

# initialize data frame 
usa.data <- data.frame()

# loop through each county 
for (index in seq_along(contUs)){
  county.name <- contUs$Area[index]
  state.name <- contUS$State[index]
  fips <- az.fips$AOI.Value[index]
  
  print(county.name)
  print(state.name)
  print(fips)
  print("\n\n")
  # get the file name 
  file.name <- paste0("C:\Users\dgm239\Downloads\Research_2025\PDSI_research_2024\countyData/USDM-", fips, ".csv")
  
  # read in the data
  drought.data <- read.csv(file.name)
  
  # processing 
  clean.data <- clean.county.data(state.name, county.name, pdsi, drought.data, TRUE)  
  
  # add all of the data together
  usa.data <- rbind(usa.data, clean.data)
}

# bin data to the nearest degree of latitude 
