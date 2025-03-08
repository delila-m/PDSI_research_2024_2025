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
      TRUE ~ Area
    )
  )

# initialize data frame 
usa.data <- data.frame()

# loop through each county 
for (index in 1:nrow(contUS)){
  county.name <- contUS$AreaClean[index]
  state.name <- contUS$State[index]
  fips <- contUS$AOI.Value[index]

  # get the file name
  file.name <- paste0("countyData/USDM-", fips, ".csv")

  # read in the data
  drought.data <- read.csv(file.name)

  # processing
  clean.data <- clean.county.data(state.name, county.name, pdsi, drought.data, TRUE)

  # add state and county column for identification 
  clean.data <- clean.data %>% mutate(State = state.name,
                                      County = county.name)
  
  # add all of the data together
  usa.data <- rbind(usa.data, clean.data)
}

# save full dataset


# bin data to the nearest degree of latitude 
