library(raster)
library(sf)
library(terra)
library(tigris)
library(dplyr)
library(tidyverse)
library(tidyterra)
library(ncdf4)
library(randomForest)
library(caret)
library(lubridate)

# This Function plots predicted PMDI data by year
plot.pmdi <- function(pmdi.predictions, predicted.col.name,
                      name.string, year, save = FALSE){
  # 
  
  # create plot for the predicted values
  pred.plot <-  ggplot(pmdi.predictions, aes(x = bin.x, y = bin.y, fill = .data[[predicted.col.name]])) +
    geom_tile() +
    scale_fill_manual(name = "USDM Category", 
                      values = c("None" = "white", 
                                 "D0" = "#E6B940", 
                                 "D1" = "#E67D2E", 
                                 "D2" = "#E5541B", 
                                 "D3" = "#C9281C", 
                                 "D4" = "#8E1C14")) +
    theme_minimal()  +
    labs(title = paste0("Actual Average USDM Across Continental US for ", year)) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "gray50"),
          legend.key.size = unit(0.8, "cm"),
          legend.key = element_rect(color = "gray50"))
  us_outline <- map_data("usa")
  
  pred.plot <- pred.plot + geom_path(data = us_outline, aes(x = long, y = lat, group = group), 
                               color = "black", linewidth = 0.7, inherit.aes = FALSE)
  pred.plot
  
  if(save){
    ggsave(pred.plot, file = paste0("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Plots/", name.string, ".png"))
  }
}

# This function creates a predicted vs. actual heatmap for categorical USDM predictions
pred.v.actual.plot.factor <- function(testing.set, actual.col.name, predicted.col.name, 
                                      save = FALSE, name.string, year){
  # First, create the basic plot without cell borders
  pred <- ggplot(testing.set, aes(x = bin.x, y = bin.y, fill = .data[[predicted.col.name]])) +
    geom_tile() +  # No borders on cells
    scale_fill_manual(name = "USDM Category", 
                      values = c("None" = "white", 
                                 "D0" = "#E6B940", 
                                 "D1" = "#E67D2E", 
                                 "D2" = "#E5541B", 
                                 "D3" = "#C9281C", 
                                 "D4" = "#8E1C14")) +
    theme_minimal() +
    labs(title = paste0("Predicted Average USDM Across Continental US for ", year)) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "gray50"),
          legend.key.size = unit(0.8, "cm"),
          legend.key = element_rect(color = "gray50"))  # Add borders just to legend keys
  
  # Now add a spatial outline for just the US
  us_outline <- map_data("usa")
  pred <- pred + geom_path(data = us_outline, aes(x = long, y = lat, group = group), 
                           color = "black", linewidth = 0.7, inherit.aes = FALSE)
  
  # create plot for the predicted values
  actual <-  ggplot(testing.set, aes(x = bin.x, y = bin.y, fill = .data[[actual.col.name]])) +
    geom_tile() +
    scale_fill_manual(name = "USDM Category", 
                      values = c("None" = "white", 
                                 "D0" = "#E6B940", 
                                 "D1" = "#E67D2E", 
                                 "D2" = "#E5541B", 
                                 "D3" = "#C9281C", 
                                 "D4" = "#8E1C14")) +
    theme_minimal()  +
    labs(title = paste0("Actual Average USDM Across Continental US for ", year)) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "gray50"),
          legend.key.size = unit(0.8, "cm"),
          legend.key = element_rect(color = "gray50"))
  
  actual <- actual + geom_path(data = us_outline, aes(x = long, y = lat, group = group), 
                               color = "black", linewidth = 0.7, inherit.aes = FALSE)
  plot <- pred + actual
  
  plot
  
  if(save){
    ggsave(plot, file = paste0(name.string, ".png"))
  }
}

library(tidyverse)

pred.v.actual.plot.factor <- function(testing.set, actual.col.name, predicted.col.name, 
                                      save = FALSE, name.string, year){
  
  # Reshape data to long format for faceting
  plot_data <- testing.set %>%
    mutate(
      Predicted = .data[[predicted.col.name]],
      Actual = .data[[actual.col.name]]
    ) %>%
    pivot_longer(cols = c(Predicted, Actual), 
                 names_to = "Type", 
                 values_to = "USDM_Category") %>%
    mutate(Type = factor(Type, levels = c("Predicted", "Actual")))
  
  # Get US outline
  us_outline <- map_data("usa")
  
  # Create faceted plot
  plot <- ggplot(plot_data, aes(x = bin.x, y = bin.y, fill = USDM_Category)) +
    geom_tile() +
    geom_path(data = us_outline, aes(x = long, y = lat, group = group), 
              color = "black", linewidth = 0.7, inherit.aes = FALSE) +
    scale_fill_manual(name = "USDM Category", 
                      values = c("None" = "white", 
                                 "D0" = "#E6B940", 
                                 "D1" = "#E67D2E", 
                                 "D2" = "#E5541B", 
                                 "D3" = "#C9281C", 
                                 "D4" = "#8E1C14")) +
    facet_wrap(~Type, ncol = 2) +
    theme_minimal() +
    labs(title = paste0("USDM Comparison Across Continental US for ", year)) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "gray50"),
          legend.key.size = unit(0.8, "cm"),
          legend.key = element_rect(color = "gray50"),
          strip.text = element_text(face = "bold", size = 11))
  
  if(save){
    ggsave(plot, filename = paste0("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Plots/", name.string, ".png"), width = 12, height = 5)
  }
  
  plot
  return(plot)
}

# This function creates a plot of variable importance based on the results from rf
plot.pdsi.importance <- function(rf.fit, train.set, save = FALSE, name.string){
  # grab the importance for each location from the rf fit object 
  importance <- data.frame(rf.fit$variable.importance.local)
  train_ungrouped <- train.set %>% ungroup()
  
  # Then add the importance columns
  train.importance <- train_ungrouped %>% 
    mutate(
      x.bin.importance = importance$bin.x, 
      y.bin.importance = importance$bin.y, 
      PDSI.importance = importance$PDSI_Avg
    )
  # bin training set by lat/long to look at importance in each grid cell
  train.binned <- train.importance %>% group_by(bin.x, bin.y) %>% 
    summarise(Mean.binx.importance = mean(x.bin.importance), 
              Mean.biny.importance = mean(y.bin.importance),
              Mean.PDSI.importance = mean(PDSI.importance), 
              .groups = 'drop')
  
  # plot the variable
  plot <- ggplot(train.binned, aes(x = bin.x, y = bin.y, fill = Mean.PDSI.importance)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Mean PDSI Importance") +
    theme_minimal() +
    labs(title = "PDSI Importance Across Continental US",
         x = "Longitude Bin", 
         y = "Latitude Bin")
  
  # save the plot if ya want 
  if(save){
    ggsave(plot, filename = paste0("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Plots/", name.string, ".png"), width = 10, height = 5)
  }
}

# This Function plots prediction error over the us for a continuous model 
plot.error <- function(test.set, error.metric,  save = FALSE, filename){
  # bin the testing set to get the MSE and MAE for each gridcell 
  test.binned <- test %>% 
    mutate(Abs_Error = abs(USDM_Avg - predicted), 
           Sq_Error = (USDM_Avg - predicted)^2) %>% group_by(bin.x, bin.y) %>% 
    summarise(MAE = mean(Abs_Error), 
              MSE = mean(Sq_Error),
              RMSE = sqrt(MSE), 
              .groups = 'drop')
  
  # select the metric to plot 
  
  
  # now actually plot the heatmap 
  ggplot(test.binned, aes(x = bin.x, y = bin.y, fill = MSE)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Mean Squared Error") +
    theme_minimal() +
    labs(title = "Prediction Error Across Continental US",
         x = "Longitude Bin", 
         y = "Latitude Bin")
}

# This function averages data from June through August to match with 
  # annual measures from the Living Blended Drought Atlas
summer.average <- function(cleaned.pdsi.data){
  # Create a column for the month/ year, select only the relevant months
  cleaned.pdsi.data <- cleaned.pdsi.data %>% 
    mutate(month = month(Date), year = year(Date)) %>% 
    filter(month == 6 | month == 7 | month == 8) %>% 
    # Average the PDSI for each year through the summer months 
    group_by(year, bin.x, bin.y) %>% 
    summarize(PDSI_Avg = mean(PDSI_Avg), 
              USDM_Avg = mean(USDM_Avg))
  return(cleaned.pdsi.data)
}

# This fucntion bins elevation data to the nearest specified degree of latitude and longitude
  # and created an average elevation for the binned area, in m
bin.elevation.data <- function(raw.elevation.data, grid.size = 0.1){
  
  # bins the lat/long into whatever decimal place you wanted 
  clean.data.binned <- raw.elevation.data %>% 
    mutate(bin.x = round.to.nearest(x,nearest = grid.size),
           bin.y = round.to.nearest(y,nearest = grid.size)) %>% 
    # groups by the new bins and averages data for the new larger, binned gridcells
    group_by(bin.x, bin.y) %>% 
    summarise(Elev_Avg = mean(Elevation_m, na.rm = TRUE), 
              x_Avg = mean(x, na.rm = TRUE), 
              y_Avg = mean(y, na.rm = TRUE))
  # retrn the cleaned, binned, data
  return(clean.data.binned)      
}


# This function reads in the continental US data, cleans it, 
  # and bins the data to the nearest degree given. 
  # It was intended to be used with the future function, to iterate over a df containing lists of county names
  # future function 
bin.state.data <- function(state.fips, LocationFactor = TRUE, nearest.degree){
  pdsi <- rast("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/agg_met_pdsi_1979_CurrentYear_CONUS.nc")
  
  # initialize data frame to hold cleaned state data
  state.data <- data.frame()
  
  # loop through all of the counties in the state
  for(index in 1:nrow(state.fips)){
    
    # grab the state, county, and fips code 
    county.name <- state.fips$AreaClean[index]
    state.name <- state.fips$State[index]
    fips <- state.fips$AOI.Value[index]
    
    # get the file name
    file.name <- paste0("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/countyData/USDM-", fips, ".csv")
    
    # read in the data
    drought.data <- read.csv(file.name)
    
    # processing
    clean.data <- clean.county.data(state.name, county.name, pdsi, drought.data, LocationFactor)
    
    # add state and county column for identification 
    clean.data <- clean.data %>% mutate(State = state.name,
                                        County = county.name)
    
    # bin data to the nearest .25 degree of lat/long
    binned.state <- bin.lat.long(clean.data, nearest.degree)
    
    # add all of the data together
    state.data <- rbind(binned.state, state.data)
  }
  saveRDS(state.data,file = paste(state.name,".RDS"))
  return(state.data)
}


# This is an internal function used in the bin.lat.long function
  # This function rounds the decimal given in the grid.size argument 
  # to make the number compatible with binning
round.to.nearest <- function(x,nearest){
  rounded <- round(x/nearest) * nearest
  return(rounded)
}

# This function takes a cleaned county dataset with columns of USDM_Avg, 
  # PDSI, x, y, and date. It bins the lat/long coordinates into the nearest 0.1
  # degree as default, but can be changes to whatever degree you want to round to. 
  # It them takes the average of PDSI, USDM, x, and y 
  # sig.digits is the number of decimal places you want to bin to
bin.lat.long <- function(clean.data, grid.size = 0.1){
  # bins the lat/long into whatever decimal place you wanted 
  clean.data.binned <- clean.data %>% 
    mutate(bin.x = round.to.nearest(x,nearest = grid.size),
           bin.y = round.to.nearest(y,nearest = grid.size)) %>% 
    # groups by the new bins and averages data for the new larger, binned gridcells
    group_by(bin.x, bin.y, Date) %>% 
    summarise(PDSI_Avg = mean(PDSI, na.rm = TRUE), 
              x_Avg = mean(x, na.rm = TRUE), 
              y_Avg = mean(y, na.rm = TRUE), 
              USDM_Avg = mean(USDM_Avg, na.rm = TRUE))
  # retrn the cleaned, binned, data
  return(clean.data.binned)      
}

# This function crops the pdsi file to the county level
  # Takes the state and county desired as strings, a spatraster object containing pdsi data
  # I added a T/F variable: Location Factor, to tell if we want to include the x, y, and cell numbers from the PDSI gridcells
  # This addition helps with plotting results later 
  ## note this raster file has to be loaded in using the rast function from the terra package
crop.county.pdsi <- function(state, county, spat.us.pdsi, LocationFactor){ # county and state need to be saved as strings 
  # # get the specific county
  # sf.county <- counties(state = state,
  #                       class = "sf")
  # # first I'll start with just selecting  county
  # county <- sf.county[sf.county$NAME == county, ]
  # 
  # # tigris is broken, so we have to read in the whole country and filter by county

  sf.county <- st_read("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/tl_2024_us_county/tl_2024_us_county.shp")
  
  county <- sf.county[sf.county$NAME == county, ]
  
  # if we want the x/y coordinates and cell numbers from PDSI measurements, w
  if(LocationFactor){
    pdsi.extracted <- terra::extract(spat.us.pdsi, county,
                                     cells = TRUE, xy = TRUE) # creates a data frame
  }
  else{
    # use extract function, cells and xy argument return the cell number and xy coordinates, respectively. 
    # added to help with mapping visualizatona later on 
    pdsi.extracted <- terra::extract(spat.us.pdsi, county) # creates a data frame
  }
  return(pdsi.extracted)
}

# This function crops the pdsi file to the state level
  # Takes the state and county desired as strings, a spatraster object containing pdsi data
  # I added a T/F variable: Location Factor, to tell if we want to include the x, y, and cell numbers from the PDSI gridcells
  # This addition helps with plotting results later 
  ## note this raster file has to be loaded in using the rast function from the terra package
crop.state.pmdi <- function(state, spat.us.pmdi, LocationFactor){ # county and state need to be saved as strings 

  # use tigris's states function to get shape files for each state 
  sf.country <- states()
  
  # filter to the specific state you want the area of 
  state.shapefile <- sf.country[sf.country$STUSPS == state, ]
  
  # if we want the x/y coordinates and cell numbers from PDSI measurements, w
  if(LocationFactor){
    pmdi.extracted <- terra::extract(spat.us.pmdi, state.shapefile,
                                     cells = TRUE, xy = TRUE) # creates a data frame
  }
  else{
    # use extract function, cells and xy argument return the cell number and xy coordinates, respectively. 
    # added to help with mapping visualizatona later on 
    pmdi.extracted <- terra::extract(spat.us.pmdi, state.shapefile) # creates a data frame
  }
  return(pmdi.extracted)
}

# This function normalizes the land area percentages for the US drought index
  # note this should be done before calculating the weighted averages
normalize.usdm <- function(data){
  # D4 stays the same
  # D3 gets D4 subtracted from it 
  data$normal.D3 <- data$D3 - data$D4
  # D2 has D3 subtracted from it 
  data$normal.D2 <- data$D2 - data$D3
  # and so on 
  data$normal.D1 <- data$D1 - data$D2
  data$normal.D0 <- data$D0 - data$D1
  # save values from temporary columns to the original to keep names consistent 
  data$D3 <- data$normal.D3
  data$D2 <- data$normal.D2
  data$D1 <- data$normal.D1
  data$D0 <- data$normal.D0

  data$total <- rowSums(data %>% select(c(D0, D1, D2, D3, D4, None)))
  # return the data frame minus the temporary columns 
  data <- subset(data, select = 
                   -c(normal.D3, normal.D2, normal.D1, normal.D0))
  return(data)
}

# This function calculates and creates a new column for the weighted average of the US drought index
  # the average is weighted by the land area percent for each index 
usdm.weighted.average <- function(data) {
  # Define the weights
  weights <- 0:5
  
  # Calculate the weighted average across columns D0 to D4
  data$USDM_Avg <- mapply(function(None, D0, D1, D2, D3, D4) {
    values <- c(None, D0, D1, D2, D3, D4)
    weighted_sum <- sum(values * weights)
    weight_sum <- sum(values)
    
    if (weight_sum == 0) {
      return(NA)  # Avoid division by zero error
    } else {
      return(weighted_sum / weight_sum)
    }
  }, data$None, data$D0, data$D1, data$D2, data$D3, data$D4)
  
  return(data)  # Return the modified data frame
}

# This function interpolates values for PDSI to find the midpoint of the 5 day cycle which the data is collected on. 
  # the midpoint date is the same for the interpolate.usdm function
interpolate.pdsi <- function(pdsi.data, LocationFactor){
  
  # assign grid cell IDs
  pdsi.data$ID <- seq_len(nrow(pdsi.data))
  
  # create test to see if the location data is included, if to processing will be a bit different
  if(LocationFactor){
    # select only the layers which have pdsi measurements, including location data
    goodData <- select(pdsi.data, ID, cell, x, y,
                       starts_with("daily_mean_palmer_drought_severity_index_day="))
  }
  else{
    # select only the layers which have pdsi measurements
    goodData <- select(pdsi.data, ID, 
                       starts_with("daily_mean_palmer_drought_severity_index_day="))
  }
  
  # get rid of of the long name
  allNames <- names(goodData)
  newNames <- str_remove(allNames,
                         "daily_mean_palmer_drought_severity_index_day=")
  names(goodData) <- newNames
  
  # check again for location data
  if(LocationFactor){
    # pivot so we can have all of the dates in one column
    goodLong <- pivot_longer(goodData, cols = -c(ID, x, y, cell),
                             names_to = "Day", values_to = "PDSI") %>%
      # deal with the date
      mutate(
        day = as.integer(Day),
        Date = as.Date(day, origin = "1900-01-01"), # starting date
        midPointDay = day - 2.5) # subtract 2.5 to find the midpoint of the 5 day cycle
  }
  else{
    # pivot so we can have all of the dates in one column 
    goodLong <- pivot_longer(goodData,cols = -ID, names_to = "Day", values_to = "PDSI")
    
    
    # fix the dating system from day numbers to actual dates
    baselineDate <- as.Date("1900-01-01") # still didn't find the start date from NOAA..
    goodLong$date <- as.numeric(goodLong$Day) + baselineDate
    
    # now we begin aligning the grid cells and finding the mid point days 
    goodLong$midPointDay <- as.numeric(goodLong$Day) - 2.5 # subtract 2.5 to find the midpoint of the 5 day cycle 
  }
  return(goodLong)
}

# This function interpolates values for the us drought monitor 
  # the interpolated dates are the same for the PDSI data 
interpolate.usdm <- function(drought.data){
  # deal with date
  baselineDate <- as.Date("1900-01-01") # still didn't find the start date from NOAA..
  
  # now we need to convert the drought monitor into days since 1900, and calculate midPointDay
  drought.data$ValidMid <- as.Date(drought.data$ValidEnd) - 3.5 # subtract 3.5 to find the midpoint of the 7 day cycle
  
  # actually apply the dates to our drought data  
  drought.data$midPointDay <- as.numeric(drought.data$ValidMid - baselineDate)
  
  return(drought.data)
}

# This function joins the interpolated USDM and PDSI data sets, and puts the data in a format ready for modeling 
join.pdsi.usdm <- function(pdsi.data, usdm.data, LocationFactor){
  # we simply have to join by date
  drought.pdsi <- left_join(pdsi.data, usdm.data)
  
  # leaving out any na values 
  drought.pdsi <- na.omit(drought.pdsi)
  
  # create another location factor check 
  if(LocationFactor){
    # select only the columns we want
    reshaped_data <- drought.pdsi %>%
      select(Date, x, y, cell, USDM_Avg, PDSI)
  }
  else{
    # reshape the data so each row contains a single day's worth of PSDI measurements, containing one column for each gridcell's pdsi values
    reshaped_data <- drought.pdsi %>%
      select(date, ID, PDSI, USDM_Avg) %>%         # Select relevant columns
      pivot_wider(
        names_from = ID,                # Pivot based on the grid cell ID
        values_from = PDSI,             # Use PDSI as values for the new columns
        names_prefix = "PDSI_"          # Add a prefix to column names
      )
  }
  return(reshaped_data)
}

# This function combines all of the above cleaning for a streamlined approach 
  # the arguments needed are the state and county names as strings, 
  # the pdsi dataset from drought.gov, loaded in as a spatraster object, 
  # and the usdm data from drought.gov for the specific county
clean.county.data <- function(state, county, spat.us.pdsi, 
                              cropped.usdm.data, LocationFactor){
  # crop the pdsi data 
  pdsi <- crop.county.pdsi(state, county, spat.us.pdsi, LocationFactor)
  
  # clean the usdm data
  normal.drought <- normalize.usdm(cropped.usdm.data)
  weighted.drought <- usdm.weighted.average(normal.drought)
  interpolated.drought <- interpolate.usdm(weighted.drought)
  
  # clean the pdsi data 
  interpolated.pdsi <- interpolate.pdsi(pdsi, LocationFactor)
  
  #Put the usdm data on the same timescale as the pdsi data
  allPDSIDays <- unique(interpolated.pdsi$midPointDay)
  
  interpolated.usdm <- approx(x = interpolated.drought$midPointDay, 
                              y = interpolated.drought$USDM_Avg,
                              xout = allPDSIDays) %>% 
    as.data.frame() %>% 
    filter(!is.na(y))
  
  names(interpolated.usdm) <- c("midPointDay","USDM_Avg")
  
  
  # join the datasets and reshaping 
  full.data <- join.pdsi.usdm(interpolated.pdsi, interpolated.usdm, LocationFactor)
  
  # add an ID column 
  # full.data$id <- 1:nrow(full.data)
  
  # return the final cleaned dataset
  return(full.data)
}

clean.pmdi.data <- function(state, spat.us.pmdi, LocationFactor = TRUE){
  # extract pmsi values for each state using terra
  pmdi <- crop.state.pmdi(state, spat.us.pmdi, LocationFactor)
  
  # assign grid cell IDs
  pmdi$ID <- seq_len(nrow(pmdi))
  
  # select only the layers with PMDI data, including the location 
  goodData <- select(pmdi, ID, cell, x, y, 
                     starts_with("PMDI_"))
  
  # get rid of of the long name
  allNames <- names(goodData)
  newNames <- str_remove(allNames,
                         "PMDI_")
  names(goodData) <- newNames
  
  str(goodData)
  
  # pivot so we can have all of the dates in one column
  goodLong <- pivot_longer(goodData, cols = -c(ID, x, y, cell),
                           names_to = "Year", values_to = "PMDI")
  
  return(goodLong)
  
}

# This function converts the weighted average USDM rating back into categorical variable 
convert.cat.USDM <- function(usdm.data){
  usdm.data <- usdm.data %>% 
    mutate(USDM_factor = case_when(
            USDM_Avg < 0.5 ~ "None",
            USDM_Avg < 1.5 ~ "D0",
            USDM_Avg < 2.5 ~ "D1",
            USDM_Avg < 3.5 ~ "D2",
            USDM_Avg < 4.5 ~ "D3",
            TRUE ~ "D4")) %>% 
    mutate(USDM_factor = factor(USDM_factor))
  return(usdm.data)
}


# This function creates a random forest model, splitting the data into a training set, 
  # tests the model, and returns A dataset with the weighted averages used to test the rf model, 
  # and the predictions made
quick.rf <- function(pdsi.usdm.data){
  
  set.seed(1)
  
  # split into training and testing 
  train <- pdsi.usdm.data %>% sample_frac(0.80)
  test <- anti_join(pdsi.usdm.data, train, by = 'id')
  
  # build model
  rf.fit <- randomForest(USDM_Avg ~ ., 
                         data = train, 
                         importance = TRUE)
  
  # predict and find RMSE
  preds <- predict(rf.fit, test)
  
  results <- data.frame(test = test$USDM_Avg, 
                        predictions = preds)
  
  return(results)
}


