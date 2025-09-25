library(httr2)
library(purrr)
library(glue)
library(dplyr)


#get county data by fips
getUsdmByCounty <- function(fips, startDate = "01/01/2000", endDate = "12/31/2024", output_file = NA) {

  fips <- as.character(fips)

  if(nchar(fips) == 4){
    fips <- paste0("0",fips)
  }


  if(nchar(fips) != 5){
    stop("The fips code must be 5 digits. Add a zero at the beginning if need be")
  }

  # Construct the URL with proper encoding
  base_url <- "https://usdmdataservices.unl.edu/api/CountyStatistics/GetDroughtSeverityStatisticsByAreaPercent"

  # Build the query parameters
  params <- list(
    aoi = fips,
    startdate = startDate,
    enddate = endDate,
    statisticsType = 1
  )

  # Create the request object with encoded URL
  req <- httr2::request(base_url) |>
    httr2::req_url_query(!!!params)

  # Perform the request
  resp <- req |> httr2::req_perform()

  if(all(is.na(output_file))){
    output_file <- glue::glue("countyData/USDM-{fips}.csv")
  }

  response_body <- httr2::resp_body_string(resp)

  # Save the raw CSV content to a file
  csv_content <- httr2::resp_body_raw(resp)
  writeBin(csv_content, con = output_file)

  message("CSV file saved to ", output_file)
  return(output_file)
}

allCountyCodes <- googlesheets4::read_sheet(ss = "1G5ZT98aSEZuZ-YFuNv8PetSfkLXlYwYOKKk9w-dvD6A")


#only look at contiguous US
contUs <- filter(allCountyCodes, ! State %in% c("AK","HI","PR"))

#download all counties
walk(contUs$`AOI Value`,getUsdmByCounty)





