library(dplyr)
library(tidyverse)
# load in drought functions
#setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Functions/") #glg
#setwd("C:/Users/delil/Desktop/NAU/Research 2024-2025/PDSI_research_2024_2025/Functions/") #laptop
source("Functions/drought_functions.R")

# file path on laptop
# load("C:/Users/delil/Desktop/NAU/Research 2024-2025/PDSI_research_2024_2025/Data/PMDIPredictionSet.Rdata")
# # file path in glg
# load("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/PMDIPredictionSet.Rdata")
load("Data/PMDIPredictionSet.Rdata")

# weekly data
load("Data/SampledWeekly_0.5.RData")

# Previous categorical model for comparison
#setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
load("Data/UpdatedCleaned_0.5/RFAnalysis0.5YearSplit.Rdata")

test_0.5 <- convert.cat.USDM(test.yearsplit.0.5, "predictions")
test_0.5_actual <- convert.cat.USDM(test.yearsplit.0.5, "USDM_Avg")
train_0.5 <- convert.cat.USDM(train.yearsplit.0.5, "USDM_Avg")

train_preds <- data.frame(preds = rf.ranger.fit.yearsplit.0.5$predictions)
train_preds_cat <- convert.cat.USDM(train_preds, "preds")
train_preds_cat$USDM_factor <- factor(
  train_preds_cat$USDM_factor,
  levels = levels(train_0.5$USDM_factor)
)

full_modern_set_0.5 <- bind_rows(train_0.5, test_0.5)

# table evaluating class metrics 
library(yardstick)
library(gt)
bal_acc_oob <- tibble(
  truth = train_0.5$USDM_factor,
  estimate = train_preds_cat$USDM_factor) %>%
  bal_accuracy(truth, estimate) %>%
  pull(.estimate)

bal_acc_test <- tibble(
  truth = test_0.5$USDM_factor,
  estimate = test_0.5_actual$USDM_factor) %>%
  bal_accuracy(truth, estimate) %>%
  pull(.estimate)

bal_acc_tbl <- tibble(
  Dataset = c("Training (2000-2020)", "Test (2020-2024)"),
  `Balanced Accuracy` = c(bal_acc_oob, bal_acc_test))

library(gt)

bal_acc_tbl %>%
  gt() %>%
  tab_header(
    title = "Random Forest Balanced Accuracy") %>%
  fmt_percent(
    columns = `Balanced Accuracy`,
    decimals = 2) %>%
  cols_align(
    align = "center",
    columns = `Balanced Accuracy`) %>%
  opt_table_outline() %>%
  tab_source_note(
    source_note = "Training performance estimated via OOB predictions" )



pmdi_one_cell <- pmdi_prediction_set %>%
  filter(bin.x == -113.75 & bin.y == 35.25) %>%
  mutate(intensity = case_when(predictions == "None" ~ 1,
                                    predictions == "D0" ~ 2,
                                    predictions == "D1" ~ 3,
                                    predictions == "D2" ~ 4,
                                    predictions == "D3" ~ 5,
                                    predictions == "D4" ~ 6,
                                    TRUE ~ -1 ))

# claude generated loop to calculate drought duration
# Define drought intensity mapping
drought_map <- c('None' = 1, 'D0' = 2, 'D1' = 3, 'D2' = 4, 'D3' = 5, 'D4' = 6)


# Identify drought events (any value >= D0)
pmdi_one_cell$is_drought <- pmdi_one_cell$intensity >= 2

# Create drought event IDs
pmdi_one_cell$drought_id <- NA
drought_event_id <- 0

for (i in 1:nrow(pmdi_one_cell)) {
  if (pmdi_one_cell$is_drought[i]) {
    if (i == 1 || !pmdi_one_cell$is_drought[i-1]) {
      # Start of new drought event
      drought_event_id <- drought_event_id + 1
    }
    pmdi_one_cell$drought_id[i] <- drought_event_id
  }
}

# Calculate drought event statistics
drought_events <- pmdi_one_cell %>%
  filter(!is.na(drought_id)) %>%
  group_by(drought_id) %>%
  summarise(
    start_year = first(year),
    end_year = last(year),
    duration = n(),
    avg_intensity = mean(intensity),
    max_intensity = max(intensity),
    min_intensity = min(intensity),
    avg_PDSI = mean(PDSI_Avg),
    drought_sequence = paste(predictions, collapse = ", "),
    .groups = 'drop'
  )
#Add max intensity label
intensity_reverse_map <- names(drought_map)
names(intensity_reverse_map) <- drought_map

drought_events$max_intensity_label <- intensity_reverse_map[as.character(drought_events$max_intensity)]
drought_events$avg_intensity_label <- intensity_reverse_map[as.character(round(drought_events$avg_intensity))]

# 2. Return interval by drought intensity
intensity_levels <- c('D0', 'D1', 'D2', 'D3', 'D4')
return_intervals <- data.frame(
  intensity = character(),
  count = integer(),
  return_interval = numeric(),
  stringsAsFactors = FALSE
)

for (level in intensity_levels) {
  count <- sum(drought_events$max_intensity_label == level, na.rm = TRUE)
  if (count > 0) {
    ri <- 2018 / count
    return_intervals <- rbind(return_intervals,
                              data.frame(intensity = level,
                                         count = count,
                                         return_interval = ri))
  }
}
severe_droughts <- drought_events %>%
  filter(max_intensity >= 3)  # D2 and above

# Calculate return intervals for different duration thresholds
duration_values <- sort(unique(severe_droughts$duration))

# Create data for two lines:
# Line 1: Number of drought events by duration
# Line 2: Return interval by duration

plot_data <- data.frame()

for (dur in duration_values) {
  count <- sum(severe_droughts$duration >= dur)
  ri <- 2018 / count
  plot_data <- rbind(plot_data,
                     data.frame(duration = dur,
                                count = count,
                                return_interval = ri))
}

plota <- ggplot(plot_data, aes(x = duration)) +
  geom_line(aes(y = count, color = "Number of Events"), size = 1.2) +
  geom_point(aes(y = count, color = "Number of Events"), size = 3) +
  geom_line(aes(y = return_interval, color = "Return Interval"), size = 1.2) +
  geom_point(aes(y = return_interval, color = "Return Interval"), size = 3) +
  scale_color_manual(values = c("Number of Events" = "blue",
                                "Return Interval" = "red")) +
  labs(
    title = "Drought Duration vs Return Interval",
    subtitle = "For drought events >= D2",
    x = "Drought Duration (years)",
    y = "Years",
    color = "Metric"
  ) +
  #scale_x_log10()+
  scale_y_log10()+
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

plota

# test functionality of above operations
  # test all separated out
cropped_test <- crop.cell(pmdi_prediction_set, xbin = -124.25, ybin = 40.25, pred_col = "predictions")
drought_events_test <- identify.drought(cropped_test, intensity_threshold = 4)

list_test <- summarize.drought.events(drought_events_test, time_col = "year", pred_col = "predictions")
                                     # duration_unit = "years")
drought_events_test2 <- list_test[[1]]
drought_intervals_test <- list_test[[2]]

plot_data <- plot.data(drought_events_test2)
slope_test <- recurrence.slope(plot_data)


  # test inside the evaluate.recurrence function
      # cell (-91.25, 33.25) only has dourght events with the duration of 1 year
full_test <- evaluate.annual.recurrence(pmdi_prediction_set, xbin = -113.75, ybin = 35.25,
                                 intensity_threshold = 4,
                                 pred_col = "predictions", time_col = "year")


full_drought_events_test2 <- full_test[[2]]
full_drought_intervals_test <- full_test[[3]]

plot_data <- plot.data(full_drought_events_test2)
slope_test <- recurrence.slope(plot_data)

# test the plotting
test_plot <- plot.annual.duration.v.return(plot_data)
test_plot


# plot weekly/annual together
weekly_df <- plot_data_weekly %>%
  mutate(scale = "Weekly")

annual_df <- plot_data %>%
  mutate(scale = "Annual")

combined_df <- bind_rows(weekly_df, annual_df)

ggplot(combined_df, aes(x = duration, y = return_interval, color = scale)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Weekly" = "#2E55A3",
      "Annual" = "#A34E2E"
    )
  ) +
  scale_y_log10() +
  labs(
    title = "Drought Duration vs Return Interval",
    subtitle = "Weekly and Annual Aggregations (Drought ≥ D2)",
    x = "Drought Duration (Years)",
    y = "Return Interval (Years)",
    color = "Aggregation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave("Plots/modern_and_paleo_duration_v_return_D2.png", width = 10, height = 6)


#create a data.frame of unique X/Y pairs to loop through
latlongpairs <- full_modern_set_0.5 %>%
  distinct(bin.x, bin.y)

# initialize data frame
us_slopes <- data.frame()

# add a progress bar
pb <- txtProgressBar(min = 0, max = nrow(latlongpairs), style = 3)


# now loop through the whole country
for(index in 1:nrow(latlongpairs)){

  # grab the cell coordinates
  current_xbin = latlongpairs$bin.x[index]
  current_ybin = latlongpairs$bin.y[index]

  # evaluate the recurrence intervals for that cell
  recurrence_list <- evaluate.annual.recurrence(full_modern_set_0.5, xbin = current_xbin,
                                         ybin = current_ybin,
                                         intensity_threshold = 3,
                                         pred_col = "USDM_factor", time_col = "Date")

  drought_events <- recurrence_list[[2]]

  # Check if valid drought events exist (handles NA, NULL, and empty)
  has_droughts <- !is.null(drought_events) &&
    is.data.frame(drought_events) &&
    nrow(drought_events) > 0

  # make sure there are drought events over our threshold in that cell
  if(has_droughts){

    # grab the plotting data and find the slope
    intervals <- plot.data(drought_events)
    #print(intervals)
    # make sure that there are drought events with at least two different durations
    if(nrow(intervals) > 0 && max(intervals$duration) >= 2){
      slope <- recurrence.slope(intervals)
     # print(slope)
      # add it all to the existing data frame
      us_slopes <- bind_rows(us_slopes,
                             data.frame(xbin = current_xbin,
                                        ybin = current_ybin,
                                        slope = slope$slope,
                                        intercept = slope$intercept,
                                          rsquared = slope$r_squared))
    }
    # if there are no valid events in that cell return null values
    else{
      us_slopes <- bind_rows(us_slopes,
                             data.frame(xbin = current_xbin,
                                        ybin = current_ybin,
                                        slope = NA,
                                        intercept = NA,
                                        rsquared = NA))
    }
  }
  # return null values
  else{
    us_slopes <- bind_rows(us_slopes,
                           data.frame(xbin = current_xbin,
                                      ybin = current_ybin,
                                      slope = NA,
                                      intercept = NA,
                                      rsquared = NA))
    }

  # Update progress bar
  setTxtProgressBar(pb, index)
}

# Close progress bar
close(pb)

# Load required library for plot arrangement
library(gridExtra)

save(us_slopes, file = "Data/Instrumental_slopes_D1.RData")
load("Data/instrumental_slopes_D1.RData")


# Slice weekly data for use on laptop
unique_cells <- unique(test_0.5[, c("bin.x", "bin.y")])
set.seed(123)
sampled_indices <- sample(1:nrow(unique_cells), 50)
sampled_cells <- unique_cells[sampled_indices, ]
weekly_PDSI_Sampled_0.5 <- test_0.5 %>%
  semi_join(sampled_cells, by = c("bin.x", "bin.y"))

save(weekly_PDSI_Sampled_0.5, file = "C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/SampledWeekly_0.5.RData")

# testing for weekly data
cropped_test_weekly <- crop.cell(test_0.5, xbin = -120.5, ybin = 40.5,
                                 pred_col = "predicted")
drought_events_test_weekly <- identify.weekly.drought(cropped_test_weekly)
list_test_weekly <- summarize.weekly.drought.events(drought_events_test_weekly,
                                             time_col = "Date", pred_col = "predicted")
drought_events_test2_weekly_fixed <- list_test_weekly[[1]]
drought_intervals_test_weekly <- list_test_weekly[[2]]

full_test_weekly <- evaluate.weekly.recurrence(full_modern_set_0.5, xbin = -113.5, ybin = 35.5,
                                 pred_col = "USDM_factor", time_col = "Date", 
                                 intensity_threshold = 4 )


full_drought_events_test2 <- full_test_weekly[[2]]
full_drought_intervals_test <- full_test_weekly[[3]]

plot_data_weekly <- plot.weekly.data(full_drought_events_test2)
weekly_slope <- recurrence.slope(plot_data_weekly)

test_plot <- plot.weekly.duration.v.return(plot_data_weekly)
test_plot
ggsave("Plots/modern_duration_v_return_D2.png", width = 10, height = 6)

save(plot_data, file = "Data/paleo_return_plot_data_D2.Rdata")


testplot <- recurrence.plot("instrumental", "D1", 2)
testplot


# difference plot and plot two lines on the histogram 


# histogram plotting 
us_slopes <- us_slopes %>% mutate(RI = slope*2 + intercept)
us_slopes$RI_yrs <- 10^us_slopes$RI

histo <- ggplot(us_slopes, aes(x=RI_yrs))+
  geom_histogram(color = "#244380", fill = "#2E55A3")+
  labs(title = "Instrumental Return Intervals of 2 Year Drought", 
       x = "Return Interval (Years)",
       y = "Count")+
  theme_minimal()+
  scale_x_log10()
  
histo

paleo_slopes <- us_slopes

plot_df <- bind_rows(
  us_slopes  %>% mutate(source = "Modern"),
  paleo_slopes %>% mutate(source = "Paleo")
)

ggplot(plot_df, aes(x = RI_yrs, fill = source)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 30, 
    color = "#4A4746"
  ) +
  scale_fill_manual(
    values = c(
      "Modern" = "#2E55A3",
      "Paleo" = "#A34E2E"
    )
  ) +
  scale_x_log10() +
  labs(
    title = "Return Intervals of 2-Year Drought Events >= D1",
    x = "Return Interval (Years)",
    y = "Count",
    fill = "Data Source"
  ) +
  theme_minimal()


ggsave("Plots/modern_v_paleo_2yr_D1_histo.png", width = 8, height = 2)
# power law
# functionalize approach http://127.0.0.1:10757/graphics/plot_zoom_png?width=1080&height=692
# weekly data making the same plots with the same gridcell

# random samples for 10 cells around the country- possibly put lines all on the same plot
# fun for everything and map out slopes

# calculate recurrence for one level of duration: ex. recurrence of 2 year long D2 level drought and plot across the whole country 

# highlight cell we are looking at 
us_slopes <- us_slopes %>% 
  mutate(correct_cell = (xbin == -113.5 & ybin == 35.5))

us_outline <- map_data("usa")

plot <- ggplot(us_slopes, aes(x = xbin, y = ybin, fill = correct_cell)) +
  geom_tile() +
  geom_path(data = us_outline, aes(x = long, y = lat, group = group), 
            color = "black", linewidth = 0.7, inherit.aes = FALSE) +
  theme_minimal()+
  scale_fill_manual(values = c("TRUE" = "#D4260F", 
                               "FALSE" = "grey"))+
  labs(title = "",
       x = "", 
       y = "")+
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "gray50"),
        legend.key.size = unit(0.8, "cm"),
        legend.key = element_rect(color = "gray50"),
        strip.text = element_text(face = "bold", size = 11))

plot

ggsave("Plots/highlighted_cell.png", width = 10, height = 6)

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

# Assume df and drought_events are already created from previous analysis
# For demonstration, using sample data structure

# =============================================================================
# PLOT 5: Severity-Duration-Frequency (SDF) Contour Plot
# =============================================================================

create_sdf_plot <- function(drought_events, study_period) {

  # Create grid of duration and intensity values
  duration_seq <- seq(1, max(drought_events$duration), length.out = 20)
  intensity_seq <- seq(1, 5, length.out = 20)

  # Calculate return intervals for each combination
  sdf_data <- expand.grid(duration = duration_seq, intensity = intensity_seq)

  sdf_data$return_interval <- sapply(1:nrow(sdf_data), function(i) {
    count <- sum(drought_events$duration >= sdf_data$duration[i] &
                 drought_events$max_intensity >= sdf_data$intensity[i])
    if (count == 0) return(NA)
    return(study_period / count)
  })

  plot <- ggplot(sdf_data, aes(x = duration, y = intensity, z = return_interval)) +
    geom_contour_filled(bins = 10) +
    geom_contour(color = "white", alpha = 0.5) +
    scale_y_continuous(
      breaks = 1:5,
      labels = c("D0", "D1", "D2", "D3", "D4")
    ) +
    labs(
      title = "Severity-Duration-Frequency (SDF) Curves",
      subtitle = "Contours show return interval in years",
      x = "Drought Duration (years)",
      y = "Drought Severity",
      fill = "Return\nInterval\n(years)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))

  ggsave("plot5_sdf_contour.png", plot, width = 10, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 6: Cumulative Drought Deficit
# =============================================================================

create_cumulative_deficit_plot <- function(pmdi_one_cell, drought_events) {

  # Calculate cumulative deficit for each drought event
  deficit_data <- pmdi_one_cell %>%
    filter(!is.na(drought_id)) %>%
    group_by(drought_id) %>%
    arrange(year) %>%
    mutate(
      # Only count negative PDSI values as deficit
      deficit = pmin(PDSI_Avg, 0),
      cumulative_deficit = cumsum(deficit)
    ) %>%
    ungroup()

  plot <- ggplot(deficit_data, aes(x = year, y = cumulative_deficit,
                                   group = drought_id, color = factor(drought_id))) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(
      title = "Cumulative Drought Deficit Over Time",
      subtitle = "Accumulated moisture deficit during each drought event",
      x = "Year",
      y = "Cumulative PDSI Deficit",
      color = "Drought\nEvent ID"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )

  ggsave("plot6_cumulative_deficit.png", plot, width = 10, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 7: Seasonal/Monthly Drought Onset Analysis
# =============================================================================

create_seasonal_onset_plot <- function(pmdi_one_cell) {

  # Identify drought onset years (first year of each drought event)
  onset_data <- pmdi_one_cell %>%
    filter(!is.na(drought_id)) %>%
    group_by(drought_id) %>%
    filter(year == min(year)) %>%
    ungroup()

  # If you have actual month data, use it. Otherwise simulate for demonstration
  # Assuming year data only - create distribution across months (random for demo)
  set.seed(123)
  onset_data$month <- sample(1:12, nrow(onset_data), replace = TRUE)

  monthly_counts <- onset_data %>%
    group_by(month) %>%
    summarise(count = n()) %>%
    complete(month = 1:12, fill = list(count = 0))

  monthly_counts$month_name <- factor(month.abb[monthly_counts$month],
                                     levels = month.abb)

  plot <- ggplot(monthly_counts, aes(x = month_name, y = count)) +
    geom_bar(stat = "identity", fill = "coral", color = "black") +
    labs(
      title = "Drought Onset by Month",
      subtitle = "When do droughts typically begin?",
      x = "Month",
      y = "Number of Drought Events",
      caption = "Note: Month data simulated for demonstration"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave("plot7_seasonal_onset.png", plot, width = 10, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 8: Return Period vs Intensity Curve
# =============================================================================

create_return_period_curve <- function(drought_events, study_period) {

  # Calculate return period for each intensity level
  intensity_levels <- 1:5
  return_data <- data.frame(
    intensity = intensity_levels,
    intensity_label = c("D0", "D1", "D2", "D3", "D4")
  )

  return_data$count <- sapply(intensity_levels, function(level) {
    sum(drought_events$max_intensity >= level)
  })

  return_data$return_period <- study_period / return_data$count
  return_data$return_period[return_data$count == 0] <- NA

  plot <- ggplot(return_data, aes(x = intensity, y = return_period)) +
    geom_line(color = "darkred", size = 1.5) +
    geom_point(size = 4, color = "darkred") +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("D0", "D1", "D2", "D3", "D4")
    ) +
    scale_y_log10() +
    labs(
      title = "Return Period vs Drought Intensity",
      subtitle = "How often do droughts of different severities occur?",
      x = "Drought Intensity",
      y = "Return Period (years, log scale)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave("plot8_return_period_curve.png", plot, width = 8, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 1: Exceedance Probability Curve
# =============================================================================

create_exceedance_probability_plot <- function(drought_events) {

  # Calculate exceedance probability for duration
  duration_sorted <- sort(drought_events$duration, decreasing = TRUE)
  n <- length(duration_sorted)
  exceedance_prob <- (1:n) / (n + 1) * 100

  duration_ep <- data.frame(
    value = duration_sorted,
    exceedance_prob = exceedance_prob,
    type = "Duration (years)"
  )

  # Calculate exceedance probability for intensity
  intensity_sorted <- sort(drought_events$max_intensity, decreasing = TRUE)
  intensity_ep <- data.frame(
    value = intensity_sorted,
    exceedance_prob = (1:length(intensity_sorted)) / (length(intensity_sorted) + 1) * 100,
    type = "Max Intensity"
  )

  ep_data <- rbind(duration_ep, intensity_ep)

  plot <- ggplot(ep_data, aes(x = value, y = exceedance_prob, color = type)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    facet_wrap(~type, scales = "free_x") +
    labs(
      title = "Exceedance Probability Curves",
      subtitle = "Probability that a drought will exceed given duration/intensity",
      x = "Value",
      y = "Exceedance Probability (%)",
      color = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold")
    )

  ggsave("plot1_exceedance_probability.png", plot, width = 10, height = 5, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 3: Drought Frequency by Decade
# =============================================================================

create_frequency_by_decade_plot <- function(drought_events) {

  # Assign each drought to a decade based on start year
  drought_events$decade <- floor(drought_events$start_year / 10) * 10

  decade_counts <- drought_events %>%
    group_by(decade) %>%
    summarise(
      total_droughts = n(),
      severe_droughts = sum(max_intensity >= 2)
    ) %>%
    pivot_longer(cols = c(total_droughts, severe_droughts),
                 names_to = "category",
                 values_to = "count")

  decade_counts$category <- factor(decade_counts$category,
                                  levels = c("total_droughts", "severe_droughts"),
                                  labels = c("All Droughts (≥D0)", "Severe Droughts (≥D2)"))

  plot <- ggplot(decade_counts, aes(x = decade, y = count, fill = category)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_fill_manual(values = c("All Droughts (≥D0)" = "lightblue",
                                 "Severe Droughts (≥D2)" = "darkred")) +
    labs(
      title = "Drought Frequency by Decade",
      subtitle = "Are droughts becoming more or less frequent over time?",
      x = "Decade",
      y = "Number of Drought Events",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )

  ggsave("plot3_frequency_by_decade.png", plot, width = 10, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# PLOT 2: Time Series with Drought Events Highlighted
# =============================================================================

create_timeseries_highlighted_plot <- function(pmdi_one_cell) {

  # Create color mapping for drought intensity
  pmdi_one_cell$color_group <- case_when(
    is.na(pmdi_one_cell$drought_id) ~ "No Drought",
    pmdi_one_cell$predictions == "D0" ~ "D0",
    pmdi_one_cell$predictions == "D1" ~ "D1",
    pmdi_one_cell$predictions == "D2" ~ "D2",
    pmdi_one_cell$predictions == "D3" ~ "D3",
    pmdi_one_cell$predictions == "D4" ~ "D4",
    TRUE ~ "No Drought"
  )

  pmdi_one_cell$color_group <- factor(pmdi_one_cell$color_group,
                          levels = c("No Drought", "D0", "D1", "D2", "D3", "D4"))

  plot <- ggplot(pmdi_one_cell, aes(x = year, y = PDSI_Avg)) +
    geom_line(color = "gray50", size = 0.5) +
    geom_point(aes(color = color_group), size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c("No Drought" = "gray70", "D0" = "yellow",
                "D1" = "orange", "D2" = "darkorange",
                "D3" = "red", "D4" = "darkred"),
      name = "Drought Category"
    ) +
    labs(
      title = "PDSI Time Series with Drought Events Highlighted",
      subtitle = "Color-coded by drought severity",
      x = "Year",
      y = "PDSI Average"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )

  ggsave("plot2_timeseries_highlighted.png", plot, width = 12, height = 6, dpi = 300)
  return(plot)
}

# =============================================================================
# EXECUTE ALL PLOTS (uncomment when you have full dataset)
# =============================================================================

cat("Drought Visualization Suite\n")
cat("============================\n\n")

# Assuming df and drought_events exist from previous analysis
# and study_period is defined (e.g., 2018)

# Uncomment these lines when running with your full dataset:
study_period <- 2018

cat("Creating Plot 5: SDF Contour...\n")
plot5 <- create_sdf_plot(drought_events, study_period)

cat("Creating Plot 6: Cumulative Deficit...\n")
plot6 <- create_cumulative_deficit_plot(pmdi_one_cell, drought_events)

cat("Creating Plot 7: Seasonal Onset...\n")
plot7 <- create_seasonal_onset_plot(pmdi_one_cell)

cat("Creating Plot 8: Return Period Curve...\n")
plot8 <- create_return_period_curve(drought_events, study_period)

cat("Creating Plot 1: Exceedance Probability...\n")
plot1 <- create_exceedance_probability_plot(drought_events)

cat("Creating Plot 3: Frequency by Decade...\n")
plot3 <- create_frequency_by_decade_plot(drought_events)

cat("Creating Plot 2: Time Series Highlighted...\n")
plot2 <- create_timeseries_highlighted_plot(pmdi_one_cell)

ggsave("C:/Users/delil/Desktop/NAU/Research 2024-2025/PDSI_research_2024_2025/Plots/pred.vs.actual_pmdi_1000.png", test.biplot)

pmdi_1000 <- pmdi_prediction_set %>% filter(year == 1000)

test.plot <- plot.pmdi(pmdi_1000, "predictions", year = 1000)

test2 <- plot.pdsi(pmdi_1000, "PDSI_Avg", year = 1000)

test.biplot <- test.plot + test2
# plot predicted vs actual
