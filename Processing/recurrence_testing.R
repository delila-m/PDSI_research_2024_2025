library(dplyr)
library(tidyverse)
# load in drought functions
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Functions/") #glg
setwd("C:/Users/delil/Desktop/NAU/Research 2024-2025/PDSI_research_2024_2025/Functions/") #laptop
source("drought_functions.R")

# file path on laptop
load("C:/Users/delil/Desktop/NAU/Research 2024-2025/PDSI_research_2024_2025/Data/PMDIPredictionSet.Rdata")
# file path in glg
load("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/PMDIPredictionSet.Rdata")

# Previous categorical model for comparison 
setwd("C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/UpdatedCleaned_0.5")
load("RFAnalysis0.5_factor_updated.Rdata")

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
cropped_test <- crop.cell(pmdi_prediction_set, xbin = -113.75, ybin = 35.25, pred_col = "predictions")
drought_events_test <- identify.drought(cropped_test)
list_test <- summarize.drought.events(drought_events_test, time_col = "year")
drought_events_test2 <- list_test[[1]]
drought_intervals_test <- list_test[[2]]

full_test <- evaluate.recurrence(pmdi_prediction_set, xbin = -113.75, ybin = 35.25, 
                                 pred_col = "predictions", time_col = "year")


full_drought_events_test2 <- full_test[[2]]
full_drought_intervals_test <- full_test[[3]]

test_plot <- plot.duration.v.return(full_test$Drought_events)

# Load required library for plot arrangement
library(gridExtra)

# Get unique cell combinations
unique_cells <- unique(pmdi_prediction_set[, c("bin.x", "bin.y")])

# Randomly sample 9 cells
set.seed(123)  
sampled_indices <- sample(1:nrow(unique_cells), 9)
sampled_cells <- unique_cells[sampled_indices, ]

# Initialize list to store plots
plot_list <- list()
all_drought_data <- list()


# Loop through each sampled cell
for(i in 1:9) {
  # Get cell coordinates
  xbin <- sampled_cells$bin.x[i]
  ybin <- sampled_cells$bin.y[i]
  
  # Run your evaluation function
  full_test <- evaluate.recurrence(pmdi_prediction_set, 
                                   xbin = xbin, 
                                   ybin = ybin, 
                                   pred_col = "predictions", 
                                   time_col = "year")
  
  # Generate plot
  test_plot <- plot.duration.v.return(full_test$Drought_events)
  
  # Store plot in list
  plot_list[[i]] <- test_plot

  # print progress
  cat("Processed cell", i, "of 9: (", xbin, ",", ybin, ")\n")
}

# Arrange all plots in a 3x3 grid
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)


# Use the function - you may want to increase n_cells to ensure you get 9 successful ones
result <- plot.duration.v.return.combined(train_0.5, n_cells = 15, 
                                          seed = 123,  
                                          pred_col = "predicted", 
                                          time_col = "Date")

# Display the plot
print(result$plot)

# Slice weekly data for use on laptop
unique_cells <- unique(test_0.5[, c("bin.x", "bin.y")])
set.seed(123)
sampled_indices <- sample(1:nrow(unique_cells), 50)
sampled_cells <- unique_cells[sampled_indices, ]
weekly_PDSI_Sampled_0.5 <- test_0.5 %>%
  semi_join(sampled_cells, by = c("bin.x", "bin.y")) 

save(weekly_PDSI_Sampled_0.5, file = "C:/Users/dgm239/Downloads/Research_2025/PDSI_research_2024/Data/SampledWeekly_0.5.RData")

# testing for weekly data 
cropped_test_weekly <- crop.cell(sampled_cells, xbin = -113.5, ybin = 35.5, pred_col = "predicted")
drought_events_test <- identify.drought(cropped_test_weekly)
list_test <- summarize.drought.events(drought_events_test, time_col = "Date", pred_col = "predicted")
drought_events_test2 <- list_test[[1]]
drought_intervals_test <- list_test[[2]]

full_test <- evaluate.recurrence(test_0.5, xbin = -113.75, ybin = 35.25, 
                                 pred_col = "predicted", time_col = "Date")


full_drought_events_test2 <- full_test[[2]]
full_drought_intervals_test <- full_test[[3]]

test_plot <- plot.duration.v.return(list_test$Droughts)


# power law 
# functionalize approach 
# weekly data making the same plots with the same gridcell 

# random samples for 10 cells around the country- possibly put lines all on the same plot 
# fun for everything and map out slopes 





ggsave("")

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
