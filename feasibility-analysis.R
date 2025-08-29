# feasibility-analysis.R
# author: Kim P.C. van der Sanden
# date: 12/08/2025

# Load required packages
library(dplyr)    # For data manipulation
library(ggplot2)  # For data visualization

# Load and inspect data
# data_full contains raw EMA data for all participants
data_full <- readRDS("C:/path/to/data.rds")

#-------------------------------------------------------------------------------
# Feasibility analysis
#-------------------------------------------------------------------------------

# Report the number of unique participants and total observations
num_participants <- length(unique(data_full$subjno))
num_observations <- nrow(data_full)

# Create a subset for compliance analysis based on Fat (physical fatigue)
# This assumes that if Fat is not NA, the prompt was filled out
data_compliance <- data_full %>%
  select(subjno, Daynumber, Fat) %>%
  mutate(
    is_filled_out = ifelse(!is.na(Fat), 1, 0)
  )

# Compliance analysis
# Aggregate filled-out responses per day for each participant
daily_compliance <- data_compliance %>%
  group_by(subjno, Daynumber) %>%
  summarise(
    filled_in_count = sum(is_filled_out, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate daily percentage based on 5 expected responses per day
    daily_percentage = (filled_in_count / 5) * 100
  )

# Calculate the mean compliance percentage for each participant over all days
participant_avg_compliance <- daily_compliance %>%
  group_by(subjno) %>%
  summarize(
    mean_percentage = mean(daily_percentage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_percentage))

# Join the average compliance back to the daily data
daily_compliance <- daily_compliance %>%
  left_join(
    participant_avg_compliance %>% rename(avg_compliance = mean_percentage),
    by = "subjno"
  )

# Create anonymous labels for each participant for plotting
anon_map <- participant_avg_compliance %>%
  mutate(anon_label = paste0("P", seq_len(n())))

daily_compliance <- daily_compliance %>%
  left_join(anon_map, by = "subjno")

# Assign participants to cohorts based on their original ID
# cohort_map <- c(
#  ... # sanitized
# )

daily_compliance <- daily_compliance %>%
  mutate(
    cohort = factor(cohort_map[as.character(subjno)], levels = c("PRO-Lung", "PRO-RCC"))
  )

# Statistical reporting
# Create compliance groups based on a 75% threshold.
daily_compliance <- daily_compliance %>%
  mutate(compliance_group = case_when(
    avg_compliance < 75 ~ "low",
    TRUE ~ "high"
  ))

# Calculate and report daily mean compliance by group.
group_compliance_stats <- daily_compliance %>%
  group_by(anon_label, compliance_group) %>%
  summarise(mean_beeps = mean(filled_in_count, na.rm = TRUE), .groups = "drop") %>%
  group_by(compliance_group) %>%
  summarise(
    mean_compliance = mean(mean_beeps, na.rm = TRUE),
    sd_compliance = sd(mean_beeps, na.rm = TRUE),
    n_participants = n()
  )
print(group_compliance_stats)

# Calculate and report compliance statistics for specific day ranges.
day1_stats <- daily_compliance %>%
  filter(Daynumber == 1) %>%
  summarise(
    Mday1 = mean(filled_in_count, na.rm = TRUE),
    SDday1 = sd(filled_in_count, na.rm = TRUE)
  )

day2_20_stats <- daily_compliance %>%
  filter(Daynumber >= 2, Daynumber <= 20) %>%
  group_by(anon_label) %>%
  summarise(mean_beeps = mean(filled_in_count, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    Mday2_20 = mean(mean_beeps, na.rm = TRUE),
    SDday2_20 = sd(mean_beeps, na.rm = TRUE)
  )

day21_stats <- daily_compliance %>%
  filter(Daynumber == 21) %>%
  summarise(
    Mday21 = mean(filled_in_count, na.rm = TRUE),
    SDday21 = sd(filled_in_count, na.rm = TRUE)
  )

cat(sprintf("Day 1 Mean: %.2f, SD: %.2f\n", day1_stats$Mday1, day1_stats$SDday1))
cat(sprintf("Days 2-20 Mean: %.2f, SD: %.2f\n", day2_20_stats$Mday2_20, day2_20_stats$SDday2_20))
cat(sprintf("Day 21 Mean: %.2f, SD: %.2f\n", day21_stats$Mday21, day21_stats$SDday21))


# Heatmap of daily compliance
# Define the levels for the y-axis to ensure consistent sorting
anon_levels <- anon_map$anon_label
p <- ggplot(daily_compliance, aes(x = Daynumber, y = factor(anon_label, levels = rev(anon_levels)), fill = daily_percentage)) +
  # Use geom_tile to create the heatmap squares.
  geom_tile(width = 0.9, height = 0.7, color = "white") +
  # Custom color scale to represent compliance percentages
  scale_fill_gradientn(
    colors = c("#FF0D0D", "#FF4E11", "#FAB733", "#ACB334", "#69B34C"),
    values = scales::rescale(c(0, 20, 40, 60, 80, 100)),
    breaks = c(0, 20, 40, 60, 80, 100),
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
    name = "Completed out of 5 (%)"
  ) +
  # Customize x-axis to show specific day numbers
  scale_x_continuous(
    limits = c(0.5, 21.5),
    expand = c(0.02, 0),
    breaks = 1:21,
    labels = ifelse(1:21 %in% c(1,5,10,15,20), 1:21, "")
  ) +
  # Customize y-axis to use anonymous labels instead of the sensitive participant labels
  scale_y_discrete(
    breaks = anon_levels,
    labels = ifelse(anon_levels %in% c("P1", "P5", "P10", "P15", "P20", "P25", "P30", "P35", "P40", "P46"), anon_levels, ""),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Day number",
    y = "Patient"
  ) +
  # Apply a clean, minimal theme
  theme_bw() +
  theme(
    legend.position = "hide",
    legend.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
print(p)


