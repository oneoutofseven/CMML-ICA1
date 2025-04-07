# README: Comparing Computational Models in the Morris Water Maze (MWM)

## Project Goals

1. **Objective 1:** Determine **if and when** the distance‐cell model can outperform the place‐cell model under **fixed‐platform** conditions.  
2. **Objective 2:** Evaluate the performance of the **combined place‐distance cell model** across various balance ratios, compared to place‐only and distance‐only models.

---

## Important Notes Before Running

- **Set** `pln <- 0` **for multiple repetitions and reliable statistical comparisons**:
  ```r
  pln <- 0
  Nruns <- 10
For single‐run trajectory visualizations, set:

r
复制
编辑
pln <- 1
Code Snippets for Each Objective
Objective 1: Distance vs. Place Model Performance
Extract performance metrics (example for latency)
r
复制
编辑
library(tidyr)
library(dplyr)

# Place model
latency_place <- PMs_place[1,,,]
df_place <- data.frame(latency_place)
df_place$Day <- 1:nrow(df_place)
df_place_long <- gather(df_place, "Trial", "Latency", -Day)
df_place_long$Model <- "place cell"

# Distance model
latency_distance <- PMs_distance[1,,,]
df_distance <- data.frame(latency_distance)
df_distance$Day <- 1:nrow(df_distance)
df_distance_long <- gather(df_distance, "Trial", "Latency", -Day)
df_distance_long$Model <- "distance cell"

# Combine both models
df_final <- rbind(df_place_long, df_distance_long)
Plot performance
r
复制
编辑
library(ggplot2)

summary_df <- df_final %>%
  group_by(Day, Model) %>%
  summarise(
    Mean = mean(Latency, na.rm = TRUE),
    SE   = sd(Latency, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_df, aes(Day, Mean, color = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = "Day", y = "Latency (s)") +
  theme_classic()
ANOVA statistical comparison
r
复制
编辑
anova_result <- aov(Latency ~ Model * factor(Day), data = df_final)
summary(anova_result)
Objective 2: Combined Model vs. Separate Models
Prepare combined models data (example ratio 0.25:0.75)
r
复制
编辑
# Combined model (25% distance, 75% place)
latency_combined25 <- PMs_combined25[1,,,,]
df_combined25 <- data.frame(latency_combined25)
df_combined25$Day <- 1:nrow(df_combined25)
df_combined25_long <- gather(df_combined25, "Trial", "Latency", -Day)
df_combined25_long$Model <- "wall:place=0.25:0.75"
(Repeat above for other ratios, e.g., 0.5:0.5, 0.75:0.25.)

Combine all data for plotting
r
复制
编辑
df_all <- rbind(
  df_place_long,
  df_distance_long,
  df_combined25_long,
  df_combined50_long,
  df_combined75_long
)
Plot comparison
r
复制
编辑
summary_combined <- df_all %>%
  group_by(Day, Model) %>%
  summarise(
    Mean = mean(Latency, na.rm=TRUE),
    SE   = sd(Latency, na.rm=TRUE)/sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_combined, aes(Day, Mean, color = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width=0.2) +
  labs(x = "Day", y = "Latency (s)") +
  theme_classic()
Statistical analysis (ANOVA and Tukey’s post‐hoc)
r
复制
编辑
df_all$Day <- factor(df_all$Day)
anova_combined <- aov(Latency ~ Model, data = df_all)
summary(anova_combined)

TukeyHSD(anova_combined)
Interpreting Results
Latency: Lower latency implies improved learning.

Target Quadrant Time: Higher occupancy indicates better spatial memory.

Wall Zone Time: Lower values suggest less reliance on boundary cues.

Final Reminder
Ensure you run multiple repetitions (Nruns >= 10) for statistically robust conclusions. Always check statistical assumptions (e.g., normality, homogeneity of variances) before analyzing results.

