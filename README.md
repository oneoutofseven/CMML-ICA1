# Comparing Computational Models in the Morris Water Maze (MWM)

## 1. Project Goals

1. **Objective 1:** Determine **if and when** the distance‐cell model can outperform the place‐cell model under **fixed‐platform** conditions.  
2. **Objective 2:** Evaluate the performance of the **combined place‐distance cell model** across various balance ratios, compared to place‐only and distance‐only models.

\bigskip

## 2. Important Notes Before Running

- **Set** `pln <- 0` **for multiple repetitions and reliable statistical comparisons**:
  ```r
  pln <- 0
  Nruns <- 50
  ```
- For single-run visualizations, pln = 1, 2, or 3
- In distance‐cell model, set the **variable_platform <- 0**. In the combined model, ensure that the platform is fixed by setting **whichplatform = 1**.


  **Note**: Ensure you run enough repetitions (`Nruns >= 50`) for robust statistical results.
---
## 3. Code Snippets for Each Objective

### Objective 1: Distance vs. Place Model Performance (Using "place-cell model.R" and "distance-cell model.R") 

#### Step 1: Run "place-cell model.R"
#### Step 2: Extracting Performance Metrics (Example: Latency)
```r
library(tidyr)
library(dplyr)

# For the place-cell-only model (Use it after runing "place-cell model.R")
PMs_place <- PMs
latency_place <- PMs_place[1,,,] # In the model, it only has PMs. So, you can use "PMs_place <- PMs" to store the data.
df_place <- data.frame(latency_place)
df_place$Day <- 1:nrow(df_place)
df_place_long <- gather(df_place, "Trial", "Latency", -Day)
df_place_long$Model <- "place cell"
```

#### Step 3: Run "distance-cell model.R"
#### Step 4: Extracting Performance Metrics (Example: Latency)
```r
# For the distance-cell-only model (Use it after runing "place-cell model.R")
PMs_distance <- PMs
latency_distance <- PMs_distance[1,,,] # In the model, it only has PMs. So, you can use "PMs_distance <- PMs" to store the data.
df_distance <- data.frame(latency_distance)
df_distance$Day <- 1:nrow(df_distance)
df_distance_long <- gather(df_distance, "Trial", "Latency", -Day)
df_distance_long$Model <- "distance cell"

# Combine the two datasets
df_final <- rbind(df_place_long, df_distance_long)
```

#### Step 5: Plotting the Performance
```r
library(ggplot2)

summary_df <- df_final %>%
  group_by(Day, Model) %>%
  summarise(
    Mean = mean(Latency, na.rm = TRUE),
    SE   = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_df, aes(x = Day, y = Mean, color = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = "Day", y = "Latency (s)") +
  theme_classic()
```

#### Step 6: Statistical Analysis
```r
anova_latency <- aov(Latency ~ as.factor(Day) * Model, data = df_final)
summary(anova_latency)

TukeyHSD(anova_latency)
```

**Note**: For "target quadrant" and "wall zone" time, change PM[1,,,] to PM[3,,,] or PM[5,,,]

---
### Objective 2: Combined Model vs. Separate Models (Using "place distance cells combined model.R") 

#### Step 1: Run "place distance cells combined model.R"
#### Step 2: Preparing Data for the Combined Model (Example: Ratio 0.25:0.75)
```r
# For the combined model with 25% distance cell signals and 75% place cell signals
PMs_combined25[1,,,,] <- PMs[1,,,,]
latency_combined25 <- PMs_combined25[1,,,,] # In the model, it only has PMs. So, you can use PMs_combined25 <- PMs to store the data.
df_combined25 <- data.frame(latency_combined25)
df_combined25$Day <- 1:nrow(df_combined25)
df_combined25_long <- gather(df_combined25, "Trial", "Latency", -Day)
df_combined25_long$Model <- "wall:place=0.25:0.75"
```
#### Step 3: Repeat similar steps for other ratios (0:1, 0.5:0.5, 0.75:0.25, 1:0) and label them accordingly.  
**Note**: The weight of distance cells can be changed by altering **weight_wall <- 0.5**

#### Step 4: Combine All Models for Comparison
```r
df_all <- rbind(
  df_place_long,
  df_distance_long,
  df_combined25_long,
  df_combined50_long,
  df_combined75_long
)
```

#### Step 5: Plotting Comparison
```r
summary_combined <- df_all %>%
  group_by(Day, Model) %>%
  summarise(
    Mean = mean(Latency, na.rm = TRUE),
    SE   = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_combined, aes(x = Day, y = Mean, color = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = "Day", y = "Latency (s)") +
  theme_classic()
```

#### Step 6: Statistical Analysis
```r
anova_latency <- aov(Latency ~ as.factor(Day) * Model, data = df_all)
summary(anova_latency)

TukeyHSD(anova_latency)
```
**Note**: For "target quadrant" and "wall zone" time, change PM[1,,,,] to PM[3,,,,] or PM[5,,,,]

---
## 4. Interpreting Results

- **Latency**: Lower latency implies improved learning.
- **Target Quadrant Time**: Higher occupancy indicates better spatial memory.
- **Wall Zone Time**: Lower values suggest less reliance on boundary cues.
