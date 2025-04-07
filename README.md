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
  ```
- For single-run trajectory visualizations, set:

  ```r
  pln <- 1

  **Note**: Ensure you run enough repetitions (`Nruns >= 10`) for robust statistical results.

## Code Snippets for Each Objective

### Objective 1: Distance vs. Place Model Performance

#### Extracting Performance Metrics (Example: Latency)

```r
library(tidyr)
library(dplyr)

# For the place-cell-only model
latency_place <- PMs_place[1,,,]
df_place <- data.frame(latency_place)
df_place$Day <- 1:nrow(df_place)
df_place_long <- gather(df_place, "Trial", "Latency", -Day)
df_place_long$Model <- "place cell"

# For the distance-cell-only model
latency_distance <- PMs_distance[1,,,]
df_distance <- data.frame(latency_distance)
df_distance$Day <- 1:nrow(df_distance)
df_distance_long <- gather(df_distance, "Trial", "Latency", -Day)
df_distance_long$Model <- "distance cell"

# Combine the two datasets
df_final <- rbind(df_place_long, df_distance_long)
