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
