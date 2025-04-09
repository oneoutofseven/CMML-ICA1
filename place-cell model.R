run_trial <- function(weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y,
                      Vdecay, ac_const, beta, etdecay, lrate, discf,
                      noise, platform_x, platform_y, starting_x, starting_y,
                      speed, wall_pun)
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  pool_diameter    <- 1.4    # Maze diameter in meters
  platform_radius <- 0.06    # Platform radius in meters

  N_pc <- 211  # Number of place cells
  N_ac <- 36   # Number of action cells

  # Initialize some counters / variables
  which     <- 0        # Index of the currently selected action
  dist      <- 0        # Total travel distance
  wall_zone <- 0        # Count of time steps spent near the wall
  quadrants <- c(0,0,0,0)  # Count of time steps spent in each of the 4 quadrants

  # Copy the initial weights
  weights <- weights0
  
  # Initialize eligibility traces
  el_tr <- matrix(0, nrow=N_pc, ncol=N_ac)

  # Initialize trajectory variables
  track_x <- starting_x
  track_y <- starting_y
  vel_x   <- 0
  vel_y   <- 0

  # NAVIGATION LOOP
  # Continue until the agent is within the platform radius
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 >
         platform_radius^2)
  {
    # Apply a bit of random noise to weights
    weights <- weights * (1 - noise) + matrix(runif(N_pc*N_ac), nrow=N_pc) * Wmult * noise

    # Compute place-cell activation
    PC_activation <- numeric(N_pc)
    for (i in 1:N_pc){
      dx <- track_x[length(track_x)] - PC_x[i]
      dy <- track_y[length(track_y)] - PC_y[i]
      PC_activation[i] <- exp(- (dx^2 + dy^2) / (2 * sigma_pc^2))
    }

    # Store the previous Q-value (for the chosen action) if this is not the first step
    if (length(track_x) > 1) {
      prevQ <- AC_activation[which] 
    }

    # Compute the Q-values (action-cell activation) by multiplying place-cell activation with the weights
    AC_activation <- as.vector(PC_activation %*% weights)

    # --- SAFETY CHECKS to avoid zero or NaN sums ---
    # Replace any NaN in AC_activation with 0, and ensure they're not negative
    AC_activation[is.na(AC_activation)] <- 0
    AC_activation <- pmax(AC_activation, 0)

    # Raise AC_activation to the power of beta
    ACsel <- AC_activation^beta
    sAC   <- sum(ACsel)

    # If the sum is zero, NaN, or Inf, revert to a uniform distribution
    if (sAC <= 0 || is.na(sAC) || is.infinite(sAC)) {
      ACsel <- rep(1/N_ac, N_ac)
    } else {
      ACsel <- ACsel / sAC
    }

    # Select an action stochastically based on ACsel
    ASrand <- runif(1)
    which  <- 1
    ASsum  <- ACsel[1]
    while (which < N_ac && ASsum < ASrand) {
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }

    # Update eligibility traces
    el_tr <- el_tr * etdecay
    for (j in 1:N_ac){
      # Gaussian-like spreading of eligibility over neighboring actions
      itmp    <- min(abs(j - which), N_ac - abs(j - which))
      actgaus <- exp(- (itmp^2) / (2 * sigma_ac^2))
      el_tr[, j] <- el_tr[, j] + actgaus * AC_activation[j] * PC_activation
    }

    # Update velocity
    new_vel_x <- (vel_x[length(vel_x)] + ac_const * cos(which / N_ac * 2 * pi)) * Vdecay
    new_vel_y <- (vel_y[length(vel_y)] + ac_const * sin(which / N_ac * 2 * pi)) * Vdecay
    vel_x <- c(vel_x, new_vel_x)
    vel_y <- c(vel_y, new_vel_y)

    # Update position
    new_x <- track_x[length(track_x)] + new_vel_x
    new_y <- track_y[length(track_y)] + new_vel_y
    track_x <- c(track_x, new_x)
    track_y <- c(track_y, new_y)

    # If out of bounds (outside the maze), push back to the edge
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter / 2)^2) {
      ratio <- (track_x[length(track_x)]^2 + track_y[length(track_y)]^2) / ((pool_diameter / 2)^2)
      track_x[length(track_x)] <- track_x[length(track_x)] / sqrt(ratio)
      track_y[length(track_y)] <- track_y[length(track_y)] / sqrt(ratio)
      vel_x[length(vel_x)]     <- track_x[length(track_x)] - track_x[length(track_x) - 1]
      vel_y[length(vel_y)]     <- track_y[length(track_y)] - track_y[length(track_y) - 1]
    }

    # Compute reward and update weights
    if (length(track_x) > 2) {
      # Check if the platform is found
      if ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 <
          platform_radius^2) {
        rew <- 10  # Found the platform
      } else if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 >
                 (0.99 * pool_diameter / 2)^2) {
        rew <- -wall_pun  # Hit the wall
      } else {
        rew <- 0
      }

      currQ  <- AC_activation[which]
      tderr  <- rew + discf * currQ - prevQ  # Temporal difference error
      weights <- pmax(weights + lrate * tderr * el_tr, 0)
    }

    # Update total traveled distance
    last_step_dist <- sqrt(
      (track_x[length(track_x)] - track_x[length(track_x) - 1])^2 +
      (track_y[length(track_y)] - track_y[length(track_y) - 1])^2
    )
    dist <- dist + last_step_dist

    # Record time step in wall zone or in a certain quadrant
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 >
        0.8 * (pool_diameter / 2)^2) {
      wall_zone <- wall_zone + 1
    } else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0) {
      quadrants[1] <- quadrants[1] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0) {
      quadrants[2] <- quadrants[2] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0) {
      quadrants[3] <- quadrants[3] + 1
    } else {
      quadrants[4] <- quadrants[4] + 1
    }

    # If the agent has been navigating for a long time (>60s), stop
    if (length(track_x) > 100) {
      # Average speed in meters / time step
      speed_ts <- mean(sqrt(vel_x[-1]^2 + vel_y[-1]^2))
      # Convert number of steps to seconds
      latency <- (length(track_x) - 1) * speed_ts / speed
      if (latency > 60) {
        break
      }
    }
  } # End of while loop

  # Compute final latency in seconds
  num_steps <- length(track_x) - 1
  speed_ts  <- mean(sqrt(vel_x[-1]^2 + vel_y[-1]^2))
  latency   <- num_steps * speed_ts / speed

  # Convert wall_zone and quadrants to fractions of total steps
  wall_zone <- wall_zone / num_steps
  quadrants <- quadrants / num_steps

  # Return key results
  return(list(
    weights,     # updated weights
    track_x,     # trajectory X
    track_y,     # trajectory Y
    vel_x,       # velocity X over time
    vel_y,       # velocity Y over time
    dist,        # total distance
    wall_zone,   # fraction of time steps near the wall
    quadrants,   # fraction of time steps in each quadrant
    latency      # total latency in seconds
  ))
}
  
N_pc <- 211 #Population of place cells [100..300]
N_ac <- 36 #Population of action cells [25..50]

plot_trajectories <- 0 #yes - 1, no - 0
plot_cognitive_maps <- 0 #yes - 1, no - 0
pln <- plot_trajectories + plot_cognitive_maps
Nruns <- 50 #how many runs to run if not plotting anything

pool_diameter <- 1.4 #Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 #Platform radius (m)
sigma_pc <- 0.1 #place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- 2 #action cell sigma (standard deviation), in action cells [1..3]

etdecay <- 0.83 #Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- 6 #Exploration-exploitation factor (\beta) [0.5..12]
alpha <- 0.01 #Learning rate (\alpha) [0.005..0.02]
gamma <- 0.85 #Discount factor (\gamma) [0.75..0.95]

Vdecay <- 0.82 #velocity decay [0.75..0.95]
ac_const <- 0.02 #acceleration const [0.01..0.03]
Wnoise <- 0.0004 #Weight noise [0.0001..0.0007]
Wmult <- 0.1 #Weight multiplier [0.05..0.15]
hitwall <- 0.5 #punishment for hitting the wall [0..1]
speed <- 0.175 #mouse speed (m/s) [0.1..0.25]

Ntrials <- 4 #number of trials per day
Ndays <- 8 #number of days

#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
if (pln > 0.5) #if any plots
{ PMs <- array(rep(0,5*Ndays*Ntrials), c(5,Ndays,Ntrials))
} else {
  PMs <- array(rep(0,5*Ndays*Ntrials*Nruns), c(5,Ndays,Ntrials,Nruns)) }
#multiple runs


#Platform coordinates:
platform_x <- cos(-pi/4)*pool_diameter/4 #x coordinate
platform_y <- sin(-pi/4)*pool_diameter/4 #y coordinate

#Starting locations of the modeled animal (4 different ones):
strad <- pool_diameter/2*0.85 #15% of maze radius to the wall
starting_xs <- strad * c(cos(pi/6), cos(pi/3), cos(7*pi/6), cos(4*pi/3)) #x coordinates
starting_ys <- strad * c(sin(pi/6), sin(pi/3), sin(7*pi/6), sin(4*pi/3)) #y coordinates

th <- (0:100)/50*pi #for plotting circles :)

if (pln > 0.5) {

#Generate initial weights
weights <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult

#Generate place cells
PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell
PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell
for (i in 1:N_pc) {
    #For each place cell:
    PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
    PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
        #Checks for out of bounds
        PC_x[i] <- (runif(1) - 0.5)*pool_diameter
        PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    }
}

for (day in 1:Ndays) {
    idxs = sample(4) #randomly choose 4 starting locations
    for (trial in 1:Ntrials){

        idx <- idxs[trial] #take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]

        modresults <- run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x, platform_y, starting_x, starting_y, speed, hitwall)
        #run trial
        weights <- modresults[[1]]
        track_x <- modresults[[2]]
        track_y <- modresults[[3]]
        vel_x <- modresults[[4]]
        vel_y <- modresults[[5]]
        
#        weights <- wres

        PMs[1,day,trial] <- modresults[[9]] #latency
        PMs[2,day,trial] <- modresults[[6]] #dist
        PMs[3,day,trial] <- modresults[[8]][4]*100 #target quadrant
        PMs[4,day,trial] <- modresults[[8]][2]*100 #opposite quadrant
        PMs[5,day,trial] <- modresults[[7]]*100 #wall zone
        #record performance measures

        if (plot_trajectories)
        {
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l", xlab = paste("day ",day,", trial ",trial), ylab = "trajectory")
          #plot the trajectory
          lines(track_x, track_y, type = "l")
          #plot the platform
          lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
        }

        if (plot_cognitive_maps){
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day ",day,", trial ",trial), ylab = "cognitive map")
          #plot the cognitive map
          for (x in (-3:3)*(pool_diameter/6)){
            for (y in (-3:3)*(pool_diameter/6)){
                if (x^2 + y^2 <= (pool_diameter/2)^2){
                    x2 = x
                    y2 = y
                    for (k in 1:N_ac){
                        PC_activation <- rep(0,N_pc)
                        for (i in 1:N_pc){
                            PC_activation[i] <- exp(-((x - PC_x[i])^2 + (y - PC_y[i])^2)/(2*sigma_pc^2))
                        }
                #Calculate AC activation (i.e. value of the action)
                        AC_activation <- rep(0,N_ac)
                        for (i in 1:N_ac){
                            for (j in 1:N_pc){
                                AC_activation[i] <- AC_activation[i] + PC_activation[j]*weights[j,i]
                            }
                        }
                         x2 <- c(x2, x + (AC_activation[k]/10)*cos(k/N_ac*2*pi))
                         y2 <- c(y2, y + (AC_activation[k]/10)*sin(k/N_ac*2*pi))
        #                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
                    }
                    lines(x2,y2,type = "l",col = "blue")
               }
            }
          }
          # plot the platform
          lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
        }
    }
}

} else {
# run multiple times without plotting!

for (reps in 1:Nruns){
#Generate initial weights for each run
weights <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult

#Generate place cells for each run
PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell
PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell
for (i in 1:N_pc) {
  #For each place cell:
  PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
  PC_y[i] <- (runif(1) - 0.5)*pool_diameter
  while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
    #Checks for out of bounds
    PC_x[i] <- (runif(1) - 0.5)*pool_diameter
    PC_y[i] <- (runif(1) - 0.5)*pool_diameter
  }
}

for (day in 1:Ndays){
  idxs = sample(4)  #randomly choose 4 starting locations
    for (trial in 1:Ntrials){

      idx <- idxs[trial] #take each location
      starting_x <- starting_xs[idx]
      starting_y <- starting_ys[idx]

      modresults <- run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x, platform_y, starting_x, starting_y, speed, hitwall)
      #run trial
      weights <- modresults[[1]]
      
      PMs[1,day,trial,reps] <- modresults[[9]] #latency
      PMs[2,day,trial,reps] <- modresults[[6]] #dist
      PMs[3,day,trial,reps] <- modresults[[8]][4]*100 #target quadrant
      PMs[4,day,trial,reps] <- modresults[[8]][2]*100 #opposite quadrant
      PMs[5,day,trial,reps] <- modresults[[7]]*100 #wall zone
      #record performance measures
     }
}
}
}
