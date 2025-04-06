# ======================= run trial function =======================
run_trial <- function(weights0_pc, weights0_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun, weight_wall, weight_place)# c
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  
  pool_diameter <- 1.4 #Maze diameter in metres (m)
  platform_radius <- 0.06 #Platform radius
  
  N_dc <- 10 #Population of place cells
  N_pc <- 211 #Population of place cells # c
  N_ac <- 36 #Population of action cells
  which_pc <- 0 # c
  which_dc <- 0 # c
  which <- 0 # c
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0,0,0,0) #Percentage spent on each quadrant
  
  weights_pc <- weights0_pc #Initialize modifiable weights # c
  weights_dc <- weights0_dc #Initialize modifiable weights # c
  
  
  el_tr_dc <- matrix(rep(0, N_dc*N_ac), nrow = N_dc) #Initialize eligibility traces matrix # c
  el_tr_pc <- matrix(rep(0, N_pc*N_ac), nrow = N_pc) #Initialize eligibility traces matrix # c
  
  #Initialize trajectories
  track_x <- starting_x #Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) # while out of platform
  {
    
    weights_dc <- weights_dc*(1-noise) + matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult*noise # c
    weights_pc <- weights_pc*(1-noise) + matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult*noise # c
    
    dist.to.wall <- pool_diameter/2 - sqrt(track_x[length(track_x)]^2+track_y[length(track_y)]^2)
    
    #Calculate DC and PC activation
    DC_activation <- rep(0, N_dc)
    for (i in 1:N_dc){
      DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
    }
    PC_activation <- rep(0, N_pc) # c
    for (i in 1:N_pc){# c
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2)/(2*sigma_pc^2))
    }# c
    
    #Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1){
      prevQ_pc <- AC_activation_pc[which_pc] #Displays the Q value before movement # c
      prevQ_dc <- AC_activation_dc[which_dc] #Displays the Q value before movement # c
    }
    #print(DC_activation)
    AC_activation_dc <- DC_activation %*% weights_dc # c
    AC_activation_pc <- PC_activation %*% weights_pc # c
    
    
    #Make an action
    ACsel_pc <- pmax(0,AC_activation_pc)^beta# c
    ACsel_pc <- ACsel_pc / sum(ACsel_pc)# c
    ACsel_dc <- pmax(0,AC_activation_dc)^beta# c
    ACsel_dc <- ACsel_dc / sum(ACsel_dc)# c
    
    
    shift = round((180-atan2(track_y[length(track_y)],track_x[length(track_x)])/pi*180)/10)
    if(shift == 36){
      shift = 0
    }
    
    
    if (shift < 0){
      ACsel_dc = c(ACsel_dc[(shift+37):36],ACsel_dc[1:(shift+36)])
    }else if(shift == 0){
      ACsel_dc = ACsel_dc
    }else{
      ACsel_dc = c(ACsel_dc[(shift+1):36],ACsel_dc[1:shift])
    }
    
    
    ACsel = ACsel_dc * weight_wall + ACsel_pc * weight_place
    
    ACsel <- ACsel / sum(ACsel)
    
    
    ASrand <- runif(1) # c
    which <- 1 # c
    ASsum <- ACsel[1]# c
    
    while (which < N_ac & ASsum <= ASrand){# c
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    if (which + shift> 36){
      which_dc = which + shift - 36
    }else if(which + shift <= 0){
      which_dc = which + shift + 36
    }else{
      which_dc = which + shift
    }
    
    
    which_pc = which
    
    
    #Eligibility traces
    el_tr_pc <- el_tr_pc * etdecay # c
    el_tr_dc <- el_tr_dc * etdecay # c
    
    for (j in 1:N_ac){# c
      itmp_pc <- min(abs(j-which_pc), N_ac-abs(j-which_pc))
      itmp_dc <- min(abs(j-which_dc), N_ac-abs(j-which_dc))
      actgaus_pc <- exp(-(itmp_pc*itmp_pc)/(2*sigma_ac*sigma_ac))
      actgaus_dc <- exp(-(itmp_dc*itmp_dc)/(2*sigma_ac*sigma_ac))
      el_tr_pc[,j] <- el_tr_pc[,j] + actgaus_pc*AC_activation_pc[j]*t(t(PC_activation))
      el_tr_dc[,j] <- el_tr_dc[,j] + actgaus_dc*AC_activation_dc[j]*t(t(DC_activation))
    }
    
    

    vel_x = c(vel_x, (vel_x[length(vel_x)]+ac_const*cos(which/N_ac*2*pi))*Vdecay)
    vel_y = c(vel_y, (vel_y[length(vel_y)]+ac_const*sin(which/N_ac*2*pi))*Vdecay)

    #velocity per time step (not second)
    track_x = c(track_x, track_x[length(track_x)]+vel_x[length(vel_x)])
    track_y = c(track_y, track_y[length(track_y)]+vel_y[length(vel_y)])
    
    #Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter/2)^2)
    {
      ratio = (track_x[length(track_x)]^2 + track_y[length(track_y)]^2)/((pool_diameter/2)^2)
      track_x[length(track_x)] = track_x[length(track_x)]/sqrt(ratio)
      track_y[length(track_y)] = track_y[length(track_y)]/sqrt(ratio)
      vel_x[length(vel_x)] = track_x[length(track_x)] - track_x[length(track_x)-1]
      vel_y[length(vel_y)] = track_y[length(track_y)] - track_y[length(track_y)-1]
    }
    
    
    if (length(track_x) > 2)
    { if ((track_x[length(track_x)]  - platform_x)^2 + (track_y[length(track_y)]  - platform_y)^2 < platform_radius^2)
    { rew = 10 } #found platform - reward
      else if (track_x[length(track_x)]^2+track_y[length(track_y)]^2 > (0.99*pool_diameter/2)^2)
      { rew = -wall_pun } #hit wall - punishment
      else
      { rew = 0 } #didn't find - no reward
      
      currQ_pc = AC_activation_pc[which_pc]# c
      currQ_dc = AC_activation_dc[which_dc]# c
      tderr_pc = rew + discf*currQ_pc - prevQ_pc #temporal difference error # c
      tderr_dc = rew + discf*currQ_dc - prevQ_dc #temporal difference error # c
      weights_pc = pmax(weights_pc + lrate*tderr_pc*el_tr_pc, 0) # c
      weights_dc = pmax(weights_dc + lrate*tderr_dc*el_tr_dc, 0) # c
    }
    
    laststep = sqrt((track_x[length(track_x)]-track_x[length(track_x)-1])^2 + (track_y[length(track_y)]-track_y[length(track_y)-1])^2)
    dist = dist + laststep
    
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8*(pool_diameter/2)^2)
    { wall_zone = wall_zone + 1 }
    else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0)
    { quadrants[1] = quadrants[1] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0)
    { quadrants[2] = quadrants[2] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0)
    { quadrants[3] = quadrants[3] + 1 }
    else
    { quadrants[4] = quadrants[4] + 1 }
    
    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
    { speed_ts = mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
    latency = (length(track_x)-1) * speed_ts / speed #convert to seconds
    if (latency > 60) # if more than a minute, stop
    { break }
    }
    
  }
  
  latency <- length(track_x)-1 # latency in time steps
  wall_zone <- wall_zone/latency
  quadrants <- quadrants/latency
  speed_ts <- mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
  # speed per action step from Hanbing
  speed_ps = (vel_x[-1]^2+vel_y[-1]^2)^0.5
  
  # time step
  time_step = speed_ts/speed
  
  # mean turning angle 
  vel = as.matrix(data.frame(vel_x,vel_y))
  angle=c()
  for (steps in 2:(length(vel_x)-1)){
    A = vel[steps,]
    B = vel[steps+1,]
    angle = c(angle, acos((A%*%B)[1,1])/norm(as.matrix(A)*norm(as.matrix(B))))
  }
  angle = angle * 180 / pi
  mean_angle = mean(angle)
  
  # speed standard deviation
  speed_std = sd((vel_x[-1]^2+vel_y[-1]^2)^0.5)
  speed_std = speed_std / time_step
  
  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights_pc, weights_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency, speed_std, speed_ps,mean_angle,time_step))
}


# =========================== main =============================== 
timestart=Sys.time()

N_dc <- 10 #Population of place cells [100..300]
N_pc <- 211 #Population of place cells # c
N_ac <- 36 #Population of action cells [25..50]

plot_trajectories <- 0 #yes - 1, no - 0
plot_cognitive_maps <- 0 #yes - 1, no - 0
plot_integrated_cognitive_map <- 0
pln <- plot_trajectories + plot_cognitive_maps + plot_integrated_cognitive_map
Nruns <- 100 #how many runs to run if not plotting anything

pool_diameter <- 1.4 #Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 #Platform radius (m)
sigma_pc <- 0.1# c
sigma_dc <- c(0.1,0.1) #place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- c(2,2) #action cell sigma (standard deviation), in action cells [1..3]

etdecay <- 0.83 #Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- c(6,6) #Exploration-exploitation factor (\beta) [0.5..12]
alpha <- c(0.006,0.006) #Learning rate (\alpha) [0.005..0.02]
gamma <- c(0.85,0.85) #Discount factor (\gamma) [0.75..0.95]

Vdecay <- c(0.82,0.82) #velocity decay [0.75..0.95]
ac_const <- c(0.02,0.02) #acceleration const [0.01..0.03]
Wnoise <- c(0.0004,0.0004) #Weight noise [0.0001..0.0007]
Wmult <- c(0.1,0.1) #Weight multiplier [0.05..0.15]
hitwall <- c(0,0)#punishment for hitting the wall [0..1]
speed <- c(0.175,0.175) #mouse speed (m/s) [0.1..0.25]
weight_wall = 0.2
weight_place = 1 - weight_wall


Ntrials <- 4 #number of trials per day
Ndays <- 8 #number of days
track_x_sum=vector(mode='list',length = Ntrials*Ndays)
track_y_sum=vector(mode='list',length = Ntrials*Ndays)

Npsets = 1
params<-array(rep(0,12*Npsets),c(12,Npsets))

params[1,] =runif(Npsets)*(sigma_dc[2]-sigma_dc[1])+sigma_dc[1]
params[2,] =runif(Npsets)*(sigma_ac[2]-sigma_ac[1])+sigma_ac[1]
params[3,] =runif(Npsets)*(beta[2]-beta[1])+beta[1]
params[4,] =runif(Npsets)*(alpha[2]-alpha[1])+alpha[1]
params[5,] =runif(Npsets)*(gamma[2]-gamma[1])+gamma[1]
params[6,] =runif(Npsets)*(Vdecay[2]-Vdecay[1])+Vdecay[1]
params[7,] =runif(Npsets)*(ac_const[2]-ac_const[1])+ac_const[1]
params[8,] =runif(Npsets)*(Wnoise[2]-Wnoise[1])+Wnoise[1]
params[9,] =runif(Npsets)*(Wmult[2]-Wmult[1])+Wmult[1]
params[10,] =runif(Npsets)*(hitwall[2]-hitwall[1])+hitwall[1]
params[11,] =runif(Npsets)*(speed[2]-speed[1])+speed[1]
params[12,] =etdecay


#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
if (pln > 0.5) #if any plots
{ PMs <- array(rep(0,9*Ndays*Ntrials), c(9,Ndays,Ntrials))
AMs <- array(rep(0,Ndays*Ntrials),c(Ndays,Ntrials))
} else {
  PMs <- array(rep(0,9*Ndays*Ntrials*Nruns), c(9,Ndays,Ntrials,Nruns,Npsets))
  AMs <- array(rep(0,Ndays*Ntrials*Nruns),c(Ndays,Ntrials,Nruns))}
#multiple runs

for (pset in 1:Npsets){
  print(paste('No.',pset))
  sigma_dc = params[1,pset]
  sigma_ac = params[2,pset]
  beta = params[3,pset]
  alpha = params[4,pset]
  gamma = params[5,pset]
  Vdecay = params[6,pset]
  ac_const = params[7,pset]
  Wnoise = params[8,pset]
  Wmult = params[9,pset]
  hitwall = params[10,pset]
  speed = params[11,pset]
  etdecay = params[12,pset]
  
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
    weights_pc <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult# c
    weights_dc <- matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult# c
    
    #Generate place and distance cells
    DC <- rep(0,N_dc)
    PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell# c
    PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell# c
    
    for (i in 1:N_dc) {
      #For each distance cell:
      DC[i] <- runif(1)*(pool_diameter/2) #Random positions of distance cells
      
    }
    for (i in 1:N_pc) {# c
      #For each place cell:
      PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
      PC_y[i] <- (runif(1) - 0.5)*pool_diameter
      while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
        #Checks for out of bounds
        PC_x[i] <- (runif(1) - 0.5)*pool_diameter
        PC_y[i] <- (runif(1) - 0.5)*pool_diameter
      }
    }
    
    par(mfrow=c(4,8))
    
    
    
    for (day in 1:Ndays) {
      idxs = sample(4) #randomly choose 4 starting locations
      for (trial in 1:Ntrials){
        
        # fix or variable platform
        # whichplatform = 1 # fix
        whichplatform=sample(c(1,-1),1) # variable
        
        idx <- idxs[trial] #take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
        
        modresults <- run_trial (weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, speed, hitwall, weight_wall, weight_place)
        #run trial
        weights_pc <- modresults[[1]]# c
        weights_dc <- modresults[[2]]# c
        track_x <- modresults[[3]]
        track_y <- modresults[[4]]
        vel_x <- modresults[[5]]
        vel_y <- modresults[[6]]
        track_x_sum[[(day-1)*4+trial]]=track_x
        track_y_sum[[(day-1)*4+trial]]=track_y
        #switch_times = modresults[[15]]
        #strategy_record = modresults[[16]]
        #findtheplatform = modresults[[17]]
        #        weights <- wres
        
        PMs[1,day,trial] <- modresults[[10]] #latency
        PMs[2,day,trial] <- modresults[[7]] #dist
        if (whichplatform==1){
          PMs[3,day,trial] <- modresults[[9]][4]*100 #target quadrant
          PMs[4,day,trial] <- modresults[[9]][2]*100 #opposite quadrant
        }else{
          PMs[3,day,trial] <- modresults[[9]][2]*100 #target quadrant
          PMs[4,day,trial] <- modresults[[9]][4]*100 #opposite quadrant
        }
        PMs[5,day,trial] <- modresults[[8]]*100 #wall zone
        PMs[6,day,trial] <- modresults[[11]]*100 # speed_std
        PMs[7,day,trial] <- modresults[[13]] # mean_angle
        PMs[8,day,trial] <- modresults[[14]] # time_step
        
        #AMs[day,trial] <- modresults[[11]] # speed_ps
        #PMs[9,day,trial] <- modresults[[15]]#which model # c
        
        #record performance measures
        
        if (plot_trajectories)
        {
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l", xlab = paste("day",day,", trial",trial), ylab = "trajectory")
          #plot the trajectory
          lines(track_x, track_y, type = "l")
          
          #plot the platform
          lines(platform_x*whichplatform+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
        }
        
        if (plot_cognitive_maps){
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day",day,", trial",trial), ylab = "cognitive map")
          #plot the cognitive map
          for (x in (-6:6)*(pool_diameter/6)){
            for (y in (-6:6)*(pool_diameter/6)){
              if (x^2 + y^2 <= (pool_diameter/2)^2){
                x2 = x
                y2 = y
                dist.to.wall = pool_diameter/2-sqrt(x2^2+y2^2)
                for (k in 1:N_ac){
                  DC_activation <- rep(0,N_dc)
                  for (i in 1:N_dc){
                    DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
                  }
                  #Calculate AC activation (i.e. value of the action)
                  AC_activation <- rep(0,N_ac)
                  for (i in 1:N_ac){
                    for (j in 1:N_dc){
                      AC_activation[i] <- AC_activation[i] + DC_activation[j]*weights_dc[j,i]
                    }
                  }
                  #direction of AC activation
                  if(x>=0){
                    central.angle = atan(y/x)
                  } else {
                    central.angle = pi+atan(y/x)
                  }
                  moving.dir = pi+central.angle+k/N_ac*2*pi
                  
                  x2 <- c(x2, x + (AC_activation[k]/10)*cos(moving.dir))
                  y2 <- c(y2, y + (AC_activation[k]/10)*sin(moving.dir))
                  lines(x2,y2,type = "l",col = "red")
                  #                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
                }
                
              }
            }
          }
          # plot the platform
          lines(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
          
          # place cell based cognitive map # c
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day",day,", trial",trial), ylab = "cognitive map")
          #plot the cognitive map
          for (x in (-6:6)*(pool_diameter/6)){
            for (y in (-6:6)*(pool_diameter/6)){
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
                      AC_activation[i] <- AC_activation[i] + PC_activation[j]*weights_pc[j,i]# c
                    }
                  }
                  x2 <- c(x2, x + (AC_activation[k]/10)*cos(k/N_ac*2*pi))
                  y2 <- c(y2, y + (AC_activation[k]/10)*sin(k/N_ac*2*pi))
                  lines(x2,y2,type = "l",col = "blue")#                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
                }
                
              }
            }
          }
          # plot the platform
          lines(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
          
          
        }
        
        if (plot_integrated_cognitive_map){
          #plot the maze
          plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day",day,", trial",trial), ylab = "cognitive map")
          #plot the cognitive map
          for (x in (-6:6)*(pool_diameter/6)){
            for (y in (-6:6)*(pool_diameter/6)){
              if (x^2 + y^2 <= (pool_diameter/2)^2){
                x2 = x
                y2 = y
                dist.to.wall = pool_diameter/2-sqrt(x2^2+y2^2)
                for (k in 1:N_ac){
                  DC_activation <- rep(0,N_dc)
                  PC_activation <- rep(0,N_pc)
                  for (i in 1:N_dc){
                    DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
                  }
                  for (i in 1:N_pc){
                    PC_activation[i] <- exp(-((x - PC_x[i])^2 + (y - PC_y[i])^2)/(2*sigma_pc^2))
                  }
                  #Calculate AC activation (i.e. value of the action)
                  AC_activation_dc <- rep(0,N_ac)
                  AC_activation_pc <- rep(0,N_ac)
                  for (i in 1:N_ac){
                    for (j in 1:N_dc){
                      AC_activation_dc[i] <- AC_activation_dc[i] + DC_activation[j]*weights_dc[j,i]
                    }
                    for (j in 1:N_pc){
                      AC_activation_pc[i] <- AC_activation_pc[i] + PC_activation[j]*weights_pc[j,i]# c
                    }
                  }
                  shift = round((180-atan2(y,x)/pi*180)/10)
                  if(shift == 36){
                    shift = 0
                  }
                  
                  
                  if (shift < 0){
                    AC_activation_dc = c(AC_activation_dc[(shift+37):36],AC_activation_dc[1:(shift+36)])
                  }else if(shift == 0){
                    AC_activation_dc = AC_activation_dc
                  }else{
                    AC_activation_dc = c(AC_activation_dc[(shift+1):36],AC_activation_dc[1:shift])
                  }
                  
                  AC_activation_dc = AC_activation_dc/sum(AC_activation_dc)*15
                  AC_activation_pc = AC_activation_pc/sum(AC_activation_pc)*15
                  
                  AC_activation = AC_activation_dc * weight_wall + AC_activation_pc * weight_place
                  
                  x2 <- c(x2, x + (AC_activation[k]/10)*cos(k/N_ac*2*pi))
                  y2 <- c(y2, y + (AC_activation[k]/10)*sin(k/N_ac*2*pi))
                  lines(x2,y2,type = "l",col = "blue")
                  
                }
                
              }
            }
          }
          # plot the platform
          lines(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
          
        }
        
      }
    }
  } else {
    # run multiple times without plotting!
    print('else')
    for (reps in 1:Nruns){
      print(reps)
      
      #Generate initial weights
      weights_pc <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult# c
      weights_dc <- matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult# c
      
      #Generate place and distance cells
      DC <- rep(0,N_dc)
      PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell# c
      PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell# c
      
      for (i in 1:N_dc) {
        #For each place cell:
        DC[i] <- runif(1)*(pool_diameter) #Random positions of distance cells
        
      }
      for (i in 1:N_pc) {# c
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
          
          # fix or variable platform
          # whichplatform = 1 # fix
          whichplatform=sample(c(1,-1),1) # variable
          
          idx <- idxs[trial] #take each location
          starting_x <- starting_xs[idx]
          starting_y <- starting_ys[idx]
          
          modresults <- run_trial (weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, speed, hitwall,weight_wall, weight_place)
          #run trial
          weights_pc <- modresults[[1]]# c
          weights_dc <- modresults[[2]]# c
          track_x <- modresults[[3]]
          track_y <- modresults[[4]]
          vel_x <- modresults[[5]]
          vel_y <- modresults[[6]]
          track_x_sum[[(day-1)*4+trial]]=track_x
          track_y_sum[[(day-1)*4+trial]]=track_y
          #        weights <- wres
          
          PMs[1,day,trial,reps,pset] <- modresults[[10]] #latency
          PMs[2,day,trial,reps,pset] <- modresults[[7]] #dist
          if (whichplatform==1){
            PMs[3,day,trial,reps,pset] <- modresults[[9]][4]*100 #target quadrant
            PMs[4,day,trial,reps,pset] <- modresults[[9]][2]*100 #opposite quadrant
          }else{
            PMs[3,day,trial,reps,pset] <- modresults[[9]][2]*100 #target quadrant
            PMs[4,day,trial,reps,pset] <- modresults[[9]][4]*100 #opposite quadrant
          }
          PMs[5,day,trial,reps,pset] <- modresults[[8]]*100 #wall zone
          PMs[6,day,trial,reps,pset] <- modresults[[11]]*100 # speed_std
          PMs[7,day,trial,reps,pset] <- modresults[[13]] # mean_angle
          PMs[8,day,trial,reps,pset] <- modresults[[14]] # time_step
        }}
    }
  }
}

#########################output file##############################
if (reps == Nruns){
  latency = PMs[1,,,,]
  dist = PMs[2,,,,]
  target_quadrant = PMs[3,,,,]
  opposite_quadrant = PMs[4,,,,]
  wall_zone = PMs[5,,,,]
  day = c(1:Ndays)
  latency = cbind(day,latency)
  dist = cbind(day,dist)
  target_quadrant = cbind(day,target_quadrant)
  opposite_quadrant = cbind(day,opposite_quadrant)
  wall_zone = cbind(day,wall_zone)
  all = cbind(latency,dist[,2],target_quadrant[,2],opposite_quadrant[,2],wall_zone[,2])
  colnames(all)=c("day","latency","dist","target_quadrant","opposite_quadrant","wall_zone")
  write.csv(all,"variable_wall_p-0.csv")
}

# p1 = read.csv("fixed_start-wall_p-0.004_1.csv")[,2:7]
# p2 = read.csv("fixed_start-wall_p-0.004_2.csv")[,2:7]
# p3 = read.csv("fixed_start-wall_p-0.004_3.csv")[,2:7]
# p4 = read.csv("fixed_start-wall_p-0.004_4.csv")[,2:7]
# p5 = read.csv("fixed_start-wall_p-0.004_5.csv")[,2:7]
# p6 = read.csv("fixed_start-wall_p-0.004_6.csv")[,2:7]
# p7 = read.csv("fixed_start-wall_p-0.004_7.csv")[,2:7]
# p8 = read.csv("fixed_start-wall_p-0.004_8.csv")[,2:7]
# p9 = read.csv("fixed_start-wall_p-0.004_9.csv")[,2:7]
# p10 = read.csv("fixed_start-wall_p-0.004_10.csv")[,2:7]
# 
# all = rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
# 
# write.csv(all,"fixed_start-wall_p-0.004.csv")

