competing_risks_algo <- function(
    r_prob=0.5, r_dist, d_prob=0.5, d_dist,
    beta, S_init, I_init, R_init, D_init, N=1000,
    dt = 0.05, duration = 75,
    return_dt = 0.5 # dt in returned data.frame
){
  # ----- Compute out transition here --------
  # assume dwell time of r >= dwell time of d
  # i_dist is the new dwell time distribution in I
  r_maxtime <- length(r_dist)
  d_maxtime <- length(d_dist)
  
  # prob of people that would end up in R and stay in I up to d_maxtime
  # equivalent of z'
  p_prime <- r_prob*sum(r_dist[1:d_maxtime])/sum(r_dist)
  
  # distribution normalizer (denominator for new dwell time computation)
  normalizer <- sapply(1:d_maxtime, \(i){
    p_prime*r_dist[i] + (1-r_prob)*d_dist[i]
  })
  
  i_dist <- rep(0, r_maxtime)
  # Bug due to sum of i_dist is greater than 1
  i_dist[1:d_maxtime] <- sapply(1:d_maxtime, \(i){
    (p_prime*r_dist[i] + d_prob*d_dist[i])/sum(normalizer)
  })
  
  # for i_dist at timestep > d_maxtime is just the probability distribution of I->R transition
  if(r_maxtime > d_maxtime){
    i_dist[(d_maxtime+1):r_maxtime] <- sapply((d_maxtime+1):r_maxtime, \(i){
      r_prob*r_dist[i]
    })
  }
  
  # transition prob
  transprob <- rep(0, r_maxtime)
  transprob[1] <- i_dist[1]
  transprob[2:r_maxtime] <- sapply(2:r_maxtime, \(i){
    i_dist[i]/(1 - sum(i_dist[1:(i-1)]))
  })
  
  # ----- Model update here ------
  no_timesteps <- duration/dt + 1
  S <- c(S_init, rep(0, no_timesteps-1))
  R <- c(R_init, rep(0, no_timesteps-1))
  D <- c(D_init, rep(0, no_timesteps-1))
  I <- matrix(0, nrow=r_maxtime, ncol=no_timesteps)
  I[,1] <- I_init
  
  for (i in 2:no_timesteps){
    S[i] <- S[i-1] - dt * S[i-1] * beta * sum(I[,i-1]) / N
    I[1,i] <- dt * S[i-1] * beta* sum(I[,i-1])/N
    
    I[2:r_maxtime, i] <- sapply(2:r_maxtime,\(subcomp){
      I[subcomp-1, i-1]*(1-transprob[subcomp-1])
    })
    
    I_to_R <- sapply(1:r_maxtime, \(subcomp){
      transprob[subcomp]*I[subcomp, i-1]*
        (p_prime*r_dist[subcomp])/
        (p_prime*r_dist[subcomp] + (1-r_prob)*d_dist[subcomp])
    })
    R[i] <- R[i-1] + sum(I_to_R)
    
    I_to_D <- sapply(1:d_maxtime, \(subcomp){
      transprob[subcomp]*I[subcomp, i-1]*
        (d_prob*d_dist[subcomp])/
        (p_prime*r_dist[subcomp] + (1-r_prob)*d_dist[subcomp])
    })
    D[i] <- D[i-1] + sum(I_to_D)
  }
  
  data.frame(
    time = seq(0, duration, dt),
    S = S, 
    I = colSums(I),
    R = R, 
    D = D
  ) %>% 
    filter(
      time %in% seq(0, duration, return_dt)
    )
}

competing_risks_algo(
  r_dist = dist_r, d_dist= dist_d,
  beta = mod2_config$beta,
  S_init = mod2_config$S,
  I_init = c(mod2_config$I, rep(0, length(dist_r)-1)),
  R_init = mod2_config$R,
  D_init = mod2_config$D,
  dt = timestep,
  duration = 75,
  return_dt = 0.2
) 