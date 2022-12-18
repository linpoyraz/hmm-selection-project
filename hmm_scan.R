# SLIM VCF editor 
library(data.table)
library(parallel)
library(pbmcapply)

# times for which we have observation 
times = seq(0,200, by = 10)
times[1] = 1

# function to read and clean up SLiM outputs 
read_slim = function(path){
  muts = read.table(path, sep = " ", fill = T)
  # strip until genomes
  until = which(muts[,1] == "Genomes:")
  muts = muts[2:(until-1),]
  muts = muts[,c(2,3,4,9)]
  names(muts) = c("id", "type", "pos", "N.A")
  return(muts)
}

# read t1 to get a list of all mutations we will be tracking 
t = read_slim("./slim_s2/t_1.txt")
names(t) = c("id", "type", "pos", "1")

# read all other times, merge on initial mutations (newly arising mutations excluded)
for (i in times[2:length(times)]){
  f = paste("./slim_s2/t_", i, ".txt", sep = "")
  muts = read_slim(f)
  names(muts)[4] = i
  # merge with t
  t = merge(x = t, y = muts, by = c("id", "type", "pos"), all.x = T)
  # replace NA with 0 
  is_na = which(is.na(t[,ncol(t)]))
  t[is_na,ncol(t)] = 0
}

# turn to list of obs data frames to prepare for mclapply 
t_list = list()

for (i in 1:nrow(t)){
  N = rep(0, 200)
  N[times] = 100
  N_A = rep(0,200)
  N_A[times] = unlist(t[i,4:ncol(t)])
  obs_t = list(data.frame("id" = rep(t$id[i], 200), "pos" = rep(t$pos[i],200), "N" = N, "N.A" = N_A))
  t_list = c(t_list, obs_t)
}

################################################################################
# THESE ARE THE FUNCTIONS FOR A SINGLE ALLELE WHICH WE WILL RUN IN PARALLEL
### PARAMETERS #################################################################

set.seed(444)
Ne <- 10000                   # effective population size
gens <- 200                  # number of generations
K <- 2                      # number of hidden selection states 
D <- 100                     # grid size of the discrete allele freq
sigma_init <- c(0, 0.03) # initial guess for selection coefficient
tolerance = 0.001

### SIMULATION #################################################################

#obs = data.frame("N" = sim$n, "N.A" = sim$a)

### HIDDEN MARKOV MODEL DEFITINION #############################################

hmm_constructor <- function(obs, D, Ne, s){
  
  ## initialize discrete frequency space
  obs_freq <- obs$N.A / obs$N
  
  f.max <- min(1, max(obs_freq + 0.1, na.rm=TRUE))
  f.min <- max(0, min(obs_freq - 0.1, na.rm=TRUE))
  
  freq <- seq(f.min, f.max, length.out = D)
  delta <- (f.max - f.min) / (D - 1)
  
  
  # define the transition probability function for all z_t
  # same indicates which index is has prob 1-u
  tf <- function(f_from, f_to, same){
    u <- 1 / NROW(obs)
    mu <- f_from + 2 * rep(s, each=D) * f_from * (1 - f_from) * (f_from + (1 - 2 * f_from) / 2)
    nu <- sqrt(f_from * (1 - f_from) / Ne)
    high <- f_to + delta/2
    low <- f_to - delta/2
    probs <- pnorm(high, mu, nu) - pnorm(low, mu, nu)
    all <- probs * ifelse(same == 1, 1-u, u)
    return(all)
  }
  
  ## define the emission probability function
  ef <- function(a, n, f){
    return(dbinom(a, n, f))
  }
  
  ## return the hmm object
  hmm <- list(obs = obs,
              total_t = NROW(obs),
              f_states = freq,
              transition_func = tf,
              emission_func = ef)
  return(hmm)
}

### HIDDEN MARKOV MODEL ALGORITHMS #############################################

forward_backward <- function(hmm){
  
  ### Compute forward probabilities
  
  F <- matrix(nrow = D*K, ncol = hmm$total_t) # rows are the hidden states, columns are t
  scaling <- rep(0,hmm$total_t)        # we use scaling factors since even log probs don't have enough precision
  
  # initialize F and scale
  F[,1] <- 1/(D*K) * rep(hmm$emission_func(hmm$obs$N.A[1], hmm$obs$N[1], hmm$f_states), times = K)
  scaling[1] <- sum(F[,1])
  F[,1] <- F[,1] / scaling[1]
  
  for(t in 2:hmm$total_t) {
    for(i in 1:(D*K)) {
      
      # mark same index region
      same <- rep(0,D*K)
      start <- D * ((i-1) %/% D)
      for (j in (start + 1):(start + D)) {
        same[j] <- 1
      }
      
      # forward algorithm recurrence 
      f_from = rep(hmm$f_states, times=K)
      f_to = hmm$f_states[(i-1) %% D + 1]
      state_probs <- F[,t-1] * hmm$transition_func(f_from, f_to, same)
      F[i,t] <- hmm$emission_func(hmm$obs$N.A[t], hmm$obs$N[t], hmm$f_states[(i-1) %% D + 1]) * sum(state_probs)    
    }
    
    # scale 
    scaling[t] <- sum(F[,t])
    F[,t] <-   F[,t] / scaling[t]
  }
  
  ### Compute backward probabilities (using the scaling factor from forward)
  
  # initialize B
  B <- matrix(nrow = D*K, ncol = hmm$total_t)
  B[,hmm$total_t] <- 1
  
  for(t in rev(1:(hmm$total_t-1))){
    for(i in 1:(D*K)){
      
      # mark same index region
      same <- rep(0,D*K)
      start <- D * ((i-1) %/% D)
      for (j in (start + 1):(start + D)) {
        same[j] <- 1
      }
      
      # backward algorithm recurrence 
      f_from = hmm$f_states[(i-1) %% D + 1]
      f_to = rep(hmm$f_states, times = K)
      state_probs <- B[,t+1] *
        hmm$transition_func(f_from, f_to, same) *
        rep(hmm$emission_func(hmm$obs$N.A[t+1], hmm$obs$N[t+1], hmm$f_states), times = K)
      B[i,t] <- sum(state_probs) / scaling[t+1]
    }
  }
  
  P = sum(F[,hmm$total_t])
  FB <- F * B / P
  return(list(F=F, B=B, FB=FB, scaling=scaling, P=P))
}

viterbi <- function(hmm){
  
  # traceback matrix
  B <- matrix(nrow=D*K, ncol=hmm$total_t)
  
  # initialization
  M <- matrix(nrow=D*K, ncol=hmm$total_t)
  M[,1] <- 1/(D*K) * rep(hmm$emission_func(hmm$obs$N.A[1], hmm$obs$N[1], hmm$f_states), times = K)
  
  for(t in 2:hmm$total_t){
    for(i in 1:(D*K)){
      
      # mark same index region   
      same <- rep(0,D*K)
      start <- D * ((i-1) %/% D)
      for (j in (start + 1):(start + D)) {
        same[j] <- 1
      }
      
      state_probs <- M[,t-1] * hmm$transition_func(rep(hmm$f_states, times=K), hmm$f_states[(i-1) %% D + 1], same)
      max_idx <- which.max(state_probs)
      
      M[i,t] <- state_probs[max_idx] * hmm$emission_func(hmm$obs$N.A[t], hmm$obs$N[t], hmm$f_states[(i-1) %% D + 1])
      B[i,t] <- max_idx
    }
  }
  
  # traceback viterbi path
  trace <- c(which.max(M[,hmm$total_t]))
  cur <- hmm$total_t - 1
  while (cur > 0) {
    trace <- c(trace, B[trace[hmm$total_t - cur], cur+1])
    cur <- cur - 1
  }
  return(rev(trace))
}

hmm_estimator <- function(obs, s){
  
  # construct hmm
  hmm <- hmm_constructor(obs, D, Ne, s)
  
  # compute viterbi path
  #path <- viterbi(hmm)
  #f_path <- rep(hmm$f_states, times=length(s))[path]
  #s_path <- rep(1:length(s), each=length(hmm$f_states))[path]
  
  # compute forward-backward probabilities
  fb_results <- forward_backward(hmm)
  
  f_FB <- matrix(0, nrow=D, ncol=NCOL(fb_results$FB))
  s_FB <- matrix(0, nrow=K, ncol=NCOL(fb_results$FB))
  for(i in 1:K){
    start <- D*(i-1) + 1
    end <- D*i
    f_FB <- f_FB + fb_results$FB[start:end,]
    s_FB[i,] <- colSums(fb_results$FB[start:end,])
  }
  
  # return results as an object
  results <- list(
    hmm = hmm,
    #f_viterbi = f_path,
    #s_viterbi = s_path,
    F = fb_results$F,
    B = fb_results$B,
    FB = fb_results$FB,
    f_FB = f_FB,
    s_FB = s_FB,
    log_likelihood = sum(log(fb_results$scaling))
  )
  
  return(results)
}

### EXPECTATION-MAXIMIZATION ALGORITHM (APPROXIMATE) ###########################

# E[f_{t+1} | f_t,z_t = state_i]
cond_exp <- function(hmm_est, i, t){
  
  hmm <- hmm_est$hmm
  
  # mark same index region
  same <- rep(0,D*K)
  start <- D * ((i-1) %/% D)
  for (j in (start + 1):(start + D)) { same[j] <- 1 }
  
  # initialize
  probs <- hmm_est$F[i,t] *
    hmm$transition_func(
      hmm$f_states[(i-1) %% D + 1], 
      rep(hmm$f_states, times=K), same
    ) *
    hmm_est$B[,t+1] 
  
  probs <- probs / hmm_est$FB[i,t]
  
  # multiply by emmision probs  
  probs <- probs * 
    rep(hmm$emission_func(hmm$obs$N.A[t+1], 
                          hmm$obs$N[t+1], hmm$f_states), 
        times=K
    )
  
  # normalize
  probs <- probs / sum(probs)             
  
  probs <- probs * rep(hmm$f_states, times=K)
  return(sum(probs))
}


em_rule_update <- function(hmm_est){
  hmm <- hmm_est$hmm
  
  # update each s_state individually
  new_s <- rep(0, K)
  for(i in 1:K) {
    FB_segment <- hmm_est$FB[(D*(i-1) + 1):(D*i),]
    exp <- rep(0, hmm$total_t-1)
    for(t in 1:(hmm$total_t-1)){
      sum <- 0
      for(f in 1:D){
        f1 <- FB_segment[f,t]
        f2 <- cond_exp(hmm_est, D*(i-1)+f, t)
        
        if(!is.nan(f1*f2)){
          sum <- sum + (f1 * f2)
        }
      }
      exp[t] <- sum
    }
    
    mean_f <- colSums(apply(FB_segment, 2, "*", hmm$f_states))[1:(hmm$total_t-1)]
    num <- sum(exp - mean_f)
    den <- sum(colSums(apply(FB_segment[,1:(hmm$total_t-1)], 2, "*", hmm$f_states * (1-hmm$f_states))))
    new_s[i] <- num / den  
  }
  
  return(new_s)
}


### INFERENCE ##################################################################

estimate_s <-  function(obs, init_s, tol = tolerance, max_iter = 100){
  
  # use EM-update rule until convergence  
  s <- init_s
  cur_s <- s + 2 * tol
  
  i = 1
  while(tol < max(abs(s - cur_s)) & i < max_iter) {
    #print(s)
    cur_s <- s
    hmm_est <- hmm_estimator(obs, s)
    s <- em_rule_update(hmm_est)
    i <- i + 1
  }
  
  return(list(s=s, hmm_est=hmm_estimator(obs, s)))
}

infer = function(obs){
  inferred = estimate_s(obs, sigma_init)
  s_est <- apply(inferred$hmm_est$s_FB, 2, which.max)
  return(list("s" = inferred$s, "s_est" = s_est))
}

#infer(t_list[[162]])
print("Start parallel run")
system.time({res = pbmclapply(t_list, infer, mc.cores = 4)})

# create data frame to plot for t = 25, t = 75, t = 125, t = 175 
generation = c(25,75,125,175)
for (i in 1:4){
  to_plot = data.frame("type" = t$type, "pos"= t$pos, "s" = rep(NA, length(t$pos)))
  for (j in 1:length(res)){
    if (sum(t_list[[j]]$N.A == 0)/200 < 0.97){
      s_idx = res[[j]]$s_est[generation[i]]
      s_coef = res[[j]]$s[s_idx]
      to_plot$s[j] = s_coef
    }
  }
  main_lab = paste("Selection Landscape at Generation", generation[i])
  plot = ggplot(data = to_plot, aes(x = pos, y = s, color = type)) + geom_point() + 
    theme_bw() + labs(x = "Position", y = "Inferred s", color = "", title = main_lab)
  ggsave(filename = paste("./results/landscape_", generation[i], ".png", sep = ""))
}

# 3D plot 
x = rep(t$pos, each = 200)
y = rep(1:200, length(t$pos))
z = rep(NA, 200*length(t$pos)) 

for (j in 1:length(res)){
  if (sum(t_list[[j]]$N.A == 0)/200 < 0.97){
    s_idx = res[[j]]$s_est
    s_coef = res[[j]]$s[s_idx]
    z[((i-1)*200+1):(i*200)] = s_coef
  }
}
  
to_plot = data.frame(x = x, y = y, z = z)
# z matrix 

p = plot_ly(to_plot, x = ~x, y = ~y, z = ~z) %>% add_markers(size = 1, color = ~z) %>% 
    layout(scene = list(xaxis = list(title = "position"),
                        yaxis = list(title = "time"),
                        zaxis = list(title = "selection coefficient")))

export(p = p, file = "./results/3d.png")


