# re-implementation of the paper
# usage: hmm_single_allele.R <simulation file> <outname>
# output: results plot of estimated and observed frequency values 

library(ggplot2)
library(plotrix)
library(data.table)
### PARAMETERS #################################################################

args = commandArgs(trailingOnly = T)
sim = read.table(args[1], sep = "\t", header = T)
Ne <- 10000                   # effective population size
gens <- 200                  # number of generations
K <- 3                       # number of hidden selection states 
D <- 100                     # grid size of the discrete allele freq
sigma_init <- c(-0.03, 0, 0.03) # initial guess for selection coefficient
tolerance = 0.001
  
### SIMULATION #################################################################

obs = data.frame("N" = sim$n, "N.A" = sim$a)

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
  path <- viterbi(hmm)
  f_path <- rep(hmm$f_states, times=length(s))[path]
  s_path <- rep(1:length(s), each=length(hmm$f_states))[path]
  
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
    f_viterbi = f_path,
    s_viterbi = s_path,
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
    print(s)
    cur_s <- s
    hmm_est <- hmm_estimator(obs, s)
    s <- em_rule_update(hmm_est)
    i <- i + 1
  }
  
  return(list(s=s, hmm_est=hmm_estimator(obs, s)))
}


### RUN INFERENCE ##############################################################

inferred<-estimate_s(obs, sigma_init)
print(inferred$s)
f_est <- inferred$hmm_est$hmm$f_states[apply(inferred$hmm_est$f_FB, 2, which.max)]
s_est <- apply(inferred$hmm_est$s_FB, 2, which.max)
print(data.frame(f = f_est, s = s_est))

## PLOT ##############################################################

f_plot = paste("./results/single_allele/", args[2], ".png", sep = "")
png(f_plot, units = "px", height = 1500, width = 1500, res = 200) # rename 
plot.new()
# for k = 3, currently hard-coded 
cols = c("blue", "green", "red")
g = gens
plot(obs$N.A/obs$N, pch=16, bty="n", ylab="Frequency", xlab="Generation", ylim=c(0, 0.5))

# segments and colors for true underlying
f = sim$freq
f_s = unique(sim$s)
f_col = rep(cols[2], g)
f_col[which(sim$s == min(f_s))] = cols[1]
f_col[which(sim$s == max(f_s))] = cols[3]
segments(1:(g-1), f[1:(g-1)], 2:g, f[2:g], col=f_col, lwd=2)

# segments and colors for inferred underlying using posterior mode 
segments(1:(g-1), f_est[1:(g-1)], 2:g, f_est[2:g], col= cols[s_est], lwd=1.5, lty=3)

legend("topleft", c("True", "Inferred", "Observed", sort(unique(sim$s))), lwd=c(2,2,NA,4,4,4), lty=c(1,3,NA,1,1,1), pch=c(NA, NA, 16,NA,NA,NA), col=c(1, 1, 1, cols), bty="n")

dev.off()



