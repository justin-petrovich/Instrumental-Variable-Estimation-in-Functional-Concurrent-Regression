# scheme: an integer: 1, 2, or 3, corresponding to the desired sampling scheme
# K: the length of the grid on which all functions are observed
# N: the sample size (number of functions)
# sparsity: a vector of integers indicating the range of possible values for the
# sparsity of each observed function. If a single value, all functions will be 
# observed that number of times. The values in sparisty must be no bigger than K
# data_full: the data frame of all fully observed functions, in long form

library(dplyr)

sampscheme <- function(scheme = 1,sparsity,data_full){
  
  N <- length(unique(data_full$subj))
  K <- length(unique(data_full$argvals))
  
  if(length(sparsity)==1){
    mi <- rep(sparsity,N)
  }else{
    mi <- sample(sparsity,size = N,replace = T)
  }
  
  #-------------------------------------------------------------------#
  #--- Scheme 1: Random sampling of grid points from each function ---#
  #-------------------------------------------------------------------#
  if(scheme==1){
    row_obs <- rep(NA,sum(mi))
    row_obs[1:mi[1]] <- sort(sample(1:K,size = mi[1],replace = F))
    for(i in 2:N){
      ri <- sort(sample(1:K,size = mi[i],replace = F))
      row_obs[(cumsum(mi)[i-1] + 1):(cumsum(mi)[i-1] + mi[i])] <- ri + K*(i-1)
    }
    
    # Ensure that the boundaries of the grid are observed
    row_obs[1] <- 1
    row_obs[sum(mi)] <- N*K
    
    # observed rows of full data frame
    data_obs <- data_full[row_obs,]
  }
  #-----------------------------------------------------------#
  #--- Scheme 2: Consecutive grid points for each function ---#
  #-----------------------------------------------------------#
  if(scheme==2){
    if(length(sparsity) > 1){
      stop("For scheme 2, sparsity should be an integer of length 1. All functions are
           expected to be sampled the same number of times.")
    }
    
    # Grid point of first observation for each function
    start <- sample(1:K,N-2,replace = T)
    
    # Ensure that the boundaries of the grid are observed
    start <- c(1,start)
    start <- c(start,K - mi[N] + 1)
    
    # mi need adjusted downward if first observation occurs too late in the grid
    if(any(start + mi - 1 > K)){
      ids <- which(start + mi - 1 > K)
      notobs <- (start + mi - 1 - K)[ids]
      mi[ids] <- mi[ids] - notobs
    }
    
    row_obs <- rep(NA,sum(mi))
    
    for(i in 1:N){
      Ki <- start[i]:(start[i] + mi[i] - 1)
      ri <- Ki + K*(i-1)
      row_obs[(cumsum(mi)[i] - mi[i] + 1):(cumsum(mi)[i])] <- ri
    }
    
    # observed rows of full data frame
    data_obs <- data_full[row_obs,]
  }
  #---------------------------------------------------------------#
  #--- Scheme 3:                                               ---#
  #---------------------------------------------------------------#
  
  
  

  return(data_obs)
}