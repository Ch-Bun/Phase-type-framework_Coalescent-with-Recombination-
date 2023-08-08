# load librarys ----
library(PhaseTypeR)
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
# plot graph
library(igraph)
library(ggplot2)
library(gtools)
library(partitions)

# Function: Generate state space ----

  # generate all possible unique states of Rate-Matrix and returns Rate-Matrix filled with 0s
  # all possible unique states: (n_ab,n_a,n_b), with  n_ab + max(n_a,n_b) <= n and n_ab + min(n_a,n_b) >= 1
  
  # For state space reduction all states with same n_a+n_ab,n_b+n_ab are saved in:
  # same_K: key value pair -> key = n_a+n_ab,n_b+n_ab, value: list of all states with this key 

stateSpace <- function(n){
  same_K <- list()
  cnt <- 0
  lst <- list()
  for (i in rev(0:n)){
    for(j in rev(0:n)){
      if((i+j > n) | (i+j < 1)){
        next
      }
      for(k in rev(0:n)){
        if((i+k > n) | (i+k < 1) ){
          next
        }
        else{
          cnt <- cnt+1
          lst[[length(lst)+1]] <- paste(c(i,j,k),collapse=",")
          key = toString(c(i+j,i+k))
          same_K[[key]] <- c(same_K[[key]], paste(c(i,j,k),collapse=","))
        }
      }
    }
  }
  return(list(StateSpace = array(0, dim = c(cnt, cnt),dimnames = list(lst, lst)), List_of_same_K = same_K))
}

# Function: Calculate intensity matrix and reward vectors for T_MRCA and T_Total ----

  # Idea: For every state max. 7 possible paths:
    # 6 Coalescense:  - 3 times seq of same type coalesces
    #                 - 3 times seq of diff. type coalesce
    # Recombination:  - if n_ab > 0 -> recombine with rate n_ab*rho

IntensityMatrix <- function(n, rho){
  # Initialize
  absorbing_state <- c("1,0,0", "0,1,1")
  stateSpace_listK <- stateSpace(n)
  stateSpace <- stateSpace_listK[["StateSpace"]]
  list_k <- stateSpace_listK[["List_of_same_K"]]
  # Remove absorbing states from list of same K's
  list_k$`1, 1`<- NULL
  
  # Initialize reward vectors 
  rho <- rho/2
  reward_1_total <- vector()
  reward_2_total <- vector()
  reward_1_MRCA <- vector()
  reward_2_MRCA <- vector()
  
  # For every row calculate all rates to leave the state
  for (row in rownames(stateSpace)){
    # If absorbing state (p+1), not possible to leave -> skip
    if(row %in% absorbing_state){
      next
    }
    
    # i = [1], j = [2], k = [3]
    row_state <- as.integer(strsplit(row, ",")[[1]])
    
    # Generate reward vector = #lineages in state i + (j or k)
    # T_total
    reward_1_total <- append(reward_1_total, (if((row_state[1]+row_state[2]) > 1)(row_state[1]+row_state[2]) else 0) )
    reward_2_total <- append(reward_2_total, (if ((row_state[1]+row_state[3]) > 1)(row_state[1]+row_state[3]) else 0) )
    # T_MRCA
    if ((row_state[1]+row_state[2]) > 1){
      reward_1_MRCA <- append(reward_1_MRCA, 1)
    }else{
      reward_1_MRCA <- append(reward_1_MRCA, 0)
    }
    if ((row_state[1]+row_state[3]) > 1){
      reward_2_MRCA <- append(reward_2_MRCA, 1)
    }else{
      reward_2_MRCA <- append(reward_2_MRCA, 0)
    }
    
    # Coalescence of two lineages of the same type
    if (row_state[1] > 1){
      col <- paste(c(row_state[1]-1, row_state[2], row_state[3]),collapse=",")
      stateSpace[row,col] <- choose(row_state[1], 2)
    }
    if (row_state[2] > 1){
      col <- paste(c(row_state[1], row_state[2]-1, row_state[3]),collapse=",")
      stateSpace[row,col] <- choose(row_state[2], 2)
    }
    if (row_state[3] > 1){
      col <- paste(c(row_state[1], row_state[2], row_state[3]-1),collapse=",")
      stateSpace[row,col] <- choose(row_state[3], 2)
    }
    
    if (row_state[1] > 0 ){
      # Recombination
      if (rho > 0){
        col <- paste(c(row_state[1]-1, row_state[2]+1, row_state[3]+1),collapse=",")
        stateSpace[row,col] <- row_state[1]*rho
      }
      # Coalescence of two lineages of the different type
      # i and j
      if(row_state[2] > 0){
        col <- paste(c(row_state[1], row_state[2]-1, row_state[3]),collapse=",")
        stateSpace[row,col] <- (row_state[1]*row_state[2])
      }
      # i and k
      if(row_state[3] > 0){
        col <- paste(c(row_state[1], row_state[2], row_state[3]-1),collapse=",")
        stateSpace[row,col] <- (row_state[1]*row_state[3])
      }
    }
    # j and k
    if ((row_state[2] > 0) & (row_state[3] > 0)){
      col <- paste(c(row_state[1]+1, row_state[2]-1, row_state[3]-1),collapse=",")
      stateSpace[row,col] <- (row_state[2]*row_state[3])
    }
  }
  
  # Add -rowSum on diagonal to have rowSum = 0
  diag(stateSpace) <- -(rowSums(stateSpace))
  
  # Remove absorbing state to get sub-intensity matrix S
  stateSpace <- stateSpace[!rownames(stateSpace) %in% absorbing_state, !colnames(stateSpace) %in% absorbing_state]
  
  
  # Initial vector alpha
  # -> Assuming that we always start from having n sequences of type i (both loci ancestral)
  alpha <- c(1, rep(0, length(reward_2_total)-1))
  
  # Return intensity matrix and rewards  
  return(list(subIntMat = stateSpace, initProbs = alpha, rewardMatTotal = cbind(reward_1_total, reward_2_total), rewardMatMRCA = cbind(reward_1_MRCA, reward_2_MRCA), List_of_same_K=list_k))
}


# Reduce state space function ----
  # requires object of return-type from IntesityMatrix()
  # Idea: If rows are the same and have same exit rate we can join the states

reduce_S <- function(obj){
  # load object
  stateSpace <- obj$subIntMat
  alpha <- obj$initProbs
  rewardMatTotal <- obj$rewardMatTotal
  rewardMatMRCA <- obj$rewardMatMRCA
  List_of_same_K <- obj$List_of_same_K
  reduced <- TRUE
  
  # To keep track on which elements to remove from reward matrix and alpha
  drop_list <- c()
  absorbing_state <- c("1,0,0", "0,1,1")
  original_stateSpace <- stateSpace
  # Keep trying to reduce, as long as something got reduced in the last iteration
  while(reduced == TRUE){
    reduced = FALSE
    # For every group with same n_ab+n_a, n_ab+n_b
    for(k_list in names(List_of_same_K)){
      # For every element/state of that group
      for(row in List_of_same_K[[k_list]]){
        # Compare element to all elements before in that list (To avoid double checking)
        for(row_against in List_of_same_K[[k_list]][0: max(0, match(row, List_of_same_K[[k_list]])-1)]){
          # Absorbing state 0,1,1 and 1,0,0 need to be treated extra as sum
          # If rows are the same and have same exit rate we can join the states
          if(all(stateSpace[row_against,!colnames(stateSpace) %in% c(row, row_against,absorbing_state)] == 
                 stateSpace[row,!colnames(stateSpace) %in% c(row, row_against,absorbing_state)]) &
             (sum(stateSpace[row_against,colnames(stateSpace) %in% absorbing_state]) ==
              sum(stateSpace[row,colnames(stateSpace) %in% absorbing_state]))) {
            # Add pos to drop_list, to remove later from rewards and alpha
            drop_list <- append(drop_list, which(rownames(original_stateSpace) == row_against)) 
            # Sum columns 
            stateSpace[,row] <- stateSpace[, row_against] + stateSpace[, row]
            # Remove the other column and row after joining
            stateSpace <- stateSpace[!rownames(stateSpace) %in% row_against,!colnames(stateSpace) %in% row_against]
            # Remove joint element from K list
            List_of_same_K[[k_list]] <- List_of_same_K[[k_list]][List_of_same_K[[k_list]] != row]
            reduced = TRUE
          }
        }
      }
    }
  }
  return(list(subIntMat = stateSpace, initProbs = alpha[-drop_list], rewardMatTotal = rewardMatTotal[-drop_list,], rewardMatMRCA = rewardMatMRCA[-drop_list,], List_of_same_K=List_of_same_K))
}

# General functions to generate S when n = 2, form Wakeley's Book ----
# Matrix combines state 4&6 and 5&7
subint_book <- function(rho){
  matrix <- matrix(c(-(rho +1),  rho,  0,  0, 0,
                     1, -(3+(rho/2)),  (rho/2),  1, 1,
                     0,  4, -6,  1, 1,
                     0,  0,  0, -1, 0,
                     0, 0, 0, 0, -1), nrow = 5, byrow = T)
  return(matrix)
}

# Function to generate the sub-intensity matrix, Reward-Matrix and alpha for standard coalescent model
basicCoalModel <- function(n){
  S <- matrix(0,n-1,n-1)
  rate <- choose(rev(2:n), 2)
  diag(S) <- rate
  S <- cbind(rep(0,n-1), S)
  # add -rowSum on diagonal to have rowSum = 0
  diag(S) <- -(rowSums(S))
  rewardMatTotal <- cbind(rev(2:n),rev(2:n))
  init_p <- c(1, rep(0, n-2))
  return(list(subIntMat = S[,-(n)], init_probs = init_p, rewardMatTotal = rewardMatTotal))
}

# Function to calculate covariance of two loci 
cov_2Loci <- function(obj) {
  if (is(obj, 'mult_cont_phase_type')) {
    # select values
    reward <- obj$reward  
    U <- -(solve(obj$subint_mat))
    init_probs <- obj$init_probs
    
    #E[Y_1,Y_2] = αU diag(R_1)U R_2 − αU diag(R_2)U R_1
    E12 <- init_probs %*% U %*% diag(reward[,1])%*% U %*% reward[,2] + init_probs %*% U %*% diag(reward[,2])%*% U %*% reward[,1]
    
    #E[Y_1] = α U R_1
    E1 <- init_probs %*% U %*% reward[,1]
    E2 <- init_probs %*% U %*% reward[,2]
    
    
    cov <- E12 - E1*E2
    return(cov)
  } 
  else {
    stop("Please provide an object of class 'mult_cont_phase_type'.")
  }
}

# Function to calculate correlation between no. of segregating sites ----
corr_SFS <- function(MPH_obj, m){
  cov <- var(MPH_obj)
  exp_1 <- mean(MPH_obj, 1)
  return(cov[1,2]/(cov[1,1]+((2/(m*2))*exp_1)))
}


# Wattersons theta: Estimate the mutation rate ----
  # theta = K/a_n, with K: #segregating sites a_n = Sum_i = 1 till n-1_(1/i) 
watterson <- function(n,data){
  K <- mean(data)
  i <- seq(1, n-1)
  sum_i <- sum(1/i)
  result <- K/sum_i
  return(result)
}

###################################################################
# Recombination rate estimation methods

# Recombination rate estimation using the method of moments
  # For optim-function:
MoM_rho_optim <- function(x, n, m_estimate, data_joint){
    # Optimization parameter rho
    rho <- x
    # Calculate the observed correlation in the sample
    corr_observed <- abs(cor(data_joint[,1], data_joint[,2]))
    EvoMdl_test <- IntensityMatrix(n,rho)
    MPH_gen <- MPH(EvoMdl_test [["subIntMat"]], EvoMdl_test [["initProbs"]], EvoMdl_test [["rewardMatTotal"]])
    corr_model<- corr_SFS(MPH_gen, m_estimate)
    # Minimize the absolute error between the theoretical correlation and the empirical correlation
    return(abs(corr_observed-corr_model))
}

# Recombination rate estimation using the maximum likelihood method 

# Functions to calculate probability of an observed sample 
# Function to converts an evolutionary model to the desired ----
ConvertEvoMdlToPrbMdl_charly <- function(EvoMdl, lam){
  SMat <- EvoMdl$subIntMat
  RMat <- EvoMdl$rewardMatTotal
  # Sub-transition probability matrix P_total
  P_Total <- solve( -SMat+lam*( diag( rowSums( RMat ) ) ) )
  # Exit vector
  nSt <- nrow( SMat )
  sVec <- -SMat %*% rep(1,nSt)                 
  exitVec <- P_Total %*% sVec
  # Marginal sub-transition probability matrices
  PrbArray <- array( 0, dim = c( ncol(RMat), rep(nSt, 2) ) )
  for (i in 1:ncol(RMat)){
    PrbArray[i,,] <- P_Total %*% ( lam*diag( RMat[,i] )) 
  }
  out <- list()
  out$initVec <- EvoMdl$initProbs
  out$PrbArray <- PrbArray
  out$exitVec <- exitVec
  return(out)
}

# Function to calculate Probability of the sample ----
SmplPrb <- function(smpl,PrbMdl){
  initVec <- PrbMdl$initVec
  PrbArray <- PrbMdl$PrbArray
  exitVec <- PrbMdl$exitVec
  # If zero mutations
  if (sum(smpl)==0) return( as.double( initVec %*% exitVec) )
  # If positive number of mutations
  nSt <- length(initVec)
  mSet <- multiset( c( rep( 1:length(smpl), smpl ) ) )  
  sum.prb <- 0
  for (k in 1:ncol(mSet)){
    ProdMat <- diag(nSt)
    for (m in 1:nrow(mSet)){
      ProdMat <-  ProdMat %*% PrbArray[mSet[m,k],,]
    }
    sum.prb <- sum.prb + initVec %*% ProdMat %*% exitVec
  }
  return( as.double( sum.prb ) )
}

# Recombination rate estimation using the maximum likelihood method ----
# For optim-function:
ML_prob <- function(x, n, m_estimate, data_joint) {
  # Optimization parameter rho
  rho <- x
  EvoMdl_test <- IntensityMatrix(n,rho)
  PrbMdl <- ConvertEvoMdlToPrbMdl_charly(EvoMdl_test, m_estimate)
  # Add up prob. to observe whole data set: call function 'SmplPrb' for each sample
  prob_sample <- vector()
  # Taking the log of the probabilities: Convert the product of small numbers into the sum of their logarithms.
  for(row in 1:nrow(data_joint)){
    prob_sample[row] <- SmplPrb(data_joint[row,], PrbMdl)%>%log
  }
  return(-(sum(prob_sample)))
}



