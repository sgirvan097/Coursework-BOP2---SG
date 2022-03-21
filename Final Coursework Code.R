
# Function to determine min per arm sample size taking account
# of type I and II error rates of 0.05 and 0.2 respectively

Min_Sample_Size <- function(theta,delta,alpha,beta){
  
  n <- (theta*(1-theta)/(delta^2))*((qnorm(1-alpha)-qnorm(beta))^2)
  
  return(n)
}

N_min <-  Min_Sample_Size(0.5,0.2,0.05,0.2)
N_min


# Note the current time.
ptm <- proc.time()
# Call the function several times.
y2 <- replicate(1e2, Min_sample_size(5,0.2,0.05,0.2))
# Print how much time has elapsed.
proc.time() - ptm
# A very short run time can be observed for the per
# arm sample size calculation

############################################################################

Error_Function <- function(x,n1,n2,theta,a0,b0) {
  
  lambda <- x[1]
  gamma  <- x[2]
  
  #Set number of simulations
  M <- 1e4
  
  #Simulate number of successes for stage 1
  y1 <- rbinom(M, n1, theta)
  
  
  # Get posterior Beta(a1, b1) parameters.
  a1 <- a0 + y1
  b1 <- b0 + n1 - y1
  
  # Probability of futility.
  fut1 <- pbeta(theta, a1, b1)
  
  # Threshold to determine progression, based on the decision rule.
  C1 <- 1 - lambda * (n1 / n2)^gamma
  
  #return number of simulations where threshold is exceeded for stage 1
  s1<- sum(fut1<C1)
  
  #Simulate number of successes for stage 2 given stage 1 is passed
  y2 <- rbinom(s1, n2 - n1, theta)
  
  # Get posterior Beta(a2, b2) parameters.
  a2 <- a0 + y2 + y1
  b2 <- b0 + n2 - y2 - y1
  
  # Probability of futility2.
  fut2 <- pbeta(theta, a2, b2)
  
  # Threshold to determine progression, based on the decision rule.
  C2 <- 1 - lambda * (n2 / n2)^gamma
  
  #return number of simulations where threshold is exceeded for stage 2
  #i.e probability of type 1 error
  s2<- sum(fut2<C2)
  
  return(sum(fut1 < C1 & fut2 < C2)  /M)
}

#############################################################################

# We want to evaluate a grid of possible values of lambda and gamma
grid_search <- expand.grid(lambda = seq(0, 2, 0.05),gamma = seq(0.1, 1, 0.1))

# Functions to determine which values of lambda and gamma satisfy
# the restrictions on error rate

a<- apply(grid_search, 1, Error_Function,n1=ceiling(N_min),n2=2*ceiling(N_min), theta=0.5,
          a0=0.5, b0=5)
b<-  apply(grid_search, 1, Error_Function,n1=ceiling(N_min),n2=2*ceiling(N_min), theta=0.7,
           a0=0.5, b0=5)

Type_I_II_Error<-cbind(a,b)

Lamb_Gam_sets <- which(Type_I_II_Error[, 1] >0 & Type_I_II_Error[, 1] <= 0.05 
                       & Type_I_II_Error[, 2] >0 & Type_I_II_Error[, 2] <= 0.2)

Type_I_II_Error[Lamb_Gam_sets,]

updated_grid_search <- grid_search[Lamb_Gam_sets,]


# Note the current time.
ptm <- proc.time()
# Call the function several times.
y <- replicate(10, apply(grid_search, 1, Error_Function,n1=ceiling(N_min),
                         n2=2*ceiling(N_min), theta=0.5,
                         a0=0.5, b0=5))
# Print how much time has elapsed.
proc.time() - ptm
# A time of 4.6 seconds can be observed
# for the lambda and gamma selection calculation

###############################################################################

set.seed(283765)

evaluate_design <- function(lambda, gamma, n1, n2) {
  
  # Estimate the expected sample size of a design defined by its 
  # decision rule parameters (lambda, gamma) and sample size 
  # parameters (n1, n2), along with its standard error.
  
  # Set the number of simulations.
  M <- 10^4
  # Create an empty vector to store simulated NS.
  Ns <- rep(NA, M)
  for (i in 1:M) {
    # Simulate theta from its prior, and then the stage 1 data conditional
    # on this theta.
    theta <- rbeta(1, 0.5, 0.5)
    y1 <- rbinom(1, n1, theta)
    
    # Get posterior Beta(a1, b1) parameters.
    a1 <- 0.5 + y1
    b1 <- 0.5 + n1 - y1
    
    # Probability of futility.
    fut1 <- pbeta(0.5, a1, b1)
    
    # Threshold to determine progression, based on the decision rule.
    C1 <- 1 - lambda * (n1 / n2)^gamma
    
    # Note the final total sample size and store in the vector Ns.
    if (fut1 > C1) {
      Ns[i] <- n1
    } else {
      Ns[i] <- n2
    }
  }
  
  
  # Return the estimated expected sample size and its estimated standard error.
  return(c(mean(Ns), sqrt(var(Ns)/M)))
}

###############################################################################

# For each pair in the grid, find the expected sample size and 
# store in a vector.
exp_n <- rep(NULL, nrow(updated_grid_search))
for(i in 1:nrow(updated_grid_search)) {
  exp_n[i] <- evaluate_design(updated_grid_search[i, 1], updated_grid_search[i, 2], 
                              n1 = ceiling(N_min), n2 = 2*ceiling(N_min))
}

# Note the current time.
ptm <- proc.time()
#run the min expected sample size function
exp_n <- rep(NULL, nrow(updated_grid_search))
for(i in 1:nrow(updated_grid_search)) {
  exp_n[i] <- evaluate_design(updated_grid_search[i, 1], updated_grid_search[i, 2], 
                              n1 = ceiling(N_min), n2 = 2*ceiling(N_min))[1]
}
# Print how much time has elapsed.
proc.time() - ptm
# A run time of around 2.2 seconds is observed


# Min sample size given it is greater than n1 
# (as that would mean no trial makes it to stage 2)
k <- min(exp_n[exp_n>ceiling(N_min)])

# Position of lambda and gamma value in vector()
p <- which(exp_n == k)

# Min lambda and gamma values along with min expected sample size
updated_grid_search[p,]
k

##############################################################################

# Condensed grid search after initial values of lambda = 0.95 and gamma = 0.1 
# were found to minimize the expected sample size
grid_search2 <- expand.grid(lambda = seq(0.8, 1, 0.01),
                            gamma = seq(0.01, 0.2, 0.01))


# Functions to determine which values of lambda and gamma satisfy
# the restrictions on error rate
a2<- apply(grid_search2, 1, Error_Function,n1=ceiling(N_min),n2=2*ceiling(N_min), theta=0.5,
           a0=0.5,b0=5)
b2<-  apply(grid_search2, 1, Error_Function,n1=ceiling(N_min),n2=2*ceiling(N_min), theta=0.7,
            a0=0.5,b0=5)


Type_I_II_Error2<-cbind(a2,b2)

Lamb_Gam_sets2 <- which(Type_I_II_Error2[, 1] >0 & Type_I_II_Error2[, 1] <= 0.05 
                        & Type_I_II_Error2[, 2] >0 & Type_I_II_Error2[, 2] <= 0.2)

Type_I_II_Error2[Lamb_Gam_sets2,]

updated_grid_search2 <- grid_search2[Lamb_Gam_sets2,]

# Note the current time.
ptm <- proc.time()
# Call the function several times.
y <- replicate(1, apply(grid_search2, 1, Error_Function,n1=ceiling(N_min),
                        n2=2*ceiling(N_min), theta=0.5,
                        a0=0.5, b0=5))
# Print how much time has elapsed.
proc.time() - ptm
# A time of 5.1 seconds can be observed
# for the lambda and gamma selection calculation

##############################################################################

# For each pair in the grid, find the expected sample size and 
# store in a vector.
exp_n2 <- rep(NULL, nrow(updated_grid_search2))
for(i in 1:nrow(updated_grid_search2)) {
  exp_n2[i] <- evaluate_design(updated_grid_search2[i, 1], updated_grid_search2[i, 2],
                               n1 = ceiling(N_min), n2 = 2*ceiling(N_min))
}

# Note the current time.
ptm <- proc.time()
# Run the min expected sample size function
exp_n2 <- rep(NULL, nrow(updated_grid_search2))
for(i in 1:nrow(updated_grid_search2)) {
  exp_n2[i] <- evaluate_design(updated_grid_search2[i, 1], updated_grid_search2[i, 2], 
                               n1 = ceiling(N_min), n2 = 2*ceiling(N_min))[1]
}
# Print how much time has elapsed.
proc.time() - ptm
# A run time of around 29 seconds is observed due to the 
# increased number of pairs to be simulated



# Min sample size given > than n1
k2 <- min(exp_n2[exp_n2>ceiling(N_min)])

# Position of lambda and gamma value in vector()
p2 <- which(exp_n2 == k2)

# Min lambda and gamma values along with min expected sample size
updated_grid_search2[p2,]
k2

