#The multivariate geometric Brownian motion
#Simulation and parameter estimation

#Simulate multivariate normal distribution
simulate.mv.normal <- function(n,mu,Sigma){
  # Simulation of multivariate normal random vectors X ~ N(mu,Sigma)
  # n: number of random vectors to be simulated
  # mu: mean vector
  # Sigma: covariance matrix
  # to be a symmetric square matrix, it needs to have equal dimensions 
  # and needs to be equal to it transpose
  if(diff(dim(Sigma))!=0 || !all(Sigma==t(Sigma)))
    return("Sigma is not a symmetric square matrix!")
  if(length(mu) != dim(Sigma)[1])
    return("Incompatible dimensions of mu and Sigma!")
  ev <- eigen(Sigma)
  #we want a positive semi-definite matrix therefore negative eigenvalues are not accepted
  if(any(ev$values<0)) 
    return("Sigma is not positive semi-definite!")
  # A is a root of Sigma
  A <- ev$vectors %*% diag(sqrt(ev$values))
  M <- matrix(numeric(n*length(mu)),length(mu),n)
  for(i in 1:n) 
    M[,i] <- A %*% rnorm(length(mu)) + mu
  return(M) 
}

#Simulate geometric brownian motion
simulate.mv.geomBB <- function(S0,mu,Sigma,t=10,m=1000){
  # Simulation of a geometric Brownian motion with starting vector S0 on [0,t]
  # m defines the number of increments to be simulated
  if(diff(dim(Sigma))!=0 || !all(Sigma==t(Sigma)))
    return("Sigma is not a symmetric square matrix!")
  if(length(mu) != dim(Sigma)[1] | length(S0) != dim(Sigma)[1])
    return("Incompatible dimensions of S0, mu, and Sigma!")
  ev <- eigen(Sigma)
  if(any(ev$values < 0)) 
    return("Sigma is not positive semi-definite!")
  # A is a root of Sigma
  A <- ev$vectors %*% diag(sqrt(ev$values))
  M <- matrix(numeric((m+1)*length(mu)),length(mu),m+1)
  for(i in 2:(m+1)) #2:m+1 because the first column is full of zeros because e^0=1
    M[,i] <- sqrt(t/m)*A %*% rnorm(length(mu)) + t/m*(mu-diag(Sigma)/2)
  return(S0*exp(t(apply(M,1,cumsum))))
}
