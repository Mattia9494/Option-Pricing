american.put.approx <- function(S,K,t,r,sigma,n){
  # Valuation of an American put option in the Black-Scholes model
  # using a binomial tree approximation
  if(S <= 0 || K < 0 || t <= 0 || r < 0 || sigma <= 0)
    return("Invalid parameters: S,t,sigma > 0, and K,r >= 0 must hold!")
  # Use the Normal scaling to ensure convergence to the geometric Brownian motion
  u <- exp(sigma*sqrt(t/n))
  d <- exp(-sigma*sqrt(t/n))
  #we use r, not mu because it is a risk-neutral measure
  q <- 1/2+1/2*(r-sigma^2/2)/sigma*sqrt(t/n)
  # possible stock prices at time t_n=t
  stock <- S*u^(n:0)*d^(0:n)
  # put option prices max(K-S_n,0) at time t_n
  option <- pmax(K-stock,0)
  # calculate the put option prices iteratively backwards
  for(i in (n-1):0){
    # discounted expected option price with respect to q at time t_i
    option <- exp(-r*t/n)*(q*option[-(i+2)]+(1-q)*option[-1])
    # possible stock prices at time t_i
    stock <- stock[-1]/d
    # intrinsic values max(K-S_i,0) at time t_i
    iv <- pmax(K-stock,0)
    # if the intrinsic value is higher than the discounted expected option price,
    # early exercise is optimal --> replace the corresponding values in the
    # option vector by the intrinsic values
    option[option<iv] <- iv[option<iv]
  }
  return(option)
}

#Although the put prices computed with the above program seem to converge fast, they can
#never be exact, so for a better overview practitioners might want to determine 
#an interval or price range in addition in which the true put price should lie.

american.put.boundaries <- function(S,K,t,r,sigma,m=500,n=1000){
  # Upper and lower bound for the value of an American put option in the Black-Scholes
  # model based on path simulation
  if(S <= 0 || K < 0 || t <= 0 || r < 0 || sigma <= 0)
    return("Invalid parameters: S,t,sigma > 0, and K,r >= 0 must hold!")
  # Simulate m possible stock price paths under the risk-neutral measure, calculate
  # the corresponding option values and store these in the vectors lower and upper
  # initialize lower and upper vectors
  lower <- numeric(m)
  upper <- numeric(m)
  for(i in 1:m){
    stock <- S*exp(cumsum(c(0,rnorm(n,mean=(r-sigma^2/2)*t/n, sd=sigma*sqrt(t/n)))))
    # compute the intrinsic values for every point in time t_i and compare them with
    # the BS-European put prices
    iv <- pmax(K-stock,0)
    europ.put <- bs.option.price(stock,K,r,seq(t,0,-t/n),sigma,type="put")
    compare <- iv-europ.put
    # If the intrinsic values are always smaller than the European option prices,
    # it is never optimal to exercise early
    if(sum(compare>0) == 0){
      lower[i] <- exp(-r*t)*max(K-stock[n+1],0)
      upper[i] <- lower[i]
    }
    else{
      # determine the smallest early exercise time
      ex.time <- seq(0,t,t/n)[compare>0][1]
      ex.value <- iv[compare>0][1]
      # one put price bound is the discounted intrinsic value at the smallest
      # exercise time
      lower[i] <- exp(-r*ex.time)*ex.value
      # determine the latest early exercise time
      len <- sum(compare>0)
      ex.time <- seq(0,t,t/n)[compare>0][len]
      ex.value <- iv[compare>0][len]
      # the other put price bound is the discounted intrinsic value at the largest
      # exercise time
      upper[i] <- exp(-r*ex.time)*ex.value
      ## if upper[i] < lower[i], then switch both values
      if(upper[i] < lower[i]){
        dummy <- lower[i]
        lower[i] <- upper[i]
        upper[i] <- dummy
      }
    }
  }
  return(c(mean(lower),mean(upper)))
}
