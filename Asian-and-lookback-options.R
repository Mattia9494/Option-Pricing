evaluate.asian.option <- function(S,K,r,t,sigma,mu=r,type="call",strike="fixed",m=1000,n=1000){
  # Simulation-based valuation of Asian Options in the Black-Scholes model
  if(S <= 0 || r < 0 || t <= 0 || sigma <= 0)
    return("Invalid parameters: S,t,sigma > 0 and r >= 0 must hold!")
  if(type != "call" && type != "put") 
    return("type must be either ’call’ or ’put’!")
  if(strike != "fixed" && strike != "float")
    return("strike must be either ’fixed’ or ’float’!")
  if(strike=="fixed" && K<0) 
    return("fixed strike K must be >= 0")
  # Simulate m possible stock price paths, calculate the corresponding option value
  # and store these in the vector option
  # initialize option vector
  option <- numeric(m)
  for(i in 1:m){
    stock <- S*exp(cumsum(c(0,rnorm(n,mean=(mu-sigma^2/2)*t/n,sd=sigma*sqrt(t/n)))))
    average <- mean(stock)
    if(type=="call"){
      if(strike=="fixed") 
        option[i] <- max(average-K,0)
      else 
        option[i] <- max(stock[n+1]-average,0)
    }
    else{
      if(strike=="fixed") 
        option[i] <- max(K-average,0)
      else 
        option[i] <- max(average-stock[n+1],0)
    }
  }
  return(exp(-r*t)*mean(option))
}



evaluate.lookback.option <- function(S,K,r,t,sigma,mu=r,type="call",strike="fixed",m=1000,n=1000){
  # Simulation-based valuation of Lookback Options in the Black-Scholes model
  if(S <= 0 || r < 0 || t <= 0 || sigma <= 0)
    return("Invalid parameters: S,t,sigma > 0 and r >= 0 must hold!")
  if(type != "call" && type != "put") 
    return("type must be either ’call’ or ’put’!")
  if(strike != "fixed" && strike != "float")
    return("strike must be either ’fixed’ or ’float’!")
  if(strike=="fixed" && K<0) 
    return("fixed strike K must be >= 0")
  # Simulate m possible stock price paths, calculate the corresponding option value
  # and store these in the vector option
  # initialize option vector
  option <- numeric(m)
  for(i in 1:m){
    stock <- S*exp(cumsum(c(0,rnorm(n,mean=(mu-sigma^2/2)*t/n,sd=sigma*sqrt(t/n)))))
    if(type=="call"){
      if(strike=="fixed") 
        option[i] <- max(max(stock)-K,0)
      else 
        option[i] <- max(stock[n+1]-min(stock),0)
    }
    else{
      if(strike=="fixed") 
        option[i] <- max(K-min(stock),0)
      else 
        option[i] <- max(max(stock)-stock[n+1],0)
    }
  }
  return(exp(-r*t)*mean(option))
}
