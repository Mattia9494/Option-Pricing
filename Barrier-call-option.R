#Simulate the possible price paths and determine the final payoff 
#of up-and-out barrier call with the Black-Scholes model

plot.barrier.call <- function(S,mu,sigma,K,Ku,n=1000,t=1){
  #Simulation of a geometric Brownian path for barrier call valuation
  if(S <= 0 || sigma <= 0 || t <= 0 || K < 0 || Ku <= S || Ku <= K)
    return("Invalid parameters: The relations S,t,sigma>0, K>=0 and Ku>S,K must hold!")
  #Simulate a stock price path
  bm <- cumsum(c(0,rnorm(n,mean = (mu-sigma^2/2)*t/n,sd = sigma*sqrt(t/n))))
  stock <- S*exp(bm)
  plot(0:n,stock,type="n",xaxt="n",xlab="t",ylab="stock price",ylim=c(min(stock),Ku))
  #Mark the barrier level (in red) and the strike price (dashed)
  abline(h=Ku,col="red")
  abline(h=K,lty=3)
  if(all(stock < Ku)){
    # stock price does not hit the upper barrier
    lines(0:n,stock,col="blue")
    title(paste("Up-and-out call with final value", round(max(stock[n+1]-K,0),2)))
  }
  else{
    # stock price hits the upper barrier
    # determine the hitting time index hti
    hti <- (1:(n+1))[Ku-stock <= 0][1]
    lines(0:(hti-1),stock[1:hti],col="blue")
    title("Up-and-out call with final value 0")
  }
}
plot.barrier.call(100,0,1,50,200)


# determine how many simulations are performed to calculate the option price.

evaluate.barrier.call <- function(S,r,sigma,K,Ku,n=1000,t=1,m=1000){
  #Simulation-based barrier call option pricing
  if(S <= 0 || sigma <= 0 || t <= 0 || K < 0 || r < 0 || Ku <= S || Ku <= K)
    return("Invalid parameters: The relations S,t,sigma>0, K,r>=0 and Ku>S,K must hold!")
  #Simulate m possible stock price paths, calculate the corresponding option value and
  #store these in the vector option
  #initialize option vector
  option <- numeric(m)
  #loop for simulations
  for(i in 1:m){
    #simulate a stock price path foe a geometric bm
    bm <- cumsum(c(0,rnorm(n,mean=(r-sigma^2/2)*t/n,sd=sigma*sqrt(t/n))))
    stock <- S*exp(bm)
    #stock price does not hit the upper barrier
    if(all(stock < Ku)) 
      option[i] <- max(stock[n+1]-K,0)
    else 
      option[i] <- 0
  }
  #return the discounted average of all obtained prices as result
  return(exp(-r*t)*mean(option))
}
evaluate.barrier.call(100,0.05,1,50,200)
