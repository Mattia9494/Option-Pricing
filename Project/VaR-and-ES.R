#VaR and ES
allianz1 <- read.table("allianz_Jan2005-Jul2009.txt",header=TRUE)
allianz2 <- read.table("allianz_July2010-July2015.txt",header=TRUE)

VaR.ES <- function(data,ES=TRUE,width=501,stepsize=1,level=0.99,conf.level=0.05,
                   summary=TRUE,plot=TRUE,plot.info=TRUE,name="",lwd=2){
  ## Moving VaR-forecasts and backtesting for the dataset data
  ## data is assumed to be a sample of a geometric Brownian motion observed at equidistant
  ## time points
  if(any(is.na(data)) || any(is.infinite(data)) || any(data <= 0))
    return("data must be a vector with components 0 < data[i] < Inf")
  if(width < 2) 
    return("width of moving window must be >= 2")
  if(stepsize < 1) 
    return("stepsize must be >= 1")
  if(width > length(data) || stepsize > length(data))
    return("width and stepsize must be smaller than length(data)")
  if(any(c(level,conf.level) <= 0) || any(c(level,conf.level) >= 1))
    return("0 < level,conf.level < 1 must hold!")
  n <- length(data)
  window.number <- floor((n-width)/stepsize)
  VaR <- numeric(window.number)
  es <- numeric(window.number)
  FoEL <- numeric(window.number)
  loss <- numeric(window.number)
  for(i in 1:window.number){
    index <- (1:width)+(i-1)*stepsize
    logreturns <- diff(log(data[index]))
    mu <- mean(logreturns) #a little confusing because he changes how it's written
    sigma <- sd(logreturns) #correct like in the formula
    qa <- qnorm(1-level,mean=mu,sd=sigma)
    ## VaR- and ES-forecast for the next data point data[i+1]
    VaR[i] <- data[index[length(index)]]*(1-exp(qa))
    es[i] <- data[index[length(index)]]*(-1+exp(mu+sigma^2/2)*pnorm(qnorm(1-level)+sigma)/(1-level))
    ## realized loss (potential gains are treated as zero)
    loss[i] <- max(data[index[length(index)]]-data[index[length(index)]+1],0)
    ## Check if excessive loss has occured
    FoEL[i] <- as.numeric(loss[i] >= VaR[i])
  }
  ## FoEL-Test (H_0: Frequency of excessive losses is 1-level)
  p0 <- 1-level
  h <- sum(FoEL)
  p.hat <- h/window.number
  ## Under H_0: p.hat=p0 the statistic stat is (asymptotically) chisq(1)-distributed
  stat <- -2*log(((1-p0)^(window.number-h))*(p0^h))+2*log(((1-p.hat)^(window.number-h))*(p.hat^h))
  ## p-value
  p.value <- 1-pchisq(stat,df=1)
  ## confidence interval for p0
  factor <-  sqrt((p.hat*(1-p.hat))/window.number)
  conf.int <- c(p.hat-qnorm(1-0.5*conf.level)*factor,p.hat+qnorm(1-0.5*conf.level)*factor)
  conf.int <- round(conf.int*100,2)
  ## test statistic for expected shortfall
  V <- sum(loss[loss >= VaR]-es[loss >= VaR])/h
  if(summary){
    if(ES){
      ## Summary for expected shortfall
      cat("\n\t\t\tResults for ES")
      cat("\n\nES-level   excessive losses   statistic V\n",sep = "")
      cat(level*100,"%        ",h,"                  ",round(V,4),"\n",sep = "")
    } 
    else{
      ## Print a summary of the FoEL-test on the screen
      cat("\n\t\t\tResult of the FoEL-test")
      cat("\n\nVaR-level   Frequency   conf. interval     LQ-statistic      p-value   H0\n",
          sep = "")
      cat(level*100,"%         ",round(h/window.number*100, 2),"%           ",sep = "")
      cat("[",conf.int[1],"%, ",conf.int[2],"%]     ",round(stat,4),"          ",
            round(p.value*100,2),"%        ",sep = "")
      if(conf.level<p.value) cat("not rejected\n")
      else cat("rejected\n")
      cat("conf.level of LQ-test: ",round(conf.level*100,1),"%    (of the conf. interval: ",
            round((1-conf.level)*100,1),"%)\n",sep = "")
  } 
}
if(plot){
 ## save previous graphical parameters
 orig.mar <- par()$mar
 ## set appropriate plot margins
 par(mar=c(5+as.numeric(plot.info)+2,4,5,1))
 ## set y-axis limits
 if(ES) 
   ylim <- c(min(c(0,es)),max(c(loss,es)))
 else 
   ylim <- c(min(c(0,VaR)),max(c(loss,VaR)))
 x.txt <- "trading day"
 x.time <- stepsize*(0:(window.number-1))+width+1
 plot(x.time,loss,ylim=ylim,type="n",xlab="",ylab="")
 lines(c(x.time[1],x.time[length(x.time)]),c(0,0))
if(ES) 
  lines(x.time,es,lwd=lwd,col="red")
else 
  lines(x.time,VaR,lwd=lwd,col="red")
lines(x.time,loss,type="h")
## Mark times of excessive losses on the upper margins
usr <- par()$usr[c(3,4)]
if(ES){
  points(x.time[loss > es],rep(usr[2]-(usr[2]-usr[1])*0.015,sum(loss > es)),
            pch="|",col="red")
  points(x.time[loss >= VaR & loss <= es],rep(usr[2]-(usr[2]-usr[1])*0.015,
            sum(loss >= VaR & loss <= es)),pch="|",col="blue")
}
else 
  points(x.time[as.logical(FoEL)],rep(usr[2]-(usr[2]-usr[1])*0.015,sum(FoEL)),pch="|")
  ## x-axis labeling
  mtext(x.txt,side=1,line=2.5)
  if(plot.info){
    ## optional additional information (printed below the x-axis)
    mtext(paste("VaR-level: ",round(100*level,1),"%",sep=""),side=1,line=4.5)
    x.txt.3 <- paste("exc.loss-frequency: ",round(p.hat,4)*100,"%",sep = "")
    x.txt.4 <- paste(round((1-conf.level)*100,1),"%-confidence intervall: [",
                     conf.int[1],"%, ",conf.int[2],"%]",sep="")
    x.txt.5 <- paste("p-value: ",round(p.value*100,2),"%",sep="")
    if(conf.level<p.value) x.txt.6 <- "(not rejected)"
    else 
      x.txt.6 <- "(rejected)"
    if(ES) mtext(paste(x.txt.3, "           ES-statistic V: ", round(V,4),sep=""),
            side=1,line=6)
    else mtext(paste(x.txt.3, "            ",x.txt.5," ",x.txt.6,sep=""),
                      side=1,line=6)
    ## detailed y-axis labeling
    if(ES) mtext("ES (red line) and actual losses (vertical lines)",side=2,line=2.5)
    else mtext("VaR (red line) and actual losses (vertical lines)",side=2,line=2.5)
    ## detailed information on top
    mtext(paste("times of excessive losses  (total: ",sum(FoEL),")",sep=""),line=0.5)
    if(ES) title.txt <- "Losses and ES-forecasts for"
    else title.txt <- "Losses and VaR-forecasts for"
    if(name=="") title.txt <- paste(title.txt, "data")
    else title.txt <- paste(title.txt, name)
    mtext(title.txt,line=2,cex=1.2)} else{
    ## sparse annotation if plot.info=FALSE
    if(ES) mtext("ES and actual losses",side=2,line=2.5)
    else mtext("VaR and actual losses",side=2,line=2.5)
    mtext("times of excessive losses",line=0.5)
    if(ES) title.txt <- "Losses and ES-forecasts for"
    else title.txt <- "Losses and VaR-forecasts for"
    if(name=="") title.txt <- paste(title.txt, "data")
    else title.txt <- paste(title.txt, name)
    mtext(title.txt,line=2,cex=1.2)
    }
## reset graphical parameters to previous values
par(mar=orig.mar)
  }
invisible(list(VaR=VaR,ES=es,losses=loss,FoEL=FoEL,VaR.level=level,conf.level=conf.level,
                     frequency=p.hat,conf.int=conf.int,stat= stat,p.value=p.value,
                     H0=c("not rejected","rejected")[c(conf.level<p.value,conf.level>=p.value)],
                     V=V))
}
    
allianz1.var.results <- VaR.ES(allianz1$Close)
