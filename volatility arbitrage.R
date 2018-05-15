library(latticeExtra)
library(lattice)
library(pbo)
library(PerformanceAnalytics)
ms<-read.csv("MS.csv")
unh<-read.csv("UNH.csv")
amzn<-read.csv("AMZN.csv")
emr<-read.csv("EMR.csv")
aep<-read.csv("AEP.csv")
xom<-read.csv("XOM.csv")
hrg<-read.csv("HRG.csv")
aapl<-read.csv("AAPL.csv")
goog<-read.csv("GOOG.csv")

##define a function that calculate the volatility signature
Freq=30
Volsig <-function(P,Freq){                #P = stock price set, freq = frequency
    annstdev<-vector()
    for (j in 1:as.numeric(Freq)){
        logret<-vector()
        for(i in 1:length(P)-j){
            logret[i]<-log(P[i+j])-log(P[i])
        }
        annstdev[j]<-sqrt(252/j)*sd(logret,na.rm=TRUE)
    }
    plot(100*annstdev,type='o',col='blue',main=expression('Volatility Signature'),xlab='Frequency',ylab='Annualized Volatility')
}

##plot the volatility signature of spx 5yr and 2yr data

#...All Volsig


##Define the daily PnL funtion
PnL<-function(P,Freq,alpha=1)#P=stock price set,freq=frequency,alpha=initial investment
{
    rev<-rep(0,length(P)-1)
    for (i in 1:length(P)-1){
        j=Freq*(ceiling(i/Freq)-1)+1
        rev[i]=alpha*(P[i+1]-P[i])*(1/P[i]-1/P[j])
        next
    }
    rev
}
accrev<-function(rev){
    acc<-vector()
    for (i in 1:length(rev)){
        acc[i]<-sum(rev[1:i])
        next
    }
    acc
}

#use Volatility Arbitrage strategy
unhret<-accrev(PnL(unh$Close,30,1))
msret<-accrev(PnL(ms$Close,30,1))
amznret<-accrev(PnL(amzn$Close,30,1))
aepret<-accrev(PnL(aep$Close,30,1))
emrret<-accrev(PnL(emr$Close,30,1))
xomret<-accrev(PnL(xom$Close,30,1))
hrgret<-accrev(PnL(hrg$Close,30,1))
aaplret<-accrev(PnL(aapl$Close,30,1))
googret<-accrev(PnL(goog$Close,30,1))

Mat$msret<-msret
Mat$unhret<-unhret
Mat$amznret<-amznret
Mat$aepret<-aepret
Mat$emrret<-emrret
Mat$xomret<-xomret
Mat$hrgret<-hrgret
Mat$aaplret<-aaplret
Mat$googret<-googret

T<-1248
N<-9
s<-8
volpbo<-pbo(Mat[3:1250,],s,Omega,threshold=1)
#Minimum Backtest Length
minBTL<- function(Emax=1.0,N) {
    
    gama<-0.5772156649        #Euler-Mascheroni constant
    minBTL<-(((1-gama)*qnorm(1-1/N)+gama*qnorm(1-1/N*exp(-1)))/Emax)**2
    minBTL
    
}

minBTLcon<- function(Emax=1.0,N) {
    minBTLcon<-2*log(N)/Emax
    minBTLcon
}
exp<-c()
max<-c()
for (i in 1:1000){
    exp[i]<-minBTL(1.0,i)
    max[i]<-minBTLcon(1.0,i)
    next
}
length<-c(1:1000)
plot(N,exp,main="Minimum Back Test Length",type="l",xlab="Number of Trials (N)",ylab="Mimunum Backtest Length (in Year)")
line(N,max)


#after construct return matrix Mattest
mattest<-mattest1
t <- nrow(mattest)             # samples per study
n <- ncol(mattest)             # studies
s<-8
cs <- combn(s,s/2)       # combinations
sn <- t / s              # partition size

test_config <- bquote(N == .(n) ~~ T == .(t) ~~ S == .(s))
cs_results <- list()

cs_compute <- function(csi) {
  # partition indices
  is_i <- cs[,csi]
  
  # in-sample indices
  is_indices <- as.vector(sapply(is_i,function(i) {
    start <- sn * i - sn + 1
    end <- start + sn - 1
    start:end
  }))
  
  # out-of-sample indices
  os_indices <- setdiff(1:t,is_indices) 
  
  # training and test sets (in sample, out of sample)
  # choosing not to reassign row names of each to 1:(T/2)
  # after R don't need to save J or J_bar so could skip this assignment
  j <- mattest[is_indices,]
  j_bar <- mattest[os_indices,]
  
  # compute performance over the N strategies in each subset
  # could use for R any summary statistic e.g. SharpeRatio or Omega
  # r <- mapply(SharpeRatio,as.ts(j))
  # r_bar <- mapply(SharpeRatio,j_bar)
  #r <- Omega(as.ts(j))
  #r_bar <- Omega(j_bar)
  r<-vector()
  for(i in 1:9){
    r[i]<-mean(j[,i])/sqrt(var(j[,i]))
  }
  r_bar<-vector()
  for(i in 1:9){
    r_bar[i]<-mean(j_bar[,i])/sqrt(var(j_bar[,i]))
  }
  
  # compute n* by argmax over R vector
  n_star <- which.max(r)
  n_max_oos <- which.max(r_bar)
  
  # rank of n*th result from OOS performance; converted to (0,1) interval
  os_rank <- rank(r_bar)[n_star]
  omega_bar_c <- os_rank / (length(r_bar)+1)
  
  # logit
  # note the value can be Inf
  lambda_c <- log(omega_bar_c / (1 - omega_bar_c))
  
  list(r,r_bar,n_star,n_max_oos,os_rank,omega_bar_c,lambda_c)
}
cs_results <- NULL


for ( csi in 1:70) {
  cs_results <- rbind(cs_results,cs_compute(csi))
}
colnames(cs_results) <- c("R","R_bar","n*","n_max_oos","os_rank","omega_bar","lambda")
rownames(cs_results) <- 1:ncol(cs)
lambda <- as.numeric(cs_results[,"lambda"])
lambda[which(lambda==Inf)] <- 6 # for plotting
phi <- sum(ifelse(lambda<=0,1,0))/ncol(cs)
rn_pairs <- as.data.frame(do.call(rbind,lapply(1:ncol(cs),function(i) {
  n <- cs_results[[i,3]]
  r <- cs_results[[i,1]]
  rb <- cs_results[[i,2]]
  return(c(r[n],rb[n]))
})))
colnames(rn_pairs) <- c("Rn","Rbn")
linear_fit <- lm(rn_pairs)

m <- signif(as.numeric(linear_fit$coefficients[1]),digits=5) # slope
b <- signif(as.numeric(linear_fit$coefficients[2]),digits=5) # intercept
ar2 <- signif(summary(linear_fit)$adj.r.squared,digits=2) # adj R-squared
p_oos_bt <- signif(length(which(rn_pairs$Rbn<0)) / 
                     nrow(rn_pairs),digits=3)
rv<-list(
  results=cs_results,
  combos=cs,
  lambda=lambda,
  phi=phi,
  rn_pairs=rn_pairs,
  func=as.character(substitute(f)),
  slope=m,
  intercept=b,
  ar2=ar2,
  threshold=0,
  below_threshold=p_oos_bt,
  test_config=test_config,
  inf_sub=6)

class(rv) <- "pbo"
summary.pbo <- function(object,...) {
  writeLines(c(paste("Performance function",
                     object$func,
                     "with threshold",
                     object$threshold),
               ""))
  results <- c(object$phi,object$slope,object$ar2,object$below_threshold)
  names(results) <- c("p_bo","slope","ar^2","p_loss")
  results
}

summary(rv)
histogram(rv)
dotplot(rv,pch=15,col=2,cex=1.5)
xyplot(rv,plotType="cscv",cex=0.8,show_rug=FALSE,osr_threshold=100)
xyplot(rv,plotType="degradation")
xyplot(rv,plotType="dominance",lwd=2)
xyplot(rv,plotType="pairs",cex=1.1,osr_threshold=75)
xyplot(rv,plotType="ranks",pch=16,cex=1.2)
xyplot(rv,plotType="selection",sel_threshold=0,cex=1.2)
