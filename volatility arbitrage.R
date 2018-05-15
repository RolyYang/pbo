library(latticeExtra)
library(lattice)
library(pbo)
library(PerformanceAnalytics)

#import the stocks/index historical data.
#ms<-read.csv("MS.csv")
#unh<-read.csv("UNH.csv")
#amzn<-read.csv("AMZN.csv")
#emr<-read.csv("EMR.csv")
#aep<-read.csv("AEP.csv")
#xom<-read.csv("XOM.csv")
#hrg<-read.csv("HRG.csv")
#aapl<-read.csv("AAPL.csv")
#goog<-read.csv("GOOG.csv")

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

