#minimum backtest length
minBTL<- function(Emax=1.0,N) {
	
	gama<-0.5772156649		#Euler-Mascheroni constant
	minBTL<-(((1-gama)*qnorm(1-1/N)+gama*qnorm(1-1/N*exp(-1)))/Emax)**2
  minBTL
	
}

#the shell of minBTL
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
plot(length,exp,main="Minimum Back Test Length",type="l",xlab="Number of Trials (N)",ylab="Mimunum Backtest Length (in Year)")
line(length,max)

N<-9
minBTLcon(1.0,9) #4.39445
