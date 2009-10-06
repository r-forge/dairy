mc<-ReproMC()

## -----------------------------------------------------------------------------
## Plot mean progesterone level based on model defined by reproMC1
# 221 states and a rate parameter beta of 0.095 (mean length of a state)
mc1<-mc$reproMC1()
states<-1:mc1$Stepno
plot(states*mc1$beta,mc$meanProgest(states,alpha=3),type="l",xlab="days",ylab="Progesterone")
s<-mc1$tider/mc$beta    # states gaps (in states) between biological events
s<-cumsum(round(s)) # states with biological events
abline(v=mc1$beta*s,lty="dashed")
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
## Simulate progesterone and make plot
mc1<-mc$reproMC1()
states<-1:mc1$Stepno
times<-states*mc1$beta
meanProgest<-mc$meanProgest(1:mc1$Stepno,alpha=3)
stdProgest  <- 0.15     # progest ~ N(meanProgest,stdProgest^2)
dat<-NULL
for (i in (0:2)) {
  progest<-rnorm(length(states),meanProgest,stdProgest)
  dat<-rbind(dat,data.frame(t=times+i*max(times),progest=progest))
}
plot(progest~t,type="l",data=dat,main="3 cycles",xlab="Day",ylab="Progesterone")
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
## Create plot with simulated progesterone and activity
# mar=c(bottom, left, top, right)
par(mar=c(3,4,2,0),oma=c(0,0,0,0),cex=0.75)  #  mar=rep(0,4)
layout(matrix(c(1,2,3,4,5), ncol=1),heights=c(17,17,19,9,9))
xlim<-c(0,45)
## show the regions that have been allocated to each plot
#layout.show(3)

# progesterone
mc1<-mc$reproMC1()
states<-1:mc1$Stepno
curve<-mc$meanProgest(states,alpha=3)
curve<-curve+rnorm(length(curve),0,0.1)
times<-states*mc1$beta
curve<-c(curve,curve,curve)
times<-c(times,max(times)+times,2*max(times)+times)
times<-times[times<=45]
curve<-curve[1:length(times)]
dat<-data.frame(t=times, progst=curve)
plot(progst~t, type="l", data=dat, xlim=xlim, frame=F, lab=c(10,10,7), main="progesterone",
	ylab="progesterone (simulated)",xlab="time")
abline(v=c(0,21,42),col="gray")

# activity
par(mar=c(3,4,1,0))
repro<-Oestus()
tHigh=c(0,1.5,c(0,1.5)+21,c(0,1.5)+21*2)
tLow=c(1.55,c(-0.05,1.55)+21,c(-0.05,1.55)+21*2,45)
dat<-repro$activity(tLow,tHigh,sigma=1,n=20)
plot(act~t,data=dat,type="l", frame=F, ylab="activity (simulated)", main="activity", lab=c(10,5,7))
abline(v=c(0,21,42),col="gray")

# Biological events
par(mar=c(14, 4, 1, 0))
plot(1,1,type="n",xlim=xlim, ylim=c(0,0),xaxt="n", yaxt="n", frame=F, xlab="",ylab="",main="biological events")
at=c(0,0.5,1.5,2,16,c(0,0.5,1.5,2,16)+21,c(0,0.5,1.5,2)+21*2,45)
labels<-c("Oestrogen peak, Heat starts","LH peak","Heat ends","Ovulation start, CL forms","CL regression")
labels<-c(labels,labels,labels,"")
labels<-paste(labels," (",at,")",sep="")
labels<-labels[1:length(at)]
labels[length(labels)]<-""
axis(1, at=at, labels=labels, las=2)

# EJO model
par(mar=c(4, 4, 1, 0))
plot(1,1,type="n",xlim=xlim, ylim=c(0,0),xaxt="n", yaxt="n", frame=F, xlab="",ylab="",main="states simi-MC for reproMC1")
at=c(0,0.5,2,16,c(0,0.5,2,16)+21,c(0,0.5,2)+21*2,45)
axis(1, at=at, las=2)
at<-at[1:(length(at)-1)]
at<-c(0,at[1:(length(at)-1)]+(at[2:length(at)]-at[1:(length(at)-1)])/2)
labels<-c("",expression(S[1]),expression(S[2]),expression(S[3]),expression(S[4]),
	expression(S[1]),expression(S[2]),expression(S[3]),expression(S[4]),
	expression(S[1]),expression(S[2]))
mtext(labels, side=1, line=-1.5, at=at, cex=0.75)

# Jared model
par(mar=c(4, 4, 1, 0))
plot(1,1,type="n",xlim=xlim, ylim=c(0,0),xaxt="n", yaxt="n", frame=F, xlab="",ylab="", main="states model J (bivariate)")
at=c(0,1.5,2,c(0,1.5,2)+21,c(0,1.5,2)+21*2,45)
axis(1, at=at, las=2)
at<-at[1:(length(at)-1)]
at<-c(0,at[1:(length(at)-1)]+(at[2:length(at)]-at[1:(length(at)-1)])/2)
labels<-c("",expression(F[2]),expression(F[3]),"L",     #expression(F[1])
	expression(F[2]),expression(F[3]),"L",
	expression(F[2]),expression(F[3]))
mtext(labels, side=1, line=-1.5, at=at, cex=0.75)
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
## 15 day ahead transition matrix

# Note: For a continuous-time homogeneous Markov process with transition intensity matrix Q,
# the probability of occupying state s at time u + t conditional on occupying state r at
# time u is given by the (r,s) entry of the matrix P(t)=exp(tQ).

inv<-15
mc1<-mc$reproMC1()
par(mfrow=c(2,1))
states<-1:mc1$Stepno

P<-mc$calcTransMat(mc1$rate,t=inv,method="Matrix")
plot(states,P[1,],type="l",xlab="state",ylab="probability")      # black line pr given state 1 at time u
lines(states,P[150,],type="l",col="red")  # red line pr given state 150 at time u

# try another method
P<-mc$calcTransMat(mc1$rate,t=inv,method="robust")
plot(states,P[1,],type="l",xlab="state",ylab="probability")      # black line pr given state 1 at time u
lines(states,P[150,],type="l",col="red")  # red line pr given state 150 at time u
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
## Prediction of ovulation time
mc1<-mc$reproMC1()
states<-1:mc1$Stepno
times<-states*mc1$beta
meanProgest<-mc$meanProgest(1:mc1$Stepno,alpha=3)
stdProgest  <- 0.15     # progest_{state} ~ N(meanProgest_{state},stdProgest^2)
progest<-rnorm(length(states),meanProgest,stdProgest)
progestPerDay <- 1       # progest measurements per day



# at time 5
library(nor1Mix)
progestState<-norMix(mu=meanProgest,sig2=rep(stdProgest,mc1$Stepno))    # Mixture density of progesteron measurements







trans<-mc$calcTransMat(mc1$rate,t=1/progestPerDay)


res[[1]]<-simDisMC(oneDayTrans,70,mc1$kOestrousState)


prog0<- min(resprog[[1]][5:30])
prog1<- max(resprog[[1]][5:30])



for (i in (5:30)){
  #i<-5
  w<-posteriorMix(ProgState,resprog[[1]][i])
  FirstPassage<-FirstPassageDist(rate,w,kOestrousState,50,1/12)


  plot(i+FirstPassage$time,FirstPassage$d,type="l",
		  xlim=c(0,50),xlab="Day no",ylab="Probability",
		  ylim=c(0,0.032))
  points(c(5:i), (resprog[[1]][5:i]-prog0)*0.03/(prog1-prog0),col="red",pch=19)
  text(40,0.030,paste("Obs day =",i))

  neww<- w %*% OneDayTrans
  ProgState[,3]<-as.matrix(neww)
}




x=c(0.3,0.7,0.5)


#str(ProgState)



 # 1-2 (CL regression to Oestrogen Peak) 53
 # 2-3 (Oestrogen peak to LH peak 4.2)    4
 # 3-4 (LH peak to Ovulation)            18
 # 4-1 (Ovulation to CL regression)     146

  tider<-c(5,0.4,1.7)
  tider<-c(tider,21-sum(tider))

 alpha<- tider/0.095
 steps<-round(alpha)


# First-passage time
absorbstate=22
timestep<-1/24
Interval<-30


  prior<- ProgState[,3]

  dist  <-FirstPassageDist(rate,prior,absorbstate,Interval,timestep)$d

 mynewportraitplot("FirstPassage")
  plot(c(1:round(Interval/timestep))*timestep,dist,type='l',xlab="Time to oestrus, days",ylab="probability")
 dummy<-dev.off()
