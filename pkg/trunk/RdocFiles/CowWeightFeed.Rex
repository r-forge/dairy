w<-CowWeightFeed()

## plot BCS
dat<-NULL
for (con in c(30,100,150,200)) {
	for (i in 0:(con+282)) {
		dat<-rbind(dat,c(dfc=i,con=con,BCS=w$BCS(i,con)))
	}
}
dat<-data.frame(dat)

xyplot(BCS~dfc,
	group=con,
	data=dat,
	type=c("l"),
	ylab="BCS",
	xlab="dfc",
	auto.key=list(lines=TRUE,points=FALSE,columns=4),
)


## Plot Gompertz weight curve
age<-0:(8*365)
xyplot(w$bodyWStd(age)~age,
	type=c("l"),
	ylab="Standardized Body weight",
	xlab="age"
)


## plot body weight change
xyplot(w$deltaBodyWStd(age)~age,
	type=c("l"),
	ylab="Body weight change",
	xlab="age"
)


## Plot body weight
dat<-NULL
for (con in c(30,100,150,200)) {
	for (age in c(400,450,500,550)) {
		for (dfc in 0:(con+282)) {
			dat<-rbind(dat,c(age=age,con=con,dfc=dfc,BW=w$bodyW(age,dfc,con)))
		}
	}
}
dat<-data.frame(dat)

# plots of BW for different start ages
con1=100
age1=500
xyplot(BW~dfc,
	data=subset(dat,con==con1),
	group=age,
	type=c("l"),
	ylab="Body weight",
	xlab="dfc given age",
	auto.key=list(lines=TRUE,points=FALSE,columns=4),
)

# plots of BW for different conception days
xyplot(BW~dfc,
	data=subset(dat,age==age1),
	group=con,
	type=c("l"),
	ylab="Body weight",
	xlab="dfc given day of conception",
	auto.key=list(lines=TRUE,points=FALSE,columns=4),
)


## Plot fetus weight
t<-0:282
xyplot(w$fetusW(t)~t,
	type=c("l"),
	ylab="Fetus weight",
	xlab="days to calving"
)


## Plot weight of BCS
plot(1:2000,w$weightBCS(1:2000),xlab="Age in days",ylab="Total weight of BCS in kg",type="l")


## Use Wilmink function to describe milk yield
wilmink<-function(t,a=20.7,b=-0.041,c=-4.00,d=0.0000442){
  a+b*t+c*exp(-0.05*t)+d*t^2
}

## Plots of energy for different parameter values
dat<-NULL
gest=282
con=100
for (age in c(400,425,450,475,500)) {
	for (dfc in 1:(con+gest)) {
		dat<-rbind(dat,c(age=age,con=con,dfc=dfc,w$energy(age,dfc,con,wilmink(dfc))))
	}
}
dat<-data.frame(dat)

dat1<-subset(dat,age==400)
par(mfrow=c(3,2))
plot(dat1$V4~dat1$dfc,type="l",ylab="maintenance",xlab="dfc")
plot(dat1$V5~dat1$dfc,type="l",ylab="fetus",xlab="dfc")
plot(dat1$V6~dat1$dfc,type="l",ylab="milk",xlab="dfc")
plot(dat1$V7~dat1$dfc,type="l",ylab="std. BW gain",xlab="dfc")
plot(dat1$V8~dat1$dfc,type="l",ylab="BCS change",xlab="dfc")
plot(dat1$V9~dat1$dfc,type="l",ylab="total",xlab="dfc")
title("energy given conception at day 100 and start age 400",outer=TRUE)
