y<-CowYield()
y$getFields()
y$meanYieldMat[1:3,1:5]     # lactation 1-3, dfc 1-5

lac<-2
dfc<-100
m<-matrix(c(2,4))

y$meanYield(lac, dfc)    # simple lookup in meanYieldMat
y$dailyYield(m, lac, dfc)
y$dailyYieldECM(m, lac, dfc)


# plot average milk yield curves
dat<-NULL
for (lac in 1:3) {
	for (dfc in 1:400) dat<-rbind(dat,c(lac=lac,dfc=dfc,y=y$meanYield(lac,dfc)))
}
dat<-as.data.frame(dat)
xyplot(y~dfc,
	data=dat,
	group=lac,
	type=c("l"),
	ylab="mean yield",
	xlab="dfc",
	auto.key=list(lines=TRUE,points=FALSE,columns=3),
)

# plot milk yield curve for cow with m = (2,-1) in lac 2
m<-matrix(c(2,-1))
lac<-2
dat<-NULL
for (dfc in 1:400) dat<-rbind(dat,c(lac=lac,dfc=dfc,y=y$meanYield(lac,dfc),g=1))
for (dfc in 1:400) dat<-rbind(dat,c(lac=lac,dfc=dfc,y=y$dailyYield(m,lac,dfc),g=2))
dat<-as.data.frame(dat)
dat$g<-ifelse(dat$g==1,"ave","cow")
xyplot(y~dfc,
	data=dat,
	group=g,
	type=c("l"),
	ylab="mean yield",
	xlab="dfc",
	auto.key=list(lines=TRUE,points=FALSE,columns=2),
)

# convert methods
y$kg2KgECM(10)
y$liters2KgECM(10)
