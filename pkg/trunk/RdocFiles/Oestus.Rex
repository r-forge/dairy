repro<-Oestus()
repro$getFields()

## make a plot of marginal pr for pos preg test
dat<-NULL
for (lac in 1:3) {
	for (dfc in 1:300) dat<-rbind(dat,c(lac,dfc,repro$probPregTest(lac,dfc)))
}
colnames(dat)<-c("lac","dfc","y")
dat<-as.data.frame(dat)
xyplot(y~dfc,
	data=dat,
	group=lac,
	type=c("l"),
	ylab="pr",
	xlab="dfc",
	auto.key=list(lines=TRUE,points=FALSE,columns=3),
)


## some plots of log(progesterone+1) curves
dat<-NULL
for (b1 in seq(-0.5,1,by=0.1)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),b1=b1,sigma=0), g=paste("b1 =",b1)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (b2 in seq(-1,1,by=0.2)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),b2=b2,sigma=0), g=paste("b2 =",b2)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (l in seq(20,30,by=4)) {
	dat<-rbind(dat,data.frame(t=seq(0,50,by=0.25), y=repro$logProgst(seq(0,50,by=0.25),lambda=l,sigma=0), g=paste("lambda =",l)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (g in seq(-1,1,by=0.5)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),gamma=g,sigma=0), g=paste("gamma =",g)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (p in seq(1,3,by=0.5)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),psi1=p,sigma=0), g=paste("psi1 =",p)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (p in seq(-1,1,by=0.2)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),psi2=p,sigma=0), g=paste("psi2 =",p)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (p in seq(1,3,by=0.5)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),psi1=p,sigma=0, psi2=-0.3), g=paste("psi1 =",p)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))

dat<-NULL
for (p in seq(-1,0,by=0.2)) {
	dat<-rbind(dat,data.frame(t=seq(0,45,by=0.25), y=repro$logProgst(seq(0,45,by=0.25),psi2=p,sigma=0, psi1=1.5, lambda=21), g=paste("psi2 =",p)))
}
xyplot(y~t, data=dat, group=g, type="l",auto.key=list(lines=T,points=F,columns=length(unique(dat$g))))
## End plots


## Plots of activity. Time 0 is where the oestrus cycle starts
tHigh=c(0,1.5,c(0,1.5)+21,c(0,1.5)+21*2)
tLow=c(1.55,c(-0.05,1.55)+21,c(-0.05,1.55)+21*2,45)
dat<-repro$activity(tLow,tHigh,n=10)
plot(act~t,data=dat,type="l", frame=F, ylab="activity (simulated)", main="activity", lab=c(10,5,7))
abline(v=c(0,21,42),col="gray")


## Create a plot of the oestrus cycle
# mar=c(bottom, left, top, right)
par(mar=c(3,4,2,0),oma=c(0,0,0,0),cex=0.75)
layout(matrix(c(1,2,3), ncol=1),heights=c(17,17,17))
xlim<-c(0,45)
# progesterone
dat<-data.frame(t=seq(0,45,by=0.25), progst=repro$logProgst(seq(0,45,by=0.25),psi2=-0.3,psi1=2.3,lambda=21,gamma=0))
plot(progst~t, type="l", data=dat, xlim=xlim, frame=F, lab=c(10,10,7), main="progesterone",
	ylab="log(progest+1) (simulated)",xlab="time")
abline(v=c(0,21,42),col="gray")
# activity
par(mar=c(3,4,1,0))
tHigh=c(0,1.5,c(0,1.5)+21,c(0,1.5)+21*2)
tLow=c(1.55,c(-0.05,1.55)+21,c(-0.05,1.55)+21*2,45)
dat<-repro$activity(tLow,tHigh,n=10)
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
## End plot
