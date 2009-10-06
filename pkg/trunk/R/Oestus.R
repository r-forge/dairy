###########################################################################/**
# @RdocClass Oestus
#
# @title "Oestus class"
#
# \description{
#  Containing biological functions related the oestus of the cow.
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{gestLth}{ Length of pregnancy. }
#   \item{insemStart}{ The day the manager starts insemination (given in heat). }
#   \item{insemFinish}{ The day the manager stop insemination (given in heat and not pregnant). }
#   \item{pregTestLth}{ Days after insemination a pregnancy test is taken. }
#   \item{dryPeriodLth}{ Length of dry period. }
#   \item{probCalf}{ Probability of a female calf. }
#   \item{probBullCalf}{ Probability of a male calf. }
#   \item{probPregTest}{ Matrix of positive pregnancy test probabilities (parities x dfcs). If e.g. the number of rows is 3 it is assumed that lactation 3, 4,
#       .. have the same lactation probability as lactation 3. If NULL default loaded from the data folder}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods ""
# }
#
# @examples "../RdocFiles/Oestus.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("Oestus", function(gestLth=282, insemStart=35, insemFinish=250,
	pregTestLth=40, dryPeriodLth=49, probCalf=0.5, probBullCalf=0.5,
	probPregTest=NULL, ...)
{
	if (is.null(probPregTest)) {
		fileN<-system.file("data", "probPregTestMat.csv", package="dairy")
		probPregTest<-read.table(fileN, header = TRUE, sep = ";")
		colnames(probPregTest)<-NULL
		probPregTest<-as.matrix(probPregTest)
	}

	extend(Object(), "Oestus",
		gestLth=gestLth,
		insemStart=insemStart,
		insemFinish=insemFinish,
		pregTestLth=pregTestLth,
		dryPeriodLth=dryPeriodLth,
		probCalf=probCalf,
		probBullCalf=probBullCalf,
		.probPregTestMat=probPregTest
	)
})


#########################################################################/**
# @RdocMethod probPregTest
#
# @title "Marginal probability of positive pregnancy test "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{lac}{ Lactation number. }
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The probabilities.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Oestus.Rex"
#
#*/#########################################################################
setMethodS3("probPregTest", "Oestus", function(this, lac, dfc, ...) {
	if (!is.null(this$.probPregTestMat)) { # we just make a lookup in a big matrix
		if (lac>nrow(this$.probPregTestMat)) return(this$.probPregTestMat[nrow(this$.probPregTestMat),dfc])
		else return(this$.probPregTestMat[lac,dfc])
	}
	stop("probPregTestMat matrix not calculated yet!")
})


#########################################################################/**
# @RdocMethod logProgst
#
# @title "Simple simulation of log(progesterone+1) based on Højsgaard et al 2009"
#
# \description{
#   Simple curve of log(progesterone+1) (to keep variance constant) based on Højsgaard et al 2009.
#   The function use a cyclic curve based on cosinus.
# }
#
# @synopsis
#
# \arguments{
#  \item{t}{ Time in days. }
#  \item{b1}{ Parameter. }
#  \item{b2}{ Parameter. }
#  \item{lambda}{ Parameter. }
#  \item{gamma}{ Parameter. }
#  \item{psi1}{ Parameter. }
#  \item{psi2}{ Parameter. }
#  \item{sigma}{ Std of the random noise. }
#  \item{tshift}{ Time shift. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{[1] Højsgaard, S.; Green, P.J and Friggens N.C. A cyclical state space model for prediction of progesterone concentration in milk from dairy cows. Dept. of Genetics and Biotechnology, Aarhus University, 2009.}
#
# @examples "../RdocFiles/Oestus.Rex"
#
#*/#########################################################################
# create a curve of log(progesterone+1)  based on Højsgaard et al 2009
setMethodS3("logProgst", "Oestus", function(this,t,b1=0.85,b2=-0.56,lambda=22.44,
	gamma=0.70,psi1=2.9,psi2=-0.7,sigma=0.05,tshift=0,...)
{
	t<-t+tshift  # make a shift in time
	h <- t*2*pi/lambda-psi1+(exp(gamma)/(1+exp(gamma)))*cos(t*2*pi/lambda - psi2)
	return(exp(b1) + exp(b2)*sin(h)+rnorm(length(t),0,sigma))
})


#########################################################################/**
# @RdocMethod activity
#
# @title "Simple curve af activity counts"
#
# \description{
#   Simple simulation of activity counts. Make a simple curve with some peaks and add some noise.
# }
#
# @synopsis
#
# \arguments{
#  \item{tLow}{ Time intervals with low activity. }
#  \item{tHigh}{ Time intervals with high activity. }
#  \item{...}{Parameter n (number of samples in the tHigh intervals) passed to \code{\link{approx}} function.}
# }
#
# \value{
#   The probabilities.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{[1] Højsgaard, S.; Green, P.J and Friggens N.C. A cyclical state space model for prediction of progesterone concentration in milk from dairy cows. Dept. of Genetics and Biotechnology, Aarhus University, 2009.}
#
# @examples "../RdocFiles/Oestus.Rex"
#
#*/#########################################################################
setMethodS3("activity", "Oestus", function(this,tLow,tHigh,sigma=0.7,...) {
	dat<-NULL
	for (i in 1:(length(tHigh)/2)) {
		start<-tHigh[2*i-1]
		stop<-tHigh[2*i]
		x<-c(start,(stop-start)/2+start,stop)
		y<-c(20,runif(1,20,45),20)
		tmp<-approx(x,y,...)
		dat<-rbind(dat,data.frame(t=tmp$x,act=tmp$y))
	}
	for (i in 1:(length(tLow)/2)) {
		start<-tLow[2*i-1]
		stop<-tLow[2*i]
		x<-seq(start,stop,by=0.5)
		y<-rep(20,length(x))
		dat<-rbind(dat,data.frame(t=x,act=y))
	}
	dat<-dat[order(dat$t),]
	res<-rnorm(nrow(dat),0,sigma)
	dat$act<-dat$act+res   #runif(nrow(dat),dat$act-0.05*dat$act,dat$act+0.05*dat$act)
	return(dat)
})
