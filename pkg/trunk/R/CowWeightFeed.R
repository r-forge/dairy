###########################################################################/**
# @RdocClass CowWeightFeed
#
# @title "CowWeightFeed class"
#
# \description{
#  Containing all biological functions related to the weight of the cow and feeding. The default values are taken from an average Danish HF herd.
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{udvegt}{ Parameter in the Gompertz weight curve used in method \code{bodyWStd} (full-grown weight). }
#   \item{gompM}{  Parameter in the Gompertz weight curve used in method \code{bodyWStd} (horisontal shift). }
#   \item{gompN}{  Parameter in the Gompertz weight curve used in method \code{bodyWStd} (growth rate). }
#   \item{foster}{ Parameter in the method \code{fosterW} (growth rate). }
#   \item{fostervaegt}{ Parameter in the method \code{fosterW} (coefficient that multiplied with \code{udvegt} gives the fetus weight at birth). }
#   \item{L0}{ 1. parameter used in method \code{BCS}. }
#   \item{T}{ 2. parameter used in method \code{BCS}. }
#   \item{LT}{ 3. parameter used in method \code{BCS}. }
#   \item{Lnext}{ 4. parameter used in method \code{BCS}. }
#   \item{maxLLoss}{ 5. parameter used in method \code{BCS}. }
#   \item{kgBCSUnitPrStdBW}{ Additional weight of one unit of BCS above the mean 3. }
#   \item{sfuFetus}{ SFU (Scandinavian feed unit) need for one kg of fetus weight. }
#   \item{sfuECM}{ SFU need for one kg of energy corrected milk.  }
#   \item{sfuStdBWG}{ SFU need for one kg of std body weight.  }
#   \item{feKgBcsGainA}{ 1. parameter in SFU need for an increase in BCS. }
#   \item{feKgBcsGainB}{ 2. parameter in SFU need for an increase in BCS. }
#   \item{feKgBcsLossA}{ 1. parameter in SFU need for an loss in BCS. }
#   \item{feKgBcsLossB}{ 2. parameter in SFU need for an loss in BCS. }
#   \item{sfuDry}{ SFU need for one day in the dry period. }
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods ""
# }
#
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("CowWeightFeed", function(udvegt=680, gompM=2.5483,
	gompN=0.00314, foster=0.02, fostervaegt=0.1133, L0=3.0, T=70,
	LT=2.6, Lnext=3.0, maxLLoss=0.0324, kgBCSUnitPrStdBW=0.090,
	sfuFetus=0.03647, sfuECM=0.4, sfuStdBWG=4.0, feKgBcsGainA=0.48570,
	feKgBcsGainB=1.38570, feKgBcsLossA=0.47140, feKgBcsLossB=1.08570,
	sfuDry=7, ...)
{
	extend(Object(), "CowWeightFeed",
		.udvegt=udvegt,
		.gompM=gompM,
		.gompN=gompN,
		.foster=foster,
		.fostervaegt=fostervaegt,
		.L0=L0,
		.T=T,
		.LT=LT,
		.Lnext=Lnext,
		.maxLLoss=maxLLoss,
		.kgBCSUnitPrStdBW=kgBCSUnitPrStdBW,
		.sfuFetus=sfuFetus,
		.sfuECM=sfuECM,
		.sfuStdBWG=sfuStdBWG,
		.feKgBcsGainA=feKgBcsGainA,
		.feKgBcsGainB=feKgBcsGainB,
		.feKgBcsLossA=feKgBcsLossA,
		.feKgBcsLossB=feKgBcsLossB,
		.sfuDry=sfuDry
	)
})


#########################################################################/**
# @RdocMethod BCS
#
# @title "Body condition score"
#
# \description{
#   Body condition score (1-5 scale) of a cow during lactation.
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving.}
#  \item{con}{ Day of conception.}
#  \item{gest}{ Expected length of pregnancy. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The BCS (numeric).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{
# The BCS model is based on N.C. Friggens, K.L. Ingvartsen, and G.C. Emmans.
# \emph{Prediction of body lipid change in pregnancy and lactation}.
# Journal of Dairy Science, 87:988-1000, 2004. }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("BCS", "CowWeightFeed", function(this, dfc, con, gest=282, ...) {
	TT<-this$.T
	if (con>=TT) {
		b<-2*(this$.LT-this$.L0)/TT
		if (b< -this$.maxLLoss) {
			b<- -this$.maxLLoss
			TT<-2*(this$.LT-this$.L0)/-this$.maxLLoss
		}
		if (dfc<TT) return((this$.L0-this$.LT)/TT^2*dfc^2+2*(this$.LT-this$.L0)/TT*dfc+this$.L0)
		if (dfc<=con) return(this$.LT)
		return((this$.Lnext-this$.LT)/gest^2*(dfc^2-2*con*dfc+con^2)+this$.LT)
	} else {
		dLT<-2*(TT-con)*(this$.Lnext-this$.LT)/gest^2
		LTarget<-this$.LT+(dLT*(TT-con))/2
		dL0<-2*(LTarget-this$.L0)/TT-dLT
		if (dL0< -this$.maxLLoss) {
			dL0<- -this$.maxLLoss
			TT<-2*(LTarget-this$.L0)/(dL0+dLT)
		}
		if (dfc<TT) {
			a<-(dLT-dL0)/TT
			b<-dL0
			c<-this$.L0
		} else {
			dLnext<-2*(this$.Lnext-LTarget)/gest
			a<-(dLnext-dLT)/(con+gest-TT)
			b<-dLT-a*TT
			c<- LTarget - (a/2*TT^2+b*TT)
		}
		return(0.5*a*dfc^2+b*dfc+c)
	}
})



#########################################################################/**
# @RdocMethod deltaBCS
#
# @title "Body condition score change"
#
# \description{
#   Body condition score (1-5 scale) change.
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving.}
#  \item{con}{ Day of conception.}
#  \item{gest}{ Expected length of pregnancy. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The BCS change (numeric).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{
# The BCS model is based on N.C. Friggens, K.L. Ingvartsen, and G.C. Emmans.
# \emph{Prediction of body lipid change in pregnancy and lactation}.
# Journal of Dairy Science, 87:988-1000, 2004. }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("deltaBCS", "CowWeightFeed", function(this, dfc, con, gest=282, ...) {
	T<-this$.T
	LT<-this$.LT
	L0<-this$.L0
	maxLLoss<-this$.maxLLoss
	Lnext<-this$.Lnext
	if (con>=T) {
		b<-2*(LT-L0)/T
		if (b< -maxLLoss) {
			b<- -maxLLoss
			T<-2*(LT-L0)/-maxLLoss
		}
		if (dfc<T) return(b*(-dfc/T+1))
		if (dfc<=con) return(0)
		return(2*(Lnext-LT)/gest^2*(dfc-con))
	} else {
		dLT<-2*(T-con)*(Lnext-LT)/gest^2
		LTarget<-LT+(dLT*(T-con))/2
		dL0<-2*(LTarget-L0)/T-dLT
		if (dL0< -maxLLoss) {
			dL0<- -maxLLoss
			T<-2*(LTarget-L0)/(dL0+dLT)
		}
		if (dfc<T) {
			a<-(dLT-dL0)/T
			b<-dL0

		} else {
			dLnext<-2*(Lnext-LTarget)/gest
			a<-(dLnext-dLT)/(con+gest-T)
			b<-dLT-a*T
		}
		return(a*dfc+b)
	}
})


#########################################################################/**
# @RdocMethod bodyW
#
# @title "Body weight of the cow (fetus not included)"
#
# \description{
#   @get "title". The body weigth follows a Gompertz weight curve.
# }
#
# @synopsis
#
# \arguments{
#  \item{age}{ Age in days at the start of the lactation. }
#  \item{dfc}{ Days from calving. }
#  \item{con}{ Day of conception. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The body weight (numeric).
# }
#
# @author
#
# \seealso{
#   @seemethod bodyWStd
#   @seeclass
# }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("bodyW", "CowWeightFeed", function(this, age, dfc, con, ...) {
	return(this$bodyWStd(age+dfc)+this$weightBCS(age+dfc)*(this$BCS(dfc,con)-3))
})


#########################################################################/**
# @RdocMethod bodyWStd
#
# @title "Standardized body weight"
#
# \description{
#   Gompertz curve of body weight at body condition score 3 (1-5 scale). The fetus not included.
# }
#
# @synopsis
#
# \arguments{
#  \item{age}{ Age in days from birth. }
#  \item{...}{Not used.}
# }
#
# \value{
#   Standardized body weight in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod bodyW
#   @seeclass
# }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("bodyWStd", "CowWeightFeed", function(this, age, ...) {
	x<-this$.udvegt*exp(-this$.gompM*exp(-this$.gompN*age))
	return(x)
})


#########################################################################/**
# @RdocMethod deltaBodyWStd
#
# @title "Daily standardized body weight change"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{age}{ Age in days from birth. }
#  \item{...}{Not used.}
# }
#
# \value{
#   Daily standardized body weight change in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod bodyW
#   @seemethod bodyWStd
#   @seeclass
# }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("deltaBodyWStd", "CowWeightFeed", function(this, age, ...) {
	return(-this$.gompN*this$bodyWStd(age)*log(this$bodyWStd(age)/this$.udvegt))
})

#########################################################################/**
# @RdocMethod fetusW
#
# @title "Weight of the fetus in the cow"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{dtc}{ Days to calving. }
#  \item{...}{Not used.}
# }
#
# \value{
#   Fetus weight in kg.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("fetusW", "CowWeightFeed", function(this, dtc, ...) { # dtc = days to calving
	return(exp(-this$.foster*dtc)*this$.fostervaegt*this$.udvegt)
})

#########################################################################/**
# @RdocMethod weightBCS
#
# @title "Weight of BCS"
#
# \description{
#   The total weight of the body lipid for a standardized cow with BCS = 3.
# }
#
# @synopsis
#
# \arguments{
#  \item{age}{ Age in days from birth. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The weight in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod BCS
#   @seeclass
# }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("weightBCS", "CowWeightFeed", function(this, age, ...)
	return(this$bodyWStd(age)*this$.kgBCSUnitPrStdBW)
)


#########################################################################/**
# @RdocMethod energy
#
# @title "Daily energy requirements in SFU units"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{age}{ Age at calving (dfc=0).  }
#  \item{dfc}{ Days from calving. }
#  \item{con}{ Day of conception (in dfc). }
#  \item{milk}{ Kg ECM milk produced. }
#  \item{gest}{ Expected length of pregnancy. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A vector containing: energy for maintaince, energy for fetus, energy for milk,
#   energy for body weight gain, energy for BCS, total energy (sum of previous).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{
# The BCS model is based on N.C. Friggens, K.L. Ingvartsen, and G.C. Emmans.
# \emph{Prediction of body lipid change in pregnancy and lactation}.
# Journal of Dairy Science, 87:988-1000, 2004. }
#
# @examples "../RdocFiles/CowWeightFeed.Rex"
#
#*/#########################################################################
setMethodS3("energy", "CowWeightFeed", function(this, age, dfc, con, milk, gest, ...) {
	eMain <- 1.1*(this$bodyWStd(age+dfc)/200 + 1.5)
	if (dfc>con) eFetus <- this$fetusW((gest+con)-dfc)*this$.sfuFetus else eFetus <- 0
	eECM <- milk*this$.sfuECM
	eBWG <-  this$.sfuStdBWG*(this$bodyWStd(age+dfc+1)-this$bodyWStd(age+dfc))                 #this$.sfuStdBWG*this$deltaBodyWStd(age+dfc)
	dKgBCS <- (this$BCS(dfc+1,con)-this$BCS(dfc,con))*this$weightBCS(age+dfc)                  #this$deltaBCS(dfc,con)*this$weightBCS(age+dfc)
	bcs <- this$BCS(dfc,con)
	if (dKgBCS<0) eBCS <- dKgBCS * (this$.feKgBcsLossA * bcs + this$.feKgBcsLossB) else eBCS <- dKgBCS * (this$.feKgBcsGainA * bcs + this$.feKgBcsGainB)
	#cat("main:",eMain," fetus:",eFetus," ECM:",eECM," BWG:",eBWG," BCS:",eBCS,"\n",sep="")
	return(c(eMain,eFetus,eECM,eBWG,eBCS,eMain+eFetus+eECM+eBWG+eBCS))
})
