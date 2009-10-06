###########################################################################/**
# @RdocClass CowYield
#
# @title "CowYield class"
#
# \description{
#  Containing all biological functions related to the yield of the cow.
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{ssmCov}{ A matrix where row i contain the parameters (sigmaA, sigmaX, rho, sigmaNu) for lactation i used by the SSM. }
#   \item{meanYield}{ A matrix of size (lactations x dfc's) containing the mean yield values
#       (lactation curves). If e.g. the number of rows is 3 it is assumed that lactation 3, 4,
#       .. have the same lactation curve as lactation 3. If NULL default loaded from the data folder}
#   \item{milkDensity}{ Weight in kg of 1 liter milk. }
#   \item{fat}{ Percentage of fat in 1 liter of milk. }
#   \item{protein}{ Percentage of protein in 1 liter of milk. }
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods ""
# }
#
# \section{To do}{
#   Add functions for Wilmink and Wood.
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("CowYield", function(
	ssmCov=matrix(c(4.9877614,4.63,4.08, 5.6342848, 4.7387595, 4.7387595, 0.979716727, 0.979716727, 0.979716727, 6.3565842, 6.3565842, 6.3565842),nrow=3),
	meanYield=NULL,
	milkDensity=1.033, fat=4.12, protein=3.38, ...)
{
	if (is.null(meanYield)) {
		#data(meanYieldMat)
		#meanYield<-as.matrix(meanYieldMat)
		fileN<-system.file("data", "meanYieldMat.csv", package="dairy")
		meanYield<-read.table(fileN, header = TRUE, sep = ";")
		colnames(meanYield)<-NULL
		meanYield<-as.matrix(meanYield)
	}
	# Create the SSM model for daily yields (list used by dlm package)
	buildSSM<-function(lac) {
		if (nrow(ssmCov)<lac) return;
		x<-rep(0,4)
		x[1]<-ssmCov[lac,3]    # rho
		x[2]<-ssmCov[lac,2]    # sigmaX
		x[3]<-ssmCov[lac,4]    # sigmaNu
		x[4]<-ssmCov[lac,1]    # sigmaA
		GG<-diag(c(1,x[1]))
		V = matrix(x[3]^2,1,1)
		W<-matrix(c(0,0,0,(1-x[1]^2)*x[2]^2),2,2)
		return(list(
		m0 = rep(0,2),
		C0 = diag(c(x[4]^2,x[2]^2)),
		FF = matrix(c(1,1),1,2),
		GG = GG,
		V = V,
		W = W))
	}

	colnames(ssmCov)<-c("sigmaA", "sigmaX", "rho", "sigmaNu")
	ssmModels<-list()    # list of ssm models
	for (i in 1:nrow(ssmCov)) ssmModels[[i]] <- buildSSM(i)

	extend(Object(), "CowYield",
		.ssmCov=ssmCov,
		.milkDensity=milkDensity,
		.fat=fat,
		.protein=protein,
		.ssmModels=ssmModels,
		meanYieldMat=meanYield
	)
})


#########################################################################/**
# @RdocMethod getSsmModel
#
# @title "Return the SSM model for a specific lactation"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{lac}{ Lactation number. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("getSsmModel", "CowYield", function(this, lac, ...) {
	if (lac>nrow(this$.ssmCov)) return(this$.ssmModels[[nrow(this$.ssmCov)]])
	return(this$.ssmModels[[lac]])
})




#########################################################################/**
# @RdocMethod meanYield
#
# @title "Mean daily milk yield in kg"
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
#   The mean yield in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod dailyYield,
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("meanYield", "CowYield", function(this, lac, dfc, ...) {
	if (lac>nrow(this$meanYieldMat)) return(this$meanYieldMat[nrow(this$meanYieldMat),dfc])
		else return(this$meanYieldMat[lac,dfc])
})



#########################################################################/**
# @RdocMethod dailyYield
#
# @title "Total daily milk yield in kg"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{m}{ Conditional mean of the latent vector (A,X) in the SSM (2x1 matrix) }
#  \item{lac}{ Lactation number. }
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The yield in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod meanYield,
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("dailyYield", "CowYield", function(this, m, lac, dfc, ...) {
	ave<-this$meanYield(lac,dfc)
	ssmModel<-this$getSsmModel(lac)
	y<-ssmModel$FF %*% ssmModel$GG %*% m
	return(as.numeric(ave+y))
})


#########################################################################/**
# @RdocMethod dailyYieldECM
#
# @title "Total daily milk yield in kg ECM"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{m}{ Conditional mean of the latent vector (A,X) in the SSM (2x1 matrix) }
#  \item{lac}{ Lactation number. }
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The yield in kg.
# }
#
# @author
#
# \seealso{
#   @seemethod meanYield,
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("dailyYieldECM", "CowYield", function(this, m, lac, dfc, ...) {
	ave<-this$meanYield(lac,dfc)
	ssmModel<-this$getSsmModel(lac)
	y<-ssmModel$FF %*% ssmModel$GG %*% m
	milk <- as.numeric(ave+y)
	return(milk* ((383 * this$.fat + 242 * this$.protein + 783.2)/3140))
})


#########################################################################/**
# @RdocMethod kg2KgECM
#
# @title "Convert the milk yield in kg to energy corrected milk in kg"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{milk}{ Milk yield in kg. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The yield in kg ECM.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("kg2KgECM", "CowYield", function(this, milk, ...) { # fat, protein and lactose are in percentage
	return(milk* ((383 * this$.fat + 242 * this$.protein + 783.2)/3140))
})


#########################################################################/**
# @RdocMethod liters2KgECM
#
# @title "Convert the milk yield in liters to energy corrected milk in kg"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{milk}{ Milk yield in liters. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The yield in kg ECM.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/CowYield.Rex"
#
#*/#########################################################################
setMethodS3("liters2KgECM", "CowYield", function(this, milk, ...) { # fat, protein and lactose are in percentage (lactose=4.61)
	#return(milk*density * ((383 * fat + 242 * protein + 163.2 * lactose)/3140))
	return(milk*this$.milkDensity * ((383 * this$.fat + 242 * this$.protein + 783.2)/3140))
})
