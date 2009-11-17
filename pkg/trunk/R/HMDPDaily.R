###########################################################################/**
# @RdocClass HMDPDaily
#
# @title "HMDP model based on daily yields"
#
# \description{
#  Class of the model described in [1] which is a 3-level hierarchial Markov decision process with an SSM embedded based on daily milk yields.
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{ssmLac}{ The lactation where the mean and variance of the SSM model is used to discretize. }
#   \item{maxLac}{ Maximum lactation number in the HMDP. }
#   \item{priceECM}{ Price of one energy corrected milk. }
#   \item{priceSFU}{ Price of one scandinavian feed unit. }
#   \item{priceHeifer}{ Price of a heifer. }
#   \item{priceCalf}{ Price of a calf. }
#   \item{priceBullCalf}{ Price of a bull calf. }
#   \item{priceCarcassKg}{ Price of one kg of the carcass. }
#   \item{oCowY}{ Object of the CowYield class. }
#   \item{oCowWF}{ Object of the CowWeightFeed class. }
#   \item{oDis}{ Object of the Discretize class. }
#   \item{oOestus}{ Object of the Oestus class. }
#   \item{csvPrM}{ File name of the csv containing probabilities of the conditional latent mean (if already generated). }
#   \item{cowAgeParam}{ Parameters of the linear function describing cow age given lactation and dfc. }
#   \item{...}{Additional parameters passed to \code{oDis$discretize2DNonunif}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods ""
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
# \references{[1] Nielsen, L.R.; Jørgensen, E.; Kristensen, A.R. & Østergaard, S. Optimal Replacement Policies for Dairy Cows Based on Daily Yield Measurements Dept. of Genetics and Biotechnology, Aarhus University, 2009.}
#
# @author
#*/###########################################################################
setConstructorS3("HMDPDaily", function(ssmLac=1, maxLac=10, priceECM=2.8,
	priceSFU=1.37, priceHeifer=8600, priceCalf=850, priceBullCalf=950,
	priceCarcassKg=14.5, oCowY = CowYield(), oCowWF=CowWeightFeed(),
	oDis=Discretize(), oOestus=Oestus(), csvPrM=NULL,
	cowAgeParam=c(493.1604,379.3949), ...)
{
	ssmModel<-oCowY$.ssmModels[[ssmLac]]
	# state variable sM
	sM<-oDis$discretize2DNonunif(as.numeric(ssmModel$m0), ssmModel$C0, modifyCenter="split2", ...)
	sizeM<-length(sM)
	# state variable sM1
	m1<-sort(unique(unlist(lapply(sM, function(x) x$center[1]))))   # number of different center points of A
	sA<-oDis$discretize1DVec(m1)
	sA<-as.matrix(sA,rownames.force = FALSE)
	sizeA<-nrow(sA)
	# state variable sDry
	sDry<-c(ceiling((oOestus$insemStart+oOestus$gestLth-oOestus$dryPeriodLth)/7):ceiling((oOestus$insemFinish+oOestus$gestLth-oOestus$dryPeriodLth)/7),-1)      # possible weeks where dry the cow (-1 indicate unknown) always dry at day dry*dryInv. Note -1 put last in vector!
	sDry<-cbind(idxDry=1:length(sDry)-1,sDry)
	sizeDry<-nrow(sDry)
	sv<-list(sM=sM,sizeM=sizeM,sA=sA,sizeA=sizeA,sDry=sDry,sizeDry=sizeDry)

	prM<-NULL
	if (!is.null(csvPrM)) {
		prMat<-as.matrix(read.csv(csvPrM,header=TRUE))
		prMat<-cbind(prMat,NA)    # to have a NA in each row
		colnames(prMat)<-NULL
		prM<-array(list(0),c(max(prMat[,1])+1,max(prMat[,2])+1,max(prMat[,3])+1))
		add<-function(x) {
			prM[[x[1]+1,x[2]+1,x[3]+1]]<<-x[4:(which(is.na(x))[1]-1)]
			invisible(NULL)
		}
		apply(prMat,1,add)
	}

	extend(Object(), "HMDPDaily",
		stateVar=sv,
		maxLac=maxLac,
		maxDfc=oOestus$insemFinish+oOestus$gestLth-oOestus$dryPeriodLth,
		priceECM=priceECM,
		priceSFU=priceSFU,
		priceHeifer=priceHeifer,
		priceCalf=priceCalf,
		priceBullCalf=priceBullCalf,
		priceCarcassKg=priceCarcassKg,
		oCowY=oCowY,
		oCowWF=oCowWF,
		oOestus=oOestus,
		oDis=oDis,
		.prInvol=NULL,
		.prPosPregTest=NULL,
		.prM=prM,
		.prM2A=NULL,
		.cowAgeParam=list(a=cowAgeParam[1],b=cowAgeParam[2],
		ptr=NULL
	)
	)
})



#########################################################################/**
# @RdocMethod prM
#
# @title "Probability of the conditional latent mean"
#
# \description{
#   Probability of the conditional latent mean \eqn{pr(m_{t+10}|m_{t})}.
# }
#
# @synopsis
#
# \arguments{
#  \item{iM}{ Index of the cube in field \code{stateVar$sM} at time t.}
#  \item{dfc}{ Days from calving. }
#  \item{lac}{ Lactation number. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A matrix where the first row contains the index of \eqn{m_{t+1}} and the second the probability.
# }
#
# \note{Method \code{calcPrM} must have been called before!}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "calcPrM"
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("prM", "HMDPDaily", function(this, lac, dfc, iM, ...){
	return(matrix(this$.prM[[lac+1,dfc+1,iM+1]],nrow=2))
	#return(matrix(this$.prM[[lac,dfc,iM+1]],ncol=2,byrow=TRUE))
})



#########################################################################/*
# @RdocMethod .covM
#
# @title "Calculate covariance matrices for (m_[t+1]|m_[t]) distribution"
#
# \description{
#   Calculate covariance matrices B_t Q_t t(B_t) during a lactation.
# }
#
# @synopsis
#
# \arguments{
#  \item{lac}{ Lactation number. }
#  \item{stages}{ Number of time instances to calculate the covariance for. }
#  \item{...}{Not used.}
# }
#
# \value{
#   Return a list with the covariance for each time instance.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".covM", "HMDPDaily", function(this, lac, stages, ...) {
	ssmModel<-this$oCowY$.ssmModels[[lac]]
	modFilt<-dlmFilter(rnorm(stages),ssmModel) # run Kalman filter to get variance matrices
	matM<-vector("list",length(modFilt$f))
	for (i in 1:length(modFilt$f)) { # store variances
		R<-modFilt$U.R[[i]] %*% diag(modFilt$D.R[i,]^2) %*% t(modFilt$U.R[[i]])
		Q<-as.vector(ssmModel$FF %*% R %*% t(ssmModel$FF) + ssmModel$V)
		B<-R %*% t(ssmModel$FF) * (1/Q)
		matM[[i]]<-Q * B %*% t(B)
	}
	return(matM)
}, private=TRUE)




#########################################################################/**
# @RdocMethod calcPrM
#
# @title "Calculate all probabilities of the conditional latent mean"
#
# \description{
#   Calculate all \eqn{pr(m_{t+10}|m_{t})} which afterwards can be accessed using the prM method.
# }
#
# @synopsis
#
# \arguments{
#  \item{precision}{ Probability must be greater than 10^-precision otherwise assumed zero.}
#  \item{...}{Not used.}
# }
#
# \value{
#   NULL
# }
#
# @author
#
# \seealso{
#   @seeclass. Method @seemethod "prM".
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("calcPrM", "HMDPDaily", function(this, precision = 6, ...){
	ptm <- proc.time()
	this$.prM<-array(list(0),c(length(this$oCowY$.ssmModels)+1,this$maxDfc,this$stateVar$sizeM))    # create 3-dim array prM[lac+1,dfc+1,iM+1] of vectors
	for (lac in 1:length(this$oCowY$.ssmModels)) {
		ssmModel<-this$oCowY$.ssmModels[[lac]]
		covar<-this$.covM(lac,this$maxDfc+1)
		for (dfc in 1:(this$maxDfc-1)) {
			sigma<-covar[[dfc+1]]
			for (iM in 1:this$stateVar$sizeM-1) {
				m<-matrix(this$stateVar$sM[[iM+1]]$center,ncol=1) # Note that since sMIdx starts from zero we have that sMIdx = i corresponds to sM[[i+1]]!!
				mu<-as.vector(ssmModel$GG %*% m)
				pr<-NULL
				for (k in 1:this$stateVar$sizeM-1) {   # for each m_t+1 idx
					lower<-this$stateVar$sM[[k+1]]$cube[1,]
					upper<-this$stateVar$sM[[k+1]]$cube[2,]
					p <- pmvnorm(lower, upper, mu, sigma=sigma)
					if (p>10^-precision) {
						pr<-c(pr,k,p)
					}
				}
				this$.prM[[lac+1,dfc+1,iM+1]]<-pr
			}
		}
	}
	cat("All pr(m_{t}|m_{t-1}) calculated which can be accessed using the prM method.\n")
	print(proc.time() - ptm)
	invisible(NULL)
})


#########################################################################/**
# @RdocMethod csvPrM
#
# @title "Create a csv file with pr(m_[t+1]|m_[t]) "
#
# \description{
#   One row for each lac x dfc x iM and afterwards altering columns containing iM probM idx probM ...
#   Can be read again using \code{read.csv(fileName,header=FALSE)}.
# }
#
# @synopsis
#
# \arguments{
#  \item{fileName}{ File name of the csv. }
#  \item{...}{Passed on to \code{calcPrM} if not already calculated.}
# }
#
# \value{
#   NULL
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("csvPrM", "HMDPDaily", function(this, fileName="prM.csv", ...) {
	if (is.null(this$.prM)) this$calcPrM(...)
	sizes<-dim(this$.prM)
	lth<-0
	for (lac in 1:sizes[1]-1) { # find number of cols
		for (dfc in 1:sizes[2]-1) {
			for (iM in 1:sizes[3]-1) {
				lth <- max(lth,length(this$.prM[[lac+1,dfc+1,iM+1]]))
			}
		}
	}
	con<-file(fileName, "w")
	cat(paste("col",1:(lth+3+1), collapse=",", sep=""),"\n",file=con)
	for (lac in 1:sizes[1]-1) {
		for (dfc in 1:sizes[2]-1) {
			for (iM in 1:sizes[3]-1) {
				cat(lac, ",", dfc, ",", iM, ",", paste(this$.prM[[lac+1,dfc+1,iM+1]],collapse=","), "\n", sep="", file=con)
			}
		}
	}
	close(con)
	cat("Csv file",fileName,"created\n")
	invisible(NULL)
})


#########################################################################/**
# @RdocMethod probInvol
#
# @title "Marginal probability of involuntary culling "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving}
#  \item{lac}{ Lactation number}
#  \item{preg}{ Boolean variable indicating whether pregnancy has been confirmed}
#  \item{dryoff}{ Boolean variable indicating if the animal is dried of on the
#     current day. NB! the function returns a probability of culling if \code{dryoff=TRUE}}
#  \item{matinginterval}{ A (3x1) vector with start and end of mating period and length of drying of period}
#  \item{parmvec}{ A (5 x 1) vector with model parameters for the different stages of the
#    lactation. See Details for further explanation}
#  \item{ageincrease}{ Relative increase in culling rates per lactation.}
#  \item{maxLac}{ Lactation where \code{ageincrease} stops }
#  \item{...}{Not used.}
# }
#
# \details{
# The parmvec consists of 5 five elements. The first 4 elements represents the different
# levels of the intensity/hazard rate, but calculated for a 350 days period, to make
# the specification more comparable to published culling levels.
#
# \tabular{rl}{
#  \tab Probability of involuntary culling within 350 days based on the start hazard rate (at calving).\cr
#  \tab Probability of involuntary culling within 350 days based on the Hazard rate at start of mating.\cr
#  \tab Probability of involuntary culling within 350 days based on the Hazard rate when pregnancy is confirmed.\cr
#  \tab Probability of involuntary culling within 350 days based on the Hazard rate after end of mating but without pregnancy.\cr
#  \tab Probability of involuntary culling after drying-off before next calving.\cr
# }
# From day 1 to matingstart (\code{matinginterval[1]}) the hazard is interpolated between
# $h_1$ and $h_2$
# }
#
# \value{
#   The function returns the probability of involuntary culling in next timestep. The length of the timestep is 1 except if
#   dryoff is TRUE. In this case the length of the timestep  is \code{matingperiod[3]}
# }
#
# \author{Erik Jørgensen.}
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("probInvol", "HMDPDaily", function(this, dfc, lac, preg=FALSE,
	dryoff=FALSE, matinginterval=c(35,250,49),
	parmvec=c(0.3,0.1,0.05,0.2,0.1), ageincrease=0.1, maxLac=5, ...)
{
   if (lac<maxLac) {parity=lac} else {parity=this$maxLac}
   kparity<- (1+ageincrease)^(parity-1)
   h<- (-log(1-parmvec[1:4])/350)*kparity
   pcull<- plogis(qlogis(parmvec[5])+log(kparity))
   #tid<-seq(0:45)
   #hx <- h0+ (h1-h0)/45 * dfc
   hx <- c(exp( log(h[1]) + (log(h[2])-log(h[1]))/matinginterval[1] * dfc), h[2], h[4])
   intervals=c(0,matinginterval[1:2],2000)
   # first ignorere preg and dryoff
   result <- 1-exp(-hx[as.numeric(cut(dfc,intervals))])
   # take dry off and preg into account
   if (dryoff) {result <- pcull+(1-pcull)*(1-exp(-hx[3]*matinginterval[3])) } else {
	  if (preg) {result<- 1-exp(-hx[3])}
   }
   result
})


#########################################################################/*
# @RdocMethod .mapMToM1
#
# @title "Find the interval m belongs to "
#
# \description{
#   Find the interval in field \code{stateVar$sA} \code{m} belongs to.
# }
#
# @synopsis
#
# \arguments{
#  \item{m}{ 2-dim vector containing the conditional mean of (A,X). }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#*/#########################################################################
setMethodS3(".mapMToM1", "HMDPDaily", function(this, m, ...) {
	for (i in 1:nrow(this$stateVar$sA)) {
		if (m[1]>=this$stateVar$sA[i,"min"] & m[1]<=this$stateVar$sA[i,"max"]) {
			return(this$stateVar$sA[i,"idxA"])
		}
	}
	return(NULL)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .mapM1ToM
#
# @title "Given mean of A find the corresponding cube "
#
# \description{
#   Given mean of A find the corresponding cube in field \code{this$stateVar$sM} with A value closest to the mean and X value zero.
# }
#
# @synopsis
#
# \arguments{
#  \item{a}{ The conditional mean of A. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the cube.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".mapM1ToM", "HMDPDaily", function(this, a, ...) { # a: mean A value
	for (i in 1:length(this$stateVar$sM)) {
		if (this$stateVar$sM[[i]]$cube[1,1]<=a & a<=this$stateVar$sM[[i]]$cube[2,1] & this$stateVar$sM[[i]]$cube[1,2]<=0 & 0<=this$stateVar$sM[[i]]$cube[2,2]) {
			return(this$stateVar$sM[[i]]$idxM)
		}
	}
	return(NULL)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .dfc2Week
#
# @title "Return week number given that we are at day dfc. "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".dfc2Week", "HMDPDaily", function(this, dfc, ...) {
	return(floor((dfc-1)/7+1))
}, private=TRUE)


#########################################################################/*
# @RdocMethod .week2Dfc
#
# @title "Return last day in week "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{w}{ Week in lactation.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".week2Dfc", "HMDPDaily", function(this, w, ...) {
	return(w*7)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .dfc2DryWeek
#
# @title "Return week number were dry given a positive preg test at day dfc "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".dfc2DryWeek", "HMDPDaily", function(this, dfc, ...) {
	return(this$.dfc2Week(dfc-this$oOestus$pregTestLth+this$oOestus$gestLth-this$oOestus$dryPeriodLth))
}, private=TRUE)


#########################################################################/*
# @RdocMethod .dfc2DryWeekIdx
#
# @title "Return the dry week index given a positive preg test at day dfc "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dfc}{ Days from calving. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".dfc2DryWeekIdx", "HMDPDaily", function(this, dfc, ...) {
	return(this$.dfc2DryWeek(dfc)-this$stateVar$sDry[1,2]);
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getMaxDryWeekIdx
#
# @title "Returns the index of the longest possible calving interval to be defined at this stage "
#
# \description{
#   Returns the index of the longest possible calving interval to be defined at this stage. Assume that the preg test is taken prior such that the observation at day insemStart+pregTestLgd could be that the cow is pregnant!
# }
#
# @synopsis
#
# \arguments{
#  \item{n}{ Stage number at level 2.}
#  \item{...}{Not used.}
# }
#
# \value{
#  The index of the longest possible calving interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getMaxDryWeekIdx", "HMDPDaily", function(this, n, ...) {
	if (n < this$oOestus$insemStart+this$oOestus$pregTestLth) # status still unknown. Return last idx
		return(this$stateVar$sizeDry-1);
	if (n < this$oOestus$insemFinish+this$oOestus$pregTestLth)
		return(this$.dfc2DryWeek(n)-this$.dfc2DryWeek(this$oOestus$insemStart+this$oOestus$pregTestLth))
	return(this$.dfc2Week(this$oOestus$insemFinish+this$oOestus$gestLth-this$oOestus$dryPeriodLth)-this$.dfc2Week(this$oOestus$insemStart+this$oOestus$gestLth-this$oOestus$dryPeriodLth)) # last possible week idx
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getMinDryWeekIdx
#
# @title "Returns the index of the shortest possible week where dry at this stage. Do not consider the unknown state. "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{n}{ Stage number at level 2.}
#  \item{...}{Not used.}
# }
#
# \value{
#  The index of the shortest possible week where dry.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getMinDryWeekIdx", "HMDPDaily", function(this, n, ...) {
	if (n < this$oOestus$insemStart+this$oOestus$pregTestLth) # status still unknown. Return last idx
		return(this$stateVar$sizeDry-1);
	if (n > this$oOestus$insemStart + this$oOestus$gestLth - this$oOestus$dryPeriodLth)
		return(this$.dfc2Week(n)-this$.dfc2Week(this$oOestus$insemStart+this$oOestus$gestLth-this$oOestus$dryPeriodLth)) # idx of current week
	return(0);
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getLev1States
#
# @title "Returns the number of states at a specified stage of the process at level 1 "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{n}{ Stage number at level 1.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getLev1States", "HMDPDaily", function(this, n, ...) {
	if (n==0) return(1);
	if (n==this$maxLac+1) return(2);
	if (n>this$maxLac+1) return(0);
	return(this$stateVar$sizeA+1)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getLev2States
#
# @title "Returns the number of states at a specified stage of a process at level 2 "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{n}{ Stage number at level 2.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getLev2States", "HMDPDaily", function(this, n, ...) {
	#cat(".getLev2States\n")
	if (n==0) return(1);
	if (n>this$maxDfc) return(0);
	if (this$.getMaxDryWeekIdx(n)==this$stateVar$sizeDry-1) { # pregnancy unknown
		return(this$stateVar$sizeM + 1);
	}
	#cat("end .getLev2States\n")
	return((this$.getMaxDryWeekIdx(n) - this$.getMinDryWeekIdx(n) + 2)*this$stateVar$sizeM + 1);
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getDryWeekIdx
#
# @title "Returns the index of the dry week "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#    \item{n}{ Stage index at level 2.}
#    \item{i}{ State index at level 2.}
#    \item{...}{Not used.}
# }
#
# \value{
#   The index of the interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getDryWeekIdx", "HMDPDaily", function(this, n, i, ...) {
	#cat(".getDryWeekIdx: n",n,"i",i,"\n")
	if (this$.getLev2States(n)<=i+1) return(NA)  # replaced state
	if (i >= this$.getLev2States(n)-1-this$stateVar$sizeM)
		return(this$stateVar$sizeDry-1);
	minI = this$.getMinDryWeekIdx(n);
	return(floor((i / this$stateVar$sizeM) + minI))
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getDryWeek
#
# @title "Return dry week number "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{iDry}{ Index of state variable \code{this$stateVar$sDry}.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
# Returns the dry week corresponding to iDry.
setMethodS3(".getDryWeek", "HMDPDaily", function(this, iDry, ...) {
	return(this$stateVar$sDry[iDry+1,2])
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getMIdx
#
# @title "Returns the index of the latent milk yield corresponding to a specified state at a given stage "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{i}{ State index at level 2.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getMIdx", "HMDPDaily", function(this, i, ...) {
	return(i %% this$stateVar$sizeM)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getIdxM
#
# @title "Find the cube containing m. Values on borders are handled random "
#
# \description{
#   Find the cube containing \code{m}. Values on borders are handled random.
# }
#
# @synopsis
#
# \arguments{
#  \item{m}{ 2-dim vector containing the conditional mean of (A,X). }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the cube.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
# return cube which fit best
setMethodS3(".getIdxM", "HMDPDaily", function(this, m, ...) { # m: mean of (A,X) matrix/df (1 x 2)
	for (j in 1:length(this$stateVar$sM)) {
		if (this$stateVar$sM[[j]]$cube[1,1]<=m[1] & m[1]<=this$stateVar$sM[[j]]$cube[2,1] & this$stateVar$sM[[j]]$cube[1,2]<=m[2] & m[2]<=this$stateVar$sM[[j]]$cube[2,2]) {
			return(this$stateVar$sM[[j]]$idxM)
		}
	}
	return(NULL)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getIdxA
#
# @title "Find the interval containing a. Values on borders are handled random  "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{a}{ The conditional mean of A. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getIdxA", "HMDPDaily", function(this, a, ...) { # a: mean of A
	for (j in 1:nrow(this$stateVar$sA)) {
		if (this$stateVar$sA[j,2]<=a & a<=this$stateVar$sA[j,3]) {
			return(this$stateVar$sA[j,4])
		}
	}
	return(NULL)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getDOC
#
# @title "Find the day of conception  "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{iDry}{ The index of the dry week. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getDOC", "HMDPDaily", function(this, iDry, ...) {
	if (iDry==this$stateVar$sizeDry-1) return(this$maxDfc + 1000)
	return(this$.week2Dfc(this$stateVar$sDry[iDry+1,2])+this$oOestus$dryPeriodLth-this$oOestus$gestLth)
}, private=TRUE)


#########################################################################/*
# @RdocMethod .getStateIdx
#
# @title "Get the index of a state in a stage at level 2  "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{iDry}{ The index of the dry week. }
#  \item{...}{Not used.}
# }
#
# \value{
#   The index of the interval.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".getStateIdx", "HMDPDaily", function(this, d2, iM, iDry, ...) {
	if (iDry==this$stateVar$sizeDry-1) # preg unknown
		return ((this$.getLev2States(d2)-1-this$stateVar$sizeM) + iM)
	minI = this$.getMinDryWeekIdx(d2);
	return(this$stateVar$sizeM * (iDry-minI) + iM)
}, private=TRUE)


#########################################################################/**
# @RdocMethod calcPrM2A
#
# @title "Calculate the transition probabilities from level 2 to level 1  "
#
# \description{
#   Create an array containing all possible values of probA. The array can be accessed using method \code{prM2A}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   NULL (invisible).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("calcPrM2A", "HMDPDaily", function(this, ...) {
	lM<-this$stateVar$sizeM
	lA<-this$stateVar$sizeA
	this$.prM2A<-array(list(),lM)
	eps<-0.0001
	for (lac in 1:1) {
		for (iM in 1:lM-1) {    # for each m_t idx
			minA<-this$stateVar$sM[[iM+1]]$cubeB[1,1]
			maxA<-maxAcube<-this$stateVar$sM[[iM+1]]$cubeB[2,1]
			lgd<-maxA-minA
			idxAmin<-this$.getIdxA(minA+eps)    # first interval
			idxAmax<-this$.getIdxA(maxA-eps)    # last interval
			#cat("iM",iM,"minA",minA,"maxA",maxA,"\n")
			tmp<-NULL
			for (i in idxAmin:idxAmax) {
				tmp<-c(tmp,i)
				ifelse(maxAcube<this$stateVar$sA[i+1,'max'],maxA<-maxAcube, maxA<-this$stateVar$sA[i+1,'max'])
				tmp<-c(tmp,(maxA-minA)/lgd)
				#cat("  idxA",i,"pr",(maxA-minA)/lgd,"\n")
				minA<-maxA
			}
			names(tmp)<-NULL
			this$.prM2A[[iM+1]]<-tmp
		}
	}
	invisible()
})


#########################################################################/**
# @RdocMethod prM2A
#
# @title "Transition probabilities from level 2 to level 1"
#
# \description{
#   Probability for transition to a state at level one \eqn{pr(s|m_{t})}.
# }
#
# @synopsis
#
# \arguments{
#  \item{iM}{ Index of the cube in field \code{stateVar$sM} at time t.}
#  \item{...}{Not used.}
# }
#
# \value{
#   A matrix where the first row contains the index of \eqn{s} and the second the probability.
# }
#
# \note{Method \code{calcPrM2A} must have been called before!}
#
# @author
#
# \seealso{
#   @seeclass
#   @seemethod "calcPrM2A"
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("prM2A", "HMDPDaily", function(this, iM, ...) {
	return(matrix(this$.prM2A[[iM+1]],nrow=2))
	#return(this$.prM2A[[iM+1]])
})


#########################################################################/**
# @RdocMethod csvPrM2A
#
# @title "Create a csv file with the transition probabilities from level 2 to level 1 "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{fileName}{ File name of the csv. }
#  \item{...}{Not used.}
# }
#
# \value{
#   NULL
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("csvPrM2A", "HMDPDaily", function(this, fileName="prM2A.csv", ...) {
	if (is.null(this$.prM2A)) this$calcPrM2A()
	sizes<-dim(this$.prM2A)
	con<-file(fileName, "w")
	for (iM in 1:sizes[1]-1) {
		cat(iM, ",", paste(this$.prM2A[[iM+1]],collapse=","), "\n", sep="", file=con)
	}
	close(con)
	cat("Csv file",fileName,"created\n")
	invisible(NULL)
})


#########################################################################/**
# @RdocMethod cowAge
#
# @title "Return the cow age in days at the start of a given lactation"
#
# \description{
#   @get "title". An integer returned.
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
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("cowAge", "HMDPDaily", function(this, lac, ...) {
	return(floor(this$.cowAgeParam$a + this$.cowAgeParam$b*lac))
})


#########################################################################/**
# @RdocMethod rewardKeep
#
# @title "Reward of keeping the cow (action keep at level 2)"
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
#  \item{iDry}{ Index of dry week. }
#  \item{iM}{ Index of m. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("rewardKeep", "HMDPDaily", function(this, lac, dfc, iDry, iM, ...) {
	#cat("rewardKeep: lac",lac,"dfc",dfc,"iDry",iDry,"iM",iM,"\n")
	m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
	reward = this$rewardMilk(lac, dfc, iM);
	#print(reward)
	milkECM = this$oCowY$dailyYieldECM(m, lac, dfc)
	#print(milkECM)
	con <- this$.getDOC(iDry)
	feedCost = this$oCowWF$energy(this$cowAge(lac), dfc, con, milkECM, this$oOestus$gestLth)[6] * this$priceSFU
	#print(feedCost)
	ic <- this$probInvol(dfc,lac,preg= iDry!=this$stateVar$sizeDry-1)
	#print(ic)
	repRew<-this$rewardReplace(lac,dfc,iDry)
	#print(repRew)
	#print((1-ic)*(reward - feedCost) + ic*repRew)
	#cat("end rewardKeep\n")
	return((1-ic)*(reward - feedCost) + ic*repRew)
})


#########################################################################/**
# @RdocMethod rewardReplace
#
# @title "Reward of replacing the cow (action replace at level 2)"
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
#  \item{iDry}{ Index of dry week. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("rewardReplace", "HMDPDaily", function(this, lac, dfc, iDry, ...) {
	con <- this$.getDOC(iDry)
	return(this$rewardCarcass(lac, dfc, con))
})


#########################################################################/**
# @RdocMethod rewardCarcass
#
# @title "Reward of carcass"
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
#  \item{con}{ Day of conception (use a con larger than dfc if not pregnant). }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("rewardCarcass", "HMDPDaily", function(this, lac, dfc, con, ...) {
	age<-this$cowAge(lac)
	return(this$oCowWF$bodyW(age, dfc, con)*0.5*this$priceCarcassKg)
})


#########################################################################/**
# @RdocMethod rewardMilk
#
# @title "Reward of milk"
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
#  \item{iM}{ Index of m. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("rewardMilk", "HMDPDaily", function(this, lac, dfc, iM, ...) {
	m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
	yield = this$oCowY$dailyYieldECM(m, lac, dfc)
	return(this$priceECM*yield)
})


#########################################################################/**
# @RdocMethod rewardDry
#
# @title "Reward of drying the cow"
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
#  \item{iDry}{ Index of dry week. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("rewardDry", "HMDPDaily", function(this, lac, dfc, iDry, ...) {
	prDryIC = this$probInvol(400,lac=lac,preg=TRUE,dryoff=TRUE)
	rewRep = prDryIC * this$rewardReplace(lac,dfc+floor(this$oOestus$dryPeriodLth/2),iDry)
	rewKeep = (1-prDryIC) * (this$oOestus$probCalf*this$priceCalf+
		 this$oOestus$probBullCalf*this$priceBullCalf - this$priceSFU*this$oCowWF$.sfuDry*this$oOestus$dryPeriodLth)
	return (rewRep+rewKeep)
})


#########################################################################/**
# @RdocMethod transPrKeep
#
# @title "Transition probabilities (action keep at level 2)"
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
#  \item{iM}{ Index of m. }
#  \item{iDry}{ Index of dry week. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A vector with values (scp1, idx1, pr1, ...).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("transPrKeep", "HMDPDaily", function(this, lac, dfc, iM, iDry, ...) {
	pr <- rbind(1,this$prM(lac,dfc,iM))
	#cat("prM:\n"); print(pr)
	ic <- this$probInvol(dfc,lac,preg= iDry!=this$stateVar$sizeDry-1)
	#cat("prIc:\n"); print(ic)
	if (iDry==this$stateVar$sizeDry-1 &
		dfc<this$oOestus$insemFinish+this$oOestus$pregTestLth &
		dfc>=this$oOestus$insemStart+this$oOestus$pregTestLth-1) {  # i.e. pregnancy unknown and possible to run preg test next day
		pPTest = this$oOestus$probPregTest(lac, dfc+1)
		dryWeekIdx = this$.dfc2DryWeekIdx(dfc+1);     # idx if positive preg test
		prPreg<-pr
		pr[2,]<-this$.getStateIdx(dfc+1, pr[2,], iDry)    # not preg idx
		pr[3,]<-(1-ic)*(1-pPTest)*pr[3,]            # not preg pr
		prPreg[2,]<-this$.getStateIdx(dfc+1, prPreg[2,], dryWeekIdx) # preg idx
		prPreg[3,]<-(1-ic)*pPTest*prPreg[3,]                   # preg pr
		pr<-cbind(pr,prPreg)
	}
	else {  # cow pregnant or not possible to take preg test
		pr[2,]<-this$.getStateIdx(dfc+1, pr[2,], iDry)
		pr[3,]<-(1-ic)*pr[3,]
	}
	pr<-cbind(pr,c(1,this$.getLev2States(dfc+1)-1,ic))
	#cat("sum:",sum(pr[3,])," ")
	return(as.numeric(pr))
})


#########################################################################/**
# @RdocMethod transPrDry
#
# @title "Transition probabilities (action dry at level 2)"
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
#  \item{iM}{ Index of m. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("transPrDry", "HMDPDaily", function(this, lac, iM, ...) {
	prDryIC = this$probInvol(400,lac=lac,preg=TRUE,dryoff=TRUE)
	if (lac>=this$maxLac) return (c(0,0,1))   # 'Replace It' state
	pr<-rbind(0,this$prM2A(iM))
	pr[3,]<-pr[3,]*(1-prDryIC)
	pr<-cbind(pr,c(0,this$stateVar$sizeA,prDryIC))
	return(as.numeric(pr))
})


#########################################################################/**
# @RdocMethod genBinaryR
#
# @title "Generate the binary files describing the HMDP model using plain R"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{saveCsv}{ Run method @seemethod{csvPrM}.  }
#  \item{...}{Arguments passed to \code{calcPrM} and \code{csvPrM}.}
# }
#
# \details{If the transition probabilities prM have not been calculated they are calculated and saved to a csv file. }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("genBinaryR", "HMDPDaily", function(this, prefix="", saveCsv=FALSE, ...) {
	idS<-matrix(NA,nrow=this$maxLac,ncol=this$stateVar$sizeM+1)      # idS[lac,j] = sId of state j-1 at dfc = 1 in lac (second level)

	icAction<-function(d1,d2) {
		#cat("icAction:\n")
		size1Next<-this$.getLev1States(d1+1)
		pr<-c(0,size1Next-1,1)
		w$action(label="Dummy",weights=c(0,0,0),prob=pr)
		w$endAction()
		#cat("end icAction\n")
	}

	replaceAction<-function(d1,d2,iDry) {
		#cat("replaceAction:\n")
		size1Next<-this$.getLev1States(d1+1)
		pr<-c(0,size1Next-1,1)
		w$action(label="Replace",weights=c(0,this$rewardReplace(lac=d1,dfc=d2,iDry=iDry),0),prob=pr)
		w$endAction()
		#cat("end replaceAction\n")
	}

	dryAction<-function(d1,d2,iM,iDry) {
		#cat("dryAction:\n")
		w$action(label="Dry",weights=c(this$oOestus$dryPeriodLth,this$rewardDry(lac=d1,dfc=d2,iDry),0),prob=this$transPrDry(lac=d1,iM=iM))
		w$endAction()
		#cat("end dryAction:\n")
	}

	keepAction<-function(d1, d2, iM, iDry) {
		#cat("keepAction: d1",d1,"d2",d2,"iM",iM,"iDry",iDry,"\n")
		m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
		milkECM = this$oCowY$dailyYieldECM(m, lac=d1, dfc=d2)
		w$action(label="Keep",weights=c(1,this$rewardKeep(lac=d1,dfc=d2,iM=iM,iDry=iDry),milkECM),prob=this$transPrKeep(lac=d1,dfc=d2,iM=iM,iDry=iDry))
		w$endAction()
		#cat("end keepAction:\n")
	}

	createLev2DummyStage<-function() {
		#cat("createLev2DummyStage:\n")
		cat(" dfc = 0 ",sep="")
		w$stage()
			w$state(label="Dummy")
				rew<-c(0,0,0)
				prob<-c(1,this$.mapM1ToM(this$stateVar$sA[1,1]),1)
				w$action(label="Dummy",weights=rew,prob)
				w$endAction()
			w$endState()
		w$endStage()
		#cat("\n")
	}

	createLev2Stage<-function(d1,d2) {
		#cat("createLev2Stage:\n")
		size1Next<-this$.getLev1States(d1+1)
		size<-this$.getLev2States(d2)
		sizeNext<-this$.getLev2States(d2+1)
		cat(" dfc = ",d2," ",sep="")
		w$stage()
			for (s2 in 0:(size-1)) {
				iM<-this$.getMIdx(s2)
				iDry<-this$.getDryWeekIdx(d2, s2)
				if (s2==size-1) lbl<-"Replaced due to IC" else lbl<-paste(iM,",",iDry,sep="")
				w$state(label=lbl)
					while (TRUE) {      # dummy loop so can use break
						if (s2==size-1) {
							icAction(d1,d2)    # last state (IC state)
							break;
						}
						if (d2==this$maxDfc) { # if last stage
							if (iDry==this$stateVar$sizeDry-1) {
								replaceAction(d1,d2,iDry)    # if not pregnant
							} else {  # dry or replace
								dryAction(d1,d2,iM,iDry)
								replaceAction(d1,d2,iDry)
							}
							break
						}
						if (iDry==this$stateVar$sizeDry-1) {    # if cow not pregnant
							keepAction(d1, d2, iM, iDry)
							replaceAction(d1,d2,iDry)
							break;
						}
						if (d2==this$stateVar$sDry[iDry+1,2]*7) { # dry the cow
							dryAction(d1,d2,iM,iDry)
							replaceAction(d1,d2,iDry)
							break;
						}
						# else a normal state with two actions
						keepAction(d1, d2, iM, iDry)
						replaceAction(d1,d2,iDry)
						break;
					}
				w$endState()
			}
		w$endStage()
		#cat("end createLev2Stage \n")
	}

	# only called when d2 = 1 so can store id's in idS
	createLev2Stage1<-function(d1) {     # do the same as createLev2Stage except store idS
		#cat("createLev2Stage1:\n")
		d2 = 1
		size1Next<-this$.getLev1States(d1+1)
		size<-this$.getLev2States(d2)
		sizeNext<-this$.getLev2States(d2+1)
		cat(" dfc = ",d2," ",sep="")
		w$stage()
			for (s2 in 0:(size-1)) {
				iM<-this$.getMIdx(s2)
				iDry<-this$.getDryWeekIdx(d2, s2)
				if (s2==size-1) lbl<-"Replaced due to IC" else lbl<-paste(iM,",",iDry,sep="")
				idS[d1,s2+1]<<-w$state(label=lbl) # save state id
					while (TRUE) {      # dummy loop so can use break
						if (s2==size-1) {
							icAction(d1,d2)    # last state (IC state)
							break;
						}
						if (iDry==this$stateVar$sizeDry-1) {    # if cow not pregnant
							keepAction(d1, d2, iM, iDry)
							replaceAction(d1,d2,iDry)
							break;
						}
					}
				w$endState()
			}
		w$endStage()
		#cat("end createLev2Stage1:\n")
	}

	# build the process for level 1
	createLev1Stage<-function(d1) {
		#cat("createLev1Stage:\n")
		cat("Lac ",d1,":\n",sep="")
		size <- this$.getLev1States(d1)
		w$stage()
			for (s1 in 0:(size-1)) {
				if (s1==size-1) {   # replaced state
					w$state(label="Replaced")
						cat("Replaced\n")
						w$action(label="Dummy",weights=c(0,0,0),prob=c(0,0,1))
						w$endAction()
					w$endState()
				} else { # idxA state
					w$state(label=paste("idxA=",s1,sep=""))
						cat(paste("idxA = ",s1,sep=""),"\n")
						if (s1==0) {    # state where define child process
							rew<-c(0,0,0)
							prob<-c(2,0,1)  # go to dummy at level 2
							w$action(label="Dummy",rew,prob)
								w$process()     # level 2
									for (d2 in 0:this$maxDfc) {   # this$maxDfc
										if (d2==0) createLev2DummyStage()
										if (d2==1) createLev2Stage1(d1)
										if (d2>1) createLev2Stage(d1,d2)
									}
								w$endProcess()
							w$endAction()
							cat("\n")
						} else {    # just make a hyperarc to stage 1 in the child process defined above (shared child process)
							rew<-c(0,0,0)
							prob<-c(3,idS[d1,this$.mapM1ToM(this$stateVar$sA[s1+1,1])+1],1)
							w$action(label="Dummy",rew,prob)
							w$endAction()
						}
					w$endState()
				}
			}
		w$endStage()
	}

	createLev1DummyStage<-function(d1) {
		#cat("createLev1DummyStage:\n")
		if (d1==0) cat("Lac 0:\n") else cat("Lac ",this$maxLac+1,":\n",sep="")
		w$stage()
			if (d1==0) {
				w$state(label="Dummy")
					cat("Dummy\n")
					rew<-c(0,-this$priceHeifer,0)
					prob<-c(1,as.integer(this$.getIdxA(0)),1)
					w$action(label="Buy cow",rew,prob)
					w$endAction()
				w$endState()
			} else {
				w$state(label="Replace it")
					cat("Replace it\n")
					rew<-c(0,this$rewardCarcass(this$maxLac+1, 1, this$maxDfc),0)
					prob<-c(0,0,1)
					w$action(label="Replace",rew,prob)
					w$endAction()
				w$endState()
				w$state(label="Replaced")
					cat("Replaced\n")
					rew<-c(0,0,0)
					prob<-c(0,0,1)
					w$action(label="Dummy",rew,prob)
					w$endAction()
				w$endState()
			}
		w$endStage()
	}

	if (is.null(this$.prM)) this$calcPrM(...)
	if (is.null(this$.prM2A)) this$calcPrM2A()
	if (saveCsv) this$csvPrM(...)

	ptm <- proc.time()
	w<-binaryMDPWriter(prefix)
	w$setWeights(c("Time","Reward","Yield"))
	w$process()     # founder level (level 0)
		w$stage()
			w$state(label="Dummy")
				w$action(label="Dummy",weights=c(0,0,0),prob=c(2,0,1))
					w$process()     # level 1
						for (d1 in 0:(this$maxLac+1)) {
							cat("\n")
							if (d1==0 | d1==this$maxLac+1) createLev1DummyStage(d1)
							else createLev1Stage(d1)
						}
					w$endProcess()
				w$endAction()
			w$endState()
		w$endStage()
	w$endProcess()
	w$close()
	cat("\nFile(s) created.\n")
	print(proc.time() - ptm)
})



#########################################################################/**
# @RdocMethod genStates2SQLite
#
# @title "Save the states of the HMDPDaily in a SQLite database table"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dsnName}{Name of the ODBC DSN.}
#  \item{tblPrefix}{Prefix to the table named <tblPrefix>_stateIdx.}
#  \item{fileName}{The csv file where the states are saved first.}
#  \item{loadIntoDB}{If true add the csv file to the DB using SQL statements. This may be slow for large csv's. Here you can create the table manually afterwards (e.g. using .import from SQLite commandline or the program TkSQLite (free)).}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("genStates2SQLite", "HMDPDaily", function(this, dsnName="mdp_models", tblPrefix,
	fileName=paste(tblPrefix,"_stateIdx.csv", sep=""), loadIntoDB=TRUE, ...) {
	nS1<-function(lac) {    # number of states at level one
		if (lac==0) return(0)
		if (lac==this$maxLac+1) return(0:1)
		return(0:this$stateVar$sizeA)
	}

	nS2<-function(dfc) {    # number of states at level two
		if (dfc==0) return(0)
		if (this$.getMaxDryWeekIdx(dfc)==this$stateVar$sizeDry-1) return(this$stateVar$sizeM)
		return((this$.getMaxDryWeekIdx(dfc) - this$.getMinDryWeekIdx(dfc) + 2)*this$stateVar$sizeM)
	}

	state<-function(idx, label=NULL) {
		tmp<-gsub("NA","",paste(idx,collapse=","))
		if (!is.null(label)) tmp<-paste(tmp,',"',label,'"\n',sep='')
		cat(tmp, file=fS) # write line to file
	}

	labelS1<-function(d1,s1) {
		if (d1==0) return("Buy cow")
		if (d1==this$maxLac+1 & s1==0) return("Replaced")
		if (d1==this$maxLac+1 & s1==1) return("Replace it")
		if (s1==this$stateVar$sizeA) return("Replaced")
		return(paste("idxA=",s1,sep=""))
	}

	labelS2<-function(d2,s2,last) {
		if (d2==0) return("Dummy")
		if (last) return("Replaced due to IC")
		iM<-this$.getMIdx(s2)
		iDry<-this$.getDryWeekIdx(d2, s2)
		return(paste(iM,",",iDry,sep=""))
	}

	ptm <- proc.time()
	# calc number of states
	vS1<-rep(NA,this$maxLac+2)  # number of level 1 states (stage i stored in index i+1)
	for (d1 in 0:(this$maxLac+1)) vS1[d1+1]<-length(nS1(d1))
	vS2<-rep(NA,this$maxDfc+1)  # number of level 2 states (stage i stored in index i+1)
	for (d2 in 0:this$maxDfc) vS2[d2+1]<-nS2(d2)+1
	# calc labels
	matL1<-matrix("",nrow=this$maxLac+2,ncol=max(vS1))    # matrix of level one labels (stage d1 is in row i+1 and state s1 is in col j+1)
	for (d1 in 0:(this$maxLac+1)) {
		for (s1 in 0:vS1[d1+1]-1) matL1[d1+1,s1+1]<-labelS1(d1,s1)
	}
	matL2<-matrix("",nrow=this$maxDfc+1,ncol=max(vS2))    # matrix of level one labels (stage d1 is in row i+1 and state s1 is in col j+1)
	for (d2 in 0:this$maxDfc) {
		lastS2<-vS2[d2+1]-1
		for (s2 in 0:lastS2) {
			matL2[d2+1,s2+1]<-labelS2(d2,s2,s2==lastS2)
		}
	}

	.Call("DAIRY_GenStates2Csv", as.integer(vS1), as.integer(vS2), matL1, matL2, fileName)
	cat("\nFile",fileName,"created.\n")
	print(proc.time() - ptm)

#     fS <- file(fileName, "w")
#     cat("d1,s1,a1,d2,s2,label\n", file=fS)
#     idx<-rep(NA,5)  # containing the stage, state or action idx's
#     state(idx,"Dummy")  # founder state
#     cat("Lac: ")
#     for (d1 in 0:(this$maxLac+1)) {
#         cat(d1,"")
#         idx[1]<-d1
#         for (s1 in 0:vS1[d1+1]) {
#             idx[2]<-s1
#             state(idx, label=matL1[d1+1,s1+1])
#             if (s1==0 & d1>0 & d1<this$maxLac+1) { # if subprocess
#                 for (a1 in 0:0) {
#                     idx[3]<-a1
#                     for (d2 in 0:this$maxDfc) {
#                         idx[4]<-d2
#                         lastS2<-vS2[d2+1]
#                         for (s2 in 0:lastS2) {
#                             idx[5]<-s2
#                             state(idx, label=matL2[d2+1,s2+1])
#                         }
#                         idx[5]<-NA
#                     }
#                     idx[4]<-NA
#                 }
#             }
#             idx[3]<-NA
#         }
#         idx[2]<-NA
#     }
#     cat("\n")
#    close(fS)
	if (loadIntoDB) this$.addStateCsvToDB(dsnName, tblPrefix, fileName)
})


#########################################################################/*
# @RdocMethod .addStateCsvToDB
#
# @title "Save the states of the HMDP in a SQLite database table"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dsnName}{Name of the ODBC DSN.}
#  \item{tblPrefix}{Prefix to the table named <tblPrefix>_stateIdx.}
#  \item{fileName}{The csv file where the states are loaded from.}
#  \item{maxRows}{Max number of rows read from the csv at once.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/#########################################################################
setMethodS3(".addStateCsvToDB", "HMDPDaily", function(this, dsnName, tblPrefix,
	fileName, maxRows=500, ...) {
	ptm <- proc.time()
	con <- odbcConnect(dsnName)
	tblName<-paste(tblPrefix,"_stateIdx",sep="")
	sqlDrop(con, tblName, errors=FALSE)
	sql<-paste("CREATE TABLE [", tblName, "] (
		[sId] INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
		[d1] INTEGER,
		[s1] INTEGER,
		[a1] INTEGER,
		[d2] INTEGER,
		[s2] INTEGER,
		[label] VARCHAR(255));", sep="")
	sqlQuery(con,sql)
	#sql<-paste("INSERT INTO", tblName, "VALUES (0,NULL,NULL,NULL,NULL,NULL,'Dummy')")
	#sqlQuery(con,sql)   # add founder state
	skip<-1 # skip col names
	sId<-0
	while (TRUE) {
		dat<-read.table(fileName, sep=",", quote='"', nrows = maxRows, skip = skip,
			col.names=c("d1","s1","a1","d2","s2","label"),
			colClasses=c("integer","integer","integer","integer","integer","character"),
			check.names = FALSE)
		if (nrow(dat)==0) break
		dat$sId<-sId:(nrow(dat)+sId-1)
		sId<-nrow(dat)+sId
		sqlSave(con,dat,tablename=tblName,rownames=FALSE,append=TRUE)
		skip<-skip+nrow(dat)
	}
	#sql<-paste('ALTER TABLE', tblName, 'ADD sId INTEGER AUTOINCREMENT')
	close(con)
	cat("Added states to DB table",tblName,".\n")
	proc.time() - ptm
})


#########################################################################/**
# @RdocMethod genActions2SQLite
#
# @title "Generate the actions to be added to the SQLite DB"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dsnName}{Name of the ODBC DSN.}
#  \item{tblPrefix}{Prefix to the table named <tblPrefix>_stateIdx.}
#  \item{maxRows}{Max number of rows read from the csv at once.}
#  \item{startSId}{The state to start generating actions for (if already generated actions for the states below startSId.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("genActions2SQLite", "HMDPDaily", function(this, dsnName="mdp_models",
	tblPrefix, maxRows=500, startSId=0,...)
{
	action<-function(label=NULL,weights,prob) {     # prob = vector with (scp,idx,pr) ...
		prob<-matrix(prob, ncol=3, byrow=TRUE)
		sql<-paste("INSERT INTO ", tblPrefix, "_actionWeight VALUES (", aId, ",", paste(weights,collapse=","), ")", sep="")
		sqlQuery(con1,sql)
		if (!is.null(label)) {
			sql<-paste("INSERT INTO ", tblPrefix, "_actionIdxLbl VALUES (", aId, ",'", label, "')", sep="")
			sqlQuery(con1,sql)
		}
		for (i in 1:nrow(prob)) {
			sql<-paste("INSERT INTO ", tblPrefix, "_actionIdx VALUES (", aId, ",", sId, ",", paste(prob[i,],collapse=","),")", sep="")
			sqlQuery(con1,sql)
		}
		aId<<-aId+1
	}

	lev0Action<-function() {
		action(weights=c(0,0,0),prob=c(2,0,1),label="Dummy")
	}

	lev1Actions<-function() {     # state = (sId d1,s1,a1,d2,s2,label)
		if (d1==0) {
			rew<-c(0,-this$priceHeifer,0)
			prob<-c(1,as.integer(this$stateVar$sizeA/2),1)
			action(label="Buy cow",rew,prob)
			return()
		}
		if (d1==this$maxLac+1) {
			if (s1==0) {    # replaced state
				action(label="Dummy",weights=c(0,0,0),prob=c(0,0,1))
			}
			if (s1==1) {
				rew<-c(0,this$rewardCarcass(this$maxLac+1, 1, this$maxDfc),0)
				prob<-c(0,0,1)
				action(label="Replace",rew,prob)
			}
			return()
		}
		if (label=="Replaced") {
			action(label="Dummy",weights=c(0,0,0),prob=c(0,0,1))
			return()
		}
		# else must be a idxA state
		if (s1==0) {    # state where define child process
			rew<-c(0,0,0)
			prob<-c(2,0,1)  # go to dummy at level 2
			action(label="Dummy",rew,prob)
		} else {    # just make a hyperarc to stage 1 in the child process (shared child process)
			iM<-this$.mapM1ToM(this$stateVar$sA[s1+1,1])
			id<-sStart[sStart$d1==d1 & sStart$s2==iM,'sId']
			action(label="Dummy",weights=c(0,0,0),prob=c(3,id,1))
		}
		return()
	}

	lev2Actions<-function() {
		if (d2==0) {  # dummy stage
			rew<-c(0,0,0)
			prob<-c(1,this$.mapM1ToM(this$stateVar$sA[1,1]),1)
			action(label="Dummy",weights=rew,prob)
			return()
		}
		if (label=="Replaced due to IC") {
			icAction()
			return()
		}
		iM<-as.numeric(unlist(strsplit(label,",")))
		iDry<-iM[2]
		iM<-iM[1]
#        size1Next<-this$.getLev1States(d1+1)
#        size<-this$.getLev2States(d2)
#        sizeNext<-this$.getLev2States(d2+1)
		while (TRUE) {      # dummy loop so can use break
			if (d2==this$maxDfc) { # if last stage
				if (iDry==this$stateVar$sizeDry-1) {
					replaceAction(iDry)    # if not pregnant
				} else {  # dry or replace
					dryAction(iM,iDry)
					replaceAction(iDry)
				}
				break;
			}
			if (iDry==this$stateVar$sizeDry-1) {    # if cow not pregnant
				keepAction(iM, iDry)
				replaceAction(iDry)
				break;
			}
			if (d2==this$stateVar$sDry[iDry+1,2]*7) { # dry the cow
				dryAction(iM,iDry)
				replaceAction(iDry)
				break;
			}
			# else a normal state with two actions
			keepAction(iM, iDry)
			replaceAction(iDry)
			break;
		}
		return()
	}

	icAction<-function() {
		size1Next<-this$.getLev1States(d1+1)
		ifelse(d1==this$maxLac, pr<-c(0,0,1), pr<-c(0,size1Next-1,1))
		action(label="Dummy",weights=c(0,0,0),prob=pr)
	}

	replaceAction<-function(iDry) {
		size1Next<-this$.getLev1States(d1+1)
		ifelse(d1==this$maxLac, pr<-c(0,0,1), pr<-c(0,size1Next-1,1))
		action(label="Replace",weights=c(0,this$rewardReplace(lac=d1,dfc=d2,iDry=iDry),0),prob=pr)
	}

	dryAction<-function(iM,iDry) {
		action(label="Dry",weights=c(this$oOestus$dryPeriodLth,this$rewardDry(lac=d1,dfc=d2,iDry),0),prob=this$transPrDry(lac=d1,iM=iM))
	}

	keepAction<-function(iM, iDry) {
		m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
		milkECM = this$oCowY$dailyYieldECM(m, lac=d1, dfc=d2)
		action(label="Keep",weights=c(1,this$rewardKeep(lac=d1,dfc=d2,iM=iM,iDry=iDry),milkECM),prob=this$transPrKeep(lac=d1,dfc=d2,iM=iM,iDry=iDry))
	}

	ptm <- proc.time()

	if (is.null(this$.prM)) {
		this$calcPrM()
		this$csvPrM(...)
	}
	if (is.null(this$.prM2A)) this$calcPrM2A()

	con <- odbcConnect(dsnName)
	con1 <- odbcConnect(dsnName)    # connection used in action function
	# global var
	aId<-0
	sId<-0
	d1<-s1<-d2<-s2<-label<-NULL
	# first find all states at stage one in level 2 processes (such that can link to shared child)
	sql<-paste("SELECT sId,d1,s2 FROM ", tblPrefix, "_stateIdx WHERE d2=1",sep="")
	sStart<-sqlQuery(con,sql)   # states in the start of the process
	# get states
	sql<-paste("SELECT * FROM ", tblPrefix, "_stateIdx", sep="")
	if (startSId>0) sql<-paste(sql,"WHERE sId >=",startSId)
	if (startSId==0) {  # if gen actions for all states
		this$.createActionTables(dsnName,tblPrefix)
		states<-sqlQuery(con,sql,max=1,stringsAsFactors=FALSE)     # get dummy state in founder
		lev0Action()    # add level 0 action
		states<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	} else {
		aId<-sqlQuery(con,paste("SELECT MAX(aId) FROM ",tblPrefix, "_actionIdx",sep=""),stringsAsFactors=FALSE)+1
		aId<-aId[1,1] # to numeric
		states<-sqlQuery(con,sql,max=maxRows,stringsAsFactors=FALSE)
	}

	cat("Added actions for states: ")
	while(!is.numeric(states)) {
		if (nrow(states)==0) break
		for (r in 1:nrow(states)) { # state = (sId d1,s1,a1,d2,s2,label) = (sId lac,s1,a1,dfc,s2,label)
			sId<-states[r,1]; d1<-states[r,2]
			s1<-states[r,3]; d2<-states[r,5]
			label<-states[r,7]
			if (is.na(d2)) {  # if level 1 state
				lev1Actions()
			} else {  # level 2
				lev2Actions()
			}
		}
		cat(nrow(states),"")
		states<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	}
	cat("\nActions added.\n")
	proc.time() - ptm
	close(con)
	close(con1)
})



#########################################################################/*
# @RdocMethod .createActionTables
#
# @title "Create the action tables necessary to method genActions2SQLite"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dsnName}{Name of the ODBC DSN.}
#  \item{tblPrefix}{Prefix to the table named <tblPrefix>_stateIdx.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/#########################################################################
setMethodS3(".createActionTables", "HMDPDaily", function(this, dsnName="mdp_models",
	tblPrefix, ...) {
	con <- odbcConnect(dsnName)
	sql<-paste("CREATE TABLE [", tblPrefix, "_actionIdx] (
	  [aId] INTEGER,
	  [sId] INTEGER,
	  [scp] INTEGER,
	  [idx] INTEGER,
	  [pr] DOUBLE);", sep="")
	sqlDrop(con, paste(tblPrefix, "_actionIdx",sep=""), errors=FALSE)
	sqlQuery(con,sql)

	sql<-paste("CREATE TABLE [", tblPrefix, "_actionIdxLbl] (
	  [aId] INTEGER,
	  [label] TEXT);", sep="")
	sqlDrop(con, paste(tblPrefix, "_actionIdxLbl",sep=""), errors=FALSE)
	sqlQuery(con,sql)

	sql<-paste("CREATE TABLE [", tblPrefix, "_actionWeight] (
	  [aId] INTEGER NOT NULL,
	  [w0] DOUBLE,
	  [w1] DOUBLE,
	  [w2] DOUBLE);", sep="")
	sqlDrop(con, paste(tblPrefix, "_actionWeight",sep=""), errors=FALSE)
	sqlQuery(con,sql)
	close(con)
})


#########################################################################/**
# @RdocMethod convertSQLite2Binary
#
# @title "Convert the HMDP stored in SQLite tables to binary files"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dsnName}{Name of the ODBC DSN.}
#  \item{tblPrefixStates}{Prefix to the table named <tblPrefix>_stateIdx.}
#  \item{tblPrefixActions}{Prefix to the tables for the actions.}
#  \item{maxRows}{Max number of rows read from the csv at once.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("convertSQLite2Binary", "HMDPDaily", function(this, dsnName="mdp_models",
	tblPrefixActions, tblPrefixStates=tblPrefixActions, binNamePrefix=tblPrefixActions,
	maxRows=500, ...)
{
	con <- odbcConnect(dsnName)
	# create/open binary files
	binNames<-c("stateIdx.bin","stateIdxLbl.bin","actionIdx.bin","actionIdxLbl.bin",
		"actionWeight.bin","actionWeightLbl.bin","transProb.bin")
	binNames<-paste(binNamePrefix,"_",binNames,sep="")
	fS <- file(binNames[1], "wb")
	fSLbl <- file(binNames[2], "wb")
	fA <- file(binNames[3], "wb")
	fALbl <- file(binNames[4], "wb")
	fACost <- file(binNames[5], "wb")
	fACostLbl <- file(binNames[6], "wb")
	fTransP <- file(binNames[7], "wb")
	# create fS
	sql<-paste("SELECT 0 AS d0,0 AS s0,0 AS a0,d1,s1,a1,d2,s2,-1 AS sep FROM ",tblPrefixStates,"_stateIdx", sep="")
	dat<-sqlQuery(con,sql,max=maxRows,stringsAsFactors=FALSE)
	while(!is.numeric(dat)) {
		if (nrow(dat)==0) break
		dat<-as.integer(t(dat))
		writeBin(dat[!is.na(dat)],fS)
		dat<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	}
	# create fSLbl
	sql<-paste("SELECT sId,label FROM ",tblPrefixStates,"_stateIdx WHERE label is not NULL", sep="")
	dat<-sqlQuery(con,sql,max=maxRows,stringsAsFactors=FALSE)
	while(!is.numeric(dat)) {
		if (nrow(dat)==0) break
		dat<-as.character(t(dat))
		writeBin(dat,fSLbl)
		dat<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	}
	# create fA
	sql<-paste("SELECT MAX(sId) FROM ",tblPrefixActions,"_actionIdx", sep="")
	sIdMax<-sqlQuery(con,sql)[1,1]
	sId<-0
	sql<-paste("SELECT sId,scp,idx FROM ",tblPrefixActions,"_actionIdx WHERE sId BETWEEN ",sId," AND ",sId+100, sep="")
	dat<-sqlQuery(con,sql,stringsAsFactors=FALSE)
	while(TRUE) {
		vec<-NULL
		daply(dat, .(sId), .fun = function(df) {    # return the number of transitions for a given sId
			tmp<-df[1,1]    # sId
			df$sId<-NULL
			tmp<-c(tmp,as.integer(t(df)),-1)
			vec<<-c(vec,tmp)
			nrow(df)
		})
		writeBin(as.integer(vec),fA)
		sId<-sId+101
		if (sId>sIdMax) break;
		sql<-paste("SELECT sId,scp,idx FROM ",tblPrefixActions,"_actionIdx WHERE sId BETWEEN ",sId," AND ",sId+100, sep="")
		dat<-sqlQuery(con,sql,stringsAsFactors=FALSE)
	}
	# create fALbl
	sql<-paste("SELECT aId,label FROM ",tblPrefixStates,"_actionIdxLbl", sep="")
	dat<-sqlQuery(con,sql,max=maxRows,stringsAsFactors=FALSE)
	while(!is.numeric(dat)) {
		if (nrow(dat)==0) break
		dat<-as.character(t(dat))
		writeBin(dat,fALbl)
		dat<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	}
	# create fACost
	sql<-paste("SELECT w0,w1,w2 FROM ",tblPrefixStates,"_actionWeight", sep="")
	dat<-sqlQuery(con,sql,max=maxRows,stringsAsFactors=FALSE)
	while(!is.numeric(dat)) {
		if (nrow(dat)==0) break
		dat<-as.numeric(t(dat))
		writeBin(dat,fACost)
		dat<-sqlGetResults(con,max=maxRows,stringsAsFactors=FALSE)
	}
	# create fACostLbl
	writeBin(as.character(c("Time","Reward","Yield")), fACostLbl)
	# create fTransP
	sql<-paste("SELECT MAX(sId) FROM ",tblPrefixActions,"_actionIdx", sep="")
	sIdMax<-sqlQuery(con,sql)[1,1]
	sId<-0
	sql<-paste("SELECT sId,pr FROM ",tblPrefixActions,"_actionIdx WHERE sId BETWEEN ",sId," AND ",sId+100, sep="")
	dat<-sqlQuery(con,sql,stringsAsFactors=FALSE)
	while(TRUE) {
		vec<-NULL
		daply(dat, .(sId), .fun = function(df) {    # return the number of transitions for a given sId
			vec<<-c(vec,df$pr,-1)
			nrow(df)
		})
		writeBin(as.numeric(vec),fTransP)
		sId<-sId+101
		if (sId>sIdMax) break;
		sql<-paste("SELECT sId,pr FROM ",tblPrefixActions,"_actionIdx WHERE sId BETWEEN ",sId," AND ",sId+100, sep="")
		dat<-sqlQuery(con,sql,stringsAsFactors=FALSE)
	}
	close(con)
	close(fS)
	close(fSLbl)
	close(fA)
	close(fALbl)
	close(fACost)
	close(fACostLbl)
	close(fTransP)
#    readBin(binNames[3], integer(), 200000)
#    readBin(binNames[6], character(), 200000)
#    readBin(binNames[7], numeric(), 200000)
	cat("Binary files created. Names:\n")
	print(binNames)
	invisible(NULL)
})


#########################################################################/*
# @RdocMethod .deleteDairyModelDaily
#
# @title "Internal function. Remove the DairyModelDaily object from memory"
#
# \description{
#   @get "title". Should not be used except you know what you are doing.
# }
#
# @synopsis
#
# \arguments{
#  \item{p}{External pointer to the model.}
#  \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @visibility private
#
#*/#########################################################################
setMethodS3(".deleteDairyModelDaily", "HMDPDaily", function(this, p, ...) {
	.Call("DAIRY_deleteDairyModelDaily", p);
	invisible()
}, private=TRUE)


#########################################################################/**
# @RdocMethod createObjectC
#
# @title "Generate the binary files describing the HMDP model using C++ code"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{prefix}{A character string with the prefix added to the binary files. Used to identify a specific model.}
#  \item{...}{Arguments passed to \code{calcPrM} and \code{csvPrM}.}
# }
#
# \details{If the transition probabilities prM have not been calculated they are calculated and saved to a csv file. }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("createObjectC", "HMDPDaily", function(this, prefix="", saveCsv=FALSE, ...) {
	ptm <- proc.time()
	# calc number of states
	vS1<-rep(NA,this$maxLac+2)  # number of level 1 states (stage i stored in index i+1)
	for (d1 in 0:(this$maxLac+1)) vS1[d1+1]<-this$.getLev1States(d1)
	vS2<-rep(NA,this$maxDfc+1)  # number of level 2 states (stage i stored in index i+1)
	for (d2 in 0:this$maxDfc) vS2[d2+1]<-this$.getLev2States(d2)
	prICDry = rep(0,this$maxLac+1)  #IC pr in dry period (lac i stored in index i+1)
	for (d1 in 1:this$maxLac) prICDry[d1+1] = this$probInvol(400,lac=d1,preg=TRUE,dryoff=TRUE)
	if (is.null(this$.prM2A)) this$calcPrM2A(...)
	if (is.null(this$.prM)) this$calcPrM(...)
	if (saveCsv) this$csvPrM(...)
	iAZero<-this$.getIdxA(0)

	# calc yield. yield given lac, dfc and iM
	yield<-array(list(0),c(length(this$oCowY$.ssmModels)+1,this$maxDfc+1,this$stateVar$sizeM))  # create 3-dim array yield[lac+1,dfc+1,iM+1]
	for (lac in 1:(length(this$oCowY$.ssmModels)))  # yield for high lactations numbers the same
		for (dfc in 1:this$maxDfc)
			for (iM in 1:this$stateVar$sizeM-1) {
				m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
				yield[[lac+1,dfc+1,iM+1]] = this$oCowY$dailyYieldECM(m, lac, dfc)
			}

	# calc array prIC. IC pr given lac, dfc and preg (0 = not preg) is stored in prIC[lac+1][dfc+1][preg+1]
	prIC<-array(list(0),c(this$maxLac+1,this$maxDfc+1,2))
	for (lac in 1:this$maxLac)
		for (dfc in 1:this$maxDfc)
			for (preg in 0:1)
				prIC[[lac+1,dfc+1,preg+1]] <- this$probInvol(dfc,lac,(preg==1))

	# calc array prPregT. pr of positive pregnancy test (prPreg[lac+1][dfc+1]).
	prPregT<-array(list(0),c(this$maxLac+1,this$maxDfc+1))
	for (lac in 1:this$maxLac)
		for (dfc in 1:(this$oOestus$insemFinish+this$oOestus$pregTestLth))
			prPregT[[lac+1,dfc+1]] = this$oOestus$probPregTest(lac, dfc+1)

	# calc vector vA2M which assign iA to iM (stored in vA2M[iA+1])
	vA2M<-rep(-1,this$stateVar$sizeA)
	for (i in 1:this$stateVar$sizeA) vA2M[i] <- this$.mapM1ToM(this$stateVar$sA[i,'center'])

	binNames<-c("stateIdx.bin","stateIdxLbl.bin","actionIdx.bin",
		"actionIdxLbl.bin","actionWeight.bin","actionWeightLbl.bin","transProb.bin")
	binNames<-paste(prefix,binNames,sep="")

	this$ptr<-.Call("DAIRY_newDairyModelDaily",
		binNames,
		as.integer(this$oOestus$dryPeriodLth),
		as.integer(this$oOestus$gestLth),
		as.numeric(this$oOestus$probCalf),
		as.numeric(this$oOestus$probBullCalf),
		as.numeric(this$priceCalf),
		as.numeric(this$priceBullCalf),
		as.numeric(this$priceHeifer),
		as.numeric(this$priceSFU),
		as.numeric(this$priceECM),
		as.numeric(this$priceCarcassKg),
		as.integer(this$oOestus$insemStart),
		as.integer(this$oOestus$insemFinish),
		as.integer(this$oOestus$pregTestLth),
		as.integer(vS1),
		as.integer(vS2),
		as.numeric(prICDry),
		as.list(prIC),
		as.list(this$.prM),
		as.list(this$.prM2A),
		as.list(prPregT),
		as.list(yield),
		as.integer(vA2M),
		as.integer(iAZero),
		this$.deleteDairyModelDaily
	)

	cat("\nObject created.\n")
	print(proc.time() - ptm)

	invisible()
})


#########################################################################/**
# @RdocMethod genStatesBinaryC
#
# @title "Generate the binary files describing the states in the HMDP model using C++ code"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{prefix}{A character string with the prefix added to the binary files. Used to identify a specific model.}
#  \item{...}{Arguments passed to \code{calcPrM} and \code{csvPrM}.}
# }
#
# \details{If the transition probabilities prM have not been calculated they are calculated and saved to a csv file. }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("genStatesBinaryC", "HMDPDaily", function(this, ...) {
	if (is.null(this$ptr)) stop("The object have not been created in memory.")
	ptm <- proc.time()
	.Call("DAIRY_genStates",this$ptr)
	cat("\nBinary files for states created.\n")
	print(proc.time() - ptm)
	invisible()
})


#########################################################################/**
# @RdocMethod genActionsBinaryC
#
# @title "Generate the binary files describing the actions in the HMDP model using C++ code"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{prefix}{A character string with the prefix added to the binary files. Used to identify a specific model.}
#  \item{...}{Arguments passed to \code{calcPrM} and \code{csvPrM}.}
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/HMDPDaily.Rex"
#
#*/#########################################################################
setMethodS3("genActionsBinaryC", "HMDPDaily", function(this, ...) {
	if (is.null(this$ptr)) stop("The object have not been created in memory.")
	ptm <- proc.time()
	.Call("DAIRY_genActions",this$ptr)
	cat("\nBinary files for actions created.\n")
	print(proc.time() - ptm)
	invisible()
})
