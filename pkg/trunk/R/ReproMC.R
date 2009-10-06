###########################################################################/**
# @RdocClass ReproMC
#
# @title "ReproMC class"
#
# \description{
#  Containing all biological functions related to the continious markov chain modelling the reproduction of the cow.
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
# @examples "../RdocFiles/ReproMC.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("ReproMC", function(gestLth=282, insemStart=35, insemFinish=250,
	pregTestLth=40, dryPeriodLth=49, probCalf=0.5, probBullCalf=0.5,
	model="reproMC2", ...)
{
	if (model=="reproMC1") mod<-reproMC1()

	extend(Object(), "ReproMC",
		gestLth=gestLth,
		insemStart=insemStart,
		insemFinish=insemFinish,
		pregTestLth=pregTestLth,
		dryPeriodLth=dryPeriodLth
	)
})


##########################################################################################/**
# @RdocMethod reproMC1
#
# @title "Constructs a Markov chain of the Oestrous cycle (insemination not included)"
#
# \description{
# Constructs a Markov Process model the Oestrous cycle, based on mean and variance
# of the length of the estrous cycle, and the mean interval between oestrogen peak and LH-peak,
# LH-peak and Ovulation and finally, between start of CL-regression to Oestrogen peak. Approximated
# a gamma-distribution by a sum of exponentially distributed substeps.
# }
# @synopsis
# \arguments{
#   \item{MeanDuration}{The mean length of the oestrous cycle}
#   \item{VarDuration}{The variance of the length of the oestrous cycle}
#   \item{OestroLHPeak}{Interval between oestroge peak and LH peak}
#   \item{LHPeakOvul}{Interval between LH-peak and ovulation time}
#   \item{CLRegOestroPeak}{Interval between start of CL-regressin and subsequent Oestrogen peak}
#   \item{...}{Not used.}
# }
# \value{
#   \item{beta}{The mean length of each substep ( the shape parameter in the gamma-distribution of oestrous cycle length)}
#   \item{scale}{The scale parameter in the gamma-distribution of oestrous cycle length. Corresponds to the
#   number of substeps in the cycle. }
#   \item{tider}{The length of each of the overall stages in the cycle}
#   \item{Stepno}{The number of substeps in the cycle. Integer valued}
#   \item{tider}{The length of each of the overall stages in the cycle}
#   \item{rate}{The marix of transition intensities for the Marko Process}
#     \item{kOestrousState}{The state number corresponding to the state where ovulation occurs}
#     \item{kImplantationState}{The state number corresponding to the state where Implantation occurs (start of CL-regression)}
#     \item{type}{type of object = 'Oestrous'}
# }
# \author Erik Jørgensen
#
#*/##########################################################################################
setMethodS3("reproMC1", "ReproMC", function(this, MeanDuration=21, VarDuration=2,
	OestroLHPeak=0.4, LHPeakOvul=1.7, CLRegOestroPeak=5.0, ...){
 # 2-3 (Oestrogen peak to LH peak 4.2)    4 1 -   4
 # 3-4 (LH peak to Ovulation)            18 5 -  22
 # 4-1 (Ovulation to CL regression)     146 23 - 168
 # 1-2 (CL regression to Oestrogen Peak) 53 169- 221

 #  The Gamma distribution with parameters shape = beta and scale = s has density
 #  f(x)= 1/(s^beta Gamma(beta)) x^(beta-1) e^-(x/s)

 #for x >= 0, beta > 0 and s > 0. (Here Gamma(beta) is the function
 #implemented by R's gamma() and defined in its help.)

 #The mean and variance are E(X) = beta*s and Var(X) = beta*s^2.

  #The cumulative hazard H(t) = - log(1 - F(t)) is -pgamma(t, ...,
  #lower = FALSE, log = TRUE).
  beta <- VarDuration/MeanDuration  # rateparameter i skridtene
  scale <- MeanDuration/beta
  tider<-c(OestroLHPeak,LHPeakOvul,CLRegOestroPeak)
  tider<-c(tider[1:2],MeanDuration-sum(tider),tider[3])
  alpha<- tider/beta    # alpha parameter == antal states i oestrus cycle
  steps<-round(alpha+.05)
  Stepno<-sum(steps)
  kOestrousState     <-sum(steps[1:2])
  kImplantationState <-sum(steps[1:3])
  rate <- matrix(rep(0,Stepno^2),nrow=Stepno)
  for (i in (1:Stepno-1)){
	  rate[i,i+1]<-  1/beta
  }
  rate[Stepno,1]<-1/beta
  diag(rate)<- -rowSums(rate)

  list(beta=beta,scale=scale,tider=tider,Stepno=Stepno,rate=rate,
	  kOestrousState=kOestrousState,kImplantationState=kImplantationState,
	  type='Oestrous')
})


##########################################################################################/**
# @RdocMethod reproMC2
#
# @title "Constructs a Markov Process model for the Oestrous cycle with insemination included"
#
# \description{Constructs a Markov Process model for the Oestrous cycle with mating included,
# based on the model reproMC1 and additional parameters concerning probability of mating per
# cycle, and probability of pregnancy per mating. Duplicates part of the reproduction cycle to
# handle states with and without mating.
# }
# @synopsis
#
# \arguments{
#   \item{mod}{The reproMC1 model to be expanded.}
#   \item{pMate}{The probability of (successful) mating per oestrus}
#   \item{pPreg}{The probability of pregnancy per (successful) mating }
#   \item{...}{Not used.}
# }
#
# \value{
#   \item{beta}{The mean length of each substep ( the shape parameter in the gamma-distribution of oestrous cycle length)}
#   \item{scale}{The scale parameter in the gamma-distribution of oestrous cycle length. Corresponds to the number of substeps in the cycle. }
#   \item{tider}{The length of each of the overall stages in the cycle}
#   \item{Stepno}{The number of substeps in the basic oestrous cycle. Integer valued}
#   \item{AugStepno}{The number of substeps in full cycle with mating included. Integer valued}
#   \item{tider}{The length of each of the overall stages in the cycle}
#   \item{rate}{The matrix of transition intensities for the Marko Process}
#   \item{kOestrousState}{The state number corresponding to the state where ovulation occurs}
#   \item{kImplantationState}{The state number corresponding to the state where Implantation occurs (start of CL-regression)}
#   \item{kPregnantState}{The state number corresponding to the (absorbing) pregnancy state}
#   \item{kMatedImplantationState}{The state number corresponding to the state where Implantation occurs in mated subset of states.}
#   \item{kMatedState}{The state number corresponding to the first state in the mated subsets of states}
#   \item{type}{type of object = 'MateOestrous'}
# }
# @author
#
#*/##########################################################################################
setMethodS3("reproMC2", "ReproMC", function(this, mod, pMate=0.85, pPreg=0.70, ...){
	AugStepno=mod$Stepno+(mod$kImplantationState-mod$kOestrousState)+1
	kMatedState<- mod$Stepno+1
	kPregnantState<- AugStepno
	kMatedImplantationState<- AugStepno-1
	rate <- matrix(rep(0,AugStepno^2),nrow=AugStepno)
	for (i in (1:AugStepno-1)){
		rate[i,i+1]<-1/mod$beta
	}
	rate[mod$Stepno,mod$Stepno+1]<-0
	rate[mod$Stepno,1]<- 1/mod$beta
	# Mating transitions
	rate[mod$kOestrousState,mod$kOestrousState+1]<- (1-pMate)/mod$beta
	rate[mod$kOestrousState,kMatedState]     <-   pMate/mod$beta
	#implantation transitions
	rate[kMatedImplantationState,kPregnantState]  <- pPreg/mod$beta
	rate[kMatedImplantationState,mod$kImplantationState+1] <- (1-pPreg)/mod$beta
	diag(rate)<- -rowSums(rate)

	list(beta=mod$beta,scale=mod$scale,tider=mod$tider,Stepno=mod$Stepno,
	AugStepno=AugStepno,
	rate=rate,
	kOestrousState     = mod$kOestrousState,
	kImplantationState = mod$kImplantationState,
	kPregnantState     = kPregnantState,
	kMatedImplantationState = kMatedImplantationState,
	kMatedState =kMatedState,
	type='MateOestrous')
})


##########################################################################################/**
# @RdocMethod reproMC3
#
# @title "Constructs a Markov Process model for the Oestrous cycle with insemination as well as early gestation included"
# \description{
# Constructs a Markov Process model for the Oestrous cycle with insemination and early gestation incuded.
# Based on the model reproMC1 and additional parameters concerning probability of mating per
# cycle, probability of pregnancy per mating, and mean and duration of interval from implantation to
# positive pregnancy test. Duplicates part of the reproduction cycle to
# handle states with and without mating.
# }
# @synopsis
# \arguments{
#   \item{mod}{The reproMC1 model to be expanded.}
#   \item{pMate}{The probability of (successful) mating per oestrus}
#   \item{pPreg}{The probability of pregnancy per (successful) mating }
#   \item{meanDurationToTest}{Mean duration of interval from implantation to positive pregnancytest}
#   \item{varDurationToTest}{Variance of duration of interval from implantation to positive pregnancytes}
#   \item{...}{Not used.}
# }
# \details{}
# \value{
#
#   \item{beta}{The mean length of each substep ( the shape parameter in the gamma-distribution of oestrous cycle length)}
#   \item{scale}{The scale parameter in the gamma-distribution of oestrous cycle length. Corresponds to the
#   number of substeps in the cycle. }
#   \item{tider}{The length of each of the overall stages in the cycle}
#   \item{Stepno}{The number of substeps in the basic oestrous cycle. Integer valued}
#   \item{AugStepno}{The number of substeps in full cycle with mating included. Integer valued}
#   \item{rate}{The matrix of transition intensities for the Markov Process}
#     \item{kOestrousState}{The state number corresponding to the state where ovulation occurs}
#     \item{kImplantationState}{The state number corresponding to the state where Implantation occurs (start of CL-regression)}
#    \item{kPregnantState}{The state number corresponding to the (absorbing) pregnancy state}
#    \item{kMatedImplantationState}{The state number corresponding to the state where Implantation occurs
#    in mated subset of states.}
#    \item{kPregtestState}{The state number corresponding to the state where
#    positive pregnancy test occurs}
#    \item{kMatedState}{The state number corresponding to the first state in the mated subsets of states}
#      \item{type}{type of object = 'PregTestMateOestrous'}
# }
# @author
#
#*/##########################################################################################
setMethodS3("reproMC3", "ReproMC", function(this, mod, pMate=0.85, pPreg=0.70,
	meanDurationToTest=21, varDurationToTest=1, ...)
{
	betaToTest<- varDurationToTest/meanDurationToTest
	StepsToTest<-meanDurationToTest/betaToTest
	PregTestStepno <- trunc(StepsToTest) # evt +1
	OestrousAndPregStepno=mod$Stepno+(mod$kImplantationState-mod$kOestrousState)+1
	AugStepno<-OestrousAndPregStepno+PregTestStepno

	kMatedState<- mod$Stepno+1
	kPregnantState<- OestrousAndPregStepno
	kMatedImplantationState<- OestrousAndPregStepno-1
	kPregtestState<-AugStepno
	rate <- Matrix(rep(0,AugStepno^2),nrow=AugStepno)

	# Fyld ud for Oestrous cycle and mated but not pregnant
	for (i in (1:(OestrousAndPregStepno-1))){
		rate[i,i+1]<-  1/mod$beta
	}
	rate[mod$Stepno,mod$Stepno+1]<-0
	rate[mod$Stepno,1]<- 1/mod$beta
	# Mating transitions
	rate[mod$kOestrousState,mod$kOestrousState+1]<- (1-pMate)/mod$beta
	rate[mod$kOestrousState,kMatedState]     <-   pMate/mod$beta
	#implantation transitions
	rate[kMatedImplantationState,kPregnantState]  <- pPreg/mod$beta
	rate[kMatedImplantationState,mod$kImplantationState+1] <- (1-pPreg)/mod$beta
	# Fyld ud for Pregnant to Positive Test
	for (i in (kPregnantState:(kPregtestState-1))){
		rate[i,i+1]<-  1/betaToTest
	}
	diag(rate)<- -rowSums(rate)

	list(beta=mod$beta,scale=mod$scale,tider=mod$tider,Stepno=mod$Stepno,
		AugStepno=AugStepno,
		rate=rate,
		kOestrousState     =mod$kOestrousState,
		kImplantationState =mod$kImplantationState,
		kPregnantState     =kPregnantState,
		kMatedImplantationState =kMatedImplantationState,
		kMatedState =kMatedState,
		kPregtestState = kPregtestState,
		OestrusAndPregStepno=OestrousAndPregStepno,
		type='PregTestMateOestrous')
})


##########################################################################################/**
# @RdocMethod daysToPregTest
# @title "Calculates the distribution of days to Postive Pregtest"
# \description{
# Calculates the distribution of days to positive pregnancy test
# based on the underlying markov process
# of the oestrous cycle augmented with states corresponding to mating as well as the first part of the
# gestation period. Based on caculation
# of the distribution of time to First Passage of the pregnancy test state.
# }
# @synopsis
# \arguments{
#   \item{pMate}{The probability of (successful) mating per oestrus}
#   \item{pPreg}{The probability of pregnancy per (successful) mating }
#   \item{MatingStart}{Start of mating (Voluntary Waiting Period), days from calving }
#   \item{MatingEnd}{End of mating period, days from calving}
#   \item{CalcEnd}{end of calculations, days from calving}
#   \item{steplength}{The steplength (in days) used in the calculations}
#   \item{meanDurationToTest}{Mean duration of interval from implantation to positive pregnancytest}
#   \item{varDurationToTest}{Variance of duration of interval from implantation to positive pregnancytes}
#   \item{...}{Arguments passed to reproMC1.}
#  }
# \details{
# }
# \value{
#   The function returns a list
#   \item{dfc}{A vector of days from calving from \code{Matingstart} to \code{CalcEnd}}
#   \item{p}{A vector of accumulated probabilities of postive pregnancy test}
#   \item{dist}{A vector of probability of positive pregnancy test for each day}
#   \item{margp}{A vector of marginal probabilities of positve pregnancy test for each day, i.e.,
#            the conditional of a positive pregnancy test at day dfc, given no positive pregnancy before dfc}
#   }
# @author
#
#*/##########################################################################################
setMethodS3("daysToPregTest", "ReproMC", function(this, pMate = 0.45, pPreg = 0.50,
	MatingStart = 35, MatingEnd = 170,
	CalcEnd = 220, steplength = 1,
	meanDurationToTest=35,varDurationToTest=1, ...)
{
	mod <- this$reproMC1(...)
	mod <- this$reproMC3(mod, pMate = pMate,
		pPreg = pPreg, meanDurationToTest=meanDurationToTest,varDurationToTest=varDurationToTest)
	steps <- c(1:(mod$AugStepno - 1))
	mating <- rep(FALSE, mod$AugStepno - 1)
	steps[mod$kMatedState:mod$kMatedImplantationState] <- c(23:168)
	mating[mod$kMatedState:mod$kMatedImplantationState] <- TRUE
	w <- rep(1/mod$Stepno, mod$AugStepno)
	w[mod$kMatedState:mod$AugStepno] <- 0

#   First period: Maing period
	FirstDuration <- MatingEnd - MatingStart + 1
	FirstPassage <- this$FirstPassageDist(mod$rate, w, mod$AugStepno,
		MatingEnd - MatingStart + 1, steplength)

#   Second Period: Mating Ended
#   Collapse transiton rates
	newrate <- mod$rate
	diag(newrate) = 0
	with(mod,newrate[kOestrousState, kOestrousState + 1] <-
	   newrate[kOestrousState, kOestrousState + 1] +
	   newrate[kOestrousState, kMatedState])
	newrate[mod$kOestrousState, mod$kMatedState] <- 0
	diag(newrate) <- -rowSums(newrate)

	SecondDuration <- CalcEnd - MatingEnd + 1

	FirstPassageAfterMating <- this$FirstPassageDist(newrate, FirstPassage$statedist,
		mod$AugStepno, CalcEnd - MatingEnd + 1, steplength)

#   Concatenate results
	tid <- c(MatingStart:(CalcEnd + 1))
	pvalues <- (c(FirstPassage$p, FirstPassageAfterMating$p *
		(1 - FirstPassage$p[FirstDuration]) + FirstPassage$p[FirstDuration]))

	dist <- pvalues[-1]-pvalues[-length(pvalues)]
   margp <- dist/(1-pvalues[-length(pvalues)])
	list(dfc = tid, p = pvalues, dist=c(0,dist),margp=c(0,margp) )
})


##########################################################################################/**
# @RdocMethod FirstPassageDist
# @title "Calculates distribution of first-passage-times in a continous Markov-Process"
# \description{
# Calculates the approximate distribution of first passage times of a continuous Markov-Process
# using time-steps of discrete length
# }
# @synopsis
# \arguments{
#   \item{rate}{Matrix of transition intensities }
#   \item{prior}{Prior distribution over states }
#   \item{absorbstate}{The state number of the first-passage}
#   \item{Interval}{Length of the time interval over which to calculate the distribution }
#   \item{timestep}{Length of discrete time interval used for the approximation }
#   \item{...}{Not used.}
# }
# \details{
#
# }
# \value{ returns a vector of length (round(Interval/timestep)-1)
#    with the probability of passage in time interval i}
# @author
#
#*/##########################################################################################
setMethodS3("FirstPassageDist", "ReproMC", function(this, rate, prior, absorbstate,
	Interval=10, timestep=1, ...)
{
	Stages<-round(Interval/timestep)
	statedist<- prior
	newrate<-rate
	dim(statedist)<-c(1,ncol(rate))
	newrate[absorbstate,]<-rep(0,ncol(rate))
	drate<-this$calcTransMat(newrate,steplength=timestep,n=20)
	statedistnew<-Matrix(statedist)
	statedistnew[1,absorbstate]<-0
	statedistnew<-statedistnew/sum(statedistnew)
	res<-rep(0,Stages)
	time<-c(1:Stages)*timestep
	for (i in (1:Stages)){
	  statedistnew <- crossprod(t(statedistnew),drate)
		res[i]<-statedistnew[1,absorbstate]
	}
	# NB problem med første (sidste) trin ! Sættes til NA overvej at beregne et ekstra trin
	list(time=time,p=res,d=c((res[-1]-res[1:(Stages-1)]),NA),statedist=statedistnew)
})


##########################################################################################/**
# @RdocMethod calcTransMat
# @title "Calculates the transition matrix P(t) for a continuous Markov Process for a finite time-interval t"
# \description{
#   @get "title".
#   Can call different packages for calculation of matrix exponentials
# }
#
# @synopsis
#
# \arguments{
#   \item{Q}{Transition rate matrix of transition intensities in the continuous Markov process. }
#   \item{t}{Length of the time interval}
#   \item{n}{Number of iterations e.g. the length of the minimum interval is \code{t*2**(-n)}, if \code{method=='robust'}}
#   \item{method}{Method used for calculation, \code{method} should be either \code{'robust','msm','Matrix','expm'} see details}
#   \item{...}{Not used.}
# }
#
# \details{
# The methods used for calculation of the matrix exponential is either
#   \item{robust}{divides \code{t} into 2^n intervals and calculates expm(rate/(2^n))^(2^n)  }
#   \item{Matrix}{current default. Uses the \code{expm} function in the \code{Matrix} package. Should be able to utilize the sparseness of the matrix.}
# }
#
# \value{Returns the transition matrix.}
#
# @author
#
#*/##########################################################################################
setMethodS3("calcTransMat", "ReproMC", function(this, Q, t=1, n=10, method='Matrix', ...){
	if (method=='robust')
	{
		t<-t/2^n
		trans<-this$.calcTransMatShort(Q,t)
		for (i in (1:n)){
			trans<-trans%*%trans
		} # end for
	} else {
		trans<-Matrix::expm(Matrix(Q*t))
	}
	trans
})


##########################################################################################/*
# @RdocMethod .calcTransMatShort
# @title "Calculates the transition matrix P(t) for a continuous Markov Process for a short time-interval t"
# \description{
#  Calculates the transition matrix based on transition intensity based on a short time-interval.
#  Note the time-interval must be short since does not calculate multiple transitions.
# }
# @synopsis
# \arguments{
#  \item{Q}{Matrix of transition intensities for continuous Markov Chains }
#  \item{t}{time interval }
#  \item{...}{Not used.}
# }
#
# \value{Returns the transition matrix.}
#
# @author
#
#*/##########################################################################################
setMethodS3(".calcTransMatShort", "ReproMC", function(this, Q, t, ...){
	# Optimalt set skal det være Overgangsmatricen for t
	tmat <- 1-exp(-Q*t)
	diag(tmat)<-0
	diag(tmat)<-1-rowSums(tmat)
	tmat
})


##########################################################################################/**
# @RdocMethod meanProgest
# @title "Calculates mean progesterone level for each state in the Markov Process"
#
# \description{
#  Calculates mean progesterone level for each state in the Markov chain defined by \code{\link{ReproMC1}}.
#  The progesterone curve is calculated as two exponential curves joining at mu[2].
# }
# @synopsis
#\arguments{
#   \item{x}{A vector of state numbers in (1,...,maxState) }
#   \item{mu}{Vector of size 3 containing c(state where obvolution, state where CL regression start, state of oestrogen peak).
#   \item{alpha}{A parameter describing the shape of the level (2 = Gaussian dist/error function). }
#   \item{p}{The relative level at the extreme stateno}
#   \item{...}{Not used.}
#}
#\value{A vector of length Maxstate with the mean progesterone level for each state.}
# @author
#
#*/##########################################################################################
setMethodS3("meanProgest", "ReproMC", function(this, x, mu=c(22,168,221), alpha=2, p=0.02, ...){
	sleft=-abs(mu[1]-mu[2])^alpha/log(p)
	sright=-abs(mu[3]-mu[2])^alpha/log(p)
	s<-rep(sleft,length(x))
	s[x>mu[2]]<-sright
	exp(-abs(x-mu[2])^alpha/s)
})


##########################################################################################/**
# @RdocMethod simDisMC
# @title "Simulates a discrete Markov chain based on the transition matrix"
#
# \description{
#  Calculates a sequence of states over a number of stages.
# }
# @synopsis
#\arguments{
#   \item{P}{The matrix of transition probabilities.}
#   \item{stages}{Number of stages to simulate. }
#   \item{initState}{The initial state.}
#   \item{...}{Not used.}
#}
#\value{Sequence of states starting with the initial state.}
#
# @author
#
#*/##########################################################################################
setMethodS3("simDisMC", "ReproMC", function(this, P, stages, initState=1, ...){
	lmat<-P
	states<-nrow(lmat)
	#diag(lmat)<-rep(0,nrow(lmat))   # hvorfor det??
	events<-rep(NA,stages)
	events[1]<-initState
	for (i in (2:stages)){
		events[i]=sample(states,1,prob=lmat[events[i-1],])
	}
	events
})


##########################################################################################/**
# @RdocMethod posteriorMix
# @title "Calculates the posterior distribution of a normal mixture model "
#
# \description{
#  Calculates the posterior distribution of a normal mixture model, based on a vector of
#  observations from observations from the same state. Based on the \code{\link{nor1Mix}}
#  package.
# }
# @synopsis
#\arguments{
#   \item{obj}{The object containing the normal mixture }
#   \item{obs}{The vector of observations }
#   \item{...}{Not used.}
#}
#\value{A vector with posterior probabilities.}
#
# @author
#
#*/##########################################################################################
setMethodS3("posteriorMix", "ReproMC", function(this, obj, obs, ...){
  mx<-mean(obs)
  nx<-length(obs)
  prior<- obj[,3]
  likely<-dnorm(x=mx,obj[,1],sqrt(obj[,2]/nx))
  posterior <-prior*likely
  posterior/sum(posterior)
})
