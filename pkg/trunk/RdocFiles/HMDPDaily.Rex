
## model from Relund09
mod<-HMDPDaily()
mod$getFields()
mod$oDis$plotHypercubes(mod$stateVar$sM)

## smaller model with worse discretization 
mod<-HMDPDaily(maxLac=2,oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
mod$oDis$plotHypercubes(mod$stateVar$sM)

\dontrun{
## calculate pr(m_{t+1}|m_{t})
mod$calcPrM()   # calc all prM which are stored internally in the private variable .prM
dim(mod$.prM)
pr<-mod$prM(lac=2,dfc=50,iM=3)    # get pr(m_{51}|m_{50}) where the idx of m_{50} is 3
pr
sum(pr[,2])
mod$csvPrM()   # save as csv file
}


## calculate transition pr from level 2 to 1
mod$calcPrM2A()   # calc all prM2A which are stored internally in the private variable .prM2A
mod$prM2A(iM=0)
mod$prM2A(iM=10)

\dontrun{
## load small model with csv file for prM
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", csvPrM2A="prM2A.csv",
	oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
mod$prM(lac=2,dfc=50,iM=3)
}

## calculate probability for day 1 to day 399 with positive pregnancy test on day 100
h<- rep(NA,399)
pregday<-100
h[1:pregday]<-sapply(c(1:pregday),mod$probInvol,lac=2)
h[(pregday+1):399]<-sapply(c((pregday+1):399),mod$probInvol,lac=2,preg=TRUE)
cull<-mod$probInvol(400,lac=2,preg=TRUE,dryoff=TRUE)
hcumsum<-cumsum(-log(1-h))  # calculate accumulated hazard
par(mfrow=c(1,2))
plot(c(1:399),h,type='l',ylim=c(0,0.0015),xlab='Days from calving',ylab='Marginal probability of involuntary culling')
plot(c(1:399),1-exp(-hcumsum),type='l',xlab='Days from calving',ylab='Accumulated probability of involuntary culling')

## cow age at start of lactation
par(mfrow=c(1,1))
plot(1:12,mod$cowAge(1:12),xlab="Lactation",ylab="Age in days")

\dontrun{
## generate a small model
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
ptm <- proc.time()
mod$genBinary()
proc.time() - ptm

## generate small model with only one state at level one (idxA) (such than can be converted to hmp)
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
# hack the model so only one state at level one
ssmModel<-mod$oCowY$.ssmModels[[1]]
mod$stateVar$sA<-mod$oDis$discretize1DUnifEqLth(ssmModel$m0[1],ssmModel$C0[1,1],1)
mod$stateVar$sizeA<-1
# generate model
ptm <- proc.time()
mod$genBinary()
proc.time() - ptm
convertBinary2HMP()     # convert to hmp format that can be read by MLHMP

## Solve the model using MDP package that load the binary files
mdp<-loadMDP()
mdp
sIdx<-stateIdxDf()
aInfo<-actionInfo()
iW<-1
iDur<-0
rate<-0.03
rateBase<-365
policyIteDiscount(mdp, iW, iDur, rate, rateBase)
policy<-getPolicy(mdp, labels=TRUE)
policy<-merge(sIdx,policy)
policyW<-getPolicyW(mdp, iW)
#head(policyW)
policy<-merge(policy,policyW)
rpo<-calcRPO(mdp, iW, iA=0, criterion="discount", iDur=iDur, rate=rate, rateBase=rateBase)
policy<-merge(policy,rpo)
head(policy)
}











