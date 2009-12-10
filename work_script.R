# R script file
# Author: Lars Relund
# Description: Preliminary script for DairyCowYield package

# To use the tree view in WinEdt keep the following structure. Then chunks/blocks will be
# visualized in the tree together with bookmarks/comments



### CHUNK: Source function for all files in the package R dir
# ---------------------------------------------------------------------------------------------

srcIt <- function(){
	f  <- list.files("pkg/trunk/R",pattern="\\.[R]",full.names=TRUE)
	#fl <- f[(1:length(f))[-grep("~",f)]]
	a<-lapply(f,source)
}
srcIt()

#srcit <-  function(){
#  flist  <- list.files("DairyCowYield/",pattern="\\.[R]",full.names=TRUE)
#  flist  <- flist[(1:length(flist))[-grep("~",flist)]]
#  a<-lapply(flist, source)
#}

# useful commands
#system("Rcmd check DairyCowYield",invisible = FALSE,wait=FALSE)
#system("Rcmd build DairyCowYield")
#system("Rcmd build --binary DairyCowYield")
#system("Rcmd INSTALL DairyCowYield")
#system("Rcmd Rd2dvi DairyCowYield")   # create a dvi file with doc
#prompt(mdpModel)   # create doc file for object
### -------------------------------------------------------------------------------------------




### -------------------------------------------------------------------------------------------
## Playing with dairy package




library(dairy)
mod<-HMDPDaily(maxLac=2,csvPrM="prM.csv")
this<-mod
lac<-1
dfc<-273
iM=0
iDry=0
age<-this$cowAge(lac)
con <- this$.getDOC(iDry)
con
this$.getLev2States(dfc)
this$.getLev2States(dfc+1)
this$.getMinDryWeekIdx(dfc)
this$transPrKeep(lac, dfc, iM, iDry)


mod$oCowWF$BCS(dfc,con)
m <- as.matrix(this$stateVar$sM[[iM+1]]$center)
milkECM = this$oCowY$dailyYieldECM(m, lac, dfc)
reward = this$rewardMilk(lac, dfc, iM);
this$oCowWF$energy(this$cowAge(lac), dfc, con, milkECM)[6]
gest=282
this$oCowWF$fetusW((gest+con)-dfc)
this$oCowWF$.sfuFetus







## Compare binary vectors
library(dairy)
binNames<-c("stateIdx.bin","stateIdxLbl.bin","actionIdx.bin",
	"actionIdxLbl.bin","actionWeight.bin","actionWeightLbl.bin","transProb.bin")
binNamesR<-paste("genR_",binNames,sep="")
binNamesC<-paste("genC_",binNames,sep="")

# compare actionIdx/stateIdx vectors (integer())
i<-1    # 1,3
vR<-readBin(binNamesR[i], integer(), 200000000)
vC<-readBin(binNamesC[i], integer(), 200000000)
idx<-which(vR!=vC)[1]
idx
vR[c((idx-20):(idx+10))]
vC[c((idx-20):(idx+10))]

i<-3    # 1,3
vR<-readBin(binNamesR[i], integer(), 200000000)
vC<-readBin(binNamesC[i], integer(), 200000000)
idx<-which(vR!=vC)[2]
idx
vR[c((idx-20):(idx+10))]
vC[c((idx-20):(idx+10))]

Fejl i action "dry" til state "0,0" (iM=0, iDry=31, dfc=273, id=135894, lac=1)



dfR1<-actionIdxDf("genR_")
dfR2<-stateIdxDf("genR_")
subset(dfR1, sId==2870)
subset(dfR2, sId==562273)
subset(dfR2, sId==7)

# compare label vectors (character())
i<-6    # 2,4,6
vR<-readBin(binNamesR[i], character(), 20000000)
vC<-readBin(binNamesC[i], character(), 20000000)
idx<-which(vR!=vC)[1]
idx
vR[c((idx-20):(idx+10))]
vC[c((idx-20):(idx+10))]

# compare flt vectors (numeric())
i<-7    # 5,7
vR<-readBin(binNamesR[i], numeric(), 20000000)
vC<-readBin(binNamesC[i], numeric(), 20000000)
idx<-which(abs(vR-vC) > 0.0000000001)[1]
idx
vR[c((max(idx-20,0)):(idx+10))]
vC[c((max(idx-20,0)):(idx+10))]

## Try to compare df's
library(dairy)

dfR<-transProbMat("genR_")
dfC<-transProbMat("genC_")
tmp<-which(abs(dfR-dfC)> 0.0000000001,arr.ind = T)

dfR<-stateIdxDf("genR_")
dfC<-stateIdxDf("genC_")
dfR==dfC

dfR<-actionIdxMat("genR_")
dfC<-actionIdxMat("genC_")
tmp<-which(dfR!=dfC,arr.ind = T)
dfR[tmp]
dfC[tmp]

dfR<-actionInfo("genR_")
dfC<-actionInfo("genC_")
tmp<-which(abs(dfR-dfC)> 0.0000000001,arr.ind = T)

tmp <- dfR==dfC

dfR<-actionWeightMat("genR_")
dfC<-actionWeightMat("genC_")
tmp<-dfR!=dfC
dim(tmp)
idx<-apply(tmp,1,any)
head(dfR[idx,])
head(dfC[idx,])
tmp<- abs(dfR-dfC) > 0.0000001
idx<-apply(tmp,1,any)
head(dfR[idx,])
head(dfC[idx,])
dfR1<-actionIdxDf("genR_")
dfR2<-stateIdxDf("genR_")
subset(dfR1, aId==4252)
subset(dfR2, sId==2152)


head(which(dfR[idx,]-dfC[idx,] > 0.0000001))

dfR[which(dfR[idx,]-dfC[idx,] > 0.0001)]








## test C function
a<-array(list(),c(2,2,2))   # array where each element contain a numeric vector
a[[1,1,1]]<-rnorm(1)
a[[1,1,2]]<-rnorm(3)
a[[2,1,1]]<-rnorm(2)
a[[2,1,2]]<-rnorm(3)
a[[1,2,1]]<-rnorm(4)
a[[1,2,2]]<-rnorm(3)
a[[2,2,1]]<-rnorm(5)
a[[2,2,2]]<-rnorm(1)
d = dim(a)
d
i<-1; j<-1; k<-2
q<-i + j*d[1] + k*(d[1]*d[2]) - (d[1] + d[1]*d[2])
a[[i,j,k]]      # access an element
a[[q]]          # access the same element using a single index
.Call("DAIRY_PrintValue", as.list(a))



library(dairy)
#mod<-HMDPDaily()
mod<-HMDPDaily(maxLac=2,oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
mod$getFields()

mod$createObjectC(prefix="test_")
mod$genStatesBinaryC()
mod$genActionsBinaryC()

binNames<-c("stateIdx.bin","stateIdxLbl.bin","actionIdx.bin",
	"actionIdxLbl.bin","actionWeight.bin","actionWeightLbl.bin","transProb.bin")
binNames<-paste("test_",binNames,sep="")

readBin(binNames[1], integer(), 20)
readBin(binNames[2], character(), 20)
readBin(binNames[3], numeric(), 20)
readBin(binNames[4], character(), 20)
readBin(binNames[5], integer(), 20)
readBin(binNames[6], character(), 20)
readBin(binNames[7], numeric(), 20)








dis<-Discretize()
dis$plotHypercubes(mod$stateVar$sM)

mod<-HMDPDaily(maxLac=2,oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)

dis<-Discretize()
dis$plotHypercubes(mod$stateVar$sM)



## calculate pr(m_{t+1}|m_{t})
mod$calcPrM()   # calc all prM which are stored internally
dim(mod$.prM)
pr<-mod$prM(lac=2,dfc=50,idxM=3)    # get pr(m_{51}|m_{50}) where the idx of m_{50} is 3
pr
sum(pr[,2])
mod$csvPrM()   # save as csv file
#unlink("prM.csv")


## calculate transition pr from level 2 to 1
mod$calcPrM2A()   # calc all prM2A which are stored internally
mod$prM2A(idxM=0)
mod$prM2A(idxM=10)
mod$csvPrM2A()   # save as csv file
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", csvPrM2A="prM2A.csv",
	oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
mod$prM2A(idxM=0)
mod$prM2A(idxM=10)


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
plot(1:12,mod$cowAge(1:12),xlab="Lactation",ylab="Age in days")


## generate a small model
source("./R/HMDPDaily.R")
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
#Rprof()
ptm <- proc.time()
mod$genHMDP(noShared=TRUE)
proc.time() - ptm
#Rprof(NULL)
#summaryRprof()

ptm <- proc.time()
mod$genHMDP(format="hmp", noShared=TRUE)
proc.time() - ptm

convertBinary2HMP()


mdp<-loadMDP()
mdp
sIdx<-stateIdxDf()
aInfo<-actionInfo()



lac=1
dfc=56
iDry=0
mod$rewardDry(lac, dfc, iDry)

## generate small model with only one idxA (such than can be converted to hmp)
mod<-HMDPDaily(maxLac=2, csvPrM="prM.csv", oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
ssmModel<-mod$oCowY$.ssmModels[[1]]
mod$stateVar$sA<-mod$oDis$discretize1DUnifEqLth(ssmModel$m0[1],ssmModel$C0[1,1],1)
mod$stateVar$sizeA<-1
ptm <- proc.time()
mod$genHMDP()
proc.time() - ptm
convertBinary2HMP()

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


states<-list()

state<-list(sId, sIdx, label)

actions<-list()
action<-list(sId, weights, probs, label)



library(dairy)
## try to generate states and actions to db
mod<-HMDPDaily(maxLac=2, oOestus=Oestus(gestLth=100,pregTestLth=10,insemStart=15, insemFinish=30), maxKL=0.9)
mod<-HMDPDaily()
mod$genStates2SQLite(tblPrefix="full")
mod$genActions2SQLite(tblPrefix="full")
mod$convertSQLite2Binary(tblPrefixActions="full")





# test binary files
fS <- file("test.bin", "wb")
con <- odbcConnect("mdp_models")
sql<-paste("SELECT 0 AS d0,0 AS s0,0 AS a0,d1,s1,a1,d2,s2,-1 AS sep FROM test_stateIdx", sep="")
dat<-sqlQuery(con,sql,max=10,stringsAsFactors=FALSE)
dat<-as.integer(t(dat))
writeBin(dat[!is.na(dat)],fS)
close(fS)
close(con)
readBin("test.bin", integer(), 200)

## calc all prM (takes approx 1 hour)
mod<-HMDPDaily()
ptm <- proc.time()
mod$calcPrM()
mod$csvPrM()
proc.time() - ptm







param<-setParameters()
ssmModel<-buildSSMYieldModel()
sv<-mdpStateVariables(ssmModel, maxKL=0.8, parameters = param)

# sv$sA indices
mapMToM1(c(-4,12))
mapMToM1(c(3,12))
mapMToM1(c(3,0))
mapMToM1(c(30,12))

# sv$sM indices
mapM1ToM(3)
srcIt()
getMIdx(3)
getIdxM(c(0,0))
getIdxM(c(0.1,15))

# conversion
dfc2Week(7)
dfc2Week(8)
dfc2DryWeek(120)
getMaxDryWeekIdx(120)
getMinDryWeekIdx(1)

getLev2States(120)
getDryWeekIdx(120,400)
getDryWeekIdx(120,4)
