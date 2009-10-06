dis<-Discretize()

## univariate 
mu<-0
sigma2<-3
n<-10
dis$discretize1DUnifEqLth(mu,sigma2,n)  # eq length intervals
dis$discretize1DUnifEqProb(mu,sigma2,n) # equal probability

## bivariate (nonuniform)
mu<-c(0,0)
sigma<-matrix(c(2,-2,-2,4),nrow=2)

par(mfrow=c(2,3))
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes)
dis$plotHypercubes(cubes,borders=TRUE)
cols<-1:length(cubes)
dis$plotHypercubes(cubes, colors=cols)
cubes<-dis$discretize2DNonunif(mu,sigma,maxKL=0.3)
dis$plotHypercubes(cubes)
dis$plotHypercubes(cubes, text="index")
text<-paste("i",1:length(cubes)-1)
dis$plotHypercubes(cubes, text=text)

par(mfrow=c(2,2))
sigma<-matrix(c(4,0,0,4),nrow=2)
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split1")
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split2")
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DNonunif(mu,sigma, modifyCenter="split2", split=0.5)
dis$plotHypercubes(cubes)

# split cubes further
par(mfrow=c(2,2))
cubes<-dis$discretize2DNonunif(mu,sigma)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=4)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=19)
dis$plotHypercubes(cubes, text="index")
cubes<-dis$splitCube2D(cubes,mu,sigma,iM=13)
dis$plotHypercubes(cubes, text="index")


## bivariate (uniform)
mu<-c(0,0)
sigma<-matrix(c(4,0,0,4),nrow=2)

par(mfrow=c(2,2))
cubes<-dis$discretize2DUnifEqInv(mu,sigma,4,6)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqInv(mu,sigma,10,10)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqProb(mu,sigma,4,6)
dis$plotHypercubes(cubes)
cubes<-dis$discretize2DUnifEqProb(mu,sigma,10,10)
dis$plotHypercubes(cubes)






