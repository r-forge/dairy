###########################################################################/**
# @RdocClass Discretize
#
# @title "Discretize class"
#
# \description{
#  Containing all methods related to discritizing an uni/bi-variate normal distribution. Both uniform and nonuniform discretization possible.
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods ""
# }
#
#
# @examples "../RdocFiles/Discretize.Rex"
#
# \references{
# [1] Nielsen, L.R.; Jørgensen, E. & Højsgaard, S. Embedding a state space model into a Markov decision process Dept. of Genetics and Biotechnology, Aarhus University, 2008. \cr
# }
#
# @author
#*/###########################################################################
setConstructorS3("Discretize", function(...)
{
	extend(Object(), "Discretize"
	)
})


#########################################################################/**
# @RdocMethod .volCube
#
# @aliasmethod .volCube
# @aliasmethod .klBound1D
# @aliasmethod .klBound2D
# @aliasmethod .bounds1D
# @aliasmethod .bounds2D
# @aliasmethod .ratio
# @aliasmethod .dir
# @aliasmethod .plotCubes
# @aliasmethod .addCube
# @aliasmethod .addCubeCol
# @aliasmethod .addPoints
# @aliasmethod .addIdx
# @aliasmethod .addText
#
# @title "Internal functions for discretizing a normal distribution"
#
# \usage{
#   \method{.volCube}{Discretize}(this, cube, ...)
#   \method{.klBound1D}{Discretize}(this, cube, mu, sigma2, ...)
#   \method{.klBound2D}{Discretize}(this, cube, mu, sigma, ...)
#   \method{.bounds1D}{Discretize}(this, cube, mu, sigma2, len=100, ...)
#   \method{.bounds2D}{Discretize}(this, cube, mu, sigma, len=100, ...)
#   \method{.ratio}{Discretize}(this, x, mu, sigma, ...)
#   \method{.dir}{Discretize}(this, cube, mu, sigma, ...)
#   \method{.plotCubes}{Discretize}(this, cubes, start, colors, ...)
#   \method{.addCube}{Discretize}(this, cube, col = "black", ...)
#   \method{.addCubeCol}{Discretize}(this, cube, color = NULL, ...)
#   \method{.addPoints}{Discretize}(this, cubes, ...)
#   \method{.addIdx}{Discretize}(this, cubes, ...)
#   \method{.addText}{Discretize}(this, cubes, text, ...)
# }
# \description{
#   \code{.volCube} returns volume/length of cube.
#   \code{.klBound1D} and \code{.klBound2D} returns upper bound on KL distance on a cube.
#   \code{.bounds1D} and \code{.bounds2D} returns min, mean and max density on a cube.
#   \code{.ratio} calc max/min function value.
#   \code{.dir} finds the optimal (approximate) direction to spilt a cube. Return the variable index to split.
#   \code{.plotCubes} plot the cubes (only bivariate distributions).
#   \code{.addCube} adds a 2D cube to the plot.
#   \code{.addCubeCol} adds a 2D cube with color to the plot.
#   \code{.addPoints} adds center points to the plot.
#   \code{.addIdx} adds cube index to the plot.
#   \code{.addText} adds text to the plot.
# }
#
# \arguments{
#  \item{cube}{ The cube under consideration which is a (2x2) matrix containing the bounds of the A variable in the first column and the bounds of X in the second (bivariate case) or vector of lenght 2 (univariate case).   }
#  \item{mu}{ The mean. }
#  \item{sigma2}{ The variance. }
#  \item{sigma}{ The std dev. }
#  \item{len}{ The number of samples of each coordinate.}
#  \item{cubes}{ The list of hypercubes. }
#  \item{start}{ An cube used to set the plot area.}
#  \item{colors}{ An integer vector of same length as the number of cubes used to give the cubes colors. The color is set by the integer value. }
#  \item{color}{ An integer setting the color. }
#  \item{text}{ Text to be added to each hypercube. Value 'center'
#    show the center point, 'index' show the index of the cube and if text i a vector of same length as the number of cube plot the text. }
#  \item{...}{Not used.}
# }
#
# @author
#
# \references{
# Based on
# \emph{Kozlov, A. & Koller, D. Nonuniform dynamic discretization in hybrid networks The Thirteenth Conference on Uncertainty in Artificial Intelligence (UAI-97), 1997, 314-325  }}
#
# @visibility "private"
#
#*/#########################################################################
setMethodS3(".volCube", "Discretize", function(this, cube, ...){
	return(prod(cube[2,]-cube[1,]))
})

setMethodS3(".klBound1D", "Discretize", function(this,cube,mu,sigma2, ...){
	if (cube[1,1]<= -Inf) cube[1,1]<-mu-2.5*sqrt(sigma2)
	if (cube[2,1]>= Inf) cube[2,1]<-mu+2.5*sqrt(sigma2)
	b<-this$.bounds1D(cube,mu,sigma2)
	return(((b$max-b$mean)/(b$max-b$min)*b$min*log(b$min/b$mean)+(b$mean-b$min)/(b$max-b$mean)*b$max*log(b$max/b$mean))*this$.volCube(cube))
})

setMethodS3(".klBound2D", "Discretize", function(this,cube,mu,sigma, ...){
	b<-this$.bounds2D(cube,mu,sigma)
	return(((b$max-b$mean)/(b$max-b$min)*b$min*log(b$min/b$mean)+(b$mean-b$min)/(b$max-b$mean)*b$max*log(b$max/b$mean))*this$.volCube(cube))
})

setMethodS3(".bounds1D", "Discretize", function(this, cube,mu,sigma2,len=100, ...){
	tmp<-matrix(NA,len,length(cube[1,]))
	for (i in 1:length(cube[1,])) {
		tmp[,i]<-seq(cube[1,i],cube[2,i],len=len)
	}
	g <- expand.grid(as.list(as.data.frame(tmp)))
	f<-dnorm(g[,1],mu,sqrt(sigma2))
	return(list(min=min(f),mean=mean(f),max=max(f)))
})

setMethodS3(".bounds2D", "Discretize", function(this, cube,mu,sigma,len=100, ...){
	tmp<-matrix(NA,len,length(cube[1,]))
	for (i in 1:length(cube[1,])) {
		tmp[,i]<-seq(cube[1,i],cube[2,i],len=len)
	}
	g <- expand.grid(as.list(as.data.frame(tmp)))
	f<-dmvnorm(g,mean=mu,sigma=sigma)
	return(list(min=min(f),mean=mean(f),max=max(f)))
})

setMethodS3(".ratio", "Discretize", function(this, x,mu,sigma, ...){
	f<-dmvnorm(x,mean=mu,sigma=sigma)
	return(max(f)/min(f))
})

setMethodS3(".dir", "Discretize", function(this, cube,mu,sigma, ...){
	l<-cube[2,]-cube[1,]    # length of variables in the cube
	center<-cube[1,]+l/2    # cube center
	idx<-0
	maxV<- -Inf
	for (i in 1:length(cube[1,])) {
		tmp<-matrix(center,100,length(center),byrow=TRUE)
		tmp[,i]<-seq(cube[1,i],cube[2,i],len=100)
		rat<-this$.ratio(tmp,mu,sigma)
		if (rat>maxV) {idx<-i; maxV=rat}
	}
	return(idx)
})

setMethodS3(".plotCubes", "Discretize", function(this, cubes, start, colors, ...) {
	plot(0,0,xlim=c(start[1,1],start[2,1]),ylim=c(start[1,2],start[2,2]),type="n",xlab="",ylab="")
	if (is.null(colors)) {
		for (i in 1:length(cubes)) {
			this$.addCube(cubes[[i]]$cubeB)
		}
	} else {
		for (i in 1:length(cubes)) {
			this$.addCubeCol(cubes[[i]]$cubeB,colors[i])
		}
		for (i in 1:length(cubes)) {
			this$.addCube(cubes[[i]]$cubeB)
		}
	}
})

setMethodS3(".addCube", "Discretize", function(this, cube,col="black", ...) {
	lines(c(cube[1,1],cube[1,1],cube[2,1],cube[2,1],cube[1,1]),c(cube[1,2],cube[2,2],cube[2,2],cube[1,2],cube[1,2]),col=col)
})

setMethodS3(".addCubeCol", "Discretize", function(this, cube,color=NULL, ...) {
	rect(cube[1,1], cube[1,2], cube[2,1], cube[2,2], col = color ,border="black")
})

setMethodS3(".addPoints", "Discretize", function(this, cubes, ...) {
	x<-y<-NULL
	for (i in 1:length(cubes)) {
		cube<-cubes[[i]]$center
		x<-c(x,cube[1])
		y<-c(y,cube[2])
	}
	points(x,y,pch=".")
})

setMethodS3(".addIdx", "Discretize", function(this, cubes, ...) {
	x<-y<-idx<-NULL
	for (i in 1:length(cubes)) {
		cube<-cubes[[i]]$center
		x<-c(x,cube[1])
		y<-c(y,cube[2])
		idx<-c(idx,i-1)
	}
	text(x,y,labels=paste(1:length(cubes)-1,sep=""))
})

setMethodS3(".addText", "Discretize", function(this, cubes,text, ...) {
	x<-y<-yield<-NULL
	for (i in 1:length(cubes)) {
		cube<-cubes[[i]]$center
		x<-c(x,cube[1])
		y<-c(y,cube[2])
	}
	text(x,y,labels=text)
})


#########################################################################/**
# @RdocMethod discretize1DUnifEqLth
#
# @title "Discretize a normal distribution such that intervals have equal length"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{mu}{ The mean. }
#  \item{sigma2}{ The variance. }
#  \item{n}{ Number of intervals. }
#  \item{asDF}{ Return result as a data frame. If false return matrix. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list of intervals (data frame if \code{asDF = TRUE}).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize1DUnifEqLth", "Discretize", function(this, mu, sigma2,
	n, asDF=TRUE, ...)
{
	lgd<-c(mu-2.5*sqrt(sigma2),mu+2.5*sqrt(sigma2))    # bounds used
	lgdInvX<-diff(lgd)/n  # length in each interval
	dat<-data.frame(center=NA,min=NA,max=NA,idxA=1:n-1)
	minX<-lgd[1]
	for (i in 1:n) {
		dat$min[i]<-minX
		dat$center[i]<-minX+lgdInvX/2
		dat$max[i]<-minX+lgdInvX
		minX<-minX+lgdInvX
	}
	dat$min[1]<- -Inf
	dat$max[nrow(dat)]<-Inf
	if (!asDF) return(as.matrix(dat))
	return(dat)
})


#########################################################################/**
# @RdocMethod discretize1DUnifEqProb
#
# @title "Discretize a normal distribution such that intervals have equal probability"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{mu}{ The mean. }
#  \item{sigma2}{ The variance. }
#  \item{n}{ Number of intervals. }
#  \item{asDF}{ Return result as a data frame. If false return matrix. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list of intervals (data frame if \code{asDF = TRUE}).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize1DUnifEqProb", "Discretize", function(this, mu, sigma2,
	n, asDF=TRUE, ...)
{
	pX<-1/n  # prob in each interval
	#xB<-c(mu-3*sqrt(sigma2),mu+3*sqrt(sigma2))    # bounds used
	q<-0; meanX<-NULL
	x<- -Inf
	for (i in 1:(n-1)) {
		x<-c(x,qnorm(pX*i,mu,sqrt(sigma2)))
		if (i==1) {
			meanX<-c(meanX,mu-sigma2*(dnorm(x[i+1],mu,sqrt(sigma2)))/pX)
		} else {
			meanX<-c(meanX,mu-sigma2*(dnorm(x[i+1],mu,sqrt(sigma2))-dnorm(x[i],mu,sqrt(sigma2)))/pX)
		}
	}
	x<-c(x,Inf)
	meanX<-c(meanX,mu-sqrt(sigma2)^2*(0-dnorm(x[i+1],mu,sqrt(sigma2)))/pX)

	elements<-vector("list", 2) # empty list of maxIte
	queue<-list(elements=elements)
	for (i in 1:(length(x)-1)) {
		cube<-matrix(c(x[i],x[i+1]),2,1)
		center<-meanX[i]
		element<-list(center=center,cube=cube)
		queue$elements[[i]]<-element
		queue$elements[[i]]<-element
	}
	for (i in 1:length(queue$elements)) {
		queue$elements[[i]]$idxA<- i-1
	}
	KL<-0
	for (i in 1:length(queue$elements)) {
		KL<-KL+this$.klBound1D(queue$elements[[i]]$cube,mu,sigma2)
	}
	cat("   KL-bound:",KL,"\n")
	if (!asDF) return(queue$elements)
	dF<-NULL
	for (i in 1:(length(x)-1)) {
		tmp1<-queue$elements[[i]]$cube
		tmp2<-queue$elements[[i]]$center
		tmp3<-queue$elements[[i]]$idxA
		dF<-rbind(dF,c(center=tmp2,min=tmp1[1,1],max=tmp1[2,1],idxA=tmp3))
	}
	rownames(dF)<-1:(length(x)-1)
	return(as.data.frame(dF))
})


#########################################################################/**
# @RdocMethod discretize1DVec
#
# @title "Discretize the real numbers according to a set of center points"
#
# \description{
#   @get "title". Create intervals with center points as given in the argument.
# }
#
# @synopsis
#
# \arguments{
#  \item{v}{ A vector of center points. }
#  \item{asDF}{ Return result as a data frame. If false return matrix. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list of intervals (data frame if \code{asDF = TRUE}).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize1DVec", "Discretize", function(this, v, asDF=TRUE, ...)
{
	v<-sort(v)
	dat<-data.frame(center=v,min=NA,max=NA,idxA=1:length(v)-1)
	for (i in 1:length(v)) {
		if (i==1) dat$min[i]<- -Inf
		else dat$min[i]<-dat$center[i]-(dat$center[i]-dat$center[i-1])/2
		if (i==length(v)) dat$max[i]<-Inf
		else dat$max[i]<-dat$center[i]+(dat$center[i+1]-dat$center[i])/2
	}
	if (!asDF) return(as.matrix(dat))
	return(dat)
})

#########################################################################/**
# @RdocMethod discretize2DNonunif
#
# @title "Discretize a bivariate normal distribution using a non-uniform discretization "
#
# \description{
#   Discretize a bivariate normal distribution into hypercubes (squares)
#   such that the approximation have a certain Kulback Libler (KL) distance.
# }
#
# @synopsis
#
# \arguments{
#  \item{mu}{ The mean (2-dim vector). }
#  \item{sigma}{ The covariance (2x2 matrix). }
#  \item{maxKL}{ Max KL distance. }
#  \item{maxIte}{ Max number of iterations. }
#  \item{modifyCenter}{ If no don't split the cubes around the mean center. If "split1" split the 4 cubes around the mean into 9 squares such that the mean is the center of a cube. If "split2" first add cubes such that the axis of the mean always in the center of the cubes. }
#  \item{split}{ Only used if modifyCenter = "split2" to set the size of the nine cubes around the mean. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list where each element describe the cube and contains: KL - an upper bound on the KL-distance, cube - the bounds, center - the center, idxM - the index, cubeB - the fixed bounds (used for plotting).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize2DNonunif", "Discretize", function(this, mu, sigma,
	maxKL=0.5, maxIte=500, modifyCenter="no", split=0.25, ...)
{
	xB<-c(mu[1]-2.5*sqrt(sigma[1,1]),mu[1]+2.5*sqrt(sigma[1,1]))
	yB<-c(mu[2]-2.5*sqrt(sigma[2,2]),mu[2]+2.5*sqrt(sigma[2,2]))
	cube0<-cube<-matrix(c(xB[1],xB[2],yB[1],yB[2]),2,2)
	if (modifyCenter!="split2") {
		KL<-this$.klBound2D(cube,mu,sigma)
		element<-list(KL=KL,cube=cube)   # the first element in the queue
		elements<-vector("list", 2) # empty list of 2 elements
		elements[[1]]<-element
		queue<-list(maxIdx=1, lastIdx=1, KL=KL,elements=elements)
	}
	if (modifyCenter=="split2") {
		# add nine cubes split around zero (numbered from topleft to bottomright)
		x<-c(mu[1]-split*sqrt(sigma[1,1]),mu[1]+split*sqrt(sigma[1,1]))
		y<-c(mu[2]-split*sqrt(sigma[2,2]),mu[2]+split*sqrt(sigma[2,2]))
		cube<-list()
		cube[[1]]<-matrix(c(xB[1],x[1],y[2],yB[2]),2,2)
		cube[[2]]<-matrix(c(x[1],x[2],y[2],yB[2]),2,2)
		cube[[3]]<-matrix(c(x[2],xB[2],y[2],yB[2]),2,2)
		cube[[4]]<-matrix(c(xB[1],x[1],y[1],y[2]),2,2)
		cube[[5]]<-matrix(c(x[1],x[2],y[1],y[2]),2,2) # the center cube
		cube[[6]]<-matrix(c(x[2],xB[2],y[1],y[2]),2,2)
		cube[[7]]<-matrix(c(xB[1],x[1],yB[1],y[1]),2,2)
		cube[[8]]<-matrix(c(x[1],x[2],yB[1],y[1]),2,2)
		cube[[9]]<-matrix(c(x[2],xB[2],yB[1],y[1]),2,2)
		elements<-list() # empty list
		KL<-maxI<-Max<-0
		for (i in 1:9) {
			cubeKL<-this$.klBound2D(cube[[i]],mu,sigma)
			if (cubeKL>Max) {
				Max<-cubeKL
				maxI<-i
			}
			KL<-KL+cubeKL
			element<-list(KL=cubeKL,cube=cube[[i]])
			elements[[i]]<-element
		}
		queue<-list(maxIdx=maxI, lastIdx=9, KL=KL,elements=elements)
	}
	ite<-1
	while (queue$KL>maxKL & ite<maxIte){
		maxIdx<-queue$maxIdx
		#cat("Total KL = ",queue$KL,"\n")
		KL<-queue$KL-queue$elements[[maxIdx]]$KL
		cube<-queue$elements[[maxIdx]]$cube
		#cat("Split cube:\n"); print(cube)
		splitIdx<-this$.dir(cube,mu,sigma)
		#cat("Split variable number ",splitIdx,"\n")
		split<-cube[1,splitIdx]+(cube[2,splitIdx]-cube[1,splitIdx])/2
		cube1<-cube2<-cube
		cube1[2,splitIdx]<-split
		cube2[1,splitIdx]<-split
		KL1<-this$.klBound2D(cube1,mu,sigma)
		KL2<-this$.klBound2D(cube2,mu,sigma)
		queue$KL<-KL+KL1+KL2
		element1<-list(KL=KL1,cube=cube1)
		element2<-list(KL=KL2,cube=cube2)
		queue$elements[[maxIdx]]<-element1
		queue$lastIdx<-queue$lastIdx+1
		queue$elements[[queue$lastIdx]]<-element2
		#cat("The two new elements:\n"); print(element1); print(element2);
		maxVal<- -Inf;
		for (i in 1:queue$lastIdx) {
			if (queue$elements[[i]]$KL>maxVal) {
				maxIdx<-i; maxVal<-queue$elements[[i]]$KL
			}
		}
		queue$maxIdx<-maxIdx; ite<-ite+1
	}
	if (modifyCenter=="split1") {
		# split the 4 cubes close to mu such that mu becomes the center of a cube
		idx<-NULL
		for (i in 1:queue$lastIdx) {    # first find cubes
			if (queue$elements[[i]]$cube[1,1]==mu[1] | queue$elements[[i]]$cube[2,1]==mu[1]) {
				if (queue$elements[[i]]$cube[1,2]==mu[2] | queue$elements[[i]]$cube[2,2]==mu[2]) {
					idx<-c(idx,i)
				}
			}
		}
		maxY=maxX=-Inf
		minY=minX=Inf
		for (i in idx) {
			maxX=max(maxX,queue$elements[[i]]$cube[2,1])
			maxY=max(maxY,queue$elements[[i]]$cube[2,2])
			minX=min(minX,queue$elements[[i]]$cube[1,1])
			minY=min(minY,queue$elements[[i]]$cube[1,2])
			queue$KL<-queue$KL-queue$elements[[i]]$KL
		}
		difX=(maxX-minX)/3
		difY=(maxY-minY)/3
		for (i in 0:2) {
			for (j in 0:2) {
				x=c(minX+i*difX,minX+(i+1)*difX)
				y=c(minY+j*difY,minY+(j+1)*difY)
				cube<-matrix(c(x[1],x[2],y[1],y[2]),2,2)
				KL<-this$.klBound2D(cube,mu,sigma)
				element<-list(KL=KL,cube=cube)
				if (!is.null(idx)) {    # if still some idx to change
					queue$elements[[idx[1]]]<-element
					if(length(idx)>1) {
						idx<-idx[2:length(idx)]
					} else {
						idx<-NULL
					}
				} else {
					queue$lastIdx<-queue$lastIdx+1
					queue$elements[[queue$lastIdx]]<-element
				}
				queue$KL<-queue$KL+KL
			}
		}
	}
	# find center
	for (i in 1:queue$lastIdx) {
		cube<-queue$elements[[i]]$cube
		x<-cube[1,1]+(cube[2,1]-cube[1,1])/2
		y<-cube[1,2]+(cube[2,2]-cube[1,2])/2
		queue$elements[[i]]$center<-c(x,y)
	}
	# set index
	for (i in 1:queue$lastIdx) {
		queue$elements[[i]]$idxM<- i-1
	}
	# remove borders (the one with borders saved in cubeB)
	cubes<-queue$elements
	m<-matrix(c(Inf,-Inf,Inf,-Inf),nrow=2,ncol=2)
	for (i in 1:length(cubes)) {    # min and max values, i.e. borders
		idx1<-cubes[[i]]$cube[1,]<m[1,]
		idx2<-cubes[[i]]$cube[2,]>m[2,]
		m[1,idx1]<-cubes[[i]]$cube[1,idx1]
		m[2,idx2]<-cubes[[i]]$cube[2,idx2]
	}
	for (i in 1:length(cubes)) {
		cubes[[i]]$cubeB<-cubes[[i]]$cube
		if (cubes[[i]]$cube[1,1]==m[1,1]) cubes[[i]]$cube[1,1]<- -Inf
		if (cubes[[i]]$cube[1,2]==m[1,2]) cubes[[i]]$cube[1,2]<- -Inf
		if (cubes[[i]]$cube[2,1]==m[2,1]) cubes[[i]]$cube[2,1]<- Inf
		if (cubes[[i]]$cube[2,2]==m[2,2]) cubes[[i]]$cube[2,2]<- Inf
	}
	cat("Total KL = ",queue$KL,"\n")
	return(cubes)
})


#########################################################################/**
# @RdocMethod discretize2DUnifEqInv
#
# @title "Discretize a bivariate normal distribution using a uniform discretization with intervals of equal length "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{mu}{ The mean (2-dim vector). }
#  \item{sigma}{ The covariance (2x2 matrix). }
#  \item{lgdX}{ Number for intervals of x coordinate. }
#  \item{lgdY}{ Number for intervals of y coordinate. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list where each element describe the cube and contains: KL - an upper bound on the KL-distance, cube - the bounds, center - the center, idxM - the index, cubeB - the fixed bounds (used for plotting).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize2DUnifEqInv", "Discretize", function(this, mu, sigma, lgdX, lgdY, ...){
	xB<-c(mu[1]-2.5*sqrt(sigma[1,1]),mu[1]+2.5*sqrt(sigma[1,1]))
	yB<-c(mu[2]-2.5*sqrt(sigma[2,2]),mu[2]+2.5*sqrt(sigma[2,2]))
	cube0<-cube<-matrix(c(xB[1],xB[2],yB[1],yB[2]),2,2)
	x<-seq(xB[1],xB[2],length=lgdX+1)
	y<-seq(yB[1],yB[2],length=lgdY+1)
	g <- expand.grid(x = x, y = y)
	z<-matrix(dmvnorm(g, mu, sigma),lgdX+1,lgdY+1)
	elements<-vector("list", 2) # empty list od two elements
	queue<-list(maxIdx=NA, lastIdx=1, KL=0,elements=elements)
	for (i in 1:(length(x)-1)) {
		for (j in 1:(length(y)-1)) {
			cube<-matrix(c(x[i],x[i+1],
						y[j],y[j+1]),2,2)
			KL<-this$.klBound2D(cube,mu,sigma)
			element<-list(KL=KL,cube=cube)
			queue$KL=queue$KL+KL
			queue$elements[[queue$lastIdx]]<-element
			queue$lastIdx<-queue$lastIdx+1
		}
	}
	queue$lastIdx<-queue$lastIdx-1
	# calc center point
	for (i in 1:queue$lastIdx) {
		cube<-queue$elements[[i]]$cube
		x<-cube[1,1]+(cube[2,1]-cube[1,1])/2
		y<-cube[1,2]+(cube[2,2]-cube[1,2])/2
		queue$elements[[i]]$center<-c(x,y)
	}
	# set index
	for (i in 1:queue$lastIdx) {
		queue$elements[[i]]$idxM <- i-1
	}
	# remove borders (the one with borders saved in cubeB)
	cubes<-queue$elements
	m<-matrix(c(Inf,-Inf,Inf,-Inf),nrow=2,ncol=2)
	for (i in 1:length(cubes)) {
		idx1<-cubes[[i]]$cube[1,]<m[1,]
		idx2<-cubes[[i]]$cube[2,]>m[2,]
		m[1,idx1]<-cubes[[i]]$cube[1,idx1]
		m[2,idx2]<-cubes[[i]]$cube[2,idx2]
	}
	for (i in 1:length(cubes)) {
		cubes[[i]]$cubeB<-cubes[[i]]$cube
		if (cubes[[i]]$cube[1,1]==m[1,1]) cubes[[i]]$cube[1,1]<- -Inf
		if (cubes[[i]]$cube[1,2]==m[1,2]) cubes[[i]]$cube[1,2]<- -Inf
		if (cubes[[i]]$cube[2,1]==m[2,1]) cubes[[i]]$cube[2,1]<- Inf
		if (cubes[[i]]$cube[2,2]==m[2,2]) cubes[[i]]$cube[2,2]<- Inf
	}
	cat("Total KL = ",queue$KL,"\n")
	return(cubes)
})



#########################################################################/**
# @RdocMethod discretize2DUnifEqProb
#
# @title "Discretize a bivariate normal distribution using a uniform discretization with intervals of equal probability "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{mu}{ The mean (2-dim vector). }
#  \item{sigma}{ The covariance (2x2 matrix). }
#  \item{lgdX}{ Number for intervals of x coordinate. }
#  \item{lgdY}{ Number for intervals of y coordinate. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list where each element describe the cube and contains: KL - an upper bound on the KL-distance, cube - the bounds, center - the center, idxM - the index, cubeB - the fixed bounds (used for plotting).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("discretize2DUnifEqProb", "Discretize", function(this,mu,sigma,lgdX,lgdY, ...){
	pX<-1/lgdX  # prob in each interval
	pY<-1/lgdY
	xB<-c(mu[1]-2.5*sqrt(sigma[1,1]),mu[1]+2.5*sqrt(sigma[1,1]))    # bounds used
	yB<-c(mu[2]-2.5*sqrt(sigma[2,2]),mu[2]+2.5*sqrt(sigma[2,2]))
	cube0<-cube<-matrix(c(xB[1],xB[2],yB[1],yB[2]),2,2)
	q<-0; meanX<-NULL
	x<-xB[1]
	for (i in 1:(lgdX-1)) {
		x<-c(x,qnorm(pX*i,mu[1],sqrt(sigma[1,1])))
		if (i==1) {
			meanX<-c(meanX,mu[1]-sqrt(sigma[1,1])^2*(dnorm(x[i+1],mu[1],sqrt(sigma[1,1])))/pX)
		} else {
			meanX<-c(meanX,mu[1]-sqrt(sigma[1,1])^2*(dnorm(x[i+1],mu[1],sqrt(sigma[1,1]))-dnorm(x[i],mu[1],sqrt(sigma[1,1])))/pX)
		}
	}
	x<-c(x,xB[2])
	i<-lgdX
	meanX<-c(meanX,mu[1]-sqrt(sigma[1,1])^2*(0-dnorm(x[i],mu[1],sqrt(sigma[1,1])))/pX)

	q<-0; meanY<-NULL
	y<-yB[1]
	for (i in 1:(lgdY-1)) {
		y<-c(y,qnorm(pY*i,mu[2],sqrt(sigma[2,2])))
		if (i==1) {
			meanY<-c(meanY,mu[2]-sqrt(sigma[2,2])^2*(dnorm(y[i+1],mu[2],sqrt(sigma[2,2])))/pY)
		} else {
			meanY<-c(meanY,mu[2]-sqrt(sigma[2,2])^2*(dnorm(y[i+1],mu[2],sqrt(sigma[2,2]))-dnorm(y[i],mu[2],sqrt(sigma[2,2])))/pY)
		}
	}
	y<-c(y,yB[2])
	i<-lgdY
	meanY<-c(meanY,mu[2]-sqrt(sigma[2,2])^2*(0-dnorm(y[i],mu[2],sqrt(sigma[2,2])))/pY)

	g <- expand.grid(x = x, y = y)
	m <- expand.grid(m1 = meanX, m2 = meanY)
	z<-matrix(dmvnorm(g, mu, sigma),lgdX+1,lgdY+1)
	elements<-vector("list", 2) # empty list of maxIte
	queue<-list(maxIdx=NA, lastIdx=1, KL=0,elements=elements)
	for (i in 1:(length(x)-1)) {
		for (j in 1:(length(y)-1)) {
			cube<-matrix(c(x[i],x[i+1],
						y[j],y[j+1]),2,2)
			KL<-this$.klBound2D(cube,mu,sigma)
			center<-c(meanX[i],meanY[j])
			element<-list(KL=KL,cube=cube,center=center)
			queue$KL=queue$KL+KL
			queue$elements[[queue$lastIdx]]<-element
			queue$lastIdx<-queue$lastIdx+1
		}
	}
	queue$lastIdx<-queue$lastIdx-1
	# calc center point
	for (i in 1:queue$lastIdx) {
		cube<-queue$elements[[i]]$cube
		x<-cube[1,1]+(cube[2,1]-cube[1,1])/2
		y<-cube[1,2]+(cube[2,2]-cube[1,2])/2
		queue$elements[[i]]$center<-c(x,y)
	}
	# set index
	for (i in 1:queue$lastIdx) {
		queue$elements[[i]]$idxM<-i-1
	}
	# remove borders (the one with borders saved in cubeB)
	cubes<-queue$elements
	m<-matrix(c(Inf,-Inf,Inf,-Inf),nrow=2,ncol=2)
	for (i in 1:length(cubes)) {
		idx1<-cubes[[i]]$cube[1,]<m[1,]
		idx2<-cubes[[i]]$cube[2,]>m[2,]
		m[1,idx1]<-cubes[[i]]$cube[1,idx1]
		m[2,idx2]<-cubes[[i]]$cube[2,idx2]
	}
	for (i in 1:length(cubes)) {
		cubes[[i]]$cubeB<-cubes[[i]]$cube
		if (cubes[[i]]$cube[1,1]==m[1,1]) cubes[[i]]$cube[1,1]<- -Inf
		if (cubes[[i]]$cube[1,2]==m[1,2]) cubes[[i]]$cube[1,2]<- -Inf
		if (cubes[[i]]$cube[2,1]==m[2,1]) cubes[[i]]$cube[2,1]<- Inf
		if (cubes[[i]]$cube[2,2]==m[2,2]) cubes[[i]]$cube[2,2]<- Inf
	}
	cat("Total KL = ",queue$KL,"\n")
	return(cubes)
})





#########################################################################/**
# @RdocMethod plotHypercubes
#
# @title "Plotting the discretization of a bivariate random normal variable  "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{cubes}{ The list of hypercubes. }
#  \item{text}{ Text to be added to each hypercube. Value 'center'
#    show the center point, 'index' show the index of the cube and if text i a vector of same length as the number of cube plot the text. }
#  \item{borders}{ Show the border of the hypercubes if true.}
#  \item{colors}{ A integer vector of same length as the number of cubes used to give the cubes colors. The color is set by the integer value. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A plot is produced.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("plotHypercubes", "Discretize", function(this, cubes,text="center",borders=FALSE, colors=NULL, ...){
	if (!is.null(colors)) {
		if (length(colors)!=length(cubes)) stop("Argument colors must have length equal to the number of cubes")
	}
	m<-matrix(c(Inf,-Inf,Inf,-Inf),nrow=2,ncol=2)
	for (i in 1:length(cubes)) {
		idx1<-cubes[[i]]$cubeB[1,]<m[1,]
		idx2<-cubes[[i]]$cubeB[2,]>m[2,]
		m[1,idx1]<-cubes[[i]]$cubeB[1,idx1]
		m[2,idx2]<-cubes[[i]]$cubeB[2,idx2]
	}
	cube0<-m
	this$.plotCubes(cubes,cube0,colors)
	if (!borders) this$.addCube(cube0,col="white")
	if (length(text)==length(cubes)) {
		this$.addText(cubes,text)
	} else {
		if (text=="index") this$.addText(cubes,1:length(cubes)-1)
		if (text=="center") this$.addPoints(cubes)
	}

	title(xlab=expression(m[1]),ylab=expression(m[2]))
	cat("   Plotted", length(cubes), "cubes.\n")
	invisible(NULL)
})




#########################################################################/**
# @RdocMethod splitCube2D
#
# @title "Split a cube further up "
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{cubes}{ The list of hypercubes. }
#  \item{mu}{ The mean (2-dim vector). }
#  \item{sigma}{ The covariance (2x2 matrix). }
#  \item{iM}{ Index of the cube that we want to split. }
#  \item{...}{Not used.}
# }
#
# \value{
#   A list where each element describe the cube and contains: KL - an upper bound on the KL-distance, cube - the bounds, center - the center, idxM - the index, cubeB - the fixed bounds (used for plotting).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @examples "../RdocFiles/Discretize.Rex"
#
#*/#########################################################################
setMethodS3("splitCube2D", "Discretize", function(this, cubes, mu, sigma, iM, ...) {
	# bounds the area
	xB<-c(mu[1]-2.5*sqrt(sigma[1,1]),mu[1]+2.5*sqrt(sigma[1,1]))
	yB<-c(mu[2]-2.5*sqrt(sigma[2,2]),mu[2]+2.5*sqrt(sigma[2,2]))

	cube<-cubes[[iM+1]]$cubeB
	#cat("Split cube:",maxIdx,"\n"); print(cube)
	#cat("KL=",queue$elements[[maxIdx]]$KL,"\n")
	splitIdx<-this$.dir(cube,mu,sigma)
	#cat("Split variable number ",splitIdx,"\n")
	split<-cube[1,splitIdx]+(cube[2,splitIdx]-cube[1,splitIdx])/2
	cube1<-cube2<-cube
	cube1[2,splitIdx]<-split
	cube2[1,splitIdx]<-split
	KL1<-this$.klBound2D(cube1,mu,sigma)
	KL2<-this$.klBound2D(cube2,mu,sigma)
	element1<-list(KL=KL1,cube=cube1,center=NA,idxM=iM,cubeB=cube1)
	element2<-list(KL=KL2,cube=cube2,center=NA,idxM=length(cubes),cubeB=cube2)
	cubes[[iM+1]]<-element1
	cubes[[length(cubes)+1]]<-element2
	# find center
	for (i in c(iM+1,length(cubes))) {
		cube<-cubes[[i]]$cube
		x<-cube[1,1]+(cube[2,1]-cube[1,1])/2
		y<-cube[1,2]+(cube[2,2]-cube[1,2])/2
		cubes[[i]]$center<-c(x,y)
	}
	# remove borders (the one with borders saved in cubeB)
	m<-matrix(c(Inf,-Inf,Inf,-Inf),nrow=2,ncol=2)
	for (i in 1:length(cubes)) {    # min and max values, i.e. borders of the cube
		idx1<-cubes[[i]]$cubeB[1,]<m[1,]
		idx2<-cubes[[i]]$cubeB[2,]>m[2,]
		m[1,idx1]<-cubes[[i]]$cubeB[1,idx1]
		m[2,idx2]<-cubes[[i]]$cubeB[2,idx2]
	}
	for (i in c(iM+1,length(cubes))) {
		if (cubes[[i]]$cube[1,1]==m[1,1]) cubes[[i]]$cube[1,1]<- -Inf
		if (cubes[[i]]$cube[1,2]==m[1,2]) cubes[[i]]$cube[1,2]<- -Inf
		if (cubes[[i]]$cube[2,1]==m[2,1]) cubes[[i]]$cube[2,1]<- Inf
		if (cubes[[i]]$cube[2,2]==m[2,2]) cubes[[i]]$cube[2,2]<- Inf
	}
	return(cubes)
})
