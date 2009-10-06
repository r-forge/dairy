#########################################################################/**
# @RdocData meanYieldMat
#
# @title "Matrix of mean milk yields"
# \docType{data}
#
# \description{
#    A matrix of size (parities x dfc's) containing the default mean yield values used by class @see "CowYield".
# }
#
# \usage{data(meanYieldMat)}
#
# \details{
#   We assume that if for instance the number of rows is 3 then parity 3+ have the
#   same mean yield curve as parity 3.
# }
# \examples{
# data(meanYieldMat)
# meanYieldMat<-as.matrix(meanYieldMat)
# colnames(meanYieldMat)<-NULL
# plot(meanYieldMat[1,], type="l")
# }
#*/#########################################################################



#########################################################################/**
# @RdocData probPregTestMat
#
# @title "Matrix of positive pregnancy test probabilities"
# \docType{data}
#
# \description{
#   A matrix of size (parities x dfcs) containing the default probabilities used by class @see "Discretize".
# }
#
# \usage{data(probPregTestMat)}
#
# \details{
#   We assume that if for instance the number of rows is 3 then parity 3+ have the
#   same probabilities.
# }
#
# \examples{
# data(probPregTestMat)
# probPregTestMat<-as.matrix(probPregTestMat)
# colnames(probPregTestMat)<-NULL
# plot(probPregTestMat[1,], type="l")
# }
#*/#########################################################################
