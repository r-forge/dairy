# R script file
# Author: Lars Relund
# Description: Preliminary script for DairyCowYield package testing





### -------------------------------------------------------------------------------------------
## Playing with dairy package

library(dairy)
#mod<-HMDPDaily()
#mod<-HMDPDaily(maxLac=10,csvPrM="prM.csv")
#mod<-HMDPDaily(maxLac=3,
#    oOestus=Oestus(gestLth=100, insemStart=10, insemFinish=30, pregTestLth=5, dryPeriodLth=10),
#    csvPrM="prM.csv",  maxKL=0.8)
#mod$getFields()
#mod$genBinaryR(prefix="genR_")

mod<-HMDPDaily(maxLac=3,csvPrM="prM.csv")
mod$genBinaryC(prefix="genCLac3_")

#mod$maxLac<-2
#mod$genBinaryC(prefix="genCLac2_")
#
#mod$maxLac<-3
#mod$genBinaryC(prefix="genCLac3_")

#mod$maxLac<-4
#mod$genBinaryC(prefix="genCLac4_")
#
#mod$maxLac<-5
#mod$genBinaryC(prefix="genCLac5_")
#
#mod$maxLac<-10
#mod$genBinaryC(prefix="genCLac10_")
#rm(mod)
#gc()

#mdp<-loadMDP(prefix="genC_")

### -------------------------------------------------------------------------------------------
