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

#mod<-HMDPDaily(maxLac=1,oOestus=Oestus(gestLth=30,pregTestLth=2,insemStart=2, insemFinish=4, dryPeriodLth=4), maxKL=0.9)
mod<-HMDPDaily(maxLac=1,csvPrM="prM.csv")
mod$genBinaryC(prefix="genCLac1_")

mod$maxLac<-2
mod$genBinaryC(prefix="genCLac2_")

mod$maxLac<-3
mod$genBinaryC(prefix="genCLac3_")

mod$maxLac<-4
mod$genBinaryC(prefix="genCLac4_")

mod$maxLac<-5
mod$genBinaryC(prefix="genCLac5_")

mod$maxLac<-6
mod$genBinaryC(prefix="genCLac6_")

mod$maxLac<-7
mod$genBinaryC(prefix="genCLac7_")

mod$maxLac<-8
mod$genBinaryC(prefix="genCLac8_")

mod$maxLac<-9
mod$genBinaryC(prefix="genCLac9_")

mod$maxLac<-10
mod$genBinaryC(prefix="genCLac10_")
rm(mod)
gc()

#mdp<-loadMDP(prefix="genC_")

### -------------------------------------------------------------------------------------------
