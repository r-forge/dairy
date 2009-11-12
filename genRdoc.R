# generate Rd files. MUST be run with working dir equal the R source folder!
library(R.oo)
library(mvtnorm)
library(lattice)
doc<-Rdoc()
f <- list.files("pkg/trunk/R",pattern="\\.[R]",full.names=TRUE)
a <- lapply(f,source)
setwd("./pkg/trunk/R")
author <- "Lars Relund, \\url{http://www.research.relund.dk}"
doc$compile(check=F)#, source=FALSE, check=TRUE, debug=TRUE,)
#setwd("../..")
#system("gen_and_install_dairy.bat", wait = FALSE, show.output.on.console = FALSE, invisible = FALSE)
#setwd("./dairy")
#q()
