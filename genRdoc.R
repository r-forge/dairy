# generate Rd files. MUST be run with working dir equal the R source folder!
library(R.oo)
library(mvtnorm)
library(lattice)
setwd("./pkg/trunk/R")
doc<-Rdoc()
doc$setNameFormat("class.method")
f <- list.files(".",pattern="\\.[R]",full.names=TRUE)
a <- lapply(f,source)
author <- "Lars Relund, \\url{http://www.research.relund.dk}"
auEJO <- "\\author{Erik Jørgensen, \\url{http://gbi.agrsci.dk/~ejo/}}"
doc$compile(check=F, source=FALSE) #check=TRUE, debug=TRUE
setwd("../../..")
#system("gen_and_install_dairy.bat", wait = FALSE, show.output.on.console = FALSE, invisible = FALSE)
#setwd("./dairy")
#q()

#doc$compile("Oestus.R")
