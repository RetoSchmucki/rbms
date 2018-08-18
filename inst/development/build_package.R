R --vanilla
library(devtools)
library(roxygen2)

setwd("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/rbms")
document()
build()
install()