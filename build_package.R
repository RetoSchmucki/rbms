R --vanilla
library(devtools)
library(roxygen2)


setwd("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/")
create('rbms')
setwd("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/rbms")
document()
build()
install()

library(rbms)

m_visit <- data.table::fread("data/m_visit.csv", header=TRUE)
m_count <- data.table::fread("data/m_count.csv", header=TRUE)
devtools::use_data(m_visit,m_count)

help(flight_curve)