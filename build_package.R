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
metzger_v3_europe <- raster::raster('C:\\Users\\RETOSCHM\\Documents\\update2016\\metzger_crop.tif')
metzger_v3_class <- data.table::fread('W:\\PYWELL_SHARED\\Pywell Projects\\BRC\\BMS\\eBMS\\reto_workfiles\\ebms_database\\data\\Metzger_Climate\\global_dataset\\GEnSv3_25012013\\GEnS_v3_classification.csv' )
names(metzger_v3_class) <- tolower(names(metzger_v3_class))

devtools::use_data(m_visit,m_count)
devtools::use_data(metzger_v3_europe)
devtools::use_data(metzger_v3_class)
help(flight_curve)