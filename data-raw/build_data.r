R --vanilla

m_visit <- data.table::fread("data-raw/m_visit.csv", header=TRUE)
m_count <- data.table::fread("data-raw/m_count.csv", header=TRUE)
metzger_v3_europe_orig <- raster::raster('data-raw/metzger_v3_eu_crop.tif')
metzger_v3_europe_values <- raster::values(metzger_v3_europe_orig)

##build empty raster
metzger_v3_europe <- raster::raster(ncol = 6480, nrow = 4440, xmn = -10.99999, xmx = 43.00001, ymn = 35, ymx = 72.00001)
raster::res(metzger_v3_europe) <- 0.008333334
raster::projection(metzger_v3_europe) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

metzger_v3_class <- data.table::fread('data-raw/GEnS_v3_classification.csv' )
names(metzger_v3_class) <- tolower(names(metzger_v3_class))

devtools::use_data(m_visit, m_count, metzger_v3_europe, metzger_v3_europe_values, metzger_v3_class, overwrite = TRUE)

library(rbms)

raster::values(metzger_v3_europe) <- metzger_v3_europe_values
raster::writeRaster(metzger_v3_europe, file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'), overwrite=TRUE)