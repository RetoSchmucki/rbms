m_visit <- data.table::fread("data-raw/m_visit.csv", header=TRUE)
m_count <- data.table::fread("data-raw/m_count.csv", header=TRUE)
metzger_v3_europe <- raster::raster('data-raw/metzger_crop.tif')
metzger_v3_class <- data.table::fread('data-raw/GEnS_v3_classification.csv' )
names(metzger_v3_class) <- tolower(names(metzger_v3_class))

devtools::use_data(m_visit, m_count, metzger_v3_europe, metzger_v3_class, overwrite = TRUE)