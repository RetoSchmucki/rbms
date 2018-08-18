

#' get_raster_value
#' Extract raster value for a set of geographic points
#' @param x data.frame with longitude and latitude in column 1 and 2 respectively
#' @param yPath complete path to raster layer from which value should be extracted
#' @param Classification a data.frame with classification for the raster values
#' @param xCrs EPSG number of the cordinate reference system (CRS) of the x coordinates, see \url{http://spatialreference.org}
#' @param BufferDist extent of the buffer to be added around the bounding box to insure, default set to 50 pixels of the raster
#' @param OutDf logical if true return a data.frame, if false, a sf object is returned
#' @param PlotRaster logical informing if points, their bounding box and the raster layer should be plotted
#' @return return a data.frame or a sf object with value extracted from the y layer, see \code{out_df}.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export get_raster_value
#' @examples
#' x <- data.frame(longitude = c(4, 4.1, 4.5), latitude = c(50, 50.45, 50.5), id = c('a','b','c'))
#' (x_value <- get_raster_value(x))
#' (x_value <- get_raster_value(x, PlotRaster = TRUE))
#'

get_raster_value <- function(x, yPath = 'metzger_v3_europe' , Classification = NULL, xCrs = 4326, BufferDist = NULL, OutDf = TRUE, PlotRaster = FALSE){

        if(yPath == 'metzger_v3_europe'){
            if(!file.exists(file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'))){
                message(paste0('On first use, we need to build the metzger_v3_europe raster, this is stored in \n',
                                    file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif')))
                raster::values(metzger_v3_europe) <- metzger_v3_europe_values
                raster::writeRaster(metzger_v3_europe, file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'), overwrite=TRUE)
            }
            y_r <- raster::raster(file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'))
        } else {
            y_r <- raster::raster(yPath)
        }

        x_sf <- sf::st_as_sf(x, coords = names(x)[1:2], crs = xCrs, agr = 'constant')
        x_sf <- sf::st_transform(x_sf, as.character(raster::crs(y_r)))
        if(is.null(BufferDist)){
            r <- raster::res(y_r)[1]*50
        }
        focal_box <- as(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(x_sf)), raster::res(y_r)[1]*50), "Spatial")
        layer_cr <- raster::crop(y_r, raster::extent(focal_box))

        if(isTRUE(PlotRaster)){
            dev.new()
            raster::plot(y_r)
            plot(sf::st_as_sf(focal_box)$geometry, col = 'orange', lwd = 1, lty = 2, add = TRUE)
            plot(x_sf$geometry, col = 'blue', pch = 19, add = TRUE)
        }

        x_sf$value <- raster::extract(layer_cr,as(x_sf,"Spatial"))
        x_sf <- sf::st_transform(x_sf, xCrs)

        if(isTRUE(OutDf)){
            x_sf$longitude <- sf::st_coordinates(x_sf)[,1]
            x_sf$latitude <- sf::st_coordinates(x_sf)[,2]
            sf::st_geometry(x_sf) <- NULL
        }

    return(x_sf)
}

#' get_bioclim
#' Assign a bioclimatic region to the value extractacted from raster layer for a set of geographic points
#' @param x data.frame, or sf object, with longitude and latitude in column 1 and 2 respectively and bioclim value extracted with get_raster_value
#' @param y data.frame with bioclimatic region classification key refering to the raster layer used by get_raster_value
#' @param byY character name of the variable name corresponding to the raster layer value in y
#' @param xCrs EPSG number of the cordinate reference system (CRS) of the x coordinates, see \url{http://spatialreference.org}
#' @param OutDf logical if true return a data.frame, if false, a sf object is returned
#' @return return a data.frame or a sf object with bioclimatic region extracted from the y layer.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export get_bioclim
#' @examples
#' x <- data.frame(longitude = c(4, 4.1, 4.5), latitude = c(50, 50.45, 50.5), id = c('a','b','c'))
#' x_value <- get_raster_value(x, OutDf = FALSE)
#' get_bioclim(x_value)
#' get_bioclim(x_value, OutDf = FALSE)
#'

get_bioclim <- function(x, y = 'metzger_v3_class', byY = 'gens_seq', xCrs = 4326, OutDf = TRUE){

        if(class(x)[1]== 'sf'){
            x$longitude <- sf::st_coordinates(x)[,1]
            x$latitude <- sf::st_coordinates(x)[,2]
            sf::st_geometry(x) <- NULL
        }

        if(y == 'metzger_v3_class'){
            x_bioclim <- merge(x, metzger_v3_class[,c('gens_seq', 'genzname', 'genz', 'gens')], by.x = 'value', by.y = byY, all.x = TRUE)
        } else {
            x_bioclim <- merge(x, y, by.x = 'value', by.y = byY, all.x = TRUE)
        }

        if(!isTRUE(OutDf)){
          x_bioclim <- sf::st_as_sf(x_bioclim, coords = c('longitude','latitude'), crs = xCrs, agr = 'constant')
        }

    return(x_bioclim)

}

#' get_raster_value
#' Extract raster value for a set of geographic points
#' @param x data.frame with longitude and latitude in column 1 and 2 respectively
#' @param yPath complete path to raster layer from which value should be extracted
#' @param Classification a data.frame with classification for the raster values
#' @param xCrs EPSG number of the cordinate reference system (CRS) of the x coordinates, see \url{http://spatialreference.org}
#' @param BufferDist extent of the buffer to be added around the bounding box to insure, default set to 50 pixels of the raster
#' @param OutDf logical if true return a data.frame, if false, a sf object is returned
#' @param PlotRaster logical informing if points, their bounding box and the raster layer should be plotted
#' @return return a data.frame or a sf object with value extracted from the y layer, see \code{out_df}.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export get_raster_value
#' @examples
#' x <- data.frame(longitude = c(4, 4.1, 4.5), latitude = c(50, 50.45, 50.5), id = c('a','b','c'))
#' (x_value <- get_raster_value(x))
#' (x_value <- get_raster_value(x, PlotRaster = TRUE))
#'

get_raster_value <- function(x, yPath = 'metzger_v3_europe' , Classification = NULL, xCrs = 4326, BufferDist = NULL, OutDf = TRUE, PlotRaster = FALSE){

        if(yPath == 'metzger_v3_europe'){
            if(!file.exists(file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'))){
                message(paste0('On first use, we need to build the metzger_v3_europe raster, this is stored in \n',
                                    file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif')))
                raster::values(metzger_v3_europe) <- metzger_v3_europe_values
                raster::writeRaster(metzger_v3_europe, file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'), overwrite=TRUE)
            }
            y_r <- raster::raster(file.path(system.file(package = 'rbms'), 'raster_data/metzger_v3_europe.tif'))
        } else {
            y_r <- raster::raster(yPath)
        }

        x_sf <- sf::st_as_sf(x, coords = names(x)[1:2], crs = xCrs, agr = 'constant')
        x_sf <- sf::st_transform(x_sf, as.character(raster::crs(y_r)))
        if(is.null(BufferDist)){
            r <- raster::res(y_r)[1]*50
        }
        focal_box <- as(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(x_sf)), raster::res(y_r)[1]*50), "Spatial")
        layer_cr <- raster::crop(y_r, raster::extent(focal_box))

        if(isTRUE(PlotRaster)){
            dev.new()
            raster::plot(y_r)
            plot(sf::st_as_sf(focal_box)$geometry, col = 'orange', lwd = 1, lty = 2, add = TRUE)
            plot(x_sf$geometry, col = 'blue', pch = 19, add = TRUE)
        }

        x_sf$value <- raster::extract(layer_cr,as(x_sf,"Spatial"))
        x_sf <- sf::st_transform(x_sf, xCrs)

        if(isTRUE(OutDf)){
            x_sf$longitude <- sf::st_coordinates(x_sf)[,1]
            x_sf$latitude <- sf::st_coordinates(x_sf)[,2]
            sf::st_geometry(x_sf) <- NULL
        }

    return(x_sf)
}

#' get_bioclim
#' Assign a bioclimatic region to the value extractacted from raster layer for a set of geographic points
#' @param x data.frame, or sf object, with longitude and latitude in column 1 and 2 respectively and bioclim value extracted with get_raster_value
#' @param y data.frame with bioclimatic region classification key refering to the raster layer used by get_raster_value
#' @param byY character name of the variable name corresponding to the raster layer value in y
#' @param xCrs EPSG number of the cordinate reference system (CRS) of the x coordinates, see \url{http://spatialreference.org}
#' @param OutDf logical if true return a data.frame, if false, a sf object is returned
#' @return return a data.frame or a sf object with bioclimatic region extracted from the y layer.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export get_bioclim
#' @examples
#' x <- data.frame(longitude = c(4, 4.1, 4.5), latitude = c(50, 50.45, 50.5), id = c('a','b','c'))
#' x_value <- get_raster_value(x, OutDf = FALSE)
#' get_bioclim(x_value)
#' get_bioclim(x_value, OutDf = FALSE)
#'

get_bioclim <- function(x, y = 'metzger_v3_class', byY = 'gens_seq', xCrs = 4326, OutDf = TRUE){

        if(class(x)[1]== 'sf'){
            x$longitude <- sf::st_coordinates(x)[,1]
            x$latitude <- sf::st_coordinates(x)[,2]
            sf::st_geometry(x) <- NULL
        }

        if(y == 'metzger_v3_class'){
            x_bioclim <- merge(x, metzger_v3_class[,c('gens_seq', 'genzname', 'genz', 'gens')], by.x = 'value', by.y = byY, all.x = TRUE)
        } else {
            x_bioclim <- merge(x, y, by.x = 'value', by.y = byY, all.x = TRUE)
        }

        if(!isTRUE(OutDf)){
          x_bioclim <- sf::st_as_sf(x_bioclim, coords = c('longitude','latitude'), crs = xCrs, agr = 'constant')
        }

    return(x_bioclim)

}
