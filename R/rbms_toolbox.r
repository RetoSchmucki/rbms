##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   02.01.2018
##
##      rbms_toolbox: useful functions for the rbms package
##
##==========================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(" Welcome to rbms, version 0.1.0 \n This package is still in its beta version and active development.")
}


#' initiate_project
#' Build the initial folder structure for a typical research project
#' @param project_name name of the project name to be used for the parent folder
#' @param project_home path to where the project folder structure to be built, default is the current working directory
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export initiate_project

initiate_project <- function(project_name, project_home = NULL){

  if(is.null(project_home)){
      project_home <- getwd()
  }

	prj_path <- file.path(project_home,project_name)

	if (dir.exists(prj_path)) {
		setwd(prj_path)
	} else {
		cat(paste0('Creating project folder: ',prj_path,'\n'))
			dir.create(prj_path,recursive=TRUE)
			setwd(prj_path)}

	if (dir.exists(file.path(prj_path,'data'))) {
	     cat(paste('You are all set to work on',project_name,'project!','\n'))
	} else {
	   	setwd(prj_path)
	   	cat(paste0('Creating the folder structure and initiate the version control for ',project_name,'\n'))

        dir.create("R")
        dir.create("data")
        dir.create("documentation")
        dir.create("figures")
        dir.create("output")
        dir.create("analysis")
        dir.create("manuscript")
        dir.create("sql_script")
     }

	setwd(prj_path)
}



#' check_package
#' Internal function to verified the required package is installed
#' @param pkgName A string with the package name
#' @param message1 A string to inform about the dependency
#' @param message2 A string to inform what happen if not installed
#' @return If package is not installed, the function ask to install the package.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export check_package

check_package <- function(pkgName=NULL, message1='you need to install the package ',message2='This version requires '){
        if (!requireNamespace(pkgName)) {
            print(paste(message1, pkgName))
            x <- readline(paste("Do you want to install",pkgName,"? Y/N"))
            if (toupper(x) == 'Y') {
                    install.packages(pkgName)
            }
            if (toupper(x) == 'N') {
                print(paste(message2,pkgName))
            }
        }
    }


#' check_names
#' Verify for the required column names in the data
#' @param x a data.table object with column names
#' @param y vector with the required variable names
#' @return Verify if column names listed in \code{y} vector are found in the data set \code{x}, if not, a message identifies the
#' missing column name and stops.
#' @details This function is not case sensitive, but it does not accept different names or spelling.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @examples
#' DF <- data.frame(DAY=c(1:5),MONTH=rep(month.name[2],5),YEAR=rep(format(Sys.Date(),'%Y'),5))
#' check_names(DF,c('DAY','month','Years'))
#' @export check_names
#'

check_names <- function(x, y){
        dt_names <- y %in% names(x)

        if(sum(dt_names) != length(y)) {
            stop(paste('You need to have a variable named -', paste(y[!dt_names], collapse = ' & '), '- in table', deparse(substitute(x)), '\n'))
        }
    }


#' ts_date_seq
#' Generate a time-series of dates (per day) from the beginning of a starting year to the end of an ending years.
#' @param InitYear start year of the time-series, 4 numbers format (e.g 1987)
#' @param LastYear end year of the time-series, if not provided, current year is used instead
#' @return return a POSIXct vector with the format 'YYYY-MM-DD HH:MM:SS'
#' @keywords time series
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export ts_date_seq
#'

ts_date_seq <- function(InitYear=1970,LastYear=format(Sys.Date(),"%Y")) {

        init_date <- as.Date(paste((InitYear-1), "01-01", sep="-"))
        last_date <- as.Date(paste((as.numeric(LastYear) + 1), "12-31", sep = "-"))

        date_series <- as.POSIXct(seq(from=init_date, to=last_date, by="day"), format="%Y-%m-%d")
        date_series <- date_series[!format(date_series,'%Y') %in% c((InitYear-1),as.numeric(LastYear)+1)]

        return(date_series)

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
