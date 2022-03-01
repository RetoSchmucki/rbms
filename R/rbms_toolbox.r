##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   31.01.2021
##
##      rbms_toolbox: useful functions for the rbms package
##
##==========================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(" Welcome to rbms, version 1.1.2 \n This package has been tested, but is still in active development and feedbacks are welcome\n https://github.com/RetoSchmucki/rbms/issues")
}


#' initiate_project
#' Build the folder structure for a generic research project
#' @param project_name Name of the project to be used for the parent folder
#' @param project_home Path where the project folder structure to be built, default is the current working directory
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export initiate_project
#'

initiate_project <- function(project_name, project_home = NULL){

    if(!is.character(project_name)){
        stop("project_name must be defined as characters")
    }

    if(is.null(project_home)){
        project_home <- getwd()
    }
    
    prj_path <- file.path(project_home, project_name)

	if (dir.exists(prj_path)){
		setwd(prj_path)
        } else {
        cat(paste0('Building project folder: ', prj_path, '\n'))
            dir.create(prj_path, recursive = TRUE)
            setwd(prj_path)
    }

	if (dir.exists(file.path(prj_path, 'data'))) {
	     cat(paste('Ready to work on', project_name,'\n'))
	    } else {
            setwd(prj_path)
            cat(paste0('Building the folder structure for ', project_name, '\n'))
            dir.create("R")
            dir.create("data")
            dir.create("documentation")
            dir.create("output")
            dir.create("manuscript")

    }

	setwd(prj_path)
}


#' check_package
#' Internal function to check if a package is installed.
#' @param pkgName The package name to be verified.
#' @param message1 Inform the user about the package dependency.
#' @param message2 Inform the user what happen if the package is not installed.
#' @return If package is not found, the user is offered the option to install the missing package.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export check_package
#'

check_package <- function(pkgName = NULL, message1 = 'You need to install ', message2 = 'This version requires'){

    if(!is.character(pkgName)){
        stop("pkgName must be defined as characters")
    }

    if (!requireNamespace(pkgName)) {
        print(paste(message1, pkgName))
        x <- readline(paste("Do you want to install", pkgName, "? Y/N"))
        if (toupper(x) == 'Y') {
            install.packages(pkgName)
            }
        if (toupper(x) == 'N') {
                print(paste(message2, pkgName))
            }
    }
}


#' check_names
#' Verify for the required column names in the data object.
#' @param x Data object in which the variables names are verified.
#' @param y Variable names to be verified.
#' @return Verifyies if column names listed in \code{y} are found in the data set \code{x}, if not, a message identifies the
#' missing column name and stops.
#' @details This function is not case sensitive, but it does not accept different names or spelling.
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export check_names
#'

check_names <- function(x, y){

        if(!is.character(y)){
            stop("pkgName must be characters")
        }

        dt_names <- y %in% names(x)

        if(sum(dt_names) != length(y)) {
            stop(paste('You need a variable -', paste(y[!dt_names], collapse = ' & '), '- in table', deparse(substitute(x)), '\n'))
        }
    }


#' ts_date_seq
#' Generate a time-series with dates starting from January of the starting year to December of the ending years.
#' @param InitYear First year of the time-series, four digits format (e.g. 1987).
#' @param LastYear Last year of the time-series, if not provided, the current year is used.
#' @return Returns a POSIXct vector with the format 'YYYY-MM-DD HH:MM:SS'
#' @keywords time-series
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export ts_date_seq
#'

ts_date_seq <- function(InitYear = 1970, LastYear = format(Sys.Date(), "%Y")) {

    if(nchar(InitYear) != 4){stop("InitialYear must be provided in a four digits format (e.g. 1987)")}

    init_date <- as.Date(paste((InitYear - 1), "01-01", sep = "-"))
    last_date <- as.Date(paste((as.numeric(LastYear) + 1), "12-31", sep = "-"))

    date_series <- as.POSIXct(seq(from = init_date, to = last_date, by = "day"), format = "%Y-%m-%d")
    date_series <- date_series[!format(date_series, '%Y') %in% c((InitYear - 1), as.numeric(LastYear) + 1)]

    return(date_series)

}
