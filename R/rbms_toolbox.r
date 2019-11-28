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
  packageStartupMessage(" Welcome to rbms, version 1.0.0 \n While this package has been tested by serveral users,\n it is still in active development and feedbacks are welcome \n https://github.com/RetoSchmucki/rbms/issues")
}


#' initiate_project
#' Build the initial folder structure for a typical research project
#' @param project_name name of the project name to be used for the parent folder
#' @param project_home path to where the project folder structure to be built, default is the current working directory
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export initiate_project
#'

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
#'

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
