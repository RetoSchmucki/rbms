#' plot.pheno_curve
#' Generic method to plot the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data_x pheno_curve object which is the outcome of the light_curve() function
#' @param year integer or vector for the year(s) to be displayed (e.g. 2015 or c(2000, 2001)), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param SiteID integer or string to ID the site to display, default the first site is displayed.
#' @param ymax maximum value for y axis
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export
#'

plot.pheno_curve <- function(data_x, year = NULL, weekday = 3, SiteID = NULL, ymax = NULL, ...) {

    if ("DAY" %in% names(data_x$pheno)) {
        data_x$pheno[, DATE := data.table::as.IDate(as.Date(paste(YEAR, MONTH, DAY, sep = "-")))]
    } else {
        f_y <- data_x$pheno[order(YEAR), YEAR][1]
        l_y <- data_x$pheno[rev(order(YEAR)), YEAR][1]
        date_seq <- ts_date_seq(f_y, l_y)
        w <- c(7, 1:6)

        iso_week <- data.table::data.table(
            DATE = data.table::as.IDate(date_seq),
            YEAR = data.table::year(date_seq),
            MONTH = data.table::month(date_seq),
            DAY = data.table::mday(date_seq),
            WEEK = data.table::isoweek(date_seq),
            WDAY = w[data.table::wday(date_seq)]
        )

        setkey(data_x$pheno, YEAR, MONTH, WEEK)
        setkey(iso_week, YEAR, MONTH, WEEK)

        data_x$pheno <- merge(data_x$pheno, iso_week[WDAY == weekday, .(YEAR, MONTH, WEEK, DATE)], all.x = TRUE)

        setkey(data_x$pheno, WEEK_SINCE)
        f_date <- which(is.na(data_x$pheno[, DATE]))
        if (f_date[1] != 1) {
            data_x$pheno[f_date, DATE2 := data_x$pheno[f_date - 1, DATE + 7]]
        } else {
            data_x$pheno[f_date, DATE2 := data_x$pheno[f_date + 1, DATE - 7]]
        }
        data_x$pheno[is.na(DATE), DATE := DATE2][, DATE2 := NULL]
    }
    
    if(!is.null(SiteID)){
        SiteID <- as.character(SiteID)
        site_list <- as.character(unlist(data_x$pheno[, unique(SITE_ID)]))
        if(!SiteID %in% site_list) stop("No flight curve available for the selected site, check site 'ts_flight_curve$pheno[, unique(SITE_ID)]'")
        if(!is.null(year)) {
            y_site_list <- as.character(unlist(data_x$pheno[YEAR %in% year, unique(SITE_ID)]))
            if(!SiteID %in% y_site_list) stop(paste("No flight curve available for the selected site and year, check site 'ts_flight_curve$pheno[YEAR ==",year,"year, unique(SITE_ID)]'"))
     }
     }else{
        SiteID <- as.character(data_x$pheno[, unique(SITE_ID)][1])
        if (!is.null(year)) {
            SiteID <- as.character(data_x$pheno[YEAR %in% year, unique(SITE_ID)][1])
        }else{
            SiteID <- as.character(data_x$pheno[, unique(SITE_ID)][1])
        }
    }
    
    x <- data_x$pheno[as.character(SITE_ID) == SiteID, .(WEEK_SINCE, DATE)]
    y <- data_x$pheno[as.character(SITE_ID) == SiteID, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data_x$pheno[as.character(SITE_ID) == SiteID & YEAR %in% year, .(WEEK_SINCE, DATE)]
        y <- data_x$pheno[as.character(SITE_ID) == SiteID & YEAR %in% year, .(NM, ANCHOR)]
    }

    plot(x$DATE, y$NM, type = "l", xlab = "Time", ylab = "Relative abundance (NM)", 
        ylim = c(0, ifelse(is.null(ymax), max(data_x$pheno$NM), ymax)), ...)

}

#' points.pheno_curve
#' Generic method to add points on a plot of the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data_x pheno_curve object which is the outcome of the light_curve() function
#' @param year integer or vector for the year(s) to be displayed (e.g. 2015 or c(2000, 2001)), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param SiteID integer or string to ID the site to display, default the first site is displayed.
#' @param BaseYear integer to identify the base year to plot the points, default is the actual year, but can be use to plot additional year on existing plot
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @importFrom graphics points
#' @export
#'

points.pheno_curve <- function(data_x, year = NULL, weekday = 3, SiteID = NULL, BaseYear = NULL, ...) {

    if ("DAY" %in% names(data_x$pheno)) {
        data_x$pheno[, DATE := data.table::as.IDate(as.Date(paste(YEAR, MONTH, DAY, sep = "-"), format = "%Y-%m-%d"))]
    } else {
        f_y <- data_x$pheno[order(YEAR), YEAR][1]
        l_y <- data_x$pheno[rev(order(YEAR)), YEAR][1]
        date_seq <- ts_date_seq(f_y, l_y)
        w <- c(7, 1:6)

        iso_week <- data.table::data.table(
            DATE = data.table::as.IDate(date_seq),
            YEAR = data.table::year(date_seq),
            MONTH = data.table::month(date_seq),
            DAY = data.table::mday(date_seq),
            WEEK = data.table::isoweek(date_seq),
            WDAY = w[data.table::wday(date_seq)]
        )

        setkey(data_x$pheno, YEAR, MONTH, WEEK)
        setkey(iso_week, YEAR, MONTH, WEEK)

        data_x$pheno <- merge(data_x$pheno, iso_week[WDAY == weekday, .(YEAR, MONTH, WEEK, DATE)], all.x = TRUE)

        setkey(data_x$pheno, WEEK_SINCE)
        f_date <- which(is.na(data_x$pheno[, DATE]))
        if (f_date[1] != 1) {
            data_x$pheno[f_date, DATE2 := data_x$pheno[f_date - 1, DATE + 7]]
        } else {
            data_x$pheno[f_date, DATE2 := data_x$pheno[f_date + 1, DATE - 7]]
        }
        data_x$pheno[is.na(DATE), DATE := DATE2][, DATE2 := NULL]
    }

    if(!is.null(SiteID)){
        SiteID <- as.character(SiteID)
        site_list <- as.character(unlist(data_x$pheno[, unique(SITE_ID)]))
        if(!SiteID %in% site_list) stop("No flight curve available for the selected site, check site 'ts_flight_curve$pheno[, unique(SITE_ID)]'")
        if(!is.null(year)) {
            y_site_list <- as.character(unlist(data_x$pheno[YEAR %in% year, unique(SITE_ID)]))
            if(!SiteID %in% y_site_list) stop(paste("No flight curve available for the selected site and year, check site 'ts_flight_curve$pheno[YEAR ==",year,"year, unique(SITE_ID)]'"))
     }
     }else{
        SiteID <- as.character(data_x$pheno[, unique(SITE_ID)][1])
        if (!is.null(year)) {
            SiteID <- as.character(data_x$pheno[YEAR %in% year, unique(SITE_ID)][1])
        }else{
            SiteID <- as.character(data_x$pheno[, unique(SITE_ID)][1])
        }
    }

    x <- data_x$pheno[as.character(SITE_ID) == SiteID, .(WEEK_SINCE, DATE)]
    y <- data_x$pheno[as.character(SITE_ID) == SiteID, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data_x$pheno[as.character(SITE_ID) == SiteID & YEAR %in% year, .(WEEK_SINCE, DATE)]
        y <- data_x$pheno[as.character(SITE_ID) == SiteID & YEAR %in% year, .(NM, ANCHOR)]
    }

    if(!is.null(BaseYear)){
       x[, DATE := data.table::as.IDate(as.Date(paste(BaseYear, month(DATE), mday(DATE), sep = "-"), format = "%Y-%m-%d"))]
    }

    points(x$DATE, y$NM, ...)
}