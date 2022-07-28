#' plot.pheno_curve
#' Generic method to plot the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data pheno_curve object which is the outcome of the light_curve() function
#' @param year integer for the year to be displayed (e.g. 2015), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param SiteID integer or string to ID the site to display, default the first site is displayed.
#' @param ymax maximum value for y axis
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export
#'

plot.pheno_curve <- function(data, year = NULL, weekday = 3, SiteID = NULL, ymax = NULL, ...) {

    if ("DAY" %in% names(data$pheno)) {
        data$pheno[, DATE := data.table::as.IDate(as.Date(paste(YEAR, MONTH, DAY, sep = "-")))]
    } else {
        f_y <- data$pheno[order(YEAR), YEAR][1]
        l_y <- data$pheno[rev(order(YEAR)), YEAR][1]
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

        setkey(data$pheno, YEAR, MONTH, WEEK)
        setkey(iso_week, YEAR, MONTH, WEEK)

        data$pheno <- merge(data$pheno, iso_week[WDAY == weekday, .(YEAR, MONTH, WEEK, DATE)], all.x = TRUE)

        setkey(data$pheno, WEEK_SINCE)
        f_date <- which(is.na(data$pheno[, DATE]))
        if (f_date[1] != 1) {
            data$pheno[f_date, DATE2 := data$pheno[f_date - 1, DATE + 7]]
        } else {
            data$pheno[f_date, DATE2 := data$pheno[f_date + 1, DATE - 7]]
        }
        data$pheno[is.na(DATE), DATE := DATE2][, DATE2 := NULL]
    }
    
    if(!is.null(SiteID)){
        SiteID <- as.character(SiteID)
        if(!SiteID %in% as.character(data$pheno[, unique(SITE_ID)])) stop("No flight curve available for the selected site, check site 'ts_flight_curve$pheno[, unique(SITE_ID)]'")
    }else{
        SiteID <- as.character(data$pheno[, unique(SITE_ID)][1])
    }
    
    x <- data$pheno[as.character(SITE_ID) == SiteID, .(WEEK_SINCE, DATE)]
    y <- data$pheno[as.character(SITE_ID) == SiteID, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data$pheno[as.character(SITE_ID) == SiteID & YEAR == year, .(WEEK_SINCE, DATE)]
        y <- data$pheno[as.character(SITE_ID) == SiteID & YEAR == year, .(NM, ANCHOR)]
    }

    plot(x$DATE, y$NM, type = "l", xlab = "Time", ylab = "Relative abundance (NM)", 
        ylim = c(0, ifelse(is.null(ymax), max(data$pheno$NM), ymax)), ...)

}

#' points.pheno_curve
#' Generic method to add points on a plot of the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data pheno_curve object which is the outcome of the light_curve() function
#' @param year integer for the year to be displayed (e.g. 2015), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param SiteID integer or string to ID the site to display, default the first site is displayed.
#' @param BaseYear integer to identify the base year to plot the points, default is the actual year, but can be use to plot additional year on existing plot
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @export
#'

points.pheno_curve <- function(data, year = NULL, weekday = 3, SiteID = NULL, BaseYear = NULL, ...) {

    if ("DAY" %in% names(data$pheno)) {
        data$pheno[, DATE := data.table::as.IDate(as.Date(paste(YEAR, MONTH, DAY, sep = "-"), format = "%Y-%m-%d"))]
    } else {
        f_y <- data$pheno[order(YEAR), YEAR][1]
        l_y <- data$pheno[rev(order(YEAR)), YEAR][1]
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

        setkey(data$pheno, YEAR, MONTH, WEEK)
        setkey(iso_week, YEAR, MONTH, WEEK)

        data$pheno <- merge(data$pheno, iso_week[WDAY == weekday, .(YEAR, MONTH, WEEK, DATE)], all.x = TRUE)

        setkey(data$pheno, WEEK_SINCE)
        f_date <- which(is.na(data$pheno[, DATE]))
        if (f_date[1] != 1) {
            data$pheno[f_date, DATE2 := data$pheno[f_date - 1, DATE + 7]]
        } else {
            data$pheno[f_date, DATE2 := data$pheno[f_date + 1, DATE - 7]]
        }
        data$pheno[is.na(DATE), DATE := DATE2][, DATE2 := NULL]
    }

    if(!is.null(SiteID)){
        SiteID <- as.character(SiteID)
        if(!SiteID %in% as.character(data$pheno[, unique(SITE_ID)])) stop("No flight curve available for the selected site, check site 'ts_flight_curve$pheno[, unique(SITE_ID)]'")
    }else{
        SiteID <- as.character(data$pheno[, unique(SITE_ID)][1])
    }

    x <- data$pheno[as.character(SITE_ID) == SiteID, .(WEEK_SINCE, DATE)]
    y <- data$pheno[as.character(SITE_ID) == SiteID, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data$pheno[as.character(SITE_ID) == SiteID & YEAR == year, .(WEEK_SINCE, DATE)]
        y <- data$pheno[as.character(SITE_ID) == SiteID & YEAR == year, .(NM, ANCHOR)]
    }

    if(!is.null(BaseYear)){
       x[, DATE := data.table::as.IDate(as.Date(paste(BaseYear, month(DATE), mday(DATE), sep = "-"), format = "%Y-%m-%d"))]
    }

    points(x$DATE, y$NM, ...)
}