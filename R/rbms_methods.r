#' plot.pheno_curve
#' Generic method to plot the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data pheno_curve object which is the outcome of the light_curve() function
#' @param year integer for the year to be displayed (e.g. 2015), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export plot.pheno_curve
#'

plot.pheno_curve <- function(data, year = NULL, weekday = 3, ...) {

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

    x <- data$pheno[, .(WEEK_SINCE, DATE)]
    y <- data$pheno[, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data$pheno[YEAR == year, .(WEEK_SINCE, DATE)]
        y <- data$pheno[YEAR == year, .(NM, ANCHOR)]
    }

    plot(x$DATE, y$NM, type = "l", ...)
}


#' points.pheno_curve
#' Generic method to add points on a plot of the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)
#' @param data pheno_curve object which is the outcome of the light_curve() function
#' @param year integer for the year to be displayed (e.g. 2015), default is NULL.
#' @param weekday weekday to be used for date in weekly count, where 1 refers to Monday, default is 3 (Wednesday)
#' @param ... additional parameters for base plot.
#' @return Returns a base plot with relative abundance (y) over time (x)
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export plot.pheno_curve
#'

points.pheno_curve <- function(data, year = NULL, weekday = 3, ...) {

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

    x <- data$pheno[, .(WEEK_SINCE, DATE)]
    y <- data$pheno[, .(NM, ANCHOR)]

    if (!is.null(year)) {
        x <- data$pheno[YEAR == year, .(WEEK_SINCE, DATE)]
        y <- data$pheno[YEAR == year, .(NM, ANCHOR)]
    }

    points(x$DATE, y$NM, ...)
}