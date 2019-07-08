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
##      organize_data: functions for organizing BMS datad for modelling flight curve and compute indices
##
##==========================================

#' ts_dwmy_table
#' Generate a time-series of dates with day, week, month and year (dwmy) from one initial to an end years
#' @param InitYear start year of the time-series, 4 numbers format (e.g 1987)
#' @param LastYear end year of the time-series, if not provided, current year is used instead
#' @param WeekDay1 to start the week on monday, use 'monday', otherwise the week start on sunday
#' @return a data.table object with the date, the day since the first date, the week since the first week, the year, the month,
#'          the day in the month, the ISO week number and the day in the week.
#' @keywords time series
#' @seealso \link{ts_date_seq}, \link[data.table]{IDate}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export ts_dwmy_table
#'

ts_dwmy_table = function(InitYear=1970, LastYear=format(Sys.Date(), "%Y"), WeekDay1='monday') {

        check_package('data.table')
        date_seq <- ts_date_seq(InitYear,LastYear)

        if(WeekDay1=='monday'){
            w <- c(7,1:6)
        }else{
            w <- c(1:7)
        }

        w_s <- c(1, diff(data.table::isoweek(date_seq)))
        w_s[w_s < 0] <- 1
        w_s <- cumsum(w_s)

        dt_iso_dwmy <- data.table::data.table(
                        DATE = data.table::as.IDate(date_seq),
                        DAY_SINCE = seq_along(date_seq),
                        YEAR = data.table::year(date_seq),
                        MONTH = data.table::month(date_seq),
                        DAY = data.table::mday(date_seq),
                        WEEK = data.table::isoweek(date_seq),
                        WEEK_SINCE = w_s,
                        WEEK_DAY = w[data.table::wday(date_seq)])

        return(dt_iso_dwmy)

    }


#' set_anchor
#' Add Anchors of "zeros" at determined distance on each side of the monitoring season with specific weight (length),
#' this function is used by ts_monit_season()
#' @param FirstObs integer defining the start of the monitoring season - correspond to the day since
#' @param LastObs integer defining the end of the monitoring season - correspond to the day since
#' @param AnchorLength integer defining the number of days used as Anchor each side of the monitoring season
#' @param AnchorLag integer defining the number of days between the Anchor and the monitoring season
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export set_anchor
#'

set_anchor <- function(FirstObs,LastObs,AnchorLength=7,AnchorLag=7){

        x <- FirstObs$V1-(AnchorLength+AnchorLag)
        y <- FirstObs$V1-(AnchorLag+1)
        before_anchor <- c(apply(as.matrix(cbind(x,y)),1,function(w) w[1]:w[2]))

        x <- LastObs$V1+(AnchorLength+AnchorLag)
        y <- LastObs$V1+(AnchorLag+1)
        after_anchor <- c(apply(as.matrix(cbind(x,y)),1,function(w) w[2]:w[1]))

        anchor_day <- c(before_anchor,after_anchor)

        return(anchor_day)

    }


#' ts_monit_season
#' Build a time-series of dates with specific detail about the monitoring season
#' @param d_series time series of dates returned by \link{ts_dwmy_table}
#' @param StartMonth integer for the first month of the monitoring season, default=4
#' @param EndMonth integer for the last month of the monitoring season, default=9
#' @param StartDay integer for the day of the month when the monitoring season start, default=1
#' @param EndDay integer for the day of the month when the monitoring season end, default=last day of the EndMonth
#' @param CompltSeason logical if only the date for a complete season should be returned, default=TRUE
#' @param Anchor logical if Anchor should be used at the beginning and end of the monitoring season, default=TRUE
#' @param AnchorLength integer for the number of day used as Anchor, default=7
#' @param AnchorLag integer for the number of days before and after where the Anchor should start, default=TRUE
#' @return A data.table with the entire time-series of date (\code{DATE}, \code{DAY_SINCE}, \code{YEAR}, \code{MONTH},
#'         \code{DAY}, \code{WEEK} (ISO), the \code{WEEK_DAY}, and details about the Monitoring Year (\code{M_YEAR}) which
#'         is the monitoring year to which that date refers to, using the year of the starting month of the monitoring season,
#'         the monitoring season (\code{M_SEASON}), if the Monitoring season is complete (\code{COMPLT_SEASON}), the
#'         location of the anchors (\code{ANCHOR}) and the count for every day used as anchor (0).
#' @details The monitoring season can start in year y and end in year y+1, so if the monitoring season goes over two
#'         years (i.e. winter, December-January) the SEASON YEAR is shifted to set the monitoring season within a
#'         continuous year named according to the year at the start. Here the monitoring season is set in the middle
#'         or year' to give some room to set ANCHORs  before and after the monitoring season.
#' @seealso \link{ts_date_seq}, \link{ts_dwmy_table}, \link{set_anchor}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export ts_monit_season
#'

ts_monit_season = function(d_series, StartMonth=4, EndMonth=9, StartDay=1, EndDay=NULL, CompltSeason=TRUE,
                            Anchor=TRUE, AnchorLength=7, AnchorLag=7){

        check_package('data.table')
        d_series <- data.table::copy(d_series)

        names(d_series) <- toupper(names(d_series))
        check_names(d_series,c("DATE"))

        ## NO leap year
        if(is.null(EndDay)) {
            EndDay <- max(data.table::mday(seq(from = as.Date(paste0('2017-', EndMonth, '-01')), by = 'day', length = 31)))
        }

        if (StartMonth < EndMonth){
            s_month <- c(StartMonth:EndMonth)
            month_days_out <- c(paste(StartMonth, c(0:(StartDay-1))), paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(M_YEAR=as.factor(data.table::year(d_series$DATE)),
                                               M_SEASON=ifelse(((data.table::month(d_series$DATE) %in% (s_month)) &
                                               (!(paste(data.table::month(d_series$DATE), data.table::mday(d_series$DATE))
                                               %in% month_days_out))), 1L, 0L))
        }

        if (StartMonth > EndMonth){
            s_month <- c(StartMonth:12, 1:EndMonth)
            month_days_out <- c(paste(StartMonth, c(0:(StartDay - 1))), paste(EndMonth, c((EndDay + 1):32)))
            y_series <- data.table::data.table(M_YEAR = as.factor(ifelse(data.table::month(d_series$DATE) >= StartMonth - floor((12 - length(s_month)) / 2), data.table::year(d_series$DATE), (data.table::year(d_series$DATE) - 1))),
                                               M_SEASON=ifelse(((data.table::month(d_series$DATE)%in%(s_month)) & (!(paste(data.table::month(d_series$DATE), data.table::mday(d_series$DATE)) %in% month_days_out))), 1L, 0L))
        }

        d_series <- d_series[, c("M_YEAR", "M_SEASON") := y_series[, .(M_YEAR, (as.numeric(M_YEAR) * M_SEASON))]]

        if(isTRUE(CompltSeason)){
            d_series[, START_END := d_series[,ifelse(data.table::month(DATE) == StartMonth & data.table::mday(DATE) == StartDay, 1L, 0L)] + d_series[, ifelse(data.table::month(DATE) == EndMonth & data.table::mday(DATE) == EndDay, 1L, 0L)]]
            d_series[, M_SEASON := ifelse(M_YEAR %in% (d_series[, sum(START_END), by = M_YEAR][V1 == 2, M_YEAR]), 1L, 0L) * M_SEASON][,START_END:=NULL]
            d_series[M_SEASON > 0L, M_SEASON := M_SEASON - (min(d_series[M_SEASON > 0L, M_SEASON]) - 1)]
            d_series[, COMPLT_SEASON := ifelse(M_YEAR %in% d_series[M_SEASON > 0L, unique(M_YEAR)], 1L, 0L)]
        }

        d_series[, ANCHOR := 0L]

        if(isTRUE(Anchor)){
            first_obs <- d_series[M_SEASON > 0L, min(DAY_SINCE), by = .(M_YEAR)]
            last_obs <- d_series[M_SEASON > 0L, max(DAY_SINCE), by = .(M_YEAR)]
            anchor_day <- set_anchor(FirstObs = first_obs, LastObs = last_obs, AnchorLength = AnchorLength, AnchorLag = AnchorLag)
            d_series <- d_series[DAY_SINCE %in% anchor_day, ANCHOR := 1L] [DAY_SINCE %in% anchor_day, COUNT := 0L]
        }

        return(d_series)
    }


#' df_visit_season
#' Link each recorded visit to a corresponding monitoring season, this function is used in \link{df_visit_season}
#' @param m_visit data.table or data.frame with Date and Site ID for each monitoring visit.
#' @param ts_season data.table returned by \link{ts_monit_season} with the detail time-series of the monitoring season.
#' @param DateFormat format used for the date in the visit data, default="\%Y-\%m-\%d".
#' @return A data.table where each visit is attributed a monitoring year, the \code{M_YEAR} as this can differ the date (see detail)
#' @details The value of \code{M_YEAR} should be used as Monitoring year as this can differ form the \code{YEAR} if the
#'         monitoring season covers two calendar years (e.g. November to June).
#' @seealso \link{ts_monit_season}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export df_visit_season
#'

df_visit_season <- function(m_visit, ts_season, DateFormat="%Y-%m-%d"){

            check_package('data.table')
            m_visit <- data.table::data.table(m_visit)
            m_visit <- data.table::copy(m_visit[, .(DATE, SITE_ID)])
            m_visit[, DATE := data.table::as.IDate(as.Date(m_visit$DATE, format = DateFormat))]

            season_year <- ts_season[, .(DATE, M_YEAR)]

            data.table::setkey(m_visit, DATE)
            data.table::setkey(season_year, DATE)

            m_visit <- merge(m_visit, season_year, all.x=FALSE)

        return(m_visit)

    }


#' ts_monit_site
#' Augment the time series in m_season with all sites and visits with "zeros", leaving all non visited day with and <NA>
#' @param m_visit data.table or data.frame with Date and Site ID for each monitoring visit.
#' @param ts_season data.table returned by \link{ts_monit_season} with the detail time-series of the monitoring season.
#' @param DateFormat format used for the date in the visit data, default="\%Y-\%m-\%d".
#' @return A data.table with a complete time-series where absences have been implemented and that is ready to receive counts
#' @seealso \link{df_visit_season}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export ts_monit_site
#'

ts_monit_site = function(m_visit, ts_season, DateFormat="%Y-%m-%d") {

            check_package('data.table')
            names(ts_season) <- toupper(names(ts_season))
            check_names(ts_season,c("DATE", "M_YEAR", "M_SEASON"))

            names(m_visit) <- toupper(names(m_visit))
            check_names(m_visit, c("DATE", "SITE_ID"))

            m_visit <- df_visit_season(m_visit, ts_season, DateFormat = DateFormat)

            data.table::setkey(ts_season, DATE)
            data.table::setkey(m_visit, DATE)

            r_year <- ts_season[, range(data.table::year(DATE))]
            m_visit <- m_visit[DATE %in% ts_season[M_SEASON > 0L, DATE], ]

            monit_syl <- m_visit[data.table::year(DATE) >= min(r_year) &
                        data.table::year(DATE) <= max(r_year),
                        .(SITE_ID = .SD[, unique(SITE_ID)]), by = M_YEAR]

            data.table::setkey(monit_syl, M_YEAR, SITE_ID)
            data.table::setkey(ts_season, M_YEAR)

            ts_season_site <- merge(ts_season, monit_syl, by.x = "M_YEAR", by.y = "M_YEAR", allow.cartesian = TRUE)

            data.table::setkey(ts_season_site, DATE, SITE_ID)
            data.table::setkey(m_visit, DATE, SITE_ID)

            ts_season_site <- ts_season_site[m_visit, COUNT := 0L] [M_SEASON == 0L & ANCHOR == 0L, COUNT := NA]

        return(ts_season_site)

    }


#' ts_monit_count_site
#' Generate a full time series of observed count, for all sites and each day since a starting and ending years of the defined time-series
#' @param m_season_visit data.table with potential absence based on visit data as returned by \code{ts_monit_site}
#' @param m_count data.table or data.frame with monitoring count data
#' @param sp integer or string for the species ID or name
#' @param DateFormat format used for the date in the count data, default="\%Y-\%m-\%d".
#' @return A data.table with a complete time-series where absences (0) and recorded counts for a specific species are implemented for each
#'         site.
#' @seealso \link{ts_monit_site}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export ts_monit_count_site
#'

ts_monit_count_site = function(m_season_visit, m_count, sp=1, DateFormat="%Y-%m-%d") {

        check_package('data.table')
        names(m_season_visit) <- toupper(names(m_season_visit))
        check_names(m_season_visit, c("DATE", "SITE_ID", "COUNT"))
        m_season_visit_copy <- data.table::copy(m_season_visit)

        m_count <- data.table::data.table(m_count)
        names(m_count) <- toupper(names(m_count))
        check_names(m_count,c("DATE", "SITE_ID", "SPECIES", "COUNT"))
        m_count[, DATE := data.table::as.IDate(as.Date(m_count$DATE, format=DateFormat))]

        if(!sp %in% m_count[, unique(SPECIES)]){
            stop(paste("Species", sp, "is not found in your dataset, check your \"sp\" argument."))
        }else{
            m_sp_count <- m_count[SPECIES %in% sp, ]
            data.table::setkey(m_sp_count, DATE, SITE_ID, SPECIES)
            data.table::setkey(m_season_visit, DATE, SITE_ID)
            spcount_site_series <- m_season_visit_copy[m_sp_count, COUNT := m_sp_count[, as.integer(COUNT)]]
            spcount_site_series[, SPECIES := sp]
        }

        return(spcount_site_series)

    }
