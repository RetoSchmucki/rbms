##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   02.01.2018
##      Update:   13.04.2020
##
##      index_modelling: functions for abundance index model
##      removed this arguments #' @param ConLikelihood Logical to use the concentrated likelihood approach without site parameters (default = TRUE).
##==========================================


#' fit_gam
#' fit a Generalized Additive Model to butterfly count data along a temporal variable and accounting for site effect when multiple are available.
#' @param dataset_y data.table Filtered butterfly counts for species x over year y over all sites.
#' @param NbrSample integer Inherited from \link{flight_curve}, default=100.
#' @param GamFamily string Inherited from \link{flight_curve}, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param MaxTrial integer Inherited from \link{flight_curve}, default=3.
#' @param SpeedGam logical to use the \link[mgcv]{bam} method instead of the \link[mgcv]{gam} method.
#' @param OptiGam logical Set if the \link[mgcv]{bam} method should be used, default instead of the default \link[mgcv]{gam} method.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week.
#' @param MultiVisit string Function to apply for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abundance, normalize to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. This is provide for one year 'y'.
#' @keywords gam
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export fit_gam
#'

fit_gam <- function(dataset_y, NbrSample = NULL, GamFamily = 'poisson', MaxTrial = 4,
                    SpeedGam = TRUE, OptiGam = TRUE, ConLikelihood = TRUE,  TimeUnit = 'd', MultiVisit = "mean"){

        check_package('data.table')

        tr <- 1
        gam_obj_site <- c()

        while((tr == 1 | class(gam_obj_site)[1] == "try-error") & tr <= MaxTrial){

          if(!is.null(NbrSample)){
            if (dataset_y[, uniqueN(SITE_ID)] > NbrSample) {
                sp_data_all <- data.table::copy(dataset_y[SITE_ID %in% sample(unique(dataset_y[, SITE_ID]), NbrSample, replace = FALSE), ])
            } else {
                sp_data_all <- data.table::copy(dataset_y)
            }
          } else {
            sp_data_all <- data.table::copy(dataset_y)
          }

            if(isTRUE(OptiGam)){
                if(sp_data_all[, uniqueN(SITE_ID)] < 100){
                    SpeedGam <- FALSE
                }
            }

            gamMethod <-'gam()'
            if(isTRUE(SpeedGam)){
                gamMethod <- 'SpeedGAM [bam()]'
            }

            print(paste("Fitting the flight curve spline for species", as.character(sp_data_all$SPECIES[1]), "and year", sp_data_all$M_YEAR[1], "with",
                        sp_data_all[, uniqueN(SITE_ID)], "sites, using", gamMethod, ":", Sys.time(), "-> trial", tr))

            if(MultiVisit == "max") {
                if(TimeUnit == 'd'){
                  tp_col <- "trimDAYNO"
                  dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)])
                  sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), ][dup, ]
                } else {
                  tp_col <- "trimWEEKNO"
                  dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)])
                  sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), ][dup, ]
                }
            } else {
                if (TimeUnit == "d") {
                  tp_col <- "trimDAYNO"
                  sp_data_all[ , meanCOUNT := ceiling(mean(COUNT, na.rm = TRUE)), by = .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)]
                  dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, meanCOUNT), .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)])
                  sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, meanCOUNT), ][dup, ]
                  sp_data_all <- sp_data_all[ , COUNT := meanCOUNT][ , meanCOUNT := NULL]
                } else {
                  tp_col <- "trimWEEKNO"
                  sp_data_all[, meanCOUNT := ceiling(mean(COUNT, na.rm = TRUE)), by = .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)]
                  dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, meanCOUNT), .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)])
                  sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, meanCOUNT), ][dup, ]
                  sp_data_all <- sp_data_all[, COUNT := meanCOUNT][, meanCOUNT := NULL]
                }
            }

           # if(isTRUE(ConLikelihood)){
           # mod_form <- as.formula(paste0("COUNT ~ s(", tp_col,", bs =\"cr\")"))
           # } else {
            mod_form <- as.formula(paste0("COUNT ~ s(", tp_col,", bs =\"cr\")", ifelse(sp_data_all[, uniqueN(SITE_ID)] > 1, "+ factor(SITE_ID)", "")))
           # }

              if(isTRUE(SpeedGam)){
                gam_obj_site <- try(mgcv::bam(mod_form, data=sp_data_all, family=GamFamily), silent = TRUE)
              } else {
                gam_obj_site <- try(mgcv::gam(mod_form, data=sp_data_all, family=GamFamily), silent = TRUE)
              }

          tr <- tr + 1

        }

        if (class(gam_obj_site)[1] == "try-error") {
            print(paste("Error in fitting the GAM for species", as.character(sp_data_all$SPECIES[1]), "and year", sp_data_all$M_YEAR[1],
                    "; Model did not converge after", tr, "trials"))
            sp_data_all[, c("FITTED", "NM") := .(NA, NA)]
        } else {
            colnames <- c(tp_col, "SITE_ID")
            pred_data <- sp_data_all[, ..colnames]
            sp_data_all[, FITTED := mgcv::predict.gam(gam_obj_site, newdata = pred_data, type = "response")]
            sp_data_all[M_SEASON == 0L, FITTED := 0]

            if(sum(is.infinite(sp_data_all[, FITTED])) > 0){
                sp_data_all[, c("FITTED", "NM") := .(NA, NA)]
            } else {
                sp_data_all[, SITE_SUM := sum(FITTED), by = SITE_ID]
                sp_data_all[, NM := round(FITTED / SITE_SUM, 5)]
            }
        }

        sid_1 <- sp_data_all[, SITE_ID][1]
        f_curve <- sp_data_all[SITE_ID == sid_1, ][ , c("COUNT", "SITE_ID","FITTED", "SITE_SUM") := NULL]
        data.table::setkey(f_curve)
        sp_data_all[ , c("FITTED", "SITE_SUM", "NM") := NULL]

        f_curve_mod <- list(f_curve = f_curve, f_model = gam_obj_site, f_data = sp_data_all)

    return(f_curve_mod)
}


#' get_nm
#' Compute the normalized flight curve by fitting a spline in a Generalized Additive Model for one year 'y' to butterfly count data.
#' @param y integer Vector of years for which to compute the flight curve.
#' @param ts_season_count data.table Time-series of count and season information returned by \link{ts_monit_count_site}
#' @param MinVisit integer The minimum number of visits required for a site to be included in the computation, default=3.
#' @param MinOccur integer The minimum number of positive records (e.g. >= 1) observed over the year in a site default=2.
#' @param MinNbrSite integer The minimum number of sites required to compute the flight curve, default=1.
#' @param NbrSample integer Value inherited from \link{flight_curve}, when set to 'NULL' (default), all site are considered in the GAM model
#' @param GamFamily string Value inherited from \link{flight_curve}, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param MaxTrial integer Value inherited from \link{flight_curve}, default=3.
#' @param SpeedGam logical Set if the \link[mgcv]{bam} method should be used, default instead of the default \link[mgcv]{gam} method.
#' @param OptiGam logical Set the use the \link[mgcv]{bam} method when data are larger than 200 and gam for smaller datasets
#' @param TimeUnit character Time-step for which the spline should be computed, 'd' day or 'w' week.
#' @param MultiVisit string Function for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @return A list of lists, each list containing three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abundance, normalised to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. This is provided for all year provided in 'y'.
#' @keywords gam, spline
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export get_nm
#'

get_nm <- function(y, ts_season_count, MinVisit, MinOccur, MinNbrSite, NbrSample, GamFamily, MaxTrial, SpeedGam, OptiGam, TimeUnit, MultiVisit){

  dataset_y <- ts_season_count[as.integer(M_YEAR) == y, ]
  visit_occ_site <- merge(dataset_y[!is.na(COUNT) & ANCHOR == 0L, .N, by=SITE_ID],
                          dataset_y[!is.na(COUNT) & ANCHOR == 0L & COUNT > 0, .N, by=SITE_ID],
                          by="SITE_ID", all=TRUE)
  dataset_y <- data.table::copy(dataset_y[SITE_ID %in%
                                visit_occ_site[N.x >= MinVisit & N.y >= MinOccur, SITE_ID],])

  if(TimeUnit == 'd'){
    tp_col <- "trimDAYNO"
  } else {
    tp_col <- "trimWEEKNO"
  }

  if(dataset_y[, uniqueN(SITE_ID)] < MinNbrSite){

    f_curve  <- unique(ts_season_count[as.integer(M_YEAR) == y, ][ , SITE_ID := NULL][, COUNT := NULL])[, NM := NA]
    f_curve_mod <- list(f_curve=f_curve[order(get(tp_col)),], f_model=list(NA), f_data=data.table(NA))
    print(paste("You have not enough sites with observations for estimating the flight curve for species", as.character(ts_season_count$SPECIES[1]), "in", unique(ts_season_count[as.integer(M_YEAR) == y, M_YEAR])))
  } else {
    f_curve_mod <- fit_gam(dataset_y, NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial,
                          SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit, MultiVisit = MultiVisit) #ConLikelihood = ConLikelihood,
  }
  return(f_curve_mod)
}


#' flight_curve
#' Compute the annual flight curve from butterfly count data collated across sites.
#' @param ts_season_count data.table Time-series of counts and season information returned by \link{ts_monit_count_site}
#' @param NbrSample integer The maximum number of sites to use for computing the flight curve, default=100.
#' @param MinVisit integer The minimum number of visits required for a site to be included, default=3.
#' @param MinOccur integer The minimum number of positive records (e.g. >= 1) required in one year for a site to be included, default=2.
#' @param MinNbrSite integer The minimum number of sites required to compute a flight curve, default=1.
#' @param MaxTrial integer The maximum number of trials for model convergence, default=3.
#' @param GamFamily string The distribution of the error term used in the GAM, default='poisson', but can be the negative-binomial 'nb' or 'quasipoisson'.
#' @param CompltSeason Logical to restrict computation of flight curve for years where the complete season has been sampled, default=TRUE.
#' @param SelectYear integer Select a specific year to compute the flight curve, default=NULL.
#' @param SpeedGam Logical Set if the \link[mgcv]{bam} method should be used, default instead of the default \link[mgcv]{gam} method.
#' @param OptiGam Logical Set the use of the \link[mgcv]{bam} method when data are larger than 200 and gam for smaller datasets
#' @param KeepModel Logical to keep model output in a list object named \code{flight_curve_model}.
#' @param KeepModelData Logical to keep the data used for the GAM.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week.
#' @param MultiVisit string Function to apply for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **pheno**: a vector with annual flight curves \code{f_pheno} with expected relative abundance, normalize to sum to one over a full season,
#'         ii) **model**: a list of the resulting gam models \code{f_model} fitted on the count data for each year and iii) **data**: a data.table with the data used to fit the GAM model.
#' @keywords gam, flight curve
#' @seealso \link{fit_gam}
#' @seealso \link[mgcv]{gam}
#' @seealso \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export flight_curve
#'

flight_curve <- function(ts_season_count, NbrSample = 100, MinVisit = 3, MinOccur = 2, MinNbrSite = 1, MaxTrial = 3,
                         GamFamily = 'poisson', CompltSeason = TRUE, SelectYear = NULL, SpeedGam = TRUE,
                         OptiGam = TRUE, KeepModel = TRUE, KeepModelData = TRUE, 
                         TimeUnit = 'd', MultiVisit = "mean", ...) {

        check_package('data.table')

        names(ts_season_count) <- toupper(names(ts_season_count))
        check_names(ts_season_count, c("COMPLT_SEASON", "M_YEAR", "SITE_ID", "SPECIES", "DATE", "WEEK", "WEEK_DAY", "DAY_SINCE", "WEEK_SINCE", "M_SEASON", "COUNT", "ANCHOR"))

        if(isTRUE(CompltSeason)){
            ts_season_count <- ts_season_count[COMPLT_SEASON == 1, ]
        }

        if(is.null(SelectYear)){
            year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, unique(as.integer(M_YEAR))]
        } else {
            year_series <- ts_season_count[M_YEAR %in% SelectYear, unique(as.integer(M_YEAR))]
        }

        if(length(year_series) == 0) stop(paste0(" No count data found for year ", SelectYear))

        if(TimeUnit == 'd'){
            tp_col <- "trimDAYNO"
            dup <- !duplicated(ts_season_count[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)])
            ts_season_count <- ts_season_count[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), ][dup, ]
            ts_season_count[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1, by = M_YEAR]
            ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, WEEK_SINCE, DAY_SINCE, trimDAYNO, M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]
        } else {
            tp_col <- "trimWEEKNO"
            dup <- !duplicated(ts_season_count[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)])
            ts_season_count <- ts_season_count[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), ][dup, ]
            ts_season_count[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
            ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, WEEK, WEEK_SINCE, trimWEEKNO, M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]
        }

        result_fc <- lapply(year_series, get_nm, ts_season_count=ts_season_count, MinVisit=MinVisit, MinOccur=MinOccur, MinNbrSite=MinNbrSite,
                                                  NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial, SpeedGam = SpeedGam,
                                                  OptiGam = OptiGam, TimeUnit = TimeUnit, MultiVisit = MultiVisit) #ConLikelihood = ConLikelihood, 

        result_fcurve <- data.table::rbindlist(lapply(result_fc, function(x) x$f_curve), fill = TRUE)
        result_fdata <- data.table::rbindlist(lapply(result_fc, function(x) x$f_data), fill = TRUE)
        result_fmodel <- lapply(result_fc, function(x) x$f_model)

        if(!isTRUE(KeepModelData) & !isTRUE(KeepModel)){
          result_fc <- list(pheno = result_fcurve)
        }

        if(!isTRUE(KeepModelData) & isTRUE(KeepModel)){
          result_fc <- list(pheno = result_fcurve, model = result_fmodel)
        }

        if(isTRUE(KeepModelData) & isTRUE(KeepModel)){
          result_fc <- list(pheno = result_fcurve, model = result_fmodel, data = result_fdata)
        }

    return(result_fc)
}


#' get_nny
#' find the nearest year with a computed flight curve
#' @param x integer Year for which to find the nearest flight curve.
#' @param y vector Years with available flight curve.
#' @return z integer The value of the nearest year with flight curve.
#' @keywords flight curve
#' @seealso \link{check_pheno}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export check_pheno
#'

get_nny <- function(x, y) {
      z <- which.min(abs(x - y))
      if(length(z) > 1)
      z <- sample(z, 1)
      return(z)
      }


#' check_pheno
#' Check for the flight curve of a specific year. If the specific year is missing, use the nearest year available within a 5-year period to impute missing count. Function used in \link{impute_count}.
#' @param sp_count_flight data.table Object with the flight curve, relative abundance (NM), for a specific year.
#' @param ts_flight_curve data.table Object with the all flight curves and relative abundance (NM) for the years available for the search as returned by \link{flight_curve}.
#' @param YearCheck integer or vector Year to check for nearest flight curve, set internally in \link{impute_count}.
#' @param YearLimit integer Define the span of years (+/- number of year) to look for a flight curve, if NULL no restriction is set.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week.
#' @return A data.table with time-series of the expected relative abundance of butterfly count per day (NM) for the year or the nearest year where
#'         phenology is available.
#' @keywords flight curve
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export check_pheno
#'

check_pheno <- function(sp_count_flight, ts_flight_curve, YearCheck, YearLimit, TimeUnit){

       sp_count_flight_y <- sp_count_flight[M_YEAR == YearCheck, ]

       if(any(!is.na(ts_flight_curve[!is.na(NM), .N, by = M_YEAR]$N)) == 0){
            warning("No flight curve available for this species!")
        }

        if(TimeUnit == 'd'){
            tp_col <- "trimDAYNO"
        } else {
            tp_col <- "trimWEEKNO"
        }

        fcol <- c("M_YEAR", tp_col, "NM")

        count_y <- as.integer(as.character(unique(sp_count_flight_y[ , M_YEAR])))
        fc_y <- as.integer(as.character(unique(ts_flight_curve[!is.na(NM), M_YEAR])))
        fc_year <- fc_y[unlist(lapply(count_y, get_nny, y = fc_y))]

        alt_flight <- unique(ts_flight_curve[as.integer(as.character(M_YEAR)) == fc_year, fcol, with = FALSE])

        if(is.null(YearLimit)){
          warning(paste("We used the flight curve of", alt_flight[1, M_YEAR], "to compute abundance indices for year", sp_count_flight_y[1, M_YEAR, ]))
          data.table::setnames(alt_flight, 'NM', 'NMnew')
          alt_flight[, M_YEAR := NULL]
          data.table::setkeyv(sp_count_flight_y, tp_col)
          data.table::setkeyv(alt_flight, tp_col)
          sp_count_flight_y <- merge(sp_count_flight_y, alt_flight, by=tp_col, all.x=TRUE)
          sp_count_flight_y[, NM := NMnew][, NMnew := NULL]
        } else {
          if(abs(count_y - fc_year) > YearLimit){
            print(paste("No reliable flight curve available within a ",YearLimit," year horizon of", sp_count_flight_y[1, M_YEAR, ]))
          } else {
            warning(paste("We used the flight curve of", alt_flight[1, M_YEAR], "to compute abundance indices for year", sp_count_flight_y[1, M_YEAR, ]))
            data.table::setnames(alt_flight, 'NM', 'NMnew')
            alt_flight[, M_YEAR := NULL]
            data.table::setkeyv(sp_count_flight_y, tp_col)
            data.table::setkeyv(alt_flight, tp_col)
            sp_count_flight_y <- merge(sp_count_flight_y, alt_flight, by=tp_col, all.x=TRUE)
            sp_count_flight_y[, NM := NMnew][, NMnew := NULL]
          }
       }

      return(sp_count_flight_y)

    }


#' impute_count
#' @param ts_season_count data.table Time-series of counts for a specific species across all sites as returned by \link{ts_monit_count_site}.
#' @param ts_flight_curve data.table Flight curves and relative abundances (NM) for a specific species as returned by \link{flight_curve}.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week.
#' @param sp integer or string Species ID or name.
#' @param YearLimit integer Define the span of years (+/- number of year) to look for a flight curve, if NULL no restriction is set.
#' @param SelectYear integer Select a specific year to compute the flight curve, default=NULL.
#' @param CompltSeason logical Restrict computation of flight curve for years where the complete season was sampled, default=TRUE.
#' @return A data.table based on the entry count data, augmented with site indices 'SINDEX' and imputed weekly count 'IMPUTED_COUNT'.
#' @details Site indices can be extracted from the data.table returned by this function. The site index is currently computed by adjusting the count by the proportion of the flight curve covered by the visits.
#' @keywords site index, flight curve
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export impute_count
#'

impute_count <- function(ts_season_count, ts_flight_curve, TimeUnit, sp = NULL, YearLimit= NULL, SelectYear = NULL, CompltSeason = TRUE){

        check_package('data.table')

        if(!is.null(SelectYear)){
          ts_season_count <- ts_season_count[M_YEAR == SelectYear, ]
        }

        if(!is.null(sp)){
          ts_season_count <- ts_season_count[SPECIES == sp, ]
        }

        if(!ts_season_count$SPECIES[1] %in% unique(ts_flight_curve$SPECIES)){
            stop('Species in count data must be in the flight curve data!')
        }

        if(TimeUnit == 'd'){
            tp_col <- "trimDAYNO"
            dup <- !duplicated(ts_season_count[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)])
            ts_season_count <- ts_season_count[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), ][dup, ]
            ts_season_count[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1, by = M_YEAR]
            ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, WEEK_SINCE, DAY_SINCE, trimDAYNO, M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]
        } else {
            tp_col <- "trimWEEKNO"
            dup <- !duplicated(ts_season_count[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)])
            ts_season_count <- ts_season_count[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), ][dup, ]
            ts_season_count[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
            ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, WEEK, WEEK_SINCE, trimWEEKNO, M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]
        }

        keycol <- c("SPECIES", "M_YEAR", tp_col)
        data.table::setkeyv(ts_season_count, keycol)
        data.table::setkeyv(ts_flight_curve, keycol)
        mcol <- c(keycol, "NM")

        sp_count_flight <- merge(ts_season_count, ts_flight_curve[!duplicated(ts_flight_curve[, mcol, with = FALSE]), mcol, with = FALSE], all.x=TRUE)

        YearCheck <- as.integer(as.character(unique(sp_count_flight[ is.na(NM), M_YEAR])))

        if(length(YearCheck) > 0){
          a <- lapply(YearCheck, check_pheno, sp_count_flight=sp_count_flight, ts_flight_curve=ts_flight_curve, YearLimit=YearLimit, TimeUnit = TimeUnit)
          sp_count_flight <- rbind(data.table::rbindlist(a), sp_count_flight[!M_YEAR %in% data.table::rbindlist(a)[!is.na(NM), .N, by= M_YEAR][N>0, M_YEAR], ])
        }

        total_count <- sp_count_flight[M_SEASON != 0 , sum(COUNT, na.rm = TRUE), by = .(M_YEAR, SITE_ID)][, TOTAL_COUNT := V1][, V1 := NULL]
        total_nm <- sp_count_flight[M_SEASON != 0 & !is.na(COUNT) , sum(NM, na.rm = TRUE), by = .(M_YEAR, SITE_ID)][, TOTAL_NM := V1][, V1 := NULL]
        setkey(total_count, M_YEAR, SITE_ID)
        setkey(total_nm, M_YEAR, SITE_ID)
        setkey(sp_count_flight, M_YEAR, SITE_ID)

        sp_count_flight <- merge(sp_count_flight, merge(total_count, total_nm, all.x=TRUE)[, SINDEX := TOTAL_COUNT / TOTAL_NM], all.x = TRUE)
        sp_count_flight[, IMPUTED_COUNT:= COUNT][M_SEASON != 0 & is.na(COUNT), IMPUTED_COUNT := as.integer(round(SINDEX*NM))]

        if(TimeUnit == 'd') {
          sp_count_flight[ , SINDEX := SINDEX / 7]
          week_count <- unique(sp_count_flight[M_SEASON != 0 & !is.na(COUNT), .(M_YEAR, SITE_ID, WEEK)])[ , WEEK_COUNT := 1]
          setkey(week_count, SITE_ID, M_YEAR, WEEK)
          setkey(sp_count_flight, SITE_ID, M_YEAR, WEEK)
          sp_count_flight <- merge(sp_count_flight, week_count, all.x = TRUE)
          sp_count_flight[ , WEEK_NM := mean(NM * WEEK_COUNT, na.rm = TRUE), by = .(M_YEAR, SITE_ID, WEEK)][ ,
                                             TOTAL_NM := sum(WEEK_NM, na.rm=TRUE), by = .(M_YEAR, SITE_ID)][ , WEEK_NM := NULL][, WEEK_COUNT := NULL]
        }

     return(sp_count_flight)
  }


#' site_index
#' Extract abundance indices per site and year based on flight curve imputation.
#' @param butterfly_count data.table Observed and imputed weekly or daily counts and the estimated total counts (SINDEX) computed by the function impute_count2().
#' @param MinFC numeric Value between 0 and 1 to define the threshold for the proportion of the flight curve covered by the visits, if NULL all site-year available indices are returned.
#' @return data.table Estimated annual abundance index and the proportion of the flight curve covered by the visit - total weekly or daily count over the entire monitoring season.
#' @keywords butterfly count
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export site_index
#'

site_index <- function(butterfly_count, MinFC = NULL){

      if(is.null(MinFC)){
        sindex <- unique(butterfly_count[!is.na(SINDEX), .(SPECIES, SITE_ID, M_YEAR, TOTAL_NM, SINDEX)])
      } else {
        if(!MinFC >= 0 & MinFC <= 1) stop('MinFC must be between 0 and 1, or NULL')
        sindex <- unique(butterfly_count[!is.na(SINDEX) & TOTAL_NM >= MinFC, .(SPECIES, SITE_ID, M_YEAR, TOTAL_NM, SINDEX)])
      }

    return(sindex)
}


#' collated_index_old
#' compute a collated index from the site indices, using a Generalized Linear Model.
#' @param site_indices data.table or data.frame with site indices per year and proportion of flight curve covered by the monitoring.
#' @param GlmWeight vector of weight used in the GLM.
#' @param GlmFamily family used for the GLM model.
#' @return a list of three objects, a vector of site, a glm model object, and a vector of collated indices per year.
#' @keywords annual index
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export collated_index_old
#'

collated_index_old <- function(site_indices, GlmWeight = NULL, GlmFamily = poisson()){
      
      site_indices <- setDT(site_indices)

      if(is.null(GlmWeight)){
        col_mod <- try(speedglm::speedglm(round(SINDEX) ~ factor(M_YEAR) + factor(SITE_ID) - 1, data = site_indices, family = GlmFamily), silent = TRUE)
      } else {
        col_mod <- try(speedglm::speedglm(round(SINDEX) ~ factor(M_YEAR) + factor(SITE_ID) - 1, data = site_indices, family = GlmFamily, weigth = GlmWeight), silent = TRUE)
      }

      if(class(col_mod)[1] == "try-error"){
        if(is.null(GlmWeight)){
          col_mod <- try(glm(round(SINDEX) ~ factor(M_YEAR) + factor(SITE_ID) - 1, data = site_indices, family = GlmFamily), silent = TRUE)
        } else {
          col_mod <- try(glm(round(SINDEX) ~ factor(M_YEAR) + factor(SITE_ID) - 1, data = site_indices, family = GlmFamily, weigth = GlmWeight), silent = TRUE)
        }
      }

      if(class(col_mod)[1] == "try-error") stop('Unable to fit GLM to estimate a collated index')

      yr_coefs = coef(summary(col_index))
      res <- yr_coefs[grepl("M_YEAR",row.names(yr_coefs)), c("Estimate","Std. Error")]

      res <- list(site = site_indices$SITE_ID, glm_model = col_mod, collated_index = res)

    return(res)
}

#' collated_index
#' compute a collated index from the site indices, using a Generalized Linear Model.
#' @param data data.table or data.frame Site indices per year and proportion of flight curve covered by the monitoring.
#' @param s_sp string Species name to be found in the data.
#' @param sindex_value string Name of the response variable to be used the computation of the collated index, default is "SINDEX", but could be other standardized values.
#' @param bootID integer Identify the n^th bootstrap.
#' @param boot_ind data.table Index of the bootstrap and the site ids of the specific bootstrap sample.
#' @param glm_weights logical Use the proportion of the flight curve sampled as weight for the collated index, if FALSE the function uses equal weights.
#' @param rm_zero logical Remove the sites where species was not observed to speed-up the fit of the GLM without altering the output, default TRUE.
#' @return a list of two objects, a vector of site, a glm model object, and a vector of collated indices per year.
#' @keywords annual index
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export collated_index
#'

collated_index <- function(data, s_sp, sindex_value = "SINDEX", bootID=NULL, boot_ind=NULL, glm_weights=TRUE, rm_zero=TRUE){

  data <- setDT(data)
  data = data[SPECIES == s_sp, ]
  y = as.numeric(as.character(unique(data$M_YEAR)))

  if(is.null(boot_ind)){
    site_id_ <- unique(data[, SITE_ID])
    boot_ind_ <- matrix(NA, nrow = 1, ncol = length(site_id_))
    boot_ind_[1, ] <- seq_len(length(site_id_))
    boot_ind <- list(site_id = site_id_, boot_ind = boot_ind_)
  } else {
    boot_ind = boot_ind
  }

  if(is.null(bootID)){
    bootID <- 0
  }else{
    booID <- bootID
  }

  if(bootID == 0){
    if(is.null(boot_ind)){
    a <- data.table(uID = as.character(paste0("s_",seq_len(data[, uniqueN(SITE_ID)]))), SITE_ID = as.character(unique(data[, SITE_ID])))
    }else{
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID+1, ]))), SITE_ID = as.character(boot_ind$site_id[seq_along(boot_ind$boot_ind[bootID+1, ])]))
    }
  } else {
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID, ]))), SITE_ID = as.character(boot_ind$site_id[boot_ind$boot_ind[bootID, ]]))
  }

  s_y <- data.table(expand.grid(a$uID, y[order(y)]))
  names(s_y) <- c("uID", "M_YEAR")
  s_y[ , uID:=as.character(uID)]
  s_y[ , M_YEAR:=as.numeric(as.character(M_YEAR))]
  data[ , SITE_ID:=as.character(SITE_ID)]
  data[ , M_YEAR:=as.numeric(as.character(M_YEAR))]

  setkey(a, uID); setkey(s_y, uID)
  s_y <- merge(s_y, a, all.x=TRUE)

  setkey(s_y, SITE_ID, M_YEAR); setkey(data, SITE_ID, M_YEAR)

  if(sindex_value != 'SINDEX'){
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)][, (sindex_value) := data[, get(sindex_value)]])
    } else {
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)])
    }

  sum_data <- data.table( BOOTi = bootID,
                          M_YEAR = y[order(y)],
                          NSITE = data_boot[, uniqueN(SITE_ID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_OBS = data_boot[, sum(SINDEX>0), by = M_YEAR][order(M_YEAR), V1])

  setkey(sum_data, M_YEAR)

  if(nrow(sum_data[NSITE_OBS > 0,]) != 0){

    if(isTRUE(rm_zero)){
      zero_site <- data_boot[ , V1:=sum(SINDEX>0), uID][V1==0, uID]
      data_boot <- data_boot[!uID %in% zero_site, ]
    }

    pred_data_boot <- expand.grid(unique(data_boot$uID), unique(data_boot$M_YEAR))
    names(pred_data_boot) <- c("uID", "M_YEAR")

    if (data_boot[, uniqueN(M_YEAR)] == 1 & data_boot[, uniqueN(SITE_ID)] == 1) warning(paste0(sindex_value, " available for only one site and one year: no collated index computed for ", s_sp, " - ", unique(data_boot$M_YEAR)), call. = FALSE)
    
    if (data_boot[, uniqueN(M_YEAR)] == 1) {
      mod_form <- as.formula(paste0("I(round(", sindex_value, ")) ~ factor(uID) - 1"))
    } else {
      mod_form <- as.formula(ifelse(data_boot[, uniqueN(SITE_ID)] > 1, paste0("I(round(", sindex_value, ")) ~ factor(uID) + factor(M_YEAR) - 1"), paste0("I(round(", sindex_value, "))) ~ factor(M_YEAR) - 1")))
    }

    if(!isTRUE(glm_weights)){
      data_boot[ , weights := 1]
    }else{
      ## make sure that column TOTAL_NM exists
      data_boot[ , weights := TOTAL_NM]
    }

    col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights), silent = TRUE)
      if (class(col_index)[1] == "try-error") {
        col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights, method = "qr"), silent = TRUE)
        if (class(col_index)[1] == "try-error") {
          col_index <- try(glm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights), silent = TRUE)
            if (class(col_index)[1] == "try-error") {
              res <- sum_data[, COL_INDEX := NA]
              return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_OBS, COL_INDEX)], site_id = a$SITE_ID ))
            }
        }
     }

    pred_data_boot$fitted <- predict(col_index, newdata = pred_data_boot, type = "response")

    co_index_res <- data.table(pred_data_boot)[order(M_YEAR), sum(fitted)/length(unique(data$SITE_ID)), by = M_YEAR]
    setnames(co_index_res, "V1", "COL_INDEX"); setkey(co_index_res, M_YEAR)
    res <- merge(sum_data, co_index_res, all.x = TRUE)[is.na(COL_INDEX), COL_INDEX := 0]

  } else {
    res <- sum_data[, COL_INDEX := 0]
  }

  return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_OBS, COL_INDEX)], site_id = a$SITE_ID ))
}


#' boot_sample
#' Generate n bootstrap sample of the monitoring sites to be used for each iteration.
#' @param data data.table or data.frame Data with all site id.
#' @param boot_n integer The number of bootstrap samples to be generated.
#' @return A list with site id and bootstrap indices for n bootstrap sample.
#' @keywords bootstrap collated index
#' @seealso \link{impute_count}, \link{collated_index}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export boot_sample
#'

boot_sample <- function(data, boot_n = 1000){

  site_id <- unique(data[, SITE_ID])
  boot_n <- boot_n
  boot_ind <- matrix(NA, nrow = boot_n, ncol = length(site_id))

    for(i in seq_len(boot_n)){
      boot_ind[i, ] <- sample(seq_len(length(site_id)), replace = TRUE)
    }

  boot_ind <- list(site_id = site_id, boot_ind = boot_ind)
  return(boot_ind)
}
