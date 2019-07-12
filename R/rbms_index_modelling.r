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
##      index_modelling: functions for abundance index model
##
##==========================================


#' fit_gam
#' fit a Generalized Additive Model to butterfly count data along a temporal variable and accounting for site efffect when multiple are available.
#' @param dataset_y data.table with filtered butterfly counts for species x over year y over all sites.
#' @param NbrSample integer inherited from \link{flight_curve}, default=100.
#' @param GamFamily string inherited from \link{flight_curve}, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param MaxTrial integer inherited from \link{flight_curve}, default=3.
#' @param SpeedGam Logical to use the \link[mgcv]{bam} method instead of the \link[mgcv]{gam} method.
#' @param OptiGam Logical to set use bam when data are larger than 100 and gam for smaller dataset
#' @param TimeUnit Character defining if the spline should be computed at the day 'd' or the week 'd'.
#' @param ... additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abudance, normalize to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. Thi is provide for one year 'y'.
#' @keywords gam
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export fit_gam
#'

fit_gam <- function(dataset_y, NbrSample = NULL, GamFamily = 'poisson', MaxTrial = 4,
                    SpeedGam = TRUE, OptiGam = TRUE, TimeUnit = 'd'){

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

            if(TimeUnit == 'd'){
              tp_col <- "trimDAYNO"
              dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)])
              sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, DAY_SINCE, -COUNT), ][dup, ]
            } else {
              tp_col <- "trimWEEKNO"
              dup <- !duplicated(sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)])
              sp_data_all <- sp_data_all[order(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE, -COUNT), ][dup, ]
            }

            mod_form <- as.formula(paste0("COUNT ~ s(",tp_col,",bs =\"cr\")",ifelse(sp_data_all[, uniqueN(SITE_ID)] > 1, "+ factor(SITE_ID)","")))

              if(isTRUE(SpeedGam)){
                gam_obj_site <- try(mgcv::bam(mod_form, data=sp_data_all, family=GamFamily), silent = TRUE)
              } else {
                gam_obj_site <- try(mgcv::gam(mod_form, data=sp_data_all, family=GamFamily), silent = TRUE)
              }

          tr <- tr + 1

        }

        if (class(gam_obj_site)[1] == "try-error") {
            print(paste("Error in fitting the RegionalGAM for species", as.character(sp_data_all$SPECIES[1]), "and year", sp_data_all$M_YEAR[1],
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
#' Compute the Normalized flight curve by fitting a spline in a Generalized Additive Model for one year 'y' to butterfly count data.
#' @param y integer of vector of years for wich to compute flight curve.
#' @param ts_season_count data.table with complete time series of count and season information returned by \link{ts_monit_count_site}
#' @param MinVisit integer setting the minimum number of visit required for a site to included in the computation, default=3.
#' @param MinOccur integer setting the minimum number of positive records (e.g. >= 1) observed over the year in a site default=2.
#' @param MinNbrSite integer setting the minimum number of site required to compute the flight curve, default=1.
#' @param NbrSample integer inherited from \link{flight_curve}, when set to 'NULL' (default), all site are considered in the GAM model
#' @param GamFamily string inherited from \link{flight_curve}, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param MaxTrial integer inherited from \link{flight_curve}, default=3.
#' @param SpeedGam Logical to use the \link[mgcv]{bam} method instead of the \link[mgcv]{gam} method.
#' @param OptiGam Logical to set use bam when data are larger than 200 and gam for smaller dataset
#' @param TimeUnit Character defining if the spline should be computed at the day 'd' or the week 'd'.
#' @return A list of lists, each containing three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abudance, normalize to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. This is provided for all year provided in 'y'.
#' @keywords gam, spline
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export get_nm
#'
#ts_season_count_d <- ts_season_count
get_nm <- function(y, ts_season_count, MinVisit, MinOccur, MinNbrSite, NbrSample, GamFamily, MaxTrial, SpeedGam, OptiGam, TimeUnit){

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
                          SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit)
  }
  return(f_curve_mod)
}

#' flight_curve
#' Compute the annual flight curve from butterfly count data collated across multiple sites.
#' @param ts_season_count data.table with complete time series of count and season information returned by \link{ts_monit_count_site}
#' @param NbrSample integer setting the maximum number of site to use to compute the flight curve, default=100.
#' @param MinVisit integer setting the minimum number of visit required for a site to included in the computation, default=3.
#' @param MinOccur integer setting the minimum number of positive records (e.g. >= 1) observed over the year in a site default=2.
#' @param MinNbrSite integer setting the minimum number of site required to compute the flight curve, default=1.
#' @param MaxTrial integer setting the maximum number of trial to reach convergence of the model, default=3.
#' @param FcMethod string defining the method to be used for computation of the flight curve, default='regionalGAM' (no other method is available yet).
#' @param GamFamily string setting the distribution of the error term in the GAM, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param CompltSeason Logical to restrict computation of flight curve for years where the complete season has been sampled, default=TRUE.
#' @param SelectYear integer to select a specific year to compute the flight curve, default=NULL.
#' @param SpeedGam Logical to use the \link[mgcv]{bam} method instead of the \link[mgcv]{gam} method.
#' @param OptiGam Logical to set use bam when data are larger than 100 and gam for smaller dataset.
#' @param KeepModel Logical to keep model output in a list object named \code{flight_curve_model}.
#' @param KeepModelData Logical to keep the data used for the GAM.
#' @param TimeUnit Character to define days 'd' or week 'w' as variable for the GAM.
#' @param ... additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **pheno**: a vector with annual flight curves \code{f_pheno} with expected relative abudance, normalize to sum to one over a full season,
#'         ii) **model**: a list of the resulting gam models \code{f_model} fitted on the count data for each year and iii) **data**: a data.table with the data used to fit the GAM model.
#' @keywords gam, flight curve
#' @seealso \link{fit_gam}
#' @seealso \link[mgcv]{gam}
#' @seealso \link[mgcv]{bam}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export flight_curve
#'

flight_curve <- function(ts_season_count, NbrSample = 100, MinVisit = 3, MinOccur = 2, MinNbrSite = 1, MaxTrial = 3, FcMethod = 'regionalGAM',
                            GamFamily = 'poisson', CompltSeason = TRUE, SelectYear = NULL, SpeedGam = TRUE, OptiGam = TRUE, KeepModel = TRUE, KeepModelData = TRUE,
                            TimeUnit = 'd', ...) {

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
                                                  OptiGam = OptiGam, TimeUnit = TimeUnit)

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
#' find nearest year with flight period
#' @param x year to find nearest flight period.
#' @param y years with availabe flight period.
#' @return z a value of the nearest year
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
#' Check for the flight curve of a specific year and if missing, impute the nearest available within a span of 5 years.Function used in \link{impute_count}.
#' @param sp_count_flight data.table with the flight curve, relative abundance (NM), for a specific year.
#' @param ts_flight_curve data.table with the flight curves, relative abundance (NM), for all years available for search as returned by \link{flight_curve}.
#' @param YearCheck integer or vector of year to check for nearest phenology, set internally in \link{impute_count2}.
#' @param YearLimit integer defining the range (+/- number of year) of year to look for a flight curve, if NULL no restriction is set.
#' @param TimeUnit Character defining if the spline should be computed at the day 'd' or the week 'd'.
#' @return A data.table with time series of the expected relative abundance of butterfly count per day (NM) for the year or the nearest year where
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


#' impute_count2
#' @param ts_season_count data.table with time series of counts for a specific species across all sites as returned by \link{ts_monit_count_site}.
#' @param ts_flight_curve data.table with the flight curves, relative abundance (NM), for a specific species as returned by \link{flight_curve}.
#' @param TimeUnit Character to define days 'd' or week 'w' as variable for the GAM.
#' @param sp integer or string for the species ID or name.
#' @param YearLimit integer defining the range (+/- number of year) of year to look for a flight curve, if NULL no restriction is set.
#' @param SelectYear integer to select a specific year to compute the flight curve, default=NULL.
#' @param CompltSeason Logical to restrict computation of flight curve for years where the complete season has been sampled, default=TRUE.
#' @return A data.table based on the entry count data, augmented with site indices 'SINDEX' and imputed weekly count 'IMPUTED_COUNT'.
#' @details Site indices can be extracted from the data.table returned from this function.The Site index is currently computed by adjusting the count by the proportion of the flight curve covered by the visits.
#' @keywords site index, flight curve
#' @seealso \link{fit_glm}, \link{fit_glm.nb}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export impute_count2
#'

impute_count2 <- function(ts_season_count, ts_flight_curve, TimeUnit, sp = NULL, YearLimit= NULL, SelectYear = NULL, CompltSeason = TRUE){

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

        if(length(YearCheck)>0){
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

        if(TimeUnit == 'd') {sp_count_flight[ , SINDEX := SINDEX / 7]}

     return(sp_count_flight)
  }


#' fit_glm
#' Fit a Generalized Linear Model (GLM) to predict daily count per site, using the regional flight curve as offset. Function used in \link{impute_count}.
#' @param sp_count_flight_y data.table with time series of the expected relative abundance of butterfly count per day (NM) for the year of interest,
#'        or the nearest available.
#' @param non_zero vector of sites with non-zero value to be included in the model (see detail)
#' @param FamilyGlm string for the distribution to be used for the error term in the GLM, inherited from \link{impute_count}, default='quasipoisson'.
#' @return A a list of objects, i) **sp_count_flight_y**: data.table with data and expected relative abundance of butterfly count per day (NM) for the year or the nearest year where
#'         a fligh curve (phenology) is available \code{sp_count_flight_y} and ii) **glm_obj_site**:a glm object for the model \code{glm_obj_site}.
#' @details GLM model is only fitted for site with observations, non-zero, as the abundance for sites where count are only zeros is expected to be null
#' @keywords flight curve
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export fit_glm
#'

fit_glm <- function(sp_count_flight_y, non_zero, FamilyGlm){

            if(sp_count_flight_y[SITE_ID %in% non_zero,uniqueN(SITE_ID)] > 1){
                glm_obj_site <- try(glm(COUNT ~ factor(SITE_ID) + offset(log(NM)) -1, data = sp_count_flight_y[SITE_ID %in% non_zero, ],
                family = FamilyGlm, control = list(maxit = 100)), silent = TRUE)
            } else {
                glm_obj_site <- try(glm(COUNT ~ offset(log(NM)) -1, data = sp_count_flight_y[SITE_ID %in% non_zero, ],
                family = FamilyGlm, control = list(maxit = 100)), silent = TRUE)
            }

            if (class(glm_obj_site)[1] == "try-error") {

                sp_count_flight_y[SITE_ID %in% non_zero, c("FITTED","COUNT_IMPUTED") := .(NA, NA)]

                print(paste("Computation of abundance indices for year",sp_count_flight_y[1, M_YEAR,],
                            "failed with the RegionalGAM, verify the data you provided for that year"))
                # next()

            } else {
                sp_count_flight_y[SITE_ID %in% non_zero, FITTED := predict.glm(glm_obj_site, newdata = sp_count_flight_y[SITE_ID %in% non_zero, ], type = "response")]
            }

            sp_count_flight_mod_y <- list(sp_count_flight_y = sp_count_flight_y, glm_obj_site = glm_obj_site)

        return(sp_count_flight_mod_y)
    }


#' fit_glm.nb
#' Fit a Generalized Linear Model (GLM) to predict daily count per site, using the regional flight curve as offset and the negative binomial error distribution
#' as implemented in the package link[MASS]{glm.nb}. Function used in \link{impute_count}.
#' @param sp_count_flight_y data.table with time series of the expected relative abundance of butterfly count per day (NM) for the year of interest,
#'        or the nearest available.
#' @param non_zero vector of sites with non-zero value to be included in the model (see detail)
#' @return A a list of two objects, i) **sp_count_flight_y**: data.table with data and expected relative abundance of butterfly count per day (NM) for the year or the nearest year where
#'         a fligh curve (phenology) is available \code{sp_count_flight_y} and ii) **glm_obj_site**: a glm object for the model \code{glm_obj_site}.
#' @details GLM model is only fitted for site with observations, non-zero, as the abundance for sites where count are only zeros is expected to be null
#' @keywords flight curve
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export fit_glm.nb
#'

fit_glm.nb <- function(sp_count_flight_y, non_zero){

            if(sp_count_flight_y[SITE_ID %in% non_zero,uniqueN(SITE_ID)] > 1){
                glm_obj_site <- try(MASS::glm.nb(COUNT ~ factor(SITE_ID) + offset(NM), data=sp_count_flight_y[SITE_ID %in% non_zero, ]), silent = TRUE)
            } else {
                glm_obj_site <- try(MASS::glm.nb(COUNT ~ offset(NM) -1, data = sp_count_flight_y[SITE_ID %in% non_zero, ]), silent = TRUE)
            }

            if (class(glm_obj_site)[1] == "try-error") {
                sp_count_flight_y[SITE_ID %in% non_zero, c("FITTED", "COUNT_IMPUTED") := .(NA,NA)]

                print(paste("Computation of abundance indices for year", sp_count_flight_y[1, M_YEAR, ],
                "failed with the RegionalGAM, verify the data you provided for that year"))

                # next()

            } else {
                sp_count_flight_y[SITE_ID %in% non_zero, FITTED := predict(glm_obj_site, newdata = sp_count_flight_y[SITE_ID %in% non_zero, ], type = "response")]
            }

            sp_count_flight_mod_y <- list(sp_count_flight_y = sp_count_flight_y, glm_obj_site = glm_obj_site)

        return(sp_count_flight_mod_y)
    }


#' impute_count
#' Fit a Generalized Linear Model (GLM) to predict daily count per site, using the regional flight curve as offset.
#' @param ts_season_count data.table with time series of counts for a specific species across all sites as returned by \link{ts_monit_count_site}.
#' @param ts_flight_curve data.table with the flight curves, relative abundance (NM), for a specific species as returned by \link{flight_curve}.
#' @param FamilyGlm string for the distribution to be used for the error term in the GLM, default='quasipoisson'.
#' @param CompltSeason logical to define if only years where a complete season is available should be modelled.
#' @param SelectYear vector of specific years of interest, can be a single value (e.g. 2015).
#' @param NearPheno Logical if flight curve from neirhest year should be used, if available within 5 years, default=TRUE
#' @param KeepModel Logical to keep model output in a list object named \code{imp_glm_model}
#' @return A a list of one or two objects, i) **impute_count**: data.table with observed and expected butterfly counts per day imputed based on the flight curve of the year or the nearest year where
#'         a phenology is available \code{sp_ts_season_count} and ii) **model**: a glm object for the GLM model \code{glm_model}.
#' @details GLM model is only fitted for site with observations, non-zero, as the abundance for sites where count are only zeros is expected to be null.
#' @keywords flight curve
#' @seealso \link{fit_glm}, \link{fit_glm.nb}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export impute_count
#'

impute_count <- function(ts_season_count, ts_flight_curve, FamilyGlm = quasipoisson(), CompltSeason = TRUE,
                                    SelectYear = NULL, KeepModel = TRUE, NearPheno = TRUE) {

        if(ts_season_count$SPECIES[1] != ts_flight_curve$SPECIES[1]){
            stop('Species in count data and flight curve must be the same!')
        }

        if(FamilyGlm[1]=='nb'){
            stop('Negative binomial distribution for GLM is not resolved, please use quasipoisson distribution instead')
        }

        check_package('data.table')

        if(isTRUE(CompltSeason)){
            ts_season_count <- ts_season_count[COMPLT_SEASON==1]
        }

        sp_ts_season_count <- data.table::copy(ts_season_count)
        data.table::setkey(sp_ts_season_count, DATE)
        data.table::setkey(ts_flight_curve, DATE)
        sp_count_flight <- merge(sp_ts_season_count, ts_flight_curve[, .(DATE, trimDAYNO, NM)], all.x=TRUE)
        data.table::setkey(sp_count_flight, M_YEAR, DATE, SITE_ID)

        glmMet <- "glm()"

        if( FamilyGlm[1]=='nb'){
            glmMet <- "glm.nb()"
        }

        if(is.null(SelectYear)){
            year_series <- ts_season_count[, unique(as.integer(M_YEAR))]
        } else {
            year_series <- ts_season_count[M_YEAR %in% SelectYear, unique(as.integer(M_YEAR))]
        }

        i_glm_model <- list()

        for(y in year_series){

            sp_count_flight_y <-  data.table::copy(sp_count_flight[as.integer(M_YEAR) == y, ])

            if(isTRUE(NearPheno)){
                sp_count_flight_y <- check_pheno(sp_count_flight_y, sp_count_flight)
            }

            if(sp_count_flight_y[is.na(NM) & trimDAYNO != 366, .N] > 0){
                print(paste("No glm model will be fitted for", sp_count_flight_y[1, M_YEAR]))
            } else {
                print(paste("Computing abundance indices for", sp_count_flight_y[1, SPECIES], "in", sp_count_flight_y[1, M_YEAR], "across",
                        sp_count_flight_y[,uniqueN(SITE_ID)], "sites, using", glmMet, ":", Sys.time()))
            }

            sp_count_flight_y[M_SEASON == 0L, COUNT := NA]
            sp_count_flight_y[M_SEASON != 0L & NM == 0, NM := 0.000001]
            non_zero <- sp_count_flight_y[, sum(COUNT, na.rm=TRUE), by=(SITE_ID)] [V1 > 0, SITE_ID]
            zero <- sp_count_flight_y[, sum(COUNT, na.rm=TRUE), by=(SITE_ID)] [V1 == 0, SITE_ID]

            if(length(non_zero) >= 1 & sp_count_flight_y[is.na(NM) & trimDAYNO != 366, .N] == 0){
                if(FamilyGlm[1] == 'nb'){
                    sp_count_flight_l <- fit_glm.nb(sp_count_flight_y, non_zero)
                } else {
                    sp_count_flight_l <- fit_glm(sp_count_flight_y, non_zero, FamilyGlm)
                }
            } else {
               sp_count_flight_l <- list(sp_count_flight_y = sp_count_flight_y,
                                    glm_obj_site = paste('No glm fitted for', sp_count_flight_y[1, M_YEAR]))
            }

            sp_count_flight_y <- sp_count_flight_l$sp_count_flight_y

            sp_count_flight_y[SITE_ID %in% zero, FITTED := 0]
            sp_count_flight_y[is.na(COUNT), COUNT_IMPUTED := FITTED][!is.na(COUNT), COUNT_IMPUTED := as.numeric(COUNT)] [M_SEASON == 0L, COUNT_IMPUTED := 0]

            data.table::setkey(sp_ts_season_count, SITE_ID, DAY_SINCE)
            data.table::setkey(sp_count_flight_y, SITE_ID, DAY_SINCE)

            if("FITTED" %in% names(sp_ts_season_count)){
                sp_ts_season_count[sp_count_flight_y, ':=' (trimDAYNO=i.trimDAYNO, NM=i.NM, FITTED=i.FITTED, COUNT_IMPUTED=i.COUNT_IMPUTED)]
            } else {
                sp_ts_season_count <- merge(sp_ts_season_count, sp_count_flight_y[, .(DAY_SINCE, SITE_ID, trimDAYNO, NM, FITTED, COUNT_IMPUTED)], all.x = TRUE)
            }

            if(isTRUE(KeepModel)){
                    glm_model <- list(sp_count_flight_l$glm_obj_site)
                    names(glm_model) <- paste0('glm_mod_', gsub(' ', '_', sp_count_flight_y[1, SPECIES]), '_', sp_count_flight_y[1, M_YEAR])
                    i_glm_model <- c(i_glm_model, glm_model)
                }
            }


        if(!is.null(SelectYear)){
            sp_ts_season_count <- sp_ts_season_count[M_YEAR %in% SelectYear, ]
        } else {
            sp_ts_season_count <- sp_ts_season_count
        }


        if(isTRUE(KeepModel)){
          return(list(impute_count = sp_ts_season_count, model = i_glm_model))
        } else {
          return(list(impute_count = sp_ts_season_count))
        }
}


#' butterfly_day
#' count cumulative butterfly count observed over one monitoring season.
#' @param sp_ts_season_count list of objects, including a data.table with observed and expected butterfly counts per day imputed for each sites based by \link{impute_count}.
#' @param WeekCount logical defining if butterfly day should be counted by week (default), using the 4th day of the week, or as a total of predicted daily counts.
#' @return data.table with total butterfly day observed per species, year, and site.
#' @keywords butterfly count
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{reto.schmucki@@mail.mcgill.ca}
#' @import data.table
#' @export butterfly_day
#'

butterfly_day <- function(sp_ts_season_count, WeekCount = TRUE){
            if (WeekCount == TRUE){
                b_day <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
            } else {
                b_day <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
            }

            data.table::setnames(b_day,"V1","BUTTERFLY_DAY")

        return(b_day)
    }
