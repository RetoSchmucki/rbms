##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   02.01.2018
##      Update:   22.04.2021
##
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
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param mod_form string with formula to be passed to the gam model, default null.
#' @param tp_col string or vector of string with additional variable used in the gam model, default null.
#' @param verbose a logical indicating if some “progress report” should be given. 
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abundance, normalize to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. This is provide for one year 'y'.
#' @keywords gam
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @importFrom stats as.formula coef fitted glm poisson predict weights
#' @export fit_gam
#'

fit_gam <- function(dataset_y, NbrSample = NULL, GamFamily = 'poisson', MaxTrial = 4,
                    SpeedGam = TRUE, OptiGam = TRUE, TimeUnit = 'd', MultiVisit = 'mean',
                    weekday = 3, mod_form = NULL, tp_col = NULL, verbose = TRUE, ...){


        check_package('data.table')
        tr <- 1
        gam_obj_site <- c()

        while((tr == 1 | inherits(gam_obj_site,"try-error")) & tr <= MaxTrial){

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

            if(verbose){
              message(paste("Fitting the flight curve spline for species", as.character(sp_data_all$SPECIES[1]), "and year", sp_data_all$M_YEAR[1], "with",
                            sp_data_all[, uniqueN(SITE_ID)], "sites, using", gamMethod, ":", Sys.time(), "-> trial", tr))
            }

            if(is.null(mod_form)){
            mod_form <- as.formula(paste0("COUNT ~ s(", tp_col,", bs =\"cr\")", ifelse(sp_data_all[, uniqueN(SITE_ID)] > 1, "+ factor(SITE_ID)", "")))
            } else {
              mod_form <- as.formula(mod_form)
            }

              if(isTRUE(SpeedGam)){
                gam_obj_site <- try(mgcv::bam(mod_form, data=sp_data_all, family=GamFamily, ...), silent = TRUE)
              } else {
                gam_obj_site <- try(mgcv::gam(mod_form, data=sp_data_all, family=GamFamily, ...), silent = TRUE)
              }

          tr <- tr + 1

        }

        if (inherits(gam_obj_site, "try-error")) {
            if(verbose) {
              message(paste("Error in fitting the GAM for species", as.character(sp_data_all$SPECIES[1]), "and year", sp_data_all$M_YEAR[1],
                            "; Model did not converge after", tr, "trials"))
            }
            sp_data_all[, c("FITTED", "NM") := .(NA, NA)]
        } else {
            colnames <- names(sp_data_all)
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
        f_curve <- data.table::copy(sp_data_all)
        
        ifelse(is.null(tp_col),
          f_curve <- f_curve[SITE_ID == sid_1, , c("COUNT", "SITE_ID","FITTED", "SITE_SUM") := NULL],
          f_curve <- f_curve[, c("SITE_SUM") := NULL])  
        # keep fitted values for testing against observed counts

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
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param mod_form string with formula to be passed to the gam model, default null.
#' @param tp_col string or vector of string with additional variable used in the gam model, default null.
#' @param verbose a logical indicating if some “progress report” should be given. 
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list of lists, each list containing three objects, i) **f_curve**: a data.table with the flight curve \code{f_curve} with expected relative abundance, normalised to sum to one over a full season,
#'         ii) **f_model**: the resulting gam model \code{f_model} fitted on the count data and iii) **f_data**: a data.table with the data used to fit the GAM model. This is provided for all year provided in 'y'.
#' @keywords gam, spline
#' @seealso \link{flight_curve}, \link[mgcv]{gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export get_nm
#'

get_nm <- function(y, ts_season_count, MinVisit, MinOccur, MinNbrSite, NbrSample, GamFamily, MaxTrial, 
                  SpeedGam, OptiGam, TimeUnit, MultiVisit, weekday, mod_form, tp_col, verbose = TRUE, ...){

  dataset_y <- ts_season_count[as.integer(M_YEAR) == y, ]
  visit_occ_site <- unique(dataset_y[!is.na(COUNT) & ANCHOR == 0L, visitN := .N, by = SITE_ID][
                                  !is.na(COUNT) & ANCHOR == 0L & COUNT > 0, occurN := .N, by=SITE_ID][
                                  !is.na(COUNT) & ANCHOR == 0L, ][order(SITE_ID), .(SITE_ID, occurN, visitN)])

  dataset_y <- data.table::copy(dataset_y[SITE_ID %in% 
                                  visit_occ_site[visitN >= MinVisit & occurN >= MinOccur, SITE_ID], ][,
                                  visitN := NULL][, occurN := NULL])
  if(is.null(tp_col)){
    if(TimeUnit == 'd'){
      tp_col <- "trimDAYNO"
    } else {
      tp_col <- "trimWEEKNO"
    }
  }else{
    tp_col <- tp_col
  }

  if(dataset_y[, uniqueN(SITE_ID)] < MinNbrSite){

    f_curve  <- unique(ts_season_count[as.integer(M_YEAR) == y, ][ , SITE_ID := NULL][, COUNT := NULL])[, NM := NA]
    f_curve_mod <- list(f_curve=f_curve[order(tp_col),], f_model=list(NA), f_data=data.table(NA))

    if(verbose){
      message(paste("You have not enough sites with observations for estimating the flight curve for species", 
                    as.character(ts_season_count$SPECIES[1]), "in", unique(ts_season_count[as.integer(M_YEAR) == y, M_YEAR])))
    }  
  } else {
    f_curve_mod <- fit_gam(dataset_y, NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial,
                          SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit, MultiVisit = MultiVisit, weekday = weekday,
                          mod_form = mod_form, tp_col = tp_col, verbose = verbose, ...)
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
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param mod_form string with formula to be passed to the gam model, default null.
#' @param tp_col string or vector of string with additional variable used in the gam model, default null.
#' @param verbose a logical indicating if some “progress report” should be given. 
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **pheno**: a vector with annual flight curves \code{f_pheno} with expected relative abundance, normalize to sum to one over a full season,
#'         ii) **model**: a list of the resulting gam models \code{f_model} fitted on the count data for each year and iii) **data**: a data.table with the data used to fit the GAM model.
#' @keywords gam, flight curve
#' @seealso \link{fit_gam}
#' @seealso \link[mgcv]{gam}
#' @seealso \link[mgcv]{bam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export flight_curve
#'

flight_curve <- function(ts_season_count, NbrSample = 100, MinVisit = 3, MinOccur = 2, MinNbrSite = 1, MaxTrial = 3,
                         GamFamily = 'poisson', CompltSeason = TRUE, SelectYear = NULL, SpeedGam = TRUE,
                         OptiGam = TRUE, KeepModel = TRUE, KeepModelData = TRUE, TimeUnit = 'd', 
                         MultiVisit = "mean", weekday = 3, mod_form = NULL, tp_col = NULL, verbose = TRUE, ...) {


        check_package('data.table')

        names(ts_season_count) <- toupper(names(ts_season_count))
        if(!is.null(tp_col)){ 
          tp_col <- toupper(tp_col)
        }

        col_name <- c("COMPLT_SEASON", "M_YEAR", "SITE_ID", "SPECIES", "DATE", "WEEK", "WEEK_DAY", "DAY_SINCE", "WEEK_SINCE", "M_SEASON", "COUNT", "ANCHOR")
        check_names(ts_season_count, ifelse(is.null(tp_col), col_name, c(col_name, tp_col)))

        if(isTRUE(CompltSeason)){
            ts_season_count <- ts_season_count[COMPLT_SEASON == 1, ]
        }

        if(is.null(SelectYear)){
            year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, unique(as.integer(M_YEAR))]
        } else {
            year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, ][M_YEAR %in% SelectYear, unique(as.integer(M_YEAR))]
        }

        if(length(year_series) == 0) stop(paste0(" No count data found for year ", SelectYear))

        if(!is.null(tp_col)){
          if(!all(tp_col %in% col_name)){
          if(TimeUnit == 'd'){
            ts_season_count_tp_col <- unique(ts_season_count[ , c(c("SPECIES", "SITE_ID", "DAY_SINCE"), tp_col), with = FALSE])
            setkey(ts_season_count_tp_col, SPECIES, SITE_ID, DAY_SINCE)
          }else{
            ## compute mean weekly value for tp_col (temporal variable column) if TimeUnit is week
            ts_season_count_m <- ts_season_count[, lapply(.SD, mean, na.rm = TRUE), by = .(SPECIES, SITE_ID, WEEK_SINCE), .SDcols = tp_col]
            ts_season_count_tp_col <- unique(ts_season_count_m[ , c(c("SPECIES", "SITE_ID", "WEEK_SINCE"), tp_col), with = FALSE])
            setkey(ts_season_count_tp_col, SPECIES, SITE_ID, WEEK_SINCE)
          }
          }
        }

        ts_season_count <- day_week_summary(ts_season_count, MultiVisit = MultiVisit, TimeUnit = TimeUnit, weekday = weekday)

        if(!is.null(tp_col)){
          if(!all(tp_col %in% col_name)){
            ifelse(TimeUnit == 'd', setkey(ts_season_count, SPECIES, SITE_ID, DAY_SINCE), setkey(ts_season_count, SPECIES, SITE_ID, WEEK_SINCE))
            ts_season_count <- merge(ts_season_count, ts_season_count_tp_col, all.x = TRUE)
          }
        }

        result_fc <- lapply(year_series, get_nm, ts_season_count = ts_season_count, MinVisit = MinVisit, MinOccur = MinOccur, MinNbrSite = MinNbrSite,
                                                  NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial, SpeedGam = SpeedGam,
                                                  OptiGam = OptiGam, TimeUnit = TimeUnit, MultiVisit = MultiVisit, weekday = weekday, mod_form = mod_form,
                                                  tp_col = tp_col, verbose = verbose, ...) 

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
    class(result_fc) <- "pheno_curve"
    return(result_fc)
}


#' get_nny
#' find the nearest year with a computed flight curve
#' @param x integer Year for which to find the nearest flight curve.
#' @param y vector Years with available flight curve.
#' @return z integer The value of the nearest year with flight curve.
#' @keywords flight curve
#' @seealso \link{check_pheno}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
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
#' @param verbose a logical indicating if some “progress report” should be given. 
#' @return A data.table with time-series of the expected relative abundance of butterfly count per day (NM) for the year or the nearest year where
#'         phenology is available.
#' @keywords flight curve
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export check_pheno
#'

check_pheno <- function(sp_count_flight, ts_flight_curve, YearCheck, YearLimit, TimeUnit, verbose = TRUE){

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
            if(verbose){
              message(paste("No reliable flight curve available within a ",YearLimit," year horizon of", sp_count_flight_y[1, M_YEAR, ]))
            }
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
#' @param ts_flight_curve pheno_curve class object or data.table with light curves and relative abundances (NM) for a specific species as returned by \link{flight_curve}.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week.
#' @param MultiVisit string Function to apply for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param sp integer or string Species ID or name.
#' @param YearLimit integer Define the span of years (+/- number of year) to look for a flight curve, if NULL no restriction is set.
#' @param SelectYear integer Select a specific year to compute the flight curve, default=NULL.
#' @param CompltSeason logical Restrict computation of flight curve for years where the complete season was sampled, default=TRUE.
#' @param SindexMethod character Method to compute the site index, 'ratio' (default) or 'glm', see details.
#' @param GlmFit character Fitting of the GLM when SindexMethod = 'glm', 'newton' (default) or 'ipf', see details.
#' @param verbose a logical indicating if some “progress report” should be given.
#' @return A data.table based on the entry count data, augmented with site indices 'SINDEX' and imputed weekly count 'IMPUTED_COUNT'.
#' @details Site indices can be extracted from the data.table returned by this function. With SindexMethod = 'ratio', the site index
#' is computed by adjusting the count by the proportion of the flight curve covered by the visits (TOTAL_COUNT / TOTAL_NM), each site-year
#' being estimated on its own. With SindexMethod = 'glm', the counts of the time-steps not visited are predicted from a Poisson GLM with site
#' and year effects and the flight curve as offset, COUNT ~ factor(SITE_ID) + factor(M_YEAR) + offset(log(NM)), and the site index is the sum
#' of the observed and predicted counts. The GLM borrows the level of the site and the year from the other years and sites, which avoids extreme
#' site indices when a small part of the flight curve is covered, but assumes no site x year interaction for the time-steps not visited.
#' Sites where the species was never observed keep the ratio (0).
#' The GLM is fitted without model matrix, on the counts summed per site-year. With GlmFit = 'newton', the site effects are
#' concentrated out of the likelihood and Newton steps are made on the year effects only (few iterations); with
#' GlmFit = 'ipf', iterative proportional fitting is used (more iterations when years are weakly linked). Both give the
#' estimates of glm(), 'newton' uses 'ipf' if it fails. With 'newton', the year effects of the GLM are returned in the
#' attribute "glm_year_effect", \code{attr(x, "glm_year_effect")}: log effect relative to the first year of each group of
#' linked years (LOG_EFFECT), with a quasi-Poisson standard error conditional on the site totals (SE). This is a model-based
#' annual trend from the counts, a quick check of the collated index, but its SE ignores the variation among sites captured
#' by the bootstrap.
#' @keywords site index, flight curve
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export impute_count
#'

impute_count <- function(ts_season_count, ts_flight_curve, TimeUnit = 'd', MultiVisit = 'mean', weekday = 3,
                         sp = NULL, YearLimit= NULL, SelectYear = NULL, CompltSeason = TRUE, SindexMethod = c("ratio", "glm"),
                         GlmFit = c("newton", "ipf"), verbose = TRUE){

        check_package('data.table')
        SindexMethod <- match.arg(SindexMethod)
        GlmFit <- match.arg(GlmFit)

        if(!is.null(SelectYear)){
          ts_season_count <- ts_season_count[M_YEAR == SelectYear, ]
        }

        if(!is.null(sp)){
          ts_season_count <- ts_season_count[SPECIES == sp, ]
        }

        if(inherits(ts_flight_curve, "pheno_curve")){
          ts_flight_curve <- ts_flight_curve$pheno
        }

        if(!ts_season_count$SPECIES[1] %in% unique(ts_flight_curve$SPECIES)){
            stop('Species in count data must be in the flight curve data!')
        }

        if(TimeUnit == 'd'){
            tp_col <- "trimDAYNO"
        } else {
            tp_col <- "trimWEEKNO"
        }

        ts_season_count <- day_week_summary(ts_season_count, MultiVisit = MultiVisit, TimeUnit = TimeUnit, weekday = weekday)

        ts_season_count$M_YEAR <- as.integer(as.character(ts_season_count$M_YEAR))
        ts_flight_curve$M_YEAR <- as.integer(as.character(ts_flight_curve$M_YEAR))

        keycol <- c("SPECIES", "M_YEAR", tp_col)
        data.table::setkeyv(ts_season_count, keycol)
        data.table::setkeyv(ts_flight_curve, keycol)
        mcol <- c(keycol, "NM")

        sp_count_flight <- merge(ts_season_count, ts_flight_curve[!duplicated(ts_flight_curve[, mcol, with = FALSE]), mcol, with = FALSE], all.x=TRUE)

        YearCheck <- as.integer(as.character(unique(sp_count_flight[ is.na(NM), M_YEAR])))

        if(length(YearCheck) > 0){
          a <- lapply(YearCheck, check_pheno, sp_count_flight=sp_count_flight, ts_flight_curve=ts_flight_curve, YearLimit=YearLimit, TimeUnit = TimeUnit,
                      verbose = verbose)
          a_rb <- data.table::rbindlist(a, use.names = TRUE, fill = TRUE)
          sp_count_flight <- rbind(a_rb, sp_count_flight[!M_YEAR %in% a_rb[!is.na(NM), .N, by= M_YEAR][N>0, M_YEAR], ], use.names = TRUE, fill=TRUE)
        }

        total_count <- sp_count_flight[M_SEASON != 0 , sum(COUNT, na.rm = TRUE), by = .(M_YEAR, SITE_ID)][, TOTAL_COUNT := V1][, V1 := NULL]
        total_nm <- sp_count_flight[M_SEASON != 0 & !is.na(COUNT) , sum(NM, na.rm = TRUE), by = .(M_YEAR, SITE_ID)][, TOTAL_NM := V1][, V1 := NULL]
        setkey(total_count, M_YEAR, SITE_ID)
        setkey(total_nm, M_YEAR, SITE_ID)
        setkey(sp_count_flight, M_YEAR, SITE_ID)

        sp_count_flight <- merge(sp_count_flight, merge(total_count, total_nm, all.x=TRUE)[, SINDEX := TOTAL_COUNT / TOTAL_NM], all.x = TRUE)
        sp_count_flight[, IMPUTED_COUNT:= COUNT][M_SEASON != 0 & is.na(COUNT), IMPUTED_COUNT := as.integer(round(SINDEX*NM))]

        glm_year_effect <- NULL
        if(SindexMethod == "glm"){
          ## site-year level from the Poisson GLM COUNT ~ site + year + offset(log(NM)), fitted on the counts of the
          ## sites where the species was observed; other site-years keep the ratio (e.g. 0 for sites never observed)
          obs <- sp_count_flight[M_SEASON != 0 & !is.na(COUNT) & !is.na(NM) & NM > 0, .(M_YEAR, SITE_ID, COUNT, NM)]
          obs <- obs[SITE_ID %in% obs[, sum(COUNT), by = SITE_ID][V1 > 0, SITE_ID], ]

          if(nrow(obs) > 0){
            glm_level <- site_year_level(obs$SITE_ID, obs$M_YEAR, obs$COUNT, obs$NM, method = GlmFit)
            glm_year_effect <- attr(glm_level, "year_effect")
            obs[, GLM_LEVEL := as.numeric(glm_level)]
            sp_count_flight <- merge(sp_count_flight, unique(obs[, .(M_YEAR, SITE_ID, GLM_LEVEL)]), by = c("M_YEAR", "SITE_ID"), all.x = TRUE)
            sp_count_flight[!is.na(GLM_LEVEL) & M_SEASON != 0 & is.na(COUNT), IMPUTED_COUNT := as.integer(round(GLM_LEVEL*NM))]
            ## SINDEX: observed counts plus counts predicted for the time-steps not visited
            glm_sindex <- sp_count_flight[!is.na(GLM_LEVEL) & M_SEASON != 0,
                                          .(SINDEX_GLM = sum(data.table::fifelse(is.na(COUNT), GLM_LEVEL*NM, as.numeric(COUNT)), na.rm = TRUE)),
                                          by = .(M_YEAR, SITE_ID)]
            sp_count_flight <- merge(sp_count_flight, glm_sindex, by = c("M_YEAR", "SITE_ID"), all.x = TRUE)
            sp_count_flight[!is.na(SINDEX_GLM), SINDEX := SINDEX_GLM][, c("GLM_LEVEL", "SINDEX_GLM") := NULL]
          }
        }

        if(TimeUnit == 'd') {
          sp_count_flight[ , SINDEX := SINDEX / 7]
          week_count <- unique(sp_count_flight[M_SEASON != 0 & !is.na(COUNT), .(M_YEAR, SITE_ID, WEEK)])[ , WEEK_COUNT := 1]
          setkey(week_count, SITE_ID, M_YEAR, WEEK)
          setkey(sp_count_flight, SITE_ID, M_YEAR, WEEK)
          sp_count_flight <- merge(sp_count_flight, week_count, all.x = TRUE)
          sp_count_flight[ , WEEK_NM := mean(NM * WEEK_COUNT, na.rm = TRUE), by = .(M_YEAR, SITE_ID, WEEK)][ ,
                                             TOTAL_NM := sum(WEEK_NM, na.rm=TRUE), by = .(M_YEAR, SITE_ID)][ , WEEK_NM := NULL][, WEEK_COUNT := NULL]
        }

        if(!is.null(glm_year_effect)) data.table::setattr(sp_count_flight, "glm_year_effect", glm_year_effect)

     return(sp_count_flight)
  }


#' site_index
#' Extract abundance indices per site and year based on flight curve imputation.
#' @param butterfly_count data.table Observed and imputed weekly or daily counts and the estimated total counts (SINDEX) computed by the function impute_count().
#' @param MinFC numeric Value between 0 and 1 to define the threshold for the proportion of the flight curve covered by the visits, if NULL all site-year available indices are returned.
#' @return data.table Estimated annual abundance index and the proportion of the flight curve covered by the visit - total weekly or daily count over the entire monitoring season.
#' @keywords butterfly count
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
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
#' Deprecated, use \link{collated_index}. The function is removed and only returns an error; it will be deleted in a
#' future version.
#' @param ... not used.
#' @keywords internal
#' @seealso \link{collated_index}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @export collated_index_old
#'

collated_index_old <- function(...){

  .Deprecated("collated_index", package = "rbms")
  stop("collated_index_old() has been removed, use collated_index()", call. = FALSE)
}

#' collated_index
#' compute a collated index from the site indices, using a Generalized Linear Model.
#' @details Site-years with NA in \code{sindex_value} (e.g. years where the flight curve failed) are treated as missing
#' and removed with a warning, as done by \link{site_index}; a site with NA in all years is therefore not part of the sample.
#' When \code{boot_ind} was generated for a species with \code{boot_sample(s_sp = )}, \code{s_sp} must be the same species.
#' A bootstrap sample that failed in \code{boot_sample()} (NA row) returns NA for the collated index and the number of sites.
#' With GlmFit = 'newton' or 'ipf', the same weighted Poisson GLM (site and year effects) is fitted without building the model
#' matrix (see \link{impute_count} for the two methods), which is much faster with many sites and gives the same collated index as
#' \code{glm()}. With 'newton', year_effect gives per year the log year effect relative to the first year of each group of linked
#' years (LOG_EFFECT), i.e. the log of the collated index relative to that year, with its quasi-Poisson standard error (SE), an
#' alternative to the bootstrap that ignores the variation among sites. The column GROUP gives the group of years linked through
#' sites with observations in the original data (see \link{n_components}), NA for years without observation. When the years form
#' more than one group, the collated index is only comparable between years of the same group (warning): with SplitGroup = TRUE,
#' each group is fitted separately on its own years and sites, its index being the mean over the sites with data in the years of
#' the group, so each group is a separate series. Fitting all years together (SplitGroup = FALSE, previous versions) gives arbitrary
#' levels between groups and can distort the trend within a group, because zero counts of sites observed only in another group
#' enter the fit.
#' @param data data.table or data.frame Site indices per year and proportion of flight curve covered by the monitoring.
#' @param s_sp string Species name to be found in the data.
#' @param sindex_value string Name of the response variable to be used the computation of the collated index, default is "SINDEX", but could be other standardized values.
#' @param bootID integer Identify the n^th bootstrap.
#' @param boot_ind data.table Index of the bootstrap and the site ids of the specific bootstrap sample.
#' @param glm_weights logical Use the proportion of the flight curve sampled as weight for the collated index, if FALSE the function uses equal weights.
#' @param rm_zero logical Remove the sites where species was not observed to speed-up the fit of the GLM without altering the output, default TRUE.
#' @param GlmFit character Fitting of the GLM, 'glm' (default, \code{glm()}), 'newton' or 'ipf', see details.
#' @param SplitGroup logical When the years form groups not linked through sites with observations, fit each group separately
#' (default TRUE, the correct method). FALSE fits all years together as in previous versions (glm only), the levels of the
#' groups are then arbitrary; 'newton' and 'ipf' always fit the groups separately.
#' @return a list with the collated indices per year (col_index), the site ids of the sample (site_id) and, with GlmFit = 'newton',
#' the year effects (year_effect, see details), NULL otherwise.
#' @keywords annual index
#' @seealso \link{impute_count}, \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export collated_index
#'

collated_index <- function(data, s_sp, sindex_value = "SINDEX", bootID=NULL,
    boot_ind=NULL, glm_weights=TRUE, rm_zero=TRUE, GlmFit = c("glm", "newton", "ipf"), SplitGroup = TRUE){

  GlmFit <- match.arg(GlmFit)
  year_effect <- NULL
  data <- setDT(data)
  data = data[SPECIES == s_sp, ]

  ## site-years without index (e.g. failed flight curve) are missing, not zero: removed as in site_index()
  n_na <- data[is.na(get(sindex_value)), .N]
  if(n_na > 0){
    warning(paste0(s_sp, ": ", n_na, " site-years with NA in ", sindex_value, " removed before computing the collated index"), call. = FALSE)
    data <- data[!is.na(get(sindex_value)), ]
  }
  if(nrow(data) == 0) stop(paste0("no ", sindex_value, " available for species ", s_sp), call. = FALSE)

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
    bootID <- bootID
  }

  ## bootstrap samples drawn per species with boot_sample(s_sp = ) must match the species modelled
  if(!is.null(boot_ind$s_sp) && !identical(as.character(boot_ind$s_sp), as.character(s_sp))){
    stop(paste0("boot_ind was generated for species ", boot_ind$s_sp, ", not for ", s_sp), call. = FALSE)
  }

  ## bootstrap sample that failed in boot_sample() after max_trial trials (NA row): no index computed
  if(bootID > 0 && anyNA(boot_ind$boot_ind[bootID, ])){
    res <- data.table(BOOTi = bootID, M_YEAR = y[order(y)], NSITE = NA_integer_, NSITE_B = NA_integer_,
                      NSITE_OBS = NA_integer_, COL_INDEX = NA_real_, GROUP = NA_integer_)
    setkey(res, M_YEAR)
    return(list(col_index = res, site_id = NA_character_, year_effect = NULL))
  }

  if(bootID == 0){
    if(is.null(boot_ind)){
    a <- data.table(uID = as.character(paste0("s_",seq_len(data[, uniqueN(SITE_ID)]))),
     SITE_ID = as.character(unique(data[, SITE_ID])))
    }else{
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID+1, ]))),
     SITE_ID = as.character(boot_ind$site_id[seq_along(boot_ind$boot_ind[bootID+1, ])]))
    }
  } else {
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID, ]))),
     SITE_ID = as.character(boot_ind$site_id[boot_ind$boot_ind[bootID, ]]))
  }

  s_y <- data.table(expand.grid(uID = a$uID, M_YEAR = y[order(y)]))
  # names(s_y) <- c("uID", "M_YEAR")
  s_y[ , uID:=as.character(uID)]
  s_y[ , M_YEAR:=as.numeric(as.character(M_YEAR))]
  data[ , SITE_ID:=as.character(SITE_ID)]
  data[ , M_YEAR:=as.numeric(as.character(M_YEAR))]

  ## groups of years linked through sites with observations in the original data (see n_components),
  ## years without observation have no group
  obs_0 <- data[get(sindex_value) > 0, ]
  year_group <- if(nrow(obs_0) > 0) year_components(obs_0$SITE_ID, obs_0$M_YEAR) else integer(0)
  grp <- data.table(M_YEAR = as.numeric(names(year_group)), GROUP = as.integer(year_group))
  if(bootID == 0 && uniqueN(grp$GROUP) > 1){
    warning(paste0(s_sp, ": the years form ", uniqueN(grp$GROUP), " groups not linked through sites with observations (",
                   paste(grp[, paste(range(M_YEAR), collapse = "-"), by = GROUP]$V1, collapse = ", "),
                   "); the collated index is only comparable between years of the same group, ",
                   if(isTRUE(SplitGroup) || GlmFit != "glm") "each group is fitted separately (SplitGroup = TRUE)"
                   else "all years are fitted together (SplitGroup = FALSE) and the levels of the groups are arbitrary, use SplitGroup = TRUE"),
            call. = FALSE)
  }

  setkey(a, uID);
  setkey(s_y, uID)
  s_y <- merge(s_y, a, all.x=TRUE)

  setkey(s_y, SITE_ID, M_YEAR); setkey(data, SITE_ID, M_YEAR)

  if(sindex_value != 'SINDEX'){
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)][, (sindex_value) := data[, get(sindex_value)]])
    } else {
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)])
    }

  sum_data <- data.table( BOOTi = bootID,
                          M_YEAR = y[order(y)],
                          NSITE = data_boot[, uniqueN(uID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_B = data_boot[, uniqueN(SITE_ID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_OBS = data_boot[, sum(get(sindex_value)>0), by = M_YEAR][order(M_YEAR), V1])

  setkey(sum_data, M_YEAR)

  if(nrow(sum_data[NSITE_OBS > 0,]) != 0){

    if(isTRUE(rm_zero)){
      zero_site <- data_boot[ , tmpV1:=sum(SINDEX>0), uID][tmpV1==0, uID]
      data_boot <- data_boot[!uID %in% zero_site, -c("tmpV1")]
    }

    if (data_boot[, uniqueN(M_YEAR)] == 1 & data_boot[, uniqueN(SITE_ID)] == 1) warning(paste0(sindex_value, " available for only one site and one year: no collated index computed for ", s_sp, " - ", unique(data_boot$M_YEAR)), call. = FALSE)

    if(!isTRUE(glm_weights)){
      data_boot[ , weights := 1]
    }else{
      ## make sure that column TOTAL_NM exists
      data_boot[ , weights := TOTAL_NM]
    }

    # col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(),
    #                   weights = data_boot$weights), silent = TRUE)
    #   if (inherits(col_index, "try-error")) {
    #     col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(),
    #                   weights = data_boot$weights, method = "qr"), silent = TRUE)
    #     if (inherits(col_index, "try-error")) {
    #       col_index <- try(glm(mod_form, data = data_boot, family = poisson(), 
    #                   weights = data_boot$weights), silent = TRUE)
    #         if (inherits(col_index, "try-error")) {
    #           res <- sum_data[, COL_INDEX := NA]
    #           return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_OBS, COL_INDEX)], site_id = a$SITE_ID ))
    #         }
    #     }
    #  }

# patch for glmmTMB issue
#    col_index <- try(glmmTMB::glmmTMB(mod_form, data = data_boot, family = "nbinom1",
#                      weights = data_boot$weights), silent = TRUE)
#      if (inherits(col_index, "try-error")) {
#        col_index <- try(glmmTMB::glmmTMB(mod_form, data = data_boot, family = poisson(),
#                      weights = data_boot$weights), silent = TRUE)
#        if (inherits(col_index, "try-error")) {
    ## years fitted together: each group of linked years separately (SplitGroup = TRUE, always with newton and ipf),
    ## or all years at once (one group, or SplitGroup = FALSE with glm, the former joint fit). Index of a group = prediction
    ## over its sites and years / number of sites with data in its years (all sites of the species when fitted at once)
    if(uniqueN(grp$GROUP) > 1 && (isTRUE(SplitGroup) || GlmFit != "glm")){
      fit_years <- split(grp$M_YEAR, grp$GROUP)
    } else {
      fit_years <- list(`1` = unique(data_boot$M_YEAR))
    }

    co_index_res <- list()
    ye <- list()
    for(g in names(fit_years)){
      yrs <- fit_years[[g]]
      d_g <- data_boot[M_YEAR %in% yrs, ]
      if(nrow(d_g) == 0 || !any(d_g[[sindex_value]] > 0)) next
      n_site_g <- if(length(fit_years) == 1) length(unique(data$SITE_ID)) else data[M_YEAR %in% yrs, uniqueN(SITE_ID)]

      if(GlmFit == "glm"){
        if (d_g[, uniqueN(M_YEAR)] == 1) {
          mod_form <- as.formula(paste0("I(round(", sindex_value, ")) ~ factor(uID) - 1"))
        } else {
          mod_form <- as.formula(ifelse(d_g[, uniqueN(SITE_ID)] > 1, paste0("I(round(", sindex_value, ")) ~ factor(uID) + factor(M_YEAR) - 1"), paste0("I(round(", sindex_value, ")) ~ factor(M_YEAR) - 1")))
        }
          col_index <- try(glm(mod_form, data = d_g, family = poisson(),
                      weights = d_g$weights), silent = TRUE)
            if (inherits(col_index, "try-error")) {
              res <- merge(sum_data[, COL_INDEX := NA_real_], grp, by = "M_YEAR", all.x = TRUE)
              return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_B, NSITE_OBS, COL_INDEX, GROUP)], site_id = a$SITE_ID, year_effect = NULL))
            }
#        }
#     }
        pred_data_boot <- data.table(expand.grid(uID = unique(d_g$uID), M_YEAR = unique(d_g$M_YEAR)))
        pred_data_boot$fitted <- predict(col_index, newdata = pred_data_boot, type = "response")
        co_g <- pred_data_boot[order(M_YEAR), .(COL_INDEX = sum(fitted)/n_site_g), by = M_YEAR]
      } else {
        ## same weighted Poisson GLM fitted without model matrix (weighted counts, weights as flight curve covered):
        ## prediction over all sites and years = year effect x sum of site effects
        fit_level <- site_year_level(d_g$uID, d_g$M_YEAR, d_g$weights * round(d_g[[sindex_value]]), d_g$weights, method = GlmFit)
        eff <- attr(fit_level, "effects")
        co_g <- data.table(M_YEAR = as.numeric(names(eff$year)), COL_INDEX = eff$year * sum(eff$site) / n_site_g)
        if(!is.null(attr(fit_level, "year_effect"))){
          ye[[length(ye) + 1]] <- attr(fit_level, "year_effect")[, M_YEAR := as.numeric(M_YEAR)][, GROUP := NULL]
        }
      }
      co_index_res[[length(co_index_res) + 1]] <- co_g
    }
    co_index_res <- if(length(co_index_res) > 0) rbindlist(co_index_res) else data.table(M_YEAR = numeric(0), COL_INDEX = numeric(0))
    if(length(ye) > 0) year_effect <- merge(rbindlist(ye), grp, by = "M_YEAR", all.x = TRUE)[, .(M_YEAR, GROUP, LOG_EFFECT, SE, DISPERSION)]
    setkey(co_index_res, M_YEAR)
    res <- merge(sum_data, co_index_res, all.x = TRUE)[is.na(COL_INDEX), COL_INDEX := 0]

  } else {
    res <- sum_data[, COL_INDEX := 0]
  }

  res <- merge(res, grp, by = "M_YEAR", all.x = TRUE)
  setkey(res, M_YEAR)

  return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_B, NSITE_OBS, COL_INDEX, GROUP)],
                site_id = a$SITE_ID, year_effect = year_effect))
}


#' boot_sample
#' Generate n bootstrap sample of the monitoring sites to be used for each iteration.
#' @details When \code{s_sp} is given, the samples are drawn for that species only, stratified by the occurrence of zeros
#' at the sites (never zero, some zeros, always zero), and a sample is drawn again, up to \code{max_trial} trials,
#' when it has fewer than \code{min_year} years or when its years are not linked through the sites with observations
#' (see \code{n_components}). A sample still failing after \code{max_trial} trials is returned as a row of NA.
#' When the years of the original data form several groups not linked through sites with observations (see
#' \link{collated_index}), the strata are crossed with the groups (a site with observations belongs to one group, a site
#' never observed to the group of most of its years), so each sample keeps the number of sites of each group, and each group
#' must remain linked in the sample (not split or lost). With one group, this is the stratification by zeros only.
#' For reproducible samples, call \code{set.seed()} before \code{boot_sample()}; when looping over species, setting
#' the seed before each species makes its samples independent of the order or the subset of species processed.
#' @param data data.table or data.frame Data with all site id, with SPECIES, M_YEAR and \code{sindex_value} when \code{s_sp} is given.
#' @param boot_n integer The number of bootstrap samples to be generated.
#' @param s_sp string Species for which the bootstrap samples are drawn, default NULL draws among all sites without checks.
#' @param sindex_value string Name of the site index variable used to define strata and observations, default is "SINDEX".
#' @param strata logical Stratify the sampling by the occurrence of zeros at the sites and by group of linked years, default TRUE.
#' @param check_link logical Draw again the samples where a group of linked years is split or lost, default TRUE.
#' @param min_year integer Minimum number of years in a bootstrap sample, default 3.
#' @param max_trial integer Maximum number of trials to draw a valid bootstrap sample, default 5.
#' @return A list with site id and bootstrap indices for n bootstrap sample. When \code{s_sp} is given, the list also
#' contains the species (s_sp), the number of groups of years of the original data (n_comp_0), the site indices per
#' stratum (strata_ind, named group_zero pattern), the group of each year (year_group) and of each site (site_group),
#' the number of trials per sample (n_trial) and the number of failed samples (n_failed).
#' @keywords bootstrap collated index
#' @seealso \link{impute_count}, \link{collated_index}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export boot_sample
#'

boot_sample <- function(data, boot_n = 1000, s_sp = NULL, sindex_value = "SINDEX",
    strata = TRUE, check_link = TRUE, min_year = 3, max_trial = 5){

  data <- data.table::as.data.table(data)

  if(is.null(s_sp)){
    site_id <- unique(data[, SITE_ID])
    boot_ind <- matrix(NA, nrow = boot_n, ncol = length(site_id))

      for(i in seq_len(boot_n)){
        boot_ind[i, ] <- sample(seq_len(length(site_id)), replace = TRUE)
      }

    boot_ind <- list(site_id = site_id, boot_ind = boot_ind)
    return(boot_ind)
  }

  dt <- data[SPECIES == s_sp, .(SITE_ID = as.character(SITE_ID), M_YEAR = as.numeric(as.character(M_YEAR)), SVALUE = get(sindex_value))]
  dt <- dt[!is.na(SVALUE), ]
  if(nrow(dt) == 0) stop(paste0("no ", sindex_value, " available for species ", s_sp), call. = FALSE)

  site_id <- unique(dt$SITE_ID)
  obs <- dt[SVALUE > 0, ]
  n_comp_0 <- n_components(obs$SITE_ID, obs$M_YEAR)

  ## groups of years linked through sites with observations (see n_components); a site with observations belongs to
  ## one group, a site never observed to the group of most of its years (years without observation have no group)
  year_group <- if(nrow(obs) > 0) year_components(obs$SITE_ID, obs$M_YEAR) else integer(0)
  year_group <- data.table(M_YEAR = as.numeric(names(year_group)), GROUP = as.integer(year_group))
  if(n_comp_0 > 1 && is.finite(n_comp_0)){
    warning(paste0(s_sp, ": the years of the original data form ", n_comp_0, " groups not linked through sites with observations (",
                   paste(year_group[, paste(range(M_YEAR), collapse = "-"), by = GROUP]$V1, collapse = ", "),
                   "); the samples are stratified by group and each group must remain linked in each sample."), call. = FALSE)
  }
  site_group <- merge(dt, year_group, by = "M_YEAR", all.x = TRUE)[, .(GROUP = {
    g <- if(any(SVALUE > 0)) GROUP[SVALUE > 0] else GROUP[!is.na(GROUP)]
    if(length(g) == 0) NA_integer_ else as.integer(names(which.max(table(g))))}), by = SITE_ID]
  if(uniqueN(year_group$GROUP) <= 1) site_group[, GROUP := 1L]

  ## strata: group x zero pattern, 1 = never zero, 2 = some zeros, 3 = always zero (keeps the number of sites of each
  ## group and the share of zero sites in each sample); with one group, the zero pattern only
  st <- dt[, .(ZERO = if(all(SVALUE > 0)) 1L else if(any(SVALUE > 0)) 2L else 3L), by = SITE_ID]
  st <- merge(st, site_group, by = "SITE_ID", sort = FALSE)
  st[, STRATA := if(isTRUE(strata)) paste(GROUP, ZERO, sep = "_") else "1"]
  strata_ind <- split(match(st$SITE_ID, site_id), st$STRATA)

  boot_ind <- matrix(NA_integer_, nrow = boot_n, ncol = length(site_id))
  n_trial <- integer(boot_n)

    for(i in seq_len(boot_n)){
      for(tr in seq_len(max_trial)){
        ind <- draw_strata(strata_ind)
        s_i <- site_id[unique(ind)] # duplicated sites do not change the years or their linkage
        ok <- dt[SITE_ID %in% s_i, uniqueN(M_YEAR)] >= min_year
        if(ok && isTRUE(check_link)){
          ## each group of the original data must remain linked (one group) in the sample, not split or lost
          o <- obs[SITE_ID %in% s_i, ]
          for(g in unique(year_group$GROUP)){
            o_g <- o[M_YEAR %in% year_group[GROUP == g, M_YEAR], ]
            if(n_components(o_g$SITE_ID, o_g$M_YEAR) != 1){ ok <- FALSE; break }
          }
        }
        if(ok) break
      }
      n_trial[i] <- tr
      if(ok) boot_ind[i, ] <- ind
    }

  n_failed <- sum(is.na(boot_ind[, 1]))
  if(n_failed > 0) warning(paste0(s_sp, ": ", n_failed, " of ", boot_n, " bootstrap samples failed after ", max_trial, " trials (NA rows)"), call. = FALSE)

  boot_ind <- list(site_id = site_id, boot_ind = boot_ind, s_sp = s_sp, n_comp_0 = n_comp_0,
                   strata_ind = strata_ind, year_group = year_group, site_group = site_group[match(site_id, SITE_ID), ],
                   n_trial = n_trial, n_failed = n_failed)
  return(boot_ind)
}


#' day_week_summary
#' Summarize the count(s) per day or week using either the maximum or the average count when multiple counts (visit) are observed.
#' @param ts_season_count data.table Time-series of butterfly counts for species x over year y over all sites.
#' @param TimeUnit character The time-step for which the spline should be computed, 'd' day or 'w' week. For TimeUnit week, 'w', we refer to week_day 3 to match a mid-week calendar date.
#' @param MultiVisit string Function to apply for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @return data.table Time-series of butterfly counts for species x over year y over all sites, but summarized with the chosen function.
#' @keywords flight_curve
#' @seealso \link{flight_curve}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export day_week_summary
#'

day_week_summary <- function(ts_season_count, MultiVisit, TimeUnit, weekday){

  if(MultiVisit == "mean"){
            if(TimeUnit == 'd'){
              ts_season_count[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1, by = M_YEAR]
              ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                      WEEK_DAY, WEEK_SINCE, DAY_SINCE, trimDAYNO, 
                                                      M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]

              ts_season_count[ , meanCOUNT := ceiling(mean(COUNT, na.rm = TRUE)),
                                                    by = .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)][
                                                        is.nan(meanCOUNT), meanCOUNT := NA]

              ts_season_count_summary <- unique(ts_season_count[, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                                    WEEK_DAY, WEEK_SINCE, DAY_SINCE, trimDAYNO, 
                                                                    M_SEASON, COMPLT_SEASON, ANCHOR, meanCOUNT)])[, 
                                                                    COUNT := meanCOUNT][, meanCOUNT := NULL]
            } else {
              ts_season_count[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
              ts_season_count <- ts_season_count[, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                        WEEK_DAY, WEEK_SINCE, trimWEEKNO, 
                                                        M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]

              ts_season_count[ , meanCOUNT := ceiling(mean(COUNT, na.rm = TRUE)),
                                                      by = .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)][
                                                        is.nan(meanCOUNT), meanCOUNT := NA]

              ts_season_count_summary <- unique(ts_season_count[WEEK_DAY == weekday, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                                            WEEK_DAY, WEEK_SINCE, trimWEEKNO, 
                                                                            M_SEASON, COMPLT_SEASON, ANCHOR, meanCOUNT)])[, 
                                                                            COUNT := meanCOUNT][, meanCOUNT := NULL]
              # make trimWEEKNO start with 1, it might start with 2 if the selected weekday was not in the first week 
              ts_season_count_summary[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
            }
  } else {
            if(TimeUnit == 'd'){

              ts_season_count[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1, by = M_YEAR]
              ts_season_count <- ts_season_count[ , .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                      WEEK_DAY, WEEK_SINCE, DAY_SINCE, trimDAYNO, 
                                                      M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]

              ts_season_count[ , maxCOUNT := suppressWarnings(max(COUNT, na.rm = TRUE)),
                                                    by = .(SPECIES, SITE_ID, M_YEAR, DAY_SINCE)][
                                                        maxCOUNT == -Inf, maxCOUNT := NA]

              ts_season_count_summary <- unique(ts_season_count[, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                                    WEEK_DAY, WEEK_SINCE, DAY_SINCE, trimDAYNO, 
                                                                    M_SEASON, COMPLT_SEASON, ANCHOR, maxCOUNT)])[, 
                                                                    COUNT := maxCOUNT][, maxCOUNT := NULL]
            } else {
                ts_season_count[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
                ts_season_count <- ts_season_count[, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                        WEEK_DAY, WEEK_SINCE, trimWEEKNO, 
                                                        M_SEASON, COMPLT_SEASON, ANCHOR, COUNT)]

                ts_season_count[ , maxCOUNT := suppressWarnings(max(COUNT, na.rm = TRUE)),
                                                      by = .(SPECIES, SITE_ID, M_YEAR, WEEK_SINCE)][
                                                        maxCOUNT == -Inf, maxCOUNT := NA]

                ts_season_count_summary <- ts_season_count[WEEK_DAY == weekday, .(SPECIES, SITE_ID, YEAR, M_YEAR, MONTH, DAY, WEEK, 
                                                                            WEEK_DAY, WEEK_SINCE, trimWEEKNO, 
                                                                            M_SEASON, COMPLT_SEASON, ANCHOR, maxCOUNT)][, 
                                                                            COUNT := maxCOUNT][, maxCOUNT := NULL]
                # make trimWEEKNO start with 1, it might start with 2 if the selected weekday was not in the first week 
                ts_season_count_summary[, trimWEEKNO := WEEK_SINCE - min(WEEK_SINCE) + 1, by = M_YEAR]
            }
  }

return(ts_season_count_summary)
}
