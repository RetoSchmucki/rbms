##==========================================
## Optimised version of flight_curve function
##
## This optimized version maintains identical results
## but improves performance through:
## 1. Reduced data copies
## 2. Pre-computed tp_col
## 3. Optimized filtering logic
## 4. Better use of data.table operations
## 5. Parallel processing option
##
##      Date:   29.12.2025
##
##==========================================


#' flight_curve_optimised
#' Optimised version: Compute the annual flight curve from butterfly count data collated across sites.
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
#' @param parallel logical Use parallel processing for fitting models across years, default FALSE.
#' @param n_cores integer Number of cores to use if parallel=TRUE, default NULL (uses all available cores - 1).
#' @param verbose a logical indicating if some "progress report" should be given.
#' @param ... Additional parameters passed to gam or bam function from the \link[mgcv]{gam} package.
#' @return A list with three objects, i) **pheno**: a vector with annual flight curves \code{f_pheno} with expected relative abundance, normalize to sum to one over a full season,
#'         ii) **model**: a list of the resulting gam models \code{f_model} fitted on the count data for each year and iii) **data**: a data.table with the data used to fit the GAM model.
#' @keywords gam, flight curve
#' @seealso \link{flight_curve}, \link{fit_gam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export flight_curve_optimised
#'

flight_curve_optimised <- function(ts_season_count, NbrSample = 100, MinVisit = 3, MinOccur = 2, MinNbrSite = 1, MaxTrial = 3,
                         GamFamily = 'poisson', CompltSeason = TRUE, SelectYear = NULL, SpeedGam = TRUE,
                         OptiGam = TRUE, KeepModel = TRUE, KeepModelData = TRUE, TimeUnit = 'd',
                         MultiVisit = "mean", weekday = 3, mod_form = NULL, tp_col = NULL,
                         parallel = FALSE, n_cores = NULL, verbose = TRUE, ...) {

    check_package('data.table')

    # OPTIMIZATION 1: Work with reference to avoid unnecessary copy
    # Only uppercase names in place
    setnames(ts_season_count, names(ts_season_count), toupper(names(ts_season_count)))
    if(!is.null(tp_col)){
      tp_col <- toupper(tp_col)
    }

    col_name <- c("COMPLT_SEASON", "M_YEAR", "SITE_ID", "SPECIES", "DATE", "WEEK", "WEEK_DAY", "DAY_SINCE", "WEEK_SINCE", "M_SEASON", "COUNT", "ANCHOR")
    check_names(ts_season_count, ifelse(is.null(tp_col), col_name, c(col_name, tp_col)))

    # OPTIMIZATION 2: Filter in-place with subset instead of creating new object
    if(isTRUE(CompltSeason)){
        ts_season_count <- ts_season_count[COMPLT_SEASON == 1, ]
    }

    # OPTIMIZATION 3: Simplified year extraction
    if(is.null(SelectYear)){
        year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, unique(as.integer(M_YEAR))]
    } else {
        year_series <- unique(as.integer(SelectYear[SelectYear %in% ts_season_count[!is.na(COUNT) & ANCHOR != 1, M_YEAR]]))
    }

    if(length(year_series) == 0) stop(paste0(" No count data found for year ", SelectYear))

    # OPTIMIZATION 4: Handle custom temporal columns BEFORE day_week_summary if provided by user
    # (trimDAYNO and trimWEEKNO will be created BY day_week_summary)
    if(!is.null(tp_col)){
      if(!all(tp_col %in% col_name)){
        if(TimeUnit == 'd'){
          ts_season_count_tp_col <- unique(ts_season_count[, c("SPECIES", "SITE_ID", "DAY_SINCE", tp_col), with = FALSE])
          setkey(ts_season_count_tp_col, SPECIES, SITE_ID, DAY_SINCE)
        } else {
          ts_season_count_m <- ts_season_count[, lapply(.SD, mean, na.rm = TRUE), by = .(SPECIES, SITE_ID, WEEK_SINCE), .SDcols = tp_col]
          ts_season_count_tp_col <- unique(ts_season_count_m[, c("SPECIES", "SITE_ID", "WEEK_SINCE", tp_col), with = FALSE])
          setkey(ts_season_count_tp_col, SPECIES, SITE_ID, WEEK_SINCE)
        }
      }
    }

    # OPTIMIZATION 5: Call day_week_summary once before loop (creates trimDAYNO/trimWEEKNO)
    ts_season_count <- day_week_summary(ts_season_count, MultiVisit = MultiVisit, TimeUnit = TimeUnit, weekday = weekday)

    # OPTIMIZATION 6: Merge custom temporal columns once before loop (if provided)
    if(!is.null(tp_col)){
      if(!all(tp_col %in% col_name)){
        if(TimeUnit == 'd'){
          setkey(ts_season_count, SPECIES, SITE_ID, DAY_SINCE)
        } else {
          setkey(ts_season_count, SPECIES, SITE_ID, WEEK_SINCE)
        }
        ts_season_count <- merge(ts_season_count, ts_season_count_tp_col, all.x = TRUE)
      }
    }

    # OPTIMIZATION 7: Set tp_col to default trim column if not provided by user
    if(is.null(tp_col)){
      tp_col <- ifelse(TimeUnit == 'd', "trimDAYNO", "trimWEEKNO")
    }

    # OPTIMIZATION 8: Use parallel processing if requested
    if(isTRUE(parallel)){
      check_package('parallel')

      if(is.null(n_cores)){
        n_cores <- max(1, parallel::detectCores() - 1)
      }

      if(verbose){
        message(paste("Using parallel processing with", n_cores, "cores"))
      }

      # Create cluster
      cl <- parallel::makeCluster(n_cores)

      # Ensure clean exit on error
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export required objects to cluster
      parallel::clusterExport(cl, c("ts_season_count", "MinVisit", "MinOccur", "MinNbrSite",
                                   "NbrSample", "GamFamily", "MaxTrial", "SpeedGam", "OptiGam",
                                   "TimeUnit", "MultiVisit", "weekday", "mod_form", "tp_col"),
                            envir = environment())

      # Load required packages on each worker first
      parallel::clusterEvalQ(cl, {
        library(data.table)
        library(mgcv)
      })

      # Load rbms package and required functions on workers
      load_success <- parallel::clusterEvalQ(cl, {
        # Try to load rbms package
        if(requireNamespace("rbms", quietly = TRUE)) {
          library(rbms)
          TRUE
        } else {
          FALSE
        }
      })

      # If rbms couldn't be loaded on workers, export functions manually
      if(!all(unlist(load_success))) {
        if(verbose) message("Exporting functions to workers manually")

        # Get functions from current environment or global environment
        get_nm_opt_func <- if(exists("get_nm_optimised")) get("get_nm_optimised") else get_nm_optimised
        fit_gam_func <- if(exists("fit_gam")) get("fit_gam") else fit_gam

        parallel::clusterExport(cl, c("get_nm_opt_func", "fit_gam_func"), envir = environment())

        # Assign to proper names on workers
        parallel::clusterEvalQ(cl, {
          get_nm_optimised <- get_nm_opt_func
          fit_gam <- fit_gam_func
        })
      }

      # Run get_nm_optimised in parallel
      result_fc <- tryCatch({
        parallel::parLapply(cl, year_series, get_nm_optimised,
                          ts_season_count = ts_season_count,
                          MinVisit = MinVisit, MinOccur = MinOccur, MinNbrSite = MinNbrSite,
                          NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial,
                          SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit,
                          MultiVisit = MultiVisit, weekday = weekday, mod_form = mod_form,
                          tp_col = tp_col, verbose = verbose, ...)
      }, error = function(e) {
        parallel::stopCluster(cl)
        on.exit()  # Remove the on.exit handler
        stop("Parallel processing failed: ", e$message)
      })

      parallel::stopCluster(cl)
      on.exit()  # Remove the on.exit handler after successful completion

    } else {
      # OPTIMIZATION 9: Use lapply with optimized get_nm function
      result_fc <- lapply(year_series, get_nm_optimised, ts_season_count = ts_season_count,
                         MinVisit = MinVisit, MinOccur = MinOccur, MinNbrSite = MinNbrSite,
                         NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial,
                         SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit,
                         MultiVisit = MultiVisit, weekday = weekday, mod_form = mod_form,
                         tp_col = tp_col, verbose = verbose, ...)
    }

    # OPTIMIZATION 10: Use rbindlist with fill=TRUE directly, avoid intermediate lapply
    result_fcurve <- data.table::rbindlist(lapply(result_fc, `[[`, "f_curve"), fill = TRUE)

    if(isTRUE(KeepModelData)){
      result_fdata <- data.table::rbindlist(lapply(result_fc, `[[`, "f_data"), fill = TRUE)
    }

    if(isTRUE(KeepModel)){
      result_fmodel <- lapply(result_fc, `[[`, "f_model")
    }

    # OPTIMIZATION 11: Simplified output construction
    result_fc <- list(pheno = result_fcurve)

    if(isTRUE(KeepModel)){
      result_fc$model <- result_fmodel
    }

    if(isTRUE(KeepModelData)){
      result_fc$data <- result_fdata
    }

    class(result_fc) <- "pheno_curve"
    return(result_fc)
}


#' get_nm_optimised
#' Optimised version: find the nearest year with a computed flight curve
#' @param y integer Year for which to compute the flight curve.
#' @param ts_season_count data.table Time-series of count and season information (pre-processed).
#' @param MinVisit integer The minimum number of visits required.
#' @param MinOccur integer The minimum number of positive records.
#' @param MinNbrSite integer The minimum number of sites required.
#' @param NbrSample integer Maximum number of sites to sample.
#' @param GamFamily string Distribution family for GAM.
#' @param MaxTrial integer Maximum number of trials.
#' @param SpeedGam logical Use bam instead of gam.
#' @param OptiGam logical Optimize choice between bam and gam.
#' @param TimeUnit character Time unit 'd' or 'w'.
#' @param MultiVisit string Aggregation function.
#' @param weekday integer Weekday selection.
#' @param mod_form string Model formula.
#' @param tp_col string Temporal column name.
#' @param verbose logical Verbose output.
#' @param ... Additional parameters.
#' @return A list with flight curve data, model, and input data.
#' @keywords gam
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @export get_nm_optimised
#'

get_nm_optimised <- function(y, ts_season_count, MinVisit, MinOccur, MinNbrSite, NbrSample, GamFamily, MaxTrial,
                            SpeedGam, OptiGam, TimeUnit, MultiVisit, weekday, mod_form, tp_col, verbose = TRUE, ...){

  # OPTIMIZATION 1: Use direct subsetting without intermediate copy
  dataset_y <- ts_season_count[as.integer(M_YEAR) == y, ]

  # OPTIMIZATION 2: Combine visit and occurrence counting in single pass
  # Calculate visitN and occurN by reference
  dataset_y[!is.na(COUNT) & ANCHOR == 0L, visitN := .N, by = SITE_ID]
  dataset_y[!is.na(COUNT) & ANCHOR == 0L & COUNT > 0, occurN := .N, by = SITE_ID]

  # Get unique sites that meet criteria
  valid_sites <- unique(dataset_y[!is.na(COUNT) & ANCHOR == 0L & visitN >= MinVisit & occurN >= MinOccur, SITE_ID])

  # Filter dataset to valid sites and remove temporary columns
  dataset_y <- dataset_y[SITE_ID %in% valid_sites, ][, c("visitN", "occurN") := NULL]

  # Check if we have enough sites
  if(dataset_y[, uniqueN(SITE_ID)] < MinNbrSite){

    # OPTIMIZATION 3: More efficient construction of empty result
    f_curve <- unique(ts_season_count[as.integer(M_YEAR) == y, !c("SITE_ID", "COUNT")])
    f_curve[, NM := NA]
    setkeyv(f_curve, tp_col)

    f_curve_mod <- list(f_curve = f_curve, f_model = list(NA), f_data = data.table(NA))

    if(verbose){
      message(paste("Not enough sites with observations for estimating the flight curve for species",
                   as.character(ts_season_count$SPECIES[1]), "in", y))
    }

  } else {
    # OPTIMIZATION 4: Call fit_gam with pre-processed data
    f_curve_mod <- fit_gam(dataset_y, NbrSample = NbrSample, GamFamily = GamFamily, MaxTrial = MaxTrial,
                          SpeedGam = SpeedGam, OptiGam = OptiGam, TimeUnit = TimeUnit, MultiVisit = MultiVisit,
                          weekday = weekday, mod_form = mod_form, tp_col = tp_col, verbose = verbose, ...)
  }

  return(f_curve_mod)
}
