##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   29.12.2025
##
##==========================================


#' fit_bam_multi
#' Fit a Generalized Additive Model using bam() with factor-smooth interactions to butterfly count data across multiple years simultaneously.
#' This function fits a single model across multiple years using factor-smooth interactions (fs), which allows year-specific smooth terms
#' while sharing information across years. This is more efficient than fitting separate models per year.
#' @param dataset_multi data.table Butterfly counts for species x over multiple years y over all sites.
#' @param NbrSample integer Maximum number of sites to sample per year, default=NULL (use all sites).
#' @param GamFamily string Family for GAM, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param MaxTrial integer Maximum number of fitting attempts, default=4.
#' @param TimeUnit character Time-step for the spline, 'd' for day or 'w' for week.
#' @param MultiVisit string Function for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param tp_col string Name of temporal variable used in the GAM model, default NULL.
#' @param smooth_basis string Basis for the smooth term, default "cr" (cubic regression spline). Other options include "tp" (thin plate), "ps" (P-splines).
#' @param k_value integer Basis dimension for the smooth, default=-1 (automatic selection).
#' @param verbose logical Indicating if progress reports should be given.
#' @param ... Additional parameters passed to bam function from the \link[mgcv]{bam} package.
#' @return A list with three objects:
#'   i) **f_curve**: a data.table with the flight curves for all years with expected relative abundance (NM), normalized to sum to one per year,
#'   ii) **f_model**: the resulting bam model fitted on the count data across all years,
#'   iii) **f_data**: a data.table with the data used to fit the BAM model.
#' @details This function uses factor-smooth interactions (s(tp_col, M_YEAR, bs="fs")) which allows:
#'   - Fitting a single model for multiple years (more efficient)
#'   - Year-specific smooths that can differ in shape
#'   - Sharing of information across years for more stable estimates
#'   - Better handling of years with sparse data
#' @keywords gam, bam, factor-smooth
#' @seealso \link{flight_curve_multi}, \link{fit_gam}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @importFrom data.table uniqueN setkey setkeyv setnames rbindlist
#' @importFrom stats as.formula predict
#' @export fit_bam_multi
#'

fit_bam_multi <- function(dataset_multi, NbrSample = NULL, GamFamily = 'poisson', MaxTrial = 4,
                          TimeUnit = 'd', MultiVisit = 'mean', weekday = 3, tp_col = NULL,
                          smooth_basis = "cr", k_value = -1, verbose = TRUE, ...){

  check_package('data.table')
  check_package('mgcv')

  tr <- 1
  gam_obj_multi <- c()

  while((tr == 1 | inherits(gam_obj_multi, "try-error")) & tr <= MaxTrial){

    # Sample sites if NbrSample is specified
    if(!is.null(NbrSample)){
      # Sample per year to ensure balanced representation
      sampled_sites <- dataset_multi[, .(SITE_ID = sample(unique(SITE_ID),
                                                          min(uniqueN(SITE_ID), NbrSample),
                                                          replace = FALSE)),
                                    by = M_YEAR]
      sp_data_all <- data.table::copy(dataset_multi[sampled_sites, on = .(M_YEAR, SITE_ID)])
    } else {
      sp_data_all <- data.table::copy(dataset_multi)
    }

    # Determine temporal column if not provided
    if(is.null(tp_col)){
      if(TimeUnit == 'd'){
        tp_col <- "trimDAYNO"
      } else {
        tp_col <- "trimWEEKNO"
      }
    }

    # Ensure M_YEAR is a factor for factor-smooth interaction
    sp_data_all[, M_YEAR_fac := as.factor(M_YEAR)]

    if(verbose){
      message(paste("Fitting multi-year BAM with factor-smooth interaction for species",
                   as.character(sp_data_all$SPECIES[1]),
                   "across", sp_data_all[, uniqueN(M_YEAR)], "years with",
                   sp_data_all[, uniqueN(SITE_ID)], "sites, using bam():",
                   Sys.time(), "-> trial", tr))
    }

    # Build formula with factor-smooth interaction
    # s(tp_col, M_YEAR_fac, bs="fs") creates year-specific smooths
    if(sp_data_all[, uniqueN(SITE_ID)] > 1){
      if(k_value == -1){
        mod_form <- as.formula(paste0("COUNT ~ s(", tp_col, ", M_YEAR_fac, bs='fs', xt=list(bs='", smooth_basis, "')) + factor(SITE_ID)"))
      } else {
        mod_form <- as.formula(paste0("COUNT ~ s(", tp_col, ", M_YEAR_fac, bs='fs', k=", k_value, ", xt=list(bs='", smooth_basis, "')) + factor(SITE_ID)"))
      }
    } else {
      if(k_value == -1){
        mod_form <- as.formula(paste0("COUNT ~ s(", tp_col, ", M_YEAR_fac, bs='fs', xt=list(bs='", smooth_basis, "'))"))
      } else {
        mod_form <- as.formula(paste0("COUNT ~ s(", tp_col, ", M_YEAR_fac, bs='fs', k=", k_value, ", xt=list(bs='", smooth_basis, "'))"))
      }
    }

    # Fit BAM model
    gam_obj_multi <- try(mgcv::bam(mod_form, data = sp_data_all, family = GamFamily, ...), silent = TRUE)

    tr <- tr + 1
  }

  # Check if model fitting failed
  if (inherits(gam_obj_multi, "try-error")) {
    if(verbose) {
      message(paste("Error in fitting the multi-year BAM for species",
                   as.character(sp_data_all$SPECIES[1]),
                   "; Model did not converge after", tr, "trials"))
    }
    sp_data_all[, c("FITTED", "NM") := .(NA, NA)]

    f_curve <- unique(sp_data_all[, -c("COUNT", "FITTED")])
    f_curve_mod <- list(f_curve = f_curve, f_model = gam_obj_multi, f_data = sp_data_all)

  } else {
    # Predict fitted values
    sp_data_all[, FITTED := mgcv::predict.bam(gam_obj_multi, newdata = sp_data_all, type = "response")]
    sp_data_all[M_SEASON == 0L, FITTED := 0]

    # Check for infinite fitted values
    if(sum(is.infinite(sp_data_all[, FITTED])) > 0){
      sp_data_all[, c("FITTED", "NM") := .(NA, NA)]
    } else {
      # Compute normalized values per site and year
      sp_data_all[, SITE_YEAR_SUM := sum(FITTED), by = .(SITE_ID, M_YEAR)]
      sp_data_all[, NM := round(FITTED / SITE_YEAR_SUM, 5)]
    }

    # Prepare flight curve output
    f_curve <- data.table::copy(sp_data_all)

    # Average across sites for each year
    cols_to_keep <- c("SPECIES", "M_YEAR", tp_col, "M_SEASON", "NM", "FITTED")
    if(TimeUnit == 'd'){
      cols_to_keep <- c(cols_to_keep, "YEAR", "MONTH", "DAY", "WEEK", "WEEK_DAY",
                       "WEEK_SINCE", "DAY_SINCE", "COMPLT_SEASON", "ANCHOR")
    } else {
      cols_to_keep <- c(cols_to_keep, "YEAR", "MONTH", "DAY", "WEEK", "WEEK_DAY",
                       "WEEK_SINCE", "COMPLT_SEASON", "ANCHOR")
    }

    # Keep only relevant columns and aggregate across sites
    available_cols <- intersect(cols_to_keep, names(f_curve))
    f_curve <- f_curve[, ..available_cols]

    # Average NM and FITTED across sites for plotting
    agg_cols <- setdiff(available_cols, c("NM", "FITTED"))
    f_curve <- f_curve[, .(NM = mean(NM, na.rm = TRUE),
                          FITTED = mean(FITTED, na.rm = TRUE)),
                      by = agg_cols]

    data.table::setkey(f_curve)

    # Clean up sp_data_all
    sp_data_all[, c("FITTED", "SITE_YEAR_SUM", "NM", "M_YEAR_fac") := NULL]

    f_curve_mod <- list(f_curve = f_curve, f_model = gam_obj_multi, f_data = sp_data_all)
  }

  return(f_curve_mod)
}


#' flight_curve_multi
#' Compute flight curves across multiple years simultaneously using a single BAM model with factor-smooth interactions.
#' This is more efficient than fitting separate models per year and can provide better estimates for years with sparse data.
#' @param ts_season_count data.table Time-series of counts and season information returned by \link{ts_monit_count_site}
#' @param NbrSample integer Maximum number of sites to use per year, default=100.
#' @param MinVisit integer Minimum number of visits required for a site to be included, default=3.
#' @param MinOccur integer Minimum number of positive records required per year in a site, default=2.
#' @param MinNbrSite integer Minimum number of sites required per year, default=1.
#' @param MaxTrial integer Maximum number of trials for model convergence, default=3.
#' @param GamFamily string Distribution for GAM, default='poisson', but can be 'nb' or 'quasipoisson'.
#' @param CompltSeason logical Restrict to years where complete season was sampled, default=TRUE.
#' @param SelectYear integer Vector of specific years to include, default=NULL (use all).
#' @param TimeUnit character Time-step for the spline, 'd' for day or 'w' for week.
#' @param MultiVisit string Function for summarising multiple counts within a time unit, 'max' or 'mean' (default).
#' @param weekday Integer for selected day of the week for weekly summary, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param tp_col string Name of temporal variable, default NULL (auto-determined from TimeUnit).
#' @param smooth_basis string Basis for smooth term, default "cr". Options: "tp", "ps", "cr".
#' @param k_value integer Basis dimension, default=-1 (automatic).
#' @param KeepModel logical Keep model output, default=TRUE.
#' @param KeepModelData logical Keep data used for fitting, default=TRUE.
#' @param verbose logical Show progress reports, default=TRUE.
#' @param ... Additional parameters passed to bam function.
#' @return A list with three objects:
#'   i) **pheno**: flight curves with expected relative abundance (NM) across all years,
#'   ii) **model**: the bam model object fitted across all years,
#'   iii) **data**: the data used to fit the model (if KeepModelData=TRUE).
#' @details This function uses mgcv::bam with factor-smooth interactions which is particularly useful for:
#'   - Large datasets with many years
#'   - Years with varying data quality/quantity
#'   - Sharing information across years while allowing year-specific patterns
#' @keywords gam, flight curve, factor-smooth
#' @seealso \link{fit_bam_multi}, \link{flight_curve}, \link[mgcv]{bam}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @importFrom data.table uniqueN setkey setkeyv rbindlist
#' @export flight_curve_multi
#'

flight_curve_multi <- function(ts_season_count, NbrSample = 100, MinVisit = 3, MinOccur = 2,
                               MinNbrSite = 1, MaxTrial = 3, GamFamily = 'poisson',
                               CompltSeason = TRUE, SelectYear = NULL, TimeUnit = 'd',
                               MultiVisit = "mean", weekday = 3, tp_col = NULL,
                               smooth_basis = "cr", k_value = -1,
                               KeepModel = TRUE, KeepModelData = TRUE, verbose = TRUE, ...) {

  check_package('data.table')
  check_package('mgcv')

  names(ts_season_count) <- toupper(names(ts_season_count))
  if(!is.null(tp_col)){
    tp_col <- toupper(tp_col)
  }

  col_name <- c("COMPLT_SEASON", "M_YEAR", "SITE_ID", "SPECIES", "DATE", "WEEK",
                "WEEK_DAY", "DAY_SINCE", "WEEK_SINCE", "M_SEASON", "COUNT", "ANCHOR")
  check_names(ts_season_count, ifelse(is.null(tp_col), col_name, c(col_name, tp_col)))

  # Filter for complete seasons if required
  if(isTRUE(CompltSeason)){
    ts_season_count <- ts_season_count[COMPLT_SEASON == 1, ]
  }

  # Select years
  if(is.null(SelectYear)){
    year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, unique(as.integer(M_YEAR))]
  } else {
    year_series <- ts_season_count[!is.na(COUNT) & ANCHOR != 1, ][M_YEAR %in% SelectYear, unique(as.integer(M_YEAR))]
  }

  if(length(year_series) == 0) stop(paste0("No count data found for selected years"))

  # Handle temporal columns
  if(!is.null(tp_col)){
    if(!all(tp_col %in% col_name)){
      if(TimeUnit == 'd'){
        ts_season_count_tp_col <- unique(ts_season_count[, c(c("SPECIES", "SITE_ID", "DAY_SINCE"), tp_col), with = FALSE])
        setkey(ts_season_count_tp_col, SPECIES, SITE_ID, DAY_SINCE)
      } else {
        ts_season_count_m <- ts_season_count[, lapply(.SD, mean, na.rm = TRUE),
                                             by = .(SPECIES, SITE_ID, WEEK_SINCE),
                                             .SDcols = tp_col]
        ts_season_count_tp_col <- unique(ts_season_count_m[, c(c("SPECIES", "SITE_ID", "WEEK_SINCE"), tp_col), with = FALSE])
        setkey(ts_season_count_tp_col, SPECIES, SITE_ID, WEEK_SINCE)
      }
    }
  }

  # Summarize by day or week
  ts_season_count <- day_week_summary(ts_season_count, MultiVisit = MultiVisit,
                                     TimeUnit = TimeUnit, weekday = weekday)

  # Merge temporal columns if needed
  if(!is.null(tp_col)){
    if(!all(tp_col %in% col_name)){
      ifelse(TimeUnit == 'd',
             setkey(ts_season_count, SPECIES, SITE_ID, DAY_SINCE),
             setkey(ts_season_count, SPECIES, SITE_ID, WEEK_SINCE))
      ts_season_count <- merge(ts_season_count, ts_season_count_tp_col, all.x = TRUE)
    }
  }

  # Filter sites based on visit and occurrence criteria
  dataset_multi_list <- lapply(year_series, function(y) {
    dataset_y <- ts_season_count[as.integer(M_YEAR) == y, ]
    visit_occ_site <- unique(dataset_y[!is.na(COUNT) & ANCHOR == 0L, visitN := .N, by = SITE_ID][
      !is.na(COUNT) & ANCHOR == 0L & COUNT > 0, occurN := .N, by = SITE_ID][
        !is.na(COUNT) & ANCHOR == 0L, ][order(SITE_ID), .(SITE_ID, occurN, visitN)])

    dataset_y_filtered <- data.table::copy(dataset_y[SITE_ID %in%
                                                       visit_occ_site[visitN >= MinVisit & occurN >= MinOccur, SITE_ID], ][,
                                                                                                                              visitN := NULL][, occurN := NULL])

    if(dataset_y_filtered[, uniqueN(SITE_ID)] < MinNbrSite){
      return(NULL)
    } else {
      return(dataset_y_filtered)
    }
  })

  # Remove NULL entries (years without enough sites)
  dataset_multi_list <- Filter(Negate(is.null), dataset_multi_list)

  if(length(dataset_multi_list) == 0){
    stop("No years have sufficient sites meeting MinVisit, MinOccur, and MinNbrSite criteria")
  }

  # Combine all years into single dataset
  dataset_multi <- data.table::rbindlist(dataset_multi_list, use.names = TRUE, fill = TRUE)

  # Fit single multi-year BAM model
  result_fc <- fit_bam_multi(dataset_multi, NbrSample = NbrSample, GamFamily = GamFamily,
                            MaxTrial = MaxTrial, TimeUnit = TimeUnit, MultiVisit = MultiVisit,
                            weekday = weekday, tp_col = tp_col, smooth_basis = smooth_basis,
                            k_value = k_value, verbose = verbose, ...)

  # Prepare output based on KeepModel and KeepModelData settings
  if(!isTRUE(KeepModelData) & !isTRUE(KeepModel)){
    result_fc <- list(pheno = result_fc$f_curve)
  }

  if(!isTRUE(KeepModelData) & isTRUE(KeepModel)){
    result_fc <- list(pheno = result_fc$f_curve, model = result_fc$f_model)
  }

  if(isTRUE(KeepModelData) & isTRUE(KeepModel)){
    result_fc <- list(pheno = result_fc$f_curve, model = result_fc$f_model, data = result_fc$f_data)
  }

  class(result_fc) <- "pheno_curve"
  return(result_fc)
}


#' plot_flight_curve_multi
#' Plot annual flight curves from multi-year BAM model with each year displayed in a different color.
#' @param ts_flight_curve_multi pheno_curve object returned by \link{flight_curve_multi} or a data.table with flight curve data.
#' @param SelectYear integer or vector of years to plot, default NULL (plot all years).
#' @param TimeUnit character Time-step used, 'd' for day or 'w' for week, default 'd'.
#' @param weekday Integer for selected day of the week for weekly data, default is 3 (Wednesday). [1-7] where 1 = Monday.
#' @param color_palette string Name of color palette to use. Options: "rainbow", "heat", "terrain", "topo", "cm", "viridis", or a vector of color names/hex codes.
#' @param add_legend logical Add legend to plot, default TRUE.
#' @param legend_pos string Position of legend: "topright", "topleft", "bottomright", "bottomleft", "top", "bottom", "left", "right", default "topright".
#' @param line_width numeric Width of lines, default 2.
#' @param ymax numeric Maximum value for y-axis, default NULL (auto-scale).
#' @param ylab string Label for y-axis, default "Relative abundance (NM)".
#' @param xlab string Label for x-axis, default "Date".
#' @param main string Main title for plot, default NULL.
#' @param ... Additional parameters passed to base plot function.
#' @return Returns a base plot with relative abundance (y) over time (x), with each year in a different color.
#' @keywords flight curve, visualization
#' @seealso \link{flight_curve_multi}, \link{fit_bam_multi}
#' @author Reto Schmucki - \email{retoshm@@ceh.ac.uk}
#' @import data.table
#' @importFrom grDevices rainbow heat.colors terrain.colors topo.colors cm.colors hcl.colors
#' @importFrom graphics lines legend
#' @export plot_flight_curve_multi
#'

plot_flight_curve_multi <- function(ts_flight_curve_multi, SelectYear = NULL, TimeUnit = 'd',
                                   weekday = 3, color_palette = "rainbow", add_legend = TRUE,
                                   legend_pos = "topright", line_width = 2, ymax = NULL,
                                   ylab = "Relative abundance (NM)", xlab = "Date",
                                   main = NULL, ...) {

  check_package('data.table')

  # Extract pheno data if input is pheno_curve object
  if (inherits(ts_flight_curve_multi, "pheno_curve")) {
    pheno_data <- data.table::copy(ts_flight_curve_multi$pheno)
  } else {
    pheno_data <- data.table::copy(ts_flight_curve_multi)
  }

  # Filter years if specified
  if (!is.null(SelectYear)) {
    pheno_data <- pheno_data[M_YEAR %in% SelectYear, ]
  }

  if (nrow(pheno_data) == 0) stop("No data available for selected years")

  # Get unique years
  years <- sort(unique(pheno_data$M_YEAR))
  n_years <- length(years)

  # Handle DATE column
  if ("DAY" %in% names(pheno_data)) {
    pheno_data[, DATE := data.table::as.IDate(as.Date(paste(YEAR, MONTH, DAY, sep = "-")))]
  } else {
    # For weekly data, construct dates
    f_y <- pheno_data[order(YEAR), YEAR][1]
    l_y <- pheno_data[rev(order(YEAR)), YEAR][1]
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

    setkey(pheno_data, YEAR, MONTH, WEEK)
    setkey(iso_week, YEAR, MONTH, WEEK)

    pheno_data <- merge(pheno_data, iso_week[WDAY == weekday, .(YEAR, MONTH, WEEK, DATE)], all.x = TRUE)

    # Determine temporal column
    if (TimeUnit == 'd') {
      tp_col <- "DAY_SINCE"
    } else {
      tp_col <- "WEEK_SINCE"
    }

    setkey(pheno_data, eval(as.name(tp_col)))
    f_date <- which(is.na(pheno_data[, DATE]))
    if (length(f_date) > 0) {
      if (f_date[1] != 1) {
        pheno_data[f_date, DATE2 := pheno_data[f_date - 1, DATE + 7]]
      } else {
        pheno_data[f_date, DATE2 := pheno_data[f_date + 1, DATE - 7]]
      }
      pheno_data[is.na(DATE), DATE := DATE2][, DATE2 := NULL]
    }
  }

  # Set up color palette
  if (length(color_palette) == 1) {
    colors <- switch(color_palette,
                    "rainbow" = rainbow(n_years),
                    "heat" = heat.colors(n_years),
                    "terrain" = terrain.colors(n_years),
                    "topo" = topo.colors(n_years),
                    "cm" = cm.colors(n_years),
                    "viridis" = hcl.colors(n_years, palette = "viridis"),
                    rainbow(n_years)  # default fallback
    )
  } else {
    # User provided custom colors
    if (length(color_palette) < n_years) {
      warning("Not enough colors provided, recycling color palette")
      colors <- rep_len(color_palette, n_years)
    } else {
      colors <- color_palette[1:n_years]
    }
  }

  # Calculate y-axis limits
  if (is.null(ymax)) {
    ymax <- max(pheno_data$NM, na.rm = TRUE) * 1.1
  }

  # Create empty plot
  plot(1, type = "n",
       xlim = range(pheno_data$DATE, na.rm = TRUE),
       ylim = c(0, ymax),
       xlab = xlab,
       ylab = ylab,
       main = main,
       ...)

  # Plot each year
  for (i in seq_along(years)) {
    year_data <- pheno_data[M_YEAR == years[i], ]
    year_data <- year_data[order(DATE), ]

    lines(year_data$DATE, year_data$NM,
          col = colors[i],
          lwd = line_width)
  }

  # Add legend
  if (add_legend) {
    legend(legend_pos,
           legend = as.character(years),
           col = colors,
           lwd = line_width,
           bty = "n",
           title = "Year")
  }

  invisible(NULL)
}
