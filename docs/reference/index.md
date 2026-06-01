# Package index

## All functions

- [`boot_sample()`](boot_sample.md) : boot_sample Generate n bootstrap
  sample of the monitoring sites to be used for each iteration.
- [`check_names()`](check_names.md) : check_names Verify for the
  required column names in the data object.
- [`check_package()`](check_package.md) : check_package Internal
  function to check if a package is installed.
- [`check_pheno()`](check_pheno.md) : check_pheno Check for the flight
  curve of a specific year. If the specific year is missing, use the
  nearest year available within a 5-year period to impute missing count.
  Function used in impute_count.
- [`collated_index()`](collated_index.md) : collated_index compute a
  collated index from the site indices, using a Generalized Linear
  Model.
- [`collated_index_old()`](collated_index_old.md) : collated_index_old
  compute a collated index from the site indices, using a Generalized
  Linear Model.
- [`day_week_summary()`](day_week_summary.md) : day_week_summary
  Summarize the count(s) per day or week using either the maximum or the
  average count when multiple counts (visit) are observed.
- [`df_visit_season()`](df_visit_season.md) : df_visit_season Link each
  recorded visit to a corresponding monitoring season, this function is
  used in ts_monit_site
- [`fit_bam_multi()`](fit_bam_multi.md) : fit_bam_multi Fit a
  Generalized Additive Model using bam() with factor-smooth interactions
  to butterfly count data across multiple years simultaneously. This
  function fits a single model across multiple years using factor-smooth
  interactions (fs), which allows year-specific smooth terms while
  sharing information across years. This is more efficient than fitting
  separate models per year.
- [`fit_gam()`](fit_gam.md) : fit_gam fit a Generalized Additive Model
  to butterfly count data along a temporal variable and accounting for
  site effect when multiple are available.
- [`flight_curve()`](flight_curve.md) : flight_curve Compute the annual
  flight curve from butterfly count data collated across sites.
- [`flight_curve_multi()`](flight_curve_multi.md) : flight_curve_multi
  Compute flight curves across multiple years simultaneously using a
  single BAM model with factor-smooth interactions. This is more
  efficient than fitting separate models per year and can provide better
  estimates for years with sparse data.
- [`flight_curve_optimised()`](flight_curve_optimised.md) :
  flight_curve_optimised Optimised version: Compute the annual flight
  curve from butterfly count data collated across sites.
- [`get_nm()`](get_nm.md) : get_nm Compute the normalized flight curve
  by fitting a spline in a Generalized Additive Model for one year 'y'
  to butterfly count data.
- [`get_nm_optimised()`](get_nm_optimised.md) : get_nm_optimised
  Optimised version: find the nearest year with a computed flight curve
- [`get_nny()`](get_nny.md) : get_nny find the nearest year with a
  computed flight curve
- [`impute_count()`](impute_count.md) : impute_count
- [`initiate_project()`](initiate_project.md) : initiate_project Build
  the folder structure for a generic research project
- [`m_count`](m_count.md) : Toy data set with butterfly count for x
  species across y sites
- [`m_visit`](m_visit.md) : Toy data set with the date when the sites
  have been visited for monitoring
- [`plot(`*`<pheno_curve>`*`)`](plot.pheno_curve.md) : plot.pheno_curve
  Generic method to plot the flight curve, where values are extracted
  from a pheno_curve object (outcome of the light_curve() function)
- [`plot_flight_curve_multi()`](plot_flight_curve_multi.md) :
  plot_flight_curve_multi Plot annual flight curves from multi-year BAM
  model with each year displayed in a different color.
- [`points(`*`<pheno_curve>`*`)`](points.pheno_curve.md) :
  points.pheno_curve Generic method to add points on a plot of the
  flight curve, where values are extracted from a pheno_curve object
  (outcome of the light_curve() function)
- [`set_anchor()`](set_anchor.md) : set_anchor Add Anchors of "zeros" at
  determined distance on each side of the monitoring season with
  specific weight (length), this function is used by ts_monit_season()
- [`site_index()`](site_index.md) : site_index Extract abundance indices
  per site and year based on flight curve imputation.
- [`ts_date_seq()`](ts_date_seq.md) : ts_date_seq Generate a time-series
  with dates starting from January of the starting year to December of
  the ending years.
- [`ts_dwmy_table()`](ts_dwmy_table.md) : ts_dwmy_table Generate a
  time-series of dates with day, week, month and year (dwmy) from one
  initial to an end years
- [`ts_monit_count_site()`](ts_monit_count_site.md) :
  ts_monit_count_site Generate a full time series of observed count, for
  all sites and each day since a starting and ending years of the
  defined time-series
- [`ts_monit_season()`](ts_monit_season.md) : ts_monit_season Build a
  time-series of dates with specific detail about the monitoring season
- [`ts_monit_site()`](ts_monit_site.md) : ts_monit_site Augment the time
  series in m_season with all sites and visits with "zeros", leaving all
  non visited day with and \<NA\>
