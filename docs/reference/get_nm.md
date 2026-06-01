# get_nm Compute the normalized flight curve by fitting a spline in a Generalized Additive Model for one year 'y' to butterfly count data.

get_nm Compute the normalized flight curve by fitting a spline in a
Generalized Additive Model for one year 'y' to butterfly count data.

## Usage

``` r
get_nm(
  y,
  ts_season_count,
  MinVisit,
  MinOccur,
  MinNbrSite,
  NbrSample,
  GamFamily,
  MaxTrial,
  SpeedGam,
  OptiGam,
  TimeUnit,
  MultiVisit,
  weekday,
  mod_form,
  tp_col,
  verbose = TRUE,
  ...
)
```

## Arguments

- y:

  integer Vector of years for which to compute the flight curve.

- ts_season_count:

  data.table Time-series of count and season information returned by
  [ts_monit_count_site](https://retoschmucki.github.io/rbms/reference/ts_monit_count_site.md)

- MinVisit:

  integer The minimum number of visits required for a site to be
  included in the computation, default=3.

- MinOccur:

  integer The minimum number of positive records (e.g. \>= 1) observed
  over the year in a site default=2.

- MinNbrSite:

  integer The minimum number of sites required to compute the flight
  curve, default=1.

- NbrSample:

  integer Value inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  when set to 'NULL' (default), all site are considered in the GAM model

- GamFamily:

  string Value inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  default='poisson', but can be 'nb' or 'quasipoisson'.

- MaxTrial:

  integer Value inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  default=3.

- SpeedGam:

  logical Set if the [bam](https://rdrr.io/pkg/mgcv/man/bam.html) method
  should be used, default instead of the default
  [gam](https://rdrr.io/pkg/mgcv/man/gam.html) method.

- OptiGam:

  logical Set the use the [bam](https://rdrr.io/pkg/mgcv/man/bam.html)
  method when data are larger than 200 and gam for smaller datasets

- TimeUnit:

  character Time-step for which the spline should be computed, 'd' day
  or 'w' week.

- MultiVisit:

  string Function for summarising multiple counts within a time unit,
  'max' or 'mean' (default).

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- mod_form:

  string with formula to be passed to the gam model, default null.

- tp_col:

  string or vector of string with additional variable used in the gam
  model, default null.

- verbose:

  a logical indicating if some “progress report” should be given.

- ...:

  Additional parameters passed to gam or bam function from the
  [gam](https://rdrr.io/pkg/mgcv/man/gam.html) package.

## Value

A list of lists, each list containing three objects, i) \*\*f_curve\*\*:
a data.table with the flight curve `f_curve` with expected relative
abundance, normalised to sum to one over a full season, ii)
\*\*f_model\*\*: the resulting gam model `f_model` fitted on the count
data and iii) \*\*f_data\*\*: a data.table with the data used to fit the
GAM model. This is provided for all year provided in 'y'.

## See also

[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
[gam](https://rdrr.io/pkg/mgcv/man/gam.html),
[bam](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
