# flight_curve Compute the annual flight curve from butterfly count data collated across sites.

flight_curve Compute the annual flight curve from butterfly count data
collated across sites.

## Usage

``` r
flight_curve(
  ts_season_count,
  NbrSample = 100,
  MinVisit = 3,
  MinOccur = 2,
  MinNbrSite = 1,
  MaxTrial = 3,
  GamFamily = "poisson",
  CompltSeason = TRUE,
  SelectYear = NULL,
  SpeedGam = TRUE,
  OptiGam = TRUE,
  KeepModel = TRUE,
  KeepModelData = TRUE,
  TimeUnit = "d",
  MultiVisit = "mean",
  weekday = 3,
  mod_form = NULL,
  tp_col = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- ts_season_count:

  data.table Time-series of counts and season information returned by
  [ts_monit_count_site](ts_monit_count_site.md)

- NbrSample:

  integer The maximum number of sites to use for computing the flight
  curve, default=100.

- MinVisit:

  integer The minimum number of visits required for a site to be
  included, default=3.

- MinOccur:

  integer The minimum number of positive records (e.g. \>= 1) required
  in one year for a site to be included, default=2.

- MinNbrSite:

  integer The minimum number of sites required to compute a flight
  curve, default=1.

- MaxTrial:

  integer The maximum number of trials for model convergence, default=3.

- GamFamily:

  string The distribution of the error term used in the GAM,
  default='poisson', but can be the negative-binomial 'nb' or
  'quasipoisson'.

- CompltSeason:

  Logical to restrict computation of flight curve for years where the
  complete season has been sampled, default=TRUE.

- SelectYear:

  integer Select a specific year to compute the flight curve,
  default=NULL.

- SpeedGam:

  Logical Set if the [bam](https://rdrr.io/pkg/mgcv/man/bam.html) method
  should be used, default instead of the default
  [gam](https://rdrr.io/pkg/mgcv/man/gam.html) method.

- OptiGam:

  Logical Set the use of the
  [bam](https://rdrr.io/pkg/mgcv/man/bam.html) method when data are
  larger than 200 and gam for smaller datasets

- KeepModel:

  Logical to keep model output in a list object named
  `flight_curve_model`.

- KeepModelData:

  Logical to keep the data used for the GAM.

- TimeUnit:

  character The time-step for which the spline should be computed, 'd'
  day or 'w' week.

- MultiVisit:

  string Function to apply for summarising multiple counts within a time
  unit, 'max' or 'mean' (default).

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

A list with three objects, i) \*\*pheno\*\*: a vector with annual flight
curves `f_pheno` with expected relative abundance, normalize to sum to
one over a full season, ii) \*\*model\*\*: a list of the resulting gam
models `f_model` fitted on the count data for each year and iii)
\*\*data\*\*: a data.table with the data used to fit the GAM model.

## See also

[fit_gam](fit_gam.md)

[gam](https://rdrr.io/pkg/mgcv/man/gam.html)

[bam](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
