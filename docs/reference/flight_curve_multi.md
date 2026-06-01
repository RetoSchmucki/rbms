# flight_curve_multi Compute flight curves across multiple years simultaneously using a single BAM model with factor-smooth interactions. This is more efficient than fitting separate models per year and can provide better estimates for years with sparse data.

flight_curve_multi Compute flight curves across multiple years
simultaneously using a single BAM model with factor-smooth interactions.
This is more efficient than fitting separate models per year and can
provide better estimates for years with sparse data.

## Usage

``` r
flight_curve_multi(
  ts_season_count,
  NbrSample = 100,
  MinVisit = 3,
  MinOccur = 2,
  MinNbrSite = 1,
  MaxTrial = 3,
  GamFamily = "poisson",
  CompltSeason = TRUE,
  SelectYear = NULL,
  TimeUnit = "d",
  MultiVisit = "mean",
  weekday = 3,
  tp_col = NULL,
  smooth_basis = "cr",
  k_value = -1,
  KeepModel = TRUE,
  KeepModelData = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- ts_season_count:

  data.table Time-series of counts and season information returned by
  [ts_monit_count_site](ts_monit_count_site.md)

- NbrSample:

  integer Maximum number of sites to use per year, default=100.

- MinVisit:

  integer Minimum number of visits required for a site to be included,
  default=3.

- MinOccur:

  integer Minimum number of positive records required per year in a
  site, default=2.

- MinNbrSite:

  integer Minimum number of sites required per year, default=1.

- MaxTrial:

  integer Maximum number of trials for model convergence, default=3.

- GamFamily:

  string Distribution for GAM, default='poisson', but can be 'nb' or
  'quasipoisson'.

- CompltSeason:

  logical Restrict to years where complete season was sampled,
  default=TRUE.

- SelectYear:

  integer Vector of specific years to include, default=NULL (use all).

- TimeUnit:

  character Time-step for the spline, 'd' for day or 'w' for week.

- MultiVisit:

  string Function for summarising multiple counts within a time unit,
  'max' or 'mean' (default).

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- tp_col:

  string Name of temporal variable, default NULL (auto-determined from
  TimeUnit).

- smooth_basis:

  string Basis for smooth term, default "cr". Options: "tp", "ps", "cr".

- k_value:

  integer Basis dimension, default=-1 (automatic).

- KeepModel:

  logical Keep model output, default=TRUE.

- KeepModelData:

  logical Keep data used for fitting, default=TRUE.

- verbose:

  logical Show progress reports, default=TRUE.

- ...:

  Additional parameters passed to bam function.

## Value

A list with three objects: i) \*\*pheno\*\*: flight curves with expected
relative abundance (NM) across all years, ii) \*\*model\*\*: the bam
model object fitted across all years, iii) \*\*data\*\*: the data used
to fit the model (if KeepModelData=TRUE).

## Details

This function uses mgcv::bam with factor-smooth interactions which is
particularly useful for: - Large datasets with many years - Years with
varying data quality/quantity - Sharing information across years while
allowing year-specific patterns

## See also

[fit_bam_multi](fit_bam_multi.md), [flight_curve](flight_curve.md),
[bam](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
