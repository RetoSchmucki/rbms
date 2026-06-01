# get_nm_optimised Optimised version: find the nearest year with a computed flight curve

get_nm_optimised Optimised version: find the nearest year with a
computed flight curve

## Usage

``` r
get_nm_optimised(
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

  integer Year for which to compute the flight curve.

- ts_season_count:

  data.table Time-series of count and season information
  (pre-processed).

- MinVisit:

  integer The minimum number of visits required.

- MinOccur:

  integer The minimum number of positive records.

- MinNbrSite:

  integer The minimum number of sites required.

- NbrSample:

  integer Maximum number of sites to sample.

- GamFamily:

  string Distribution family for GAM.

- MaxTrial:

  integer Maximum number of trials.

- SpeedGam:

  logical Use bam instead of gam.

- OptiGam:

  logical Optimize choice between bam and gam.

- TimeUnit:

  character Time unit 'd' or 'w'.

- MultiVisit:

  string Aggregation function.

- weekday:

  integer Weekday selection.

- mod_form:

  string Model formula.

- tp_col:

  string Temporal column name.

- verbose:

  logical Verbose output.

- ...:

  Additional parameters.

## Value

A list with flight curve data, model, and input data.

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
