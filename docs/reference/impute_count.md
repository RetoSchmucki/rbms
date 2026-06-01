# impute_count

impute_count

## Usage

``` r
impute_count(
  ts_season_count,
  ts_flight_curve,
  TimeUnit = "d",
  MultiVisit = "mean",
  weekday = 3,
  sp = NULL,
  YearLimit = NULL,
  SelectYear = NULL,
  CompltSeason = TRUE,
  verbose = TRUE
)
```

## Arguments

- ts_season_count:

  data.table Time-series of counts for a specific species across all
  sites as returned by [ts_monit_count_site](ts_monit_count_site.md).

- ts_flight_curve:

  pheno_curve class object or data.table with light curves and relative
  abundances (NM) for a specific species as returned by
  [flight_curve](flight_curve.md).

- TimeUnit:

  character The time-step for which the spline should be computed, 'd'
  day or 'w' week.

- MultiVisit:

  string Function to apply for summarising multiple counts within a time
  unit, 'max' or 'mean' (default).

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- sp:

  integer or string Species ID or name.

- YearLimit:

  integer Define the span of years (+/- number of year) to look for a
  flight curve, if NULL no restriction is set.

- SelectYear:

  integer Select a specific year to compute the flight curve,
  default=NULL.

- CompltSeason:

  logical Restrict computation of flight curve for years where the
  complete season was sampled, default=TRUE.

- verbose:

  a logical indicating if some “progress report” should be given.

## Value

A data.table based on the entry count data, augmented with site indices
'SINDEX' and imputed weekly count 'IMPUTED_COUNT'.

## Details

Site indices can be extracted from the data.table returned by this
function. The site index is currently computed by adjusting the count by
the proportion of the flight curve covered by the visits.

## See also

[flight_curve](flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
