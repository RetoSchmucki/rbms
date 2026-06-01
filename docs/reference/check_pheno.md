# check_pheno Check for the flight curve of a specific year. If the specific year is missing, use the nearest year available within a 5-year period to impute missing count. Function used in [impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md).

check_pheno Check for the flight curve of a specific year. If the
specific year is missing, use the nearest year available within a 5-year
period to impute missing count. Function used in
[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md).

## Usage

``` r
check_pheno(
  sp_count_flight,
  ts_flight_curve,
  YearCheck,
  YearLimit,
  TimeUnit,
  verbose = TRUE
)
```

## Arguments

- sp_count_flight:

  data.table Object with the flight curve, relative abundance (NM), for
  a specific year.

- ts_flight_curve:

  data.table Object with the all flight curves and relative abundance
  (NM) for the years available for the search as returned by
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md).

- YearCheck:

  integer or vector Year to check for nearest flight curve, set
  internally in
  [impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md).

- YearLimit:

  integer Define the span of years (+/- number of year) to look for a
  flight curve, if NULL no restriction is set.

- TimeUnit:

  character The time-step for which the spline should be computed, 'd'
  day or 'w' week.

- verbose:

  a logical indicating if some “progress report” should be given.

## Value

A data.table with time-series of the expected relative abundance of
butterfly count per day (NM) for the year or the nearest year where
phenology is available.

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
