# site_index Extract abundance indices per site and year based on flight curve imputation.

site_index Extract abundance indices per site and year based on flight
curve imputation.

## Usage

``` r
site_index(butterfly_count, MinFC = NULL)
```

## Arguments

- butterfly_count:

  data.table Observed and imputed weekly or daily counts and the
  estimated total counts (SINDEX) computed by the function
  impute_count().

- MinFC:

  numeric Value between 0 and 1 to define the threshold for the
  proportion of the flight curve covered by the visits, if NULL all
  site-year available indices are returned.

## Value

data.table Estimated annual abundance index and the proportion of the
flight curve covered by the visit - total weekly or daily count over the
entire monitoring season.

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
