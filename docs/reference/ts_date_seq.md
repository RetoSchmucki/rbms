# ts_date_seq Generate a time-series with dates starting from January of the starting year to December of the ending years.

ts_date_seq Generate a time-series with dates starting from January of
the starting year to December of the ending years.

## Usage

``` r
ts_date_seq(InitYear = 1970, LastYear = format(Sys.Date(), "%Y"))
```

## Arguments

- InitYear:

  First year of the time-series, four digits format (e.g. 1987).

- LastYear:

  Last year of the time-series, if not provided, the current year is
  used.

## Value

Returns a POSIXct vector with the format 'YYYY-MM-DD HH:MM:SS'

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
