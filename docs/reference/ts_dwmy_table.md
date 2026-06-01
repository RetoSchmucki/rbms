# ts_dwmy_table Generate a time-series of dates with day, week, month and year (dwmy) from one initial to an end years

ts_dwmy_table Generate a time-series of dates with day, week, month and
year (dwmy) from one initial to an end years

## Usage

``` r
ts_dwmy_table(
  InitYear = 1970,
  LastYear = format(Sys.Date(), "%Y"),
  WeekDay1 = "monday"
)
```

## Arguments

- InitYear:

  start year of the time-series, 4 numbers format (e.g 1987)

- LastYear:

  end year of the time-series, if not provided, current year is used
  instead

- WeekDay1:

  to start the week on Monday, use 'monday', otherwise the week start on
  Sunday

## Value

a data.table object with the date, the day since the first date, the
week since the first week, the year, the month, the day in the month,
the ISO week number and the day in the week.

## See also

[ts_date_seq](ts_date_seq.md),
[IDateTime](https://rdrr.io/pkg/data.table/man/IDateTime.html)

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
