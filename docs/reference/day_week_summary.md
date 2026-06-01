# day_week_summary Summarize the count(s) per day or week using either the maximum or the average count when multiple counts (visit) are observed.

day_week_summary Summarize the count(s) per day or week using either the
maximum or the average count when multiple counts (visit) are observed.

## Usage

``` r
day_week_summary(ts_season_count, MultiVisit, TimeUnit, weekday)
```

## Arguments

- ts_season_count:

  data.table Time-series of butterfly counts for species x over year y
  over all sites.

- MultiVisit:

  string Function to apply for summarising multiple counts within a time
  unit, 'max' or 'mean' (default).

- TimeUnit:

  character The time-step for which the spline should be computed, 'd'
  day or 'w' week. For TimeUnit week, 'w', we refer to week_day 3 to
  match a mid-week calendar date.

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

## Value

data.table Time-series of butterfly counts for species x over year y
over all sites, but summarized with the chosen function.

## See also

[flight_curve](flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
