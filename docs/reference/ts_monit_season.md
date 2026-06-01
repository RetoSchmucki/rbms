# ts_monit_season Build a time-series of dates with specific detail about the monitoring season

ts_monit_season Build a time-series of dates with specific detail about
the monitoring season

## Usage

``` r
ts_monit_season(
  d_series,
  StartMonth = 4,
  EndMonth = 9,
  StartDay = 1,
  EndDay = NULL,
  CompltSeason = TRUE,
  Anchor = TRUE,
  AnchorLength = 7,
  AnchorLag = 7,
  TimeUnit = "d"
)
```

## Arguments

- d_series:

  time series of dates returned by [ts_dwmy_table](ts_dwmy_table.md)

- StartMonth:

  integer for the first month of the monitoring season, default=4

- EndMonth:

  integer for the last month of the monitoring season, default=9

- StartDay:

  integer for the day of the month when the monitoring season start,
  default=1

- EndDay:

  integer for the day of the month when the monitoring season end,
  default=last day of the EndMonth

- CompltSeason:

  logical if only the date for a complete season should be returned,
  default=TRUE

- Anchor:

  logical if Anchor should be used at the beginning and end of the
  monitoring season, default=TRUE

- AnchorLength:

  integer for the number of day used as Anchor, default=7

- AnchorLag:

  integer for the number of days before and after where the Anchor
  should start, default=TRUE

- TimeUnit:

  character defining the time unit used for the Anchor, can be 'd' for
  day or 'w' for week

## Value

A data.table with the entire time-series of date (`DATE`, `DAY_SINCE`,
`YEAR`, `MONTH`, `DAY`, `WEEK` (ISO), the `WEEK_DAY`, and details about
the Monitoring Year (`M_YEAR`) which is the monitoring year to which
that date refers to, using the year of the starting month of the
monitoring season, the monitoring season (`M_SEASON`), if the Monitoring
season is complete (`COMPLT_SEASON`), the location of the anchors
(`ANCHOR`) and the count for every day used as anchor (0).

## Details

The monitoring season can start in year y and end in year y+1, so if the
monitoring season goes over two years (i.e. winter, December-January)
the SEASON YEAR is shifted to set the monitoring season within a
continuous year named according to the year at the start. Here the
monitoring season is set in the middle or year' to give some room to set
ANCHORs before and after the monitoring season.

## See also

[ts_date_seq](ts_date_seq.md), [ts_dwmy_table](ts_dwmy_table.md),
[set_anchor](set_anchor.md)

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
