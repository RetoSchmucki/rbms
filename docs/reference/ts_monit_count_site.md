# ts_monit_count_site Generate a full time series of observed count, for all sites and each day since a starting and ending years of the defined time-series

ts_monit_count_site Generate a full time series of observed count, for
all sites and each day since a starting and ending years of the defined
time-series

## Usage

``` r
ts_monit_count_site(
  m_season_visit,
  m_count,
  sp = 1,
  DateFormat = "%Y-%m-%d"
)
```

## Arguments

- m_season_visit:

  data.table with potential absence based on visit data as returned by
  `ts_monit_site`

- m_count:

  data.table or data.frame with monitoring count data

- sp:

  integer or string for the species ID or name

- DateFormat:

  format used for the date in the count data, default="%Y-%m-%d".

## Value

A data.table with a complete time-series where absences (0) and recorded
counts for a specific species are implemented for each site.

## See also

[ts_monit_site](ts_monit_site.md)

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
