# df_visit_season Link each recorded visit to a corresponding monitoring season, this function is used in [ts_monit_site](https://retoschmucki.github.io/rbms/reference/ts_monit_site.md)

df_visit_season Link each recorded visit to a corresponding monitoring
season, this function is used in
[ts_monit_site](https://retoschmucki.github.io/rbms/reference/ts_monit_site.md)

## Usage

``` r
df_visit_season(ts_season, m_visit, DateFormat = "%Y-%m-%d")
```

## Arguments

- ts_season:

  data.table returned by
  [ts_monit_season](https://retoschmucki.github.io/rbms/reference/ts_monit_season.md)
  with the detail time-series of the monitoring season.

- m_visit:

  data.table or data.frame with Date and Site ID for each monitoring
  visit.

- DateFormat:

  format used for the date in the visit data, default="%Y-%m-%d".

## Value

A data.table where each visit is attributed a monitoring year, the
`M_YEAR` as this can differ the date (see detail)

## Details

The value of `M_YEAR` should be used as Monitoring year as this can
differ form the `YEAR` if the monitoring season covers two calendar
years (e.g. November to June).

## See also

[ts_monit_season](https://retoschmucki.github.io/rbms/reference/ts_monit_season.md)

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
