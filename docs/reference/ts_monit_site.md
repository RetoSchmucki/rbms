# ts_monit_site Augment the time series in m_season with all sites and visits with "zeros", leaving all non visited day with and \<NA\>

ts_monit_site Augment the time series in m_season with all sites and
visits with "zeros", leaving all non visited day with and \<NA\>

## Usage

``` r
ts_monit_site(
  ts_season,
  m_visit,
  DateFormat = "%Y-%m-%d",
  expand_sy = "FALSE"
)
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

- expand_sy:

  logical for all site-year combination to be kept, if FALSE, only
  sampled site-year combination are kept.

## Value

A data.table with a complete time-series where absences have been
implemented and that is ready to receive counts

## See also

[df_visit_season](https://retoschmucki.github.io/rbms/reference/df_visit_season.md)

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
