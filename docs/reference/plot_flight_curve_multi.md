# plot_flight_curve_multi Plot annual flight curves from multi-year BAM model with each year displayed in a different color.

plot_flight_curve_multi Plot annual flight curves from multi-year BAM
model with each year displayed in a different color.

## Usage

``` r
plot_flight_curve_multi(
  ts_flight_curve_multi,
  SelectYear = NULL,
  TimeUnit = "d",
  weekday = 3,
  color_palette = "rainbow",
  add_legend = TRUE,
  legend_pos = "topright",
  line_width = 2,
  ymax = NULL,
  ylab = "Relative abundance (NM)",
  xlab = "Date",
  main = NULL,
  ...
)
```

## Arguments

- ts_flight_curve_multi:

  pheno_curve object returned by
  [flight_curve_multi](https://retoschmucki.github.io/rbms/reference/flight_curve_multi.md)
  or a data.table with flight curve data.

- SelectYear:

  integer or vector of years to plot, default NULL (plot all years).

- TimeUnit:

  character Time-step used, 'd' for day or 'w' for week, default 'd'.

- weekday:

  Integer for selected day of the week for weekly data, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- color_palette:

  string Name of color palette to use. Options: "rainbow", "heat",
  "terrain", "topo", "cm", "viridis", or a vector of color names/hex
  codes.

- add_legend:

  logical Add legend to plot, default TRUE.

- legend_pos:

  string Position of legend: "topright", "topleft", "bottomright",
  "bottomleft", "top", "bottom", "left", "right", default "topright".

- line_width:

  numeric Width of lines, default 2.

- ymax:

  numeric Maximum value for y-axis, default NULL (auto-scale).

- ylab:

  string Label for y-axis, default "Relative abundance (NM)".

- xlab:

  string Label for x-axis, default "Date".

- main:

  string Main title for plot, default NULL.

- ...:

  Additional parameters passed to base plot function.

## Value

Returns a base plot with relative abundance (y) over time (x), with each
year in a different color.

## See also

[flight_curve_multi](https://retoschmucki.github.io/rbms/reference/flight_curve_multi.md),
[fit_bam_multi](https://retoschmucki.github.io/rbms/reference/fit_bam_multi.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
