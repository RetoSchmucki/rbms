# plot.pheno_curve Generic method to plot the flight curve, where values are extracted from a pheno_curve object (outcome of the light_curve() function)

plot.pheno_curve Generic method to plot the flight curve, where values
are extracted from a pheno_curve object (outcome of the light_curve()
function)

## Usage

``` r
# S3 method for class 'pheno_curve'
plot(data_x, year = NULL, weekday = 3, SiteID = NULL, ymax = NULL, ...)
```

## Arguments

- data_x:

  pheno_curve object which is the outcome of the light_curve() function

- year:

  integer or vector for the year(s) to be displayed (e.g. 2015 or
  c(2000, 2001)), default is NULL.

- weekday:

  weekday to be used for date in weekly count, where 1 refers to Monday,
  default is 3 (Wednesday)

- SiteID:

  integer or string to ID the site to display, default the first site is
  displayed.

- ymax:

  maximum value for y axis

- ...:

  additional parameters for base plot.

## Value

Returns a base plot with relative abundance (y) over time (x)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
