# fit_gam fit a Generalized Additive Model to butterfly count data along a temporal variable and accounting for site effect when multiple are available.

fit_gam fit a Generalized Additive Model to butterfly count data along a
temporal variable and accounting for site effect when multiple are
available.

## Usage

``` r
fit_gam(
  dataset_y,
  NbrSample = NULL,
  GamFamily = "poisson",
  MaxTrial = 4,
  SpeedGam = TRUE,
  OptiGam = TRUE,
  TimeUnit = "d",
  MultiVisit = "mean",
  weekday = 3,
  mod_form = NULL,
  tp_col = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- dataset_y:

  data.table Filtered butterfly counts for species x over year y over
  all sites.

- NbrSample:

  integer Inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  default=100.

- GamFamily:

  string Inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  default='poisson', but can be 'nb' or 'quasipoisson'.

- MaxTrial:

  integer Inherited from
  [flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
  default=3.

- SpeedGam:

  logical to use the [bam](https://rdrr.io/pkg/mgcv/man/bam.html) method
  instead of the [gam](https://rdrr.io/pkg/mgcv/man/gam.html) method.

- OptiGam:

  logical Set if the [bam](https://rdrr.io/pkg/mgcv/man/bam.html) method
  should be used, default instead of the default
  [gam](https://rdrr.io/pkg/mgcv/man/gam.html) method.

- TimeUnit:

  character The time-step for which the spline should be computed, 'd'
  day or 'w' week.

- MultiVisit:

  string Function to apply for summarising multiple counts within a time
  unit, 'max' or 'mean' (default).

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- mod_form:

  string with formula to be passed to the gam model, default null.

- tp_col:

  string or vector of string with additional variable used in the gam
  model, default null.

- verbose:

  a logical indicating if some “progress report” should be given.

- ...:

  Additional parameters passed to gam or bam function from the
  [gam](https://rdrr.io/pkg/mgcv/man/gam.html) package.

## Value

A list with three objects, i) \*\*f_curve\*\*: a data.table with the
flight curve `f_curve` with expected relative abundance, normalize to
sum to one over a full season, ii) \*\*f_model\*\*: the resulting gam
model `f_model` fitted on the count data and iii) \*\*f_data\*\*: a
data.table with the data used to fit the GAM model. This is provide for
one year 'y'.

## See also

[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md),
[gam](https://rdrr.io/pkg/mgcv/man/gam.html),
[bam](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
