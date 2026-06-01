# fit_bam_multi Fit a Generalized Additive Model using bam() with factor-smooth interactions to butterfly count data across multiple years simultaneously. This function fits a single model across multiple years using factor-smooth interactions (fs), which allows year-specific smooth terms while sharing information across years. This is more efficient than fitting separate models per year.

fit_bam_multi Fit a Generalized Additive Model using bam() with
factor-smooth interactions to butterfly count data across multiple years
simultaneously. This function fits a single model across multiple years
using factor-smooth interactions (fs), which allows year-specific smooth
terms while sharing information across years. This is more efficient
than fitting separate models per year.

## Usage

``` r
fit_bam_multi(
  dataset_multi,
  NbrSample = NULL,
  GamFamily = "poisson",
  MaxTrial = 4,
  TimeUnit = "d",
  MultiVisit = "mean",
  weekday = 3,
  tp_col = NULL,
  smooth_basis = "cr",
  k_value = -1,
  verbose = TRUE,
  ...
)
```

## Arguments

- dataset_multi:

  data.table Butterfly counts for species x over multiple years y over
  all sites.

- NbrSample:

  integer Maximum number of sites to sample per year, default=NULL (use
  all sites).

- GamFamily:

  string Family for GAM, default='poisson', but can be 'nb' or
  'quasipoisson'.

- MaxTrial:

  integer Maximum number of fitting attempts, default=4.

- TimeUnit:

  character Time-step for the spline, 'd' for day or 'w' for week.

- MultiVisit:

  string Function for summarising multiple counts within a time unit,
  'max' or 'mean' (default).

- weekday:

  Integer for selected day of the week for weekly summary, default is 3
  (Wednesday). \[1-7\] where 1 = Monday.

- tp_col:

  string Name of temporal variable used in the GAM model, default NULL.

- smooth_basis:

  string Basis for the smooth term, default "cr" (cubic regression
  spline). Other options include "tp" (thin plate), "ps" (P-splines).

- k_value:

  integer Basis dimension for the smooth, default=-1 (automatic
  selection).

- verbose:

  logical Indicating if progress reports should be given.

- ...:

  Additional parameters passed to bam function from the
  [bam](https://rdrr.io/pkg/mgcv/man/bam.html) package.

## Value

A list with three objects: i) \*\*f_curve\*\*: a data.table with the
flight curves for all years with expected relative abundance (NM),
normalized to sum to one per year, ii) \*\*f_model\*\*: the resulting
bam model fitted on the count data across all years, iii)
\*\*f_data\*\*: a data.table with the data used to fit the BAM model.

## Details

This function uses factor-smooth interactions (s(tp_col, M_YEAR,
bs="fs")) which allows: - Fitting a single model for multiple years
(more efficient) - Year-specific smooths that can differ in shape -
Sharing of information across years for more stable estimates - Better
handling of years with sparse data

## See also

[flight_curve_multi](https://retoschmucki.github.io/rbms/reference/flight_curve_multi.md),
[fit_gam](https://retoschmucki.github.io/rbms/reference/fit_gam.md),
[bam](https://rdrr.io/pkg/mgcv/man/bam.html)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
