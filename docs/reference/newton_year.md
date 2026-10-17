# newton_year Newton steps on the concentrated likelihood of the year effects, see [site_year_level](https://retoschmucki.github.io/rbms/reference/site_year_level.md).

newton_year Newton steps on the concentrated likelihood of the year
effects, see
[site_year_level](https://retoschmucki.github.io/rbms/reference/site_year_level.md).

## Usage

``` r
newton_year(sy, C_s, C_y, n_s, n_y, tol = 1e-08, maxit = 100)
```

## Arguments

- sy:

  data.table Site-year table with site (S_IDX), year (Y_IDX), count
  (C_SY) and flight curve covered (N_SY).

- C_s, C_y:

  numeric Total count per site and per year.

- n_s, n_y:

  integer Number of sites and years.

- tol:

  numeric Tolerance on the gradient, relative to the largest year total.

- maxit:

  integer Maximum number of Newton steps.

## Value

list with the year effects b (exp of the log effects, 0 for years
without count) and a data.table with GROUP, LOG_EFFECT and SE (not
scaled by the dispersion) per year.

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
