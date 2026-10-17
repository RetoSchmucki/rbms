# site_year_level Expected count per unit of flight curve (level) for each site and year, from the Poisson GLM COUNT ~ factor(SITE_ID) + factor(M_YEAR) + offset(log(NM)), fitted without building the model matrix (no limit on the number of sites). The counts are first summed per site-year, with the flight curve covered (sum of NM). With method = "newton", the site effects are concentrated out of the likelihood (closed form given the year effects) and the concentrated likelihood, a multinomial per site, is maximised by Newton steps on the year effects only. One year per group of linked years (see [n_components](https://retoschmucki.github.io/rbms/reference/n_components.md)) is the reference and years without count have a level of 0. If Newton fails, iterative proportional fitting is used. With method = "ipf", iterative proportional fitting matches the fitted total of each site and each year to the observed total. Both give the maximum likelihood estimates of glm().

site_year_level Expected count per unit of flight curve (level) for each
site and year, from the Poisson GLM COUNT ~ factor(SITE_ID) +
factor(M_YEAR) + offset(log(NM)), fitted without building the model
matrix (no limit on the number of sites). The counts are first summed
per site-year, with the flight curve covered (sum of NM). With method =
"newton", the site effects are concentrated out of the likelihood
(closed form given the year effects) and the concentrated likelihood, a
multinomial per site, is maximised by Newton steps on the year effects
only. One year per group of linked years (see
[n_components](https://retoschmucki.github.io/rbms/reference/n_components.md))
is the reference and years without count have a level of 0. If Newton
fails, iterative proportional fitting is used. With method = "ipf",
iterative proportional fitting matches the fitted total of each site and
each year to the observed total. Both give the maximum likelihood
estimates of glm().

## Usage

``` r
site_year_level(
  site,
  year,
  count,
  nm,
  method = c("newton", "ipf"),
  tol = 1e-08,
  maxit = NULL
)
```

## Arguments

- site:

  vector Site id of each count.

- year:

  vector Year of each count.

- count:

  numeric Observed counts.

- nm:

  numeric Proportion of the flight curve (NM) of each count, \> 0.

- method:

  character Fitting method, "newton" (default) or "ipf".

- tol:

  numeric Tolerance for convergence, default 1e-8.

- maxit:

  integer Maximum number of iterations, default 100 for Newton and 10000
  for IPF.

## Value

numeric vector The level of the site-year of each count, the fitted
count being `level * nm`. The attribute "effects" is a list with the
site (site) and year (year) multiplicative effects, named by site and
year, level = site x year (the scale between them is arbitrary, their
product is not). With Newton, the attribute "year_effect" is a
data.table with, per year, the group of linked years (GROUP), the log
year effect relative to the first year of the group (LOG_EFFECT, -Inf
for years without count), its standard error (SE) scaled by the
quasi-Poisson dispersion and conditional on the site totals, and the
dispersion (DISPERSION).

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
