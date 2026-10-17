# collated_index_mixed compute a collated index from the site indices, using a mixed model with a random site effect.

collated_index_mixed compute a collated index from the site indices,
using a mixed model with a random site effect.

## Usage

``` r
collated_index_mixed(
  data,
  s_sp,
  sindex_value = "SINDEX",
  bootID = NULL,
  boot_ind = NULL,
  glm_weights = TRUE,
  rm_zero = TRUE,
  family = c("poisson", "nbinom2"),
  habitat = NULL,
  habitat_weights = NULL,
  site_length = NULL,
  index_sites = c("all", "observed")
)
```

## Arguments

- data:

  data.table or data.frame Site indices per year and proportion of
  flight curve covered by the monitoring.

- s_sp:

  string Species name to be found in the data.

- sindex_value:

  string Name of the response variable, default is "SINDEX".

- bootID:

  integer Identify the n^th bootstrap.

- boot_ind:

  list Bootstrap samples from
  [boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md).

- glm_weights:

  logical Use the proportion of the flight curve sampled as weight, if
  FALSE the function uses equal weights.

- rm_zero:

  logical Remove the sites where the species was not observed, default
  TRUE (they contribute 0 to the index, but are not used to estimate the
  distribution of the site effects).

- family:

  character Distribution of the counts, 'poisson' (default) or 'nbinom2'
  (negative binomial). With 'poisson', the index matches
  [collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)
  for linked years and its level (mean count per site) is not biased by
  overdispersion; 'nbinom2' weights the counts differently, which can
  give a more precise trend but a lower level of the index.

- habitat:

  string Name of a column of `data` with the habitat (or stratum) of
  each site, constant within site, default NULL (no habitat). The
  habitat is a fixed effect of the model, the site effects being the
  deviation of each site from the mean of its habitat (sites
  representative of their habitat, e.g. sampled at random within habitat
  strata).

- habitat_weights:

  Weights of the habitats in the collated index: NULL (default)
  proportional to the number of sites of each habitat (the mean over all
  sites, as in
  [collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)),
  "equal" (each habitat counts equally), or a numeric vector named by
  habitat, e.g. the share of each habitat in the area of the region for
  a design stratified by habitat (normalised to sum 1).

- site_length:

  string Name of a column of `data` with the length of each site
  (transect) in km, constant within site, default NULL. With a length,
  log(length) is an offset of the model and the collated index is the
  expected count on 1 km of transect, otherwise the expected count per
  site.

- index_sites:

  character Sites the collated index is averaged over: "all" (default)
  all sites monitored for the species, the sites where it was never
  observed counting as zero, or "observed" only the sites where the
  species was observed at least once (habitats without such site are
  then left out of the index and the weights of the other habitats
  normalised).

## Value

a list with the collated indices per year (col_index, same columns as
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)),
the site ids of the sample (site_id), the year effects (year_effect: log
effect relative to the first year, LOG_EFFECT, and its standard error,
SE), the model parameters (model: family, standard deviation of the site
effects, size of the negative binomial, convergence) and, with
`habitat`, the index of each habitat (habitat_index: number of sites,
weight and index per habitat and year; the collated index is the
weighted sum of the habitat indices).

## Details

The model is the GLM of
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)
with the site effect as a random effect, round(SINDEX) ~ 0 +
factor(M_YEAR) + (1 \| site), weighted by the proportion of the flight
curve sampled (TOTAL_NM), fitted with
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
(package glmmTMB required). The collated index is defined as in
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md),
the expected index summed over the sites of the sample and divided by
the number of sites of the species, exp(year effect) x sum(exp(site
effects)) / number of sites, using the predicted site effects. The
random site effect assumes that the sites are a sample of one population
of sites. It links years that are not linked through sites with
observations (column GROUP), which is valid when the sites of all groups
are representative of the same population, e.g. sites assigned to years
at random in a rotating panel. With linked years, the index is close to
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md),
the site effects being shrunk towards their mean. The bootstrap is used
as in
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md),
with `bootID` and `boot_ind` from
[boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md).
When the index is used as a density, see the section "Collated index as
a density" of
[`vignette("vignette", package = "rbms")`](https://retoschmucki.github.io/rbms/articles/vignette.md).

## See also

[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md),
[boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
