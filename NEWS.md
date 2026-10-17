# rbms 1.3.0


## New features

* `collated_index()` gains `GlmFit = c("glm", "newton", "ipf")`. With `"newton"` or `"ipf"`, the same weighted Poisson GLM
  (site and year effects) is fitted without building the model matrix and gives the same collated index as `glm()`, much
  faster with many sites (e.g. 0.2 s instead of 3 minutes with 2000 sites and 20 years). With `"newton"`, the year effects
  and their standard error are returned (`year_effect`).
* `collated_index()` gains `site_length` (index per 1 km of transect) and `index_sites = c("all", "observed")` (index
  averaged over all monitored sites or over the sites where the species was observed).
* New `collated_index_mixed()`: the collated index from a mixed model with a random site effect (package glmmTMB,
  Poisson or negative binomial), with an optional habitat (stratum) effect and habitat weights (by sites, equal or user
  defined, e.g. area), `site_length` as offset and `index_sites`. It links years that are not linked through common sites,
  for randomised designs such as rotating panels.
* `boot_sample()` gains `s_sp` to draw the bootstrap samples per species, stratified by the occurrence of zeros at the
  sites and by group of linked years, with `habitat` as additional stratum. A sample with fewer than `min_year` years or
  where a group of linked years is split or lost is drawn again, up to `max_trial` times; failed samples are returned as
  NA rows. Without `s_sp`, the samples are drawn as before.
* `impute_count()` gains `SindexMethod = c("ratio", "glm")`. The default ratio is unchanged; with `"glm"`, the counts of the
  time steps not visited are predicted from a Poisson GLM with site and year effects and the flight curve as offset
  (`GlmFit = c("newton", "ipf")`), which avoids extreme site indices at low coverage. The year effects are returned as
  attribute `glm_year_effect`.

## Changes

* `collated_index()`: when the years form groups not linked through sites with observations, each group is fitted
  separately (`SplitGroup = TRUE`, default) and identified in the new column `GROUP`; fitting all years together gave
  arbitrary levels between groups. `SplitGroup = FALSE` reproduces the previous fit.
* Site-years with NA or infinite site index are treated as missing: `impute_count()` returns NA instead of `Inf` or `NaN`
  when the visits cover no part of the flight curve (`TOTAL_NM = 0`), and `site_index()`, `collated_index()`,
  `collated_index_mixed()` and `boot_sample()` remove them (previously a site with NA in all years set the collated index
  to 0 in all years).
* A bootstrap sample that failed in `boot_sample()` returns NA in `collated_index()`; `boot_ind` drawn for another species
  is an error.
* `collated_index_old()` is deprecated and returns an error, use `collated_index()`.
* `plot()` and `points()` methods for `pheno_curve` objects: the first argument is now `x`, as in the generics.
* speedglm is no longer needed; glmmTMB and ggplot2 are suggested.

## Potential breaking changes

With the default arguments, the results of the standard workflow (flight curve, imputed counts, site indices, collated
index and bootstrap) are identical to version 1.2.1. A pipeline is only affected in the following, unlikely cases:

* `collated_index()` returns an additional column `GROUP` in `col_index` and an additional element `year_effect`; code that
  selects columns by position or binds outputs with `rbind()` may need adjusting (selection by name and
  `data.table::rbindlist()` are not affected).
* Data where the years form groups not linked through common sites: the collated index now differs, each group being
  fitted separately (`SplitGroup = FALSE` gives the previous result, whose levels between groups were arbitrary).
* Data with NA or infinite site indices: these site-years are now removed; previously an NA site set the collated index
  to 0 in all years, and `site_index(MinFC = NULL)` kept infinite values.
* speedglm is no longer installed with rbms; scripts calling `speedglm()` directly need it installed.
* `collated_index_old()` now returns an error (it did not run in previous versions).

## Documentation

* Vignette: per-species stratified bootstrap, linked years and groups of years, when separate groups can be compared,
  the collated index as a density (per km, sites used) and references with DOIs.
