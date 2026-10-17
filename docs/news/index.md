# Changelog

## rbms 1.3.1

### New features

- [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  gains `parallel` and `n_cores` to fit the flight curves of the years
  in parallel. The random streams of the workers are seeded from the
  session, so the result is reproducible with
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

### Changes

- [`flight_curve_optimised()`](https://retoschmucki.github.io/rbms/reference/flight_curve_optimised.md)
  and
  [`get_nm_optimised()`](https://retoschmucki.github.io/rbms/reference/get_nm_optimised.md)
  are deprecated: they gave the same results as
  [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  and
  [`get_nm()`](https://retoschmucki.github.io/rbms/reference/get_nm.md),
  which they now call with a warning.
- [`flight_curve_multi()`](https://retoschmucki.github.io/rbms/reference/flight_curve_multi.md),
  [`fit_bam_multi()`](https://retoschmucki.github.io/rbms/reference/fit_bam_multi.md)
  and
  [`plot_flight_curve_multi()`](https://retoschmucki.github.io/rbms/reference/plot_flight_curve_multi.md)
  (one `bam()` model for all years) are removed from the package for
  further development; they return an error, the previous code is
  available in rbms 1.3.0.
- [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  checks all the required columns of `ts_season_count` (only the first
  one was checked).
- [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  with `KeepModel = FALSE` and `KeepModelData = TRUE` returns the flight
  curves and the data (it returned the unformatted list).

### Documentation

- Vignette: the installation code is shown but not run, the messages of
  [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  are not printed, `linewidth` replaces the deprecated `size` of ggplot2
  lines; author and date updated.
- LICENSE year 2022-2026.

### Potential breaking changes

- Code using
  [`flight_curve_multi()`](https://retoschmucki.github.io/rbms/reference/flight_curve_multi.md),
  [`fit_bam_multi()`](https://retoschmucki.github.io/rbms/reference/fit_bam_multi.md)
  or
  [`plot_flight_curve_multi()`](https://retoschmucki.github.io/rbms/reference/plot_flight_curve_multi.md)
  (exported in 1.3.0) needs rbms 1.3.0, or
  [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md).
- [`flight_curve()`](https://retoschmucki.github.io/rbms/reference/flight_curve.md)
  now stops with a clear message when a required column of
  `ts_season_count` is missing.

## rbms 1.3.0

### New features

- [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md)
  gains `GlmFit = c("glm", "newton", "ipf")`. With `"newton"` or
  `"ipf"`, the same weighted Poisson GLM (site and year effects) is
  fitted without building the model matrix and gives the same collated
  index as [`glm()`](https://rdrr.io/r/stats/glm.html), much faster with
  many sites (e.g. 0.2 s instead of 3 minutes with 2000 sites and 20
  years). With `"newton"`, the year effects and their standard error are
  returned (`year_effect`).
- [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md)
  gains `site_length` (index per 1 km of transect) and
  `index_sites = c("all", "observed")` (index averaged over all
  monitored sites or over the sites where the species was observed).
- New
  [`collated_index_mixed()`](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md):
  the collated index from a mixed model with a random site effect
  (package glmmTMB, Poisson or negative binomial), with an optional
  habitat (stratum) effect and habitat weights (by sites, equal or user
  defined, e.g. area), `site_length` as offset and `index_sites`. It
  links years that are not linked through common sites, for randomised
  designs such as rotating panels.
- [`boot_sample()`](https://retoschmucki.github.io/rbms/reference/boot_sample.md)
  gains `s_sp` to draw the bootstrap samples per species, stratified by
  the occurrence of zeros at the sites and by group of linked years,
  with `habitat` as additional stratum. A sample with fewer than
  `min_year` years or where a group of linked years is split or lost is
  drawn again, up to `max_trial` times; failed samples are returned as
  NA rows. Without `s_sp`, the samples are drawn as before.
- [`impute_count()`](https://retoschmucki.github.io/rbms/reference/impute_count.md)
  gains `SindexMethod = c("ratio", "glm")`. The default ratio is
  unchanged; with `"glm"`, the counts of the time steps not visited are
  predicted from a Poisson GLM with site and year effects and the flight
  curve as offset (`GlmFit = c("newton", "ipf")`), which avoids extreme
  site indices at low coverage. The year effects are returned as
  attribute `glm_year_effect`.

### Changes

- [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md):
  when the years form groups not linked through sites with observations,
  each group is fitted separately (`SplitGroup = TRUE`, default) and
  identified in the new column `GROUP`; fitting all years together gave
  arbitrary levels between groups. `SplitGroup = FALSE` reproduces the
  previous fit.
- Site-years with NA or infinite site index are treated as missing:
  [`impute_count()`](https://retoschmucki.github.io/rbms/reference/impute_count.md)
  returns NA instead of `Inf` or `NaN` when the visits cover no part of
  the flight curve (`TOTAL_NM = 0`), and
  [`site_index()`](https://retoschmucki.github.io/rbms/reference/site_index.md),
  [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md),
  [`collated_index_mixed()`](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md)
  and
  [`boot_sample()`](https://retoschmucki.github.io/rbms/reference/boot_sample.md)
  remove them (previously a site with NA in all years set the collated
  index to 0 in all years).
- A bootstrap sample that failed in
  [`boot_sample()`](https://retoschmucki.github.io/rbms/reference/boot_sample.md)
  returns NA in
  [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md);
  `boot_ind` drawn for another species is an error.
- [`collated_index_old()`](https://retoschmucki.github.io/rbms/reference/collated_index_old.md)
  is deprecated and returns an error, use
  [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`points()`](https://rdrr.io/r/graphics/points.html) methods for
  `pheno_curve` objects: the first argument is now `x`, as in the
  generics.
- speedglm is no longer needed; glmmTMB and ggplot2 are suggested.

### Potential breaking changes

With the default arguments, the results of the standard workflow (flight
curve, imputed counts, site indices, collated index and bootstrap) are
identical to version 1.2.1. A pipeline is only affected in the
following, unlikely cases:

- [`collated_index()`](https://retoschmucki.github.io/rbms/reference/collated_index.md)
  returns an additional column `GROUP` in `col_index` and an additional
  element `year_effect`; code that selects columns by position or binds
  outputs with [`rbind()`](https://rdrr.io/r/base/cbind.html) may need
  adjusting (selection by name and
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  are not affected).
- Data where the years form groups not linked through common sites: the
  collated index now differs, each group being fitted separately
  (`SplitGroup = FALSE` gives the previous result, whose levels between
  groups were arbitrary).
- Data with NA or infinite site indices: these site-years are now
  removed; previously an NA site set the collated index to 0 in all
  years, and `site_index(MinFC = NULL)` kept infinite values.
- speedglm is no longer installed with rbms; scripts calling
  `speedglm()` directly need it installed.
- [`collated_index_old()`](https://retoschmucki.github.io/rbms/reference/collated_index_old.md)
  now returns an error (it did not run in previous versions).

### Documentation

- Vignette: per-species stratified bootstrap, linked years and groups of
  years, when separate groups can be compared, the collated index as a
  density (per km, sites used) and references with DOIs.
