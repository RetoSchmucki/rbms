# collated_index compute a collated index from the site indices, using a Generalized Linear Model.

collated_index compute a collated index from the site indices, using a
Generalized Linear Model.

## Usage

``` r
collated_index(
  data,
  s_sp,
  sindex_value = "SINDEX",
  bootID = NULL,
  boot_ind = NULL,
  glm_weights = TRUE,
  rm_zero = TRUE,
  GlmFit = c("glm", "newton", "ipf"),
  SplitGroup = TRUE,
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

  string Name of the response variable to be used the computation of the
  collated index, default is "SINDEX", but could be other standardized
  values.

- bootID:

  integer Identify the n^th bootstrap.

- boot_ind:

  data.table Index of the bootstrap and the site ids of the specific
  bootstrap sample.

- glm_weights:

  logical Use the proportion of the flight curve sampled as weight for
  the collated index, if FALSE the function uses equal weights.

- rm_zero:

  logical Remove the sites where species was not observed to speed-up
  the fit of the GLM without altering the output, default TRUE.

- GlmFit:

  character Fitting of the GLM, 'glm' (default,
  [`glm()`](https://rdrr.io/r/stats/glm.html)), 'newton' or 'ipf', see
  details.

- SplitGroup:

  logical When the years form groups not linked through sites with
  observations, fit each group separately (default TRUE, the correct
  method). FALSE fits all years together as in previous versions (glm
  only), the levels of the groups are then arbitrary; 'newton' and 'ipf'
  always fit the groups separately.

- site_length:

  string Name of a column of `data` with the length of each site
  (transect) in km, constant within site, default NULL. With a length,
  the collated index is the expected count on 1 km of transect (the site
  effects are divided by the length of the site), otherwise the expected
  count per site.

- index_sites:

  character Sites the collated index is averaged over: "all" (default)
  all sites monitored for the species, the sites where it was never
  observed counting as zero (mean over the monitored area), or
  "observed" only the sites where the species was observed at least once
  (mean over the occupied sites).

## Value

a list with the collated indices per year (col_index), the site ids of
the sample (site_id) and, with GlmFit = 'newton', the year effects
(year_effect, see details), NULL otherwise.

## Details

Site-years with NA or infinite `sindex_value` (e.g. years where the
flight curve failed) are treated as missing and removed with a warning,
as done by
[site_index](https://retoschmucki.github.io/rbms/reference/site_index.md);
a site with NA in all years is therefore not part of the sample. When
`boot_ind` was generated for a species with `boot_sample(s_sp = )`,
`s_sp` must be the same species. A bootstrap sample that failed in
[`boot_sample()`](https://retoschmucki.github.io/rbms/reference/boot_sample.md)
(NA row) returns NA for the collated index and the number of sites. With
GlmFit = 'newton' or 'ipf', the same weighted Poisson GLM (site and year
effects) is fitted without building the model matrix (see
[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md)
for the two methods), which is much faster with many sites and gives the
same collated index as [`glm()`](https://rdrr.io/r/stats/glm.html). With
'newton', year_effect gives per year the log year effect relative to the
first year of each group of linked years (LOG_EFFECT), i.e. the log of
the collated index relative to that year, with its quasi-Poisson
standard error (SE), an alternative to the bootstrap that ignores the
variation among sites. The column GROUP gives the group of years linked
through sites with observations in the original data (see
[n_components](https://retoschmucki.github.io/rbms/reference/n_components.md)),
NA for years without observation. When the years form more than one
group, the collated index is only comparable between years of the same
group (warning): with SplitGroup = TRUE, each group is fitted separately
on its own years and sites, its index being the mean over the sites with
data in the years of the group, so each group is a separate series.
Fitting all years together (SplitGroup = FALSE, previous versions) gives
arbitrary levels between groups and can distort the trend within a
group, because zero counts of sites observed only in another group enter
the fit. The data cannot relate the levels of separate groups. Comparing
them rests on the sites of each group being representative of the same
population. This holds by design when sites are assigned to years at
random, e.g. a rotating panel where a random set of sites is surveyed
each year and each site is revisited every few years (each panel then
forms a group), but not when the groups result from site turnover, a
change in the scheme or different regions. In that case, the groups are
separate series. For randomised designs,
[collated_index_mixed](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md)
links the groups through a random site effect. The level of the collated
index depends on the sites it is averaged over (`index_sites`) and on
the length of the transects (`site_length`), its changes between years
do not. When the level is used as a density, e.g. to weight schemes by
their population in a multi-scheme index, see the section "Collated
index as a density" of
[`vignette("vignette", package = "rbms")`](https://retoschmucki.github.io/rbms/articles/vignette.md).

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
