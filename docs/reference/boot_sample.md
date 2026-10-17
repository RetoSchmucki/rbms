# boot_sample Generate n bootstrap sample of the monitoring sites to be used for each iteration.

boot_sample Generate n bootstrap sample of the monitoring sites to be
used for each iteration.

## Usage

``` r
boot_sample(
  data,
  boot_n = 1000,
  s_sp = NULL,
  sindex_value = "SINDEX",
  strata = TRUE,
  check_link = TRUE,
  min_year = 3,
  max_trial = 5,
  habitat = NULL
)
```

## Arguments

- data:

  data.table or data.frame Data with all site id, with SPECIES, M_YEAR
  and `sindex_value` when `s_sp` is given.

- boot_n:

  integer The number of bootstrap samples to be generated.

- s_sp:

  string Species for which the bootstrap samples are drawn, default NULL
  draws among all sites without checks.

- sindex_value:

  string Name of the site index variable used to define strata and
  observations, default is "SINDEX".

- strata:

  logical Stratify the sampling by the occurrence of zeros at the sites
  and by group of linked years, default TRUE.

- check_link:

  logical Draw again the samples where a group of linked years is split
  or lost, default TRUE.

- min_year:

  integer Minimum number of years in a bootstrap sample, default 3.

- max_trial:

  integer Maximum number of trials to draw a valid bootstrap sample,
  default 5.

- habitat:

  string Name of a column of `data` with the habitat (or stratum) of
  each site, constant within site, default NULL. The habitat is crossed
  with the strata, so each sample keeps the number of sites of each
  habitat (design stratified by habitat, see
  [collated_index_mixed](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md)).

## Value

A list with site id and bootstrap indices for n bootstrap sample. When
`s_sp` is given, the list also contains the species (s_sp), the number
of groups of years of the original data (n_comp_0), the site indices per
stratum (strata_ind, named group_zero pattern), the group of each year
(year_group) and of each site (site_group), the number of trials per
sample (n_trial) and the number of failed samples (n_failed).

## Details

When `s_sp` is given, the samples are drawn for that species only,
stratified by the occurrence of zeros at the sites (never zero, some
zeros, always zero), and a sample is drawn again, up to `max_trial`
trials, when it has fewer than `min_year` years or when its years are
not linked through the sites with observations (see `n_components`). A
sample still failing after `max_trial` trials is returned as a row of
NA. When the years of the original data form several groups not linked
through sites with observations (see
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)),
the strata are crossed with the groups (a site with observations belongs
to one group, a site never observed to the group of most of its years),
so each sample keeps the number of sites of each group, and each group
must remain linked in the sample (not split or lost). With one group,
this is the stratification by zeros only. For reproducible samples, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
`boot_sample()`; when looping over species, setting the seed before each
species makes its samples independent of the order or the subset of
species processed.

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
