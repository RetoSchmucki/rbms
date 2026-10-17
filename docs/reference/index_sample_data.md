# index_sample_data Data of a species for a collated index: site-years with NA or infinite `sindex_value` removed (warning), the sites of the original data or of the bootstrap sample `bootID` (each draw of a site with its own uID), the groups of years linked through sites with observations in the original data, and the number of sites per year. Used by [collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md) and [collated_index_mixed](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md).

index_sample_data Data of a species for a collated index: site-years
with NA or infinite `sindex_value` removed (warning), the sites of the
original data or of the bootstrap sample `bootID` (each draw of a site
with its own uID), the groups of years linked through sites with
observations in the original data, and the number of sites per year.
Used by
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)
and
[collated_index_mixed](https://retoschmucki.github.io/rbms/reference/collated_index_mixed.md).

## Usage

``` r
index_sample_data(
  data,
  s_sp,
  sindex_value = "SINDEX",
  bootID = NULL,
  boot_ind = NULL
)
```

## Arguments

- data, s_sp, sindex_value, bootID, boot_ind:

  see
  [collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md).

## Value

list with the species data (data), bootID, the sites of the sample (a:
uID, SITE_ID), the group of each year (grp), the data of the sample
(data_boot) and the number of sites per year (sum_data); or, for a
bootstrap sample that failed in
[boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md)
(NA row), the output of the collated index with NA (failed).

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
