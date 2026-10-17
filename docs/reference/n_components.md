# n_components Number of groups of years linked through the sites with an observation. Two years are linked when a site has an observation in both, directly or through other years and sites. With one group, the level of all years can be related in the collated index model; with more, the site and year effects of the groups cannot be related and the collated index of the bootstrap sample can run away. Years without any observation are not counted.

n_components Number of groups of years linked through the sites with an
observation. Two years are linked when a site has an observation in
both, directly or through other years and sites. With one group, the
level of all years can be related in the collated index model; with
more, the site and year effects of the groups cannot be related and the
collated index of the bootstrap sample can run away. Years without any
observation are not counted.

## Usage

``` r
n_components(site, year)
```

## Arguments

- site:

  vector Site id of each observation (e.g. SINDEX \> 0).

- year:

  vector Year of each observation, same length as `site`.

## Value

integer The number of groups of linked years, `Inf` if there is no
observation.

## See also

[boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md),
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
