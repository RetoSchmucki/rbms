# year_components Group of linked years for each year, see [n_components](https://retoschmucki.github.io/rbms/reference/n_components.md).

year_components Group of linked years for each year, see
[n_components](https://retoschmucki.github.io/rbms/reference/n_components.md).

## Usage

``` r
year_components(site, year)
```

## Arguments

- site:

  vector Site id of each observation (e.g. count \> 0).

- year:

  vector Year of each observation, same length as `site`.

## Value

integer vector Group label of each year, named by year (sorted as
`levels(factor(year))`).

## See also

[n_components](https://retoschmucki.github.io/rbms/reference/n_components.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
