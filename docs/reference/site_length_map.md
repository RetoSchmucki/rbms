# site_length_map Length of each site (transect) from a column of the data, checked to be positive, not NA and constant within site.

site_length_map Length of each site (transect) from a column of the
data, checked to be positive, not NA and constant within site.

## Usage

``` r
site_length_map(data, site_length)
```

## Arguments

- data:

  data.table Data of a species with SITE_ID.

- site_length:

  string Name of the column with the length of the site (km).

## Value

data.table with SITE_ID and SITE_LENGTH.

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
