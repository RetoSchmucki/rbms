# Toy data set with butterfly count for x species across y sites

Toy data set with butterfly count for x species across y sites

## Usage

``` r
m_count
```

## Format

data.table object with butterfly count in long format

- SITE_ID:

  number used as site id, the can be and integer or a string

- DATE:

  date when the observation/count was recorded, format YEAR-MM-DD

- SPECIES:

  number used as species id, this can be an integer or a string (e.g.
  "Aglais io")

- DAY:

  integer, day within the month

- MONTH:

  integer, month within the year

- YEAR:

  integer, 4 digit year

- COUNT:

  integer, actual butterfly count
