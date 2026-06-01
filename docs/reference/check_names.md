# check_names Verify for the required column names in the data object.

check_names Verify for the required column names in the data object.

## Usage

``` r
check_names(x, y)
```

## Arguments

- x:

  Data object in which the variables names are verified.

- y:

  Variable names to be verified.

## Value

Verifyies if column names listed in `y` are found in the data set `x`,
if not, a message identifies the missing column name and stops.

## Details

This function is not case sensitive, but it does not accept different
names or spelling.

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
