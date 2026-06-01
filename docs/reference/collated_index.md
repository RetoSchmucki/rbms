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
  rm_zero = TRUE
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

## Value

a list of two objects, a vector of site, a glm model object, and a
vector of collated indices per year.

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[flight_curve](https://retoschmucki.github.io/rbms/reference/flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
