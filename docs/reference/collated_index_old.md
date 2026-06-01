# collated_index_old compute a collated index from the site indices, using a Generalized Linear Model.

collated_index_old compute a collated index from the site indices, using
a Generalized Linear Model.

## Usage

``` r
collated_index_old(site_indices, GlmWeight = NULL, GlmFamily = poisson())
```

## Arguments

- site_indices:

  data.table or data.frame with site indices per year and proportion of
  flight curve covered by the monitoring.

- GlmWeight:

  vector of weight used in the GLM.

- GlmFamily:

  family used for the GLM model.

## Value

a list of three objects, a vector of site, a glm model object, and a
vector of collated indices per year.

## See also

[impute_count](impute_count.md), [flight_curve](flight_curve.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
