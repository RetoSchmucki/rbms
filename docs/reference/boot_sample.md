# boot_sample Generate n bootstrap sample of the monitoring sites to be used for each iteration.

boot_sample Generate n bootstrap sample of the monitoring sites to be
used for each iteration.

## Usage

``` r
boot_sample(data, boot_n = 1000)
```

## Arguments

- data:

  data.table or data.frame Data with all site id.

- boot_n:

  integer The number of bootstrap samples to be generated.

## Value

A list with site id and bootstrap indices for n bootstrap sample.

## See also

[impute_count](https://retoschmucki.github.io/rbms/reference/impute_count.md),
[collated_index](https://retoschmucki.github.io/rbms/reference/collated_index.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
