# draw_strata Draw one stratified bootstrap sample of site indices, sampling with replacement within each stratum the number of sites in that stratum. Uses `sample.int()` on positions so that a stratum with a single site is not expanded to `1:x` as `sample(x)` would do.

draw_strata Draw one stratified bootstrap sample of site indices,
sampling with replacement within each stratum the number of sites in
that stratum. Uses [`sample.int()`](https://rdrr.io/r/base/sample.html)
on positions so that a stratum with a single site is not expanded to
`1:x` as `sample(x)` would do.

## Usage

``` r
draw_strata(strata_ind)
```

## Arguments

- strata_ind:

  list Integer vectors of site indices, one per stratum.

## Value

integer vector The site indices of the bootstrap sample.

## See also

[boot_sample](https://retoschmucki.github.io/rbms/reference/boot_sample.md)

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
