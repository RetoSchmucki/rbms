# check_package Internal function to check if a package is installed.

check_package Internal function to check if a package is installed.

## Usage

``` r
check_package(
  pkgName = NULL,
  message1 = "You need to install ",
  message2 = "This version requires"
)
```

## Arguments

- pkgName:

  The package name to be verified.

- message1:

  Inform the user about the package dependency.

- message2:

  Inform the user what happen if the package is not installed.

## Value

If package is not found, the user is offered the option to install the
missing package.

## Author

Reto Schmucki - <retoshm@ceh.ac.uk>
