# rbms

[![Build Status](https://travis-ci.org/RetoSchmucki/rbms.png?branch=master)](https://travis-ci.org/RetoSchmucki/rbms)

<img style="float: right;" src="rbmshexOR200.png" hspace="20">

With `rbms`, our aim is to facilitate the implementation of statistical and mathematical methods developed for computing relative abundance indices from yearly time-series of butterfly counts. These data are characterized by a strong phenology that must be accounted for when deriving abundance from counts.  As a toolbox, we plan to implement more methods to compute and visualize metrics as they develop. The rbms package will provide the option of being coupled and work in line with other tools available and developed by the community (e.g. [rtrim](https://cran.r-project.org/web/packages/rtrim/), [BRCindicators](https://github.com/BiologicalRecordsCentre/BRCindicators)). Together with the development of the 'rbms' R package, we also aim to provide tutorial to facilitate its use and understanding.

Although `rbms` implements methods that have been develop independently and for which the original source should be cited, user should also citing the `rbms` package and it version to ensure appropriate referencing and therefore improve transparency and repeatability of the work.

#### Suggested citation for the rbms package

Schmucki R., Harrower C. Dennis E. (2019) rbms: Computing generalised abundance indices for butterfly monitoring count data. R package version 1.0.0. https://github.com/RetoSchmucki/rbms


#### Installation

To install this package from GitHub, you need to install the package `devtools` available on CRAN. Once installed, use the function `devtools::install_github()` to install the `rbms` package on your system.

> Note that `rbms` has been build with R 3.6.0, so you might have to update your R system before installation.

```R
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms")
```

#### Reporting Issues

For reporting issues related to this package, please visit the issue and see if your problem has not yet been reported before opening a new [issue here](https://github.com/RetoSchmucki/rbms/issues)
