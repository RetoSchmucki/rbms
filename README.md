# rbms

[![Build Status](https://travis-ci.org/RetoSchmucki/rbms.png?branch=master)](https://travis-ci.org/RetoSchmucki/rbms)

<img style="float: right;" src="rbmshexOR200.png" hspace="20">

With `rbms`, we aim to facilitate the implementation of statistical and mathematical methods developed for computing relative abundance indices from yearly time-series of butterfly counts. These data are characterised by temporal patterns (phenology) that must be accounted for when deriving abundance from a time-series of counts.  As a toolbox, we plan to implement more methods to compute and visualise metrics as they develop. The rbms package will provide the option of being coupled and work in line with other tools available and developed by the community (e.g. [rtrim](https://cran.r-project.org/web/packages/rtrim/), [BRCindicators](https://github.com/BiologicalRecordsCentre/BRCindicators)). With the development of the 'rbms' R package, we also provide a tutorial to facilitate its usage and understanding.

Although `rbms` implements methods that have been developed independently and for which the source should be cited. Users should also be citing the `rbms` package and its version to ensure appropriate referencing improve transparency as well as the reproducibility of the work.

#### Suggested citation for the rbms package

Schmucki R., Harrower C.A.,  Dennis E.B. (2019) rbms: Computing generalised abundance indices for butterfly monitoring count data. R package version 1.0.0. https://github.com/RetoSchmucki/rbms

#### Installation

Once you have successfully installed the [R programming system](https://cran.r-project.org/), you can install the `rbms package` from GitHub, you need to install the package `devtools` available on CRAN. Once installed, use the function `devtools::install_github()` to install the `rbms` package on your system.

> Note that `rbms` is build with R > 3.6.0, so you might need to update your R system before installation.
> To install devtools on Windows system, you will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) - note Rtools40 for R 4.x.x
> New to `R` or want to refresh your R coding skills, try the excellent tutorials available at [ourcodingclub](https://ourcodingclub.github.io/) :thumbsup:

```R
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms")
```

#### Get Started

Consult the [tutorial vignettes / Articles](https://retoschmucki.github.io/rbms/articles/Get_Started_1.html) to learn about `rbms` data format and the workflow how to compute flight curves, abundance indices, and collated index with a bootstrap confidence interval

Further documentation is also available through the help function in R or from the [rbms online references](https://retoschmucki.github.io/rbms/reference/index.html)

#### Reporting Issues

For reporting issues related to this package, please visit the issue and see if your problem has not yet been reported before opening a new [issue here](https://github.com/RetoSchmucki/rbms/issues)
