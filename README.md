# rbms

[![Build Status](https://travis-ci.org/RetoSchmucki/rbms.png?branch=master)](https://travis-ci.org/RetoSchmucki/rbms)

<img style="float: right;" src="rbmshexOR200.png" hspace="20">

With `rbms`, we aim to facilitate the implementation of statistical and mathematical methods developed for computing relative abundance indices from yearly time-series of butterfly counts. These data are characterised by temporal patterns (phenology) that must be accounted for when deriving abundance from a time-series of counts.  As a toolbox, we plan to implement more methods to compute and visualise metrics as they develop. The rbms package will provide the option of being coupled and work in line with other tools available and developed by the community (e.g. [rtrim](https://cran.r-project.org/web/packages/rtrim/), [BRCindicators](https://github.com/BiologicalRecordsCentre/BRCindicators)). With the development of the 'rbms' R package, we also provide a tutorial to facilitate its usage and understanding.

Although `rbms` implements methods that were developed independently and for which the source should be cited, users should also cite the `rbms` package and its version to ensure appropriate referencing and improve the reproducibility of the work.

#### Suggested citation for the rbms package

Schmucki R., Harrower C.A.,  Dennis E.B. (2020) rbms: Computing generalised abundance indices for butterfly monitoring count data. R package version 1.0.2. https://github.com/RetoSchmucki/rbms

#### The rbms package implements methods from:

- Dennis, E.B., Freeman, S.N., Brereton, T., Roy, D.B., 2013. Indexing butterfly abundance whilst accounting for missing counts and variability in seasonal pattern. Methods Ecol Evol 4, 637–645. https://doi.org/10.1111/2041-210X.12053
- Dennis, E.B., Morgan, B.J.T., Freeman, S.N., Brereton, T.M., Roy, D.B., 2016. A generalized abundance index for seasonal invertebrates. Biom 72, 1305–1314. https://doi.org/10.1111/biom.12506
- Schmucki, R., Pe’er, G., Roy, D.B., Stefanescu, C., Van Swaay, C.A.M., Oliver, T.H., Kuussaari, M., Van Strien, A.J., Ries, L., Settele, J., Musche, M., Carnicer, J., Schweiger, O., Brereton, T.M., Harpke, A., Heliölä, J., Kühn, E., Julliard, R., 2016. A regionally informed abundance index for supporting integrative analyses across butterfly monitoring schemes. J Appl Ecol 53, 501–510. https://doi.org/10.1111/1365-2664.12561

#### Installation

Once the [R programming system](https://cran.r-project.org/) is successfully installed, you can install the `rbms package` from GitHub. For this you will need to install the package `devtools` or `remotes`, both available on CRAN. Once installed, use the function `devtools::install_github()` to install the `rbms` package on your system.

> Note that `rbms` is build with R > 3.6.0, so you might need to update your R system before installation.
> To install devtools on Windows system, you will also need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) - note Rtools40 for R 4.x.x
> New to `R` or want to refresh your R coding skills, try the excellent tutorials available at [ourcodingclub](https://ourcodingclub.github.io/) :thumbsup:

```R
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms")
```

#### Get Started

Consult the [tutorial vignettes / Articles](https://retoschmucki.github.io/rbms/articles/Get_Started_1.html) to learn about `rbms` data format and the workflow how to compute flight curves, abundance indices, and collated index with a bootstrap confidence interval

Further documentation is also available through the help function in R or from the [rbms online references](https://retoschmucki.github.io/rbms/reference/index.html)

### Development and the beta version

The development version, rbms v.1.2.0, address an issue that was affecting some datasets and resulted in not computing flight curves and indices, although enough visits, sites and occurrence were recorded.
This version can be installed from the specific issue branch

```R
devtools::install_github("RetoSchmucki/rbms", ref = "dev_issue_no8")
```

The rbms v.1.2.0 (dev) work like the previous one but should have solved the issue noted previously, resulting in potentially more flight curve and indices being computed.
Version 1.2.0 implements a basic plot method for "pheno_curve" object that is produced by the `flight_curve()` function, use:

 ```R
 # usage or plot method
 ts_flight_curve <- flight_curve(ts_season_count, ...)
 plot(ts_flight_curve)
 points(ts_flight_curve, col = 'magenta', pch = 19)

 # for a single year
 plot(ts_flight_curve, year = 2017)
 points(ts_flight_curve, year = 2017, col = 'magenta', pch = 19)
 ```
#### Reporting Issues

For reporting issues related to this package or workflow, please visit the issue and before opening a new, see if your problem is not already in the list reported [issues here](https://github.com/RetoSchmucki/rbms/issues)
