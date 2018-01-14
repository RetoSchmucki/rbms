# rbms

### - A new package that replace the former 'regionalGAM' -

With the rapid expansion of monitoring efforts and the usefulness of conducting integrative analyses to inform conservation initiatives, the choice of a robust abundance index is crucial to adequately assess the species status. Butterfly Monitoring Schemes (BMS) operate in increasing number of countries with broadly the same methodology, yet they differ in their observation frequencies and often in the method used to compute annual abundance indices.

Here we implemented the method for computing an abundance index with the *regional GAM* approach, an extension of the two-stages model introduced by [Dennis et al. (2013)](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12053/abstract). This index offers the best performance for a wide range of sampling frequency, providing greater robustness and unbiased estimates compared to the popular linear interpolation approach [(Schmucki et al. 2015)](http://onlinelibrary.wiley.com/doi/10.1111/1365-2664.12561/abstract).

#### Suggested citation

Schmucki R., Pe’er G., Roy D.B., Stefanescu C., Van Swaay C.A.M., Oliver T.H., Kuussaari M., Van Strien A.J., Ries L., Settele J., Musche M., Carnicer J., Schweiger O., Brereton T. M., Harpke A., Heliölä J., Kühn E. & Julliard R. (2016) A regionally informed abundance index for supporting integrative analyses across butterfly monitoring schemes. Journal of Applied Ecology. Vol. 53 (2) 501–510. DOI: 10.1111/1365-2664.12561


#### Installation

To install this package from GitHub, you will fist need to install the package `devtools` that is available from CRAN. From there, simply use the the function `install_github()` to install the `rbms` package on your system. Note that this package was build with R 3.4.2, so you might have to update your R installation. If you are unable to install this package, you might consider sourcing the R script that can be found here: [rbms source code] (https://github.com/RetoSchmucki/rbms/blob/master/R/)

```R
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms")
```

#### LINUX note:
You might have to install to following library to install the `sf package`.

```
sudo apt-get install libudunits2-dev
```

For reporting errors and issues related to this package and its functions, please open a [issue here](https://github.com/RetoSchmucki/rbms/issues)
