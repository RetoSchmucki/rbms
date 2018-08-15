# rbms

[![Build Status](https://travis-ci.org/RetoSchmucki/rbms.png?branch=master)](https://travis-ci.org/RetoSchmucki/rbms)


### - A new package that replaces the former 'regionalGAM' -

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


> LINUX note:
> You might have to install to following library to install the `sf package`.
>```
>sudo apt-get install libudunits2-dev
>```

#### Usage example
##### 1. load package and data included in the package

```r
library(rbms)
data(m_visit)
data(m_count)
```

##### 2. organize the data to cover the time period and monitoring season of the BMS

This create a full time series for the period of interest with days and weeks

```r
ts_date <- ts_dwmy_table(InitYear = 2000, LastYear = 2003, WeekDay1 = 'monday')
```

You can then define use the time-series to define your monitoring season with the StartMonth and EndMonth arguments. You can refine this with more arguments (e.g. StartDay, AnchorLength, and so on).

```r
ts_season <- ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)
```

Now that you have defined the monitoring season for the period of interest, you can add your data, starting with the visit date to inform about the sites and the visits

```r
ts_season_visit <- ts_monit_site(m_visit, ts_season)
```

Once the visit data have been integrated, you can add the count for the species of interest, here we use the count data provided with the package for species "2".

```r
ts_season_count <- ts_monit_count_site(ts_season_visit, m_count, sp = 2)
```

##### 3. Compute the yearly flight curve for the data you just created

So far, I have only implemented the regionalGAM method. Here you can filter for the data used in your model by setting the Minimum number of visit, the minimum number of occurrence and number of site.
It is important to remember that the resulting flight curve will depend on the quality of information provided to the model. In other words, "garbage in, garbage out". If you have only occurrence, sites, and few visits, you might have an issue to get a reliable flight curve.

```r
 ts_flight_curve <- flight_curve(ts_season_count, NbrSample = 200, MinVisit = 5, MinOccur = 3, MinNbrSite = 1,
                            MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
```

The output of the flight_curve() function is a list of 3 objects, the phenology, the GAM model, and the data used to fit the GAM. Each of these object can be accessed by specifying the name within the list (e.g. ts_flight_curve$pheno).
So for example, you can produce a simple plot for the flight curves with the following script.

```r

plot(ts_flight_curve$pheno[M_YEAR == 2000, trimDAYNO], ts_flight_curve$pheno[M_YEAR == 2000, NM], type = 'l',
      ylim = c(0, max(ts_flight_curve$pheno[, NM])), xlab = 'Monitoring Year Day', ylab = 'Relative Abundance')
c <- 2
for(y in 2001:2003){
  points(ts_flight_curve$pheno[M_YEAR == y, trimDAYNO], ts_flight_curve$pheno[M_YEAR == y, NM], type = 'l', col = c)
  c <- c + 1
}
legend('topright', legend = c(2000:2003), col = c(seq_along(c(2000:2003))), lty = 1, bty = 'n')
```

##### 4. Impute predicted counts for missing monitoring dates

Using the shape of the flight curve computed in 3, you can now impute expected counts for each site, using a GLM with the site, the observations and the phenology as predictors of daily counts.
This is done with the impute_count function, using the count data (ts_season_count) and the calculated flight curve (ts_flight_curve$pheno).

```r
site_year_sp_count <- impute_count(ts_season_count, ts_flight_curve$pheno, SpeedGlm = FALSE, FamilyGlm = 'quasipoisson')
```

The output is again a list of object, including the imputed counts and the model used. To access and visualize the imputed data, you can
use the object "site_year_sp_count$impute_count".

For example, you can plot the imputed and observed counts for a specific site (e.g. 2) and a year (e.g. 2003).

```r
s <- 2
y <- 2003

plot(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, FITTED],
    ylim=c(0, max(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, COUNT_IMPUTED])),
    col = 'blue', type='l',
    main = paste0('Site ', s, ', Season ', y),
    xlab='Monitoring Month', ylab='Fitted Count')
points(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, COUNT],
       col='red')
```

##### 5. Compute annual site indices

There is no function implemented yet, but you can count the total of weekly butterfly count for each year and site, using one count a week during the defined monitoring period (e.g. Thursday - WEEK_DAY 4).

```r
week_count <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)]
site_total_week_count <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
data.frame(site_total_week_count)

plot(week_count[SITE_ID == 1 & M_YEAR == 2000, .(WEEK, FITTED)], type='l')
points(week_count[SITE_ID == 1 & M_YEAR == 2000, .(WEEK, FITTED)], col = 'red')
```

For reporting errors and issues related to this package and its functions, please open a [issue here](https://github.com/RetoSchmucki/rbms/issues)
