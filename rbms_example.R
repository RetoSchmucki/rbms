R --vanilla

## add a comment that is very important!

if(!requireNamespace("devtools")) install.packages("devtools")
if(!requireNamespace('curl')) install.packages('curl')
devtools::install_github("RetoSchmucki/rbms", force = TRUE)
library(rbms)


## 1.
## load data included in the package
##==================================
data(m_visit)
data(m_count)

# raster::plot(metzger_v3_europe)
# rv <- raster::values(metzger_v3_europe)
# save(rv,file='rv.rds')
# nr <- metzger_v3_europe
# raster::values(nr) <- rv
# raster::plot(nr)

## 2.
## organize the data to cover the time
## period and monitoring season of the BMS
##=========================================
ts_date <- ts_dwmy_table(InitYear = 2000, LastYear = 2003, WeekDay1 = 'monday')

ts_season <- ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)

ts_season_visit <- ts_monit_site(m_visit, ts_season)

ts_season_count <- ts_monit_count_site(ts_season_visit, m_count, sp = 2)


## 3.
## compute the flight curve, using the
## regionalGAM method
##=========================================
system.time(
  ts_flight_curve <- flight_curve(ts_season_count, NbrSample = 200, MinVisit = 5, MinOccur = 3, MinNbrSite = 1,
                            MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
                            )


## 4.
## plot the flight curves
##=========================================
ts_flight_curve

plot(ts_flight_curve$pheno[M_YEAR == 2000, trimDAYNO], ts_flight_curve$pheno[M_YEAR == 2000, NM], type = 'l',
      ylim = c(0, max(ts_flight_curve$pheno[, NM])), xlab = 'Monitoring Year Day', ylab = 'Relative Abundance')
c <- 2
for(y in 2001:2003){
  points(ts_flight_curve$pheno[M_YEAR == y, trimDAYNO], ts_flight_curve$pheno[M_YEAR == y, NM], type = 'l', col = c)
  c <- c + 1
}
legend('topright', legend = c(2000:2003), col = c(seq_along(c(2000:2003))), lty = 1, bty = 'n')


## 5.
## retrieve GAM models
##=========================================
names(ts_flight_curve$model)
summary(ts_flight_curve$model$FlightModel_2_2003)


## 6.
## impute the count for the missing day, using
## the flight curve computed with the
## regionalGAM method
##=========================================
site_year_sp_count <- impute_count(ts_season_count, ts_flight_curve$pheno, SpeedGlm = FALSE, FamilyGlm = 'quasipoisson')


## 7.
## plot the imputed and observed counts for one
## site (e.g. 2) and one year (e.g 2003)
##==============================================


unique(site_year_sp_count[, SITE_ID])
s <- 2
y <- 2003

plot(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, FITTED],
    ylim=c(0, max(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, COUNT_IMPUTED])),
    col = 'blue', type='l',
    main = paste0('Site ', s, ', Season ', y),
    xlab='Monitoring Month', ylab='Fitted Count')
points(site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$impute_count[SITE_ID == s & M_YEAR == y, COUNT],
       col='red')

## 8.
## compute annual site indices
## using the fitted count for the fourth day in the week.
##========================================================

week_count <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)]
site_total_week_count <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
data.frame(site_total_week_count)

plot(week_count[SITE_ID == 1 & M_YEAR == 2000, .(WEEK, FITTED)], type='l')
points(week_count[SITE_ID == 1 & M_YEAR == 2000, .(WEEK, FITTED)], col = 'red')


## 9.
## retrieve GLM models
##==============================================
names(site_year_sp_count$model)
summary(site_year_sp_count$model$glm_mod_2_2003)


# x <- data.frame(longitude = c(4, 4.1, 4.5), latitude = c(50, 50.45, 50.5), id = c('a','b','c'))
# x_value <- get_raster_value(x, OutDf = FALSE, PlotRaster = TRUE)
# get_bioclim(x_value)
# get_bioclim(x_value, OutDf = FALSE)
