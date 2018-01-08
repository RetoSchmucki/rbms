R --vanilla

install.packages("devtools")
devtools::install_github("RetoSchmucki/rbms", force=TRUE)
library(rbms)

data(m_visit)
data(m_count)

ts_date <- ts_dwmy_table(InitYear = 2000, LastYear = 2003, WeekDay1 = 'monday')

ts_season <- ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)

ts_season_visit <- ts_monit_site(m_visit, ts_season)

ts_season_count <- ts_monit_count_site(ts_season_visit, m_count, sp = 2)

system.time(ts_flight_curve <- flight_curve(ts_season_count, NbrSample = 200, MinVisit = 4, MinOccur = 3, MinNbrSite = 1,
                            MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'quasipoisson', SpeedGam = FALSE, CompltSeason = TRUE))

site_year_sp_count <- impute_count(ts_season_count, ts_flight_curve, SpeedGlm = FALSE, FamilyGlm = 'quasipoisson')

## plot the flight curves
ts_flight_curve

plot(ts_flight_curve[M_YEAR == 2000, trimDAYNO], ts_flight_curve[M_YEAR == 2000, NM], type = 'l', 
      ylim = c(0, max(ts_flight_curve[, NM])), xlab = 'Monitoring Year Day', ylab = 'Relative Abundance')
c <- 2
for(y in 2001:2003){
  points(ts_flight_curve[M_YEAR == y, trimDAYNO], ts_flight_curve[M_YEAR == y, NM], type = 'l', col = c)
  c <- c + 1
}
legend('topright', legend = c(2000:2003), col = c(seq_along(c(2000:2003))), lty = 1, bty = 'n')

unique(site_year_sp_count[, SITE_ID])

s <- 2
y <- 2003

plot(site_year_sp_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count[SITE_ID == s & M_YEAR == y, FITTED],
    ylim=c(0, max(site_year_sp_count[SITE_ID == s & M_YEAR == y, COUNT_IMPUTED])), 
    col = 'blue', type='l', 
    main = paste0('Site ', s, ', Season ', y), 
    xlab='Monitoring Month', ylab='Fitted Count')
points(site_year_sp_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count[SITE_ID == s & M_YEAR == y, COUNT],
       col='red')

site_year_sp_count

## retrieve GLM models 
names(impute_glm_model)
summary(impute_glm_model$glm_mod_2_2003)

## retrieve GAM models
names(flight_curve_model)
summary(flight_curve_model$FlightModel_2_2003)