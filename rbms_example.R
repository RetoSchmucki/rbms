R --vanilla

install.packages("devtools")
library(devtools)
install_github("RetoSchmucki/rbms", force=TRUE)

library(rbms)

?check_names

data(m_visit)
data(m_count)

ts_date <- ts_dwmy_table(InitYear=2000, LastYear=2003, WeekDay1='monday')

ts_season <- ts_monit_season(ts_date, StartMonth=4, EndMonth=9, StartDay=1, EndDay=NULL,
                      CompltSeason=TRUE, Anchor=TRUE, AnchorLength=7, AnchorLag=7)

ts_season_visit <- ts_monit_site(m_visit, ts_season)

ts_season_count <- ts_monit_count_site(ts_season_visit, m_count, sp=2)

system.time(ts_flight_curve <- flight_curve(ts_season_count, NbrSample=100, MinVisit=3, MinOccur=2, MinNbrSite=1,
                            MaxTrial=3, FcMethod='regionalGAM', GamFamily='quasipoisson', SpeedGam=FALSE, CompltSeason=TRUE))

## plot the flight curves
ts_flight_curve$f_pheno

plot(ts_flight_curve$f_pheno[M_YEAR==2000,trimDAYNO],ts_flight_curve$f_pheno[M_YEAR==2000,NM],type='l',ylim=c(0,max(ts_flight_curve$f_pheno[,NM])),xlab='Monitoring Year Day',ylab='Relative Abundance')
c <- 2
for(y in 2001:2003){
  points(ts_flight_curve$f_pheno[M_YEAR==y,trimDAYNO],ts_flight_curve$f_pheno[M_YEAR==y,NM],type='l',col=c)
  c <- c + 1
}
legend('topright',legend=c(2000:2003),col=c(seq_along(c(2000:2003))),lty=1,bty='n')

site_year_sp_count <- impute_count(ts_season_count, ts_flight_curve, SpeedGlm=FALSE, FamilyGlm = 'poisson')

unique(site_year_sp_count$sp_ts_season_count[,SITE_ID])
s <- 2
y <- 2003

plot(site_year_sp_count$sp_ts_season_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$sp_ts_season_count[SITE_ID == s & M_YEAR == y, FITTED],
    ylim=c(0,max(site_year_sp_count$sp_ts_season_count[SITE_ID == s & M_YEAR == y,COUNT_IMPUTED])), col='blue', type='l', main=paste0('Site ', s, ', Season ', y), 
    xlab='Monitoring Month', ylab='Fitted Count')
points(site_year_sp_count$sp_ts_season_count[SITE_ID == s & M_YEAR == y, DATE], site_year_sp_count$sp_ts_season_count[SITE_ID == s & M_YEAR == y, COUNT],
       col='red')

site_year_sp_count$sp_ts_season_count
