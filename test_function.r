R --vanilla

source('R/rbms_index_modelling.r')
source('R/rbms_organize_data.r')
source('R/rbms_toolbox.r')

load('data/m_count.rda')
load('data/m_visit.rda')


ts_date <- ts_dwmy_table(InitYear = 2000, LastYear = 2003, WeekDay1 = 'monday')

ts_season <- ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)

ts_season_visit <- ts_monit_site(m_visit, ts_season)

ts_season_count <- ts_monit_count_site(ts_season_visit, m_count, sp = 4)


## 3.
## compute the flight curve, using the
## regionalGAM method
##=========================================
library(data.table)
ts_flight_curve <- flight_curve(ts_season_count, NbrSample = 200, MinVisit = 5, MinOccur = 3, MinNbrSite = 1,
                            MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)

f_curve_pheno <- ts_flight_curve$pheno

plot(f_curve_pheno[M_YEAR == 2000, trimDAYNO], f_curve_pheno[M_YEAR == 2000, NM], type = 'l', ylim = c(0, max(f_curve_pheno[, NM])),
    xlab = 'Monitoring Year Day', ylab = 'Relative Abundance')

                            c <- 2
                            for(y in 2001:2003){
                              points(f_curve_pheno[M_YEAR == y, trimDAYNO], f_curve_pheno[M_YEAR == y, NM], type = 'l', col = c)
                              c <- c + 1
                            }
                            legend('topright', legend = c(2000:2003), col = c(seq_along(c(2000:2003))), lty = 1, bty = 'n')

# impute counts
site_year_sp_count <- impute_count(ts_season_count, f_curve_pheno, SpeedGlm = FALSE, FamilyGlm = 'quasipoisson')

site_year_count <- site_year_sp_count$impute_count

res <- site_year_count[,.(SITE_ID, M_YEAR, COUNT_IMPUTED)][, mean(COUNT_IMPUTED), by = .(SITE_ID, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]
res1 <- site_year_count[,.(SITE_ID, M_YEAR, COUNT)][, mean(COUNT,na.rm=TRUE), by = .(SITE_ID, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]

data.table::setkey(res,"SITE_ID","M_YEAR")
data.table::setkey(res1,"SITE_ID","M_YEAR")
plot(merge(res,res1)[,.(V1.x,V1.y)])

plot(

hist(res$V1)


unique(site_year_count[, SITE_ID])
s <- 2
y <- 2003

plot(site_year_count[SITE_ID == s & M_YEAR == y, DATE], site_year_count[SITE_ID == s & M_YEAR == y, FITTED],
    ylim=c(0, max(site_year_count[SITE_ID == s & M_YEAR == y, COUNT_IMPUTED])),
    col = 'blue', type='l',
    main = paste0('Site ', s, ', Season ', y),
    xlab='Monitoring Month', ylab='Fitted Count')


for (i in seq_along(unique(site_year_count[, SITE_ID]))){
    s <- unique(site_year_count[, SITE_ID])[i]
    points(site_year_count[SITE_ID == s & M_YEAR == y, DATE], site_year_count[SITE_ID == s & M_YEAR == y, COUNT], col=i)
    }
