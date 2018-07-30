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

# compute some indices

res <- site_year_count[,.(SITE_ID, M_YEAR, COUNT_IMPUTED)][, mean(COUNT_IMPUTED), by = .(SITE_ID, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]
res1 <- site_year_count[,.(SITE_ID, M_YEAR, COUNT)][, mean(COUNT,na.rm=TRUE), by = .(SITE_ID, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]

data.table::setkey(res,"SITE_ID","M_YEAR")
data.table::setkey(res1,"SITE_ID","M_YEAR")
plot(merge(res,res1)[,.(V1.x,V1.y)])

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

## real data
r --vanilla
data(package="RegionalGAM")
data(gatekeeper_CM, package='RegionalGAM')

dataset1 <- data.table::data.table(gatekeeper_CM[,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")])
m_visit_gatekeeper <- dataset1[,DATE:= as.Date(paste0(YEAR,'-',MONTH,'-',DAY),format="%Y-%m-%d")][,.(SITE,DATE)]
m_visit_gatekeeper <- unique(m_visit_gatekeeper)

dataset2 <- data.table::data.table(gatekeeper_CM[gatekeeper_CM$TREND==1,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")])
m_count_gatekeeper <- dataset2[,DATE:= as.Date(paste0(YEAR,'-',MONTH,'-',DAY),format="%Y-%m-%d")]

data.table::setnames(m_visit_gatekeeper,"SITE","SITE_ID")
data.table::setnames(m_count_gatekeeper,"SITE","SITE_ID")


ts_date <- rbms::ts_dwmy_table(InitYear = 2003, LastYear = 2012, WeekDay1 = 'monday')

ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 4, EndMonth = 9, StartDay = 1, EndDay = NULL,
                      CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)


ts_season_visit <- rbms::ts_monit_site(m_visit_gatekeeper, ts_season)

ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, m_count_gatekeeper, sp = 'Pyronia tithonus')

ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 200, MinVisit = 5, MinOccur = 3, MinNbrSite = 1,
                            MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)

f_curve_pheno <- ts_flight_curve$pheno

plot(f_curve_pheno[M_YEAR == 2005, trimDAYNO], f_curve_pheno[M_YEAR == 2005, NM], type = 'l', ylim = c(0, max(f_curve_pheno[, NM])),
    xlab = 'Monitoring Year Day', ylab = 'Relative Abundance')

                            c <- 2
                            for(y in 2004:2012){
                              points(f_curve_pheno[M_YEAR == y, trimDAYNO], f_curve_pheno[M_YEAR == y, NM], type = 'l', col = c)
                              c <- c + 1
                            }
                            legend('topright', legend = c(2004:2012), col = c(seq_along(c(2004:2012))), lty = 1, bty = 'n')


site_year_sp_count <- rbms::impute_count(ts_season_count, f_curve_pheno, SpeedGlm = FALSE, FamilyGlm = 'quasipoisson')

site_year_count <- site_year_sp_count$impute_count
res <- site_year_count[,.(SITE_ID, M_YEAR, WEEK, COUNT_IMPUTED)][, mean(COUNT_IMPUTED), by = .(SITE_ID, WEEK, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]
res1 <- site_year_count[!is.na(COUNT),.(SITE_ID, M_YEAR, WEEK, COUNT)][, mean(COUNT,na.rm = TRUE), by = .(SITE_ID, WEEK, M_YEAR)][, sum(V1), by = .(SITE_ID, M_YEAR)]

hist(res[V1>0,V1])


data.table::setkey(res,"SITE_ID","M_YEAR")
data.table::setkey(res1,"SITE_ID","M_YEAR")
plot(merge(res,res1)[,.(V1.x,V1.y)])


r --vanilla
library(RegionalGAM)
data("gatekeeper_CM")
dataset1 <- gatekeeper_CM[,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]
pheno <- flight_curve(dataset1)
dev.new()
plot(pheno$DAYNO[pheno$year==2005],pheno$nm[pheno$year==2005],pch=19,cex=0.7,type='o',col='red',xlab="day",ylab="relative abundance")

dataset2 <- gatekeeper_CM[gatekeeper_CM$TREND==1,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]

data.index <- abundance_index(dataset2, pheno)
dev.new()
hist(data.index$regional_gam)
