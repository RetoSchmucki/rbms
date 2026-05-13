source('R/fit_bam_multi.r')
source('R/rbms_index_modelling_optimised.r')

library(rbms)
data(m_visit)
data(m_count)
library(data.table)

ts_date <- rbms::ts_dwmy_table(InitYear = 2000, LastYear = 2003, WeekDay1 = 'monday')

ts_season <- rbms::ts_monit_season(ts_date,
                        StartMonth = 4,
                        EndMonth = 9, 
                        StartDay = 1,
                        EndDay = NULL,
                        CompltSeason = TRUE,
                        Anchor = TRUE,
                        AnchorLength = 2,
                        AnchorLag = 2,
                        TimeUnit = 'w')

ts_season_visit <- rbms::ts_monit_site(ts_season, m_visit)

ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, m_count, sp = 2)


st.time <- Sys.time()
ts_flight_curve_multi <-flight_curve_multi(ts_season_count, 
                        NbrSample = 300,
                        MinVisit = 3,
                        MinOccur = 2,
                        MinNbrSite = 1,
                        MaxTrial = 4,
                        GamFamily = 'nb',
                        SpeedGam = FALSE,
                        CompltSeason = TRUE,
                        SelectYear = NULL,
                        TimeUnit = 'w')
multi.time <- Sys.time()-st.time

st.time <- Sys.time()
st_flight_curve <- rbms::flight_curve(ts_season_count, 
                        NbrSample = 300,
                        MinVisit = 3,
                        MinOccur = 2,
                        MinNbrSite = 1,
                        MaxTrial = 4,
                        GamFamily = 'nb',
                        SpeedGam = FALSE,
                        CompltSeason = TRUE,
                        SelectYear = NULL,
                        TimeUnit = 'w')
old.time <- Sys.time()-st.time

st.time <- Sys.time()
ts_flight_curve_opti <-flight_curve_optimised(ts_season_count, 
                        NbrSample = 300,
                        MinVisit = 3,
                        MinOccur = 2,
                        MinNbrSite = 1,
                        MaxTrial = 4,
                        GamFamily = 'nb',
                        SpeedGam = FALSE,
                        CompltSeason = TRUE,
                        SelectYear = NULL,
                        TimeUnit = 'w')
opt.time <- Sys.time()-st.time

message(paste("multi:", multi.time, "- old:", old.time, "- optim:", opt.time))

plot_flight_curve_multi(ts_flight_curve_multi, 
                       legend_pos = "topleft",
                       line_width = 3)
dev.new()
plot(st_flight_curve)
dev.new()

plot(ts_flight_curve_opti)
