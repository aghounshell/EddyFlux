### Script to conduct ARIMA models on Eddy flux (CO2 and CH4 data) from 2020-2021

### Corresponds to Table 1

### Following MEL script: https://github.com/melofton/FCR-phytos/blob/master/2_Data_analysis/2B_FP_ARIMA_analyses.R
### A. Hounshell, 03 January 2022

###############################################################################

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

### Load in Eddy Flux data ----
# Data downloaded from EDI: https://doi.org/10.6073/pasta/a1324bcf3e1415268996ba867c636489
# Processed following scripts associated with EDI data publication
eddy_flux <- read_csv("./Data/20211210_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Separate into: Hourly, Daily, and Weekly data
# Hourly
fcr_hourly <- eddy_flux %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  summarise(NEE = mean(NEE_uStar_f, na.rm = TRUE),
            NEE05 = mean(NEE_U05_f, na.rm = TRUE),
            NEE50 = mean(NEE_U50_f, na.rm = TRUE),
            NEE95 = mean(NEE_U95_f, na.rm = TRUE),
            NEE_sd = sd(NEE_uStar_f, na.rm = TRUE),
            CH4 = mean(ch4_flux_uStar_f, na.rm = TRUE),
            CH405 = mean(ch4_flux_U05_f, na.rm = TRUE),
            CH450 = mean(ch4_flux_U50_f, na.rm = TRUE),
            CH495 = mean(ch4_flux_U95_f, na.rm = TRUE),
            CH4_sd = sd(ch4_flux_uStar_f, na.rm = TRUE),
            Tmean = mean(Tair_f, na.rm = TRUE),
            Tmax = max(Tair_f, na.rm = TRUE),
            Tmin = min(Tair_f, na.rm = TRUE),
            H = mean(H_f, na.rm = TRUE),
            LE = mean(LE_f, na.rm = TRUE),
            VPD = mean(VPD, na.rm = TRUE),
            RH = mean(rH, na.rm = TRUE),
            umean = mean(u, na.rm = TRUE),
            umax = max(u),
            umin = min(u),
            pressure = mean(airP, na.rm = TRUE),
            minpress = min(airP, na.rm = TRUE),
            maxpress = max(airP, na.rm = TRUE),
            PAR_tot = mean(PAR_f, na.rm = TRUE),
            precip_sum = sum(precip, na.rm = TRUE),
            Rg = mean(Rg_f, na.rm = TRUE),
            SW_out = mean(SW_out, na.rm = TRUE),
            Rn = mean(Rn_f, na.rm = TRUE),
            LW_in = mean(LW_in, na.rm = TRUE),
            LW_out = mean(LW_out, na.rm = TRUE),
            albedo = mean(albedo, na.rm = TRUE))

# Eddy Flux Daily
fcr_daily <- eddy_flux %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_f, na.rm = TRUE),
                   NEE05 = mean(NEE_U05_f, na.rm = TRUE),
                   NEE50 = mean(NEE_U50_f, na.rm = TRUE),
                   NEE95 = mean(NEE_U95_f, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_f, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_f, na.rm = TRUE),
                   CH405 = mean(ch4_flux_U05_f, na.rm = TRUE),
                   CH450 = mean(ch4_flux_U50_f, na.rm = TRUE),
                   CH495 = mean(ch4_flux_U95_f, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_f, na.rm = TRUE),
                   Tmean = mean(Tair_f, na.rm = TRUE),
                   Tmax = max(Tair_f, na.rm = TRUE),
                   Tmin = min(Tair_f, na.rm = TRUE),
                   H = mean(H_f, na.rm = TRUE),
                   LE = mean(LE_f, na.rm = TRUE),
                   VPD = mean(VPD, na.rm = TRUE),
                   RH = mean(rH, na.rm = TRUE),
                   umean = mean(u, na.rm = TRUE),
                   umax = max(u),
                   umin = min(u),
                   pressure = mean(airP, na.rm = TRUE),
                   minpress = min(airP, na.rm = TRUE),
                   maxpress = max(airP, na.rm = TRUE),
                   PAR_tot = mean(PAR_f, na.rm = TRUE),
                   precip_sum = sum(precip, na.rm = TRUE),
                   Rg = mean(Rg_f, na.rm = TRUE),
                   SW_out = mean(SW_out, na.rm = TRUE),
                   Rn = mean(Rn_f, na.rm = TRUE),
                   LW_in = mean(LW_in, na.rm = TRUE),
                   LW_out = mean(LW_out, na.rm = TRUE),
                   albedo = mean(albedo, na.rm = TRUE))

fcr_daily$DateTime <- as.POSIXct(paste(fcr_daily$Year, fcr_daily$Month, fcr_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# Eddy Flux Weekly (7 days)
fcr_weekly <- eddy_flux %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_f, na.rm = TRUE),
                   NEE05 = mean(NEE_U05_f, na.rm = TRUE),
                   NEE50 = mean(NEE_U50_f, na.rm = TRUE),
                   NEE95 = mean(NEE_U95_f, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_f, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_f, na.rm = TRUE),
                   CH405 = mean(ch4_flux_U05_f, na.rm = TRUE),
                   CH450 = mean(ch4_flux_U50_f, na.rm = TRUE),
                   CH495 = mean(ch4_flux_U95_f, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_f, na.rm = TRUE),
                   Tmean = mean(Tair_f, na.rm = TRUE),
                   Tmax = max(Tair_f, na.rm = TRUE),
                   Tmin = min(Tair_f, na.rm = TRUE),
                   H = mean(H_f, na.rm = TRUE),
                   LE = mean(LE_f, na.rm = TRUE),
                   VPD = mean(VPD, na.rm = TRUE),
                   RH = mean(rH, na.rm = TRUE),
                   umean = mean(u, na.rm = TRUE),
                   umax = max(u),
                   umin = min(u),
                   pressure = mean(airP, na.rm = TRUE),
                   minpress = min(airP, na.rm = TRUE),
                   maxpress = max(airP, na.rm = TRUE),
                   PAR_tot = mean(PAR_f, na.rm = TRUE),
                   precip_sum = sum(precip, na.rm = TRUE),
                   Rg = mean(Rg_f, na.rm = TRUE),
                   SW_out = mean(SW_out, na.rm = TRUE),
                   Rn = mean(Rn_f, na.rm = TRUE),
                   LW_in = mean(LW_in, na.rm = TRUE),
                   LW_out = mean(LW_out, na.rm = TRUE),
                   albedo = mean(albedo, na.rm = TRUE))

# Aggregate to monthly
fcr_monthly <- eddy_flux %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_f, na.rm = TRUE),
                   NEE05 = mean(NEE_U05_f, na.rm = TRUE),
                   NEE50 = mean(NEE_U50_f, na.rm = TRUE),
                   NEE95 = mean(NEE_U95_f, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_f, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_f, na.rm = TRUE),
                   CH405 = mean(ch4_flux_U05_f, na.rm = TRUE),
                   CH450 = mean(ch4_flux_U50_f, na.rm = TRUE),
                   CH495 = mean(ch4_flux_U95_f, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_f, na.rm = TRUE),
                   Tmean = mean(Tair_f, na.rm = TRUE),
                   Tmax = max(Tair_f, na.rm = TRUE),
                   Tmin = min(Tair_f, na.rm = TRUE),
                   H = mean(H_f, na.rm = TRUE),
                   LE = mean(LE_f, na.rm = TRUE),
                   VPD = mean(VPD, na.rm = TRUE),
                   RH = mean(rH, na.rm = TRUE),
                   umean = mean(u, na.rm = TRUE),
                   umax = max(u),
                   umin = min(u),
                   pressure = mean(airP, na.rm = TRUE),
                   minpress = min(airP, na.rm = TRUE),
                   maxpress = max(airP, na.rm = TRUE),
                   PAR_tot = mean(PAR_f, na.rm = TRUE),
                   precip_sum = sum(precip, na.rm = TRUE),
                   Rg = mean(Rg_f, na.rm = TRUE),
                   SW_out = mean(SW_out, na.rm = TRUE),
                   Rn = mean(Rn_f, na.rm = TRUE),
                   LW_in = mean(LW_in, na.rm = TRUE),
                   LW_out = mean(LW_out, na.rm = TRUE),
                   albedo = mean(albedo, na.rm = TRUE))

# Remove 2021 Apr when we don't have much data!!
fcr_monthly <- fcr_monthly[-13,]

### Load in and aggregate environmental data: Exo Sonde (Temp, DO, Chla, FDOM), N2 (from LakeAnalyzer), Temp Diff, Discharge
# Using Catwalk data: updated to EDI on 2 November 2021
# Load in data
catwalk <- read.csv("./Data/Catwalk_EDI_2018-2021.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Calculate difference between surface and bottom temps
catwalk <- catwalk %>% 
  mutate(ThermistorTemp_C_surface = na.fill(na.approx(ThermistorTemp_C_surface,na.rm=FALSE),"extend"),
         EXODO_mgL_1 = na.fill(na.approx(EXODO_mgL_1,na.rm=TRUE),"extend"),
         EXODOsat_percent_1 = na.fill(na.approx(EXODOsat_percent_1,na.rm=TRUE),"extend"),
         EXOChla_ugL_1 = na.fill(na.approx(EXOChla_ugL_1,na.rm=TRUE),"extend"),
         EXOfDOM_RFU_1 = na.fill(na.approx(EXOfDOM_RFU_1,na.rm=TRUE),"extend")) %>% 
  mutate(Temp_diff = ThermistorTemp_C_surface - ThermistorTemp_C_9)

# Load in buoyancy frequency:
la <- read_csv("./Data/FCR_results_LA.csv") %>% 
  mutate(datetime = as.POSIXct(datetime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(datetime >= as.POSIXct("2020-04-05") & datetime < as.POSIXct("2021-04-05")) %>% 
  rename(DateTime = datetime)

# Load in VT discharge:
# From EDI: https://doi.org/10.6073/pasta/c65755d4c0102dde6e3140c1c91b77d6

#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/522/1/537b4cab40700d5fd990a019684212c0" 
#infile1 <- paste0(getwd(),"/Data/inflow_for_EDI_2013_22Oct2021.csv")
#download.file(inUrl1,infile1,method="curl")

q <- read.csv("./Data/inflow_for_EDI_2013_22Oct2021.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-06 01:00:00") & DateTime < as.POSIXct("2021-04-06 01:00:00"))

# Aggregate data by hourly, daily, weekly, and monthly
# Aggregate to hourly: catwalk, discharge
catwalk_hourly <- catwalk %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
            Temp_C_surface_sd = sd(ThermistorTemp_C_surface,na.rm=TRUE),
            DO_mgL = mean(EXODO_mgL_1,na.rm=TRUE),
            DO_mgL_sd = sd(EXODO_mgL_1,na.rm=TRUE),
            DO_sat = mean(EXODOsat_percent_1,na.rm=TRUE),
            DO_sat_sd = sd(EXODOsat_percent_1,na.rm=TRUE),
            Chla_ugL = mean(EXOChla_ugL_1,na.rm=TRUE),
            Chla_ugL_sd = sd(EXOChla_ugL_1,na.rm=TRUE),
            fdom_rfu = mean(EXOfDOM_RFU_1,na.rm=TRUE),
            fdmo_rfu_sd = sd(EXOfDOM_RFU_1,na.rm=TRUE),
            Temp_diff = mean(Temp_diff,na.rm=TRUE),
            Temp_diff_sd = sd(Temp_diff,na.rm=TRUE))
  
q_hourly <- q %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

la_hourly <- la %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

catwalk_hourly_2 <- left_join(catwalk_hourly,q_hourly,by=c("DateTime"))
env_hourly <- left_join(catwalk_hourly_2,la_hourly,by="DateTime")

# Aggregate to daily: catwalk, discharge, N2
catwalk_daily <- catwalk %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
            Temp_C_surface_sd = sd(ThermistorTemp_C_surface,na.rm=TRUE),
            DO_mgL = mean(EXODO_mgL_1,na.rm=TRUE),
            DO_mgL_sd = sd(EXODO_mgL_1,na.rm=TRUE),
            DO_sat = mean(EXODOsat_percent_1,na.rm=TRUE),
            DO_sat_sd = sd(EXODOsat_percent_1,na.rm=TRUE),
            Chla_ugL = mean(EXOChla_ugL_1,na.rm=TRUE),
            Chla_ugL_sd = sd(EXOChla_ugL_1,na.rm=TRUE),
            fdom_rfu = mean(EXOfDOM_RFU_1,na.rm=TRUE),
            fdmo_rfu_sd = sd(EXOfDOM_RFU_1,na.rm=TRUE),
            Temp_diff = mean(Temp_diff,na.rm=TRUE),
            Temp_diff_sd = sd(Temp_diff,na.rm=TRUE))

catwalk_daily$DateTime <- as.POSIXct(paste(catwalk_daily$Year, catwalk_daily$Month, catwalk_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

q_daily <- q %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

q_daily$DateTime <- as.POSIXct(paste(q_daily$Year, q_daily$Month, q_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

la_daily <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

la_daily$DateTime <- as.POSIXct(paste(la_daily$Year, la_daily$Month, la_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

catwalk_daily_2 <- left_join(catwalk_daily,q_daily,by=c("DateTime","Year","Month","Day"))
env_daily <- left_join(catwalk_daily_2,la_daily,by=c("DateTime","Year","Month","Day"))

# Weekly: catwalk, discharge, N2
catwalk_weekly <- catwalk %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
            Temp_C_surface_sd = sd(ThermistorTemp_C_surface,na.rm=TRUE),
            DO_mgL = mean(EXODO_mgL_1,na.rm=TRUE),
            DO_mgL_sd = sd(EXODO_mgL_1,na.rm=TRUE),
            DO_sat = mean(EXODOsat_percent_1,na.rm=TRUE),
            DO_sat_sd = sd(EXODOsat_percent_1,na.rm=TRUE),
            Chla_ugL = mean(EXOChla_ugL_1,na.rm=TRUE),
            Chla_ugL_sd = sd(EXOChla_ugL_1,na.rm=TRUE),
            fdom_rfu = mean(EXOfDOM_RFU_1,na.rm=TRUE),
            fdmo_rfu_sd = sd(EXOfDOM_RFU_1,na.rm=TRUE),
            Temp_diff = mean(Temp_diff,na.rm=TRUE),
            Temp_diff_sd = sd(Temp_diff,na.rm=TRUE))

q_weekly <- q %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>%
  summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

la_weekly <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

catwalk_weekly_2 <- left_join(catwalk_weekly,q_weekly,by=c("Year","Week"))
env_weekly <- left_join(catwalk_weekly_2,la_weekly,by=c("Year","Week"))

# Monthly: catwalk, discharge, N2
catwalk_monthly <- catwalk %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
            Temp_C_surface_sd = sd(ThermistorTemp_C_surface,na.rm=TRUE),
            DO_mgL = mean(EXODO_mgL_1,na.rm=TRUE),
            DO_mgL_sd = sd(EXODO_mgL_1,na.rm=TRUE),
            DO_sat = mean(EXODOsat_percent_1,na.rm=TRUE),
            DO_sat_sd = sd(EXODOsat_percent_1,na.rm=TRUE),
            Chla_ugL = mean(EXOChla_ugL_1,na.rm=TRUE),
            Chla_ugL_sd = sd(EXOChla_ugL_1,na.rm=TRUE),
            fdom_rfu = mean(EXOfDOM_RFU_1,na.rm=TRUE),
            fdmo_rfu_sd = sd(EXOfDOM_RFU_1,na.rm=TRUE),
            Temp_diff = mean(Temp_diff,na.rm=TRUE),
            Temp_diff_sd = sd(Temp_diff,na.rm=TRUE))

q_monthly <- q %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

la_monthly <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

catwalk_monthly_2 <- left_join(catwalk_monthly,q_monthly,by=c("Year","Month"))
env_monthly <- left_join(catwalk_monthly_2,la_monthly,by=c("Year","Month"))

# Remove 2021 Apr when we don't have much data!!
env_monthly <- env_monthly[-13,]

### Plot Environmental Parameters for Manuscript and SI -----
# Figure S2
temp_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Temp_C_surface))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(catwalk,mapping=aes(x=DateTime,y=ThermistorTemp_C_surface),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(Temp~(C^o)))+
  theme_classic(base_size=15)

#temp_time

dosat_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=DO_sat))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(catwalk,mapping=aes(x=DateTime,y=EXODOsat_percent_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab("DO (%)")+
  theme_classic(base_size=15)

#dosat_time

chla_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Chla_ugL))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(catwalk,mapping=aes(x=DateTime,y=EXOChla_ugL_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(Chla~(mu~g~L^-1)))+
  theme_classic(base_size=15)

#chla_time

fdom_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=fdom_rfu))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(catwalk,mapping=aes(x=DateTime,y=EXOfDOM_RFU_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylim(0,11)+
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(FDOM~(R.F.U.)))+
  theme_classic(base_size=15)

#fdom_time

q_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Flow_cms))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(q,mapping=aes(x=DateTime,y=VT_Flow_cms),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(Inflow~(m^3~s^-1)))+
  theme_classic(base_size=15)

#q_time

n2_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=n2))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(N^2))+
  theme_classic(base_size=15)

#n2_time

diff_temp_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Temp_diff))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_point(catwalk,mapping=aes(x=DateTime,y=Temp_diff),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab(expression(Temp~Diff.~(C^o)))+
  theme_classic(base_size=15)

#diff_temp_time

thermo_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=-thermo.depth))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red")+
  geom_line(size=1)+
  xlab("") + 
  scale_x_datetime(labels = date_format("%b"))+
  ylab("Depth (m)")+
  theme_classic(base_size=15)

#thermo_time

ggarrange(temp_time,dosat_time,chla_time,fdom_time,q_time,n2_time,diff_temp_time,thermo_time,ncol=3,nrow=3,
          labels=c("A.","B.","C.","D.","E.","F.","G.","H."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Fig6_Env_data_all.jpg",width = 8, height=8, units="in",dpi=320)

### Organize all data for hourly, daily, and weekly ----
# Used DO_Sat instead of DO_mgL (as they are redundant!)
# Remove Temp_diff and N2 due to correlations
co2_hourly <- left_join(fcr_hourly,env_hourly,by="DateTime") %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(NEE = na.fill(na.approx(NEE,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

hourly_cor = as.data.frame(cor(co2_hourly))

# Remove Temp_diff and N2 due to correlations
ch4_hourly <- left_join(fcr_hourly,env_hourly,by="DateTime") %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(CH4 = na.fill(na.approx(CH4,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

# Remove Temp_diff and N2 due to correlations
co2_daily <- left_join(fcr_daily,env_daily,by=c("DateTime","Year","Month","Day")) %>% 
  ungroup(Year,Month) %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(NEE = na.fill(na.approx(NEE,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

daily_cor = as.data.frame(cor(co2_daily))

# Remove Temp_diff and N2 due to correlations
ch4_daily <- left_join(fcr_daily,env_daily,by=c("DateTime","Year","Month","Day")) %>% 
  ungroup(Year,Month) %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(CH4 = na.fill(na.approx(CH4,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

# Remove Temp_diff and N2 due to correlations
co2_weekly <- left_join(fcr_weekly,env_weekly,by=c("Year","Week")) %>% 
  ungroup(Year,Week) %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(NEE = na.fill(na.approx(NEE,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

weekly_cor = as.data.frame(cor(co2_weekly))

# Remove Temp_diff and N2 due to correlations
ch4_weekly <- left_join(fcr_weekly,env_weekly,by=c("Year","Week")) %>% 
  ungroup(Year,Week) %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(CH4 = na.fill(na.approx(CH4,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

# Remove Temp_diff, N2, fDOM, and Flow_cms due to correlations
co2_monthly <- left_join(fcr_monthly,env_monthly,by=c("Year","Month")) %>% 
  ungroup(Year,Month) %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,
         #fdom_rfu,
         #Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(NEE = na.fill(na.approx(NEE,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         #fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         #Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

monthly_cor = as.data.frame(cor(co2_monthly))

all_cor <- rbind(hourly_cor,daily_cor,weekly_cor,monthly_cor)

# For Table S4
write_csv(all_cor,"./Fig_output/20210822_env_cor.csv")

# Remove Temp_diff, N2, fDOM, and Flow_cms due to correlations
ch4_monthly <- left_join(fcr_monthly,env_monthly,by=c("Year","Month")) %>% 
  ungroup(Year,Month) %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,
         #fdom_rfu,
         #Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth) %>% 
  mutate(CH4 = na.fill(na.approx(CH4,na.rm=FALSE),"extend"),
         Temp_C_surface = na.fill(na.approx(Temp_C_surface,na.rm=FALSE),"extend"),
         DO_sat = na.fill(na.approx(DO_sat,na.rm=FALSE),"extend"),
         Chla_ugL = na.fill(na.approx(Chla_ugL,na.rm = FALSE),"extend"),
         #fdom_rfu = na.fill(na.approx(fdom_rfu,na.rm=FALSE),"extend"),
         #Flow_cms = na.fill(na.approx(Flow_cms,na.rm = FALSE), "extend"),
         #Temp_diff = na.fill(na.approx(Temp_diff,na.rm=FALSE),"extend"),
         #n2 = na.fill(na.approx(n2,na.rm=FALSE),"extend"),
         thermo.depth = na.fill(na.approx(thermo.depth,na.rm=FALSE),"extend"))

### Check for colinearity among variables ----
# CO2 Hourly
chart.Correlation(co2_hourly,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff and n2 - highly correlated with Temp_surface

# Ch4 Hourly
chart.Correlation(ch4_hourly,histogram=TRUE,method=c("pearson"))
# Remove Temp_diff and n2 - highly correlated with Temp_surface

# CO2 daily
chart.Correlation(co2_daily,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_Surface

# CH4 daily
chart.Correlation(ch4_daily,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_Surface

# CO2_weekly
chart.Correlation(co2_weekly,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_Surface

# CH4 weekly
chart.Correlation(ch4_weekly,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_surface

# CO2 Monthly
chart.Correlation(co2_monthly,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_surface
# Remove Flow_cms - correlated with DO_Sat
# Remove FDOM - correlated wtih Chla

# CH4 Monthly
chart.Correlation(ch4_monthly,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_surface
# Remove Flow_cms - correlated with DO_Sat
# Remove FDOM - correlated wtih Chla

### Check for skewness; following MEL script ----
# Function to calculate the cube root w/ complex numbers
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# CO2_hourly
for (i in 1:7){
  print(colnames(co2_hourly)[i])
  var <- co2_hourly[,i]
  hist(as.matrix(var), main = colnames(co2_hourly)[i])
  print(skewness(co2_hourly[,i], na.rm = TRUE))
  print(skewness(log(co2_hourly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(co2_hourly[,i]), na.rm = TRUE))
  print(skewness(co2_hourly[,i]^2),na.rm=TRUE)
  var <- log(co2_hourly[,i])
  hist(as.matrix(var), main = c("Log",colnames(co2_hourly)[i]))
  var <- Math.cbrt(co2_hourly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(co2_hourly)[i]))
  var <- (co2_hourly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(co2_hourly)[i]))
}

# Nothing: NEE, Temp_c_surface, DO_sat, thermo.depth
# Log: Chla
# Cube_rt: fdom_rfu, Flow_cms

# Scale necessary data
co2_hourly_scale <- co2_hourly %>% 
  mutate(Chla_ugL = log(Chla_ugL),
         fdom_rfu = (fdom_rfu)^(1/3),
         Flow_cms = (Flow_cms)^(1/3)) %>% 
  scale()

# ch4_hourly
for (i in 1:7){
  print(colnames(ch4_hourly)[i])
  var <- ch4_hourly[,i]
  hist(as.matrix(var), main = colnames(ch4_hourly)[i])
  print(skewness(ch4_hourly[,i], na.rm = TRUE))
  print(skewness(log(ch4_hourly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(ch4_hourly[,i]), na.rm = TRUE))
  print(skewness(ch4_hourly[,i]^2),na.rm=TRUE)
  var <- log(ch4_hourly[,i])
  hist(as.matrix(var), main = c("Log",colnames(ch4_hourly)[i]))
  var <- Math.cbrt(ch4_hourly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(ch4_hourly)[i]))
  var <- (ch4_hourly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(ch4_hourly)[i]))
}
# Nothing: CH4
# Environmental parameters are the same for CO2
# Scale necessary data
ch4_hourly_scale <- ch4_hourly %>% 
  mutate(Chla_ugL = log(Chla_ugL),
         fdom_rfu = (fdom_rfu)^(1/3),
         Flow_cms = (Flow_cms)^(1/3)) %>% 
  scale()

# CO2 daily
for (i in 1:7){
  print(colnames(co2_daily)[i])
  var <- co2_daily[,i]
  hist(as.matrix(var), main = colnames(co2_daily)[i])
  print(skewness(co2_daily[,i], na.rm = TRUE))
  print(skewness(log(co2_daily[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(co2_daily[,i]), na.rm = TRUE))
  print(skewness(co2_daily[,i]^2),na.rm=TRUE)
  var <- log(co2_daily[,i])
  hist(as.matrix(var), main = c("Log",colnames(co2_daily)[i]))
  var <- Math.cbrt(co2_daily[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(co2_daily)[i]))
  var <- (co2_daily[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(co2_daily)[i]))
}

# Nothing: NEE, Temp,
# Log: DO_sat, Chla,
# Cube Rt: fdom_rfu, Flow_cms, thermo.depth
co2_daily_scale <- co2_daily %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = (fdom_rfu)^(1/3),
         Flow_cms = (Flow_cms)^(1/3),
         thermo.depth = (thermo.depth)^(1/3)) %>% 
  scale()

# ch4 daily
for (i in 1:7){
  print(colnames(ch4_daily)[i])
  var <- ch4_daily[,i]
  hist(as.matrix(var), main = colnames(ch4_daily)[i])
  print(skewness(ch4_daily[,i], na.rm = TRUE))
  print(skewness(log(ch4_daily[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(ch4_daily[,i]), na.rm = TRUE))
  print(skewness(ch4_daily[,i]^2),na.rm=TRUE)
  var <- log(ch4_daily[,i])
  hist(as.matrix(var), main = c("Log",colnames(ch4_daily)[i]))
  var <- Math.cbrt(ch4_daily[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(ch4_daily)[i]))
  var <- (ch4_daily[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(ch4_daily)[i]))
}

# Nothing: CH4
# Everything else the same as CO2
ch4_daily_scale <- ch4_daily %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = (fdom_rfu)^(1/3),
         Flow_cms = (Flow_cms)^(1/3),
         thermo.depth = (thermo.depth)^(1/3)) %>% 
  scale()

# CO2 weekly
for (i in 1:7){
  print(colnames(co2_weekly)[i])
  var <- co2_weekly[,i]
  hist(as.matrix(var), main = colnames(co2_weekly)[i])
  print(skewness(co2_weekly[,i], na.rm = TRUE))
  print(skewness(log(co2_weekly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(co2_weekly[,i]), na.rm = TRUE))
  print(skewness(co2_weekly[,i]^2),na.rm=TRUE)
  var <- log(co2_weekly[,i])
  hist(as.matrix(var), main = c("Log",colnames(co2_weekly)[i]))
  var <- Math.cbrt(co2_weekly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(co2_weekly)[i]))
  var <- (co2_weekly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(co2_weekly)[i]))
}

# Nothing: NEE, Temp, FDOM, thermo.depth
# Log: DO_sat, Chla
# Cube rt: Flow_cms
co2_weekly_scale <- co2_weekly %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         Flow_cms = (Flow_cms)^(1/3)) %>% 
  scale()

# ch4 weekly
for (i in 1:7){
  print(colnames(ch4_weekly)[i])
  var <- ch4_weekly[,i]
  hist(as.matrix(var), main = colnames(ch4_weekly)[i])
  print(skewness(ch4_weekly[,i], na.rm = TRUE))
  print(skewness(log(ch4_weekly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(ch4_weekly[,i]), na.rm = TRUE))
  print(skewness(ch4_weekly[,i]^2),na.rm=TRUE)
  var <- log(ch4_weekly[,i])
  hist(as.matrix(var), main = c("Log",colnames(ch4_weekly)[i]))
  var <- Math.cbrt(ch4_weekly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(ch4_weekly)[i]))
  var <- (ch4_weekly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(ch4_weekly)[i]))
}

# Nothing: CH4
# Environmental parameters the same
ch4_weekly_scale <- ch4_weekly %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         Flow_cms = (Flow_cms)^(1/3)) %>% 
  scale()

# Co2 monthly
for (i in 1:5){
  print(colnames(co2_monthly)[i])
  var <- co2_monthly[,i]
  hist(as.matrix(var), main = colnames(co2_monthly)[i])
  print(skewness(co2_monthly[,i], na.rm = TRUE))
  print(skewness(log(co2_monthly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(co2_monthly[,i]), na.rm = TRUE))
  print(skewness(co2_monthly[,i]^2),na.rm=TRUE)
  var <- log(co2_monthly[,i])
  hist(as.matrix(var), main = c("Log",colnames(co2_monthly)[i]))
  var <- Math.cbrt(co2_monthly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(co2_monthly)[i]))
  var <- (co2_monthly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(co2_monthly)[i]))
}

# Nothing: NEE, Temp, thermo.depth
# Log: DO_Sat, Chla
# Cube rt: 
co2_monthly_scale <- co2_monthly %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL)) %>% 
  scale()

# Ch4 monthly
for (i in 1:5){
  print(colnames(ch4_monthly)[i])
  var <- ch4_monthly[,i]
  hist(as.matrix(var), main = colnames(ch4_monthly)[i])
  print(skewness(ch4_monthly[,i], na.rm = TRUE))
  print(skewness(log(ch4_monthly[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(ch4_monthly[,i]), na.rm = TRUE))
  print(skewness(ch4_monthly[,i]^2),na.rm=TRUE)
  var <- log(ch4_monthly[,i])
  hist(as.matrix(var), main = c("Log",colnames(ch4_monthly)[i]))
  var <- Math.cbrt(ch4_monthly[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(ch4_monthly)[i]))
  var <- (ch4_monthly[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(ch4_monthly)[i]))
}

# Nothing: CH4
# Environmental parameters the same
ch4_monthly_scale <- ch4_monthly %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL)) %>% 
  scale()

### ARIMA models ----
# Following MEL code : )
# THINGS TO CHANGE: 'cols' (change to the environmental variables); best fit!
# CO2 hourly
colnames(co2_hourly_scale)

cols <- c(2:7) # Change this to the environmental variables you want to select!!
sub.final <- NULL
final <- NULL

y <- co2_hourly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(co2_hourly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "co2_hourly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "co2_hourly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for co2_hourly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "co2_hourly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(co2_hourly_scale)[combn(cols,3)[,2]] # Change to reflect the best fit parameters as identified with 'best'
best.vars.cols <- combn(cols,3)[,2] # Change to reflect the best fit parameters as identified with 'best'

best.fit <- auto.arima(y,xreg = as.matrix(co2_hourly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(co2_hourly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(co2_hourly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(co2_hourly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(co2_hourly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(co2_hourly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### Ch4 Hourly
colnames(ch4_hourly_scale)

cols <- c(2:7) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- ch4_hourly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(ch4_hourly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "ch4_hourly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_hourly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for ch4_hourly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "ch4_hourly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(ch4_hourly_scale)[combn(cols,2)[,1]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,2)[,1] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(ch4_hourly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(ch4_hourly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(ch4_hourly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(ch4_hourly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(ch4_hourly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(ch4_hourly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### Co2 Daily
colnames(co2_daily_scale)

cols <- c(2:7) # Update this to the environmental predictor columns!
sub.final <- NULL
final <- NULL

y <- co2_daily_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(co2_daily_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "co2_daily"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "co2_daily"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for co2_daily",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "co2_daily"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(co2_daily_scale)[combn(cols,5)[,1]] # Change to follow 'best'
best.vars.cols <- combn(cols,5)[,1] # Change to follow 'best'

best.fit <- auto.arima(y,xreg = as.matrix(co2_daily_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(co2_daily_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(co2_daily_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(co2_daily_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(co2_daily_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(co2_daily_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### CH4 Daily
colnames(ch4_daily_scale)

cols <- c(2:7) # Change to reflect columns of environmental predictors
sub.final <- NULL
final <- NULL

y <- ch4_daily_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(ch4_daily_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "ch4_daily"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_daily"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for ch4_daily",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "ch4_daily"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(ch4_daily_scale)[combn(cols,3)[,6]] # Change this to reflect the 'best'
best.vars.cols <- combn(cols,3)[,6] # Change this to reflect the 'best'

best.fit <- auto.arima(y,xreg = as.matrix(ch4_daily_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(ch4_daily_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(ch4_daily_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(ch4_daily_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(ch4_daily_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(ch4_daily_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### CO2 Weekly
colnames(co2_weekly_scale)

cols <- c(2:7) # Change to reflect columns of environmental predictors
sub.final <- NULL
final <- NULL

y <- co2_weekly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(co2_weekly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "co2_weekly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "co2_weekly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for co2_weekly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "co2_weekly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(co2_weekly_scale)[combn(cols,3)[,2]] # Change this to reflect 'best'
best.vars.cols <- combn(cols,3)[,2] # Change this to reflect 'best'

best.fit <- auto.arima(y,xreg = as.matrix(co2_weekly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(co2_weekly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(co2_weekly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(co2_weekly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(co2_weekly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(co2_weekly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### CH4 weekly
colnames(ch4_weekly_scale)

cols <- c(2:7) # Change to reflect columns of environmental predictors
sub.final <- NULL
final <- NULL

y <- ch4_weekly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(ch4_weekly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "ch4_weekly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_weekly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for ch4_weekly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "ch4_weekly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(ch4_weekly_scale)[combn(cols,3)[,6]] # Change this to reflect 'best'
best.vars.cols <- combn(cols,3)[,6] # Change this to reflect 'best'

best.fit <- auto.arima(y,xreg = as.matrix(ch4_weekly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(ch4_weekly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(ch4_weekly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(ch4_weekly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  mutate(AICc = as.numeric(AICc)) %>% 
  filter(AICc >= as.numeric(min(AICc)) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(ch4_weekly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(ch4_weekly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### CO2 monthly
colnames(co2_monthly_scale)

cols <- c(2:5) # Change to reflect columns of environmental predictors
sub.final <- NULL
final <- NULL

y <- co2_monthly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(co2_monthly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "co2_monthly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "co2_monthly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for co2_monthly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "co2_monthly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(co2_monthly_scale)[combn(cols,2)[,4]] # Change this to reflect 'best'
best.vars.cols <- combn(cols,2)[,4] # Change this to reflect 'best'

best.fit <- auto.arima(y,xreg = as.matrix(co2_monthly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(co2_monthly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(co2_monthly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(co2_monthly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(co2_monthly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(co2_monthly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

### ch4 monthly
colnames(ch4_monthly_scale)

cols <- c(2:5) # Change to reflect columns of environmental predictors
sub.final <- NULL
final <- NULL

y <- ch4_monthly_scale[,1]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(ch4_monthly_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "ch4_monthly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_monthly"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for ch4_monthly",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "ch4_monthly"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(ch4_monthly_scale)[combn(cols,3)[,3]] # Change this to reflect 'best'
best.vars.cols <- combn(cols,3)[,3] # Change this to reflect 'best'

best.fit <- auto.arima(y,xreg = as.matrix(ch4_monthly_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(ch4_monthly_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(ch4_monthly_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(ch4_monthly_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(ch4_monthly_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(ch4_monthly_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}
