### Script to conduct ARIMA models on Eddy flux (CO2 and CH4 data) from 2020-2021
# Following MEL script: https://github.com/melofton/FCR-phytos/blob/master/2_Data_analysis/2B_FP_ARIMA_analyses.R
# A Hounshell, 19 July 2021

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,lognorm,forecast,utils,igraph,RColorBrewer)

install.packages('PerformanceAnalytics')
library(PerformanceAnalytics)

### Load in Eddy Flux data ----
eddy_flux <- read_csv("./Data/20210615_EC_processed.csv") %>% 
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

fcr_daily$Date <- as.POSIXct(paste(fcr_daily$Year, fcr_daily$Month, fcr_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

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

### Load in and aggregate environmental data: Exo Sonde (Temp, DO, Chla, FDOM), N2 (from LakeAnalyzer), Temp Diff, Discharge

