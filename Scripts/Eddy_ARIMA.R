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

#install.packages('PerformanceAnalytics')
#library(PerformanceAnalytics)

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

### Load in and aggregate environmental data: Exo Sonde (Temp, DO, Chla, FDOM), N2 (from LakeAnalyzer), Temp Diff, Discharge
# Using Catwalk data: aggregated from EDI and GitHub following 20210715_ColeFlux
# Load in data
catwalk <- read_csv("./Data/Catwalk_all.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00")) %>% 
  select(-X1)

# Calculate difference between surface and bottom temps
catwalk <- catwalk %>% 
  mutate(Temp_diff = ThermistorTemp_C_surface - ThermistorTemp_C_9)

# Load in buoyancy frequency: currently daily; will need to update if we want to include for hourly!
la <- read_csv("./Data/FCR_results_LA.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%m/%d/%Y", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05") & DateTime < as.POSIXct("2021-04-05"))

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

catwalk_daily_2 <- left_join(catwalk_daily,la,by="DateTime")

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

la_weekly <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  summarise(N2 = mean(N2,na.rm=TRUE),
            N2_sd = sd(N2,na.rm=TRUE))

catwalk_weekly_2 <- left_join(catwalk_weekly,la_weekly,by=c("Year","Week"))

### Organize all data for hourly, daily, and weekly ----

######### NOTE: NEED TO FILL IN NA VALUES FOR CATWALK DATA IN THE WINTER!!!
co2_hourly <- left_join(fcr_hourly,catwalk_hourly,by="DateTime") %>% 
  select(NEE,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

ch4_hourly <- left_join(fcr_hourly,catwalk_hourly,by="DateTime") %>% 
  select(CH4,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

co2_daily <- left_join(fcr_daily,catwalk_daily_2,by=c("DateTime","Year","Month","Day")) %>% 
  ungroup(Year,Month) %>% 
  select(NEE,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

ch4_daily <- left_join(fcr_daily,catwalk_daily_2,by=c("DateTime","Year","Month","Day")) %>% 
  ungroup(Year,Month) %>% 
  select(CH4,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

co2_weekly <- left_join(fcr_weekly,catwalk_weekly_2,by=c("Year","Week")) %>% 
  ungroup(Year,Week) %>% 
  select(NEE,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

ch4_weekly <- left_join(fcr_weekly,catwalk_weekly_2,by=c("Year","Week")) %>% 
  ungroup(Year,Week) %>% 
  select(CH4,Temp_C_surface,DO_mgL,DO_sat,Chla_ugL,fdom_rfu)

### Check for colinearity among variables ----
# CO2 Hourly
chart.Correlation(co2_hourly,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff - highly correlated with Temp_surface

# Ch4 Hourly
chart.Correlation(ch4_hourly,histogram=TRUE,method=c("pearson"))
# Remove Temp_diff - highly correlated with Temp_surface

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
chart.CaptureRatios(ch4_weekly,histogram=TRUE,method=c("pearson"))
# Remove N2 and Temp_diff - highly correlated with Temp_surface

### Check for skewness; following MEL script ----
# Function to calculate the cube root w/ complex numbers
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# CO2_hourly
for (i in 1:6){
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

# Nothing: NEE, Temp_c_surface, DO_sat
# Log: Chla
# Cube_rt: fdom_rfu
# Sqrt: DO_mgL

# Scale necessary data
co2_hourly_scale <- co2_hourly %>% 
  mutate(DO_mgL = (DO_mgL^2),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = Math.cbrt(fdom_rfu)) %>% 
  scale()

# ch4_hourly
for (i in 1:6){
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
  mutate(DO_mgL = (DO_mgL^2),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = Math.cbrt(fdom_rfu)) %>% 
  scale()

# CO2 daily
for (i in 1:6){
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
# Cube Rt: fdom_rfu
# Sq: DO
co2_daily_scale <- co2_daily %>% 
  mutate(DO_mgL = (DO_mgL^2),
         DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = Math.cbrt(fdom_rfu)) %>% 
  scale()

# ch4 daily
for (i in 1:6){
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
  mutate(DO_mgL = (DO_mgL^2),
         DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL),
         fdom_rfu = Math.cbrt(fdom_rfu)) %>% 
  scale()

# CO2 weekly
for (i in 1:6){
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

# Nothing: NEE, Temp, FDOM
# Log: DO_sat, Chla
# Cube rt:
# Sq: DO
co2_weekly_scale <- co2_weekly %>% 
  mutate(DO_mgL = (DO_mgL^2),
         DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL)) %>% 
  scale()

# ch4 weekly
for (i in 1:6){
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
  mutate(DO_mgL = (DO_mgL^2),
         DO_sat = log(DO_sat),
         Chla_ugL = log(Chla_ugL)) %>% 
  scale()
