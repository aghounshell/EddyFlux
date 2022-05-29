### Script to calculate climatolgical averages from met data
### Following JGR-Revisions

### 29 May 2022, A. Hounshell

##############################################################################

# Clear workspace
rm(list = ls())

# Download/load libraries
pacman::p_load(lubridate,readr,ggpubr,ggplot2,dplyr,openair)

# Set working directory - up-date for your specific working directory
wd <- getwd()
setwd(wd)

###############################################################################

# Set new dataframe with list of dates+times:
# every 30 minutes
# Constrain to study time period: 2020-04-05 (time series start date) to
# last 30 minute period - UPDATE WITH EACH NEW SET OF DATA!
ts <- seq.POSIXt(as.POSIXct("2020-04-05 00:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), 
                 as.POSIXct("2022-05-02 13:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), by = "30 min")
ts2 <- data.frame(datetime = ts)

###############################################################################

# Reading in data from the Met Station for gap-filling purposes
# Load data Meteorological data from EDI

# Downloaded from EDI: 06 May 2022
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

met_edi <- read.csv("./Data/Met_final_2015_2021.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

met_2022 <- read.csv("./Data/FCR_Met_final_2022.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

met_all <- rbind(met_edi,met_2022)

# Start timeseries on the 00:15:00 to facilitate 30-min averages
met_all <- met_all %>% 
  filter(DateTime >= as.POSIXct("2019-12-31 00:15:00"))

# Select data every 30 minutes from Jan 2020 to end of met data
met_all$Breaks <- cut(met_all$DateTime,breaks = "30 mins",right=FALSE)
met_all$Breaks <- ymd_hms(as.character(met_all$Breaks))

# Average met data to the 30 min mark (excluding Total Rain and Total PAR)
met_30 <- met_all %>% 
  select(DateTime,BP_Average_kPa,AirTemp_Average_C,RH_percent,ShortwaveRadiationUp_Average_W_m2,ShortwaveRadiationDown_Average_W_m2,
         InfraredRadiationUp_Average_W_m2,InfraredRadiationDown_Average_W_m2,Albedo_Average_W_m2,WindSpeed_Average_m_s,WindDir_degrees,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(mean,na.rm=TRUE)

# Sum met data to the 30 min mark (for Total Rain and Total PAR)
met_30_rain <- met_all %>% 
  select(Rain_Total_mm,PAR_Total_mmol_m2,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(sum,na.rm=TRUE)

# Combine averaged and summed data together
met_30_2 <- cbind.data.frame(met_30,met_30_rain)

# Adjust datetime to 30 minute intervals, select relevant parameters, and rename
# following Brenda's conventions
met_30_2 <- met_30_2 %>% 
  select(-Breaks) %>% 
  mutate(DateTime_Adj = DateTime + 30) %>% 
  select(-DateTime) %>% 
  rename(datetime = DateTime_Adj, AirTC_Avg = AirTemp_Average_C, RH = RH_percent, Pressure = BP_Average_kPa, 
         Rain_sum = Rain_Total_mm, WS_ms_Avg = WindSpeed_Average_m_s, WindDir = WindDir_degrees,SW_in = ShortwaveRadiationUp_Average_W_m2,
         SW_out = ShortwaveRadiationDown_Average_W_m2,LW_in = InfraredRadiationUp_Average_W_m2,LW_out = InfraredRadiationDown_Average_W_m2,
         PAR_Tot_Tot = PAR_Total_mmol_m2,albedo = Albedo_Average_W_m2)

# Join with 30 minute time series
met2 <- left_join(ts2, met_30_2, by = 'datetime')

# Select the study period
met2 <- met2 %>% 
  filter(datetime >= as.POSIXct("2020-05-01 00:00:00") & datetime <= as.POSIXct("2022-04-31 24:00:00"))

###############################################################################

### Calculate climatological information for FCR - during study period
# Include: Avg Temp, Min Temp, Max Temp, Yearly rainfall (both years), AVg wind speed, Max wind speed,
# predominant wind direction
met2 <- met2 %>% 
  mutate(year = ifelse(datetime < as.POSIXct("2021-05-01 00:00:00"), "year1", "year2"),
         WindDir = round(WindDir, digits = 0))

met2 %>% 
  group_by(year) %>% 
  summarise(avg_temp = mean(AirTC_Avg, na.rm = TRUE),
            min_temp = min(AirTC_Avg, na.rm = TRUE),
            max_temp = max(AirTC_Avg, na.rm = TRUE),
            rainfall = sum(Rain_sum, na.rm = TRUE),
            avg_wind = mean(WS_ms_Avg, na.rm = TRUE),
            max_wind = max(WS_ms_Avg, na.rm = TRUE))

met2 %>% 
  summarise(avg_temp = mean(AirTC_Avg, na.rm = TRUE),
            min_temp = min(AirTC_Avg, na.rm = TRUE),
            max_temp = max(AirTC_Avg, na.rm = TRUE),
            rainfall = sum(Rain_sum, na.rm = TRUE),
            avg_wind = mean(WS_ms_Avg, na.rm = TRUE),
            max_wind = max(WS_ms_Avg, na.rm = TRUE))

winddir <- met2 %>% 
  group_by(year) %>% 
  count(WindDir)

winddir_tot <- met2 %>% 
  count(WindDir)

###############################################################################

## Visualize the dominant wind direction
# Save as: 800 x 800
windRose(mydata = met2, ws = "WS_ms_Avg", wd = "WindDir", 
         width = 3, key.position = 'bottom', 
         offset = 3, paddle = FALSE, key.header = 'Wind speed (m/s)', 
         key.footer = ' ', dig.lab = 2, annotate = FALSE,
         angle.scale = 45, ws.int = 1, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
         par.settings = list(fontsize=list(text=25)))

###############################################################################

## Plot daily surface water temp and wind speed throughout the study period

## Load in Catwalk data first
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/271/6/23a191c1870a5b18cbc17f2779f719cf" 
#infile1 <- paste0(getwd(),"/Data/FCR_Catwalk_2018_2021.csv")
#download.file(inUrl1,infile1,method="curl")

catwalk_2021 <- read.csv("./Data/FCR_Catwalk_2018_2021.csv",header = T)%>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= "2020-01-01") %>% 
  select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9, EXOTemp_C_1, EXOSpCond_uScm_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1,
         EXODOsat_percent_1)

catwalk_2022 <- read.csv("./Data/Catwalk_first_QAQC_2018_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= "2022-01-01") %>% 
  select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9, EXOTemp_C_1, EXOSpCond_uScm_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1,
         EXODOsat_percent_1)

catwalk_all <- rbind(catwalk_2021,catwalk_2022) %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00")) %>% 
  mutate(Temp_diff = ThermistorTemp_C_surface - ThermistorTemp_C_9)

## Aggregate to hourly
catwalk_hourly <- catwalk_all %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  ungroup() %>% 
  mutate(Hour = hour(DateTime)) %>% 
  group_by(Hour) %>% 
  dplyr::summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
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

## Aggregate Met data to hourly
met2_hourly <- met2 %>% 
  mutate(datetime = format(as.POSIXct(datetime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(datetime = as.POSIXct(datetime, "%Y-%m-%d %H", tz = "EST")) %>% 
  mutate(Hour = hour(datetime)) %>% 
  ungroup() %>% 
  group_by(Hour) %>% 
  dplyr::summarise(Wind_speed_avg = mean(WS_ms_Avg,na.rm=TRUE),
                   Wind_speed_sd = sd(WS_ms_Avg,na.rm=TRUE))

## Plot
hour_wind <- ggplot(met2_hourly,mapping=aes(x=Hour,y=Wind_speed_avg))+
  geom_ribbon(mapping=aes(x=Hour,y=Wind_speed_avg,ymin=Wind_speed_avg-Wind_speed_sd,ymax=Wind_speed_avg+Wind_speed_sd),fill="grey",alpha=0.5)+
  geom_vline(xintercept = 12, linetype = "dashed")+
  geom_line(size=1)+
  xlab("Hour") + 
  ylab(expression(Wind~Speed~(m~s^-2)))+
  theme_classic(base_size=15)

hour_temp <- ggplot(catwalk_hourly,mapping=aes(x=Hour,y=Temp_C_surface))+
  geom_vline(xintercept = 12, linetype = "dashed")+
  geom_line(size=1)+
  xlab("Hour") + 
  ylab(expression(Temp~(C^o)))+
  theme_classic(base_size=15)

ggarrange(hour_temp,hour_wind,nrow=1,ncol=2,
          labels=c("A.","B."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_DielTempWnd.jpg",width = 8, height=4, units="in",dpi=320)
           