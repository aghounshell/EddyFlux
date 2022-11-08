### Script to calculate climatolgical averages from met data
### Following JGR-Revisions

### 29 May 2022, A. Hounshell

##############################################################################

# Clear workspace
rm(list = ls())

# Download/load libraries
pacman::p_load(lubridate,readr,ggpubr,ggplot2,dplyr,openair,clifro,tidyverse,viridis)

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

## Load in preliminary Met 2022 data: following 1a_Rev_2022_MetData.R
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
# Table S6
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

## Look at differences in winter temps
met2 %>% 
  filter(datetime >= "2021-01-01 00:00:00" & datetime <= "2021-02-28 00:00:00") %>% 
  summarise(avg_temp = mean(AirTC_Avg, na.rm = TRUE),
            min_temp = min(AirTC_Avg, na.rm = TRUE),
            max_temp = max(AirTC_Avg, na.rm = TRUE),
            rainfall = sum(Rain_sum, na.rm = TRUE),
            avg_wind = mean(WS_ms_Avg, na.rm = TRUE),
            max_wind = max(WS_ms_Avg, na.rm = TRUE))

met2 %>% 
  filter(datetime >= "2022-01-01 00:00:00" & datetime <= "2022-02-28 00:00:00") %>% 
  summarise(avg_temp = mean(AirTC_Avg, na.rm = TRUE),
            min_temp = min(AirTC_Avg, na.rm = TRUE),
            max_temp = max(AirTC_Avg, na.rm = TRUE),
            rainfall = sum(Rain_sum, na.rm = TRUE),
            avg_wind = mean(WS_ms_Avg, na.rm = TRUE),
            max_wind = max(WS_ms_Avg, na.rm = TRUE))

###############################################################################

## Visualize day vs. night wind
# Fig. S5
met2 <- met2 %>% 
  mutate(daynight = ifelse(SW_in > 0, "Day", "Night"))

with(met2, windrose(WS_ms_Avg, WindDir, daynight, n_col=2,
                    legend_title = 'Wind Speed (m/s)',
                    ggtheme = "bw"))

ggsave("./Fig_Output/SI_DayNight_Wind.jpg",width = 8, height = 4, units="in",dpi=320)

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
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/518/11/5255c0af78097b91d6b32267c32ee7be" 
#infile1 <- paste0(getwd(),"/Data/FCR_Catwalk_2018_2022.csv")
#download.file(inUrl1,infile1,method="curl")

catwalk_all <- read.csv("./Data/FCR_Catwalk_2018_2022.csv",header = T)%>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= "2020-01-01") %>% 
  select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9, EXOTemp_C_1, EXOSpCond_uScm_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1,
         EXODOsat_percent_1)

catwalk_all <- catwalk_all %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00")) %>% 
  mutate(Temp_diff = ThermistorTemp_C_surface - ThermistorTemp_C_9)

catwalk_all <- catwalk_all %>% 
  mutate(year = ifelse(DateTime < as.POSIXct("2021-05-01 00:00:00"), "year1", "year2"))

## Calculate differences between years for key environmental variables
catwalk_all %>% 
  group_by(year) %>% 
  summarise(avg_surf_temp = mean(ThermistorTemp_C_surface, na.rm = TRUE),
            min_surf_temp = min(ThermistorTemp_C_surface, na.rm = TRUE),
            max_surf_temp = max(ThermistorTemp_C_surface, na.rm = TRUE),
            avg_chla = mean(EXOChla_ugL_1, na.rm=TRUE),
            min_chla = min(EXOChla_ugL_1, na.rm = TRUE),
            max_chla = max(EXOChla_ugL_1, na.rm=TRUE),
            avg_fdom = mean(EXOfDOM_RFU_1, na.rm=TRUE),
            min_fdom = min(EXOfDOM_RFU_1, na.rm=TRUE),
            max_fdom = max(EXOfDOM_RFU_1, na.rm=TRUE),
            avg_DO = mean(EXODOsat_percent_1, na.rm=TRUE),
            min_DO = min(EXODOsat_percent_1, na.rm=TRUE),
            max_DO = max(EXODOsat_percent_1, na.rm = TRUE))

catwalk_all %>% 
  summarise(avg_surf_temp = mean(ThermistorTemp_C_surface, na.rm = TRUE),
            avg_chla = mean(EXOChla_ugL_1, na.rm=TRUE),
            avg_fdom = mean(EXOfDOM_RFU_1, na.rm=TRUE),
            avg_DO = mean(EXODOsat_percent_1, na.rm=TRUE))

## Find spring mixing for 2020 and 2021:
## Defined as: first day after ice-off when temp at 1m is similar to 8m
spring_21 <- catwalk_all %>% 
  filter(DateTime >= "2021-02-09 00:00:00" & DateTime <= "2021-03-01 00:00:00") %>% 
  mutate(mixing = ThermistorTemp_C_1-ThermistorTemp_C_8) %>% 
  filter(mixing < abs(1))

spring_22 <- catwalk_all %>% 
  filter(DateTime >= "2022-02-10 00:00:00" & DateTime <= "2022-03-01 00:00:00") %>% 
  mutate(mixing = ThermistorTemp_C_1-ThermistorTemp_C_8) %>% 
  filter(mixing < abs(1))

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

################################################################################
## Load in inflow and calculate mean, min, and max for each year
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/8/cc045f9fe32501138d5f4e1e7f40d492" 
#infile1 <- paste0(getwd(),"/Data/inflow_for_EDI_2013_15May2022.csv")
#download.file(inUrl1,infile1,method="curl")

q_all <- read.csv("./Data/inflow_2013_May2022.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST"))

q_all <- q_all %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00")) %>% 
  select(DateTime,VT_Flow_cms) %>% 
  mutate(year = ifelse(DateTime < as.POSIXct("2021-05-01 00:00:00"), "year1", "year2"))

q_all %>% 
  group_by(year) %>%
  summarise(avg_inflow = mean(VT_Flow_cms, na.rm = TRUE),
            min_inflow = min(VT_Flow_cms, na.rm = TRUE),
            max_inflow = max(VT_Flow_cms, na.rm = TRUE))
 
q_all %>% 
  summarise(avg_inflow = mean(VT_Flow_cms, na.rm = TRUE))

## Calculate daily water residence time
q_daily <- q_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
                   Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

q_daily$DateTime <- as.POSIXct(paste(q_daily$Year, q_daily$Month, q_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Calculate water residence time assuming full pond (3.1 x 10^5 m3 following Howard et al. 2021 and Gerling et al. 2014)
q_daily <- q_daily %>% 
  mutate(wtr_d = (310000)/Flow_cms/60/60/24) %>% 
  mutate(year = ifelse(DateTime < as.POSIXct("2021-05-01 00:00:00"), "year1", "year2"))

## Calculate median and mean water residence time
q_daily %>% 
  group_by(year) %>% 
  summarise(avg_wtr = mean(wtr_d,na.rm=TRUE),
            st_wtr = sd(wtr_d,na.rm=TRUE),
            med_wtr = median(wtr_d,na.rm=TRUE),
            min_wtr = min(wtr_d,na.rm=TRUE),
            max_wtr = max(wtr_d,na.rm=TRUE))

q_daily %>% 
  ungroup() %>% 
  summarise(avg_wtr = mean(wtr_d,na.rm=TRUE),
            st_wtr = sd(wtr_d,na.rm=TRUE),
            med_wtr = median(wtr_d,na.rm=TRUE),
            min_wtr = min(wtr_d,na.rm=TRUE),
            max_wtr = max(wtr_d,na.rm=TRUE))

# Plot daily water residence time
wtr_year1 <- ggplot(q_daily, mapping=aes(x=DateTime,y=wtr_d))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_hline(yintercept = 148,linetype="dashed",color="darkgrey")+ # Mean WTR
  annotate(geom="text", x=as.POSIXct("2021-04-15"), y=170, label="Mean")+
  geom_hline(yintercept = 71.1,linetype="dashed",color="darkgrey")+ # Median WTR
  annotate(geom="text", x=as.POSIXct("2021-04-15"), y=90, label="Median")+
  geom_line(size=1)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-30"))+
  xlab("") + 
  ylab(expression(Water~Residence~Time~(d)))+
  theme_classic(base_size=15)

wtr_year2 <- ggplot(q_daily, mapping=aes(x=DateTime,y=wtr_d))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_hline(yintercept = 347,linetype="dashed",color="darkgrey")+ # Mean WTR
  annotate(geom="text", x=as.POSIXct("2021-05-01"), y=330, label="Mean")+
  geom_hline(yintercept = 383,linetype="dashed",color="darkgrey")+ # Median WTR
  annotate(geom="text", x=as.POSIXct("2021-05-01"), y=400, label="Median")+
  geom_line(size=1)+
  xlim(as.POSIXct("2021-05-01"),as.POSIXct("2022-04-30"))+
  xlab("") + 
  ylab(expression(Water~Residence~Time~(d)))+
  theme_classic(base_size=15)

ggarrange(wtr_year1,wtr_year2,ncol=1,nrow=2,
          labels=c("A.","B."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Daily_wtr.jpg",width = 8, height=7.5, units="in",dpi=320)

################################################################################

## Plot dissolved GHG concentrationst through the water column
## Following comments from reviewers

## Will also want to include GHG data from 2022 - forthcoming!
ghg <- read.csv("./Data/final_GHG_2015-June2022.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
  filter(Reservoir == "FCR" & Site == 50) %>% 
  filter(DateTime >= "2020-01-01") %>% 
  mutate(DateTime = round_date(DateTime, "30 mins")) 

ghg_fluxes <- ghg %>% 
  group_by(Reservoir,Site,DateTime,Depth_m) %>% 
  summarize(ch4_mean = mean(ch4_umolL,na.rm=TRUE),
            ch4_sd = sd(ch4_umolL,na.rm = TRUE),
            co2_mean = mean(co2_umolL,na.rm=TRUE),
            co2_sd = sd(co2_umolL, na.rm=TRUE))

## Plot
ch4_all <- ggplot(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_mean,color=as.factor(Depth_m)))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(mapping=aes(ymin=ch4_mean-ch4_sd,ymax=ch4_mean+ch4_sd),size=1)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  xlab("")+
  ylab(expression(~CH[4]~(mu~mol~L^-1)))+
  scale_color_viridis(discrete=TRUE, name = "Depth (m)")+
  theme_classic(base_size = 15)

ch4_surf <- ghg_fluxes %>% 
  filter(Depth_m == 0.1) %>% 
  ggplot(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_mean,color=as.factor(Depth_m)))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(mapping=aes(ymin=ch4_mean-ch4_sd,ymax=ch4_mean+ch4_sd),size=1)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  xlab("")+
  ylab(expression(~CH[4]~(mu~mol~L^-1)))+
  scale_color_viridis(discrete=TRUE, name = "Depth (m)")+
  theme_classic(base_size = 15)

co2_all <- ggplot(ghg_fluxes,mapping=aes(x=DateTime,y=co2_mean,color=as.factor(Depth_m)))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(mapping=aes(ymin=co2_mean-co2_sd,ymax=co2_mean+co2_sd),size=1)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  xlab("")+
  ylab(expression(~CO[2]~(mu~mol~L^-1)))+
  scale_color_viridis(discrete=TRUE, name = "Depth (m)")+
  theme_classic(base_size = 15)

co2_surf <- ghg_fluxes %>% 
  filter(Depth_m == 0.1) %>% 
  ggplot(ghg_fluxes,mapping=aes(x=DateTime,y=co2_mean,color=as.factor(Depth_m)))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(mapping=aes(ymin=co2_mean-co2_sd,ymax=co2_mean+co2_sd),size=1)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  xlab("")+
  ylab(expression(~CO[2]~(mu~mol~L^-1)))+
  scale_color_viridis(discrete=TRUE, name = "Depth (m)")+
  theme_classic(base_size = 15)

ggarrange(ch4_all,ch4_surf,co2_all,co2_surf,ncol=1,nrow=4,
          labels=c("A.","B.","C.","D."), font.label = list(face="plain",size=15),common.legend=TRUE)

ggsave("./Fig_Output/SI_dGHG.jpg",width = 8, height=12, units="in",dpi=320)
