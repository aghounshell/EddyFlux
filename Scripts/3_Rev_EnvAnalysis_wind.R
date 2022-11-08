### Script to conduct environmental correlation analysis of EC fluxes from FCR
### Apr 2020 - Apr 2022

### ARIMA modeling on measured data ONLY!
### Currently using 'roughly' QA/QC'd environmental data for 2022
### Catwalk, VT discharge, LakeAnalyzer results (via Catwalk data)

### Following revisions from initial submission to JGR-Biogeosciences
### 17 May 2022, A. Hounshell

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,zoo,scales,plyr,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

###############################################################################

## First load in QA/QC'd EC fluxes 

## Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
## Data corrected following FCR_Process_BD
## Data downloaded from: https://pasta-s.lternet.edu/package/data/eml/edi/920/2/9e658ef44de05303dbc496fc25e8c49a

ec <- read.csv("./Data/20220506_EC_processed.csv") 

ec2 <- ec %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00"))

## Visualize measured fluxes at the hourly, daily, weekly, and monthly timescale
fcr_hourly <- ec2 %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
            NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
            CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
            CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

# Aggregate to daily and calculate the variability (SD) - following script for figures_BD
# data in umolm2s
fcr_daily <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_daily$DateTime <- as.POSIXct(paste(fcr_daily$Year, fcr_daily$Month, fcr_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# Aggregate into weekly
fcr_weekly <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_weekly$Date <- make_datetime(year = fcr_weekly$Year) + weeks(fcr_weekly$Week)

# Aggregate to Monthly
fcr_monthly <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_monthly$yearmon <- with(fcr_monthly, sprintf("%d-%02d", Year, Month))

###############################################################################

## Load in environmental variables and aggregate to various timescales

## Start with Catwalk data

## Load in catwalk data - updated on 17 May 2022
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

## Aggregate to hourly
catwalk_hourly <- catwalk_all %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
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

## Catwalk daily
catwalk_daily <- catwalk_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
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

catwalk_daily$DateTime <- as.POSIXct(paste(catwalk_daily$Year, catwalk_daily$Month, catwalk_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Catwalk weekly
catwalk_weekly <- catwalk_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
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

## Catwalk monthly
catwalk_monthly <- catwalk_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
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

###############################################################################

## Load in VT Q data

## Load in data from 2020-2021 via EDI
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/923/1/9e438aa8bcca18bf0ba70c8307aafed9" 
#infile1 <- paste0(getwd(),"/Data/inflow_for_EDI_2013_15May2022.csv")
#download.file(inUrl1,infile1,method="curl")

q_all <- read.csv("./Data/Inflow_2013_May2022.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST"))

q_all <- q_all %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00")) %>% 
  select(DateTime,VT_Flow_cms)
  
## Q Hourly
q_hourly <- q_all %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

## Q daily
q_daily <- q_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

q_daily$DateTime <- as.POSIXct(paste(q_daily$Year, q_daily$Month, q_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Q weekly
q_weekly <- q_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>%
  dplyr::summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

## Q monthly
q_monthly <- q_all %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(Flow_cms = mean(VT_Flow_cms,na.rm=TRUE),
            Flow_cms_sd = sd(VT_Flow_cms,na.rm=TRUE))

###############################################################################

## Load in and format LakeAnalyzer results (thermocline depth and N2)
la <- read_csv("./Data/rev_FCR_results_LA.csv") %>% 
  mutate(datetime = as.POSIXct(datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>% 
  mutate(datetime = with_tz(datetime, tz="EST")) %>% 
  filter(datetime >= as.POSIXct("2020-05-01 01:00:00") & datetime < as.POSIXct("2022-05-01 01:00:00")) %>% 
  select(datetime,thermo.depth,n2)

colnames(la) <- c("DateTime","thermo.depth","n2")

## LA results hourly
la_hourly <- la %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  dplyr::group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

## LA results daily
la_daily <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

la_daily$DateTime <- as.POSIXct(paste(la_daily$Year, la_daily$Month, la_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## LA results weekly
la_weekly <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  dplyr::summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

## LA results monthly
la_monthly <- la %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(thermo.depth = mean(thermo.depth,na.rm=TRUE),
            thermo.depth_sd = sd(thermo.depth,na.rm=TRUE),
            n2 = mean(n2,na.rm=TRUE),
            n2_sd = sd(n2,na.rm=TRUE))

###############################################################################

## Load in wind data - following suggestions from reviewer
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/143/17/02d36541de9088f2dd99d79dc3a7a853" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2022.csv")
#download.file(inUrl1,infile1,method="curl")

met_all <- read.csv("./Data/Met_final_2015_2022.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

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
  filter(DateTime_Adj >= as.POSIXct("2020-05-01 01:00:00") & DateTime_Adj < as.POSIXct("2022-05-01 01:00:00"))

names(met_30_2)[names(met_30_2) == 'DateTime_Adj'] <- 'DateTime'

## Aggregate to hourly and select wind data
met_hourly <- met_30_2 %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  dplyr::group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::summarise(wndspd_m_s_mean = mean(WindSpeed_Average_m_s,na.rm=TRUE),
                   wndspd_m_s_sd = sd(WindSpeed_Average_m_s,na.rm=TRUE))

## Met results daily
met_daily <- met_30_2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(wndspd_m_s_mean = mean(WindSpeed_Average_m_s,na.rm=TRUE),
                   wndspd_m_s_sd = sd(WindSpeed_Average_m_s,na.rm=TRUE))

met_daily$DateTime <- as.POSIXct(paste(met_daily$Year, met_daily$Month, la_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Met results weekly
met_weekly <- met_30_2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime),
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  dplyr::summarise(wndspd_m_s_mean = mean(WindSpeed_Average_m_s,na.rm=TRUE),
                   wndspd_m_s_sd = sd(WindSpeed_Average_m_s,na.rm=TRUE))

## Met results monthly
met_monthly <- met_30_2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(wndspd_m_s_mean = mean(WindSpeed_Average_m_s,na.rm=TRUE),
                   wndspd_m_s_sd = sd(WindSpeed_Average_m_s,na.rm=TRUE))

###############################################################################

## Aggregate all Env data by timescale
env_hourly <- join_all(list(fcr_hourly,catwalk_hourly,q_hourly,la_hourly,met_hourly),by="DateTime",type="left")

env_daily <- join_all(list(fcr_daily,catwalk_daily,q_daily,la_daily,met_daily),by=c("DateTime","Year","Month","Day"),type="left")

env_weekly <- join_all(list(fcr_weekly,catwalk_weekly,q_weekly,la_weekly,met_weekly),by=c("Year","Week"))

env_monthly <- join_all(list(fcr_monthly,catwalk_monthly,q_monthly,la_monthly,met_monthly),by=c("Year","Month"))

###############################################################################

## Plot Environmental Variables for SI
## Fig. S2 and S3
temp_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Temp_C_surface))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(catwalk_all,mapping=aes(x=DateTime,y=ThermistorTemp_C_surface),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab(expression(Temp~(C^o)))+
  theme_classic(base_size=15)

#temp_time

dosat_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=DO_sat))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(catwalk_all,mapping=aes(x=DateTime,y=EXODOsat_percent_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab("DO (%)")+
  theme_classic(base_size=15)

#dosat_time

chla_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Chla_ugL))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(catwalk_all,mapping=aes(x=DateTime,y=EXOChla_ugL_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab(expression(Chla~(mu~g~L^-1)))+
  theme_classic(base_size=15)

#chla_time

fdom_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=fdom_rfu))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(catwalk_all,mapping=aes(x=DateTime,y=EXOfDOM_RFU_1),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylim(0,11)+
  ylab(expression(FDOM~(R.F.U.)))+
  theme_classic(base_size=15)

#fdom_time

q_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Flow_cms))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(q_all,mapping=aes(x=DateTime,y=VT_Flow_cms),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab(expression(Inflow~(m^3~s^-1)))+
  theme_classic(base_size=15)

#q_time

n2_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=n2))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(la,mapping=aes(x=DateTime,y=n2),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab(expression(N^2))+
  theme_classic(base_size=15)

#n2_time

diff_temp_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=Temp_diff))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(catwalk_all,mapping=aes(x=DateTime,y=Temp_diff),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab(expression(Temp~Diff.~(C^o)))+
  theme_classic(base_size=15)

#diff_temp_time

thermo_time <- ggplot(env_daily,mapping=aes(x=DateTime,y=-thermo.depth))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted")+
  geom_point(la,mapping=aes(x=DateTime,y=-thermo.depth),color="lightgrey",alpha=0.1)+
  geom_line(size=1)+
  xlab("") + 
  ylab("Depth (m)")+
  theme_classic(base_size=15)

#thermo_time

ggarrange(temp_time,dosat_time,chla_time,fdom_time,ncol=1,nrow=4,
          labels=c("A.","B.","C.","D."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Env_data_a.jpg",width = 8, height=10, units="in",dpi=320)

ggarrange(q_time,n2_time,diff_temp_time,thermo_time,ncol=1,nrow=4,
          labels=c("A.","B.","C.","D."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Env_data_b.jpg",width = 8, height=10, units="in",dpi=320)

###############################################################################

## Check correlations among environmental variables - what needs to be removed?
## Hourly
co2_hourly <- env_hourly %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

hourly_cor = as.data.frame(cor(co2_hourly,use = "complete.obs"),method=c("pearson"))
chart.Correlation(co2_hourly,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff and N2 (correlated with surface temp)

ch4_hourly <- env_hourly %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

## Daily
co2_daily <- env_daily %>% 
  ungroup() %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

daily_cor = as.data.frame(cor(co2_daily,use = "complete.obs"),method=c("pearson"))
chart.Correlation(co2_daily,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff and N2 (correlated with surface temp)

ch4_daily <- env_daily %>% 
  ungroup() %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

## Weekly
co2_weekly <- env_weekly %>% 
  ungroup() %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

weekly_cor = as.data.frame(cor(co2_weekly,use = "complete.obs"),method=c("pearson"))
chart.Correlation(co2_weekly,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff and N2 (correlated with surface temp)

ch4_weekly <- env_weekly %>% 
  ungroup() %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,wndspd_m_s_mean)

## Monthly
co2_monthly <- env_monthly %>% 
  ungroup() %>% 
  select(NEE,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,
         #wndspd_m_s_mean
         )

monthly_cor = as.data.frame(cor(co2_monthly,use = "complete.obs"),method=c("pearson"))
chart.Correlation(co2_monthly,histogram = TRUE,method=c("pearson"))
# Remove Temp_diff, N2, and wnd_spd (correlated with surface temp)

ch4_monthly <- env_monthly %>% 
  ungroup() %>% 
  select(CH4,Temp_C_surface,DO_sat,Chla_ugL,fdom_rfu,Flow_cms,
         #Temp_diff,
         #n2,
         thermo.depth,
         #wndspd_m_s_mean
         )

all_cor <- rbind(hourly_cor,daily_cor,weekly_cor,monthly_cor)

# For Table S2
write_csv(all_cor,"./Fig_output/20221108_env_cor_wwind.csv")

###############################################################################

## Check for skewness; following MEL script
# Function to calculate the cube root w/ complex numbers
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# CO2_hourly
for (i in 1:8){
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

# Nothing: NEE, Temp, DO_sat, fDOM, thermo.depth
# Log: Flow
# Cube rt: Chla, Wind speed

# Scale necessary data
co2_hourly_scale <- co2_hourly %>% 
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         wndspd_m_s_mean = (wndspd_m_s_mean)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

# ch4_hourly
for (i in 1:8){
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

# Env parameters - same as CO2
# Nothing for CH4

# Scale necessary data
ch4_hourly_scale <- ch4_hourly %>% 
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         wndspd_m_s_mean = (wndspd_m_s_mean)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

## Daily
# CO2 daily
for (i in 1:8){
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

# Nothing: NEE, Temp, DO_sat, fDOM, thermo.depth
# Cube rt: Chla
# Log: Flow, Wnd spd
co2_daily_scale <- co2_daily %>% 
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms),
         wndspd_m_s_mean = log(wndspd_m_s_mean)) %>% 
  scale()

# CH4 daily
for (i in 1:8){
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
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms),
         wndspd_m_s_mean = log(wndspd_m_s_mean)) %>% 
  scale()

## Weekly
for (i in 1:8){
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

# Nothing: NEE, Temp, DO_sat, fDOM, thermo_depth, wnd_spd
# Log: flow
# Cube rt: Chla
co2_weekly_scale <- co2_weekly %>% 
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

# ch4 weekly
for (i in 1:8){
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
  mutate(Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

## Monthly
# Co2 monthly
for (i in 1:7){
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

# Nothing: NEE, Temp, FDOM, thermo_depth
# Log: DO_sat, Flow
# Cube rt: Chla
co2_monthly_scale <- co2_monthly %>% 
  mutate(DO_sat = log(DO_sat),
         Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

# Ch4 monthly
for (i in 1:7){
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
         Chla_ugL = (Chla_ugL)^(1/3),
         Flow_cms = log(Flow_cms)) %>% 
  scale()

###############################################################################

## ARIMA modeling
# Following MEL code : )
# THINGS TO CHANGE: 'cols' (change to the environmental variables); best fit!

###############################################################################

## CO2 Hourly
colnames(co2_hourly_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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

best.vars <- colnames(co2_hourly_scale)[combn(cols,6)[,6]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,6)[,6] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## CH4 Hourly
colnames(ch4_hourly_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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

best.vars <- colnames(ch4_hourly_scale)[combn(cols,4)[,8]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,8] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## CO2 daily
colnames(co2_daily_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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
      sub.sub.final[j,1] <- "ch4_daily"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_daily"
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

best.vars <- colnames(co2_daily_scale)[combn(cols,5)[,11]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,5)[,11] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## ch4 daily
colnames(ch4_daily_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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

best.vars <- colnames(ch4_daily_scale)[combn(cols,5)[,10]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,5)[,10] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## CO2 weekly
colnames(co2_weekly_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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
      sub.sub.final[j,1] <- "ch4_weekly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_weekly"
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

best.vars <- colnames(co2_weekly_scale)[combn(cols,5)[,7]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,5)[,7] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## CH4 weekly
colnames(ch4_weekly_scale)

cols <- c(2:8) # UPDATE THIS TO THE ENV. VARIABLES
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

best.vars <- colnames(ch4_weekly_scale)[combn(cols,4)[,17]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,17] # UPDATE THIS FOLLOWING 'BEST'

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
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
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

###############################################################################

## CO2 monthly
colnames(co2_monthly_scale)

cols <- c(2:7) # UPDATE THIS TO THE ENV. VARIABLES
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
      sub.sub.final[j,1] <- "ch4_monthly"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "ch4_monthly"
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

best.vars <- colnames(co2_monthly_scale)[combn(cols,3)[,9]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,3)[,9] # UPDATE THIS FOLLOWING 'BEST'

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

###############################################################################

## CH4 monthly
colnames(ch4_monthly_scale)

cols <- c(2:7) # UPDATE THIS TO THE ENV. VARIABLES
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

best.vars <- colnames(ch4_monthly_scale)[combn(cols,3)[,10]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,3)[,10] # UPDATE THIS FOLLOWING 'BEST'

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
