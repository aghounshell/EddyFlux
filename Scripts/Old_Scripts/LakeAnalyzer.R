### Script to load in and format catwalk data for use in Lake Analyzer
### 2 November 2021, A. Hounshell

# Trying to calculate Thermocline depth and N2 for hourly data
# Following script: Cole.R

## Updated 17 May 2022 with additional Catwalk data

###############################################################################

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples,rLakeAnalyzer)

###############################################################################

# Pull catwalk data from EDI - from 2018-2021
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/518/3/295d91f59699995b80bf551d431d75c5" 
#infile1 <- paste0(getwd(),"/Data/Catwalk_EDI_2018-2021.csv")
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

catwalk_all <- rbind(catwalk_2021,catwalk_2022)

catwalk_mean <- catwalk_all %>% 
  group_by(by30 = cut(DateTime, "30 min")) %>% 
  summarise_all(mean,na.rm=TRUE) %>% 
  select(-DateTime) %>% 
  rename(DateTime = by30)

# Format catwalk temp data for use in LakeAnalyzer in Matlab
catwalk_temp <- catwalk_all %>% 
  select(DateTime,ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9) %>% 
  group_by(by30 = cut(DateTime, "30 min")) %>% 
  summarise_all(mean,na.rm=TRUE) %>% 
  select(-DateTime) %>% 
  rename(dateTime = by30, wtr_0.1 = ThermistorTemp_C_surface, wtr_1.0 = ThermistorTemp_C_1, 
         wtr_2.0 = ThermistorTemp_C_2, wtr_3.0 = ThermistorTemp_C_3, wtr_4.0 = ThermistorTemp_C_4,
         wtr_5.0 = ThermistorTemp_C_5, wtr_6.0 = ThermistorTemp_C_6, wtr_7.0 = ThermistorTemp_C_7,
         wtr_8.0 = ThermistorTemp_C_8, wtr_9.0 = ThermistorTemp_C_9)

catwalk_temp$dateTime <- ymd_hms(catwalk_temp$dateTime)

catwalk_temp = data.frame(catwalk_temp)

# Export out for LakeAnalyzer
write.table(catwalk_temp, "./Data/fcr.wtr", sep='\t', row.names=FALSE)

## Load in data for Lake Analyzer - move path to rLakeAnalyzer folder
wtr.path <- system.file('extdata', 'fcr.wtr', package="rLakeAnalyzer")

fcr_temp = load.ts(wtr.path)

thermo_depth = ts.thermo.depth(fcr_temp,seasonal=TRUE)

# Replace NaNs with NA
thermo_depth <- thermo_depth %>% 
  mutate(thermo.depth = ifelse(thermo.depth == "NaN", NA, thermo.depth))

# Plot to look at results
ggplot(thermo_depth,mapping=aes(x=datetime,y=-thermo.depth)) +
  geom_line()+
  xlim(as.POSIXct("2020-04-01"), as.POSIXct("2022-05-05"))

## Calculate N2
n2_freq <- ts.buoyancy.freq(fcr_temp,seasonal=TRUE)

fcr_results_la <- full_join(thermo_depth,n2_freq,by="datetime")

write.csv(fcr_results_la,"./Data/rev_FCR_results_LA.csv")
