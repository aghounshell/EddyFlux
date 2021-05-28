### Calculating diffusive fluxes from dissolved GHGs using the Cole method
### Following McClure et al. 2018
### 21 May 2021, A Hounshell

# Following code from 2019
# https://github.com/aghounshell/GHG/blob/master/Scripts/GHG_Flux.R
# https://github.com/aghounshell/GHG/blob/master/Scripts/WindComps.R
# https://github.com/aghounshell/GHG/blob/master/Scripts/Atm_CH4.R

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales)

# First load in wind data from Met station at FCR ----
# Download 2020 Met data from EDI
# Downloaded from EDI: 21 May 2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2020.csv")
#download.file(inUrl1,infile1,method="curl")

met_edi <- read.csv("./Data/Met_final_2015_2020.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

# Load met data from 2021 (from GitHub, cleaned w/ script: MET_QAQC_2020.R)
met_21 <- read.csv("./Data/Met_GitHub_2021.csv", header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

# Combine into one data frame
met_all <- rbind(met_edi,met_21)

# Select data every 30 minutes from Jan 2020 to end of met data
met_30 <- as.data.frame(seq(as.POSIXct("2020-01-01", tz = "EST"), as.POSIXct("2021-05-19", tz = "EST"), by = "30 min"))

met_30 <- met_30 %>% 
  rename(DateTime = `seq(as.POSIXct("2020-01-01", tz = "EST"), as.POSIXct("2021-05-19", tz = "EST"), by = "30 min")`)

met_30 <- left_join(met_30,met_all,by="DateTime")

# Gut-check plot: wind
ggplot(met_all,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s))+
  geom_line()+
  ylim(0,11)+
  theme_classic(base_size = 17)

# FOR NOW: Assume atmospheric CH4 = 1893.4 ppb; CO2 = 419.05 ppb
# https://gml.noaa.gov/ccgg/trends/

### Aggregate CO2 and CH4 dissolved data ----
# From EDI (2015-2020) on 28 May 2021
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/5/38d72673295864956cccd6bbba99a1a3" 
infile1 <- paste0(getwd(),"/Data/Dissolved_CO2_CH4_Virginia_Reservoirs.csv")
download.file(inUrl1,infile1,method="curl")

ghg_edi <- read.csv("./Data/Dissolved_CO2_CH4_Virginia_Reservoirs.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(Reservoir == "FCR" & Site == 50 & Depth_m == 0.1) %>% 
  filter(DateTime >= "2020-01-01")

# Then will need to get preliminary GHG data from Reservoir Drive


### OLD CODE ----
### Load in atmospheric CH4 and CO2 data from Beech Island, SC, USA ----
# https://gml.noaa.gov/dv/data/index.php?site=SCT
# Following: https://github.com/aghounshell/GHG/blob/master/Scripts/Atm_CH4.R amd
# https://github.com/aghounshell/GHG/blob/master/Scripts/Atm_CO2.R

# CH4
ch4_atm <- nc_open('./Data/ch4_sct_tower-insitu_1_ccgg_HourlyData.nc')

# Export metadata out as text file
{sink('./Data/ch4_sct.txt')
  print(ch4_atm)
  sink()}

# Extract date/time, elevation, and CH4 (ppb)
datetime <- ncvar_get(ch4_atm,"time")
datetime <- as.POSIXct(datetime, origin = "1970-01-01", tz="UTC")
elevation <- ncvar_get(ch4_atm,"elevation")
ch4 <- ncvar_get(ch4_atm,"value")

# Combine datetime and CH4 into one data frame
methane <- cbind.data.frame(datetime,ch4)
methane <- methane %>% 
  mutate(datetime = format(as.POSIXct(datetime,"%Y-%m-%d %H:%M:%S",tz="UTC"),format="%Y-%m-%d"))

# Average to daily and select for 2020-01-01 to 2021-05-19
ch4_avg <- methane %>% 
  group_by(datetime) %>% 
  summarize_all(funs(mean)) %>% 
  arrange(datetime)


