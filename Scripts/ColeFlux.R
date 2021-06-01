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
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate)

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
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/5/38d72673295864956cccd6bbba99a1a3" 
#infile1 <- paste0(getwd(),"/Data/Dissolved_CO2_CH4_Virginia_Reservoirs.csv")
#download.file(inUrl1,infile1,method="curl")

ghg_edi <- read.csv("./Data/Dissolved_CO2_CH4_Virginia_Reservoirs.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(Reservoir == "FCR" & Site == 50 & Depth_m == 0.1) %>% 
  filter(DateTime >= "2020-01-01") %>% 
  select(-c("flag_ch4","flag_co2")) %>% 
  mutate(DateTime = DateTime + hours(12) + minutes(00) + seconds(00))

# Then will need to get preliminary GHG data from Reservoir Drive
ghg_git <- read.csv("./Data/20210601_GHGs.csv", header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(Reservoir == "fcr" & Depth_m == 0.1) %>% 
  mutate(ch4_umolL = ifelse(ch4_umolL <= 0.00103, NA, ch4_umolL)) %>% # Remove samples below MDL
  mutate(co2_umolL = ifelse(co2_umolL <= 3.75, NA, co2_umolL)) %>% # Remove samples below MDL
  mutate(Site = 50) %>% 
  relocate(Reservoir, .before = DateTime) %>% 
  relocate(Site, .after = Reservoir) %>% 
  mutate(Reservoir = "FCR")

# Merge edi and git data
ghg <- rbind(ghg_edi,ghg_git)

# Plot to check
ggplot(ghg,mapping=aes(x=DateTime,y=ch4_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

ggplot(ghg,mapping=aes(x=DateTime,y=co2_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

# Separate by rep
ghg_1 <- ghg %>% 
  filter(Rep == 1)

ghg_2 <- ghg %>% 
  filter(Rep == 2)

### Calculate Fluxes ----
# Calculate U10 for wind data
wind <- met_30 %>% 
  select(WindSpeed_Average_m_s) %>% 
  rename(wnd = WindSpeed_Average_m_s)

height <- 2
u10 <- wind.scale.base(wind,height)

# First calculate k600 using LakeMetabolizer
# in m/d
k600_2 <- k.cole.base(u10)

k600 <- cbind(met_30,k600_2) %>% 
  select(c(DateTime:Reservoir,PAR_Average_umol_s_m2:Albedo_Average_W_m2,wnd)) %>% 
  rename(k600_md = wnd)

# Merge and extrapolate GHG data
fluxes <- left_join(k600,ghg_1,by=c("DateTime","Reservoir","Site")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep1 = ch4_umolL, co2_umolL_rep1 = co2_umolL)

fluxes_2 <- left_join(fluxes,ghg_2,by=c("DateTime","Reservoir","Site")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep2 = ch4_umolL, co2_umolL_rep2 = co2_umolL)

# Calculate fluxes: in umol/m2/d
fluxes_2 <- fluxes_2 %>% 
  mutate(nv = (BP_Average_kPa*0.00986923/(0.0820573660809596*(AirTemp_Average_C + 273.15)))) %>% # units = mols/L
  mutate(ch4_umolL_atm = 1893.4/1e9*nv*1e6) %>% # units = umol/L
  mutate(ch4_flux_1 = 1000*k600_md*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_2 = 1000*k600_md*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_umolL_atm = 419.05/1e9*nv*1e6) %>%  # units = umol/L
  mutate(co2_flux_1 = 1000*k600_md*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_2 = 1000*k600_md*(co2_umolL_rep2-co2_umolL_atm)/24/60/60)  # units = umol C/m2/s

fluxes_rep1 <- fluxes_2 %>% 
  select(c(DateTime,ch4_flux_1,co2_flux_1)) %>% 
  mutate(Rep = 1) %>% 
  rename(ch4_flux_umolm2s = ch4_flux_1, co2_flux_umolm2s = co2_flux_1)

fluxes_rep2 <- fluxes_2 %>% 
  select(c(DateTime,ch4_flux_2,co2_flux_2)) %>% 
  mutate(Rep = 2) %>% 
  rename(ch4_flux_umolm2s = ch4_flux_2, co2_flux_umolm2s = co2_flux_2)

fluxes_all <- rbind(fluxes_rep1,fluxes_rep2)

fluxes_all <- fluxes_all %>% 
  arrange(DateTime,Rep) %>% 
  group_by(DateTime) %>% 
  summarise_all(list(mean,sd)) %>% 
  select(-c(Rep_fn1,Rep_fn2)) %>% 
  rename(ch4_umolm2s = ch4_flux_umolm2s_fn1, co2_umolm2s = co2_flux_umolm2s_fn1, ch4_umolm2s_sd = ch4_flux_umolm2s_fn2,
         co2_umolm2s_sd = co2_flux_umolm2s_fn2)

# Plot to check?
ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,ch4_umolm2s))+
  #geom_ribbon(fluxes_all,mapping=aes(DateTime,ymin = ch4_umolm2s-ch4_umolm2s_sd,ymax = ch4_umolm2s+ch4_umolm2s_sd),fill="grey")+
  theme_classic(base_size =15 )

ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,co2_umolm2s))+
  theme_classic(base_size=15)

# Select dates where we actually have GHG concentrations
ghg_fluxes <- left_join(ghg_1,fluxes_all,by="DateTime")

### Load in Eddy Flux data ----
# Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
eddy_flux <- read_csv("./Data/processed_data_upto_2021-05-06_BD.csv")

eddy_flux <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f)

# Plot Eddy flux data and GHG flux data
ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(ghg_fluxes,mapping=aes(x=DateTime,y=co2_umolm2s,color="GHGs"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-29"))+
  theme_classic(base_size = 15)

ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_umolm2s,color="GHGs"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Time")+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-29"))+
  theme_classic(base_size = 15)

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


