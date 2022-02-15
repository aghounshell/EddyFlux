### Calculating diffusive fluxes from dissolved GHGs using the Vachon method,
### compare to eddy covariance fluxes, calculate diel differences, calculate 
### differences in intermittent ice-on and full ice-off and calculative
### annual emissions

### Corresponds to Figures: 2, 3, 4, and 6

### Diffusive fluxes calculated following McClure et al. 2018

### Updated on 05 January 2022 to include multiple methods for calculating
### diffusive fluxes following LakeMetabolizer

### A Hounshell, 05 January 2022

##############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples)

### First load in wind data from Met station at FCR ----
## Download 2020 Met data from EDI
## Downloaded from EDI: 21 May 2021
## Updated on 22 December 2021 w/ staged data
## WILL NEED TO UPDATE WITH FINAL PUBLISHED MET DATA!!!
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/143/8/a5524c686e2154ec0fd0459d46a7d1eb" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

met_all <- read.csv("./Data/Met_final_2015_2021.csv", header=T) %>%
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
         InfaredRadiationUp_Average_W_m2,InfaredRadiationDown_Average_W_m2,Albedo_Average_W_m2,WindSpeed_Average_m_s,WindDir_degrees,Breaks) %>% 
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
  rename(DateTime = DateTime_Adj) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Gut-check plot: wind
ggplot(met_all,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s))+
  geom_line()+
  ylim(0,11)+
  theme_classic(base_size = 17)

### Aggregate CO2 and CH4 dissolved data ----
## Downloaded 2022 EDI data
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/6/38d72673295864956cccd6bbba99a1a3" 
#infile1 <- paste0(getwd(),"/Data/final_GHG_2015-2021.csv")
#download.file(inUrl1,infile1,method="curl")

ghg <- read.csv("./Data/final_GHG_2015-2021.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
  mutate(DateTime = round_date(DateTime, "30 mins")) %>% 
  filter(Reservoir == "FCR" & Site == 50 & Depth_m == 0.1) %>% 
  filter(DateTime >= "2020-01-01")

## Remove samples with Flag = 2 (below detection; set to zero)
ghg <- ghg %>% 
  mutate(ch4_umolL = ifelse(flag_ch4 == 2, 0, ch4_umolL),
         co2_umolL = ifelse(flag_co2 == 2, 0, co2_umolL)) %>% 
  select(-flag_ch4,-flag_co2)

## Check for time range when samples were generally collected
ghg_time <- ghg %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00")) %>% 
  mutate(hour = hour(DateTime))

## Plot to check
## Methane
ggplot(ghg,mapping=aes(x=DateTime,y=ch4_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Carbon dioxide
ggplot(ghg,mapping=aes(x=DateTime,y=co2_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Determine time range when samples were generally taken:
ghg_time <- ghg %>%
  filter(DateTime >= as.POSIXct("2020-04-05") & DateTime <= as.POSIXct("2021-04-06")) %>% 
  mutate(hour = hour(DateTime))

## Separate by rep
ghg_1 <- ghg %>% 
  filter(Rep == 1)

ghg_2 <- ghg %>% 
  filter(Rep == 2)

### Calculate Fluxes ----
# Following Lake Metabolizer (Winslow et al. 2016)
# Load in catwalk data for surface water temperature
## Updated with 2022 data
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/271/6/23a191c1870a5b18cbc17f2779f719cf" 
#infile1 <- paste0(getwd(),"/Data/FCR_Catwalk_2018_2021.csv")
#download.file(inUrl1,infile1,method="curl")

catwalk <- read.csv("./Data/Catwalk_EDI_2018-2021.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Calculate average half-hourly values
# Select data every 30 minutes from Jan 2020 to end of met data
catwalk$Breaks <- cut(catwalk$DateTime,breaks = "30 mins",right=FALSE)
catwalk$Breaks <- ymd_hms(as.character(catwalk$Breaks))

# Average to half-hourly catwalk measurements
catwalk_30 <- catwalk %>% 
  select(DateTime,ThermistorTemp_C_surface,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(mean,na.rm=TRUE)

catwalk_30 <- catwalk_30 %>% 
  select(Breaks,ThermistorTemp_C_surface) %>% 
  rename(DateTime = Breaks) %>% 
  mutate(DateTime = force_tz(DateTime,"EST"))

# Combine w/ met data to ensure even timeseries
all_data <- left_join(met_30_2,catwalk_30,by="DateTime")

# Use PAR profiles to estimate Kd for the entire year: from 2013-2020
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/11/d771f5e9956304424c3bc0a39298a5ce" 
#infile1 <- paste0(getwd(),"/Data/CTD_final_2013_2020.csv")
#download.file(inUrl1,infile1,method="curl")

ctd <- read.csv("./Data/CTD_final_2013_2020.csv") %>% 
  mutate(DateTime = as.POSIXct(Date, "%Y-%m-%d %H:%M:%S", tz = "EST"))

ctd <- ctd %>% 
  filter(Reservoir == "FCR", Site == 50) %>% 
  select(DateTime,Depth_m,PAR_umolm2s)

# Find PAR closest to the surface
ctd_surf <- ctd %>% 
  filter(Depth_m >= 0.45 & Depth_m <=0.55) %>% 
  group_by(DateTime) %>% 
  summarise_all(mean,na.rm=TRUE)

# Find PAR closest to 3.5 m
ctd_3 <- ctd %>% 
  filter(Depth_m >= 3.45 & Depth_m <=3.55) %>% 
  group_by(DateTime) %>% 
  summarise_all(mean,na.rm=TRUE)

# PAR total
ctd_par <- full_join(ctd_surf,ctd_3,by="DateTime")

ctd_par <- ctd_par %>% 
  drop_na()

# Calculate kdPAR for each day
ctd_par <- ctd_par %>% 
  mutate(kd = (1/3.5)*log10(PAR_umolm2s.x/PAR_umolm2s.y))

ctd_par_mean <- ctd_par %>% 
  select(kd) %>% 
  summarise_all(mean,na.rm=TRUE)

### Format for analysis in LakeMetabolizer ----
# Extend for continuous timeseries
all_data <- all_data %>% 
  mutate(BP_Average_kPa = na.fill(na.approx(BP_Average_kPa,na.rm=FALSE),"extend")) %>%
  mutate(AirTemp_Average_C = na.fill(na.approx(AirTemp_Average_C,na.rm=FALSE),"extend")) %>%
  mutate(WindSpeed_Average_m_s = na.fill(na.approx(WindSpeed_Average_m_s,na.rm=FALSE),"extend")) %>%
  mutate(RH_percent = na.fill(na.approx(RH_percent,na.rm=FALSE),"extend")) %>%
  mutate(ShortwaveRadiationDown_Average_W_m2 = na.fill(na.approx(ShortwaveRadiationDown_Average_W_m2,na.rm=FALSE),"extend")) %>%
  mutate(ThermistorTemp_C_surface = na.fill(na.approx(ThermistorTemp_C_surface,na.rm=FALSE),"extend")) %>% 
  group_by(DateTime) %>% 
  summarise_all(mean,na.rm=TRUE)
  
# Define variables following convention used in LakeMetabolizer
wnd.z <- rep(3,17520) # in m
Kd <- rep(0.41,17520) # units 1/m
lat <- rep(37.30,17520) # in degrees N
lake.area <- rep(0.119*1000000,17520) # in m2
atm.press <- all_data$BP_Average_kPa*10 # in millibar
dateTime <- all_data$DateTime
Ts <- all_data$ThermistorTemp_C_surface # in deg C
z.aml <- rep(0.15,17520) # set to 0.15 following Heiskanen method
airT <- all_data$AirTemp_Average_C # in deg C
wnd <- all_data$WindSpeed_Average_m_s # in m/s
RH <- all_data$RH_percent # percent
sw <- all_data$ShortwaveRadiationDown_Average_W_m2 # downwelling shortwave radiation in w/m2

# Calculate longwave radiation balance
lwnet <- calc.lw.net.base(dateTime, sw, Ts, lat, atm.press, airT, RH)

# Calculate U10
U10 <- wind.scale.base(wnd,wnd.z)

## First calculate k600 (m/d) using LakeMetabolizer
# Calculate k600 using multiple methods (similar to Erikkilla et al. 2018)
k600_cole <- k.cole.base(U10)

k600_crusius <- k.crusius.base(U10)

k600_vachon <- k.vachon.base(U10, lake.area, params=c(2.51,1.48,0.39))

# Will need to complete these with loops
k600_read <- rep(0,17520)

for (i in 1:length(wnd)){
  k600_read[i] <- k.read.base(wnd.z[i], Kd[i], lat[i], lake.area[i], atm.press[i], dateTime[i], Ts[i], z.aml[i],
                              airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

k600_soloviev <- rep(0,17520)

for (i in 1:length(wnd)){
  k600_soloviev[i] <- k.read.soloviev.base(wnd.z[i], Kd[i], lat[i], lake.area[i], atm.press[i], dateTime[i], Ts[i], z.aml[i],
                                        airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

k600_macIntyre <- rep(0,17520)

for (i in 1:length(wnd)){
  k600_macIntyre[i] <- k.macIntyre.base(wnd.z[i], Kd[i], atm.press[i], dateTime[i], Ts[i], z.aml[i], airT[i], wnd[i], RH[i], sw[i],
                                     lwnet[i], params=c(1.2,0.4872,1.4784))
}

k600_heiskanen <- rep(0,17520)

for (i in 1:length(wnd)){
  k600_heiskanen[i] <- k.heiskanen.base(wnd.z[i], Kd[i], atm.press[i], dateTime[i], Ts[i], z.aml[i], airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

## Then convert all k600 values by temperature and gas following Lake Metabolizer
# Set water temperature <4 deg C to 4 deg C (for K600 conversion)
Ts_qa <- Ts

for (i in 1:length(Ts_qa)){
  if (Ts_qa[i] <= 4){
    Ts_qa[i] = 4
  }
  Ts_qa[i] = Ts_qa[i]
}

# Correct k600 for CO2
k600_cole_co2 <- rep(0,17520)
k600_crusius_co2 <- rep(0,17520)
k600_vachon_co2 <- rep(0,17520)
k600_read_co2 <- rep(0,17520)
k600_soloviev_co2 <- rep(0,17520)
k600_macIntyre_co2 <- rep(0,17520)
k600_heiskanen_co2 <- rep(0,17520)

k600_cole_co2 <- k600.2.kGAS.base(k600_cole,Ts_qa,gas="CO2")
k600_crusius_co2 <- k600.2.kGAS.base(k600_crusius,Ts_qa,gas="CO2")
k600_vachon_co2 <- k600.2.kGAS.base(k600_vachon,Ts_qa,gas="CO2")
k600_read_co2 <- k600.2.kGAS.base(k600_read,Ts_qa,gas="CO2")
k600_soloviev_co2 <- k600.2.kGAS.base(k600_soloviev,Ts_qa,gas="CO2")
k600_macIntyre_co2 <- k600.2.kGAS.base(k600_macIntyre,Ts_qa,gas="CO2")
k600_heiskanen_co2 <- k600.2.kGAS.base(k600_heiskanen,Ts_qa,gas="CO2")

# Correct k600 for CH4
k600_cole_ch4 <- rep(0,17520)
k600_crusius_ch4 <- rep(0,17520)
k600_vachon_ch4 <- rep(0,17520)
k600_read_ch4 <- rep(0,17520)
k600_soloviev_ch4 <- rep(0,17520)
k600_macIntyre_ch4 <- rep(0,17520)
k600_heiskanen_ch4 <- rep(0,17520)

k600_cole_ch4 <- k600.2.kGAS.base(k600_cole,Ts_qa,gas="CH4")
k600_crusius_ch4 <- k600.2.kGAS.base(k600_crusius,Ts_qa,gas="CH4")
k600_vachon_ch4 <- k600.2.kGAS.base(k600_vachon,Ts_qa,gas="CH4")
k600_read_ch4 <- k600.2.kGAS.base(k600_read,Ts_qa,gas="CH4")
k600_soloviev_ch4 <- k600.2.kGAS.base(k600_soloviev,Ts_qa,gas="CH4")
k600_macIntyre_ch4 <- k600.2.kGAS.base(k600_macIntyre,Ts_qa,gas="CH4")
k600_heiskanen_ch4 <- k600.2.kGAS.base(k600_heiskanen,Ts_qa,gas="CH4")

# Aggregate CO2 and CH4 k600 into a single data frame
k600_corr <- data.frame(matrix(ncol = 15, nrow = 17520))

colnames(k600_corr) <- c('DateTime', 'k600_Cole_co2', 'k600_Crusius_co2','k600_Vachon_co2',
                         'k600_Read_co2','k600_Soloviev_co2','k600_MacIntyre_co2','k600_Heiskanen_co2',
                         'k600_Cole_ch4', 'k600_Crusius_ch4','k600_Vachon_ch4',
                         'k600_Read_ch4','k600_Soloviev_ch4','k600_MacIntyre_ch4','k600_Heiskanen_ch4')

k600_corr <- k600_corr %>% 
  mutate(DateTime = dateTime,
         k600_Cole_co2 = k600_cole_co2,
         k600_Crusius_co2 = k600_crusius_co2,
         k600_Vachon_co2 = k600_vachon_co2,
         k600_Read_co2 = k600_read_co2,
         k600_Soloviev_co2 = k600_soloviev_co2,
         k600_MacIntyre_co2 = k600_macIntyre_co2,
         k600_Heiskanen_co2 = k600_heiskanen_co2,
         k600_Cole_ch4 = k600_cole_ch4,
         k600_Crusius_ch4 = k600_crusius_ch4,
         k600_Vachon_ch4 = k600_vachon_ch4,
         k600_Read_ch4 = k600_read_ch4,
         k600_Soloviev_ch4 = k600_soloviev_ch4,
         k600_MacIntyre_ch4 = k600_macIntyre_ch4,
         k600_Heiskanen_ch4 = k600_heiskanen_ch4,
         BP_Average_kPa = atm.press/10,
         AirTemp_Average_C = Ts)

k600_corr <- k600_corr %>% 
  mutate(k600_Read_co2 = na.fill(na.approx(k600_Read_co2,na.rm=FALSE),"extend"),
         k600_Soloviev_co2 = na.fill(na.approx(k600_Soloviev_co2,na.rm=FALSE),"extend"),
         k600_MacIntyre_co2 = na.fill(na.approx(k600_MacIntyre_co2,na.rm=FALSE),"extend"),
         k600_Heiskanen_co2 = na.fill(na.approx(k600_Heiskanen_co2,na.rm=FALSE),"extend"),
         k600_Read_ch4 = na.fill(na.approx(k600_Read_ch4,na.rm=FALSE),"extend"),
         k600_Soloviev_ch4 = na.fill(na.approx(k600_Soloviev_ch4,na.rm=FALSE),"extend"),
         k600_MacIntyre_ch4 = na.fill(na.approx(k600_MacIntyre_ch4,na.rm=FALSE),"extend"),
         k600_Heiskanen_ch4 = na.fill(na.approx(k600_Heiskanen_ch4,na.rm=FALSE),"extend"))

## Save k600 corrected values
write.csv(k600_corr,"./Data/k600_corr.csv", row.names = FALSE)

## Merge and extrapolate GHG data
fluxes <- left_join(k600_corr,ghg_1,by=c("DateTime")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep1 = ch4_umolL, co2_umolL_rep1 = co2_umolL)

fluxes_2 <- left_join(fluxes,ghg_2,by=c("DateTime")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep2 = ch4_umolL, co2_umolL_rep2 = co2_umolL)

## Use EC concentrations (CO2 and CH4) for atmospheric concentrations
## Load in EC data from EDI: https://doi.org/10.6073/pasta/a1324bcf3e1415268996ba867c636489
ec <- read_csv("./Data/FCR_Eddy_up_to_2022-02-11.csv") %>% 
  mutate(DateTime = as.POSIXct(paste(date, time), format="%m/%d/%Y %H:%M:%S", tz="EST"))

## Plot EC concentrations (just to see!)
## Methane
ggplot(ec,mapping=aes(DateTime, co2_mole_fraction_umolmol))+
  geom_line()

## Carbon dioxide
ggplot(ec,mapping=aes(DateTime, ch4_mole_fraction_umolmol))+
  geom_line()

## Do some light QA/QC'ing to remove large spikes in the data
ec_2 <- ec %>% 
  mutate(co2_mole_fraction = ifelse(co2_mole_fraction_umolmol < 430 & co2_mole_fraction_umolmol > 330, co2_mole_fraction_umolmol, NA)) %>% 
  mutate(ch4_mole_fraction = ifelse(ch4_mole_fraction_umolmol < 2.3, ch4_mole_fraction_umolmol, NA)) %>% 
  select(DateTime,co2_mole_fraction,ch4_mole_fraction)

ggplot(ec_2,mapping=aes(DateTime, co2_mole_fraction))+
  geom_point()+
  geom_line()

ggplot(ec_2,mapping=aes(DateTime, ch4_mole_fraction))+
  geom_point()+
  geom_line()

fluxes_2 <- left_join(fluxes_2,ec_2,by="DateTime")

fluxes_2 <- fluxes_2 %>% 
  mutate(co2_mole_fraction = na.fill(na.approx(co2_mole_fraction,na.rm=FALSE),"extend")) %>% 
  mutate(ch4_mole_fraction = na.fill(na.approx(ch4_mole_fraction,na.rm=FALSE),"extend"))

# Calculate fluxes: in umol/m2/s using each k600 method
fluxes_2 <- fluxes_2 %>% 
  mutate(nv = (BP_Average_kPa*0.00986923/(0.0820573660809596*(AirTemp_Average_C + 273.15)))) %>% # units = mols/L
  mutate(ch4_umolL_atm = ch4_mole_fraction/1e6*nv*1e6) %>% # units = umol/L
  mutate(ch4_flux_k600_Cole_1 = 1000*k600_Cole_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Cole_2 = 1000*k600_Cole_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Crusius_1 = 1000*k600_Crusius_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Crusius_2 = 1000*k600_Crusius_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Vachon_1 = 1000*k600_Vachon_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Vachon_2 = 1000*k600_Vachon_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Read_1 = 1000*k600_Read_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Read_2 = 1000*k600_Read_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Soloviev_1 = 1000*k600_Soloviev_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Soloviev_2 = 1000*k600_Soloviev_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_MacIntyre_1 = 1000*k600_MacIntyre_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_MacIntyre_2 = 1000*k600_MacIntyre_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Heiskanen_1 = 1000*k600_Heiskanen_co2*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Heiskanen_2 = 1000*k600_Heiskanen_co2*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_umolL_atm = co2_mole_fraction/1e6*nv*1e6) %>%  # units = umol/L
  mutate(co2_flux_k600_Cole_1 = 1000*k600_Cole_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Cole_2 = 1000*k600_Cole_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Crusius_1 = 1000*k600_Crusius_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Crusius_2 = 1000*k600_Crusius_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Vachon_1 = 1000*k600_Vachon_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Vachon_2 = 1000*k600_Vachon_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Read_1 = 1000*k600_Read_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Read_2 = 1000*k600_Read_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Soloviev_1 = 1000*k600_Soloviev_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Soloviev_2 = 1000*k600_Soloviev_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_MacIntyre_1 = 1000*k600_MacIntyre_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_MacIntyre_2 = 1000*k600_MacIntyre_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Heiskanen_1 = 1000*k600_Heiskanen_co2*(co2_umolL_rep1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Heiskanen_2 = 1000*k600_Heiskanen_co2*(co2_umolL_rep2-co2_umolL_atm)/24/60/60)  # units = umol C/m2/s

fluxes_rep1 <- fluxes_2 %>% 
  select(c(DateTime,ch4_flux_k600_Cole_1,ch4_flux_k600_Crusius_1,ch4_flux_k600_Heiskanen_1,
           ch4_flux_k600_MacIntyre_1,ch4_flux_k600_Read_1,ch4_flux_k600_Soloviev_1,
           ch4_flux_k600_Vachon_1,co2_flux_k600_Cole_1,co2_flux_k600_Crusius_1,co2_flux_k600_Heiskanen_1,
           co2_flux_k600_MacIntyre_1,co2_flux_k600_Read_1,co2_flux_k600_Soloviev_1,
           co2_flux_k600_Vachon_1)) %>% 
  mutate(Rep = 1) %>% 
  rename(ch4_flux_k600_Cole = ch4_flux_k600_Cole_1,
         ch4_flux_k600_Crusius = ch4_flux_k600_Crusius_1,
         ch4_flux_k600_Heiskanen = ch4_flux_k600_Heiskanen_1,
         ch4_flux_k600_MacIntyre = ch4_flux_k600_MacIntyre_1,
         ch4_flux_k600_Read = ch4_flux_k600_Read_1,
         ch4_flux_k600_Soloviev = ch4_flux_k600_Soloviev_1,
         ch4_flux_k600_Vachon = ch4_flux_k600_Vachon_1,
         co2_flux_k600_Cole = co2_flux_k600_Cole_1,
         co2_flux_k600_Crusius = co2_flux_k600_Crusius_1,
         co2_flux_k600_Heiskanen = co2_flux_k600_Heiskanen_1,
         co2_flux_k600_MacIntyre = co2_flux_k600_MacIntyre_1,
         co2_flux_k600_Read = co2_flux_k600_Read_1,
         co2_flux_k600_Soloviev = co2_flux_k600_Soloviev_1,
         co2_flux_k600_Vachon = co2_flux_k600_Vachon_1)

fluxes_rep2 <- fluxes_2 %>% 
  select(c(DateTime,ch4_flux_k600_Cole_2,ch4_flux_k600_Crusius_2,ch4_flux_k600_Heiskanen_2,
           ch4_flux_k600_MacIntyre_2,ch4_flux_k600_Read_2,ch4_flux_k600_Soloviev_2,
           ch4_flux_k600_Vachon_2,co2_flux_k600_Cole_2,co2_flux_k600_Crusius_2,co2_flux_k600_Heiskanen_2,
           co2_flux_k600_MacIntyre_2,co2_flux_k600_Read_2,co2_flux_k600_Soloviev_2,
           co2_flux_k600_Vachon_2)) %>% 
  mutate(Rep = 2)%>% 
  rename(ch4_flux_k600_Cole = ch4_flux_k600_Cole_2,
         ch4_flux_k600_Crusius = ch4_flux_k600_Crusius_2,
         ch4_flux_k600_Heiskanen = ch4_flux_k600_Heiskanen_2,
         ch4_flux_k600_MacIntyre = ch4_flux_k600_MacIntyre_2,
         ch4_flux_k600_Read = ch4_flux_k600_Read_2,
         ch4_flux_k600_Soloviev = ch4_flux_k600_Soloviev_2,
         ch4_flux_k600_Vachon = ch4_flux_k600_Vachon_2,
         co2_flux_k600_Cole = co2_flux_k600_Cole_2,
         co2_flux_k600_Crusius = co2_flux_k600_Crusius_2,
         co2_flux_k600_Heiskanen = co2_flux_k600_Heiskanen_2,
         co2_flux_k600_MacIntyre = co2_flux_k600_MacIntyre_2,
         co2_flux_k600_Read = co2_flux_k600_Read_2,
         co2_flux_k600_Soloviev = co2_flux_k600_Soloviev_2,
         co2_flux_k600_Vachon = co2_flux_k600_Vachon_2)

fluxes_all <- rbind(fluxes_rep1,fluxes_rep2)

fluxes_long_co2 <- fluxes_all %>% 
  select(DateTime,Rep,co2_flux_k600_Cole,co2_flux_k600_Crusius,co2_flux_k600_Heiskanen,co2_flux_k600_MacIntyre,
         co2_flux_k600_Read,co2_flux_k600_Soloviev,co2_flux_k600_Vachon) %>% 
  pivot_longer(!c(DateTime,Rep),names_to = "flux_type",values_to = "co2_flux_umolm2s")

fluxes_long_ch4 <- fluxes_all %>% 
  select(DateTime,Rep,ch4_flux_k600_Cole,ch4_flux_k600_Crusius,ch4_flux_k600_Heiskanen,ch4_flux_k600_MacIntyre,
         ch4_flux_k600_Read,ch4_flux_k600_Soloviev,ch4_flux_k600_Vachon) %>% 
  pivot_longer(!c(DateTime,Rep),names_to = "flux_type",values_to = "ch4_flux_umolm2s")

flux_avg_co2 <- fluxes_long_co2 %>% 
  select(DateTime,co2_flux_umolm2s) %>% 
  group_by(DateTime) %>% 
  summarise_all(funs(mean,sd),na.rm=TRUE)

flux_avg_ch4 <- fluxes_long_ch4 %>% 
  select(DateTime,ch4_flux_umolm2s) %>% 
  group_by(DateTime) %>% 
  summarise_all(funs(mean,sd),na.rm=TRUE)

fluxes_all <- fluxes_all %>% 
  arrange(DateTime,Rep) %>% 
  group_by(DateTime) %>% 
  summarise_all(funs(mean,sd),na.rm=TRUE) %>% 
  select(-c(Rep_mean,Rep_sd))

## Save diffusive GHG fluxes
write.csv(fluxes_all,"diffusive_fluxes_all.csv",row.names = FALSE)

## Select dates where we actually have GHG concentrations
# Make sure dates are on the half hour for all the ghg_fluxes
fluxes_all <- fluxes_all %>% 
  mutate(DateTime = round_date(DateTime, "30 mins"))

ghg_fluxes <- left_join(ghg_1,fluxes_all,by="DateTime")

ghg_fluxes_avg_co2 <- left_join(ghg_1,flux_avg_co2,by="DateTime")

ghg_fluxes_avg_ch4 <- left_join(ghg_1,flux_avg_ch4,by="DateTime")

## Calculate stats for GHG diffusive fluxes (used in Table S5 and S6, aggregated below w/ EC data)
ghg_stats_Cole <- fluxes_all %>% 
  summarise(min_co2_Cole = min(co2_flux_k600_Cole_mean,na.rm = TRUE),
            max_co2_Cole = max(co2_flux_k600_Cole_mean,na.rm=TRUE),
            med_co2_Cole = median(co2_flux_k600_Cole_mean,na.rm=TRUE),
            mean_co2_Cole = mean(co2_flux_k600_Cole_mean,na.rm=TRUE),
            sd_co2_Cole = sd(co2_flux_k600_Cole_mean,na.rm=TRUE),
            cv_co2_Cole = sd(co2_flux_k600_Cole_mean,na.rm=TRUE)/mean(co2_flux_k600_Cole_mean,na.rm=TRUE)*100,
            min_ch4_Cole = min(ch4_flux_k600_Cole_mean,na.rm = TRUE),
            max_ch4_Cole = max(ch4_flux_k600_Cole_mean,na.rm=TRUE),
            med_ch4_Cole = median(ch4_flux_k600_Cole_mean,na.rm=TRUE),
            mean_ch4_Cole = mean(ch4_flux_k600_Cole_mean,na.rm=TRUE),
            sd_ch4_Cole = sd(ch4_flux_k600_Cole_mean,na.rm=TRUE),
            cv_ch4_Cole = sd(ch4_flux_k600_Cole_mean,na.rm=TRUE)/mean(ch4_flux_k600_Cole_mean,na.rm=TRUE)*100)

ghg_stats_Crusius <- fluxes_all %>% 
  summarise(min_co2_Crusius = min(co2_flux_k600_Crusius_mean,na.rm = TRUE),
            max_co2_Crusius = max(co2_flux_k600_Crusius_mean,na.rm=TRUE),
            med_co2_Crusius = median(co2_flux_k600_Crusius_mean,na.rm=TRUE),
            mean_co2_Crusius = mean(co2_flux_k600_Crusius_mean,na.rm=TRUE),
            sd_co2_Crusius = sd(co2_flux_k600_Crusius_mean,na.rm=TRUE),
            cv_co2_Crusius = sd(co2_flux_k600_Crusius_mean,na.rm=TRUE)/mean(co2_flux_k600_Crusius_mean,na.rm=TRUE)*100,
            min_ch4_Crusius = min(ch4_flux_k600_Crusius_mean,na.rm = TRUE),
            max_ch4_Crusius = max(ch4_flux_k600_Crusius_mean,na.rm=TRUE),
            med_ch4_Crusius = median(ch4_flux_k600_Crusius_mean,na.rm=TRUE),
            mean_ch4_Crusius = mean(ch4_flux_k600_Crusius_mean,na.rm=TRUE),
            sd_ch4_Crusius = sd(ch4_flux_k600_Crusius_mean,na.rm=TRUE),
            cv_ch4_Crusius = sd(ch4_flux_k600_Crusius_mean,na.rm=TRUE)/mean(ch4_flux_k600_Crusius_mean,na.rm=TRUE)*100)

ghg_stats_heiskanen <- fluxes_all %>% 
  summarise(min_co2_Heiskanen = min(co2_flux_k600_Heiskanen_mean,na.rm = TRUE),
            max_co2_Heiskanen = max(co2_flux_k600_Heiskanen_mean,na.rm=TRUE),
            med_co2_Heiskanen = median(co2_flux_k600_Heiskanen_mean,na.rm=TRUE),
            mean_co2_Heiskanen = mean(co2_flux_k600_Heiskanen_mean,na.rm=TRUE),
            sd_co2_Heiskanen = sd(co2_flux_k600_Heiskanen_mean,na.rm=TRUE),
            cv_co2_Heiskanen = sd(co2_flux_k600_Heiskanen_mean,na.rm=TRUE)/mean(co2_flux_k600_Heiskanen_mean,na.rm=TRUE)*100,
            min_ch4_Heiskanen = min(ch4_flux_k600_Heiskanen_mean,na.rm = TRUE),
            max_ch4_Heiskanen = max(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE),
            med_ch4_Heiskanen = median(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE),
            mean_ch4_Heiskanen = mean(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE),
            sd_ch4_Heiskanen = sd(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE),
            cv_ch4_Heiskanen = sd(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE)/mean(ch4_flux_k600_Heiskanen_mean,na.rm=TRUE)*100)

ghg_stats_MacIntyre <- fluxes_all %>% 
  summarise(min_co2_MacIntyre = min(co2_flux_k600_MacIntyre_mean,na.rm = TRUE),
            max_co2_MacIntyre = max(co2_flux_k600_MacIntyre_mean,na.rm=TRUE),
            med_co2_MacIntyre = median(co2_flux_k600_MacIntyre_mean,na.rm=TRUE),
            mean_co2_MacIntyre = mean(co2_flux_k600_MacIntyre_mean,na.rm=TRUE),
            sd_co2_MacIntyre = sd(co2_flux_k600_MacIntyre_mean,na.rm=TRUE),
            cv_co2_MacIntyre = sd(co2_flux_k600_MacIntyre_mean,na.rm=TRUE)/mean(co2_flux_k600_MacIntyre_mean,na.rm=TRUE)*100,
            min_ch4_MacIntyre = min(ch4_flux_k600_MacIntyre_mean,na.rm = TRUE),
            max_ch4_MacIntyre = max(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE),
            med_ch4_MacIntyre = median(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE),
            mean_ch4_MacIntyre = mean(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE),
            sd_ch4_MacIntyre = sd(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE),
            cv_ch4_MacIntyre = sd(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE)/mean(ch4_flux_k600_MacIntyre_mean,na.rm=TRUE)*100)

ghg_stats_Read <- fluxes_all %>% 
  summarise(min_co2_Read = min(co2_flux_k600_Read_mean,na.rm = TRUE),
            max_co2_Read = max(co2_flux_k600_Read_mean,na.rm=TRUE),
            med_co2_Read = median(co2_flux_k600_Read_mean,na.rm=TRUE),
            mean_co2_Read = mean(co2_flux_k600_Read_mean,na.rm=TRUE),
            sd_co2_Read = sd(co2_flux_k600_Read_mean,na.rm=TRUE),
            cv_co2_Read = sd(co2_flux_k600_Read_mean,na.rm=TRUE)/mean(co2_flux_k600_Read_mean,na.rm=TRUE)*100,
            min_ch4_Read = min(ch4_flux_k600_Read_mean,na.rm = TRUE),
            max_ch4_Read = max(ch4_flux_k600_Read_mean,na.rm=TRUE),
            med_ch4_Read = median(ch4_flux_k600_Read_mean,na.rm=TRUE),
            mean_ch4_Read = mean(ch4_flux_k600_Read_mean,na.rm=TRUE),
            sd_ch4_Read = sd(ch4_flux_k600_Read_mean,na.rm=TRUE),
            cv_ch4_Read = sd(ch4_flux_k600_Read_mean,na.rm=TRUE)/mean(ch4_flux_k600_Read_mean,na.rm=TRUE)*100)

ghg_stats_Soloviev <- fluxes_all %>% 
  summarise(min_co2_Soloviev = min(co2_flux_k600_Soloviev_mean,na.rm = TRUE),
            max_co2_Soloviev = max(co2_flux_k600_Soloviev_mean,na.rm=TRUE),
            med_co2_Soloviev = median(co2_flux_k600_Soloviev_mean,na.rm=TRUE),
            mean_co2_Soloviev = mean(co2_flux_k600_Soloviev_mean,na.rm=TRUE),
            sd_co2_Soloviev = sd(co2_flux_k600_Soloviev_mean,na.rm=TRUE),
            cv_co2_Soloviev = sd(co2_flux_k600_Soloviev_mean,na.rm=TRUE)/mean(co2_flux_k600_Soloviev_mean,na.rm=TRUE)*100,
            min_ch4_Soloviev = min(ch4_flux_k600_Soloviev_mean,na.rm = TRUE),
            max_ch4_Soloviev = max(ch4_flux_k600_Soloviev_mean,na.rm=TRUE),
            med_ch4_Soloviev = median(ch4_flux_k600_Soloviev_mean,na.rm=TRUE),
            mean_ch4_Soloviev = mean(ch4_flux_k600_Soloviev_mean,na.rm=TRUE),
            sd_ch4_Soloviev = sd(ch4_flux_k600_Soloviev_mean,na.rm=TRUE),
            cv_ch4_Soloviev = sd(ch4_flux_k600_Soloviev_mean,na.rm=TRUE)/mean(ch4_flux_k600_Soloviev_mean,na.rm=TRUE)*100)

ghg_stats_Vachon <- fluxes_all %>% 
  summarise(min_co2_Vachon = min(co2_flux_k600_Vachon_mean,na.rm = TRUE),
            max_co2_Vachon = max(co2_flux_k600_Vachon_mean,na.rm=TRUE),
            med_co2_Vachon = median(co2_flux_k600_Vachon_mean,na.rm=TRUE),
            mean_co2_Vachon = mean(co2_flux_k600_Vachon_mean,na.rm=TRUE),
            sd_co2_Vachon = sd(co2_flux_k600_Vachon_mean,na.rm=TRUE),
            cv_co2_Vachon = sd(co2_flux_k600_Vachon_mean,na.rm=TRUE)/mean(co2_flux_k600_Vachon_mean,na.rm=TRUE)*100,
            min_ch4_Vachon = min(ch4_flux_k600_Vachon_mean,na.rm = TRUE),
            max_ch4_Vachon = max(ch4_flux_k600_Vachon_mean,na.rm=TRUE),
            med_ch4_Vachon = median(ch4_flux_k600_Vachon_mean,na.rm=TRUE),
            mean_ch4_Vachon = mean(ch4_flux_k600_Vachon_mean,na.rm=TRUE),
            sd_ch4_Vachon = sd(ch4_flux_k600_Vachon_mean,na.rm=TRUE),
            cv_ch4_Vachon = sd(ch4_flux_k600_Vachon_mean,na.rm=TRUE)/mean(ch4_flux_k600_Vachon_mean,na.rm=TRUE)*100)

ghg_stats_all_co2 <- fluxes_long_co2 %>% 
  summarise(min_co2 = min(co2_flux_umolm2s,na.rm = TRUE),
            max_co2 = max(co2_flux_umolm2s,na.rm=TRUE),
            med_co2 = median(co2_flux_umolm2s,na.rm=TRUE),
            mean_co2 = mean(co2_flux_umolm2s,na.rm=TRUE),
            sd_co2 = sd(co2_flux_umolm2s,na.rm=TRUE),
            cv_co2 = sd(co2_flux_umolm2s,na.rm=TRUE)/mean(co2_flux_umolm2s,na.rm=TRUE)*100)

ghg_stats_all_ch4 <- fluxes_long_ch4 %>% 
  summarise(min_ch4 = min(ch4_flux_umolm2s,na.rm = TRUE),
            max_ch4 = max(ch4_flux_umolm2s,na.rm=TRUE),
            med_ch4 = median(ch4_flux_umolm2s,na.rm=TRUE),
            mean_ch4 = mean(ch4_flux_umolm2s,na.rm=TRUE),
            sd_ch4 = sd(ch4_flux_umolm2s,na.rm=TRUE),
            cv_ch4 = sd(ch4_flux_umolm2s,na.rm=TRUE)/mean(ch4_flux_umolm2s,na.rm=TRUE)*100)

## Plot Diffusive fluxes for Supplementary (Figure S5 and S6)
ch4_diff_Cole <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Cole",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Cole_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Cole_mean-ch4_flux_k600_Cole_sd,ymax=ch4_flux_k600_Cole_mean+ch4_flux_k600_Cole_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Cole_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Cole_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_crusius <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Crusius",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Crusius_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Crusius_mean-ch4_flux_k600_Crusius_sd,ymax=ch4_flux_k600_Crusius_mean+ch4_flux_k600_Crusius_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Crusius_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Crusius_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_Heiskanen <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Heiskanen",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Heiskanen_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Heiskanen_mean-ch4_flux_k600_Heiskanen_sd,ymax=ch4_flux_k600_Heiskanen_mean+ch4_flux_k600_Heiskanen_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Heiskanen_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Heiskanen_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_MacIntyre <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="MacIntyre",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_MacIntyre_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_MacIntyre_mean-ch4_flux_k600_MacIntyre_sd,ymax=ch4_flux_k600_MacIntyre_mean+ch4_flux_k600_MacIntyre_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_MacIntyre_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_MacIntyre_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_Read <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Read",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Read_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Read_mean-ch4_flux_k600_Read_sd,ymax=ch4_flux_k600_Read_mean+ch4_flux_k600_Read_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Read_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Read_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_Soloviev <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Soloviev",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Soloviev_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Soloviev_mean-ch4_flux_k600_Soloviev_sd,ymax=ch4_flux_k600_Soloviev_mean+ch4_flux_k600_Soloviev_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Soloviev_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Soloviev_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

ch4_diff_Vachon <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=0.095,label="Vachon",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=ch4_flux_k600_Vachon_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Vachon_mean-ch4_flux_k600_Vachon_sd,ymax=ch4_flux_k600_Vachon_mean+ch4_flux_k600_Vachon_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Vachon_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Vachon_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-0.005,0.1)+
  theme_classic(base_size = 15)

co2_diff_Cole <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Cole",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Cole_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Cole_mean-co2_flux_k600_Cole_sd,ymax=co2_flux_k600_Cole_mean+co2_flux_k600_Cole_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Cole_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Cole_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_Crusius <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Crusius",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Crusius_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Crusius_mean-co2_flux_k600_Crusius_sd,ymax=co2_flux_k600_Crusius_mean+co2_flux_k600_Crusius_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Crusius_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Crusius_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_Heiskanen <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Heiskanen",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Heiskanen_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Heiskanen_mean-co2_flux_k600_Heiskanen_sd,ymax=co2_flux_k600_Heiskanen_mean+co2_flux_k600_Heiskanen_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Heiskanen_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Heiskanen_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_MacIntyre <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="MacIntyre",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_MacIntyre_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_MacIntyre_mean-co2_flux_k600_MacIntyre_sd,ymax=co2_flux_k600_MacIntyre_mean+co2_flux_k600_MacIntyre_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_MacIntyre_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_MacIntyre_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_Read <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Read",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Read_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Read_mean-co2_flux_k600_Read_sd,ymax=co2_flux_k600_Read_mean+co2_flux_k600_Read_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Read_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Read_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_Soloviev <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Soloviev",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Soloviev_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Soloviev_mean-co2_flux_k600_Soloviev_sd,ymax=co2_flux_k600_Soloviev_mean+co2_flux_k600_Soloviev_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Soloviev_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Soloviev_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

co2_diff_Vachon <- ggplot(ghg_fluxes)+
  annotate(geom="text",x=as.POSIXct("2020-05-15"),y=35,label="Vachon",size=5)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fluxes_all,mapping=aes(x=DateTime,y=co2_flux_k600_Vachon_mean),color="grey",alpha=0.1)+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Vachon_mean-co2_flux_k600_Vachon_sd,ymax=co2_flux_k600_Vachon_mean+co2_flux_k600_Vachon_sd),size=1,color="#4c8bfe")+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Vachon_mean),size = 1,color="#4c8bfe") +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Vachon_mean),size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  ylim(-2,36)+
  theme_classic(base_size = 15)

ggarrange(co2_diff_Cole,co2_diff_Crusius,co2_diff_Heiskanen,co2_diff_MacIntyre,
          co2_diff_Read, co2_diff_Soloviev,co2_diff_Vachon,
          ncol=2,nrow=4,labels=c("A.","B.","C.","D.","E.","F.","G."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI5_Diff_fluxes_CO2.jpg",width = 9, height=12, units="in",dpi=320)

ggarrange(ch4_diff_Soloviev,ch4_diff_crusius,ch4_diff_Heiskanen,
          ch4_diff_MacIntyre, ch4_diff_Read,ch4_diff_Soloviev,ch4_diff_Vachon,
          ncol=2,nrow=4,labels=c("A.","B.","C.","D.","E.","F.","G."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI5_Diff_fluxes_CH4.jpg",width = 9, height=12, units="in",dpi=320)

### Load in Eddy Flux data ----
## Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
## Data corrected following FCR_Process_BD
## Data downloaded from: https://doi.org/10.6073/pasta/a1324bcf3e1415268996ba867c636489
## And corrected following associated scripts in EDI
eddy_flux <- read_csv("./Data/20211210_EC_processed") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

## Assess differences between measured and gap-filled data
## Figure S1
co2_comp <- ggplot(eddy_flux,mapping=aes(x=NEE_uStar_orig, y=NEE_uStar_fall))+
  geom_abline(intercept = 0, color="darkgrey")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point()+
  stat_smooth(method = "lm")+
  xlab(expression(~Measured~CO[2]~(mu~mol~m^-2~s^-1))) +
  ylab(expression(~Gap-filled~CO[2]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

ch4_comp <- ggplot(eddy_flux,mapping=aes(x=ch4_flux_uStar_orig, y=ch4_flux_uStar_fall))+
  geom_abline(intercept = 0, color="darkgrey")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point()+
  stat_smooth(method = "lm")+
  xlab(expression(~Measured~CH[4]~(mu~mol~m^-2~s^-1))) +
  ylab(expression(~Gap-filled~CH[4]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

co2_time_comp <- ggplot(eddy_flux)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(mapping=aes(x=DateTime, y=NEE_uStar_f, color="Gap-filled"),shape=1)+
  geom_point(mapping=aes(x=DateTime, y=NEE_uStar_orig, color="Measured"))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_time_comp <- ggplot(eddy_flux)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(mapping=aes(x=DateTime, y=ch4_flux_uStar_f, color="Gap-filled"),shape=1)+
  geom_point(mapping=aes(x=DateTime, y=ch4_flux_uStar_orig, color="Measured"))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(co2_time_comp,ch4_time_comp,co2_comp,ch4_comp,ncol=2,nrow=2,
          labels=c("A.","B.","C.","D."),font.label = list(face="plain",size=15),
          common.legend = TRUE)

ggsave("./Fig_Output/SI1_GapFilledComps.jpg",width = 10, height = 9, units="in",dpi=320)

## Aggregate to hourly and calculate the variability (SD) - 
## following script for figures_BD
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

# Calculate min, max, median, mean, standard deviation, and coefficient of variation
fcr_hourly_stats <- fcr_hourly %>% 
  summarise(min_co2 = min(NEE,na.rm = TRUE),
                max_co2 = max(NEE,na.rm=TRUE),
                med_co2 = median(NEE,na.rm=TRUE),
                mean_co2 = mean(NEE,na.rm=TRUE),
                sd_co2 = sd(NEE,na.rm=TRUE),
                cv_co2 = sd(NEE,na.rm=TRUE)/mean(NEE,na.rm=TRUE)*100,
            min_ch4 = min(CH4,na.rm = TRUE),
            max_ch4 = max(CH4,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            mean_ch4 = mean(CH4,na.rm=TRUE),
            sd_ch4 = sd(CH4,na.rm=TRUE),
            cv_ch4 = sd(CH4,na.rm=TRUE)/mean(CH4,na.rm=TRUE)*100)

# Aggregate to daily and calculate the variability (SD) - following script for figures_BD
# data in umolm2s
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

fcr_daily_stats <- fcr_daily %>% 
  ungroup() %>% 
  summarise(min_co2 = min(NEE,na.rm = TRUE),
            max_co2 = max(NEE,na.rm=TRUE),
            med_co2 = median(NEE,na.rm=TRUE),
            mean_co2 = mean(NEE,na.rm=TRUE),
            sd_co2 = sd(NEE,na.rm=TRUE),
            cv_co2 = sd(NEE,na.rm=TRUE)/mean(NEE,na.rm=TRUE)*100,
            min_ch4 = min(CH4,na.rm = TRUE),
            max_ch4 = max(CH4,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            mean_ch4 = mean(CH4,na.rm=TRUE),
            sd_ch4 = sd(CH4,na.rm=TRUE),
            cv_ch4 = sd(CH4,na.rm=TRUE)/mean(CH4,na.rm=TRUE)*100)

# Aggregate into weekly
fcr_weekly <- eddy_flux %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime)) %>% 
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

fcr_weekly$Date <- make_datetime(year = fcr_weekly$Year) + weeks(fcr_weekly$Week)

fcr_weekly_stats <- fcr_weekly %>% 
  ungroup() %>% 
  summarise(min_co2 = min(NEE,na.rm = TRUE),
            max_co2 = max(NEE,na.rm=TRUE),
            med_co2 = median(NEE,na.rm=TRUE),
            mean_co2 = mean(NEE,na.rm=TRUE),
            sd_co2 = sd(NEE,na.rm=TRUE),
            cv_co2 = sd(NEE,na.rm=TRUE)/mean(NEE,na.rm=TRUE)*100,
            min_ch4 = min(CH4,na.rm = TRUE),
            max_ch4 = max(CH4,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            mean_ch4 = mean(CH4,na.rm=TRUE),
            sd_ch4 = sd(CH4,na.rm=TRUE),
            cv_ch4 = sd(CH4,na.rm=TRUE)/mean(CH4,na.rm=TRUE)*100)


# Aggregate to Monthly
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

fcr_monthly$yearmon <- with(fcr_monthly, sprintf("%d-%02d", Year, Month))

fcr_monthly_stats <- fcr_monthly %>% 
  ungroup() %>% 
  summarise(min_co2 = min(NEE,na.rm = TRUE),
            max_co2 = max(NEE,na.rm=TRUE),
            med_co2 = median(NEE,na.rm=TRUE),
            mean_co2 = mean(NEE,na.rm=TRUE),
            sd_co2 = sd(NEE,na.rm=TRUE),
            cv_co2 = sd(NEE,na.rm=TRUE)/mean(NEE,na.rm=TRUE)*100,
            min_ch4 = min(CH4,na.rm = TRUE),
            max_ch4 = max(CH4,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            mean_ch4 = mean(CH4,na.rm=TRUE),
            sd_ch4 = sd(CH4,na.rm=TRUE),
            cv_ch4 = sd(CH4,na.rm=TRUE)/mean(CH4,na.rm=TRUE)*100)

## Combine all stats into one (EC and Diffusive fluxes) - Table S5 and S6
fcr_stats <- rbind(fcr_hourly_stats,fcr_daily_stats,fcr_weekly_stats,fcr_monthly_stats)
write_csv(fcr_stats,"./Data/20211222_SI5_ECStats.csv")

## Plot Daily means (CO2, CH4; EC + Diffusive) and Monthly means
## Figure 2
dailynee <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="EC"),size = 1) +
  geom_errorbar(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  geom_point(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  scale_color_manual(breaks=c("EC","Diff"),
                      values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

dailych4 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="EC"),size = 1) +
  geom_errorbar(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  geom_point(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  scale_color_manual(breaks=c("EC","Diff"),
                     values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_month <- ggplot(fcr_monthly) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(yearmon, ymin = NEE - NEE_sd, ymax = NEE + NEE_sd), width = 0.2) +
  geom_point(aes(yearmon, NEE), color = "#E63946", size=4) +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("") +
  theme_classic(base_size=15) +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(colour = 'black'))

ch4_month <- ggplot(fcr_monthly) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(yearmon, ymin = CH4 - CH4_sd, ymax = CH4 + CH4_sd), width = 0.2) +
  geom_point(aes(yearmon, CH4), color = "#E63946", size=4) +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("") +
  theme_classic(base_size=15) +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(colour = 'black'))

ggarrange(dailynee,dailych4,
          ggarrange(co2_month, ch4_month,nrow=1,ncol=2, labels = c("C.","D."), font.label = list(face="plain",size=15)),
          nrow=3,ncol=1,labels = c("A.","B."),
          font.label = list(face="plain",size=15),common.legend = TRUE)

ggsave("./Fig_Output/Figure2.png",width = 8, height=11, units="in",dpi=320)

## Plot Hourly and Weekly time series for SI
## Hourly: Figure S3
co2_hourly <- ggplot(fcr_hourly) +
  geom_ribbon(mapping=aes(x=DateTime,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(DateTime, NEE),color="#E63946",size = 1) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
  theme_classic(base_size = 15)

ch4_hourly <- ggplot(fcr_hourly) +
  geom_ribbon(mapping=aes(x=DateTime,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(DateTime, CH4),color="#E63946",size = 1) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
  theme_classic(base_size = 15)

ggarrange(co2_hourly,ch4_hourly,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI3_Hourly_time.jpg",width = 8, height=8, units="in",dpi=320)

## Weekly: Figure S4
co2_week <- ggplot(fcr_weekly)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=Date,y=NEE),color="#E63946",size = 1)+
  geom_point(mapping=aes(x=Date,y=NEE),color="#E63946",size=2)+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size=15)

ch4_week <- ggplot(fcr_weekly)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=Date,y=CH4),color="#E63946",size = 1)+
  geom_point(mapping=aes(x=Date,y=CH4),color="#E63946",size=2)+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size=15)

ggarrange(co2_week,ch4_week,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI4_Weekly_time.jpg",width = 8, height=8, units="in",dpi=320)

## Plot Diffusive vs. EC fluxes for each time point
# First combine EC fluxes and diffusive fluxes by time: combine by the hour
diff_hour_co2 <- ghg_fluxes_avg_co2 %>% 
  mutate(DateTime = format(round(DateTime, units="hours"), format="%Y-%m-%d %H:%M")) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
  rename(co2_mean = mean, co2_sd = sd) %>% 
  select(DateTime,co2_mean,co2_sd)

diff_hour_ch4 <- ghg_fluxes_avg_ch4 %>% 
  mutate(DateTime = format(round(DateTime, units="hours"), format="%Y-%m-%d %H:%M")) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
  rename(ch4_mean = mean, ch4_sd = sd) %>% 
  select(DateTime,ch4_mean,ch4_sd)

diff_hour_all <- left_join(diff_hour_co2,diff_hour_ch4,by="DateTime")

ec_diff <- left_join(diff_hour_all,fcr_hourly,by="DateTime")

# Select unique rows
ec_diff <- ec_diff %>% 
  distinct()

# Plot EC vs. Diffusive fluxes for each half-hour
# Figure S7
diff_ec_co2 <- ggplot(ec_diff,mapping=aes(x=co2_mean,y=NEE))+
  geom_abline(intercept = 0)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_errorbar(mapping=aes(x = co2_mean, y = NEE, xmin = co2_mean-co2_sd, xmax = co2_mean+co2_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_errorbar(mapping=aes(x=co2_mean,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_point(mapping=aes(x=co2_mean,y=NEE))+
  ylab(expression(~EC~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab(expression(~Mean~Diff~CO[2]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

diff_ec_ch4 <- ggplot(ec_diff)+
  geom_abline(intercept = 0)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_errorbar(mapping=aes(x = ch4_mean, y = CH4, xmin = ch4_mean-ch4_sd, xmax = ch4_mean+ch4_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_errorbar(mapping=aes(x=ch4_mean,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_point(mapping=aes(x=ch4_mean,y=CH4))+
  ylab(expression(~EC~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab(expression(~Mean~Diff~CH[4]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

ggarrange(diff_ec_co2,diff_ec_ch4,ncol=2,nrow=1,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Diff_EC.jpg",width = 8, height=4, units="in",dpi=320)

### Check high wind days prior to turnover (for CH4 release from the metalimnion)
## Aggregate met data to daily
met_daily <- met_30_2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(wind = mean(WindSpeed_Average_m_s,na.rm=TRUE),
                   wind_sd = sd(WindSpeed_Average_m_s,na.rm=TRUE))

met_daily$Date <- as.POSIXct(paste(met_daily$Year, met_daily$Month, met_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Plot wind + CH4 fluxes around this time period
ch4_turnover <- ggplot(fcr_daily) +
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4),color="#E63946",size = 1) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-10-01"),as.POSIXct("2020-11-30"))+
  theme_classic(base_size = 15)

wind_turnover <- ggplot(met_daily) +
  geom_ribbon(mapping=aes(x=Date,y=wind,ymin=wind-wind_sd,ymax=wind+wind_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date,wind),color="#E63946",size = 1) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(met_30_2,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s),alpha = 0.1)+
  xlab("") +
  ylab(expression(Wind~Speed~(m~s^-1))) +
  xlim(as.POSIXct("2020-10-01"),as.POSIXct("2020-11-30"))+
  theme_classic(base_size = 15)

ggarrange(ch4_turnover,wind_turnover,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI10_Turnover.jpg",width = 8, height=8, units="in",dpi=320)

### Plotting diel variations: select data from noon and midnight ----
## Noon = 11-1 pm; Midnight = 23-1 am
diel_flux <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f,NEE_uStar_orig,ch4_flux_uStar_orig,u) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  filter(Hour %in% c(11,12,13)|Hour %in% c(23,0,1)) %>% 
  mutate(diel = ifelse(Hour == 11 | Hour == 12 | Hour == 13, "Day","Night"))

diel_agg <- diel_flux %>% 
  ungroup() %>% 
  group_by(Year,Month,Day,diel) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE),
                   u = mean(u, na.rm = TRUE),
                   u_sd = sd(u, na.rm = TRUE),
                   Hour = mean(Hour)) %>% 
  mutate(season = ifelse(Month == 12 | Month == 1 | Month == 2, "Winter",
                         ifelse(Month == 3 | Month == 4 | Month == 5, "Spring",
                                ifelse(Month == 6 | Month == 7 | Month == 8 | Month == 9, "Summer",
                                       "Fall"))))

diel_agg$Date <- as.POSIXct(paste(diel_agg$Year, diel_agg$Month, diel_agg$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Calculate mean difference between day and night fluxes
diel_agg %>% 
  select(diel,NEE, CH4, u) %>% 
  group_by(diel) %>% 
  summarise_all(median,na.rm=TRUE)

diel_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

## Look at dawn vs. dusk differences in fluxes
dawn_flux <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f,NEE_uStar_orig,ch4_flux_uStar_orig,u) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  filter(Hour %in% c(5,6,7)|Hour %in% c(17,18,19)) %>% 
  mutate(diel = ifelse(Hour == 5 | Hour == 6 | Hour == 7, "Dawn","Dusk"))

dawn_agg <- dawn_flux %>% 
  ungroup() %>% 
  group_by(Year,Month,Day,diel) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE),
                   u = mean(u, na.rm = TRUE),
                   u_sd = sd(u, na.rm = TRUE),
                   Hour = mean(Hour)) %>% 
  mutate(season = ifelse(Month == 12 | Month == 1 | Month == 2, "Winter",
                         ifelse(Month == 3 | Month == 4 | Month == 5, "Spring",
                                ifelse(Month == 6 | Month == 7 | Month == 8 | Month == 9, "Summer",
                                       "Fall"))))

dawn_agg$Date <- as.POSIXct(paste(dawn_agg$Year, dawn_agg$Month, dawn_agg$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Calculate mean difference between dawn and dusk fluxes
dawn_agg %>% 
  select(diel,NEE, CH4, u) %>% 
  group_by(diel) %>% 
  summarise_all(median,na.rm=TRUE)

dawn_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

## Calculate statistics for Supplementary Table (Fig S7) - Day vs. Night
diel_stats <- diel_agg %>% 
  ungroup() %>% 
  group_by(diel) %>% 
  summarise(p25_nee = quantile(NEE,0.25,na.rm=TRUE),
            med_nee = median(NEE,na.rm=TRUE),
            p75_nee = quantile(NEE,0.75,na.rm=TRUE),
            p25_ch4 = quantile(CH4,0.25,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            p75_ch4 = quantile(CH4,0.75,na.rm=TRUE),
            p25_wind = quantile(u,0.25,na.rm=TRUE),
            med_wind = median(u,na.rm=TRUE),
            p75_wind = quantile(u,0.75,na.rm=TRUE))

## Save for supplementary table (Table S7)
write.csv(diel_stats,"./Fig_output/20211223_TableSx_DielStats.csv")

## Calculate statistics for Supplementary Table (Fig S7) - Dawn vs. Dusk
dawn_stats <- dawn_agg %>% 
  ungroup() %>% 
  group_by(diel) %>% 
  summarise(p25_nee = quantile(NEE,0.25,na.rm=TRUE),
            med_nee = median(NEE,na.rm=TRUE),
            p75_nee = quantile(NEE,0.75,na.rm=TRUE),
            p25_ch4 = quantile(CH4,0.25,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            p75_ch4 = quantile(CH4,0.75,na.rm=TRUE),
            p25_wind = quantile(u,0.25,na.rm=TRUE),
            med_wind = median(u,na.rm=TRUE),
            p75_wind = quantile(u,0.75,na.rm=TRUE))

## Save for supplementary table
write.csv(dawn_stats,"./Fig_output/20220106_TableSx_DawnStats.csv")

## Calculate Paired-Wilcoxon signed rank tests
daynight_wide_co2_f <- diel_agg %>% 
  ungroup() %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(daynight_wide_co2_f$Day,daynight_wide_co2_f$Night,paired=TRUE)

daynight_wide_ch4_f <- diel_agg %>% 
  ungroup() %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na()

wilcox.test(daynight_wide_ch4_f$Day,daynight_wide_ch4_f$Night,paired=TRUE)

daynight_wide_u <- diel_agg %>% 
  ungroup() %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na()

wilcox.test(daynight_wide_u$Day,daynight_wide_u$Night,paired=TRUE)

## Calculate Dawn/Dusk differences
dawndusk_wide_co2_f <- dawn_agg %>% 
  ungroup() %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(dawndusk_wide_co2_f$Dawn,dawndusk_wide_co2_f$Dusk,paired=TRUE)

dawndusk_wide_ch4_f <- dawn_agg %>% 
  ungroup() %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na()

wilcox.test(dawndusk_wide_ch4_f$Dawn,dawndusk_wide_ch4_f$Dusk,paired=TRUE)

dawndusk_wide_u <- dawn_agg %>% 
  ungroup() %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na()

wilcox.test(dawndusk_wide_u$Dawn,dawndusk_wide_u$Dusk,paired=TRUE)

## Plot diel (day/night) and dawn/dusk comparisons: Figure 3
diel_co2 <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label ="p = 0.02*
n = 162",
           x = "Day", hjust = 0, y = 18, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  ylim(-20,20)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

diel_ch4 <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "p = 0.43
n = 123",
           x = "Day", hjust = 0, y = 0.07, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 

diel_wind <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=u,color=diel))+
  annotate("text",label = "p < 0.005*
n = 245",
           x = "Day", hjust = 0, y = 4.5, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("Wind speed (m s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

dawn_co2 <- dawn_agg %>% 
  ggplot(mapping=aes(x=diel,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label ="p < 0.001*
n = 125",
           x = "Dawn", hjust = 0, y = 18, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  ylim(-20,20)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

dawn_ch4 <- dawn_agg %>% 
  ggplot(mapping=aes(x=diel,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "p = 0.20
n = 88",
           x = "Dawn", hjust = 0, y = 0.045, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

dawn_wind <- dawn_agg %>% 
  ggplot(mapping=aes(x=diel,y=u,color=diel))+
  annotate("text",label = "p < 0.001*
n = 172",
           x = "Dawn", hjust = 0, y = 3.75, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("Wind speed (m s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

ggarrange(diel_co2,diel_ch4,diel_wind,dawn_co2,dawn_ch4,dawn_wind,nrow=2,ncol=3,
          labels=c("A.","B.","C.","D.","E.","F."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Figure3.png",width = 9, height=7.5, units="in",dpi=320)

# Calculate and plot diel temperature and wind - Figure S10
hour_temp <- catwalk_30 %>% 
  mutate(Hour = hour(DateTime),
         Temp_C_surface_2 = ThermistorTemp_C_surface) %>% 
  group_by(Hour) %>% 
  select(ThermistorTemp_C_surface, Temp_C_surface_2, Hour) %>% 
  summarise(Temp_C_surface = mean(ThermistorTemp_C_surface,na.rm=TRUE),
            Temp_C_surface_sd = sd(Temp_C_surface_2,na.rm=TRUE)) %>% 
  ggplot(mapping=aes(x=Hour,y=Temp_C_surface))+
  geom_vline(xintercept = 12, linetype = "dashed")+
  geom_line(size=1)+
  xlab("Hour") + 
  ylab(expression(Temp~(C^o)))+
  theme_classic(base_size=15)

hour_wind <- met_30 %>% 
  mutate(Hour = hour(DateTime),
         wnd_speed_2 = WindSpeed_Average_m_s) %>% 
  group_by(Hour) %>% 
  select(Hour,WindSpeed_Average_m_s,wnd_speed_2) %>% 
  summarise(mean_wnd = mean(WindSpeed_Average_m_s,na.rm=TRUE),
            sd_wnd = sd(wnd_speed_2,na.rm=TRUE)) %>% 
  ggplot(mapping=aes(x=Hour,y=mean_wnd))+
  geom_ribbon(mapping=aes(x=Hour,y=mean_wnd,ymin=mean_wnd-sd_wnd,ymax=mean_wnd+sd_wnd),fill="grey",alpha=0.5)+
  geom_vline(xintercept = 12, linetype = "dashed")+
  geom_line(size=1)+
  xlab("Hour") + 
  ylab(expression(Wind~Speed~(m~s^-2)))+
  theme_classic(base_size=15)

ggarrange(hour_temp,hour_wind,nrow=1,ncol=2,
          labels=c("A.","B."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI9_DielTempWnd.jpg",width = 8, height=4, units="in",dpi=320)

### Think about winter variability (especially with ice!)
## Download Ice data from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/456/3/ddb604975b844125b3bd542862ca0a3b" 
#infile1 <- paste0(getwd(),"/Data/Ice_Data.csv")
#download.file(inUrl1,infile1,method="curl")

ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))
ice_on <- ice %>% 
  filter(Date>"2020-01-01" & IceOn == 1)
ice_off <- ice %>% 
  filter(Date>"2020-01-01" & IceOff == 1)

## Calculate average flux for each ice on/off period
ice_off_1 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-19 19:00:00" & DateTime < "2020-12-26 19:00:00") %>% 
  select(DateTime,NEE_uStar_orig,NEE_uStar_f,ch4_flux_uStar_orig,ch4_flux_uStar_f) %>% 
  mutate(ice_period = "Off_1") %>% 
  mutate(ice = "off") %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 7 & Hour <= 19, "Day",
                       ifelse(Hour >= 19 | Hour <= 7, "Night", NA)))

wilcox.test(NEE_uStar_orig ~ diel, data=ice_off_1)
wilcox.test(ch4_flux_uStar_orig ~ diel, data=ice_off_1, correct = FALSE)

ice_on_1 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-26 19:00:00" & DateTime < "2020-12-29 19:00:00")%>% 
  select(DateTime,NEE_uStar_orig,NEE_uStar_f,ch4_flux_uStar_orig,ch4_flux_uStar_f) %>% 
  mutate(ice_period = "On_1") %>% 
  mutate(ice = "on") %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 7 & Hour <= 19, "Day",
                       ifelse(Hour >= 19 | Hour <= 7, "Night", NA)))

wilcox.test(NEE_uStar_orig ~ diel, data=ice_on_1)
wilcox.test(ch4_flux_uStar_orig ~ diel, data=ice_on_1)

ice_off_2 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-29 19:00:00" & DateTime < "2021-01-09 19:00:00")%>% 
  select(DateTime,NEE_uStar_orig,NEE_uStar_f,ch4_flux_uStar_orig,ch4_flux_uStar_f) %>% 
  mutate(ice_period = "Off_2") %>% 
  mutate(ice = "off") %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 7 & Hour <= 19, "Day",
                       ifelse(Hour >= 19 | Hour <= 7, "Night", NA)))

wilcox.test(NEE_uStar_orig ~ diel, data=ice_off_2)
wilcox.test(ch4_flux_uStar_orig ~ diel, data=ice_off_2)

ice_on_2 <- eddy_flux %>% 
  filter(DateTime >= "2021-01-09 19:00:00" & DateTime < "2021-02-09 19:00:00")%>% 
  select(DateTime,NEE_uStar_orig,NEE_uStar_f,ch4_flux_uStar_orig,ch4_flux_uStar_f) %>% 
  mutate(ice_period = "On_2") %>% 
  mutate(ice = "on") %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 7 & Hour <= 19, "Day",
                       ifelse(Hour >= 19 | Hour <= 7, "Night", NA)))

wilcox.test(NEE_uStar_orig ~ diel, data=ice_on_2)
wilcox.test(ch4_flux_uStar_orig ~ diel, data=ice_on_2)

ice_all <- rbind(ice_off_1,ice_on_1,ice_off_2,ice_on_2)

# Test for significance for all ice on vs. ice off
wilcox.test(NEE_uStar_orig ~ ice, data=ice_all)
# Not statistically significant (p = 0.6202)
wilcox.test(ch4_flux_uStar_orig ~ ice, data=ice_all)
# Not statistically significant (p = 0.1026)

# Ice_1
ice_1 <- ice_all %>% 
  filter(ice_period %in% c("Off_1","On_1"))

wilcox.test(NEE_uStar_orig ~ ice, data=ice_1)
wilcox.test(ch4_flux_uStar_orig ~ ice, data=ice_1)

# Ice_2
ice_2 <- ice_all %>% 
  filter(ice_period %in% c("Off_2","On_2"))

wilcox.test(NEE_uStar_orig ~ ice, data=ice_2)
wilcox.test(ch4_flux_uStar_orig ~ ice, data=ice_2)

# Ice stats
ice_stats <- ice_all %>% 
  group_by(ice_period,diel) %>% 
  summarise(p25_nee = quantile(NEE_uStar_orig,0.25,na.rm=TRUE),
            med_nee = median(NEE_uStar_orig,na.rm=TRUE),
            p75_nee = quantile(NEE_uStar_orig,0.75,na.rm=TRUE),
            p25_ch4 = quantile(ch4_flux_uStar_orig,0.25,na.rm=TRUE),
            med_ch4 = median(ch4_flux_uStar_orig,na.rm=TRUE),
            p75_ch4 = quantile(ch4_flux_uStar_orig,0.75,na.rm=TRUE))

ice_stats_2 <- ice_all %>% 
  ungroup() %>% 
  group_by(ice_period) %>% 
  summarise(p25_nee = quantile(NEE_uStar_orig,0.25,na.rm=TRUE),
            med_nee = median(NEE_uStar_orig,na.rm=TRUE),
            p75_nee = quantile(NEE_uStar_orig,0.75,na.rm=TRUE),
            p25_ch4 = quantile(ch4_flux_uStar_orig,0.25,na.rm=TRUE),
            med_ch4 = median(ch4_flux_uStar_orig,na.rm=TRUE),
            p75_ch4 = quantile(ch4_flux_uStar_orig,0.75,na.rm=TRUE)) %>% 
  mutate(diel = "All") %>% 
  relocate(diel,.before = p25_nee)

ice_stats_all <- ice_all %>% 
  ungroup() %>% 
  group_by(ice) %>% 
  summarise(p25_nee = quantile(NEE_uStar_orig,0.25,na.rm=TRUE),
            med_nee = median(NEE_uStar_orig,na.rm=TRUE),
            p75_nee = quantile(NEE_uStar_orig,0.75,na.rm=TRUE),
            p25_ch4 = quantile(ch4_flux_uStar_orig,0.25,na.rm=TRUE),
            med_ch4 = median(ch4_flux_uStar_orig,na.rm=TRUE),
            p75_ch4 = quantile(ch4_flux_uStar_orig,0.75,na.rm=TRUE)) %>% 
  rename(ice_period = ice) %>% 
  mutate(diel = "All") %>% 
  relocate(diel,.before = p25_nee)

ice_stats_final <- rbind(ice_stats,ice_stats_2,ice_stats_all)

write.csv(ice_stats_final,"./Fig_output/20211222_TableS6_Ice_stats.csv")

## Daily means in winter (ice on/off): Figure 4
winter_co2 <- ggplot(fcr_daily)+
  annotate(geom="text",x = as.POSIXct("2020-12-24"),y = 9,label = "Off")+
  annotate(geom="text",x=as.POSIXct("2020-12-28 12:00"),y=9,label="On")+
  annotate(geom="text",x=as.POSIXct("2021-01-04"),y=9,label="Off")+
  annotate(geom="text",x=as.POSIXct("2021-01-25"),y=9,label="On")+
  annotate(geom="text",x=as.POSIXct("2021-02-11"),y=9,label="Off")+
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_point(eddy_flux, mapping = aes(DateTime, NEE_uStar_orig),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping = aes(x = Date, y = NEE, ymin = NEE-NEE_sd,ymax = NEE+NEE_sd),fill="#E63946",alpha=0.4)+
  geom_line(mapping = aes(Date, NEE),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-02-11"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-10,10)+
  theme_classic(base_size = 15)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

ice_all$ice_period <- factor(ice_all$ice_period, levels=c("Off_1", "On_1", "Off_2", "On_2"))

ice_co2 <- ggplot(ice_all,mapping=aes(x=ice_period,y=NEE_uStar_orig,color=ice_period))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA,size=0.7)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Off_1','On_1','Off_2','On_2'),values=c("#008080","#477CC2","#008080","#477CC2"))+
  theme_classic(base_size = 15)+
  ylim(-10,10)+
  theme(legend.position="none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

winter_ch4 <- ggplot(fcr_daily)+
  annotate(geom="text",x = as.POSIXct("2020-12-24"),y = 0.025,label = "Off")+
  annotate(geom="text",x=as.POSIXct("2020-12-28 12:00"),y=0.025,label="On")+
  annotate(geom="text",x=as.POSIXct("2021-01-04"),y=0.025,label="Off")+
  annotate(geom="text",x=as.POSIXct("2021-01-25"),y=0.025,label="On")+
  annotate(geom="text",x=as.POSIXct("2021-02-11"),y=0.025,label="Off")+
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dotted", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_point(eddy_flux, mapping = aes(DateTime, ch4_flux_uStar_orig),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping = aes(x = Date, y = CH4, ymin = CH4-CH4_sd,ymax = CH4+CH4_sd),fill="#E63946",alpha=0.4)+
  geom_line(mapping = aes(Date, CH4),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-02-11"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-0.03,0.03)+
  theme_classic(base_size = 15)

ice_all$ice_period <- factor(ice_all$ice_period, levels=c("Off_1", "On_1", "Off_2", "On_2"))

ice_ch4 <- ggplot(ice_all,mapping=aes(x=ice_period,y=ch4_flux_uStar_orig,color=ice_period))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA,size=0.7)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Off_1','On_1','Off_2','On_2'),values=c("#008080","#477CC2","#008080","#477CC2"))+
  ylim(-0.03,0.03)+
  theme_classic(base_size = 15)+
  theme(legend.position="none")+
  theme(legend.title=element_blank())

ggarrange(winter_co2,ice_co2,winter_ch4,ice_ch4,nrow=2,ncol=2,common.legend = FALSE,
          labels=c("A.","B.","C.","D."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Figure4_Ice_on_off.png",width = 8, height=6, units="in",dpi=320)

### Plotting cumulative Co2 and Ch4 throughout the study period ----
## From 2020-04-04 to 2021-04-03
eddy_flux_sum <- eddy_flux %>% 
  select(DateTime,ch4_flux_uStar_orig,ch4_flux_uStar_f,ch4_flux_uStar_fsd,ch4_flux_uStar_fall,ch4_flux_uStar_fqc,
         ch4_flux_U97.5_f,NEE_uStar_orig,NEE_uStar_f,NEE_uStar_fsd,NEE_uStar_fall,NEE_uStar_fqc,NEE_U97.5_f) %>%
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v)))

## Compare cumulative fluxes: Figure 7
sum_co2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(eddy_flux_sum,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,color="EC"),size=1, color="#E63946")+
  geom_ribbon(eddy_flux_sum,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,ymin=co2_sum_g_m2_d-co2_sum_sd,ymax=co2_sum_g_m2_d+co2_sum_sd),fill="#E63946",alpha=0.3)+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)+
  labs(color="")

sum_ch4 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(eddy_flux_sum,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,color="EC"),size=1, color="#E63946")+
  geom_ribbon(eddy_flux_sum,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,ymin=ch4_sum_g_m2_d-ch4_sum_sd,ymax=ch4_sum_g_m2_d+ch4_sum_sd),fill="#E63946",alpha=0.3)+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)+
  labs(color="")

ggarrange(sum_co2,sum_ch4,nrow=1,ncol=2,common.legend = TRUE,labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Figure7.png",width = 9, height=4.5, units="in",dpi=320)

## Calculate cumulative fluxes in summer (Jun, Jul, Aug, Sep) vs. winter (Dec, Jan, Feb, Mar)
## Summer (Jun - Sept) = 376.1 g C-CO2 m2 d; 0.7274 g C-CH4 m2 d
eddy_flux_sum %>% 
  filter(DateTime == as.POSIXct("2020-05-31 20:00:00") | DateTime == as.POSIXct("2020-09-30 20:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d)

# Winter (Dec - Mar) = 60 g C-CO2 m2 d; 0.023 g C-CH4 m2 d
eddy_flux_sum %>% 
  filter(DateTime == as.POSIXct("2020-12-31 20:00:00") | DateTime == as.POSIXct("2021-03-31 20:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d)