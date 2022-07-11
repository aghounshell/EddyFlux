### Scrip to calculate diffusive fluxes
### Following: 1_DiffusiveFluxes_ECGraphs.R

### 19 May 2022, A. Hounshell

### For re-submission of EC MS

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples,viridis)

###############################################################################

### First load in wind data from Met station at FCR ----
## Downloaded from EDI: 06 May 2022
## Included preliminary met data for 2022 downloaded from the Gateway

#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

met_edi <- read.csv("./Data/Met_final_2015_2021.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

# Load in initial 2022 Met data (from 1a_Rev_2022_MetData.R)
met_2022 <- read.csv("./Data/FCR_Met_final_2022.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

# Combine Met Data together
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
  filter(DateTime_Adj >= as.POSIXct("2020-05-01 20:00:00") & DateTime_Adj < as.POSIXct("2022-04-30 20:00:00"))

names(met_30_2)[names(met_30_2) == 'DateTime_Adj'] <- 'DateTime'

###############################################################################

## Load in GHG data
## Downloaded 2022 EDI staged data
#inUrl1  <- "https://pasta-s.lternet.edu/package/data/eml/edi/928/3/38d72673295864956cccd6bbba99a1a3" 
#infile1 <- paste0(getwd(),"/Data/final_GHG_2015_May2022.csv")
#download.file(inUrl1,infile1,method="curl")

## Will also want to include GHG data from 2022 - forthcoming!
ghg <- read.csv("./Data/final_GHG_2015_May2022.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) %>% 
  filter(Reservoir == "FCR" & Site == 50 & Depth_m == 0.1) %>% 
  filter(DateTime >= "2020-01-01") %>% 
  mutate(DateTime = round_date(DateTime, "30 mins")) 

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

## Pivot wider
ghg_fluxes <- ghg %>% 
  select(DateTime,Rep,ch4_umolL,co2_umolL) %>% 
  pivot_wider(names_from = Rep, values_from = c('ch4_umolL','co2_umolL'), values_fn = mean)

###############################################################################

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

# Calculate average half-hourly values
# Select data every 30 minutes from Jan 2020 to end of met data
catwalk_all$Breaks <- cut(catwalk_all$DateTime,breaks = "30 mins",right=FALSE)
catwalk_all$Breaks <- ymd_hms(as.character(catwalk_all$Breaks))

# Average to half-hourly catwalk measurements
catwalk_30 <- catwalk_all %>% 
  select(DateTime,ThermistorTemp_C_surface,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(mean,na.rm=TRUE)

catwalk_30 <- catwalk_30 %>% 
  select(Breaks,ThermistorTemp_C_surface) %>% 
  mutate(Breaks = force_tz(Breaks,"EST"))

names(catwalk_30)[names(catwalk_30) == 'Breaks'] <- 'DateTime'

# Combine w/ met data to ensure even timeseries
all_data <- left_join(met_30_2,catwalk_30,by="DateTime")

###############################################################################

# Use PAR profiles to estimate Kd for the entire year: from 2013-2021

#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
#infile1 <- paste0(getwd(),"/Data/CTD_final_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

ctd <- read.csv("./Data/CTD_final_2013_2021.csv") %>% 
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

###############################################################################

### Format for analysis in LakeMetabolizer
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

## Then match with GHG data and calculate K for each time period
ghg_match <- left_join(ghg_fluxes,all_data,by="DateTime") %>% 
  select(-ch4_umolL_1,-ch4_umolL_2,-co2_umolL_1,-co2_umolL_2)

## Then calculate K
# Define variables following convention used in LakeMetabolizer
wnd.z <- rep(3,66) # in m
Kd <- rep(0.402,66) # units 1/m
lat <- rep(37.30,66) # in degrees N
lake.area <- rep(0.119*1000000,66) # in m2
atm.press <- ghg_match$BP_Average_kPa*10 # in millibar
dateTime <- ghg_match$DateTime
Ts <- ghg_match$ThermistorTemp_C_surface # in deg C
z.aml <- rep(0.15,66) # set to 0.15 following Heiskanen method
airT <- ghg_match$AirTemp_Average_C # in deg C
wnd <- ghg_match$WindSpeed_Average_m_s # in m/s
RH <- ghg_match$RH_percent # percent
sw <- ghg_match$ShortwaveRadiationDown_Average_W_m2 # downwelling shortwave radiation in w/m2

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
k600_read <- rep(0,66)

for (i in 1:length(wnd)){
  k600_read[i] <- k.read.base(wnd.z[i], Kd[i], lat[i], lake.area[i], atm.press[i], dateTime[i], Ts[i], z.aml[i],
                              airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

k600_soloviev <- rep(0,66)

for (i in 1:length(wnd)){
  k600_soloviev[i] <- k.read.soloviev.base(wnd.z[i], Kd[i], lat[i], lake.area[i], atm.press[i], dateTime[i], Ts[i], z.aml[i],
                                           airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

k600_macIntyre <- rep(0,66)

for (i in 1:length(wnd)){
  k600_macIntyre[i] <- k.macIntyre.base(wnd.z[i], Kd[i], atm.press[i], dateTime[i], Ts[i], z.aml[i], airT[i], wnd[i], RH[i], sw[i],
                                        lwnet[i], params=c(1.2,0.4872,1.4784))
}

k600_heiskanen <- rep(0,66)

for (i in 1:length(wnd)){
  k600_heiskanen[i] <- k.heiskanen.base(wnd.z[i], Kd[i], atm.press[i], dateTime[i], Ts[i], z.aml[i], airT[i], wnd[i], RH[i], sw[i], lwnet[i])
}

## Then convert all k600 values by temperature and gas following Lake Metabolizer
# Correct k600 for CO2
k600_cole_co2 <- rep(0,66)
k600_crusius_co2 <- rep(0,66)
k600_vachon_co2 <- rep(0,66)
k600_read_co2 <- rep(0,66)
k600_soloviev_co2 <- rep(0,66)
k600_macIntyre_co2 <- rep(0,66)
k600_heiskanen_co2 <- rep(0,66)

k600_cole_co2 <- k600.2.kGAS.base(k600_cole,Ts,gas="CO2")
k600_crusius_co2 <- k600.2.kGAS.base(k600_crusius,Ts,gas="CO2")
k600_vachon_co2 <- k600.2.kGAS.base(k600_vachon,Ts,gas="CO2")
k600_read_co2 <- k600.2.kGAS.base(k600_read,Ts,gas="CO2")
k600_soloviev_co2 <- k600.2.kGAS.base(k600_soloviev,Ts,gas="CO2")
k600_macIntyre_co2 <- k600.2.kGAS.base(k600_macIntyre,Ts,gas="CO2")
k600_heiskanen_co2 <- k600.2.kGAS.base(k600_heiskanen,Ts,gas="CO2")

# Correct k600 for CH4
k600_cole_ch4 <- rep(0,66)
k600_crusius_ch4 <- rep(0,66)
k600_vachon_ch4 <- rep(0,66)
k600_read_ch4 <- rep(0,66)
k600_soloviev_ch4 <- rep(0,66)
k600_macIntyre_ch4 <- rep(0,66)
k600_heiskanen_ch4 <- rep(0,66)

k600_cole_ch4 <- k600.2.kGAS.base(k600_cole,Ts,gas="CH4")
k600_crusius_ch4 <- k600.2.kGAS.base(k600_crusius,Ts,gas="CH4")
k600_vachon_ch4 <- k600.2.kGAS.base(k600_vachon,Ts,gas="CH4")
k600_read_ch4 <- k600.2.kGAS.base(k600_read,Ts,gas="CH4")
k600_soloviev_ch4 <- k600.2.kGAS.base(k600_soloviev,Ts,gas="CH4")
k600_macIntyre_ch4 <- k600.2.kGAS.base(k600_macIntyre,Ts,gas="CH4")
k600_heiskanen_ch4 <- k600.2.kGAS.base(k600_heiskanen,Ts,gas="CH4")

# Aggregate CO2 and CH4 k600 into a single data frame
k600_corr <- data.frame(matrix(ncol = 15, nrow = 66))

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

## Save k600 corrected values
write.csv(k600_corr,"./Data/20220617_k600_corr.csv", row.names = FALSE)

###############################################################################

## Merge k and GHG data
ghg_fluxes <- left_join(ghg_fluxes,k600_corr,by="DateTime")

## Use EC concentrations (CO2 and CH4) for atmospheric concentrations
## Load in EC data from EDI: https://pasta-s.lternet.edu/package/data/eml/edi/920/2/9e658ef44de05303dbc496fc25e8c49a
ec <- read.csv("./Data/20220506_EddyPro_cleaned.csv")

ec <- ec %>% 
  mutate(DateTime = as.POSIXct(paste(date, time), format="%Y-%m-%d %H:%M:%S", tz="EST"))

## Plot EC concentrations (just to see!)
## Methane
ggplot(ec,mapping=aes(DateTime, co2_mole_fraction_umolmol))+
  geom_line()

## Carbon dioxide
ggplot(ec,mapping=aes(DateTime, ch4_mole_fraction_umolmol))+
  geom_line()

## Do some light QA/QC'ing
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

# Make sure there is a continous timeseries
ts <- seq.POSIXt(as.POSIXct("2020-04-05 00:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), 
                 as.POSIXct("2022-05-02 13:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), by = "30 min")
ts2 <- data.frame(DateTime = ts)

# Join Eddy Flux data with list of dates+time
ec_2 <- left_join(ts2, ec_2, by = 'DateTime')

ec_2 <- ec_2 %>% 
  mutate(co2_mole_fraction = na.fill(na.approx(co2_mole_fraction,na.rm=FALSE),"extend")) %>% 
  mutate(ch4_mole_fraction = na.fill(na.approx(ch4_mole_fraction,na.rm=FALSE),"extend"))

ghg_fluxes_2 <- left_join(ghg_fluxes,ec_2,by="DateTime")

# Calculate fluxes: in umol/m2/s using each k600 method
ghg_fluxes_2 <- ghg_fluxes_2 %>% 
  mutate(nv = (BP_Average_kPa*0.00986923/(0.0820573660809596*(AirTemp_Average_C + 273.15)))) %>% # units = mols/L
  mutate(ch4_umolL_atm = ch4_mole_fraction/1e6*nv*1e6) %>% # units = umol/L
  mutate(ch4_flux_k600_Cole_1 = 1000*k600_Cole_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Cole_2 = 1000*k600_Cole_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Crusius_1 = 1000*k600_Crusius_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Crusius_2 = 1000*k600_Crusius_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Vachon_1 = 1000*k600_Vachon_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Vachon_2 = 1000*k600_Vachon_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Read_1 = 1000*k600_Read_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Read_2 = 1000*k600_Read_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Soloviev_1 = 1000*k600_Soloviev_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Soloviev_2 = 1000*k600_Soloviev_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_MacIntyre_1 = 1000*k600_MacIntyre_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_MacIntyre_2 = 1000*k600_MacIntyre_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Heiskanen_1 = 1000*k600_Heiskanen_co2*(ch4_umolL_1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_k600_Heiskanen_2 = 1000*k600_Heiskanen_co2*(ch4_umolL_2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_umolL_atm = co2_mole_fraction/1e6*nv*1e6) %>%  # units = umol/L
  mutate(co2_flux_k600_Cole_1 = 1000*k600_Cole_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Cole_2 = 1000*k600_Cole_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Crusius_1 = 1000*k600_Crusius_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Crusius_2 = 1000*k600_Crusius_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Vachon_1 = 1000*k600_Vachon_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Vachon_2 = 1000*k600_Vachon_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Read_1 = 1000*k600_Read_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Read_2 = 1000*k600_Read_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Soloviev_1 = 1000*k600_Soloviev_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Soloviev_2 = 1000*k600_Soloviev_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_MacIntyre_1 = 1000*k600_MacIntyre_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_MacIntyre_2 = 1000*k600_MacIntyre_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Heiskanen_1 = 1000*k600_Heiskanen_co2*(co2_umolL_1-co2_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_flux_k600_Heiskanen_2 = 1000*k600_Heiskanen_co2*(co2_umolL_2-co2_umolL_atm)/24/60/60)  # units = umol C/m2/s

fluxes_rep1 <- ghg_fluxes_2 %>% 
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

fluxes_rep2 <- ghg_fluxes_2 %>% 
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

names(flux_avg_co2)[names(flux_avg_co2) == 'mean'] <- 'co2_mean_umol_m2_s'
names(flux_avg_co2)[names(flux_avg_co2) == 'sd'] <- 'co2_sd_umol_m2_s'

flux_avg_ch4 <- fluxes_long_ch4 %>% 
  select(DateTime,ch4_flux_umolm2s) %>% 
  group_by(DateTime) %>% 
  summarise_all(funs(mean,sd),na.rm=TRUE)

names(flux_avg_ch4)[names(flux_avg_ch4) == 'mean'] <- 'ch4_mean_umol_m2_s'
names(flux_avg_ch4)[names(flux_avg_ch4) == 'sd'] <- 'ch4_sd_umol_m2_s'

fluxes_avg <- left_join(flux_avg_co2,flux_avg_ch4,by="DateTime")

fluxes_all <- fluxes_all %>% 
  arrange(DateTime,Rep) %>% 
  group_by(DateTime) %>% 
  summarise_all(funs(mean,sd),na.rm=TRUE) %>% 
  select(-c(Rep_mean,Rep_sd))

## Save diffusive GHG fluxes
write.csv(fluxes_all,"./Data/20220617_diffusive_fluxes_all.csv",row.names = FALSE)

write.csv(fluxes_avg,"./Data/20220617_diffusive_fluxes_avg.csv",row.names = FALSE)

###############################################################################

## Calculate stats for GHG diffusive fluxes (used in Table S4, aggregated below w/ EC data)
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

###############################################################################

fluxes_all <- fluxes_all %>% 
  drop_na()

## Plot Diffusive fluxes for Supplementary (Figure S10)
ch4_diff <- ggplot(fluxes_all)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-11-03"), col='black', size = 1, linetype = "dotted")+
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Cole_mean-ch4_flux_k600_Cole_sd,ymax=ch4_flux_k600_Cole_mean+ch4_flux_k600_Cole_sd,color="Cole"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Cole_mean,color="Cole"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Cole_mean,color="Cole"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Crusius_mean-ch4_flux_k600_Crusius_sd,ymax=ch4_flux_k600_Crusius_mean+ch4_flux_k600_Crusius_sd,color="Crusius"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Crusius_mean,color="Crusius"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Crusius_mean,color="Crusius"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Heiskanen_mean-ch4_flux_k600_Heiskanen_sd,ymax=ch4_flux_k600_Heiskanen_mean+ch4_flux_k600_Heiskanen_sd,color="Heiskanen"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Heiskanen_mean,color="Heiskanen"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Heiskanen_mean,color="Heiskanen"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_MacIntyre_mean-ch4_flux_k600_MacIntyre_sd,ymax=ch4_flux_k600_MacIntyre_mean+ch4_flux_k600_MacIntyre_sd,color="MacIntyre"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_MacIntyre_mean,color="MacIntyre"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_MacIntyre_mean,color="MacIntyre"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Read_mean-ch4_flux_k600_Read_sd,ymax=ch4_flux_k600_Read_mean+ch4_flux_k600_Read_sd,color="Read"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Read_mean,color="Read"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Read_mean,color="Read"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Soloviev_mean-ch4_flux_k600_Soloviev_sd,ymax=ch4_flux_k600_Soloviev_mean+ch4_flux_k600_Soloviev_sd,color="Soloviev"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Soloviev_mean,color="Soloviev"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Soloviev_mean,color="Soloviev"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_k600_Vachon_mean-ch4_flux_k600_Vachon_sd,ymax=ch4_flux_k600_Vachon_mean+ch4_flux_k600_Vachon_sd,color="Vachon"),size=1)+
  geom_line(mapping=aes(DateTime,ch4_flux_k600_Vachon_mean,color="Vachon"),size = 1) +
  geom_point(mapping=aes(DateTime,ch4_flux_k600_Vachon_mean,color="Vachon"),size=3) +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_viridis(discrete=TRUE)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

co2_diff <- ggplot(fluxes_all)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-11-03"), col='black', size = 1, linetype = "dotted")+
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Cole_mean-co2_flux_k600_Cole_sd,ymax=co2_flux_k600_Cole_mean+co2_flux_k600_Cole_sd,color="Cole"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Cole_mean,color="Cole"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Cole_mean,color="Cole"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Crusius_mean-co2_flux_k600_Crusius_sd,ymax=co2_flux_k600_Crusius_mean+co2_flux_k600_Crusius_sd,color="Crusius"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Crusius_mean,color="Crusius"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Crusius_mean,color="Crusius"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Heiskanen_mean-co2_flux_k600_Heiskanen_sd,ymax=co2_flux_k600_Heiskanen_mean+co2_flux_k600_Heiskanen_sd,color="Heiskanen"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Heiskanen_mean,color="Heiskanen"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Heiskanen_mean,color="Heiskanen"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_MacIntyre_mean-co2_flux_k600_MacIntyre_sd,ymax=co2_flux_k600_MacIntyre_mean+co2_flux_k600_MacIntyre_sd,color="MacIntyre"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_MacIntyre_mean,color="MacIntyre"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_MacIntyre_mean,color="MacIntyre"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Read_mean-co2_flux_k600_Read_sd,ymax=co2_flux_k600_Read_mean+co2_flux_k600_Read_sd,color="Read"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Read_mean,color="Read"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Read_mean,color="Read"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Soloviev_mean-co2_flux_k600_Soloviev_sd,ymax=co2_flux_k600_Soloviev_mean+co2_flux_k600_Soloviev_sd,color="Soloviev"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Soloviev_mean,color="Soloviev"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Soloviev_mean,color="Soloviev"),size=3) +
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_k600_Vachon_mean-co2_flux_k600_Vachon_sd,ymax=co2_flux_k600_Vachon_mean+co2_flux_k600_Vachon_sd,color="Vachon"),size=1)+
  geom_line(mapping=aes(DateTime,co2_flux_k600_Vachon_mean,color="Vachon"),size = 1) +
  geom_point(mapping=aes(DateTime,co2_flux_k600_Vachon_mean,color="Vachon"),size=3) +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_viridis(discrete=TRUE)+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

ggarrange(co2_diff,ch4_diff,
          ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Diff_fluxes_All.jpg",width = 8, height=7, units="in",dpi=320)
