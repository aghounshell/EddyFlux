# Processing Eddy Flux data from FCR
# Following initial data corrections in Eddy Pro (using standard processing)
# Specifically for data collected from 4 April 2020 (beginning of EC record) to 6 April 2021

# Original code from Brenda D'Achuna, 21 May 2021
# Modified by Alex Hounshell on 21 May 2021
# Cleaned-up for EDI publication, 6 November 2021

####################################################################

# Clear workspace
rm(list = ls())

# Download/load libraries
pacman::p_load(lubridate,readr,ggpubr,openair,REddyProc,ggplot2,dplyr)

# Set working directory
wd <- getwd()
setwd(wd)

# Read compiled file: From Eddy Pro using basic processing
# Original file from Brenda on 11 May 2021
ec <- read_csv("./EDI/20211008_EddyPro_cleaned.csv")

# Format time
ec$datetime <- as.POSIXct(paste(ec$date, ec$time), format="%m/%d/%Y %H:%M:%S", tz="EST")
ec$datetime <- as_datetime(ec$datetime)

# Set new dataframe with list of dates+times:
# every 30 minutes
# Constrain to study time period: 2020-04-05 and 2021-04-06
ts <- seq.POSIXt(as.POSIXct("2020-04-05 00:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), 
                 as.POSIXct("2021-04-05 23:30:00",'%Y-%m-%d %H:%M:%S', tz="EST"), by = "30 min")
ts2 <- data.frame(datetime = ts)

# Join Eddy Flux data with list of dates+time
ec2 <- left_join(ts2, ec, by = 'datetime')

# Make sure time stamps are okay!
ec2 %>% group_by(year = year(datetime), month = factor(month.abb[month(datetime)], levels = c("Apr", "May", "Jun",
                                                                                                 "Jul", "Aug", "Sep", "Oct", 'Nov', 
                                                                                                 'Dec', 'Jan', 'Feb', 'Mar')), 
                    hour = hour(datetime)) %>% 
  summarise(air_temperature_k = mean(air_temperature_k, na.rm = TRUE)) %>% 
  ggplot(aes(hour, air_temperature_k, col = factor(year))) + geom_point() + 
  facet_wrap(~month) + theme_bw() + ylab('Air Temp') + xlab("") +
  scale_color_brewer(palette = "Dark2")

#################################################################
# Count how many initial NAs are in CO2 and CH4 data
# Without any data processing!
#################################################################
ec2 %>% select(datetime, co2_flux_umolm2s, ch4_flux_umolm2s) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux_umolm2s))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_umolm2s))/n()*100)
# 79% data for CO2; 59% data for CH4

# Check data availability by month
ec2 %>% group_by(year(datetime), month(datetime)) %>% select(datetime, co2_flux_umolm2s, ch4_flux_umolm2s) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux_umolm2s))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_umolm2s))/n()*100)
#################################################################

# Reading in data from the Met Station for gapfilling purposes
# Include data from: EDI (2020) and from GitHub (2021 - cleaned following MET_QAQC_2020.R)
## THIS WILL NEED TO BE UPDATED!!!

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
## THIS WILL NEED TO BE UPDATED!!!

# Combine into one data frame: EDI + GitHub data
met_all <- rbind(met_edi,met_21)

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
  rename(datetime = DateTime_Adj, AirTC_Avg = AirTemp_Average_C, RH = RH_percent, Pressure = BP_Average_kPa, 
         Rain_sum = Rain_Total_mm, WS_ms_Avg = WindSpeed_Average_m_s, WindDir = WindDir_degrees,SW_in = ShortwaveRadiationUp_Average_W_m2,
         SW_out = ShortwaveRadiationDown_Average_W_m2,LW_in = InfaredRadiationUp_Average_W_m2,LW_out = InfaredRadiationDown_Average_W_m2,
         PAR_Tot_Tot = PAR_Total_mmol_m2,albedo = Albedo_Average_W_m2)

# Join with 30 minute time series
met2 <- left_join(ts2, met_30_2, by = 'datetime')

# Compare wind speeds from EC and Met
ggplot()+
  geom_point(aes(x=ec2$wind_speed,y=met2$WS_ms_Avg))+
  theme_classic(base_size = 15)

# Calculate percent of missing wind data
ec2 %>% select(datetime, wind_speed) %>% 
  summarise(wnd_na = sum(is.na(wind_speed))/n()*100)

# Use linear model to convert from Met to EC for missing time points
linearMod <- lm(ec2$wind_speed ~ met2$WS_ms_Avg)
summary(linearMod)
# For data period: EC = Met*0.533689 + 0.134411

# Check conversion
plot(ec2$wind_speed)
points(met2$WS_ms_Avg*0.533689+0.134411, col = 'red')

###########################################
# Use converted wind speed from the Met data to fill in time points with missing
# data (EC)
ec2$wind_speed <- ifelse(is.na(ec2$wind_speed),
                         met2$WS_ms_Avg*0.533689+0.134411, ec2$wind_speed)

ec2$wind_dir <- ifelse(is.na(ec2$wind_dir),
                       met2$WindDir, ec2$wind_dir)

# Visualize wind directions that are IN FRONT of the catwalk
ec2 %>% filter(wind_dir >= 250 | wind_dir <= 80) %>% 
  ggplot(aes(wind_dir, wind_speed)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 360, 45)) +
  coord_polar() + theme_bw() + xlab('Wind direction') + ylab('Wind speed')

# Filter out wind directions that are BEHIND the catwalk
# I.e., only keep data that is IN FRONT of the catwalk for both EC and Met data
ec_filt <- ec2 %>% dplyr::filter(wind_dir < 80 | wind_dir > 250)
ec_filt <- left_join(ts2, ec_filt)

met3 <- met2 %>% dplyr::filter(WindDir < 80 | WindDir > 250)
met4 <- left_join(ts2, met3)

################################################################
# count NA after filtering for wind direction BEHIND catwalk
ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100- sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)
# Now have 56% CO2 data and 42% CH4 data

# Count number of timepoints that have data
ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))
################################################################

# Remove large CO2 values
# Visualize data that is above/below abs(100)
plot(ec_filt$co2_flux)
abline(h=100)
abline(h=-100)

# Remove values that are greater than abs(100)
# NOTE: Updated from Brenda's code to use abs(100); instead of -70 to 100 filtering
# Waldo et al. 2021 used: values greater than abs(15000)
ec_filt$co2_flux <- ifelse(ec_filt$co2_flux > 100 | ec_filt$co2_flux < -100, NA, ec_filt$co2_flux)

# Remove CO2 data if QC >= 2 (aka: data that has been flagged by Eddy Pro)
ec_filt$co2_flux <- ifelse(ec_filt$qc_co2_flux >= 2, NA, ec_filt$co2_flux)

# Additionally remove CO2 data when H and LE > 2 (following CH4 filtering)
ec_filt$co2_flux <- ifelse(ec_filt$qc_co2_flux==1 & ec_filt$qc_LE>=2, NA, ec_filt$co2_flux)
ec_filt$co2_flux <- ifelse(ec_filt$qc_co2_flux==1 & ec_filt$qc_H>=2, NA, ec_filt$co2_flux)

# Remove large CH4 values
# Visualize data that is above/below abs(0.25)
plot(ec_filt$ch4_flux)
abline(h=0.25)
abline(h=-0.25)

# Remove values that are greater than abs(0.25)
# NOTE: Updated from Brenda's code to use abs(0.25)
# Waldo et al. 2021 used: values greater than abs(500)
ec_filt$ch4_flux <- ifelse(ec_filt$ch4_flux >= 0.25 | ec_filt$ch4_flux <= -0.25, NA, ec_filt$ch4_flux)

# Remove ch4 values when signal strength < 20
ec_filt$ch4_flux <- ifelse(ec_filt$rssi_77_mean < 20, NA, ec_filt$ch4_flux)

# Remove CH4 data if QC >= 2
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux >=2, NA, ec_filt$ch4_flux)

# Additionally, remove CH4 when other parameters are QA/QC'd 
# Following Waldo et al. 2021: Remove additional ch4 flux data 
# (aka: anytime ch4_qc flag = 1 & another qc_flag =2, remove)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_co2_flux>=2, NA, ec_filt$ch4_flux)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_LE>=2, NA, ec_filt$ch4_flux)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_H>=2, NA, ec_filt$ch4_flux)

# Removing qc >= 2 for H and LE
ec_filt$H <- ifelse(ec_filt$qc_H >= 2, NA, ec_filt$H)
ec_filt$LE <- ifelse(ec_filt$qc_LE >= 2, NA, ec_filt$LE)

# Remove high H values: greater than abs(200)
# NOTE: Updated to have same upper and lower magnitude bound
# Waldo et al. 2021 used abs of 200 for H
plot(ec_filt$H)
abline(h=200)
abline(h=-200)

ec_filt$H <- ifelse(ec_filt$H >= 200 | ec_filt$H <= -200, NA, ec_filt$H)

# Remove high LE values: greter than abs(300)
# NOTE: Updated to have same upper and lower magnitude bounds
# Waldo et al. 2021 used abs of 1000 for LE
plot(ec_filt$LE)
abline(h=500)
abline(h=-500)

ec_filt$LE <- ifelse(ec_filt$LE >= 500 | ec_filt$LE <= -500, NA, ec_filt$LE)

# Plotting co2 and ch4 to see if we can filter implausible values
plot(ec_filt$co2_flux, type = 'o')
plot(ec_filt$ch4_flux, type = 'o')

################################################################
# count NA after filtering for bad fluxes
ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100- sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)
# Now have 37% CO2 data and 28% CH4 data

# Count number of timepoints that have data
ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))
################################################################

# Remove CH4 when it rains
ec_filt$precip <- met2$Rain_sum
ec_filt$ch4_flux <- ifelse(ec_filt$precip > 0, NA, ec_filt$ch4_flux)

# Remove CH4 data when thermocouple was not working (apr 05 - apr 25)
ec_filt$ch4_flux <- ifelse(ec_filt$datetime >= '2021-04-05' & ec_filt$datetime <= '2021-04-25', 
                           NA, ec_filt$ch4_flux)

# Merge with timeseries
eddy_fcr <- left_join(ts2, ec_filt, by = 'datetime')

#######################################################################
# counting data again after filtering by:
# wind direction, qc, rain, unreasonable values, signal strength
eddy_fcr %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)
# 37% CO2 data; 27% CH4 data

# Number of timepoints available
eddy_fcr %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))

# Data availability by month
eddy_fcr %>% group_by(year(datetime), month(datetime)) %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

########################################################################
# Despike data using dspike.R function

# Despike NEE (CO2 flux)
source("./Scripts/despike.R")

# Calculate low, medium, and high data flags
flag <- spike_flag(eddy_fcr$co2_flux,z = 7)
NEE_low <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)
flag <- spike_flag(eddy_fcr$co2_flux,z = 5.5)
NEE_medium <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)
flag <- spike_flag(eddy_fcr$co2_flux,z = 4)
NEE_high <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)

# Plot flagged data:
plot(eddy_fcr$datetime,eddy_fcr$co2_flux,xlab = "Date", ylab = "NEE (umol m-2s-1)", col = "gray70")
points(eddy_fcr$datetime,NEE_low,col = "gray10")
points(eddy_fcr$datetime,NEE_medium,col = "blue")
points(eddy_fcr$datetime,NEE_high,col = "red")
abline(h=0)

# Combine all flagged data
eddy_fcr$NEE.low <- NEE_low
eddy_fcr$NEE.med <- NEE_medium
eddy_fcr$NEE.high <- NEE_high

#Despike CH4 flux
flag <- spike_flag(eddy_fcr$ch4_flux,z = 7)
CH4_low <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)
flag <- spike_flag(eddy_fcr$ch4_flux,z = 5.5)
CH4_medium <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)
flag <- spike_flag(eddy_fcr$ch4_flux,z = 4)
CH4_high <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)

# Plot flagged data:
plot(eddy_fcr$datetime,eddy_fcr$ch4_flux,xlab = "Date", ylab = "CH4 (umol m-2s-1)", col = "gray70")
points(eddy_fcr$datetime,CH4_low,col = "gray10")
points(eddy_fcr$datetime,CH4_medium,col = "blue")
points(eddy_fcr$datetime,CH4_high,col = "red")
abline(h=0)

# Combine all flagged data
eddy_fcr$ch4.low <- CH4_low
eddy_fcr$ch4.med <- CH4_medium
eddy_fcr$ch4.high <- CH4_high

##########################################################################
# Convert EC temperature to celsius
eddy_fcr$air_temp_celsius <- eddy_fcr$air_temperature - 273.15
eddy_fcr$sonic_temp_celsius <- eddy_fcr$sonic_temperature - 273.15

# Plot temperature
ggplot(eddy_fcr,mapping=aes(x=datetime,y=air_temp_celsius))+
  geom_point()+
  theme_classic(base_size=15)

ggplot(eddy_fcr,mapping=aes(x=datetime,y=sonic_temp_celsius))+
  geom_point()+
  theme_classic(base_size=15)

# Remove bad air temps on 10 Feb to 14 Feb
eddy_fcr$air_temp_celsius <- ifelse(eddy_fcr$datetime >= '2021-02-10' & eddy_fcr$datetime <='2021-02-14', NA, 
                                    eddy_fcr$air_temp_celsius)
eddy_fcr$sonic_temp_celsius <- ifelse(eddy_fcr$datetime >= '2021-02-10' & eddy_fcr$datetime <='2021-02-14', NA, 
                                      eddy_fcr$sonic_temp_celsius)

# Remove high temperatures on 11 November 2020
eddy_fcr$air_temp_celsius <- ifelse(eddy_fcr$datetime >= '2020-11-11' & eddy_fcr$datetime <'2020-11-12', NA, 
                                    eddy_fcr$air_temp_celsius)

eddy_fcr$sonic_temp_celsius <- ifelse(eddy_fcr$datetime >= '2020-11-11' & eddy_fcr$datetime <'2020-11-12', NA, 
                                    eddy_fcr$sonic_temp_celsius)

# Replacing Eddy Flux air temp with Met air temp
# Check correlations between EC temp and Met Temp
ggplot()+
  geom_point(aes(x=eddy_fcr$air_temp_celsius,y=met2$AirTC_Avg))+
  theme_classic(base_size = 15)

# Calculate percent of missing air temp data
eddy_fcr %>% select(datetime, air_temp_celsius) %>% 
  summarise(temp_na = sum(is.na(air_temp_celsius))/n()*100)

# Use linear model to estimate EC temp from Met temp
linearMod <- lm(eddy_fcr$air_temp_celsius ~ met2$AirTC_Avg)
summary(linearMod)
# EC = Met*0.93201-0.58405

eddy_fcr$air_temp_celsius <- ifelse(is.na(eddy_fcr$air_temp_celsius),
                                    met2$AirTC_Avg*0.93201-0.58405, eddy_fcr$air_temp_celsius)

# Check Air Temp
ggplot(eddy_fcr,mapping=aes(x=datetime,y=air_temp_celsius))+
  geom_point()+
  theme_classic(base_size=15)

# Check correlations between EC temp (sonic) and Met Temp
ggplot()+
  geom_point(aes(x=eddy_fcr$sonic_temp_celsius,y=met2$AirTC_Avg))+
  theme_classic(base_size = 15)

# Calculate percent of missing sonic temp data
eddy_fcr %>% select(datetime, sonic_temp_celsius) %>% 
  summarise(temp_na = sum(is.na(sonic_temp_celsius))/n()*100)

# Use linear model to estimate EC temp and Met temp
linearMod <- lm(eddy_fcr$sonic_temp_celsius ~ met2$AirTC_Avg)
summary(linearMod)
# EC = Met*1.017914 - 0.260521

eddy_fcr$sonic_temp_celsius <- ifelse(is.na(eddy_fcr$sonic_temp_celsius),
                                      eddy_fcr$air_temp_celsius*1.017914-0.260521, eddy_fcr$sonic_temp_celsius)

# Check sonic temp
ggplot(eddy_fcr,mapping=aes(x=datetime,y=sonic_temp_celsius))+
  geom_point()+
  theme_classic(base_size=15)

# Replacing Eddy Flux RH with Met RH
# Check correlations between EC and Met RH
ggplot()+
  geom_point(aes(x=eddy_fcr$RH,y=met2$RH))+
  theme_classic(base_size = 15)

# Calculate percent of missing sonic temp data
eddy_fcr %>% select(datetime, RH) %>% 
  summarise(RH_na = sum(is.na(RH))/n()*100)

# Use linear model to estimate EC RH from Met RH
linearMod <- lm(eddy_fcr$RH ~ met2$RH)
summary(linearMod)
# EC = Met*0.812051 + 6.626110

eddy_fcr$RH <- ifelse(is.na(eddy_fcr$RH),
                      met2$RH*0.812051+6.626110, eddy_fcr$RH)

# Add Met data to gapfill fluxes
eddy_fcr$SW_in <- met2$SW_in
eddy_fcr$SW_out <- met2$SW_out
eddy_fcr$par_tot <- met2$PAR_Tot_Tot
eddy_fcr$air_pressure <- met2$Pressure
eddy_fcr$LW_in <- met2$LW_in
eddy_fcr$LW_out <- met2$LW_out
eddy_fcr$albedo <- met2$albedo

# Calculate vapor pressure deficit from RH and Air Temp using Met data
eddy_fcr$VPD <- ifelse(is.na(eddy_fcr$VPD), 
                       fCalcVPDfromRHandTair(rH = eddy_fcr$RH, Tair = eddy_fcr$air_temp_celsius)*100, 
                       eddy_fcr$VPD)

# Remove funky PAR_Tot values
eddy_fcr$par_tot <- ifelse(eddy_fcr$datetime >='2020-07-03' & eddy_fcr$datetime <= '2020-07-22', NA, eddy_fcr$par_tot)

# Replace wind_direction in EC data with Met data if missing
eddy_fcr$wind_dir <- ifelse(is.na(eddy_fcr$wind_dir), met4$WindDir, eddy_fcr$wind_dir)

# QA/QC'ing LW_out data
eddy_fcr$LW_out <- ifelse(eddy_fcr$LW_out <= 360, NA, eddy_fcr$LW_out)
eddy_fcr$LW_out <- ifelse(eddy_fcr$datetime >= '2020-06-22' & eddy_fcr$datetime <= '2020-07-13' & eddy_fcr$LW_out <= 420, NA, eddy_fcr$LW_out)

# Calculating Net Radiation
eddy_fcr$Rn <- eddy_fcr$SW_in - eddy_fcr$SW_out + eddy_fcr$LW_in - eddy_fcr$LW_out

plot(eddy_fcr$Rn, type = 'o')

plot(eddy_fcr$SW_in)

plot(eddy_fcr$VPD/1000)  # in kpa

###############################################################################
# Filter out all the values (x_peak) that are out of the reservoir
eddy_fcr$footprint_flag <- ifelse(eddy_fcr$wind_dir >= 15 & eddy_fcr$wind_dir <= 90 & eddy_fcr$x_peak >= 40, 1, 
                                  ifelse(eddy_fcr$wind_dir < 15 & eddy_fcr$wind_dir > 327 & eddy_fcr$x_peak > 120, 1,
                                         ifelse(eddy_fcr$wind_dir < 302 & eddy_fcr$wind_dir >= 250 & eddy_fcr$x_peak > 50, 1, 0)))

# Remove flagged data
eddy_fcr_footprint <- eddy_fcr %>% filter(footprint_flag == 0)

# Visualize wind directions that were kept
eddy_fcr_footprint %>% ggplot(aes(wind_dir, x_peak)) + 
  geom_hline(yintercept = 40, col = 'goldenrod2', lwd = 2) +
  geom_hline(yintercept = 50, col = 'green', lwd = 1.4) +
  geom_hline(yintercept = 100, col = 'blue', lwd = 1.4) +
  geom_hline(yintercept = 120, col = 'gray2', lwd = 1.4) +
  geom_hline(yintercept = 150, col = 'red',lwd = 1.4) +
  geom_point() +
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 360, 45)) +
  theme_bw() + 
  coord_polar()

# merge with timeseries
eddy_fcr_footprint_full <- left_join(ts2, eddy_fcr_footprint)

######################################################################
# FILTERING BY USTAR AND GAPFILLING
######################################################################

# Setting up a new process on REddyProc
eddy_fcr3 <- eddy_fcr_footprint_full %>% 
  select(DateTime = datetime, daytime, NEE = NEE.med, ch4_flux = ch4.med, VPD, 
         H, LE, Tair = sonic_temp_celsius, rH = RH, Ustar = `u*`, u = wind_speed, 
         pressure = air_pressure, L, z_d_L = `(z-d)/L`, sigma_v = v_var, 
         precip, Rn, SW_in, SW_out, LW_out, LW_in, albedo, par_tot, wind_dir, 
         airP = air_pressure) %>% 
  mutate(VPD = VPD/100,
         z_d = z_d_L*L,
         ln_z_d = log(z_d)) %>% 
  rename(Rg = SW_in,
         PAR = par_tot)

########################################################################
# Count available data before gapfilling
eddy_fcr3 %>% select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(NEE))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)
# 28% CO2 data; 21% CH4 data

# Total number of available time points
eddy_fcr3 %>% select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = n()-sum(is.na(NEE)),
            ch4_available = n()-sum(is.na(ch4_flux)))

# Counting data availability by month
eddy_fcr3 %>% group_by(year(DateTime), month(DateTime)) %>% 
  select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(NEE))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

# Visualizing data by wind direction/wind speed
windRose(mydata = eddy_fcr3, ws = "u", wd = "wind_dir", 
         width = 3, key.position = 'bottom', 
         offset = 3, paddle = FALSE, key.header = 'Wind speed (m/s)', 
         key.footer = ' ', dig.lab = 2, annotate = FALSE,
         angle.scale = 45, ws.int = 1, breaks = c(0, 2, 4, 6, 8))

###########################################################################
# get ustar distribution and filter by ustar
###########################################################################
# Gapfilling - using periods of similar met.
# Following Eddy Proc
Eproc <- sEddyProc$new('FCR', eddy_fcr3, c('NEE','Tair', 'VPD',
                                           'rH','H', 'LE', 'Ustar', 
                                           'ch4_flux', 'u', 'PAR', 
                                           'SW_out', 'Rg', 
                                           'Rn', 'LW_out', 'LW_in'))

# Look at available data for each year/flux
Eproc$sPlotFingerprintY('NEE',Year = '2020')
Eproc$sPlotFingerprintY('NEE',Year = '2021')

Eproc$sPlotFingerprintY('ch4_flux',Year = '2020')
Eproc$sPlotFingerprintY('ch4_flux',Year = '2021')

# gapfill air temperature, solar radiation, par, H and LE: include uncertainty calculations
Eproc$sMDSGapFill('Tair', V1 = 'Rg', V2 = 'VPD', FillAll = TRUE)
Eproc$sMDSGapFill('Rg',  V2 = 'VPD', V3 = 'Tair', FillAll = TRUE)
Eproc$sMDSGapFill('PAR', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair', FillAll = TRUE)
Eproc$sMDSGapFill('Rn', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair', FillAll = TRUE)
Eproc$sMDSGapFill('H', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair', FillAll = TRUE)
Eproc$sMDSGapFill('LE', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair', FillAll = TRUE)

# Estimate ustar threshold distribution by bootstrapping the data
Eproc$sEstimateUstarScenarios(UstarColName = 'Ustar', NEEColName = 'NEE', RgColName= 'Rg',
                              nSample = 200L, probs = c(0.05, 0.7, 0.95))

# Beef up uncertainty following: https://github.com/bgctw/REddyProc/blob/master/vignettes/useCase.md
Eproc$sEstimateUstarScenarios(UstarColName = 'Ustar', NEEColName = 'NEE', RgColName= 'Rg', 
  nSample = 1000, probs = seq(0.025,0.975,length.out = 39))

Eproc$sGetEstimatedUstarThresholdDistribution()

Eproc$sGetUstarScenarios()

Eproc$sMDSGapFillUStarScens(fluxVar = 'NEE', FillAll = TRUE)
Eproc$sMDSGapFillUStarScens(fluxVar = 'ch4_flux', FillAll = TRUE)

# Try plotting daily sums
Eproc$sPlotDailySums('NEE_uStar_f','NEE_uStar_fsd')
Eproc$sPlotDailySums('ch4_flux_uStar_f')

Eproc$sPlotDiurnalCycle('NEE_uStar_f')
Eproc$sPlotDiurnalCycle('ch4_flux_uStar_f')

Eproc$sPlotHHFluxes('NEE_uStar_f')
Eproc$sPlotHHFluxes('ch4_flux_uStar_f')

# Export results from Eddy Proc
filled_fcr <- Eproc$sExportResults()
fcr_gf <- cbind(eddy_fcr3, filled_fcr)

### Quick plots as a gut-check/initial cut

fcr_gf %>% ggplot() + 
  geom_point(aes(DateTime, ch4_flux_uStar_orig), alpha = 0.3) +
  geom_line(aes(DateTime, ch4_flux_uStar_f - ch4_flux_uStar_fsd), alpha = 0.3)+
  geom_line(aes(DateTime, ch4_flux_uStar_fsd + ch4_flux_uStar_fsd), alpha = 0.3)+
  geom_line(aes(DateTime, ch4_flux_uStar_f), col = 'red', alpha = 0.5) +
  theme_bw() +
  xlab("") + ylab(expression(~CH[4]~flux~(mu~mol~m^-2~s^-1)))

# Plot as points + daily average
fcr_gf_mean <- fcr_gf %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d")) %>% 
  group_by(DateTime) %>% 
  summarise_all(mean,na.rm=TRUE)

fcr_gf %>% ggplot() +
  #geom_point(aes(DateTime, NEE),alpha=0.1) +
  geom_line(aes(DateTime, NEE_uStar_f), col='red', alpha = 0.4) +
  geom_line(aes(DateTime, NEE_U50_f), col = 'blue', alpha = 0.4)+
  theme_bw() +
  xlab("") + ylab(expression(~CO[2]~flux~(mu~mol~m^-2~s^-1)))

ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_point(fcr_gf, mapping = aes(DateTime, NEE_uStar_orig*60*60*24*44.01/1e6),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(fcr_gf_mean, mapping = aes(DateTime, NEE_uStar_f*60*60*24*44.01/1e6),color="#E63946",size = 1)+
  theme_classic(base_size = 15)+
  xlab("") + ylab(expression(~CO[2]~flux~(g~m^-2~d^-1)))

ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_point(fcr_gf, mapping = aes(DateTime, ch4_flux_uStar_orig*60*60*24*16.04/1e6),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(fcr_gf_mean, mapping = aes(DateTime, ch4_flux_uStar_f*60*60*24*16.04/1e6),color="#E63946",size = 1)+
  theme_classic(base_size = 15)+
  xlab("") + ylab(expression(~CH[4]~flux~(g~m^-2~d^-1)))

# Save the exported data
write_csv(fcr_gf, "./Data/20210813_EC_processed.csv")
