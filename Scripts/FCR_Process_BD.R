# working with FCR data - processing
# Original code from Brenda D'Achuna, 21 May 2021
# Modifed by AGH on 21 May 2021

rm(list = ls())

pacman::p_load(lubridate,readr,ggpubr,openair,REddyProc,ggplot2,dplyr)

# Set working directory
wd <- getwd()
setwd(wd)

# read compiled file: From Eddy Pro using basic processing
# Original file from Brenda on 11 May 2021

ec <- read_csv("./Data/FCR_2021-05-06_upto.csv")

ec$datetime <- as.POSIXct(paste(ec$date, ec$time), format="%m/%d/%Y %H:%M:%S", tz="EST")
ec$datetime <- as_datetime(ec$datetime)

# convert -9999 to NA
ec[ec == -9999] <- NA


# setting up a new dataframe with complete dates

# change end date

tail(ec$datetime)

ts <- seq.POSIXt(as.POSIXct("2020-04-04 11:30:00",'%Y-%m-%d %H:%M:%S', tz="EST"), 
                 as.POSIXct("2021-05-06 13:30:00",'%Y-%m-%d %H:%M:%S', tz="EST"), by = "30 min")
ts2 <- data.frame(datetime = ts)


# join dataframes to gapfill missing dates on ec data
ec2 <- left_join(ts2, ec, by = 'datetime')

#################################################################

# count how many initial NAs are in CO2 and CH4 data
#################################################################

ec2 %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

ec2 %>% group_by(year(datetime), month(datetime)) %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)


#################################################################

# reading data from catwalk and from meteorological station at FCR

# From Brenda: updated to include EDI met data + cleaned Met data from GitHub (following MET_QAQC_2020.R)
#met <- read_csv('./Data/met_data.csv')
#met <- met %>% dplyr::rename(datetime = date) %>% 
#  mutate(datetime = as.POSIXct(datetime, format="%m/%d/%Y %H:%M:%S", tz="EST"))

# Loading in Met data from EDI and Github (cleaned following MET_QAQC_2020.R)
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

met_all <- met_all %>% 
  filter(DateTime >= as.POSIXct("2019-12-31 00:15:00"))

# Select data every 30 minutes from Jan 2020 to end of met data
met_all$Breaks <- cut(met_all$DateTime,breaks = "30 mins",right=FALSE)
met_all$Breaks <- ymd_hms(as.character(met_all$Breaks))

met_30 <- met_all %>% 
  select(DateTime,BP_Average_kPa,AirTemp_Average_C,RH_percent,ShortwaveRadiationUp_Average_W_m2,ShortwaveRadiationDown_Average_W_m2,
         InfaredRadiationUp_Average_W_m2,InfaredRadiationDown_Average_W_m2,Albedo_Average_W_m2,WindSpeed_Average_m_s,WindDir_degrees,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(mean,na.rm=TRUE)

met_30_rain <- met_all %>% 
  select(Rain_Total_mm,PAR_Total_mmol_m2,Breaks) %>% 
  group_by(Breaks) %>% 
  summarise_all(sum,na.rm=TRUE)

met_30_2 <- cbind.data.frame(met_30,met_30_rain)

met_30_2 <- met_30_2 %>% 
  select(-Breaks) %>% 
  mutate(DateTime_Adj = DateTime + 30) %>% 
  select(-DateTime) %>% 
  rename(datetime = DateTime_Adj, AirTC_Avg = AirTemp_Average_C, RH = RH_percent, Pressure = BP_Average_kPa, 
         Rain_sum = Rain_Total_mm, WS_ms_Avg = WindSpeed_Average_m_s, WindDir = WindDir_degrees,SW_in = ShortwaveRadiationUp_Average_W_m2,
         SW_out = ShortwaveRadiationDown_Average_W_m2,LW_in = InfaredRadiationUp_Average_W_m2,LW_out = InfaredRadiationDown_Average_W_m2,
         PAR_Tot_Tot = PAR_Total_mmol_m2,albedo = Albedo_Average_W_m2)

# just in case gapfilling dates on met data

met2 <- left_join(ts2, met_30_2, by = 'datetime')

# reading catwalk data - AGH: I don't think we need catwalk data? Removing for now!

#cat <- read_csv('./Data/catwalk_data.csv')
#cat <- cat %>% dplyr::rename(datetime = date) %>% 
#  mutate(datetime = as.POSIXct(datetime, format="%m/%d/%Y %H:%M:%S", tz="EST"))

## gapfilling dates on met data

#cat2 <- left_join(ts2, cat, by = 'datetime')

# compare wind speeds from met and ec

ggplot()+
  geom_point(aes(x=ec2$wind_speed,y=met2$WS_ms_Avg))+
  theme_classic(base_size = 15)

linearMod <- lm(ec2$wind_speed ~ met2$WS_ms_Avg)

plot(ec2$wind_speed)
points(met2$WS_ms_Avg*0.5238+0.1489, col = 'red')

###########################################

# if there is no wind speed or wind dir data, then enter value from met

ec2$wind_speed <- ifelse(is.na(ec2$wind_speed),
                         met2$WS_ms_Avg*0.5238+0.1489, ec2$wind_speed)

ec2$wind_dir <- ifelse(is.na(ec2$wind_dir),
                       met2$WindDir, ec2$wind_dir)


ec2 %>% filter(wind_dir >= 250 | wind_dir <= 80) %>% 
  ggplot(aes(wind_dir, wind_speed)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 360, 45)) +
  coord_polar() + theme_bw() + xlab('Wind direction') + ylab('Wind speed')


# filtering by wind direction: filtering out time points when the wind is coming from BEHIND the catwalk
# Serves as an initial exculsion - additional wind filtering below
ec_filt <- ec2 %>% dplyr::filter(wind_dir < 80 | wind_dir > 250)

met3 <- met2 %>% dplyr::filter(WindDir < 80 | WindDir > 250)

met4 <- left_join(ts2, met3)

ec_filt <- left_join(ts2, ec_filt)


################################################################
# count NA after filtering for wind direction

ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100- sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

ec_filt %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))

################################################################

# removing large values from co2

plot(ec_filt$co2_flux)
abline(h=100)
abline(h=-100)

# Absolute value of 100: AGH updated to -100 to 100
ec_filt$co2_flux <- ifelse(ec_filt$co2_flux > 100 | ec_filt$co2_flux < -100, NA, ec_filt$co2_flux)

# removing qc = 2 for co2

ec_filt$co2_flux <- ifelse(ec_filt$qc_co2_flux >= 2, NA, ec_filt$co2_flux)

# removing large values from ch4

plot(ec_filt$ch4_flux)
abline(h=0.25)
abline(h=-0.25)

# Constrain to -0.25 to 0.25
ec_filt$ch4_flux <- ifelse(ec_filt$ch4_flux >= 0.25 | ec_filt$ch4_flux <= -0.25, NA, ec_filt$ch4_flux)

# removing ch4 values when signal strength < 20

ec_filt$ch4_flux <- ifelse(ec_filt$rssi_77_mean < 20, NA, ec_filt$ch4_flux)

# removing qc = 2 for ch4

ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux >=2, NA, ec_filt$ch4_flux)

# removing CH4 data when it rains

ec_filt$precip <- met2$Rain_sum

ec_filt$ch4_flux <- ifelse(ec_filt$precip > 0, NA, ec_filt$ch4_flux)

# remove CH4 data when thermocouple was not working (apr 05 - apr 25)

ec_filt$ch4_flux <- ifelse(ec_filt$datetime >= '2021-04-05' & ec_filt$datetime <= '2021-04-25', 
                           NA, ec_filt$ch4_flux)

# Following Waldo et al. 2021: Remove additional ch4 flux data (aka: anytime ch4_qc flag = 1 & another qc_flag =2, remove)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_co2_flux>=2, NA, ec_filt$ch4_flux)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_LE>=2, NA, ec_filt$ch4_flux)
ec_filt$ch4_flux <- ifelse(ec_filt$qc_ch4_flux==1 & ec_filt$qc_H>=2, NA, ec_filt$ch4_flux)

# removing qc = 2 for H and LE

ec_filt$H <- ifelse(ec_filt$qc_H >= 2, NA, ec_filt$H)
ec_filt$LE <- ifelse(ec_filt$qc_LE >= 2, NA, ec_filt$LE)

# Waldo et al. 2021 used abs of 200 for H and 1000 LE
plot(ec_filt$H)
abline(h=200)
abline(h=-200)

ec_filt$H <- ifelse(ec_filt$H >= 200 | ec_filt$H <= -200, NA, ec_filt$H)

plot(ec_filt$LE)
abline(h=300)
abline(h=-300)

ec_filt$LE <- ifelse(ec_filt$LE >= 300 | ec_filt$LE <= -300, NA, ec_filt$LE)


# Plotting co2 and ch4 to see if we can filter implausible values

plot(ec_filt$co2_flux, type = 'o')
plot(ec_filt$ch4_flux, type = 'o')


# adding missing dates
head(ec_filt$datetime)
tail(ec_filt$datetime)


eddy_fcr <- left_join(ts2, ec_filt, by = 'datetime')

#######################################################################
# counting data again after filtering by:
# wind speed, qc, rain, unreasonable values, signal strength


eddy_fcr %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

eddy_fcr %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))


eddy_fcr %>% group_by(year(datetime), month(datetime)) %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

########################################################################

# despike co2

# Despike NEE
source("./Scripts/despike.R")

flag <- spike_flag(eddy_fcr$co2_flux,z = 7)
NEE_low <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)
flag <- spike_flag(eddy_fcr$co2_flux,z = 5.5)
NEE_medium <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)
flag <- spike_flag(eddy_fcr$co2_flux,z = 4)
NEE_high <- ifelse(flag == 1, NA, eddy_fcr$co2_flux)

plot(eddy_fcr$datetime,eddy_fcr$co2_flux,xlab = "Date", ylab = "NEE (umol m-2s-1)", col = "gray70")
points(eddy_fcr$datetime,NEE_low,col = "gray10")
points(eddy_fcr$datetime,NEE_medium,col = "blue")
points(eddy_fcr$datetime,NEE_high,col = "red")
abline(h=0)

#use medium despiking to remove outliers
eddy_fcr$NEE.med <- NEE_medium


#Despike CH4 flux
flag <- spike_flag(eddy_fcr$ch4_flux,z = 7)
CH4_low <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)
flag <- spike_flag(eddy_fcr$ch4_flux,z = 5.5)
CH4_medium <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)
flag <- spike_flag(eddy_fcr$ch4_flux,z = 4)
CH4_high <- ifelse(flag == 1, NA, eddy_fcr$ch4_flux)


plot(eddy_fcr$datetime,eddy_fcr$ch4_flux,xlab = "Date", ylab = "CH4 (umol m-2s-1)", col = "gray70")
points(eddy_fcr$datetime,CH4_low,col = "gray10")
points(eddy_fcr$datetime,CH4_medium,col = "blue")
points(eddy_fcr$datetime,CH4_high,col = "red")
abline(h=0)


eddy_fcr$ch4.low <- CH4_low
eddy_fcr$ch4.med <- CH4_medium
eddy_fcr$ch4.hig <- CH4_high


##########################################################################


head(eddy_fcr$datetime)
tail(eddy_fcr$datetime)

eddy_fcr$air_temp_celsius <- eddy_fcr$air_temperature - 273.15
eddy_fcr$sonic_temp_celsius <- eddy_fcr$sonic_temperature - 273.15

eddy_fcr$air_temp_celsius <- ifelse(eddy_fcr$datetime >= '2021-02-10' & eddy_fcr$datetime <='2021-02-14' & 
                                      eddy_fcr$air_temp_celsius >= 15, NA, 
                                    eddy_fcr$air_temp_celsius)


eddy_fcr$sonic_temp_celsius <- ifelse(eddy_fcr$datetime >= '2021-02-10' & eddy_fcr$datetime <='2021-02-14' & 
                                        eddy_fcr$sonic_temp_celsius >= 15, NA, 
                                      eddy_fcr$sonic_temp_celsius)


# Replacing Eddy Flux air temp with Met air temp
ggplot()+
  geom_point(aes(x=eddy_fcr$air_temp_celsius,y=met2$AirTC_Avg))+
  theme_classic(base_size = 15)

linearMod <- lm(eddy_fcr$air_temp_celsius ~ met2$AirTC_Avg)

eddy_fcr$air_temp_celsius <- ifelse(is.na(eddy_fcr$air_temp_celsius),
                                    met2$AirTC_Avg*0.9724-0.4088, eddy_fcr$air_temp_celsius)

ggplot()+
  geom_point(aes(x=eddy_fcr$sonic_temp_celsius,y=met2$AirTC_Avg))+
  theme_classic(base_size = 15)

linearMod <- lm(eddy_fcr$sonic_temp_celsius ~ met2$AirTC_Avg)


eddy_fcr$sonic_temp_celsius <- ifelse(is.na(eddy_fcr$sonic_temp_celsius),
                                      eddy_fcr$air_temp_celsius*1.0082-0.1038, eddy_fcr$sonic_temp_celsius)

# Replacing Eddy Flux RH with Met RH
ggplot()+
  geom_point(aes(x=eddy_fcr$RH,y=met2$RH))+
  theme_classic(base_size = 15)

linearMod <- lm(eddy_fcr$RH ~ met2$RH)

eddy_fcr$RH <- ifelse(is.na(eddy_fcr$RH),
                      met2$RH*0.8116+6.8434, eddy_fcr$RH)

# adding meteorological variables to gapfill data

eddy_fcr$SW_in <- met2$SW_in
eddy_fcr$SW_out <- met2$SW_out
eddy_fcr$par_tot <- met2$PAR_Tot_Tot
eddy_fcr$air_pressure <- met2$Pressure
eddy_fcr$LW_in <- met2$LW_in
eddy_fcr$LW_out <- met2$LW_out
eddy_fcr$albedo <- met2$albedo
eddy_fcr$air_pressure <- ifelse(is.na(eddy_fcr$air_pressure), 
                                met2$Pressure, eddy_fcr$air_pressure)

eddy_fcr$VPD <- ifelse(is.na(eddy_fcr$VPD), 
                       fCalcVPDfromRHandTair(rH = eddy_fcr$RH, Tair = eddy_fcr$air_temp_celsius)*100, 
                       eddy_fcr$VPD)

eddy_fcr$par_tot <- ifelse(eddy_fcr$datetime >='2020-07-03' & eddy_fcr$datetime <= '2020-07-22', NA, eddy_fcr$par_tot)

eddy_fcr$wind_dir <- ifelse(is.na(eddy_fcr$wind_dir), met4$WindDir, eddy_fcr$wind_dir)

eddy_fcr$LW_out <- ifelse(eddy_fcr$LW_out <= 360, NA, eddy_fcr$LW_out)


eddy_fcr$LW_out <- ifelse(eddy_fcr$datetime >= '2020-06-22' & eddy_fcr$datetime <= '2020-07-13' & eddy_fcr$LW_out <= 420, NA, eddy_fcr$LW_out)

# get net radiation
eddy_fcr$Rn <- eddy_fcr$SW_in - eddy_fcr$SW_out + eddy_fcr$LW_in - eddy_fcr$LW_out

plot(eddy_fcr$Rn, type = 'o')

plot(eddy_fcr$SW_in)

plot(eddy_fcr$VPD/1000)  # in kpa


###############################################################################
# filter out all the values (x_peak) that are out of the reservoir
# Filtering for length and direction?

# Use percentiles: 70th or 80th instead of median (50/50) - check x-peak; distance (how sensitive is this to 70?)

eddy_fcr$footprint_flag <- ifelse(eddy_fcr$wind_dir >= 15 & eddy_fcr$wind_dir <= 90 & eddy_fcr$x_peak >= 40, 1, 
                                  ifelse(eddy_fcr$wind_dir < 15 & eddy_fcr$wind_dir > 327 & eddy_fcr$x_peak > 120, 1,
                                         ifelse(eddy_fcr$wind_dir < 302 & eddy_fcr$wind_dir >= 250 & eddy_fcr$x_peak > 50, 1, 0)))



eddy_fcr_footprint <- eddy_fcr %>% filter(footprint_flag == 0)

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


# merge with full df

eddy_fcr_footprint_full <- left_join(ts2, eddy_fcr_footprint)

# counting data

eddy_fcr_footprint_full %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

eddy_fcr_footprint_full %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = n() - sum(is.na(co2_flux)),
            ch4_available = n() -sum(is.na(ch4_flux)))


eddy_fcr_footprint_full %>% group_by(year(datetime), month(datetime)) %>% select(datetime, co2_flux, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(co2_flux))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)


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
# count available data before gapfilling


eddy_fcr3 %>% select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(NEE))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

eddy_fcr3 %>% select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = n()-sum(is.na(NEE)),
            ch4_available = n()-sum(is.na(ch4_flux)))

eddy_fcr3 %>% group_by(year(DateTime), month(DateTime)) %>% 
  select(DateTime, NEE, ch4_flux) %>% 
  summarise(co2_available = 100-sum(is.na(NEE))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux))/n()*100)

windRose(mydata = eddy_fcr3, ws = "u", wd = "wind_dir", 
         width = 3, key.position = 'bottom', 
         offset = 3, paddle = FALSE, key.header = 'Wind speed (m/s)', 
         key.footer = ' ', dig.lab = 2, annotate = FALSE,
         angle.scale = 45, ws.int = 1, breaks = c(0, 2, 4, 6, 8))


###########################################################################
# get ustar distribution and filter by ustar
###########################################################################
# Gapfilling - using periods of similar met.
# Eddy proc - Max Planck for biogeochem

Eproc <- sEddyProc$new('FCR', eddy_fcr3, c('NEE','Tair', 'VPD',
                                           'rH','H', 'LE', 'Ustar', 
                                           'ch4_flux', 'u', 'PAR', 
                                           'SW_out', 'Rg', 
                                           'Rn', 'LW_out', 'LW_in'))

# gapfill air temperature, solar radiation, par, H and LE
Eproc$sMDSGapFill('Tair', V1 = 'Rg', V2 = 'VPD')
Eproc$sMDSGapFill('Rg',  V2 = 'VPD', V3 = 'Tair')
Eproc$sMDSGapFill('PAR', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair')
Eproc$sMDSGapFill('Rn', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair')
Eproc$sMDSGapFill('H', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair')
Eproc$sMDSGapFill('LE', V1 = 'Rg', V2 = 'VPD', V3 = 'Tair')

# estimate ustar threshold distribution by bootstrapping the data

Eproc$sEstimateUstarScenarios(UstarColName = 'Ustar', NEEColName = 'NEE', RgColName= 'Rg',
                              nSample = 200L, probs = c(0.05, 0.5, 0.95))

Eproc$sGetUstarScenarios()

Eproc$sMDSGapFillUStarScens(fluxVar = 'NEE')
Eproc$sMDSGapFillUStarScens(fluxVar = 'ch4_flux')


# exporting results

filled_fcr <- Eproc$sExportResults()
fcr_gf <- cbind(eddy_fcr3, filled_fcr)


fcr_gf %>% ggplot() + 
  geom_line(aes(DateTime, ch4_flux_uStar_orig)) +
  geom_line(aes(DateTime, ch4_flux_uStar_f), col = 'red', alpha = 0.3) +
  theme_bw() +
  xlab("") + ylab(expression(~CH[4]~flux~(mu~mol~m^-2~s^-1)))

fcr_gf %>% ggplot() +
  geom_line(aes(DateTime, NEE)) +
  geom_line(aes(DateTime, NEE_uStar_f), col='red', alpha = 0.3) +
  theme_bw() +
  xlab("") + ylab(expression(~CO[2]~flux~(mu~mol~m^-2~s^-1)))


# saving the data 

write_csv(fcr_gf, "./Data/20210607_EC_processed.csv")
