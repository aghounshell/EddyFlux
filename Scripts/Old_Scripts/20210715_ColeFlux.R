### Calculating diffusive fluxes from dissolved GHGs using the Cole method
### Following McClure et al. 2018
### 21 May 2021, A Hounshell

# Following code from 2019: for calculating fluxes
# https://github.com/aghounshell/GHG/blob/master/Scripts/GHG_Flux.R
# https://github.com/aghounshell/GHG/blob/master/Scripts/WindComps.R
# https://github.com/aghounshell/GHG/blob/master/Scripts/Atm_CH4.R

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples)

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
  rename(DateTime = DateTime_Adj)

# Gut-check plot: wind
ggplot(met_all,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s))+
  geom_line()+
  ylim(0,11)+
  theme_classic(base_size = 17)

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
wind <- met_30_2 %>% 
  select(WindSpeed_Average_m_s) %>% 
  rename(wnd = WindSpeed_Average_m_s)

height <- 3
u10 <- wind.scale.base(wind,height)

# First calculate k600 using LakeMetabolizer
# in m/d
k600_2 <- k.cole.base(u10)

k600 <- cbind(met_30_2,k600_2) %>% 
  select(c(DateTime,BP_Average_kPa:Albedo_Average_W_m2,wnd)) %>% 
  rename(k600_md = wnd)

# Merge and extrapolate GHG data
fluxes <- left_join(k600,ghg_1,by=c("DateTime")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep1 = ch4_umolL, co2_umolL_rep1 = co2_umolL)

fluxes_2 <- left_join(fluxes,ghg_2,by=c("DateTime")) %>% 
  select(-c(Depth_m,Rep)) %>% 
  mutate(ch4_umolL = na.fill(na.approx(ch4_umolL,na.rm=FALSE),"extend")) %>% 
  mutate(co2_umolL = na.fill(na.approx(co2_umolL,na.rm=FALSE),"extend")) %>% 
  rename(ch4_umolL_rep2 = ch4_umolL, co2_umolL_rep2 = co2_umolL)

# Why don't we use EC concentrations (CO2 and CH4 for atmospheric concentrations?)
# Load in EC data from Eddy Pro
ec <- read_csv("./Data/FCR_2021-05-06_upto.csv") %>% 
  mutate(DateTime = as.POSIXct(paste(date, time), format="%m/%d/%Y %H:%M:%S", tz="EST"))

# convert -9999 to NA
ec[ec == -9999] <- NA

# Plot EC concentrations (just to see!)
ggplot(ec,mapping=aes(DateTime, co2_mole_fraction))+
  geom_line()

ggplot(ec,mapping=aes(DateTime, ch4_mole_fraction))+
  geom_line()

# Do some light QA/QC'ing
ec_2 <- ec %>% 
  mutate(co2_mole_fraction = ifelse(co2_mole_fraction < 430 & co2_mole_fraction > 330, co2_mole_fraction, NA)) %>% 
  mutate(ch4_mole_fraction = ifelse(ch4_mole_fraction < 2.3, ch4_mole_fraction, NA)) %>% 
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

# Calculate fluxes: in umol/m2/s
fluxes_2 <- fluxes_2 %>% 
  mutate(nv = (BP_Average_kPa*0.00986923/(0.0820573660809596*(AirTemp_Average_C + 273.15)))) %>% # units = mols/L
  mutate(ch4_umolL_atm = ch4_mole_fraction/1e6*nv*1e6) %>% # units = umol/L
  mutate(ch4_flux_1 = 1000*k600_md*(ch4_umolL_rep1-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(ch4_flux_2 = 1000*k600_md*(ch4_umolL_rep2-ch4_umolL_atm)/24/60/60) %>%  # units = umol C/m2/s
  mutate(co2_umolL_atm = co2_mole_fraction/1e6*nv*1e6) %>%  # units = umol/L
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
  summarise_all(funs(mean,sd),na.rm=TRUE) %>% 
  select(-c(Rep_mean,Rep_sd))

# Plot to check?
ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,ch4_flux_umolm2s_mean/1000*60*60*24))+
  #geom_ribbon(fluxes_all,mapping=aes(DateTime,ymin = ch4_umolm2s-ch4_umolm2s_sd,ymax = ch4_umolm2s+ch4_umolm2s_sd),fill="grey")+
  theme_classic(base_size =15)

ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,co2_flux_umolm2s_mean/1000*60*60*24))+
  theme_classic(base_size=15)

# Select dates where we actually have GHG concentrations
ghg_fluxes <- left_join(ghg_1,fluxes_all,by="DateTime")

# Calculate stats for GHG diffusive fluxes
ghg_stats <- ghg_fluxes %>% 
  summarise(min_co2 = min(co2_flux_umolm2s_mean,na.rm = TRUE),
            max_co2 = max(co2_flux_umolm2s_mean,na.rm=TRUE),
            med_co2 = median(co2_flux_umolm2s_mean,na.rm=TRUE),
            mean_co2 = mean(co2_flux_umolm2s_mean,na.rm=TRUE),
            sd_co2 = sd(co2_flux_umolm2s_mean,na.rm=TRUE),
            cv_co2 = sd(co2_flux_umolm2s_mean,na.rm=TRUE)/mean(co2_flux_umolm2s_mean,na.rm=TRUE)*100,
            min_ch4 = min(ch4_flux_umolm2s_mean,na.rm = TRUE),
            max_ch4 = max(ch4_flux_umolm2s_mean,na.rm=TRUE),
            med_ch4 = median(ch4_flux_umolm2s_mean,na.rm=TRUE),
            mean_ch4 = mean(ch4_flux_umolm2s_mean,na.rm=TRUE),
            sd_ch4 = sd(ch4_flux_umolm2s_mean,na.rm=TRUE),
            cv_ch4 = sd(ch4_flux_umolm2s_mean,na.rm=TRUE)/mean(ch4_flux_umolm2s_mean,na.rm=TRUE)*100)

ch4_diff <- ggplot(ghg_fluxes,mapping=aes(DateTime,ch4_flux_umolm2s_mean))+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_errorbar(mapping=aes(x=DateTime,ymin=ch4_flux_umolm2s_mean-ch4_flux_umolm2s_sd,ymax=ch4_flux_umolm2s_mean+ch4_flux_umolm2s_sd),size=1,color="#4c8bfe")+
  geom_line(size = 1,color="#4c8bfe") +
  geom_point(size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  theme_classic(base_size = 15)

co2_diff <- ggplot(ghg_fluxes,mapping=aes(DateTime,co2_flux_umolm2s_mean))+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_errorbar(mapping=aes(x=DateTime,ymin=co2_flux_umolm2s_mean-co2_flux_umolm2s_sd,ymax=co2_flux_umolm2s_mean+co2_flux_umolm2s_sd),size=1,color="#4c8bfe")+
  geom_line(size = 1,color="#4c8bfe") +
  geom_point(size=3,color="#4c8bfe") +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-04-06"))+
  theme_classic(base_size = 15)

ggarrange(co2_diff,ch4_diff,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Diff_fluxes.jpg",width = 8, height=6, units="in",dpi=320)

### Load in Eddy Flux data ----
# Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
# Data corrected following FCR_Process_BD
eddy_flux <- read_csv("./Data/20210615_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Plot different percentiles
ggplot(eddy_flux)+
  geom_point(mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="Orig"))+
  #geom_line(mapping=aes(x=DateTime,y=ch4_flux_uStar_f,color="f"))+
  #geom_line(mapping=aes(x=DateTime,y=ch4_flux_uStar_fall,color="fall"))+
  #geom_line(mapping=aes(x=DateTime,y=ch4_flux_U80_f,color="U80"))+
  #geom_line(mapping=aes(x=DateTime,y=ch4_flux_U97.5_f,color="U97.5"))+
  #geom_line(mapping=aes(x=DateTime,y=ch4_flux_U2.5_f,color="U2.5"))
  geom_ribbon(mapping=aes(x=DateTime,y=ch4_flux_uStar_f,ymin=ch4_flux_uStar_f-ch4_flux_uStar_fsd,ymax=ch4_flux_uStar_f+ch4_flux_uStar_fsd,color="f"),alpha=0.3)

# Assess differences between measured and gap-filled data
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

ggsave("./Fig_Output/GapFilledComps.jpg",width = 10, height = 9, units="in",dpi=320)

# Aggregate to hourly and calculate the variability (SD) - following script for figures_BD
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

# Combine all stats into one
fcr_stats <- rbind(fcr_hourly_stats,fcr_daily_stats,fcr_weekly_stats,fcr_monthly_stats,ghg_stats)

write_csv(fcr_stats,"./Data/20210815_Stats.csv")

# daily means

dailynee <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="EC"),size = 1) +
  geom_errorbar(ghg_fluxes,mapping=aes(x=DateTime,y=co2_flux_umolm2s_mean,ymin=co2_flux_umolm2s_mean-co2_flux_umolm2s_sd,ymax=co2_flux_umolm2s_mean+co2_flux_umolm2s_sd,color="Diff"),size=1)+
  geom_point(ghg_fluxes,mapping=aes(DateTime,co2_flux_umolm2s_mean,color="Diff"),size=2) +
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
  geom_errorbar(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_flux_umolm2s_mean,ymin=ch4_flux_umolm2s_mean-ch4_flux_umolm2s_sd,ymax=ch4_flux_umolm2s_mean+ch4_flux_umolm2s_sd,color="Diff"),size=1)+
  geom_point(ghg_fluxes,mapping=aes(DateTime,ch4_flux_umolm2s_mean,color="Diff"),size=2) +
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

ggsave("./Fig_Output/Day_Month_time.jpg",width = 8, height=11, units="in",dpi=320)

## Plot weekly and monthly time series
co2_daily <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE),color="#E63946",size = 1) +
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
  theme_classic(base_size = 15)

ch4_daily <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4),color="#E63946",size = 1) +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
  theme_classic(base_size = 15)

ggarrange(co2_daily,ch4_daily,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Hourly_time.jpg",width = 8, height=8, units="in",dpi=320)

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

ggsave("./Fig_Output/Weekly_time.jpg",width = 8, height=8, units="in",dpi=320)

### Plotting diel variations: select data from noon and midnight ----
diel_flux <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f,NEE_uStar_orig,ch4_flux_uStar_orig) %>% 
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
                   Hour = mean(Hour)) %>% 
  mutate(season = ifelse(Month == 12 | Month == 1 | Month == 2, "Winter",
                         ifelse(Month == 3 | Month == 4 | Month == 5, "Spring",
                                ifelse(Month == 6 | Month == 7 | Month == 8 | Month == 9, "Summer",
                                       "Fall"))))

diel_agg$Date <- as.POSIXct(paste(diel_agg$Year, diel_agg$Month, diel_agg$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# Plot diel separation by date
co2_diel_time <- diel_agg %>% 
  drop_na(NEE_sd) %>% 
  ggplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(mapping=aes(x=Date,y=NEE,color=diel),size=1.5)+
  geom_line(mapping=aes(x=Date,y=NEE,color=diel),size=0.8)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd,fill=diel),alpha=0.3)+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  scale_fill_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_diel_time <- diel_agg %>% 
  drop_na(CH4_sd) %>% 
  ggplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(mapping=aes(x=Date,y=CH4,color=diel),size=1.5)+
  geom_line(mapping=aes(x=Date,y=CH4,color=diel),size=0.8)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd,fill=diel),alpha=0.3)+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  scale_fill_manual(breaks=c("Day","Night"),
                    values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_box_all <- ggplot(diel_agg,mapping=aes(x=diel,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.3)+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ch4_box_all <- ggplot(diel_agg,mapping=aes(x=diel,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(alpha=0.3)+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())

ggarrange(co2_diel_time,co2_box_all,ch4_diel_time,ch4_box_all,nrow=2,ncol=2,common.legend=TRUE,
          labels=c("A.","B.","C.","D."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Diel_Time.jpg",width=9,height=8,units="in",dpi=320)

# Calculate difference between medians for each season using Kruskal Wallis test
diel_agg$NEE[is.nan(diel_agg$NEE)] <-NA

diel_agg_stats <- diel_agg %>% 
  group_by(season,diel) %>% 
  summarise(p25_nee = quantile(NEE,0.25,na.rm=TRUE),
            med_nee = median(NEE,na.rm=TRUE),
            p75_nee = quantile(NEE,0.75,na.rm=TRUE),
            p25_ch4 = quantile(CH4,0.25,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            p75_ch4 = quantile(CH4,0.75,na.rm=TRUE))

diel_agg_stats_all <- diel_agg %>% 
  ungroup() %>% 
  group_by(diel) %>% 
  summarise(p25_nee = quantile(NEE,0.25,na.rm=TRUE),
            med_nee = median(NEE,na.rm=TRUE),
            p75_nee = quantile(NEE,0.75,na.rm=TRUE),
            p25_ch4 = quantile(CH4,0.25,na.rm=TRUE),
            med_ch4 = median(CH4,na.rm=TRUE),
            p75_ch4 = quantile(CH4,0.75,na.rm=TRUE)) %>% 
  mutate(season = "All") %>% 
  relocate(season,.before = p25_nee)

diel_agg_stats_comb <- rbind(diel_agg_stats,diel_agg_stats_all)

write_csv(diel_agg_stats_comb,"./Fig_output/20210815_diel_stats.csv")

diel_agg_summer <- diel_agg %>% 
  filter(season == "Summer")

diel_agg_fall <- diel_agg %>% 
  filter(season == "Fall")

diel_agg_winter <- diel_agg %>% 
  filter(season == "Winter")

diel_agg_spring <- diel_agg %>% 
  filter(season == "Spring")

kruskal.test(NEE ~ diel,data=diel_agg_summer)
kruskal.test(CH4 ~ diel,data=diel_agg_summer)

kruskal.test(NEE ~ diel,data=diel_agg_fall)
kruskal.test(CH4 ~ diel,data=diel_agg_fall)

kruskal.test(NEE ~ diel,data=diel_agg_winter)
kruskal.test(CH4 ~ diel,data=diel_agg_winter)

kruskal.test(NEE ~ diel,data=diel_agg_spring)
kruskal.test(CH4 ~ diel,data=diel_agg_spring) # Statistically significant

kruskal.test(NEE ~ diel,data=diel_agg)
kruskal.test(CH4 ~ diel,data=diel_agg)

# Plot diel fluxes by season
diel_agg$season <- factor(diel_agg$season, levels=c("Spring", "Summer", "Fall", "Winter"))

co2_diel <- ggplot(diel_agg,mapping=aes(x=season,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position = position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ch4_diel <- ggplot(diel_agg,mapping=aes(x=season,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position = position_jitterdodge(),alpha=0.3)+
  annotate("text",label = "*", x = 1, y = 0.075, size = 7)+
  scale_color_manual(breaks=c("Day","Night"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size=17)+
  theme(legend.title=element_blank())

ggarrange(co2_diel,ch4_diel,nrow=1,ncol=2,common.legend = TRUE, labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/GHG_Diel_all.jpg",width=9,height=4,units="in",dpi=320)

### Thinking about other visualizations ----
# Plotting cumulative Co2 and Ch4 throughout the study period
# From 2020-04-04 to 2021-04-03
eddy_flux_sum <- eddy_flux %>% 
  select(DateTime,ch4_flux_uStar_orig,ch4_flux_uStar_f,ch4_flux_uStar_fsd,ch4_flux_uStar_fall,ch4_flux_uStar_fqc,
         ch4_flux_U97.5_f,NEE_uStar_orig,NEE_uStar_f,NEE_uStar_fsd,NEE_uStar_fall,NEE_uStar_fqc,NEE_U97.5_f) %>%
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v)))

# Calculate cumulative fluxes and uncertainty for ghg data
cumulative_fluxes <- fluxes_all %>% 
  filter(DateTime >= as.POSIXct("2020-04-06 01:00:00") & DateTime < as.POSIXct("2021-04-06 01:00:00")) %>% 
  mutate(co2_sum_g_m2_d = cumsum(co2_flux_umolm2s_mean*1800*12.01/1000000)) %>% 
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_umolm2s_mean*1800*12.01/1000000)) %>% 
  mutate(co2_v = (co2_flux_umolm2s_sd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_umolm2s_sd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v)))

# Compare cumulative fluxes
sum_co2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(eddy_flux_sum,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,color="EC"),size=1)+
  geom_ribbon(eddy_flux_sum,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,ymin=co2_sum_g_m2_d-co2_sum_sd,ymax=co2_sum_g_m2_d+co2_sum_sd),fill="#E63946",alpha=0.3)+
  geom_line(cumulative_fluxes,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,color="Diff"),size=1)+
  geom_ribbon(cumulative_fluxes,mapping=aes(x=DateTime,y=co2_sum_g_m2_d,ymin=co2_sum_g_m2_d-co2_sum_sd,ymax=co2_sum_g_m2_d+co2_sum_sd),fill="#4c8bfe",alpha=0.3)+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  scale_color_manual(breaks=c("EC","Diff"), 
                     values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  labs(color="")

sum_ch4 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(eddy_flux_sum,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,color="EC"),size=1)+
  geom_ribbon(eddy_flux_sum,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,ymin=ch4_sum_g_m2_d-ch4_sum_sd,ymax=ch4_sum_g_m2_d+ch4_sum_sd),fill="#E63946",alpha=0.3)+
  geom_line(cumulative_fluxes,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,color="Diff"),size=1)+
  geom_ribbon(cumulative_fluxes,mapping=aes(x=DateTime,y=ch4_sum_g_m2_d,ymin=ch4_sum_g_m2_d-ch4_sum_sd,ymax=ch4_sum_g_m2_d+ch4_sum_sd),fill="#4c8bfe",alpha=0.3)+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  scale_color_manual(breaks=c("EC","Diff"), 
                     values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  labs(color="")

ggarrange(sum_co2,sum_ch4,nrow=1,ncol=2,common.legend = TRUE,labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Summed_Fluxes_v2.jpg",width = 9, height=4.5, units="in",dpi=320)

# Calculate cumulative fluxes in summer (Jun, Jul, Aug, Sep) vs. winter (Dec, Jan, Feb, Mar)
# Summer (Jun - Sept) = 384.1 g C-CO2 m2 d; 0.7275 g C-CH4 m2 d
eddy_flux_sum %>% 
  filter(DateTime == as.POSIXct("2020-05-31 20:00:00") | DateTime == as.POSIXct("2020-09-30 20:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d)

# Winter (Dec - Mar) = 59 g C-CO2 m2 d; 0.025 g C-CH4 m2 d
eddy_flux_sum %>% 
  filter(DateTime == as.POSIXct("2020-12-31 20:00:00") | DateTime == as.POSIXct("2021-03-31 20:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d)

# Check winter fluxes vs. summer fluxes
eddy_flux_sum %>% 
  filter(DateTime <= as.POSIXct("2020-11-01")) %>% 
  summarise(co2 = mean(NEE_uStar_f,na.rm=TRUE),
            ch4 = mean(ch4_flux_uStar_f,na.rm=TRUE))

eddy_flux_sum %>% 
  filter(DateTime > as.POSIXct("2020-11-01")) %>% 
  summarise(co2 = mean(NEE_uStar_f,na.rm=TRUE),
            ch4 = mean(ch4_flux_uStar_f,na.rm=TRUE))

eddy_flux_sum %>% 
  filter(DateTime < as.POSIXct("2020-10-01")) %>% 
  summarise(co2 = mean(NEE_uStar_f,na.rm=TRUE),
            ch4 = mean(ch4_flux_uStar_f,na.rm=TRUE))

cumulative_fluxes %>% 
  filter(DateTime < as.POSIXct("2020-10-01")) %>% 
  summarise(co2 = mean(co2_flux_umolm2s_mean,na.rm=TRUE),
            ch4 = mean(ch4_flux_umolm2s_mean,na.rm=TRUE))

eddy_flux_sum %>% 
  filter(DateTime > as.POSIXct("2020-10-01")) %>% 
  summarise(co2 = mean(NEE_uStar_f,na.rm=TRUE),
            ch4 = mean(ch4_flux_uStar_f,na.rm=TRUE))

cumulative_fluxes %>% 
  filter(DateTime > as.POSIXct("2020-10-01")) %>% 
  summarise(co2 = mean(co2_flux_umolm2s_mean,na.rm=TRUE),
            ch4 = mean(ch4_flux_umolm2s_mean,na.rm=TRUE))


### Think about winter variability (especially with ice!)
ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))
ice_on <- ice %>% 
  filter(Date>"2020-01-01" & IceOn == 1)
ice_off <- ice %>% 
  filter(Date>"2020-01-01" & IceOff == 1)

# Calculate average flux for each ice on/off period
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

kruskal.test(NEE_uStar_orig ~ diel, data=ice_off_1)
kruskal.test(ch4_flux_uStar_orig ~ diel, data=ice_off_1)

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

kruskal.test(NEE_uStar_orig ~ diel, data=ice_on_1)
kruskal.test(ch4_flux_uStar_orig ~ diel, data=ice_on_1)

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

kruskal.test(NEE_uStar_orig ~ diel, data=ice_off_2)
kruskal.test(ch4_flux_uStar_orig ~ diel, data=ice_off_2)

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

kruskal.test(NEE_uStar_orig ~ diel, data=ice_on_2)
kruskal.test(ch4_flux_uStar_orig ~ diel, data=ice_on_2)

ice_all <- rbind(ice_off_1,ice_on_1,ice_off_2,ice_on_2)

# Test for significance for all ice on vs. ice off
kruskal.test(NEE_uStar_orig ~ ice, data=ice_all)
# Not statistically significant (p = 0.592)
kruskal.test(ch4_flux_uStar_orig ~ ice, data=ice_all)
# Not statistically significant (p = 0.1085)

# Ice_1
ice_1 <- ice_all %>% 
  filter(ice_period %in% c("Off_1","On_1"))

kruskal.test(NEE_uStar_orig ~ ice, data=ice_1)
kruskal.test(ch4_flux_uStar_orig ~ ice, data=ice_1)

# Ice_2
ice_2 <- ice_all %>% 
  filter(ice_period %in% c("Off_2","On_2"))

kruskal.test(NEE_uStar_orig ~ ice, data=ice_2)
kruskal.test(ch4_flux_uStar_orig ~ ice, data=ice_2)

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

write.csv(ice_stats_final,"./Fig_output/20210816_Ice_stats.csv")

# Daily means in winter (ice on/off)
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

ggsave("./Fig_Output/Ice_on_off_v3.jpg",width = 8, height=6, units="in",dpi=320)

# Plot seprated by diel
ice_all$ice_period <- factor(ice_all$ice_period, levels=c("Off_1", "On_1", "Off_2", "On_2"))

ice_co2_d <- ggplot(ice_all,mapping=aes(x=ice_period,y=NEE_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  theme_classic(base_size = 15)+
  ylim(-10,10)+
  theme(legend.title=element_blank())+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

ice_ch4_d <- ggplot(ice_all,mapping=aes(x=ice_period,y=ch4_flux_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "*", x = "On_1", y = 0.025, size = 5)+
  annotate("text",label = "*", x = "Off_2", y = 0.025, size = 5)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylim(-0.03,0.03)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(ice_co2_d,ice_ch4_d,nrow=1,ncol=2,common.legend = TRUE,
          labels=c("A.","B."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Ice_on_off_diel.jpg",width = 10, height=5, units="in",dpi=320)

### Let's get catwalk data in hand ----
# To start thinking about environmental variables
# From EDI (up to 2020) on 03 June 2021
# NOTE: Data is AHEAD by 4 hours (in GMT!)
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/271/5/c1b1f16b8e3edbbff15444824b65fe8f" 
#infile1 <- paste0(getwd(),"/Data/Catwalk_EDI_2020.csv")
#download.file(inUrl1,infile1,method="curl")

catwalk_edi <- read.csv("./Data/Catwalk_EDI_2020.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="GMT")))

catwalk_edi_est <- lubridate::with_tz(catwalk_edi,"EST")

catwalk_edi_est <- catwalk_edi_est %>% 
  filter(DateTime >= "2020-01-01") %>% 
  select(DateTime,ThermistorTemp_C_surface:EXOfDOM_RFU_1)

# Pull most up-to-date catwalk data from github
pacman::p_load("RCurl","tidyverse","lubridate", "plotly", "magrittr")
folder <- "./Data"

# Load in data from Git
# Downloaded: 03 June 2021
#download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data/CAT_MaintenanceLog.txt",paste0(folder, "/CAT_MaintenanceLog_2020.txt"))
#download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data/Catwalk.csv",paste0(folder, "/Catwalk_2020.csv"))

# Running QA/QC from catwalk_EDI_QAQC_all_variables
data_file <- paste0(folder, '/Catwalk_2020.csv')
maintenance_file <- paste0(folder, "/CAT_MaintenanceLog_2020.txt")

catdata <- read.csv(data_file,skip=1) 
catdata$DateTime<-as.POSIXct(catdata$TIMESTAMP,format = "%Y-%m-%d %H:%M:%S")
catdata <- catdata[!duplicated(catdata$TIMESTAMP), ]
colnames(catdata)[colnames(catdata)=="Lvl_psi"] <- "Lvl_psi_9"
colnames(catdata)[colnames(catdata)=="LvlTemp_c_9"] <- "LvlTemp_C_9"

# subset file to only unpublished data
catdata_flag <- catdata[catdata$TIMESTAMP>"2021-01-01",]

catdata_flag <- catdata_flag[-1,]

catdata_flag <- catdata_flag %>% 
  rename(EXODO_mgL_1 = doobs_1,EXODOsat_percent_1 = dosat_1)

# now fix the negative DO values
catdata_flag <- catdata_flag %>%
  mutate(Flag_DO_1 = NA) %>% 
  mutate(Flag_DO_1 = ifelse(EXODO_mgL_1 < 0 | EXODOsat_percent_1 <0, 3, Flag_DO_1), #and for 1m
         EXODO_mgL_1 = ifelse(EXODO_mgL_1 < 0, 0, EXODO_mgL_1),
         EXODOsat_percent_1 = ifelse(EXODOsat_percent_1 <0, 0, EXODOsat_percent_1),
         Flag_DO_1 = ifelse(is.na(EXODO_mgL_1),7,Flag_DO_1))

# chl and phyco qaqc ----
# perform qaqc on the entire dataset for chl and phyco

catdata_flag <- catdata_flag %>% 
  rename(EXOChla_ugL_1 = Chla_1,EXOBGAPC_ugL_1 = BGAPC_1)

catdata_flag <- catdata_flag %>% 
  rename(EXOChla_RFU_1 = Chla_RFU_1)

# assign standard deviation thresholds
sd_4 <- 4*sd(catdata_flag$EXOChla_ugL_1, na.rm = TRUE)
threshold <- sd_4
sd_4_phyco <- 4*sd(catdata_flag$EXOBGAPC_ugL_1, na.rm = TRUE)
threshold_phyco <- sd_4_phyco

#chl_ugl <- ggplot(data = catdata_all, aes(x = DateTime, y = EXOChla_ugL_1)) +
#  geom_point() +
#  geom_hline(yintercept = sd_4)
#ggplotly(chl_ugl)

# QAQC on major chl outliers using DWH's method: datapoint set to NA if data is greater than 4*sd different from both previous and following datapoint
catdata_flag <- catdata_flag %>% 
  mutate(Chla = lag(EXOChla_ugL_1, 0),
         Chla_lag1 = lag(EXOChla_ugL_1, 1),
         Chla_lead1 = lead(EXOChla_ugL_1, 1)) %>%   #These mutates create columns for current fDOM, fDOM before and fDOM after. These are used to run ifelse QAQC loops
  mutate(EXOChla_ugL_1 = ifelse(Chla < 0 & !is.na(Chla), 0, EXOChla_ugL_1)) %>% 
  mutate(EXOChla_RFU_1 = ifelse(Chla < 0 & !is.na(Chla), 0, EXOChla_RFU_1)) %>%
  mutate(Chla = ifelse(Chla == "NAN", NA, Chla)) %>% 
  mutate(Chla_lag1 = ifelse(Chla_lag1 == "NAN", NA, Chla_lag1)) %>% 
  mutate(Chla_lead1 = ifelse(Chla_lead1 == "NAN", NA, Chla_lead1)) %>% 
  mutate(Chla = as.numeric(Chla)) %>% 
  mutate(Chla_lag1 = as.numeric(Chla_lag1)) %>% 
  mutate(Chla_lead1 = as.numeric(Chla_lead1)) %>% 
  mutate(EXOChla_ugL_1 = ifelse((abs(Chla_lag1 - Chla) > (threshold))  & (abs(Chla_lead1 - Chla) > (threshold) & !is.na(Chla)), 
                                NA, EXOChla_ugL_1)) %>%   
  mutate(EXOChla_RFU_1 = ifelse((abs(Chla_lag1 - Chla) > (threshold))  & (abs(Chla_lead1 - Chla) > (threshold) & !is.na(Chla)), 
                                NA, EXOChla_RFU_1)) %>% 
  select(-Chla, -Chla_lag1, -Chla_lead1)

# fdom qaqc
# QAQC from DWH to remove major outliers from fDOM data that are 2 sd's greater than the previous and following datapoint
# QAQC done on 2018-2020 dataset
catdata_flag <- catdata_flag %>% 
  rename(EXOfDOM_QSU_1 = fDOM_QSU_1, EXOfDOM_RFU_1 = fDOM_RFU_1)

sd_fDOM <- sd(catdata_flag$EXOfDOM_QSU_1, na.rm = TRUE) #deteriming the standard deviation of fDOM data 

#fDOM_pre_QAQC <- ggplot(data = catdata_all, aes(x = DateTime, y = EXOfDOM_QSU_1)) +
#  geom_point()+
#  ggtitle("fDOM (QSU) pre QAQC")
#ggplotly(fDOM_pre_QAQC)

catdata_flag <- catdata_flag %>% 
  mutate(fDOM = lag(EXOfDOM_QSU_1, 0),
         fDOM_lag1 = lag(EXOfDOM_QSU_1, 1),
         fDOM_lead1 = lead(EXOfDOM_QSU_1, 1)) %>%  #These mutates create columns for current fDOM, fDOM before and fDOM after. These are used to run ifelse QAQC loops
  mutate(EXOfDOM_QSU_1 = ifelse(fDOM < 0 & !is.na(fDOM), NA, EXOfDOM_QSU_1),
         EXOfDOM_RFU_1 = ifelse(fDOM < 0 & !is.na(fDOM), NA, EXOfDOM_RFU_1)) %>% #These mutates are QAQCing for negative fDOM QSU values and setting these to NA and making a flag for these. This was done outside of the 2 sd deviation rule because there were two negative points in a row and one was not removed with the follwoing if else statements. 
  mutate(fDOM = ifelse(fDOM == "NAN", NA, fDOM)) %>% 
  mutate(fDOM_lag1 = ifelse(fDOM_lag1 == "NAN", NA, fDOM_lag1)) %>% 
  mutate(fDOM_lead1 = ifelse(fDOM_lead1 == "NAN", NA, fDOM_lead1)) %>% 
  mutate(fDOM = as.numeric(fDOM)) %>% 
  mutate(fDOM_lag1 = as.numeric(fDOM_lag1)) %>% 
  mutate(fDOM_lead1 = as.numeric(fDOM_lead1)) %>% 
  mutate(EXOfDOM_QSU_1 = ifelse(
    (abs(fDOM_lag1 - fDOM) > (2*sd_fDOM)) & (abs(fDOM_lead1 - fDOM) > (2*sd_fDOM)  & !is.na(fDOM)), NA, EXOfDOM_QSU_1
  )) %>%  #QAQC to remove outliers for QSU fDOM data 
  mutate(EXOfDOM_RFU_1 = ifelse(
    (abs(fDOM_lag1 - fDOM) > (2*sd_fDOM)) & (abs(fDOM_lead1 - fDOM) > (2*sd_fDOM)  & !is.na(fDOM)), NA, EXOfDOM_RFU_1
  )) %>% #QAQC to remove outliers for RFU fDOM data
  select(-fDOM, -fDOM_lag1, -fDOM_lead1)  #This removes the columns used to run ifelse statements since they are no longer needed. 

#Deal with when the sensors were up
maint = read.csv(paste0(folder, "/CAT_MaintenanceLog_2020.txt"))
maint = maint[!grepl("EXO",maint$parameter),] #creating file "maint" with all sensor string maintenance
maint <- maint %>% 
  filter(parameter == " All_Cat")
maint = maint%>%
  filter(!colnumber %in% c(" c(24:26)"," 40"," 41"))
clean_start<-as.POSIXct(maint$TIMESTAMP_start, tz = "EST")
clean_end <- as.POSIXct(maint$TIMESTAMP_end, tz = "EST")

ADJ_PERIOD = 2*60*60 #amount of time to stabilization after cleaning in seconds

for (i in 1:length(clean_start)){ #Set all data during cleaning and for ADJ_PERIOD after to NA
  catdata_flag$EXODO_mgL_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<(clean_end[i]+ADJ_PERIOD)] <- NA
  catdata_flag$EXODOsat_percent_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- NA
  catdata_flag$EXO_wtr_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- NA
  catdata_flag$SpCond_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- NA
  catdata_flag$EXOChla_ugL_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- NA
  catdata_flag$EXOfDOM_RFU_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- NA
  catdata_flag$Flag_DO_1[catdata_flag$DateTime>clean_start[i]&catdata_flag$DateTime<clean_end[i]+ADJ_PERIOD] <- 1
}

# Select columns of interest (Exo columns)
catwalk_git <- catdata_flag %>% 
  select(DateTime,wtr_surface,wtr_1,wtr_2,wtr_3,wtr_4,wtr_5,wtr_6,wtr_7,wtr_8,wtr_9,
         EXO_wtr_1,SpCond_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1,EXODOsat_percent_1) %>% 
  rename(ThermistorTemp_C_surface = wtr_surface,ThermistorTemp_C_1 = wtr_1,ThermistorTemp_C_2 = wtr_2,
         ThermistorTemp_C_3 = wtr_3, ThermistorTemp_C_4 = wtr_4, ThermistorTemp_C_5 = wtr_5,
         ThermistorTemp_C_6 = wtr_6, ThermistorTemp_C_7 = wtr_7, ThermistorTemp_C_8 = wtr_8,
         ThermistorTemp_C_9 = wtr_9, EXOTemp_C_1 = EXO_wtr_1,EXOSpCond_uScm_1 = SpCond_1)

catwalk_edi_est_sel <- catwalk_edi_est %>% 
  select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9, EXOTemp_C_1, EXOSpCond_uScm_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1,
         EXODOsat_percent_1)

# Combine catwalk data from EDI and from GitHub
catwalk_all <- rbind(catwalk_edi_est_sel,catwalk_git)

catwalk_all$EXOTemp_C_1 <- as.numeric(catwalk_all$EXOTemp_C_1)
catwalk_all$EXODO_mgL_1 <- as.numeric(catwalk_all$EXODO_mgL_1)
catwalk_all$EXODOsat_percent_1 <- as.numeric(catwalk_all$EXODOsat_percent_1)
catwalk_all$EXOSpCond_uScm_1 <- as.numeric(catwalk_all$EXOSpCond_uScm_1)
catwalk_all$EXOChla_ugL_1 <- as.numeric(catwalk_all$EXOChla_ugL_1)
catwalk_all$EXOfDOM_RFU_1 <- as.numeric(catwalk_all$EXOfDOM_RFU_1)
catwalk_all$ThermistorTemp_C_surface <- as.numeric(catwalk_all$ThermistorTemp_C_surface)
catwalk_all$ThermistorTemp_C_1 <- as.numeric(catwalk_all$ThermistorTemp_C_1)
catwalk_all$ThermistorTemp_C_2 <- as.numeric(catwalk_all$ThermistorTemp_C_2)
catwalk_all$ThermistorTemp_C_3 <- as.numeric(catwalk_all$ThermistorTemp_C_3)
catwalk_all$ThermistorTemp_C_4 <- as.numeric(catwalk_all$ThermistorTemp_C_4)
catwalk_all$ThermistorTemp_C_5 <- as.numeric(catwalk_all$ThermistorTemp_C_5)
catwalk_all$ThermistorTemp_C_6 <- as.numeric(catwalk_all$ThermistorTemp_C_6)
catwalk_all$ThermistorTemp_C_7 <- as.numeric(catwalk_all$ThermistorTemp_C_7)
catwalk_all$ThermistorTemp_C_8 <- as.numeric(catwalk_all$ThermistorTemp_C_8)
catwalk_all$ThermistorTemp_C_9 <- as.numeric(catwalk_all$ThermistorTemp_C_9)

# Save catwalk data:
write.csv(catwalk_all,"./Data/Catwalk_all.csv")

catwalk_mean <- catwalk_all %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  summarise_all(mean,na.rm=TRUE)

