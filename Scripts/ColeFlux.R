### Calculating diffusive fluxes from dissolved GHGs using the Cole method
### Following McClure et al. 2018
### 21 May 2021, A Hounshell

# Following code from 2019: for calculating fluxes
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

# Calculate fluxes: in umol/m2/s
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
  summarise_all(funs(mean,sd),na.rm=TRUE) %>% 
  select(-c(Rep_mean,Rep_sd))

# Plot to check?
ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,ch4_flux_umolm2s_mean))+
  #geom_ribbon(fluxes_all,mapping=aes(DateTime,ymin = ch4_umolm2s-ch4_umolm2s_sd,ymax = ch4_umolm2s+ch4_umolm2s_sd),fill="grey")+
  theme_classic(base_size =15 )

ggplot()+
  geom_line(fluxes_all,mapping=aes(DateTime,co2_flux_umolm2s_mean))+
  theme_classic(base_size=15)

# Select dates where we actually have GHG concentrations
ghg_fluxes <- left_join(ghg_1,fluxes_all,by="DateTime")

### Load in Eddy Flux data ----
# Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
eddy_flux <- read_csv("./Data/20210607_EC_processed.csv")

eddy_flux <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f)

# Plot Eddy flux data and GHG flux data
flux_co2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(ghg_fluxes,mapping=aes(x=DateTime,y=co2_flux_umolm2s_mean,color="GHGs"))+
  geom_line(ghg_fluxes,mapping=aes(x=DateTime,y=co2_flux_umolm2s_mean,color="GHGs"))+
  geom_errorbar(ghg_fluxes,mapping=aes(x=DateTime,ymin=co2_flux_umolm2s_mean - co2_flux_umolm2s_sd, ymax = co2_flux_umolm2s_mean+co2_flux_umolm2s_sd,color="GHGs"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-29"))+
  theme_classic(base_size = 15)

flux_ch4 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_flux_umolm2s_mean,color="GHGs"))+
  geom_line(ghg_fluxes,mapping=aes(x=DateTime,y=ch4_flux_umolm2s_mean,color="GHGs"))+
  geom_errorbar(ghg_fluxes,mapping=aes(x=DateTime,ymin=ch4_flux_umolm2s_mean - ch4_flux_umolm2s_sd, ymax = ch4_flux_umolm2s_mean+ch4_flux_umolm2s_sd,color="GHGs"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Time")+
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-29"))+
  theme_classic(base_size = 15)

ggarrange(flux_co2,flux_ch4,nrow=2,ncol=1)

ggsave("./Fig_Output/Fluxes.jpg",width = 10, height=7, units="in",dpi=320)

### Thinking about other visualizations ----
# Plotting cumulative Co2 and Ch4 throughout the study period
# From 2020-04-04 to 2021-04-03
eddy_flux <- eddy_flux %>% 
  mutate(co2_sum = NA) %>% 
  mutate(ch4_sum = NA)

for (i in 1:length(eddy_flux$DateTime)){
  eddy_flux$co2_sum[i] <- sum(eddy_flux$NEE_uStar_f[1:i])
}

for (i in 1:length(eddy_flux$DateTime)){
  eddy_flux$ch4_sum[i] <- sum(eddy_flux$ch4_flux_uStar_f[1:i])
}

cumulative_fluxes <- fluxes_all %>% 
  filter(DateTime >= as.POSIXct("2020-04-04 01:00:00")) %>% 
  mutate(co2_sum = NA) %>% 
  mutate(ch4_sum = NA)

for (i in 1:length(cumulative_fluxes$DateTime)){
  cumulative_fluxes$co2_sum[i] <- sum(cumulative_fluxes$co2_flux_umolm2s_mean[1:i],na.rm=TRUE)
}

for (i in 1:length(cumulative_fluxes$DateTime)){
  cumulative_fluxes$ch4_sum[i] <- sum(cumulative_fluxes$ch4_flux_umolm2s_mean[1:i],na.rm=TRUE)
}

# Compare cumulative fluxes
sum_co2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=co2_sum/1e6*44.009*60*30,color="EC"),size=1)+
  geom_line(cumulative_fluxes,mapping=aes(x=DateTime,y=co2_sum/1e6*44.009*60*30,color="GHGs"),size=1)+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)

sum_ch4 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=ch4_sum/1e6*16.04*60*30,color="EC"),size=1)+
  geom_line(cumulative_fluxes,mapping=aes(x=DateTime,y=ch4_sum/1e6*16.04*60*30,color="GHGs"),size=1)+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  xlim(as.POSIXct("2020-04-04"),as.POSIXct("2021-04-05"))+
  theme_classic(base_size = 15)

ggarrange(sum_co2,sum_ch4,nrow=1,ncol=2,common.legend = TRUE)

ggsave("./Fig_Output/Summed_Fluxes.jpg",width = 10, height=5, units="in",dpi=320)

### Think about winter variability (especially with ice!)
ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))
ice_on <- ice %>% 
  filter(Date>"2020-01-01" & IceOn == 1)
ice_off <- ice %>% 
  filter(Date>"2020-01-01" & IceOff == 1)

# Create graph to look at ice on/off
winter_ch4 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dashed", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dashed", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dashed", color="red")+
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  ylim(-0.02,0.03)+
  theme_classic(base_size = 15)

winter_co2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2020-12-27"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2020-12-30"), linetype = "dashed", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dashed", color="red")+
  geom_vline(xintercept = as.POSIXct("2021-02-11"), linetype = "dashed", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dashed", color="red")+
  geom_line(eddy_flux,mapping=aes(x=DateTime,y=NEE_uStar_f))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-03-01"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  ylim(-7,7)+
  theme_classic(base_size = 15)

ggarrange(winter_co2,winter_ch4,nrow=2,ncol=1)

ggsave("./Fig_Output/Winter_Fluxes.jpg",width = 10, height=7, units="in",dpi=320)

# Calculate average flux for each ice on/off period
ice_off_1 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-19 19:00:00" & DateTime < "2020-12-26 19:00:00") %>% 
  mutate(ice_period = "1") %>% 
  mutate(ice = "off")

ice_on_1 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-26 19:00:00" & DateTime < "2020-12-29 19:00:00")%>% 
  mutate(ice_period = "2") %>% 
  mutate(ice = "on")

ice_off_2 <- eddy_flux %>% 
  filter(DateTime >= "2020-12-29 19:00:00" & DateTime < "2021-01-09 19:00:00")%>% 
  mutate(ice_period = "3") %>% 
  mutate(ice = "off")

ice_on_2 <- eddy_flux %>% 
  filter(DateTime >= "2021-01-09 19:00:00" & DateTime < "2021-02-09 19:00:00")%>% 
  mutate(ice_period = "4") %>% 
  mutate(ice = "on")

ice_off_3 <- eddy_flux %>% 
  filter(DateTime >= "2021-02-09 19:00:00" & DateTime < "2021-02-11 19:00:00")%>% 
  mutate(ice_period = "5") %>% 
  mutate(ice = "off")

ice_on_3 <- eddy_flux %>% 
  filter(DateTime >= "2021-02-11 19:00:00" & DateTime < "2021-02-23 19:00:00")%>% 
  mutate(ice_period = "6") %>% 
  mutate(ice = "on")

#ice_off_4 <- eddy_flux %>% 
#  filter(DateTime >= "2021-02-23 19:00:00" & DateTime < "2021-03-01 19:00:00")%>% 
#  mutate(ice_period = "7") %>% 
#  mutate(ice = "off")

ice_all <- rbind(ice_off_1,ice_on_1,ice_off_2,ice_on_2,ice_off_3,ice_on_3)

# Plot?
ice_co2 <- ggplot(ice_all,mapping=aes(x=ice_period,y=NEE_uStar_f,color=ice))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot()+
  theme_classic(base_size = 15)

ice_ch4 <- ggplot(ice_all,mapping=aes(x=ice_period,y=ch4_flux_uStar_f,color=ice))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot()+
  theme_classic(base_size = 15)

ggarrange(ice_co2,ice_ch4,nrow=1,ncol=2,common.legend = TRUE)

ggsave("./Fig_Output/Ice_Fluxes.jpg",width = 10, height=5, units="in",dpi=320)

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
  select(DateTime,EXOTemp_C_1,EXOSpCond_uScm_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1)

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
  select(DateTime,EXO_wtr_1,SpCond_1,EXODO_mgL_1,EXOChla_ugL_1,EXOfDOM_RFU_1) %>% 
  rename(EXOTemp_C_1 = EXO_wtr_1,EXOSpCond_uScm_1 = SpCond_1)

# Combine catwalk data from EDI and from GitHub
catwalk_all <- rbind(catwalk_edi_est,catwalk_git)

catwalk_all$EXOTemp_C_1 <- as.numeric(catwalk_all$EXOTemp_C_1)
catwalk_all$EXOSpCond_uScm_1 <- as.numeric(catwalk_all$EXOSpCond_uScm_1)
catwalk_all$EXOChla_ugL_1 <- as.numeric(catwalk_all$EXOChla_ugL_1)
catwalk_all$EXOfDOM_RFU_1 <- as.numeric(catwalk_all$EXOfDOM_RFU_1)

# Plot catwalk temp to check?
ggplot(catwalk_all,mapping=aes(x=DateTime,y=EXOTemp_C_1))+
  geom_line()+
  theme_classic(base_size=15)

ggplot(catwalk_all,mapping=aes(x=DateTime,y=EXOChla_ugL_1))+
  geom_line()+
  theme_classic(base_size=15)

ggplot(catwalk_all,mapping=aes(x=DateTime,y=EXOfDOM_RFU_1))+
  geom_line()+
  theme_classic(base_size=15)

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

# CO2
co2_atm <- nc_open('./Data/co2_sct_tower-insitu_1_ccgg_HourlyData.nc')

# Export metadata out as text file
{sink('./Data/co2_sct.txt')
  print(ch4_atm)
  sink()}

# Extract date/time, elevation, and CH4 (ppb)
datetime <- ncvar_get(co2_atm,"time")
datetime <- as.POSIXct(datetime, origin = "1970-01-01", tz="UTC")
elevation <- ncvar_get(co2_atm,"elevation")
co2 <- ncvar_get(co2_atm,"value")

# Combine datetime and CH4 into one data frame
carbondioxide <- cbind.data.frame(datetime,co2)
carbondioxide <- carbondioxide %>% 
  mutate(datetime = format(as.POSIXct(datetime,"%Y-%m-%d %H:%M:%S",tz="UTC"),format="%Y-%m-%d"))

