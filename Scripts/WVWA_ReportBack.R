### Script to visualize/plot Eddy Flux data for WVWA report-back on 26 May 2021
# Using data corrected by Brenda D.
# 18 May 2021, A Hounshell

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate)

# Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
flux <- read_csv("./Data/processed_data_upto_2021-05-06_BD.csv")

flux <- flux %>% 
  select(DateTime,NEE_uStar_f,ch4_flux_uStar_f)

# From Brenda: There are a lot of variables, but the ones that show the final 
# fluxes are NEE_uStar_f and ch4_flux_uStar_f (ustar filtered and gapfilled co2 and ch4 fluxes, respectively)

# Let's just start by plotting the data + turn-over schedule?
flux_year <- flux %>% 
  filter(DateTime>=as.POSIXct("2020-05-06")&DateTime<as.POSIXct("2021-05-06")) %>% 
  select(NEE_uStar_f,ch4_flux_uStar_f) %>% 
  summarize_all(sum)*85000/1000000*60*30*12.111/1000 # convert to kg C/year (summed each 30 min interval*area) converted to 30 min intervals

flux_co2 <- ggplot(flux,mapping=aes(x=DateTime,y=NEE_uStar_f))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  annotate(geom="text",x=as.POSIXct("2020-10-12"),y=50, label = "Turnover")+
  annotate(geom="text",x=as.POSIXct("2021-03-01"),y=50, label = expression(paste("52,626 kg C")))+
  theme_classic(base_size = 15)

flux_ch4 <- ggplot(flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_f))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Time")+
  annotate(geom="text",x=as.POSIXct("2020-10-12"),y=0.10, label = "Turnover")+
  annotate(geom="text",x=as.POSIXct("2021-03-01"),y=0.10, label = expression(paste("69 kg C")))+
  theme_classic(base_size = 15)

ggarrange(flux_co2,flux_ch4,nrow=2,ncol=1)

ggsave("./Fig_Output/WVWA_Fluxes.jpg",width = 10, height=7, units="in",dpi=320)


# For CCC: plot UGGA results :) ----
# Load in *preliminary* UGGA data
ugga <- read_csv("./Data/20210219_Flux_Output.csv")

ugga <- ugga %>% 
  mutate(Date = as.POSIXct(strptime(Date, "%Y-%m-%d", tz="EST"))) %>% 
  mutate(DateTime = paste(Date, Start_time)) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

flux_co2 <- ggplot(flux,mapping=aes(x=DateTime,y=NEE_uStar_f))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line()+
  geom_point(ugga,mapping=aes(x=DateTime,y=co2_flux_umolCm2s,color="UGGA"),size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  xlim(as.POSIXct("2020-07-05"),as.POSIXct("2020-11-02"))+
  theme_classic(base_size = 15)

flux_ch4 <- ggplot(flux,mapping=aes(x=DateTime,y=ch4_flux_uStar_f))+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="navyblue")+ #Turnover FCR; operationally defined
  geom_line()+
  geom_point(ugga,mapping=aes(x=DateTime,y=ch4_flux_umolCm2s,color="UGGA"),size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Time")+
  annotate(geom="text",x=as.POSIXct("2020-10-12"),y=0.10, label = "Turnover")+
  xlim(as.POSIXct("2020-07-05"),as.POSIXct("2020-11-02"))+
  theme_classic(base_size = 15)

ggarrange(flux_co2,flux_ch4,nrow=2,ncol=1)

ggsave("./Fig_Output/Fluxes_wUGGA.jpg",width = 10, height=7, units="in",dpi=320)
