# Script to plot voltage and flow pump rate for Eddy Flux system
# To check potential power issues
# Data from Brenda D. on 16 Nov 20
# A Hounshell, 16 Nov 20

setwd("C:/Users/ahoun/OneDrive/Desktop/EddyFlux")

pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load data from Brenda
data <- read.csv("./Data/FCR_2020-11-10_upto.csv") %>% 
  mutate(datetime = paste(date,time)) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%m/%d/%Y %H:%M:%S", tz="EST")))

# Plot vin_sf_mean
ggplot(data,mapping=aes(x=datetime,y=vin_sf_mean))+
  geom_line()+
  ylim(11,11.6)+
  geom_vline(xintercept = as.POSIXct("2020-05-11"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-18"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-25"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-03"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-04"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-12"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-06-22"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-06"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-20"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-27"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-03"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-10"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-17"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-16"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-21"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-25"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-30"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-06"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-09"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-19"))+
  geom_vline(xintercept = as.POSIXct("2020-10-26"),linetype="dashed")+ #  summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-10-28"))+
  geom_vline(xintercept = as.POSIXct("2020-10-30"),linetype="dashed")+ # summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-11-02"))+
  geom_vline(xintercept = as.POSIXct("2020-11-09"))+
  geom_vline(xintercept = as.POSIXct("2020-11-10"))+
  xlim(as.POSIXct("2020-10-10"),as.POSIXct("2020-11-10"))+
  theme_classic(base_size = 15)

# Plot flowrate_mean_Lmin
ggplot(data,mapping=aes(x=datetime,y=flowrate_mean_Lmin))+
  geom_line()+
  ylim(0,21)+
  geom_vline(xintercept = as.POSIXct("2020-10-09"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-19"))+
  geom_vline(xintercept = as.POSIXct("2020-10-26"),linetype="dashed")+ #  summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-10-28"))+
  geom_vline(xintercept = as.POSIXct("2020-10-30"),linetype="dashed")+ # summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-11-02"))+
  geom_vline(xintercept = as.POSIXct("2020-11-09"))+
  geom_vline(xintercept = as.POSIXct("2020-11-10"))+
  xlim(as.POSIXct("2020-10-07"),as.POSIXct("2020-11-10"))+
  theme_classic(base_size = 15)

# Plot co2_flux
ggplot(data,mapping=aes(x=datetime,y=co2_flux))+
  geom_line()+
  geom_vline(xintercept = as.POSIXct("2020-05-11"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-18"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-25"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-03"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-04"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-12"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-06-22"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-06"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-20"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-27"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-03"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-10"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-17"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-16"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-21"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-25"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-30"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-06"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-09"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-19"))+
  geom_vline(xintercept = as.POSIXct("2020-10-26"),linetype="dashed")+ #  summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-10-28"))+
  geom_vline(xintercept = as.POSIXct("2020-10-30"),linetype="dashed")+ # summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-11-02"))+
  geom_vline(xintercept = as.POSIXct("2020-11-09"))+
  geom_vline(xintercept = as.POSIXct("2020-11-10"))+
  xlim(as.POSIXct("2020-05-02"),as.POSIXct("2020-11-10"))+
  theme_classic(base_size = 15)

# Plot ch4_flux
ggplot(data,mapping=aes(x=datetime,y=ch4_flux))+
  geom_line()+
  geom_vline(xintercept = as.POSIXct("2020-05-11"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-18"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-25"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-03"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-04"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-12"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-06-22"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-06"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-20"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-27"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-03"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-10"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-17"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-16"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-21"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-25"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-30"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-06"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-09"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-19"))+
  geom_vline(xintercept = as.POSIXct("2020-10-26"),linetype="dashed")+ #  summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-10-28"))+
  geom_vline(xintercept = as.POSIXct("2020-10-30"),linetype="dashed")+ # summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-11-02"))+
  geom_vline(xintercept = as.POSIXct("2020-11-09"))+
  geom_vline(xintercept = as.POSIXct("2020-11-10"))+
  xlim(as.POSIXct("2020-05-02"),as.POSIXct("2020-11-10"))+
  theme_classic(base_size = 15)

# Load in met data
download.file('https://github.com/FLARE-forecast/FCRE-data/raw/fcre-metstation-data/FCRmet.csv','FCRmet.csv')
metheader<-read.csv("FCRmet.csv", skip=1, as.is=T) #get header minus wonky Campbell rows
metdata<-read.csv("FCRmet.csv", skip=4, header=F) #get data minus wonky Campbell rows
names(metdata)<-names(metheader) #combine the names to deal with Campbell logger formatting

metdata$TIMESTAMP <- as.POSIXct(strptime(metdata$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz="EST"))

# Plot rain data
ggplot(metdata,mapping=aes(x=TIMESTAMP,y=Rain_mm_Tot))+
  geom_line()+
  geom_vline(xintercept = as.POSIXct("2020-05-11"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-18"))+ 
  geom_vline(xintercept = as.POSIXct("2020-05-25"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-03"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-04"),linetype="dashed")+ # no summary files
  geom_vline(xintercept = as.POSIXct("2020-06-12"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-06-22"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-06"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-20"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-07-27"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-03"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-10"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-08-17"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-16"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-21"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-25"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-09-30"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-06"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-09"))+ # summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-13"),linetype="dashed")+ # no summary files generated
  geom_vline(xintercept = as.POSIXct("2020-10-19"))+
  geom_vline(xintercept = as.POSIXct("2020-10-26"),linetype="dashed")+ #  summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-10-28"))+
  geom_vline(xintercept = as.POSIXct("2020-10-30"),linetype="dashed")+ # summary files not generated
  geom_vline(xintercept = as.POSIXct("2020-11-02"))+
  geom_vline(xintercept = as.POSIXct("2020-11-09"))+
  geom_vline(xintercept = as.POSIXct("2020-11-10"))+
  xlim(as.POSIXct("2020-10-10"),as.POSIXct("2020-11-10"))+
  theme_classic(base_size = 15)
