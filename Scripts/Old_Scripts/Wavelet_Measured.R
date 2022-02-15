### Script to conduct wavelet analysis on Eddy Flux Data
# Following Wavelets_for_fDOM_final_plots_lr_DWH.R
# A Hounshell, 16 July 2021
# UPDATED ON: 1 Nov 2021
# To test influence of gap-filled data on results
# Find a time period with minimal gap-filled data to see how gap-filling
# influences results

################################################

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,lubridate)

# Load in Eddy Flux data (cleaned using FCR_Process_BD)
eddy_flux <- read_csv("./Data/20210615_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Format for wavelet analysis: Following M. Johnson's comment, conduct wavelet on half-hourly data
# Find time periods with relatively consistent data availability
# For CO2
co2_data <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,NEE_uStar_orig) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

# Count number of NAs for each 2 day time period in the summer (May to August)
co2_data_summer <- co2_data %>% 
  filter(DateTime >= as.POSIXct("2020-05-01") & DateTime <= as.POSIXct("2020-08-05")) %>% 
  mutate(num_nas = NA)

for (i in 1:length(co2_data_summer$DateTime)){
  j = i+96
  co2_data_summer$num_nas[i] = sum(is.na(co2_data_summer$NEE_uStar_orig[i:j]))
}

# Highest amount of data in a 2-day period: 2020-06-28 23:30:00
co2_data_sub <- co2_data %>% 
  filter(DateTime >= as.POSIXct("2020-06-28 20:00:00") & DateTime < as.POSIXct("2020-07-01 03:00:00"))

(96-sum(is.na(co2_data_sub$NEE_uStar_orig)))/96*100 # Results in 61% data coverage

# Linearly interpolate between missing time points
co2_data_sub <- co2_data_sub %>% 
  mutate(NEE_uStar_orig = na.approx(NEE_uStar_orig)) %>% 
  select(DateTime,NEE_uStar_f,NEE_uStar_orig) %>% 
  filter(DateTime >= as.POSIXct("2020-06-28 20:00:00") & DateTime < as.POSIXct("2020-06-30 20:00:00"))

# Separate into NEE_f and NEE_orig
co2_f <- co2_data_sub %>% 
  select(DateTime,NEE_uStar_f) %>% 
  mutate(Time = 1:length(co2_data_sub$DateTime))

co2_orig <- co2_data_sub %>% 
  select(DateTime,NEE_uStar_orig) %>% 
  mutate(Time = 1:length(co2_data_sub$DateTime))

# Format data for wavelet: start with gap-filled data (_f)
headers<-names(co2_f)
all<-headers[2]
temp<-matrix(-99,length(co2_f$Time),length(all),byrow=F)
head(temp)

temp <- scale(co2_f$NEE_uStar_f)

snt_co2_f<-as.data.frame(temp)
snt_co2_f<-setNames(snt_co2_f, all)#new dataframe with standard normal transformed data
head(snt_co2_f)

# Format data for wavelet: measured, interpolated data
headers<-names(co2_orig)
all<-headers[2]
temp<-matrix(-99,length(co2_orig$Time),length(all),byrow=F)
head(temp)

temp <- scale(co2_orig$NEE_uStar_orig)

snt_co2_orig<-as.data.frame(temp)
snt_co2_orig<-setNames(snt_co2_orig, all)#new dataframe with standard normal transformed data
head(snt_co2_orig)

#### Getting Wavelet co2_output and plotting ####
Time<-co2_f$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_f<-cbind(Time, snt_co2_f)
head(snt_co2_f)
snt_co2_f <- snt_co2_f %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_f))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_f)

co2_output_f<-morlet(snt_co2_f$NEE_uStar_f, snt_co2_f$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

co2_output_f$period <- round((co2_output_f$period/48), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
co2_dummy_f <- c(1:157)
for (j in 1:ncol(co2_output_f$Power)) {
  co2_dummy_f[j] <- mean(co2_output_f$Power[,j])
}

co2_powerplot_f <- as.data.frame(cbind(co2_dummy_f, co2_output_f$period))
co2_powerplot_f <- rename(co2_powerplot_f, c( "mean_power" = "co2_dummy_f"))
co2_powerplot_f <- rename(co2_powerplot_f, c( "period" = "V2"))

# Month mean power plot 
co2_powerplot_f <- co2_powerplot_f %>% 
  filter(period < 182)

#######
# Do the same for original Co2 data
Time<-co2_orig$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_orig<-cbind(Time, snt_co2_orig)
head(snt_co2_orig)
snt_co2_orig <- snt_co2_orig %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_orig))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_orig)

co2_output_orig<-morlet(snt_co2_orig$NEE_uStar_orig, snt_co2_orig$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

co2_output_orig$period <- round((co2_output_orig$period/48), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
co2_dummy_orig <- c(1:157)
for (j in 1:ncol(co2_output_orig$Power)) {
  co2_dummy_orig[j] <- mean(co2_output_orig$Power[,j])
}

co2_powerplot_orig <- as.data.frame(cbind(co2_dummy_orig, co2_output_orig$period))
co2_powerplot_orig <- rename(co2_powerplot_orig, c( "mean_power" = "co2_dummy_orig"))
co2_powerplot_orig <- rename(co2_powerplot_orig, c( "period" = "V2"))

# Month mean power plot 
co2_powerplot_orig <- co2_powerplot_orig %>% 
  filter(period < 182)

###################
# Methane!
# Find time periods with relatively consistent data availability
# For CO2
ch4_data <- eddy_flux %>% 
  select(DateTime,ch4_flux_uStar_f,ch4_flux_uStar_orig) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

# Count number of NAs for each 2 day time period in the summer (May to August)
ch4_data_summer <- ch4_data %>% 
  filter(DateTime >= as.POSIXct("2020-05-01") & DateTime <= as.POSIXct("2020-08-05")) %>% 
  mutate(num_nas = NA)

for (i in 1:length(ch4_data_summer$DateTime)){
  j = i+96
  ch4_data_summer$num_nas[i] = sum(is.na(ch4_data_summer$ch4_flux_uStar_orig[i:j]))
}

# Highest amount of data in a 2-day period: 2020-06-29 05:00:00
ch4_data_sub <- ch4_data %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 01:00:00") & DateTime < as.POSIXct("2020-07-01 01:00:00"))

(96-sum(is.na(ch4_data_sub$ch4_flux_uStar_orig)))/96*100 # Results in 51% data coverage

# Linearly interpolate between missing time points
ch4_data_sub <- ch4_data_sub %>% 
  mutate(ch4_flux_uStar_orig = na.approx(ch4_flux_uStar_orig)) %>% 
  select(DateTime,ch4_flux_uStar_f,ch4_flux_uStar_orig)

# Separate into NEE_f and NEE_orig
ch4_f <- ch4_data_sub %>% 
  select(DateTime,ch4_flux_uStar_f) %>% 
  mutate(Time = 1:length(ch4_data_sub$DateTime))

ch4_orig <- ch4_data_sub %>% 
  select(DateTime,ch4_flux_uStar_orig) %>% 
  mutate(Time = 1:length(ch4_data_sub$DateTime))

# Format data for wavelet: start with gap-filled data (_f)
headers<-names(ch4_f)
all<-headers[2]
temp<-matrix(-99,length(ch4_f$Time),length(all),byrow=F)
head(temp)

temp <- scale(ch4_f$ch4_flux_uStar_f)

snt_ch4_f<-as.data.frame(temp)
snt_ch4_f<-setNames(snt_ch4_f, all)#new dataframe with standard normal transformed data
head(snt_ch4_f)

# Format data for wavelet: measured, interpolated data
headers<-names(ch4_orig)
all<-headers[2]
temp<-matrix(-99,length(ch4_orig$Time),length(all),byrow=F)
head(temp)

temp <- scale(ch4_orig$ch4_flux_uStar_orig)

snt_ch4_orig<-as.data.frame(temp)
snt_ch4_orig<-setNames(snt_ch4_orig, all)#new dataframe with standard normal transformed data
head(snt_ch4_orig)

#### Getting Wavelet co2_output and plotting ####
Time<-ch4_f$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_f<-cbind(Time, snt_ch4_f)
head(snt_ch4_f)
snt_ch4_f <- snt_ch4_f %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_f))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_f)

ch4_output_f<-morlet(snt_ch4_f$ch4_flux_uStar_f, snt_ch4_f$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

ch4_output_f$period <- round((ch4_output_f$period/48), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
ch4_dummy_f <- c(1:157)
for (j in 1:ncol(ch4_output_f$Power)) {
  ch4_dummy_f[j] <- mean(ch4_output_f$Power[,j])
}

ch4_powerplot_f <- as.data.frame(cbind(ch4_dummy_f, ch4_output_f$period))
ch4_powerplot_f <- rename(ch4_powerplot_f, c( "mean_power" = "ch4_dummy_f"))
ch4_powerplot_f <- rename(ch4_powerplot_f, c( "period" = "V2"))

# Month mean power plot 
ch4_powerplot_f <- ch4_powerplot_f %>% 
  filter(period < 182)

#######
# Do the same for original ch4 data
Time<-ch4_orig$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_orig<-cbind(Time, snt_ch4_orig)
head(snt_ch4_orig)
snt_ch4_orig <- snt_ch4_orig %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_orig))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_orig)

ch4_output_orig<-morlet(snt_ch4_orig$ch4_flux_uStar_orig, snt_ch4_orig$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

ch4_output_orig$period <- round((ch4_output_orig$period/48), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
ch4_dummy_orig <- c(1:157)
for (j in 1:ncol(ch4_output_orig$Power)) {
  ch4_dummy_orig[j] <- mean(ch4_output_orig$Power[,j])
}

ch4_powerplot_orig <- as.data.frame(cbind(ch4_dummy_orig, ch4_output_orig$period))
ch4_powerplot_orig <- rename(ch4_powerplot_orig, c( "mean_power" = "ch4_dummy_orig"))
ch4_powerplot_orig <- rename(ch4_powerplot_orig, c( "period" = "V2"))

# Month mean power plot 
ch4_powerplot_orig <- ch4_powerplot_orig %>% 
  filter(period < 182)

## Make Figures/look at results!
# Plot to see differences
co2_interp_comp <- co2_data_sub %>% 
  ggplot()+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(mapping=aes(x=DateTime,y=NEE_uStar_f,color="Gap-filled"))+
  geom_point(mapping=aes(x=DateTime,y=NEE_uStar_f,color="Gap-filled"))+
  geom_line(mapping=aes(x=DateTime,y=NEE_uStar_orig,color="Measured"))+
  geom_point(mapping=aes(x=DateTime,y=NEE_uStar_orig,color="Measured"))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Global power spectra
co2_power <- ggplot()+
  geom_line(data = co2_powerplot_f, mapping=aes(x=period,y=mean_power,color="Gap-filled"))+
  geom_point(data = co2_powerplot_f, mapping=aes(x=period,y=mean_power,color="Gap-filled"))+
  geom_line(data = co2_powerplot_orig, mapping=aes(x=period,y=mean_power,color="Measured"))+
  geom_point(data = co2_powerplot_orig, mapping=aes(x=period,y=mean_power,color="Measured"))+
  labs(x = "Period (Days)",
       y = (expression(paste("Mean ",Power^2,))))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  xlim(0,1)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())


ch4_interp_comp <- ch4_data_sub %>% 
  ggplot()+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(mapping=aes(x=DateTime,y=ch4_flux_uStar_f,color="Gap-filled"))+
  geom_point(mapping=aes(x=DateTime,y=ch4_flux_uStar_f,color="Gap-filled"))+
  geom_line(mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="Measured"))+
  geom_point(mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="Measured"))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Global power spectra
ch4_power <- ggplot()+
  geom_line(data = ch4_powerplot_f, mapping=aes(x=period,y=mean_power,color="Gap-filled"))+
  geom_point(data = ch4_powerplot_f, mapping=aes(x=period,y=mean_power,color="Gap-filled"))+
  geom_line(data = ch4_powerplot_orig, mapping=aes(x=period,y=mean_power,color="Measured"))+
  geom_point(data = ch4_powerplot_orig, mapping=aes(x=period,y=mean_power,color="Measured"))+
  labs(x = "Period (Days)",
       y = (expression(paste("Mean ",Power^2,))))+
  scale_color_manual(breaks=c("Gap-filled","Measured"),
                     values=c("darkgrey","black"))+
  xlim(0,1)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(co2_interp_comp,co2_power,ch4_interp_comp,ch4_power,nrow=2,ncol=2,
          labels=c("A.","B.","C.","D."),font.label = list(face="plain",size=15),
          common.legend = TRUE)

ggsave("./Fig_Output/Wavelet_comps.jpg",width = 10, height = 9, units="in",dpi=320)


##### Look at 24 hour diel patterns #####
# Use the same two week time period as the wavelet analysis
# Separate into day (0900 to 1500) and night (2100 to 0300)
daynight <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,NEE_uStar_orig,ch4_flux_uStar_f,ch4_flux_uStar_orig) %>% 
  filter(DateTime >= as.POSIXct("2020-06-28 20:00:00") & DateTime < as.POSIXct("2020-07-01 03:00:00")) %>% 
  mutate(Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 9 & Hour <= 14, "Day",
                       ifelse(Hour >= 21 | Hour <= 2, "Night", NA))) %>% 
  filter(diel == "Day" | diel == "Night")

# Conduct statistical analysis
kruskal.test(NEE_uStar_orig ~ diel, data=daynight)
kruskal.test(NEE_uStar_f ~ diel, data=daynight)

kruskal.test(ch4_flux_uStar_orig ~ diel, data=daynight)
kruskal.test(ch4_flux_uStar_f ~ diel, data=daynight)

# Aggregate data for plotting
daynight_long <- daynight %>% 
  pivot_longer(cols = -c(DateTime,Hour,diel), names_to = "flux", values_to = "concentration")

# Plot?!
co2_diel <- daynight_long %>% 
  filter(flux == "NEE_uStar_f"|flux == "NEE_uStar_orig") %>% 
  ggplot(daynight_long,mapping=aes(x=flux,y=concentration,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  scale_x_discrete(breaks=c("NEE_uStar_f","NEE_uStar_orig"),
                   labels=c("Gap-filled", "Measured"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_diel <- daynight_long %>% 
  filter(flux == "ch4_flux_uStar_f"|flux == "ch4_flux_uStar_orig") %>% 
  ggplot(daynight_long,mapping=aes(x=flux,y=concentration,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  scale_x_discrete(breaks=c("ch4_flux_uStar_f","ch4_flux_uStar_orig"),
                   labels=c("Gap-filled", "Measured"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(co2_diel,ch4_diel,nrow=1,ncol=2,common.legend = TRUE,
          labels=c("A.","B."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/2day_diel.jpg",width = 10, height=5, units="in",dpi=320)