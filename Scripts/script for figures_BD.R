### Code originally from B. D'Acunha for plotting EC data
# Adapted by A. Hounshell, 08 July 2021

# Load in libraries
pacman::p_load(readr,dplyr,lubridate,zoo,ggplot2,ggpubr,openair,REddyProc)

# Set working directory
wd <- getwd()
setwd(wd)

# read data: Processed using EReddy Proc
fcr_gf <- read_csv("./Data/20210615_EC_processed.csv")

# set DateTime to EST
fcr_gf$DateTime <- as_datetime(fcr_gf$DateTime, tz = 'EST')

# if solar radiation more than 12, it's day, otherwise, it's night
fcr_gf$Day_time <- ifelse(fcr_gf$Rg_f >= 12, 1, 0)

#############################################################
#COUNTING NA
nas_percent <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
																 NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
						ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)


nas_count <- fcr_gf %>% select(DateTime, NEE, ch4_flux,
															 NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	summarise(co2_available = n() - sum(is.na(NEE_uStar_orig)),
						ch4_available = n()-sum(is.na(ch4_flux_uStar_orig)))

fcr_gf %>% select(DateTime, NEE, ch4_flux, 
									NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	group_by(year(DateTime), month(DateTime)) %>% 
	summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
						ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)
#############################################################

#############################################################
# Check where fluxes are coming from: using NON gap-filled data
ggplot(fcr_gf,mapping=aes(x=wind_dir,y=NEE_uStar_orig,color=Ustar))+
  geom_point()

ggplot(fcr_gf,mapping=aes(x=wind_dir,y=ch4_flux_uStar_orig,color=Ustar))+
  geom_point()

fcr_gf %>% ggplot(aes(wind_dir, NEE_uStar_orig,color=Ustar)) + 
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

############# Thinking about time periods were we have 'good' coverage of actual measured fluxes ########
# For diel comparisons in summer (stratified period) vs. winter
# Summer stratified period: June 29 to July 4, 2020
# Day: 5 am - 7 pm
ggplot(fcr_gf)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(mapping=aes(x=DateTime,y=NEE_uStar_orig),size=2)+
  geom_line(mapping=aes(x=DateTime,y=NEE_uStar_f),size=0.7)+
  xlim(as.POSIXct("2020-06-29 05:00:00"),as.POSIXct("2020-07-06 04:30:00"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-10,40)+
  theme_classic(base_size = 15)

ggplot(fcr_gf)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 04:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),size=2)+
  geom_line(mapping=aes(x=DateTime,y=ch4_flux_uStar_f),size=0.7)+
  xlim(as.POSIXct("2020-06-29 05:00:00"),as.POSIXct("2020-07-06 04:30:00"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)

# Calculate number of nas during this time period
summer_day <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                 NEE_uStar_orig, ch4_flux_uStar_orig,daytime) %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00") & DateTime <= as.POSIXct("2020-07-06 07:00:00")) %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour >= 7 & hour <= 18) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

summer_night <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                NEE_uStar_orig, ch4_flux_uStar_orig,daytime) %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00") & DateTime <= as.POSIXct("2020-07-06 06:30:00")) %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour < 7 | hour > 18) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

summer_total <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                  NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00") & DateTime <= as.POSIXct("2020-07-06 06:30:00")) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

## And then look at winter data!
# Day time: 7am to 4:30 pm
ggplot(fcr_gf)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-04 15:00:00"),xmax = as.POSIXct("2021-01-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-05 15:00:00"),xmax = as.POSIXct("2021-01-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-06 15:00:00"),xmax = as.POSIXct("2021-01-07 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-07 15:00:00"),xmax = as.POSIXct("2021-01-08 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-08 15:00:00"),xmax = as.POSIXct("2021-01-09 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-09 15:00:00"),xmax = as.POSIXct("2021-01-10 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-10 15:00:00"),xmax = as.POSIXct("2021-01-11 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(mapping=aes(x=DateTime,y=NEE_uStar_orig,color="NEE"))+
  geom_point(mapping=aes(x=DateTime,y=ch4_flux_uStar_orig*1000,color="ch4"))+
  xlim(as.POSIXct("2021-01-04 07:00:00"),as.POSIXct("2021-01-11 06:30:00"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)

# Calculate NA's
winter_day <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                NEE_uStar_orig, ch4_flux_uStar_orig,daytime) %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00") & DateTime <= as.POSIXct("2021-01-11 06:30:00")) %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour >= 7 & hour <= 18) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

winter_night <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                  NEE_uStar_orig, ch4_flux_uStar_orig,daytime) %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00") & DateTime <= as.POSIXct("2021-01-11 06:30:00")) %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour < 7 | hour > 18) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

winter_total <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
                                  NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00") & DateTime <= as.POSIXct("2021-01-11 06:30:00")) %>% 
  summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)

# plotting 30min values
#air temperature

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, Tair_f)) + 
	geom_line(aes(DateTime, Tair_f)) +
	theme_bw()

# wind speed

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, u)) + 
	geom_line(aes(DateTime, u)) +
	theme_bw()

# NEE

fcr_gf %>% ggplot() + 
	geom_point(aes(DateTime, NEE_uStar_f), color = 'red') +
	geom_line(aes(DateTime, NEE_uStar_orig), color = 'gray1') + 
	theme_bw()


# CH4

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, ch4_flux_uStar_f), col = 'red') + 
	geom_line(aes(DateTime, ch4_flux_uStar_orig), color = 'gray1') +
	theme_bw()


# air pressure

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, airP)) + 
	geom_line(aes(DateTime, airP)) +
	theme_bw()


#####
# Group by hour, plot hourly data
fcr_hourly <- fcr_gf %>% 
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

# Plot hourly data for diel comparisons
summer_hourly <- fcr_hourly %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 06:00:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-29 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2020-06-30 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-30 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2020-07-01 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-01 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2020-07-02 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-02 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2020-07-03 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-03 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2020-07-04 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-04 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2020-07-05 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2020-07-05 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2020-07-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

colnames(summer_hourly) <- paste(colnames(summer_hourly), "summer", sep="_")

summer_hourly <- summer_hourly %>% 
  rename(day = day_summer, hour = hour_summer)

winter_hourly <- fcr_hourly %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 06:00:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-04 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2021-01-05 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-05 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2021-01-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-06 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2021-01-07 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-07 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2021-01-08 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-08 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2021-01-09 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-09 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2021-01-10 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2021-01-10 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2021-01-11 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

colnames(winter_hourly) <- paste(colnames(winter_hourly), "winter", sep="_")

winter_hourly <- winter_hourly %>% 
  rename(day = day_winter, hour = hour_winter)

diel_hourly <- left_join(summer_hourly,winter_hourly,by=c("day","hour"))

summer_raw <- fcr_gf %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,NEE_uStar_f,ch4_flux_uStar_f) %>%   
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 06:30:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-29 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2020-06-30 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-30 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2020-07-01 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-01 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2020-07-02 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-02 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2020-07-03 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-03 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2020-07-04 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-04 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2020-07-05 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2020-07-05 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2020-07-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

colnames(summer_raw) <- paste(colnames(summer_raw), "summer", sep="_")

summer_raw <- summer_raw %>% 
  rename(day = day_summer, hour = hour_summer)
  
winter_raw <- fcr_gf %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,NEE_uStar_f,ch4_flux_uStar_f) %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 06:30:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-04 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2021-01-05 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-05 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2021-01-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-06 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2021-01-07 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-07 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2021-01-08 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-08 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2021-01-09 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-09 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2021-01-10 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2021-01-10 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2021-01-11 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

colnames(winter_raw) <- paste(colnames(winter_raw), "winter", sep="_")

winter_raw <- winter_raw %>% 
  rename(day = day_winter, hour = hour_winter)

diel_raw <- left_join(summer_raw,winter_raw,by=c("day","hour"))

# Aggregate data to a 24 hour period over summer and winter
diel_agg <- diel_raw %>% 
  group_by(hour) %>% 
  summarise(NEE_summer = mean(NEE_uStar_f_summer, na.rm = TRUE),
            NEE_sd_summer = sd(NEE_uStar_f_summer, na.rm = TRUE),
            CH4_summer = mean(ch4_flux_uStar_f_summer, na.rm = TRUE),
            CH4_sd_summer = sd(ch4_flux_uStar_f_summer, na.rm = TRUE),
            NEE_winter = mean(NEE_uStar_f_winter, na.rm = TRUE),
            NEE_sd_winter = sd(NEE_uStar_f_winter, na.rm = TRUE),
            CH4_winter = mean(ch4_flux_uStar_f_winter, na.rm = TRUE),
            CH4_sd_winter = sd(ch4_flux_uStar_f_winter, na.rm = TRUE))
    

# Plot together - winter and summer periods
labs <- c(1,3,5,7)

co2_week <- ggplot(diel_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(diel_raw,mapping=aes(x=DateTime_summer,y=NEE_uStar_orig_summer,color="summer"),alpha=0.2)+
  geom_line(mapping=aes(x=DateTime_summer,y=NEE_summer,color="summer"),size=1)+
  geom_ribbon(mapping=aes(x=DateTime_summer,y=NEE_summer,ymin=NEE_summer-NEE_sd_summer,ymax=NEE_summer+NEE_sd_summer),fill="#E63946",alpha=0.3)+
  geom_point(diel_raw,mapping=aes(x=DateTime_summer,y=NEE_uStar_orig_winter,color="winter"),alpha=0.2)+
  geom_line(mapping=aes(x=DateTime_summer,y=NEE_winter,color="winter"),size=1)+
  geom_ribbon(mapping=aes(x=DateTime_summer,y=NEE_winter,ymin=NEE_winter-NEE_sd_winter,ymax=NEE_winter+NEE_sd_winter),fill="#4c8bfe",alpha=0.3)+
  xlim(as.POSIXct("2020-06-29 07:00:00"),as.POSIXct("2020-07-06 06:30:00"))+
  scale_x_discrete(labels=labs)+
  scale_color_manual(breaks=c("summer","winter"), 
                     labels=c("Summer","Winter"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  labs(color="")

# Plot aggregated daily CO2 for winter and summer
co2_day <- ggplot(diel_agg)+
  annotate(geom="rect",xmin = 0,xmax = 6.5,ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin=19,xmax = 24,ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(mapping=aes(x=hour,y=NEE_summer,color="summer"),size=1)+
  geom_point(mapping=aes(x=hour,y=NEE_summer,color="summer"))+
  geom_ribbon(mapping=aes(x=hour,y=NEE_summer,ymin=NEE_summer-NEE_sd_summer,ymax=NEE_summer+NEE_sd_summer),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=hour,y=NEE_winter,color="winter"),size=1)+
  geom_point(mapping=aes(x=hour,y=NEE_winter,color="winter"))+
  geom_ribbon(mapping=aes(x=hour,y=NEE_winter,ymin=NEE_winter-NEE_sd_winter,ymax=NEE_winter+NEE_sd_winter),fill="#4c8bfe",alpha=0.3)+
  scale_color_manual(breaks=c("summer","winter"), 
                     labels=c("Summer","Winter"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("Hour")+
  theme_classic(base_size = 15)+
  labs(color="")

# Plot winter and summer for CH4
ch4_week <- ggplot(diel_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(diel_raw,mapping=aes(x=DateTime_summer,y=ch4_flux_uStar_orig_summer,color="summer"),alpha=0.2)+
  geom_line(mapping=aes(x=DateTime_summer,y=CH4_summer,color="summer"),size=1)+
  geom_ribbon(mapping=aes(x=DateTime_summer,y=CH4_summer,ymin=CH4_summer-CH4_sd_summer,ymax=CH4_summer+CH4_sd_summer),fill="#E63946",alpha=0.3)+
  geom_point(diel_raw,mapping=aes(x=DateTime_summer,y=ch4_flux_uStar_orig_winter,color="winter"),alpha=0.2)+
  geom_line(mapping=aes(x=DateTime_summer,y=CH4_winter,color="winter"),size=1)+
  geom_ribbon(mapping=aes(x=DateTime_summer,y=CH4_winter,ymin=CH4_winter-CH4_sd_winter,ymax=CH4_winter+CH4_sd_winter),fill="#4c8bfe",alpha=0.3)+
  xlim(as.POSIXct("2020-06-29 07:00:00"),as.POSIXct("2020-07-06 06:30:00"))+
  scale_color_manual(breaks=c("summer","winter"), 
                     labels=c("Summer","Winter"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic(base_size = 15)+
  labs(color="")

# Plot aggregated daily CH4 for winter and summer
ch4_day <- ggplot(diel_agg)+
  annotate(geom="rect",xmin = 0,xmax = 6.5,ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin=19,xmax = 24,ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(mapping=aes(x=hour,y=CH4_summer,color="summer"),size=1)+
  geom_point(mapping=aes(x=hour,y=CH4_summer,color="summer"))+
  geom_ribbon(mapping=aes(x=hour,y=CH4_summer,ymin=CH4_summer-CH4_sd_summer,ymax=CH4_summer+CH4_sd_summer),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=hour,y=CH4_winter,color="winter"),size=1)+
  geom_point(mapping=aes(x=hour,y=CH4_winter,color="winter"))+
  geom_ribbon(mapping=aes(x=hour,y=CH4_winter,ymin=CH4_winter-CH4_sd_winter,ymax=CH4_winter+CH4_sd_winter),fill="#4c8bfe",alpha=0.3)+
  scale_color_manual(breaks=c("summer","winter"), 
                     labels=c("Summer","Winter"),
                     values=c("#E63946","#4c8bfe"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("Hour")+
  theme_classic(base_size = 15)+
  labs(color="")

ggarrange(co2_week,co2_day,ch4_week,ch4_day,nrow=2,ncol=2,common.legend=TRUE,labels = c("A.", "B.","C.","D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fluxes_Diel.jpg",width = 8, height=7, units="in",dpi=320)

# Old code for diel variation ----

# Summer stratified period: June 29 to July 4, 2020
# Day: 5 am - 6:30 pm
nee_summer <- ggplot(fcr_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.2)+
  geom_ribbon(mapping=aes(x=DateTime,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=DateTime,y=NEE),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-06-29 07:00:00"),as.POSIXct("2020-07-06 06:30:00"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-30,50)+
  theme_classic(base_size = 15)

ch4_summer <- ggplot(fcr_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-29 19:00:00"),xmax = as.POSIXct("2020-06-30 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-06-30 19:00:00"),xmax = as.POSIXct("2020-07-01 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-01 19:00:00"),xmax = as.POSIXct("2020-07-02 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-02 19:00:00"),xmax = as.POSIXct("2020-07-03 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-03 19:00:00"),xmax = as.POSIXct("2020-07-04 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-04 19:00:00"),xmax = as.POSIXct("2020-07-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2020-07-05 19:00:00"),xmax = as.POSIXct("2020-07-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha=0.2)+
  geom_ribbon(mapping=aes(x=DateTime,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=DateTime,y=CH4),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-06-29 07:00:00"),as.POSIXct("2020-07-06 06:30:00"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-0.07,0.07)+
  theme_classic(base_size = 15)

# Winter period
nee_winter <- ggplot(fcr_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-04 19:00:00"),xmax = as.POSIXct("2021-01-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-05 19:00:00"),xmax = as.POSIXct("2021-01-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-06 19:00:00"),xmax = as.POSIXct("2021-01-07 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-07 19:00:00"),xmax = as.POSIXct("2021-01-08 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-08 19:00:00"),xmax = as.POSIXct("2021-01-09 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-09 19:00:00"),xmax = as.POSIXct("2021-01-10 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-10 19:00:00"),xmax = as.POSIXct("2021-01-11 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.2)+
  geom_ribbon(mapping=aes(x=DateTime,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=DateTime,y=NEE),color="#E63946",size = 1)+
  xlim(as.POSIXct("2021-01-04 07:00:00"),as.POSIXct("2021-01-11 06:30:00"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-3,8)+
  theme_classic(base_size = 15)

ch4_winter <- ggplot(fcr_hourly)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-04 19:00:00"),xmax = as.POSIXct("2021-01-05 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-05 19:00:00"),xmax = as.POSIXct("2021-01-06 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-06 19:00:00"),xmax = as.POSIXct("2021-01-07 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-07 19:00:00"),xmax = as.POSIXct("2021-01-08 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-08 19:00:00"),xmax = as.POSIXct("2021-01-09 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-09 19:00:00"),xmax = as.POSIXct("2021-01-10 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  annotate(geom="rect",xmin = as.POSIXct("2021-01-10 19:00:00"),xmax = as.POSIXct("2021-01-11 06:30:00"),ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha=0.2)+
  geom_ribbon(mapping=aes(x=DateTime,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=DateTime,y=CH4),color="#E63946",size = 1)+
  xlim(as.POSIXct("2021-01-04 07:00:00"),as.POSIXct("2021-01-11 06:30:00"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-0.02,0.02)+
  theme_classic(base_size = 15)

ggarrange(nee_summer,ch4_summer, nee_winter, ch4_winter, nrow = 2, ncol = 2, align = "v")

ggsave("./Fig_Output/HourlyFluxes_Season_Avg.jpg",width = 8, height=6, units="in",dpi=320)

#####

# plot co2 hourly data

fcr_gf  %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																	"Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																	'Dec', 'Jan', 'Feb', 'Mar')), 
										 hour = hour(DateTime)) %>% 
	summarise(NEE_m = mean(NEE_uStar_f, na.rm = TRUE),
						NEEsd = sd(NEE_uStar_f, na.rm = TRUE)) %>% 
	ggplot(aes(hour, NEE_m, col = factor(year))) + 
	geom_point(size = 3) + geom_line(size = 1) +
	facet_wrap(~month) + ylab(expression(CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(" ") +
	geom_errorbar(aes(ymin = NEE_m - NEEsd, ymax = NEE_m + NEEsd), 
								width = 0.5, size = 0.8) +
	geom_hline(yintercept = 0, lty = 2, col = 'red') +
	theme_bw() + labs(col='Year') + scale_color_brewer(palette = "Dark2")

# plot ch4 hourly data

fcr_gf %>% 
	dplyr::group_by(year = year(DateTime), 
									month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																				"Jul", "Aug", "Sep", "Oct", 'Nov', 
																																				'Dec', 'Jan', 'Feb', 'Mar')), 
									hour = hour(DateTime)) %>% 
	dplyr::summarise(ch4_m = mean(ch4_flux_uStar_f, na.rm = TRUE), 
									 ch4_sd = sd(ch4_flux_uStar_f)) %>% ungroup() %>% 
	ggplot(aes(hour, ch4_m, col = factor(year))) + 
	geom_point(size = 3) + geom_line(size = 1) +
	facet_wrap(~month) + ylab(expression(CH[4]~(mu~mol~m^-2~s^-1))) +
	xlab(" ") +
	geom_errorbar(aes(ymin = ch4_m - ch4_sd, ymax = ch4_m + ch4_sd), 
								width = 0.5, size = 0.8) +
	geom_hline(yintercept = 0, lty = 2, col = 'red') +
	theme_bw() + labs(col='Year') + scale_color_brewer(palette = "Dark2")

# plot radiation hourly

fcr_gf %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																 "Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																 'Dec', 'Jan', 'Feb', 'Mar')), 
										hour = hour(DateTime)) %>% 
	summarise(Rg = mean(Rg, na.rm = TRUE)) %>% 
	ggplot(aes(hour, Rg, col = factor(year))) + geom_point() + 
	facet_wrap(~month) + theme_bw() + ylab('Incoming radiation') + xlab("") +
	scale_color_brewer(palette = "Dark2")


# air temperature 

fcr_gf %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																 "Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																 'Dec', 'Jan', 'Feb', 'Mar')), 
										hour = hour(DateTime)) %>% 
	summarise(airt = mean(Tair_f, na.rm = TRUE)) %>% 
	ggplot(aes(hour, airt, col = factor(year))) + geom_point() + 
	facet_wrap(~month) + theme_bw() + ylab('Air temperature') + xlab("") +
	scale_color_brewer(palette = "Dark2")

# group by hour 

fcr_hourly$datetime <- with(fcr_hourly, ymd_h(paste(Year, Month, Day, Hour, sep= ' ')))

# plot data averaged hourly

fcr_hourly %>%
	ggplot() +
	geom_point(aes(datetime, CH4)) + 
	geom_line(aes(datetime, CH4)) +
	ylab(expression(CH[4]~(mu~mol~m^-2~s^-1))) +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), lty =2) +
	scale_x_datetime(date_breaks = '2 week', date_labels = '%b-%d') +
	xlab(" ") +
	theme_bw()


fcr_hourly %>%
	ggplot() +
	geom_point(aes(datetime, NEE)) + 
	geom_line(aes(datetime, NEE)) +
	ylab(expression(CO[2]~(mu~mol~m^-2~s^-1))) +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), lty =2) +
	scale_x_datetime(date_breaks = '2 weeks', date_labels = '%b-%d') +
	xlab(" ") +
	theme_bw()

#####################################################################################

# Group by day -> plot

# data in umolm2s
fcr_daily <- fcr_gf %>% 
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


# nee value in gCm2day, CH4 values = mgCm2day

fcr_dailysum <- fcr_hourly %>% 
  mutate(Year = year(DateTime), 
           Month = month(DateTime), 
           Day = day(DateTime), 
           Hour = hour(DateTime)) %>% dplyr::group_by(Year, Month, Day) %>% 
	dplyr::summarise(
		NEE_sum = sum(NEE*3600*12/1000000, na.rm = TRUE),
		NEE_sum05 = sum(NEE05*3600*12/1000000, na.rm = TRUE),
		NEE_sum50 = sum(NEE50*3600*12/1000000, na.rm = TRUE),
		NEE_sum95 = sum(NEE95*3600*12/1000000, na.rm = TRUE),
		CH4_sum = sum(CH4*3600*12/1000, na.rm = TRUE),
		CH4_sum05 = sum(CH405*3600*12/1000, na.rm = TRUE),
		CH4_sum50 = sum(CH450*3600*12/1000, na.rm = TRUE),
		CH4_sum95 = sum(CH495*3600*12/1000, na.rm = TRUE),
		Tmin = min(Tmin, na.rm = TRUE),
		Tmax = max(Tmax, na.rm = TRUE),
		Tmean = mean(Tmean, na.rm = TRUE),
		H = sum(H, na.rm = TRUE),
		LE = sum(LE, na.rm = TRUE),
		pressure = mean(pressure, na.rm = TRUE),
		VPD = mean(VPD, na.rm = TRUE),
		RH = mean(RH, na.rm = TRUE),
		umean = mean(umean, na.rm = TRUE),
		umax = max(umax, na.rm = TRUE),
		PAR = sum(PAR_tot, na.rm = TRUE),
		PPT_tot = sum(precip_sum, na.rm = TRUE),
		LW_in = mean(LW_in, na.rm = TRUE),
		LW_out = mean(LW_out, na.rm = TRUE),
		Rg = mean(Rg, na.rm = TRUE),
		SW_out = mean(SW_out, na.rm = TRUE),
		Rn = mean(Rn , na.rm = TRUE))


fcr_dailysum$Date <- as.POSIXct(paste(fcr_dailysum$Year, fcr_dailysum$Month, fcr_dailysum$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# daily means
dailynee <- fcr_daily %>% 
	ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE_mean,ymin=NEE_mean-NEE_sd,ymax=NEE_mean+NEE_sd),fill="#E63946",alpha=0.5)+
	geom_line(aes(Date, NEE_mean),color="#E63946",size = 1) +
	xlab("") +
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
	theme_classic(base_size = 15)


dailych4 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(fcr_gf,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4_mean,ymin=CH4_mean-CH4_sd,ymax=CH4_mean+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4_mean),color="#E63946",size = 1) +
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-04-05"),as.POSIXct("2021-05-05"))+
  theme_classic(base_size = 15)


dailynee
dailych4

ggarrange(dailynee, dailych4, nrow = 2, ncol = 1, align = "v")

ggsave("./Fig_Output/DailyFluxes_Avg.jpg",width = 8, height=8, units="in",dpi=320)

### Think about winter variability (especially with ice!)
ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))
ice_on <- ice %>% 
  filter(Date>"2020-01-01" & IceOn == 1)
ice_off <- ice %>% 
  filter(Date>"2020-01-01" & IceOff == 1)

# Calculate average flux for each ice on/off period
ice_off_1 <- fcr_gf %>% 
  filter(DateTime >= "2020-12-19 19:00:00" & DateTime < "2020-12-26 19:00:00") %>% 
  mutate(ice_period = "1") %>% 
  mutate(ice = "off")

ice_on_1 <- fcr_gf %>% 
  filter(DateTime >= "2020-12-26 19:00:00" & DateTime < "2020-12-29 19:00:00")%>% 
  mutate(ice_period = "2") %>% 
  mutate(ice = "on")

ice_off_2 <- fcr_gf %>% 
  filter(DateTime >= "2020-12-29 19:00:00" & DateTime < "2021-01-09 19:00:00")%>% 
  mutate(ice_period = "3") %>% 
  mutate(ice = "off")

ice_on_2 <- fcr_gf %>% 
  filter(DateTime >= "2021-01-09 19:00:00" & DateTime < "2021-02-09 19:00:00")%>% 
  mutate(ice_period = "4") %>% 
  mutate(ice = "on")

ice_all <- rbind(ice_off_1,ice_on_1,ice_off_2,ice_on_2)

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
  geom_point(fcr_gf, mapping = aes(DateTime, NEE_uStar_orig),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping = aes(x = Date, y = NEE, ymin = NEE-NEE_sd,ymax = NEE+NEE_sd),fill="#E63946",alpha=0.4)+
  geom_line(mapping = aes(Date, NEE),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-02-11"))+
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-10,10)+
  theme_classic(base_size = 15)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

ice_co2 <- ggplot(ice_all,mapping=aes(x=ice_period,y=NEE_uStar_orig,color=ice))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('on','off'),labels=c('On','Off'),values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())+
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
  geom_point(fcr_gf, mapping = aes(DateTime, ch4_flux_uStar_orig),alpha = 0.1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping = aes(x = Date, y = CH4, ymin = CH4-CH4_sd,ymax = CH4+CH4_sd),fill="#E63946",alpha=0.4)+
  geom_line(mapping = aes(Date, CH4),color="#E63946",size = 1)+
  xlim(as.POSIXct("2020-12-20"),as.POSIXct("2021-02-11"))+
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  ylim(-0.03,0.03)+
  theme_classic(base_size = 15)

ice_ch4 <- ggplot(ice_all,mapping=aes(x=ice_period,y=ch4_flux_uStar_orig,color=ice))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Ice Period")+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('on','off'),labels=c('On','Off'),values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(winter_co2,ice_co2,winter_ch4,ice_ch4,nrow=2,ncol=2,common.legend = TRUE)

ggsave("./Fig_Output/Ice_on_off.jpg",width = 8, height=8, units="in",dpi=320)

########################################3#
# plot meteorological variables


tfig <- fcr_daily %>% ggplot(aes(Date, Tmean)) + 
	geom_line(col = 'red', size = 1) +
	geom_ribbon(aes(ymin = Tmin, ymax = Tmax), fill = 'darkred', alpha = 0.2) +
	ylab(expression(T[a]~(degree~C))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw()  +theme(axis.text = element_text(size = 12, colour = 'black'),
										 axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
										 axis.title = element_text(size = 16, colour = 'black'))


pptfig <- fcr_daily %>% ggplot(aes(Date, PPT_tot)) + 
	geom_bar(stat= 'identity', fill = 'darkblue') +
	ylab(expression(Ppt~(mm~day^-1))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

parfig <- fcr_hourly %>% ggplot(aes(datetime, PAR_avg)) + 
	geom_line(col = 'goldenrod3') +
	ylab(expression(PAR~(mu~mol~m^-2~s^-1))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

rnfig <- fcr_hourly %>% ggplot(aes(datetime, Rn)) + 
	geom_line(col = 'goldenrod3') +
	ylab(expression(Rn~(W~m^-2))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


pressfig <- fcr_daily %>% ggplot(aes(Date, pressure)) + 
	geom_line(col = 'black', size = 1) +
	ylab(expression(P~(KPa))) +
	geom_ribbon(aes(ymin = minpress, ymax = maxpress), fill = 'black', alpha = 0.2) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ufig <- fcr_daily %>% ggplot(aes(Date, umean)) + 
	geom_line(col = 'black', size = 1) +
	ylab(expression(u~(m~s^-1))) +
	geom_ribbon(aes(ymin = umin, ymax = umax), fill = 'darkblue', alpha = 0.2) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


ggarrange(tfig, pptfig, 
					nrow = 2, align = "v")

ggarrange(ufig, pressfig, 
					nrow = 2, align = "v")



# daily totals 

fcr_dailysum %>% 
	ggplot(aes(Date, NEE_sum)) +
	geom_point() + geom_line() +
	#  geom_vline(xintercept = as.POSIXct('2020-06-29 -05'), lty = 2, col = 'red') +
	theme_bw() + ylab(expression(CO[2]~flux~(g~C~m^-2~d^-1))) +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-11-20 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1, lty =2) + 
	geom_hline(yintercept = 0, lty =2) +
	scale_x_datetime(date_breaks = '14 days', date_labels = '%b/%d') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


fcr_dailysum %>% 
	ggplot(aes(Date, CH4_sum)) +
	geom_point() + geom_line() +
	ylab(expression(CH[4]~(mg~C~m^-2~d^-1))) +
	#  geom_vline(xintercept = as.POSIXct('2020-06-29 -05'), lty = 2, col = 'red') +
	geom_hline(yintercept = 0, lty = 2) +
	theme_bw() +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#          alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-11-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1, lty =2) + 
	scale_x_datetime(date_breaks = '14 days', date_labels = '%b/%d') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


###################################################################
# CUMULATIVE VALUES
###################################################################


fcr_dailysum <- fcr_dailysum %>% ungroup() %>% 
	mutate(cumnee = cumsum(NEE_sum), 
				 cumch4 = cumsum(CH4_sum))


#

fcr_daily$cumnee <- fcr_dailysum$cumnee
fcr_daily$cumch4 <- fcr_dailysum$cumch4

coeff = 40

dailynee <- fcr_daily %>% 
	ggplot() +
	#geom_line(aes(Date, NEE_mean*3600*24*12/1000000)) +
	#geom_line(aes(Date, cumnee/coeff), col = 'darkgreen', size =1) +
  geom_line(aes(Date, cumnee), col = 'darkgreen', size =1) +
	xlab("") +
	# annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#         alpha = .2, fill = 'gray70') +
	#annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#        alpha = .2, fill = 'gray70') +
	#geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%d %b') +
	scale_y_continuous(name = expression(CO[2]~(g~C~m^-2~d^-1)), 
										 sec.axis = sec_axis(trans= ~.*coeff, 
										 										name = expression(Cumulative~flux~(g~C~m^-2~d^-1)))) +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 13, colour = 'black'))

coeff2 <- 30

dailych4 <- fcr_daily %>% 
	ggplot() +
	#geom_line(aes(Date, CH4_mean*3600*24*12/1000)) +
	#geom_line(aes(Date, cumch4/coeff2), col = 'goldenrod3', size = 1) +
  geom_line(aes(Date, cumch4/1000), col = 'goldenrod3', size = 1) +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	# annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#         alpha = .2, fill = 'gray70') +
	scale_y_continuous(name = expression(CH[4]~(g~C~m^-2~d^-1)), 
										 sec.axis = sec_axis(trans= ~.*coeff2, 
										 										name = expression(Cumulative~flux~(g~C~m^-2~d^-1)))) +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%d %b') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 13, colour = 'black'))


dailynee
dailych4

ggarrange(dailynee, dailych4, nrow = 2, ncol = 1, align = 'v')

########################################################

# montly mean

fcr_monthly <- fcr_daily %>% group_by(Year, Month) %>% 
	summarise(NEE = mean(NEE_mean, na.rm = TRUE),
						NEE05 = mean(NEE_U5, na.rm = TRUE),
						NEE50 = mean(NEE_U50, na.rm = TRUE),
						NEE95 = mean(NEE_U95, na.rm = TRUE),
						NEE_sd = sd(NEE_mean, na.rm = TRUE),
						CH4 = mean(CH4_mean, na.rm = TRUE),
						CH405 = mean(CH4_U5, na.rm = TRUE),
						CH450 = mean(CH4_U50, na.rm = TRUE),
						CH495 = mean(CH4_U95, na.rm = TRUE),
						CH4_sd = sd(CH4_mean, na.rm = TRUE),
						Tair = mean(Tmean, na.rm = TRUE),
						Tmin = min(Tmin, na.rm =TRUE),
						Tmax = max(Tmax, na.rm =TRUE),
						pressure = mean(pressure, na.rm =TRUE),
						maxpress = max(maxpress, na.rm = TRUE),
						minpress = min(minpress, na.rm = TRUE),
						H = mean(H, na.rm = TRUE),
						LE = mean(LE, na.rm = TRUE),
						VPD = mean(VPD, na.rm = TRUE),
						RH = mean(RH, na.rm = TRUE),
						u = mean(umean, na.rm = TRUE),
						umax = max(umax, na.rm = TRUE),
						umin = min(umin, na.rm = TRUE),
						PAR = mean(PAR, na.rm = TRUE),
						PPT = sum(PPT_tot, na.rm = TRUE),
						LW_in = mean(LW_in, na.rm = TRUE),
						LW_out = mean(LW_out, na.rm = TRUE),
						Rg = mean(Rg, na.rm = TRUE),
						SW_out = mean(SW_out, na.rm = TRUE),
						Rn = mean(Rn, na.rm = TRUE))



# monthly sum NEE in gCm2month, CH4 in mgCm2month

fcr_monthlysum <- fcr_dailysum %>% 
	dplyr::group_by(Year, Month) %>% 
	dplyr::summarise(NEE = sum(NEE_sum, na.rm = TRUE),
									 NEEU5 = sum(NEE_sum05, na.rm = TRUE),
									 NEEU50 = sum(NEE_sum50, na.rm = TRUE),
									 NEEU95 = sum(NEE_sum95, na.rm = TRUE),
									 NEE_sd = sd(NEE_sum, na.rm = TRUE),
									 CH4 = sum(CH4_sum, na.rm = TRUE),
									 CH4U5 = sum(CH4_sum05, na.rm = TRUE),
									 CH4U50 = sum(CH4_sum50, na.rm = TRUE),
									 CH4U95 = sum(CH4_sum95, na.rm = TRUE),
									 CH4_sd = sd(CH4_sum, na.rm = TRUE),
									 Tair = mean(Tmean, na.rm = TRUE),
									 Tmax = max(Tmax, na.rm = TRUE),
									 Tmin = min(Tmin, na.rm = TRUE),
									 H = sum(H, na.rm = TRUE),
									 LE = sum(LE, na.rm = TRUE),
									 Rn = sum(Rn, na.rm = TRUE),
									 VPD = mean(VPD, na.rm = TRUE),
									 RH = mean(RH, na.rm = TRUE),
									 u = mean(umean, na.rm = TRUE),
									 umax= max(umax, na.rm = TRUE),
									 PAR = sum(PAR, na.rm = TRUE),
									 PPT = sum(PPT_tot, na.rm = TRUE),
									 LW_in = mean(LW_in, na.rm = TRUE),
									 LW_out = mean(LW_out, na.rm = TRUE),
									 Rg = mean(Rg, na.rm = TRUE),
									 SW_out = mean(SW_out, na.rm = TRUE))


# monthly bar plot with totals NEE and CH4

fcr_monthlysum$yearmon <- with(fcr_monthlysum, sprintf("%d-%02d", Year, Month))

neemon <- fcr_monthlysum %>% ggplot() +
	geom_bar(aes(yearmon, NEE), stat = 'identity', fill = 'darkgreen', alpha = 0.8, width = 0.3) +
	geom_errorbar(aes(yearmon, ymin = NEE - NEE_sd, ymax = NEE + NEE_sd), width = 0.2) +
	ylab(expression(~CO[2]~(g~C~m^-2~month^-1))) +
	xlab("") +
	theme_bw() +
	theme(axis.text = element_text(size = 14, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ch4mon <- fcr_monthlysum %>% ggplot() +
	geom_bar(aes(yearmon, CH4), stat = 'identity', fill = 'goldenrod3', alpha = 0.6, width = 0.3) +
	geom_errorbar(aes(yearmon, ymin = CH4 - CH4_sd, ymax = CH4 + CH4_sd), width = 0.2) +
	xlab('') + ylab(expression(~CH[4]~(mg~C~m^-2~month^-1))) +
	theme_bw() +
	theme(axis.text = element_text(size = 14, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ggarrange(neemon, ch4mon, align = 'v')


# get some numbers:

# hourly means

# CO2 in mmolm2h
mean(fcr_hourly$NEE*3600/1000)

# in mgCO2m2h
mean(fcr_hourly$NEE*3600*44/1000)

# CH4 in umolm2h
mean(fcr_hourly$CH4*3600)

# CH4 in ugCH4m2h
mean(fcr_hourly$CH4*3600*16)


# daily means

# mean of the sum - CO2 in gCm2d
mean(fcr_dailysum$NEE_sum)

# mean of the mean
mean(fcr_daily$NEE_mean*3600*24*12/1000000)

# mean of the sum - CH4 in mgCm2d
mean(fcr_dailysum$CH4_sum)

# mean of the mean mgCm2d
mean(fcr_daily$CH4_mean*3600*24*12/1000)

# monthly means
mean(fcr_monthly$NEE*3600*24*12*30/1000000) #gCm2month
mean(fcr_monthly$CH4*3600*24*12*30/1000) #mgCm2month

mean(fcr_monthlysum$NEE) #gCm2month
mean(fcr_monthlysum$CH4) #mgCm2month


# total values (more than 1 year of data)

fcr_monthlysum %>% select(yearmon, NEE, NEEU5, NEEU50, NEEU95, CH4, CH4U5, CH4U50, CH4U95)


annualco2 <- sum(fcr_monthlysum$NEE) # gCm2
sdco2 <- sqrt(sum(fcr_monthlysum$NEE_sd^2, na.rm = TRUE))
annualch4 <- sum(fcr_monthlysum$CH4) # mgCm2
sdch4 <- sqrt(sum(fcr_monthlysum$CH4_sd^2, na.rm = TRUE))

annualco2
sdco2
annualch4
sdch4



#########

# Looking for some relationships
# temp
fcr_dailysum %>% ggplot(aes(Tmean, NEE_sum)) + 
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~Temperature~(degree~C))) +
	geom_point() + theme_bw() +
	geom_smooth(span = 0.9)

templm <- lm(log(fcr_dailysum$NEE_sum)~fcr_dailysum$Tmean)
summary(templm)

# LW

fcr_daily %>% ggplot(aes(LW_in, NEE_mean)) + 
	geom_point() + theme_bw() +
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~LW[incoming]~(W~m^-2))) +
	geom_smooth(span = 0.9) 

lwco2lm <- lm(log(fcr_daily$NEE_mean)~fcr_daily$LW_in)
summary(lwco2lm)



fcr_dailysum %>% ggplot(aes(Rn, CH4_sum)) + 
	geom_point() + theme_bw() +
	ylab(expression(~CH4~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~Rn~(W~m^-2~s^-1))) +
	geom_smooth(method= 'lm')

rnlm <- lm(fcr_daily$CH4_mean~fcr_daily$Rn)
summary(rnlm)
