### Script conduct analysis of EC fluxes from FCR
### Apr 2020 - Apr 2022

### Following revisions from initial submission to JGR-Biogeosciences
### 12 May 2022, A. Hounshell

### Updated to add in diffusive fluxes, 20 May 2022, A. Hounshell
### Updated with all diffusive fluxes (through Spring 2022), 17 June 2022

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate,openair)

###############################################################################

## First load in QA/QC'd EC fluxes 

## Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
## Data corrected following FCR_Process_BD
## Data downloaded from: https://pasta-s.lternet.edu/package/data/eml/edi/920/2/9e658ef44de05303dbc496fc25e8c49a

ec <- read.csv("./Data/20220506_EC_processed.csv") 

ec2 <- ec %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00"))

###############################################################################

## Visualizations of missing data - following all QA/QC

## Amount of data retained following filtering due to low u*
ec2 %>% select(DateTime, NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  summarise(co2_available = 100-sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)
# 23% CO2 fluxes; 19% CH4 fluxes

## Separate by year
ec2 %>% select(DateTime, NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  filter(DateTime < "2021-05-01") %>% 
  summarise(co2_available = 100-sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)
# 22% CO2 fluxes; 16% Ch4 fluxes for year 1

ec2 %>% select(DateTime, NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  filter(DateTime >= "2021-05-01") %>% 
  summarise(co2_available = 100-sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)
# 24% CO2 fluxes; 23% CH4 fluxes in year 2

## Distribution of missing data by season, day vs. night, etc.
## Fig S4
day_co2_data <- ec2 %>% 
  mutate(time = format(as.POSIXct(DateTime,format='%Y-%m-%d %H:%M:%S'),format='%H:%M:%S')) %>% 
  select(time,NEE_uStar_orig) %>% 
  drop_na(NEE_uStar_orig) %>% 
  count(time)

day_ch4_data <- ec2 %>% 
  mutate(time = format(as.POSIXct(DateTime,format='%Y-%m-%d %H:%M:%S'),format='%H:%M:%S')) %>% 
  select(time,ch4_flux_uStar_orig) %>% 
  drop_na(ch4_flux_uStar_orig) %>% 
  count(time)

co2_time <- ggplot(day_co2_data,mapping=aes(x=as.factor(time),y=n/730*100))+
  geom_hline(yintercept = 10, color = "lightgrey")+
  geom_hline(yintercept = 20, color = "lightgrey")+
  geom_hline(yintercept = 30, color = "lightgrey")+
  geom_bar(stat="identity")+
  ylab(expression(~Percent~CO[2]~Data)) +
  xlab("Time")+
  ylim(0,35)+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90))

ch4_time <- ggplot(day_ch4_data,mapping=aes(x=as.factor(time),y=n/730*100))+
  geom_hline(yintercept = 10, color = "lightgrey")+
  geom_hline(yintercept = 20, color = "lightgrey")+
  geom_hline(yintercept = 30, color = "lightgrey")+
  geom_bar(stat="identity")+
  ylab(expression(~Percent~CH[4]~Data)) +
  xlab("Time")+
  ylim(0,35)+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90))

ggarrange(co2_time,ch4_time,ncol=1,nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Data_Time.jpg",width = 8, height=8, units="in",dpi=320)

## Look at data availability seasonally
# Fig S6
season_data_co2 <- ec2 %>% 
  group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
                                                                                        "Jul", "Aug", "Sep", "Oct", 'Nov', 
                                                                                        'Dec', 'Jan', 'Feb', 'Mar'))) %>% 
  count(month,!is.na(NEE_uStar_orig))

season_data_ch4 <- ec2 %>% 
  group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
                                                                                        "Jul", "Aug", "Sep", "Oct", 'Nov', 
                                                                                        'Dec', 'Jan', 'Feb', 'Mar'))) %>% 
  count(month,!is.na(ch4_flux_uStar_orig))

season_data_total <- season_data_co2 %>% 
  group_by(year,month) %>% 
  summarise(sum(n))

season_co2 <- season_data_co2 %>% 
  filter(`!is.na(NEE_uStar_orig)` == "TRUE") %>% 
  ggplot(mapping=aes(x=month,y=n/season_data_total$`sum(n)`*100,fill=as.factor(year)))+
  geom_hline(yintercept = 10, color = "lightgrey")+
  geom_hline(yintercept = 20, color = "lightgrey")+
  geom_hline(yintercept = 30, color = "lightgrey")+
  geom_bar(stat="identity",position=position_dodge())+
  ylab(expression(~Percent~CO[2]~Data)) +
  xlab("Month")+
  scale_fill_brewer(palette = "Dark2")+
  ylim(0,35)+
  guides(fill=guide_legend(title="Year"))+
  theme_classic(base_size = 15)

season_ch4 <- season_data_ch4 %>% 
  filter(`!is.na(ch4_flux_uStar_orig)` == "TRUE") %>% 
  ggplot(mapping=aes(x=month,y=n/season_data_total$`sum(n)`*100,fill=as.factor(year)))+
  geom_hline(yintercept = 10, color = "lightgrey")+
  geom_hline(yintercept = 20, color = "lightgrey")+
  geom_hline(yintercept = 30, color = "lightgrey")+
  geom_bar(stat="identity",position=position_dodge())+
  ylab(expression(~Percent~CH[4]~Data)) +
  xlab("Month")+
  scale_fill_brewer(palette = "Dark2")+
  ylim(0,35)+
  guides(fill=guide_legend(title="Year"))+
  theme_classic(base_size = 15)

ggarrange(season_co2,season_ch4,ncol=2,nrow=1,labels=c("A.","B."),font.label = list(face="plain",size=15),common.legend = TRUE)

ggsave("./Fig_Output/SI_Data_Season.jpg",width = 9, height=4, units="in",dpi=320)

###############################################################################

## Visualize measured fluxes at the hourly, daily, weekly, and monthly timescale
fcr_hourly <- ec2 %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d %H"),"%Y-%m-%d %H")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
            NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
            CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
            CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

# Calculate min, max, median, mean, standard deviation, and coefficient of variation
# Table S4
fcr_stats <- fcr_hourly %>% 
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

write.csv(fcr_stats,"./Fig_output/20220506_TableSx_ECStats.csv")

# Aggregate to daily and calculate the variability (SD) - following script for figures_BD
# data in umolm2s
fcr_daily <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_daily$Date <- as.POSIXct(paste(fcr_daily$Year, fcr_daily$Month, fcr_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# Aggregate into weekly
fcr_weekly <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Week = week(DateTime)) %>% 
  dplyr::group_by(Year, Week) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_weekly$Date <- make_datetime(year = fcr_weekly$Year) + weeks(fcr_weekly$Week)

# Aggregate to Monthly
fcr_monthly <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime)) %>% 
  dplyr::group_by(Year, Month) %>% 
  dplyr::summarise(NEE = mean(NEE_uStar_orig, na.rm = TRUE),
                   NEE_sd = sd(NEE_uStar_orig, na.rm = TRUE),
                   CH4 = mean(ch4_flux_uStar_orig, na.rm = TRUE),
                   CH4_sd = sd(ch4_flux_uStar_orig, na.rm = TRUE))

fcr_monthly$yearmon <- with(fcr_monthly, sprintf("%d-%02d", Year, Month))

###############################################################################

## Load in diffusive fluxes:
## Calculated using: 1b_Diffusive_Fluxes.R
diff_flux <- read.csv("./Data/20220617_diffusive_fluxes_avg.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST"))

diff_flux <- diff_flux %>% 
  drop_na()

## Compare EC and Diff fluxes at turnover
ec_turnover_year1 <- ec2 %>% 
  filter(DateTime >= "2020-10-31 00:00:00" & DateTime < "2020-11-03 00:00:00") %>% 
  dplyr::summarise(ch4_mean = mean(ch4_flux_uStar_orig,na.rm=TRUE),
                   co2_mean = mean(NEE_uStar_orig, na.rm=TRUE))

ec_turnover_year2 <- ec2 %>% 
  filter(DateTime >= "2021-11-02 00:00:00" & DateTime < "2021-11-05 00:00:00") %>% 
  dplyr::summarise(ch4_mean = mean(ch4_flux_uStar_orig,na.rm=TRUE),
                   co2_mean = mean(NEE_uStar_orig, na.rm=TRUE))

###############################################################################

## Plot daily and monthly for both years
# Fig 2
co2_daily_year1 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-02-26"), col='black', size=1,linetype="dotted")+
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), col='black', size=1,linetype="dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,ymin=co2_mean_umol_m2_s-co2_sd_umol_m2_s,ymax=co2_mean_umol_m2_s+co2_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_daily_year2 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-02-26"), col='black', size=1,linetype="dotted")+
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), col='black', size=1,linetype="dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,ymin=co2_mean_umol_m2_s-co2_sd_umol_m2_s,ymax=co2_mean_umol_m2_s+co2_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_month <- ggplot(fcr_monthly) +
  geom_vline(xintercept = "2021-01", color="lightgrey")+
  geom_vline(xintercept = "2022-01", color="lightgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(yearmon, ymin = NEE - NEE_sd, ymax = NEE + NEE_sd), width = 0.2) +
  geom_point(aes(yearmon, NEE), color = "#E63946", size=4) +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("") +
  theme_classic(base_size=15) +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(colour = 'black'))

ggarrange(co2_daily_year1,co2_daily_year2,co2_month,ncol=1,
          nrow=3,labels=c("A.","B.","C."),font.label = list(face="plain",size=15), common.legend = TRUE)

ggsave("./Fig_Output/CO2_Daily_Monthly.jpg",width = 9, height=12, units="in",dpi=320)

## Plot CH4 daily and monthly for MS
# Fig. 3
ch4_daily_year1 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-02-26"), col='black', size=1,linetype="dotted")+
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), col='black', size=1,linetype="dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,ymin=ch4_mean_umol_m2_s-ch4_sd_umol_m2_s,ymax=ch4_mean_umol_m2_s+ch4_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_daily_year2 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-02-26"), col='black', size=1,linetype="dotted")+
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), col='black', size=1,linetype="dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,ymin=ch4_mean_umol_m2_s-ch4_sd_umol_m2_s,ymax=ch4_mean_umol_m2_s+ch4_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_month <- ggplot(fcr_monthly) +
  geom_vline(xintercept = "2021-01", color="lightgrey")+
  geom_vline(xintercept = "2022-01", color="lightgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(yearmon, ymin = CH4 - CH4_sd, ymax = CH4 + CH4_sd), width = 0.2) +
  geom_point(aes(yearmon, CH4), color = "#E63946", size=4) +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("") +
  theme_classic(base_size=15) +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(colour = 'black'))

ggarrange(ch4_daily_year1,ch4_daily_year2,ch4_month,ncol=1,
          nrow=3,labels=c("A.","B.","C."),font.label = list(face="plain",size=15), common.legend = TRUE)

ggsave("./Fig_Output/CH4_Daily_Monthly.jpg",width = 9, height=12, units="in",dpi=320)

###############################################################################

## Compare discrete diffusive fluxes and EC fluxes by timepoint
# First combine by the hour
diff_flux_hr <- diff_flux %>% 
  mutate(DateTime = format(round(DateTime, units="hours"), format="%Y-%m-%d %H:%M")) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST")))

diff_ec_hr <- left_join(diff_flux_hr,fcr_hourly,by="DateTime")

diff_ec_hr %>% 
  drop_na(NEE) %>% 
  count()

diff_ec_hr %>% 
  drop_na(CH4) %>% 
  count()

diff_ec_hr_long_co2 <- diff_ec_hr %>% 
  select(DateTime,co2_mean_umol_m2_s,NEE) %>% 
  pivot_longer(!DateTime, names_to = "Type", values_to = "flux_umol_m2_s") %>% 
  mutate(Type = ifelse(Type == "co2_mean_umol_m2_s","Diff",
                       ifelse(Type == "NEE", "EC", NA)),
         Flux = "co2")

diff_ec_hr_long_ch4 <- diff_ec_hr %>% 
  select(DateTime,ch4_mean_umol_m2_s,CH4) %>% 
  pivot_longer(!DateTime, names_to = "Type", values_to = "flux_umol_m2_s") %>% 
  mutate(Type = ifelse(Type == "ch4_mean_umol_m2_s","Diff",
                       ifelse(Type == "CH4", "EC", NA)),
         Flux = "ch4")

## Plot
# Fig S11
diff_ec_co2 <- ggplot(diff_ec_hr,mapping=aes(x=co2_mean_umol_m2_s,y=NEE))+
  geom_abline(intercept = 0)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_errorbar(mapping=aes(x = co2_mean_umol_m2_s, y = NEE, xmin = co2_mean_umol_m2_s-co2_sd_umol_m2_s, xmax = co2_mean_umol_m2_s+co2_sd_umol_m2_s),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_errorbar(mapping=aes(x=co2_mean_umol_m2_s,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_point(mapping=aes(x=co2_mean_umol_m2_s,y=NEE))+
  ylab(expression(~EC~CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab(expression(~Mean~Diff~CO[2]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

co2_box <- ggplot(diff_ec_hr_long_co2,mapping=aes(x=Type,y=flux_umol_m2_s))+
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(alpha=0.3)+
  ylab(expression(CO[2]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic()

diff_ec_ch4 <- ggplot(diff_ec_hr,mapping=aes(x=ch4_mean_umol_m2_s,y=CH4))+
  geom_abline(intercept = 0)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_errorbar(mapping=aes(x = ch4_mean_umol_m2_s, y = CH4, xmin = ch4_mean_umol_m2_s-ch4_sd_umol_m2_s, xmax = ch4_mean_umol_m2_s+ch4_sd_umol_m2_s),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_errorbar(mapping=aes(x=ch4_mean_umol_m2_s,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),position = position_dodge(0.3), width = 0.2,color="lightgrey")+
  geom_point(mapping=aes(x=ch4_mean_umol_m2_s,y=CH4))+
  ylab(expression(~EC~CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab(expression(~Mean~Diff~CH[4]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size = 15)

ch4_box <- ggplot(diff_ec_hr_long_ch4,mapping=aes(x=Type,y=flux_umol_m2_s))+
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(alpha=0.3)+
  ylab(expression(CH[4]~(mu~mol~m^-2~s^-1))) +
  xlab("")+
  theme_classic()

ggarrange(diff_ec_co2,co2_box,diff_ec_ch4,ch4_box,ncol=2,nrow=2,labels=c("A.","B.","C.","D."),
          font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Diff_EC.jpg",width = 8, height=7, units="in",dpi=320)

###############################################################################

## Plot weekly for both years
# Fig. S9
co2_week <- ggplot(fcr_weekly)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=Date,y=NEE),color="#E63946",size = 1)+
  geom_point(mapping=aes(x=Date,y=NEE),color="#E63946",size=2)+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size=15)

ch4_week <- ggplot(fcr_weekly)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.3)+
  geom_line(mapping=aes(x=Date,y=CH4),color="#E63946",size = 1)+
  geom_point(mapping=aes(x=Date,y=CH4),color="#E63946",size=2)+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  theme_classic(base_size=15)

ggarrange(co2_week,ch4_week,ncol=1,
          nrow=2,labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/EC_Weekly.jpg",width = 8, height = 8, units="in",dpi=320)

###############################################################################

## Look at fluxes around turnover - both EC and diffusive

## Load in wind speed data from Met first to plot with the fluxes
## Downloaded from EDI: 06 May 2022
## Included preliminary met data for 2022 downloaded from the Gateway

#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

met_edi <- read.csv("./Data/Met_final_2015_2021.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

# Load in preliminary Met 2022 data (from 1a_Rev_2022_MetData.R)
met_2022 <- read.csv("./Data/FCR_Met_final_2022.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

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

## Aggregate to daily
met_daily <- met_30_2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime), 
         Hour = hour(DateTime)) %>% 
  dplyr::group_by(Year, Month, Day) %>% 
  dplyr::summarise(wind_ms = mean(WindSpeed_Average_m_s, na.rm = TRUE),
                   wind_ms_sd = sd(WindSpeed_Average_m_s, na.rm = TRUE))

met_daily$DateTime <- as.POSIXct(paste(met_daily$Year, met_daily$Month, met_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

## Plot fluxes and wind speed
# Fig. S9
co2_2020 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,ymin=co2_mean_umol_m2_s-co2_sd_umol_m2_s,ymax=co2_mean_umol_m2_s+co2_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("2020") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-10-01"),as.POSIXct("2020-11-30"))+
  ylim(-30,30)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_2021 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,ymin=co2_mean_umol_m2_s-co2_sd_umol_m2_s,ymax=co2_mean_umol_m2_s+co2_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=co2_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("2021") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-10-01"),as.POSIXct("2021-11-30"))+
  ylim(-30,30)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_2020 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,ymin=ch4_mean_umol_m2_s-ch4_sd_umol_m2_s,ymax=ch4_mean_umol_m2_s+ch4_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("2020") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-10-01"),as.POSIXct("2020-11-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_2021 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig,color="30 min EC fluxes"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="Daily Mean EC"),size = 1) +
  geom_errorbar(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,ymin=ch4_mean_umol_m2_s-ch4_sd_umol_m2_s,ymax=ch4_mean_umol_m2_s+ch4_sd_umol_m2_s,color="Diff"),size=1)+
  geom_point(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"),size=2)+
  geom_line(diff_flux,mapping=aes(x=DateTime,y=ch4_mean_umol_m2_s,color="Diff"))+
  scale_color_manual(breaks=c("30 min EC fluxes","Daily Mean EC","Diff"),
                     values=c("black","#E63946","#4c8bfe"))+
  xlab("2021") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-10-01"),as.POSIXct("2021-11-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

wnd_2020 <- met_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(met_30_2,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s,color="30 min Wind"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=DateTime,y=wind_ms,ymin=wind_ms-wind_ms_sd,ymax=wind_ms+wind_ms_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(DateTime, wind_ms,color="Daily Mean Wind"),size = 1) +
  scale_color_manual(breaks=c("30 min Wind","Daily Mean Wind"),
                     values=c("black","#E63946"))+
  xlab("2020") +
  ylab(expression(~Wind~Speed~(m~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-10-01"),as.POSIXct("2020-11-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

wnd_2021 <- met_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct("2021-11-03 12:00"), col = "black", size = 1, linetype = "dotted")+
  geom_point(met_30_2,mapping=aes(x=DateTime,y=WindSpeed_Average_m_s,color="30 min Wind"),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=DateTime,y=wind_ms,ymin=wind_ms-wind_ms_sd,ymax=wind_ms+wind_ms_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(DateTime, wind_ms,color="Daily Mean Wind"),size = 1) +
  scale_color_manual(breaks=c("30 min Wind","Daily Mean Wind"),
                     values=c("black","#E63946"))+
  xlab("2021") +
  ylab(expression(~Wind~Speed~(m~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-10-01"),as.POSIXct("2021-11-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(co2_2020,co2_2021,ch4_2020,ch4_2021,wnd_2020,wnd_2021,ncol=2,
          nrow=3,labels=c("A.","B.","C.","D.","E.","F."),font.label = list(face="plain",size=15),common.legend = TRUE)

ggsave("./Fig_Output/SI_Fluxes_Turnover.jpg",width = 11, height = 10, units="in",dpi=320)

###############################################################################

## Compare day to night
## Noon = 11-1 pm; Midnight = 23-1 am
diel_flux <- ec2 %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,u) %>% 
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

season_diel <- diel_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

write.csv(season_diel,"./Fig_output/20221108_TableSx_DielSeason.csv")

## Look at dawn vs. dusk differences in fluxes
dawn_flux <- ec2 %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,u) %>% 
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

season_dawn <- dawn_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

write.csv(season_dawn,"./Fig_output/20221108_TableSx_DawnSeason.csv")

## Calculate statistics for Supplementary Table (Table S5) - Day vs. Night
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
write.csv(diel_stats,"./Fig_output/20220513_TableSx_DielStats.csv")

## Calculate statistics for Supplementary Table (Table S5) - Dawn vs. Dusk
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
write.csv(dawn_stats,"./Fig_output/20220513_TableSx_DawnStats.csv")

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

## Calculate paired-wilcoxon signed rank tests for each season
daynight_wide_co2_f_spring <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(daynight_wide_co2_f_spring$Day,daynight_wide_co2_f_spring$Night,paired=TRUE)

daynight_wide_co2_f_summer <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(daynight_wide_co2_f_summer$Day,daynight_wide_co2_f_summer$Night,paired=TRUE)

daynight_wide_co2_f_fall <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(daynight_wide_co2_f_fall$Day,daynight_wide_co2_f_fall$Night,paired=TRUE)

daynight_wide_co2_f_winter <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(daynight_wide_co2_f_winter$Day,daynight_wide_co2_f_winter$Night,paired=TRUE)

daynight_wide_ch4_f_spring <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(daynight_wide_ch4_f_spring$Day,daynight_wide_ch4_f_spring$Night,paired=TRUE)

daynight_wide_ch4_f_summer <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(daynight_wide_ch4_f_summer$Day,daynight_wide_ch4_f_summer$Night,paired=TRUE)

daynight_wide_ch4_f_fall <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(daynight_wide_ch4_f_fall$Day,daynight_wide_ch4_f_fall$Night,paired=TRUE)

daynight_wide_ch4_f_winter <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(daynight_wide_ch4_f_winter$Day,daynight_wide_ch4_f_winter$Night,paired=TRUE)

daynight_wide_u_f_spring <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(daynight_wide_u_f_spring$Day,daynight_wide_u_f_spring$Night,paired=TRUE)

daynight_wide_u_f_summer <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(daynight_wide_u_f_summer$Day,daynight_wide_u_f_summer$Night,paired=TRUE)

daynight_wide_u_f_fall <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(daynight_wide_u_f_fall$Day,daynight_wide_u_f_fall$Night,paired=TRUE)

daynight_wide_u_f_winter <- diel_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(daynight_wide_u_f_winter$Day,daynight_wide_u_f_winter$Night,paired=TRUE)

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

## Calculate paired-wilcoxon signed rank tests for each season
dawndusk_wide_co2_f_spring <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(dawndusk_wide_co2_f_spring$Dawn,dawndusk_wide_co2_f_spring$Dusk,paired=TRUE)

dawndusk_wide_co2_f_summer <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(dawndusk_wide_co2_f_summer$Dawn,dawndusk_wide_co2_f_summer$Dusk,paired=TRUE)

dawndusk_wide_co2_f_fall <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(dawndusk_wide_co2_f_fall$Dawn,dawndusk_wide_co2_f_fall$Dusk,paired=TRUE)

dawndusk_wide_co2_f_winter <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,NEE) %>% 
  pivot_wider(names_from = diel, values_from = NEE) %>% 
  drop_na

wilcox.test(dawndusk_wide_co2_f_winter$Dawn,dawndusk_wide_co2_f_winter$Dusk,paired=TRUE)

dawndusk_wide_ch4_f_spring <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(dawndusk_wide_ch4_f_spring$Dawn,dawndusk_wide_ch4_f_spring$Dusk,paired=TRUE)

dawndusk_wide_ch4_f_summer <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(dawndusk_wide_ch4_f_summer$Dawn,dawndusk_wide_ch4_f_summer$Dusk,paired=TRUE)

dawndusk_wide_ch4_f_fall <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(dawndusk_wide_ch4_f_fall$Dawn,dawndusk_wide_ch4_f_fall$Dusk,paired=TRUE)

dawndusk_wide_ch4_f_winter <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,CH4) %>% 
  pivot_wider(names_from = diel, values_from = CH4) %>% 
  drop_na

wilcox.test(dawndusk_wide_ch4_f_winter$Dawn,dawndusk_wide_ch4_f_winter$Dusk,paired=TRUE)

dawndusk_wide_u_f_spring <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Spring") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(dawndusk_wide_u_f_spring$Dawn,dawndusk_wide_u_f_spring$Dusk,paired=TRUE)

dawndusk_wide_u_f_summer <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Summer") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(dawndusk_wide_u_f_summer$Dawn,dawndusk_wide_u_f_summer$Dusk,paired=TRUE)

dawndusk_wide_u_f_fall <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Fall") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(dawndusk_wide_u_f_fall$Dawn,dawndusk_wide_u_f_fall$Dusk,paired=TRUE)

dawndusk_wide_u_f_winter <- dawn_agg %>% 
  ungroup() %>% 
  filter(season == "Winter") %>% 
  select(Date,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na

wilcox.test(dawndusk_wide_u_f_winter$Dawn,dawndusk_wide_u_f_winter$Dusk,paired=TRUE)

## Plot diel (day/night) and dawn/dusk comparisons
# Fig. 4
diel_co2 <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label ="p = 0.09
n = 295",
           x = "Day", hjust = 0, y = 40, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol m"^-2*" s"^-1*")")))+
  xlab("")+
  #ylim(-20,20)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

diel_ch4 <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "p = 0.16
n = 246",
           x = "Day", hjust = 0, y = 0.07, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 

diel_wind <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=u,color=diel))+
  annotate("text",label = "p < 0.001*
n = 502",
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
n = 189",
           x = "Dawn", hjust = 0, y = 40, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol m"^-2*" s"^-1*")")))+
  xlab("")+
  #ylim(-20,20)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

dawn_ch4 <- dawn_agg %>% 
  ggplot(mapping=aes(x=diel,y=CH4,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "p = 0.36
n = 152",
           x = "Dawn", hjust = 0, y = 0.045, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

dawn_wind <- dawn_agg %>% 
  ggplot(mapping=aes(x=diel,y=u,color=diel))+
  annotate("text",label = "p < 0.005*
n = 312",
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

ggsave("./Fig_Output/Rev2_Figure3.png",width = 9, height=7.5, units="in",dpi=320)

###############################################################################

## Ice comparisons between years: 2020-2021 = intermittent ice; 2021-2022 = ice-on

## Download Ice data from EDI - downloaded on 13 May 2022
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/456/4/8454a757d0203bd8913d8adc8607f6c4" 
#infile1 <- paste0(getwd(),"/Data/Ice_Data.csv")
#download.file(inUrl1,infile1,method="curl")

ice <- read_csv("./Data/Ice_Data.csv")
ice$Date <- as.POSIXct(strptime(ice$Date,"%Y-%m-%d"))

ice_fcr <- ice %>% 
  filter(Reservoir == "FCR" & Date>"2020-05-01")

## Plot to visualize and find ice on vs. ice off
ggplot(ice_fcr)+
  geom_point(mapping=aes(x=Date,y=IceOn),color="blue")+
  geom_point(mapping=aes(x=Date,y=IceOff),color="red")
# Use the time period from: 16 Jan to 10 Feb for each year ("ice-on")
# Include a week before and after, too: 09 Jan to 17 Feb

ice_fluxes_daily <- fcr_daily %>% 
  mutate(Year = year(Date), 
         Month = month(Date), 
         Day = day(Date)) %>% 
  filter(Month %in% c(1,2))

ice_fluxes_30min <- ec2 %>% 
  mutate(Year = year(DateTime), 
         Month = month(DateTime), 
         Day = day(DateTime),
         Ice = ifelse(DateTime >= as.POSIXct("2021-01-10") & DateTime < as.POSIXct("2021-02-09"), "Ice_2021", 
                      ifelse(DateTime >= as.POSIXct("2022-01-16") & DateTime < as.POSIXct("2022-02-10"),"Ice_2022", NA))) %>% 
  filter(Month %in% c(1,2)) %>% 
  select(DateTime,Year,Month,Day,Ice,NEE_uStar_orig,ch4_flux_uStar_orig)

## Compare ice on to ice off for 2021 and 2022
ice_fluxes_30min_comps <- ice_fluxes_30min %>% 
  filter(Ice %in% c("Ice_2021","Ice_2022"))

# Table S9
ice_stats <- ice_fluxes_30min_comps %>% 
  group_by(Ice) %>% 
  summarise(p25_nee = quantile(NEE_uStar_orig,0.25,na.rm=TRUE),
            med_nee = median(NEE_uStar_orig,na.rm=TRUE),
            p75_nee = quantile(NEE_uStar_orig,0.75,na.rm=TRUE),
            p25_ch4 = quantile(ch4_flux_uStar_orig,0.25,na.rm=TRUE),
            med_ch4 = median(ch4_flux_uStar_orig,na.rm=TRUE),
            p75_ch4 = quantile(ch4_flux_uStar_orig,0.75,na.rm=TRUE))

write.csv(ice_stats,"./Fig_output/Ice_Stats.csv")

wilcox.test(NEE_uStar_orig ~ Ice, data=ice_fluxes_30min_comps)
wilcox.test(ch4_flux_uStar_orig ~ Ice, data=ice_fluxes_30min_comps)

## Plot ice on/ice off for 2021 and 2022
# Fig. 6
co2_year1 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2021-02-26"), linetype = "dotted", color="black",size=1.3)+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=NEE),color="#E63946",size=1)+
  xlim(as.POSIXct("2021-01-01"),as.POSIXct("2021-02-28"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

ch4_year1 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2021-02-23"), linetype = "dotted", color="red",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2021-02-26"), linetype = "dotted", color="black",size=1.3)+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=CH4),color="#E63946",size=1)+
  xlim(as.POSIXct("2021-01-01"),as.POSIXct("2021-02-28"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

co2_year2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2022-01-11"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-01-14"), linetype = "dotted", color="red",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-01-16"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), linetype = "dotted", color="red",size=1.3)+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=NEE),color="#E63946",size=1)+
  xlim(as.POSIXct("2022-01-01"),as.POSIXct("2022-02-28"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

ch4_year2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2022-01-11"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-01-14"), linetype = "dotted", color="red",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-01-16"), linetype = "dotted", color="blue",size=1.3)+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), linetype = "dotted", color="red",size=1.3)+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=CH4),color="#E63946",size=1)+
  xlim(as.POSIXct("2022-01-01"),as.POSIXct("2022-02-28"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

co2_comps <- ggplot(ice_fluxes_30min_comps,mapping=aes(x=Ice,y=NEE_uStar_orig,color=Ice))+
  annotate("text",label = "p < 0.001*",
           x = "Ice_2021", hjust = 0, y = 9, size = 5)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Year")+
  geom_boxplot(outlier.shape = NA,size=0.7)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Ice_2021','Ice_2022'),values=c("#008080","#477CC2"))+
  scale_x_discrete(labels=c("Ice_2021" = "2021", "Ice_2022" = "2022"))+
  theme_classic(base_size = 15)+
  theme(legend.position="none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

ch4_comps <- ggplot(ice_fluxes_30min_comps,mapping=aes(x=Ice,y=ch4_flux_uStar_orig,color=Ice))+
  annotate("text",label = "p < 0.001*",
           x = "Ice_2021", hjust = 0, y = 0.025, size = 5)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("Year")+
  geom_boxplot(outlier.shape = NA,size=0.7)+
  geom_point(position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(breaks=c('Ice_2021','Ice_2022'),values=c("#008080","#477CC2"))+
  scale_x_discrete(labels=c("Ice_2021" = "2021", "Ice_2022" = "2022"))+
  theme_classic(base_size = 15)+
  theme(legend.position="none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

ggarrange(co2_year1,ch4_year1,co2_year2,ch4_year2,co2_comps,ch4_comps,ncol=2,nrow=3,
          labels=c("A.","B.","C.","D.","E.","F."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/Figure_WinterIce.png",width = 10, height=9, units="in",dpi=320)

###############################################################################

## Estimate and plot annual fluxes for Year 1 and Year 2, CO2 and CH4
ec2_annual_fluxes <- ec2 %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,NEE_uStar_f,ch4_flux_uStar_f,NEE_uStar_fsd,ch4_flux_uStar_fsd) %>% 
  mutate(NEE_uStar_orig = ifelse(is.na(NEE_uStar_orig), 0, NEE_uStar_orig),
         ch4_flux_uStar_orig = ifelse(is.na(ch4_flux_uStar_orig), 0, ch4_flux_uStar_orig))

ec2_annual_fluxes_year1 <- ec2_annual_fluxes %>% 
  filter(DateTime < as.POSIXct("2021-05-01")) %>% 
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d_gapfilled = cumsum(NEE_uStar_f*1800*12.01/1000000)) %>% 
  mutate(ch4_sum_g_m2_d_gapfilled = cumsum(ch4_flux_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v))) %>% 
  mutate(NumCount = 1:length(ec2_annual_fluxes_year1$DateTime))
  
ec2_annual_fluxes_year2 <- ec2_annual_fluxes %>% 
  filter(DateTime >= as.POSIXct("2021-05-01")) %>% 
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d_gapfilled = cumsum(NEE_uStar_f*1800*12.01/1000000)) %>% 
  mutate(ch4_sum_g_m2_d_gapfilled = cumsum(ch4_flux_uStar_f*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v))) %>% 
  mutate(NumCount = 1:length(ec2_annual_fluxes_year2$DateTime))

## Non gap-filled
## Fig. S12
co2_annual <- ggplot()+
  geom_vline(xintercept = 8854,linetype="dotted",color="#ED6E78",size=0.8)+ #Turnover FCR (11-01-2022); operationally defined
  geom_vline(xintercept = 8955, linetype="dotted",color="#A31420",size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=co2_sum_g_m2_d,color="Year1"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=co2_sum_g_m2_d,ymin=co2_sum_g_m2_d-co2_sum_sd,ymax=co2_sum_g_m2_d+co2_sum_sd,fill="Year1"),alpha=0.3)+
  geom_line(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=co2_sum_g_m2_d,color="Year2"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=co2_sum_g_m2_d,ymin=co2_sum_g_m2_d-co2_sum_sd,ymax=co2_sum_g_m2_d+co2_sum_sd,fill="Year2"),alpha=0.3)+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_fill_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_x_continuous(breaks=c(0,2929,5905,8833,11763,14593,17518), 
                   labels=c("May","Jul.","Sep.","Nov.","Jan.","Mar.","May"))+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  xlab("")+
  theme_classic(base_size = 15)

ch4_annual <- ggplot()+
  geom_vline(xintercept = 8854,linetype="dotted",color="#ED6E78",size=0.8)+ #Turnover FCR (11-01-2022); operationally defined
  geom_vline(xintercept = 8955, linetype="dotted",color="#A31420",size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d,color="Year1"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d,ymin=ch4_sum_g_m2_d-ch4_sum_sd,ymax=ch4_sum_g_m2_d+ch4_sum_sd,fill="Year1"),alpha=0.3)+
  geom_line(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d,color="Year2"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d,ymin=ch4_sum_g_m2_d-ch4_sum_sd,ymax=ch4_sum_g_m2_d+ch4_sum_sd,fill="Year2"),alpha=0.3)+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_fill_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_x_continuous(breaks=c(0,2929,5905,8833,11763,14593,17518), 
                     labels=c("May","Jul.","Sep.","Nov.","Jan.","Mar.","May"))+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  xlab("")+
  theme_classic(base_size = 15)

ggarrange(co2_annual,ch4_annual,nrow=1,ncol=2,common.legend = TRUE,labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Rev_AnnualFluxes.png",width = 9, height=4.5, units="in",dpi=320)

## Plot gap-filled annual estimates
## Fig. 5
co2_sum_gapfilled <- ggplot()+
  geom_vline(xintercept = 8854,linetype="dotted",color="#ED6E78",size=0.8)+ #Turnover FCR (11-01-2022); operationally defined
  geom_vline(xintercept = 8955, linetype="dotted",color="#A31420",size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=co2_sum_g_m2_d_gapfilled,color="Year1"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=co2_sum_g_m2_d_gapfilled,ymin=co2_sum_g_m2_d_gapfilled-co2_sum_sd,ymax=co2_sum_g_m2_d_gapfilled+co2_sum_sd,fill="Year1"),alpha=0.3)+
  geom_line(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=co2_sum_g_m2_d_gapfilled,color="Year2"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=co2_sum_g_m2_d_gapfilled,ymin=co2_sum_g_m2_d_gapfilled-co2_sum_sd,ymax=co2_sum_g_m2_d_gapfilled+co2_sum_sd,fill="Year2"),alpha=0.3)+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_fill_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_x_continuous(breaks=c(0,2929,5905,8833,11763,14593,17518), 
                     labels=c("May","Jul.","Sep.","Nov.","Jan.","Mar.","May"))+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  xlab("")+
  theme_classic(base_size = 15)

ch4_sum_gapfilled <- ggplot()+
  geom_vline(xintercept = 8854,linetype="dotted",color="#ED6E78",size=0.8)+ #Turnover FCR (11-01-2022); operationally defined
  geom_vline(xintercept = 8955, linetype="dotted",color="#A31420",size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d_gapfilled,color="Year1"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year1,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d_gapfilled,ymin=ch4_sum_g_m2_d_gapfilled-ch4_sum_sd,ymax=ch4_sum_g_m2_d_gapfilled+ch4_sum_sd,fill="Year1"),alpha=0.3)+
  geom_line(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d_gapfilled,color="Year2"),size=1)+
  geom_ribbon(ec2_annual_fluxes_year2,mapping=aes(x=NumCount,y=ch4_sum_g_m2_d_gapfilled,ymin=ch4_sum_g_m2_d_gapfilled-ch4_sum_sd,ymax=ch4_sum_g_m2_d_gapfilled+ch4_sum_sd,fill="Year2"),alpha=0.3)+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_fill_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  scale_x_continuous(breaks=c(0,2929,5905,8833,11763,14593,17518), 
                     labels=c("May","Jul.","Sep.","Nov.","Jan.","Mar.","May"))+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  xlab("")+
  theme_classic(base_size = 15)

ggarrange(co2_sum_gapfilled,ch4_sum_gapfilled,nrow=1,ncol=2,common.legend = TRUE,labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/SI_Rev_AnnualFluxes_GapFilled.png",width = 9, height=4.5, units="in",dpi=320)


## Calculate cumulative fluxes in summer (Jun, Jul, Aug, Sep) vs. winter (Dec, Jan, Feb, Mar)
## Summer (May - Oct) = 132 g C-CO2 m2 d; 0.221 g C-CH4 m2 d; Gap-filled: 0.91 CH4; 474 CO2
ec2_annual_fluxes_year1 %>% 
  filter(DateTime == as.POSIXct("2020-05-01 01:00:00") | DateTime == as.POSIXct("2020-10-31 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

## Summer (May - Oct) = 124 g C-CO2 m2 d; 0.234 g C-CH4 m2 d; Gap-filled: 619 CO2; 1.21 CH4
ec2_annual_fluxes_year2 %>% 
  filter(DateTime == as.POSIXct("2021-05-01 01:00:00") | DateTime == as.POSIXct("2021-10-31 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled,ch4_sum_g_m2_d,co2_sum_g_m2_d)

# Winter (Nov - Apr) = 38.1 g C-CO2 m2 d; 0.012 g C-CH4 m2 d; Gap-filled: 159 CO2; 0.11 CH4
ec2_annual_fluxes_year1 %>% 
  filter(DateTime == as.POSIXct("2020-10-31 24:00:00") | DateTime == as.POSIXct("2021-04-30 23:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

# Winter (Nov - Apr) = 33.3 g C-CO2 m2 d; 0.003 g C-CH4 m2 d; Gap-filled: 111 CO2; 0.08 CH4
ec2_annual_fluxes_year2 %>% 
  filter(DateTime == as.POSIXct("2021-10-31 24:00:00") | DateTime == as.POSIXct("2022-04-30 23:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

## Plot Summer vs.Winter
Season <- c("Summer","Winter","Summer","Winter")
year <- c("Year 1", "Year 1", "Year 2", "Year 2")
co2 <- c(474, 159, 619, 111)
ch4 <- c(0.91, 0.11, 1.21, 0.08)

ghg_comps <- data.frame(Season,year,co2,ch4)

co2_comp <- ggplot(ghg_comps,mapping=aes(x=year,y=co2,fill=Season))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(breaks=c('Summer','Winter'),values=c("#ED6E78","#A31420"),"")+
  xlab("")+
  ylab(expression(paste("CO"[2]*" (g C m"^-2*")")))+
  theme_classic(base_size = 15)

ch4_comp <- ggplot(ghg_comps,mapping=aes(x=year,y=ch4,fill=Season))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(breaks=c('Summer','Winter'),values=c("#ED6E78","#A31420"),"")+
  xlab("")+
  ylab(expression(paste("CH"[4]*" (g C m"^-2*")")))+
  theme_classic(base_size = 15)

ggarrange(co2_comp,ch4_comp,nrow=1,ncol=2,common.legend = TRUE,labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/SI_GHG_Annual_Season.png",width = 9, height=4.5, units="in",dpi=320)

# June-September
ec2_annual_fluxes_year1 %>% 
  filter(DateTime == as.POSIXct("2020-06-01 01:00:00") | DateTime == as.POSIXct("2020-09-30 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

ec2_annual_fluxes_year2 %>% 
  filter(DateTime == as.POSIXct("2021-06-01 01:00:00") | DateTime == as.POSIXct("2021-09-30 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

# December - March
ec2_annual_fluxes_year1 %>% 
  filter(DateTime == as.POSIXct("2020-12-01 01:00:00") | DateTime == as.POSIXct("2021-03-31 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)

ec2_annual_fluxes_year2 %>% 
  filter(DateTime == as.POSIXct("2021-12-01 01:00:00") | DateTime == as.POSIXct("2022-03-31 24:00:00")) %>% 
  select(DateTime,ch4_sum_g_m2_d,co2_sum_g_m2_d,ch4_sum_g_m2_d_gapfilled,co2_sum_g_m2_d_gapfilled)
