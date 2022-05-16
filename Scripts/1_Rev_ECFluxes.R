### Script conduct analysis of EC fluxes from FCR
### Apr 2020 - Apr 2022

### Following revisions from inital submission to JGR-Biogeosciences
### 12 May 2022, A. Hounshell

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales,lubridate,
               lognorm,MuMIn,rsq,Metrics,astsa,DescTools,kSamples,openair)

###############################################################################

## First load in QA/QC'd EC fluxes 

## Load in data from Brenda - 30 minute fluxes from 2020-04-04 to 2021-05-06
## Data corrected following FCR_Process_BD
## Data downloaded from: https://doi.org/10.6073/pasta/a1324bcf3e1415268996ba867c636489

ec <- read.csv("./Data/20220506_EC_processed.csv") 

ec2 <- ec %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00"))

###############################################################################

## Load in 'raw' EC data following EddyPro processing
ec_raw <- read.csv("./Data/20220506_EddyPro_cleaned.csv")

# Format time
ec_raw$datetime <- as.POSIXct(paste(ec_raw$date, ec$time), format="%Y-%m-%d %H:%M", tz="EST")
ec_raw$datetime <- as_datetime(ec_raw$datetime)

# Set new dataframe with list of dates+times:
# every 30 minutes
# Constrain to study time period: 2020-04-05 (time series start date) to
# last 30 minute period - UPDATE WITH EACH NEW SET OF DATA!
ts <- seq.POSIXt(as.POSIXct("2020-05-01 00:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), 
                 as.POSIXct("2022-04-30 24:00:00",'%Y-%m-%d %H:%M:%S', tz="EST"), by = "30 min")
ts2 <- data.frame(datetime = ts)

# Join Eddy Flux data with list of dates+time
ec2_raw <- left_join(ts2, ec_raw, by = 'datetime')

###############################################################################

## Visualizations of missing data - following all QA/QC

## Amount of data retained following filtering due to low u*
ec2 %>% select(DateTime, NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
  summarise(co2_available = 100-sum(is.na(NEE_uStar_orig))/n()*100,
            ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)
# 23% CO2 fluxes; 19% CH4 fluxes

## Distribution of wind speed and direction were removed?

### COME BACK TO THIS!!!!

## Distribution of missing data by season, day vs. night, etc.
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

season_data_total <- season_data %>% 
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

ggsave("./Fig_Output/SI_Data_Season.jpg",width = 8, height=4, units="in",dpi=320)

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

## Plot daily for both years

co2_daily_year1 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="EC"),size = 1) +
  #geom_errorbar(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  #geom_point(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  #scale_color_manual(breaks=c("EC","Diff"),
  #                   values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

co2_daily_year2 <- ggplot(fcr_daily) +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, NEE,color="EC"),size = 1) +
  #geom_errorbar(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  #geom_point(ghg_fluxes_avg_co2,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  #scale_color_manual(breaks=c("EC","Diff"),
  #                   values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_daily_year1 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="EC"),size = 1) +
  #geom_errorbar(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  #geom_point(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  #scale_color_manual(breaks=c("EC","Diff"),
  #                   values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2020-05-01"),as.POSIXct("2021-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ch4_daily_year2 <- fcr_daily %>% 
  ggplot() +
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
  geom_point(ec2,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha = 0.1)+
  geom_ribbon(mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(aes(Date, CH4,color="EC"),size = 1) +
  #geom_errorbar(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,ymin=mean-sd,ymax=mean+sd,color="Diff"),size=1)+
  #geom_point(ghg_fluxes_avg_ch4,mapping=aes(x=DateTime,y=mean,color="Diff"),size=2)+
  #scale_color_manual(breaks=c("EC","Diff"),
  #                   values=c("#E63946","#4c8bfe"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  xlim(as.POSIXct("2021-05-01"),as.POSIXct("2022-04-30"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(co2_daily_year1,co2_daily_year2,ch4_daily_year1,ch4_daily_year2,ncol=1,
          nrow=4,labels=c("A.","B.","C.","D."),font.label = list(face="plain",size=15), common.legend = TRUE)

ggsave("./Fig_Output/EC_Daily.jpg",width = 9, height=12, units="in",dpi=320)

###############################################################################

## Plot weekly for both years

co2_week <- ggplot(fcr_weekly)+
  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1,linetype="dotted") + 
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

## Plot monthly for all years
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

## Plot by half-hourly fluxes by season
# Define as: Winter = Dec, Jan, Feb; Spring = Mar, Apr, May; Summer = Jun, Jul, Aug;
# Fall = Sep, Oct, Nov
ec2_season <- ec2 %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig) %>% 
  mutate(Month = month(DateTime),
         Year = ifelse(DateTime < as.POSIXct("2021-05-01"), "Year1", "Year2")) %>% 
  mutate(Season = ifelse(Month == 12 | Month == 1 | Month == 2, "Winter",
                         ifelse(Month == 3 | Month == 4 | Month == 5, "Spring",
                                ifelse(Month == 6 | Month == 7 | Month == 8, "Summer",
                                       ifelse(Month == 9 | Month == 10 | Month == 11, "Fall", NA)))))

# Compare seasons by Year1 and Year2
ec2_fall <- ec2_season %>% 
  filter(Season == "Fall")

ec2_winter <- ec2_season %>% 
  filter(Season == "Winter")

ec2_spring <- ec2_season %>% 
  filter(Season == "Spring")

ec2_summer <- ec2_season %>% 
  filter(Season == "Summer")

# Calculate statistics between seasons of each year (Year 1 vs. Year 2)
wilcox.test(NEE_uStar_orig ~ Year, data=ec2_fall)
wilcox.test(NEE_uStar_orig ~ Year, data=ec2_winter)
wilcox.test(NEE_uStar_orig ~ Year, data=ec2_spring)
wilcox.test(NEE_uStar_orig ~ Year, data=ec2_summer)

wilcox.test(ch4_flux_uStar_orig ~ Year, data=ec2_fall)
wilcox.test(ch4_flux_uStar_orig ~ Year, data=ec2_winter)
wilcox.test(ch4_flux_uStar_orig ~ Year, data=ec2_spring)
wilcox.test(ch4_flux_uStar_orig ~ Year, data=ec2_summer)

co2_season <- ggplot(ec2_season,mapping=aes(x=Season,y=NEE_uStar_orig,color=as.factor(Year)))+
  annotate("text",label ="*",x = "Fall", hjust = 0, y = 18, size = 5)+
  annotate("text",label ="***",x="Winter",hjust=0,y=18,size=5)+
  annotate("text",label="***",x="Spring",hjust=0,y=18,size=5)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  geom_boxplot(outlier.shape = NA,size=1,alpha=0.3)+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  xlab("")+
  ylim(-20,20)+
  theme_classic(base_size = 15)

ch4_season <- ggplot(ec2_season,mapping=aes(x=Season,y=ch4_flux_uStar_orig,color=as.factor(Year)))+
  annotate("text",label ="***",x = "Fall", hjust = 0, y = 0.045, size = 5)+
  annotate("text",label ="***",x="Winter",hjust=0,y=0.045,size=5)+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  geom_boxplot(outlier.shape = NA,size=1,alpha=0.3)+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  scale_color_manual(breaks=c('Year1','Year2'),values=c("#ED6E78","#A31420"),labels=c("Year 1","Year 2"),"")+
  xlab("")+
  ylim(-0.05,0.05)+
  theme_classic(base_size = 15)

ggarrange(co2_month,ch4_month,
          ggarrange(co2_season, ch4_season,nrow=1,ncol=2, labels = c("C.","D."), font.label = list(face="plain",size=15), common.legend = TRUE),
          nrow=3,ncol=1,labels = c("A.","B."),
          font.label = list(face="plain",size=15),common.legend = TRUE)

ggsave("./Fig_Output/EC_Monthly.jpg",width = 8, height = 11, units="in",dpi=320)

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

diel_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

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

dawn_agg %>% 
  select(diel,season,NEE,CH4,u) %>% 
  group_by(diel,season) %>% 
  summarise_all(median,na.rm=TRUE)

## Calculate statistics for Supplementary Table (Fig S7) - Day vs. Night
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

## Calculate statistics for Supplementary Table (Fig S7) - Dawn vs. Dusk
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

## Plot diel (day/night) and dawn/dusk comparisons: Figure XX
diel_co2 <- diel_agg %>% 
  ggplot(mapping=aes(x=diel,y=NEE,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label ="p = 0.09
n = 295",
           x = "Day", hjust = 0, y = 18, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  ylim(-20,20)+
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
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
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
           x = "Dawn", hjust = 0, y = 18, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Dusk','Dawn'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  ylim(-20,20)+
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
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
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

ggsave("./Fig_Output/Rev_Figure3.png",width = 9, height=7.5, units="in",dpi=320)

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

wilcox.test(NEE_uStar_orig ~ Ice, data=ice_fluxes_30min_comps)
wilcox.test(ch4_flux_uStar_orig ~ Ice, data=ice_fluxes_30min_comps)

## Plot ice on/ice off for 2021 and 2022
co2_year1 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=NEE),color="#E63946",size=1)+
  xlim(as.POSIXct("2021-01-01"),as.POSIXct("2021-02-28"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

ch4_year1 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2021-01-10"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2021-02-09"), linetype = "dotted", color="red")+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=ch4_flux_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=CH4,ymin=CH4-CH4_sd,ymax=CH4+CH4_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=CH4),color="#E63946",size=1)+
  xlim(as.POSIXct("2021-01-01"),as.POSIXct("2021-02-28"))+
  xlab("") +
  ylab(expression(~CH[4]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

co2_year2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2022-01-16"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), linetype = "dotted", color="red")+
  geom_point(ice_fluxes_30min,mapping=aes(x=DateTime,y=NEE_uStar_orig),alpha=0.1)+
  geom_ribbon(ice_fluxes_daily,mapping=aes(x=Date,y=NEE,ymin=NEE-NEE_sd,ymax=NEE+NEE_sd),fill="#E63946",alpha=0.5)+
  geom_line(ice_fluxes_daily,mapping=aes(x=Date,y=NEE),color="#E63946",size=1)+
  xlim(as.POSIXct("2022-01-01"),as.POSIXct("2022-02-28"))+
  xlab("") +
  ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic(base_size = 15)

ch4_year2 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2022-01-16"), linetype = "dotted", color="blue")+
  geom_vline(xintercept = as.POSIXct("2022-02-10"), linetype = "dotted", color="red")+
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

ggsave("./Fig_Output/Figure_WinterIce.png",width = 8, height=9, units="in",dpi=320)

###############################################################################

## Estimate and plot annual fluxes for Year 1 and Year 2, CO2 and CH4
ec2_annual_fluxes <- ec2 %>% 
  select(DateTime,NEE_uStar_orig,ch4_flux_uStar_orig,NEE_uStar_fsd,ch4_flux_uStar_fsd) %>% 
  mutate(NEE_uStar_orig = ifelse(is.na(NEE_uStar_orig), 0, NEE_uStar_orig),
         ch4_flux_uStar_orig = ifelse(is.na(ch4_flux_uStar_orig), 0, ch4_flux_uStar_orig))

ec2_annual_fluxes_year1 <- ec2_annual_fluxes %>% 
  filter(DateTime < as.POSIXct("2021-05-01")) %>% 
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v))) %>% 
  mutate(NumCount = 1:length(ec2_annual_fluxes_year1$DateTime))
  
ec2_annual_fluxes_year2 <- ec2_annual_fluxes %>% 
  filter(DateTime >= as.POSIXct("2021-05-01")) %>% 
  mutate(ch4_sum_g_m2_d = cumsum(ch4_flux_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_sum_g_m2_d = cumsum(NEE_uStar_orig*1800*12.01/1000000)) %>% 
  mutate(co2_v = (NEE_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(ch4_v = (ch4_flux_uStar_fsd*1800*12.01/1000000)^2) %>% 
  mutate(co2_sum_sd = sqrt(cumsum(co2_v))) %>% 
  mutate(ch4_sum_sd = sqrt(cumsum(ch4_v))) %>% 
  mutate(NumCount = 1:length(ec2_annual_fluxes_year2$DateTime))

co2_annual <- ggplot()+
  geom_vline(xintercept = 8854,linetype="dotted",color="#ED6E78",size=0.8)+ #Turnover FCR (11-01-2022); operationally defined
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
