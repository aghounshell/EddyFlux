# Master QAQC script in prep for publishing catwalk sensor string data to EDI
# this script combines code from other QAQC scripts found in misc_QAQC_scipts folder
# as well as other data files from misc_data_files
# final EDI-ready file outputs directly to MakeEMLCatwalk/2020 folder
# Set up ----
pacman::p_load("RCurl","tidyverse","lubridate", "plotly", "magrittr")

# Set working directory
wd <- getwd()
setwd(wd)

folder <- "./Scripts/Old_Scripts/"
source(paste0(folder, "FCR_catwalk_QAQC_function_2018_2021.R"))

# download most up to date catwalk data and maintenance log
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data/CAT_MaintenanceLog.txt",paste0(wd, "./Data/FCR_CAT_MaintenanceLog_2018_2021.txt"))
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data/Catwalk.csv",paste0(wd, "./Data/Catwalk.csv"))
download.file('https://raw.githubusercontent.com/CareyLabVT/ManualDownloadsSCCData/master/CR6_Files/FCRcatwalk_manual_2021.csv', paste0(wd, "./Data/CAT_2.csv"))

# run standard qaqc function from FCR_catwalk_QAQC_function_2021.R
data_file <- paste0(wd, './Data/Catwalk.csv')#current file from the data logger
data2_file <- paste0(wd, './Data/CAT_2.csv')#manual downloads to add missing data 
maintenance_file <- paste0(wd, "./Data/FCR_CAT_MaintenanceLog_2018_2021.txt")#maintenance file
output_file <- paste0(wd, "./Data/Catwalk_first_QAQC_2018_2021.csv")#name of the output file
qaqc(data_file,data2_file, maintenance_file, output_file)#function to do the main qaqc


# read in qaqc function output

catdata <- read.csv(output_file)

#current time of QAQC for graphing
start_time="2021-01-01 00:00"
end_time="2021-12-31 23:59"


# subset file to only unpublished data

catdata_flag =catdata%>%
  filter(DateTime<"2022-01-01 00:00")

#make sure no NAS in the Flag columns
Flags=catdata_flag%>%
  select(DateTime, starts_with("Flag"))

RowsNA=Flags[!complete.cases(Flags), ] # Keep only the complete rows

# Flag values
# 0: no flag
# 1: value removed due to maintenance and set to NA
# 2: negative or outlier value removed and set to NA, see Methods section for more detail on QAQC process
# 3: negative values set to 0
# 4: value removed due to fouling and set to NA
# 5: questionable value due to potential fouling
# 6: very questionable value due to potential fouling. Values adjusted using a linear or square root function to match high-resolution CTD profiles are given in RDO_mgL_5 and RDO_sat_percent_5
# 7: missing data
# 8: Value corrected using a constant offset due to two thermistor malfunctions in Fall 2020


###########################################################################################################################################################################
# temp qaqc ----

#Two of the thermistors started to read higher than the one above them in fall 2020. Fixed this using a constant offset. 
#methods described in metadata

#start time for 1m is 30 Oct 2020 13:00EST
#start time for 4m is 31 Oct 2020 5:00EST
catdata_flag <- catdata_flag %>%
  mutate(Flag_Temp_1 = ifelse(DateTime >= "2020-10-30 13:00" & DateTime < "2020-12-31 23:50" & (! is.na(ThermistorTemp_C_1)) ,8, Flag_Temp_1))%>%
  mutate(Flag_Temp_4 = ifelse(DateTime >= "2020-10-31 5:00" & DateTime < "2020-12-31 23:50" &
                                (! is.na(ThermistorTemp_C_4)),8, Flag_Temp_4))%>%
  mutate(ThermistorTemp_C_1 = ifelse(DateTime >= "2020-10-30 13:00" & DateTime < "2020-12-31 23:50" &
                                        (! is.na(ThermistorTemp_C_1)), (ThermistorTemp_C_1-0.22617), ThermistorTemp_C_1 )) %>%
  mutate(ThermistorTemp_C_4 = ifelse(DateTime >= "2020-10-31 5:00" & DateTime < "2020-12-31 23:50" &
                                        (! is.na(ThermistorTemp_C_4)), (ThermistorTemp_C_4-0.18122), ThermistorTemp_C_4 )) 
#thermistors were replaced in 2021
catdata_flag <- catdata_flag %>%
  mutate(Flag_Temp_1 = ifelse(DateTime >= "2021-02-08 14:30" & DateTime < "2021-02-26 12:00" & (! is.na(ThermistorTemp_C_2)) ,7, Flag_Temp_1))%>%
  mutate(Flag_Temp_4 = ifelse(DateTime >= "2021-02-08 14:30" & DateTime < "2021-02-26 12:00" &
                                (! is.na(ThermistorTemp_C_2)),7, Flag_Temp_4))

#Checktime=catdata_flag%>%
#  select(c(DateTime, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_4, Flag_Temp_1, Flag_Temp_4))
# check surface temp data



################
#Correcting some DO values 

catdata_flag <- catdata_flag%>%
  mutate(
    #DateTime=as.character(DateTime),#change DateTime to as.character so they line up when making changes
    #For 5m READ NOTE: These are the sections I noticed apparent following in TS and have tried to correct with linear adjustments
    Flag_DO_5_obs = ifelse(DateTime<"2019-08-12 12:35:00" & DateTime>"2019-08-11 00:00:00", 6, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime<"2019-08-12 12:35:00" & DateTime> "2019-08-11 00:00:00", 6, Flag_DO_5_sat),
    RDO_mgL_5_adjusted = ifelse(DateTime<"2019-08-12 12:35:00" & DateTime> "2019-08-11 00:00:00", 
                                RDO_mgL_5 + sqrt(as.numeric(difftime(DateTime,"2019-08-11 00:00:00", units = "mins")))/30,
                                RDO_mgL_5),
    RDOsat_percent_5_adjusted = ifelse(DateTime<"2019-08-12 12:35:00" & DateTime> "2019-08-11 00:00:00", 
                                       RDOsat_percent_5 + sqrt(as.numeric(difftime(DateTime,"2019-08-11 00:00:00", units = "mins")))/30/11.3*100,
                                       RDOsat_percent_5),
    Flag_DO_5_obs = ifelse(DateTime<"2019-07-17 11:40:00" & DateTime>"2019-07-13 00:00:00", 5, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime<"2019-07-17 11:40:00" & DateTime>"2019-07-13 00:00:00", 5, Flag_DO_5_sat),
    #9 meters #### READ NOTE: These are the sections I noticed apparent following in TS and have tried to correct with linear adjustments
    #first section fixed 11aug to 17aug
    Flag_DO_9_obs = ifelse(DateTime<"2020-08-17 12:50:00" & DateTime> "2020-08-11 7:00:00",6, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime<"2020-08-17 12:50:00" & DateTime> "2020-08-11 15:00:00",6, Flag_DO_9_sat),
    RDO_mgL_9_adjusted = ifelse(DateTime<"2020-08-17 12:50:00" & DateTime> "2020-08-11 7:00:00", 
                                RDO_mgL_9 + as.numeric(difftime(DateTime,"2020-08-11 7:00:00", units = "mins"))/6500, #A linear adjustment here
                                RDO_mgL_9),
    RDOsat_percent_9_adjusted = ifelse(DateTime<"2020-08-17 12:50:00" & DateTime> "2020-08-11 15:00:00", 
                                       RDOsat_percent_9 + (as.numeric(difftime("2020-08-17 12:50:00","2020-08-12 15:00:00", units = "mins")))/6500/11.3*100,
                                       RDOsat_percent_9),
    
    # 9 meters 19aug to 24aug
    Flag_DO_9_obs = ifelse(DateTime<"2020-08-24 10:40:00" & DateTime> "2020-08-19 20:00:00",6, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime<"2020-08-24 10:40:00" & DateTime> "2020-08-19 20:00:00",6, Flag_DO_9_sat),
    RDO_mgL_9_adjusted = ifelse(DateTime<"2020-08-24 10:40:00" & DateTime> "2020-08-19 20:00:00", 
                                RDO_mgL_9 + (as.numeric(difftime(DateTime,"2020-08-19 20:00:00", units = "mins")))/6500, #A linear adjustment here
                                RDO_mgL_9_adjusted),
    RDOsat_percent_9_adjusted = ifelse(DateTime<"2020-08-24 10:40:00" & DateTime> "2020-08-19 20:00:00", 
                                       RDOsat_percent_9 + as.numeric(difftime(DateTime,"2020-08-19 20:00:00", units = "mins"))/6500/11.3*100,
                                       RDOsat_percent_9_adjusted),
    
    # 9 meters 26aug to 02sep
    Flag_DO_9_obs = ifelse(DateTime<"2020-09-02 10:50:00" & DateTime> "2020-08-26 12:00:00" ,6, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime<"2020-09-02 10:50:00" & DateTime> "2020-08-26 12:00:00" ,6, Flag_DO_9_sat),
    RDO_mgL_9_adjusted = ifelse(DateTime< "2020-09-02 10:50:00" & DateTime> "2020-08-26 12:00:00", 
                                RDO_mgL_9 + (as.numeric(difftime(DateTime, "2020-08-26 12:00:00", units = "mins")))/10000, #A linear adjustment here
                                RDO_mgL_9_adjusted),
    RDOsat_percent_9_adjusted = ifelse(DateTime< "2020-09-02 10:50:00" & DateTime> "2020-08-26 12:00:00", 
                                       RDOsat_percent_9 + as.numeric(difftime(DateTime, "2020-08-26 12:00:00", units = "mins"))/10000/11.3*100,
                                       RDOsat_percent_9_adjusted),
    
    # 9 meters 02sep to 11sep
    Flag_DO_9_obs = ifelse(DateTime<"2020-09-09 17:50:00" & DateTime> "2020-09-05 06:00:00" ,6, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime<"2020-09-09 17:50:00" & DateTime> "2020-09-05 06:00:00" ,6, Flag_DO_9_sat),
    RDO_mgL_9_adjusted = ifelse(DateTime< "2020-09-09 17:50:00" & DateTime> "2020-09-05 06:00:00", 
                                RDO_mgL_9 + (as.numeric(difftime(DateTime, "2020-09-05 06:00:00", units = "mins")))/3000, #A linear adjustment here
                                RDO_mgL_9_adjusted),
    RDOsat_percent_9_adjusted = ifelse(DateTime< "2020-09-09 17:00:00" & DateTime> "2020-09-05 06:00:00", 
                                       RDOsat_percent_9 + as.numeric(difftime(DateTime, "2020-09-05 06:00:00", units = "mins"))/3000/11.3*100,
                                       RDOsat_percent_9_adjusted))

#take out questionable values for the 5m DO from 2018-2021
catdata_flag <- catdata_flag %>%
  mutate(
    Flag_DO_5_obs = ifelse(DateTime >= "2018-09-25 10:40" & DateTime<"2018-09-25 10:55",2, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime >= "2018-09-25 10:40" & DateTime<"2018-09-25 10:55",2, Flag_DO_5_sat),
    RDOsat_percent_5 = ifelse(DateTime >= "2018-09-25 10:40"& DateTime<"2018-09-25 10:55",NA, RDOsat_percent_5), 
    RDO_mgL_5 = ifelse(DateTime >= "2018-09-25 10:40"& DateTime<"2018-09-25 10:55",NA, RDO_mgL_5),
    Flag_DO_5_obs = ifelse(DateTime >= "2019-06-17 13:50" & DateTime < "2019-06-17 15:40" & 
                             (! is.na(RDO_mgL_5)) ,2, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime >= "2019-06-17 13:50" & DateTime < "2019-06-17 15:40" & 
                             (! is.na(RDOsat_percent_5)) ,2, Flag_DO_5_sat),
    RDOsat_percent_5 = ifelse(DateTime >= "2019-06-17 13:50" & DateTime < "2019-06-17 15:40" & 
                                (! is.na(RDOsat_percent_5)) ,NA, RDOsat_percent_5),
    RDO_mgL_5 = ifelse(DateTime >= "2019-06-17 13:50" & DateTime < "2019-06-17 15:40" &
                         (! is.na(RDO_mgL_5)),NA, RDO_mgL_5),
    Flag_DO_5_sat = ifelse(DateTime >= "2020-07-06 14:10" & DateTime < "2020-07-06 14:40" & 
                             (! is.na(RDOsat_percent_5)) ,2, Flag_DO_5_sat),
    Flag_DO_5_obs = ifelse(DateTime >= "2020-07-06 14:10" & DateTime < "2020-07-06 14:40" & 
                             (! is.na(RDO_mgL_5)) ,2, Flag_DO_5_obs),
    RDOsat_percent_5 = ifelse(DateTime >= "2020-07-06 14:10" & DateTime < "2020-07-06 14:40" & 
                                (! is.na(RDOsat_percent_5)) ,NA, RDOsat_percent_5), 
    RDO_mgL_5 = ifelse(DateTime >= "2020-07-06 14:10" & DateTime < "2020-07-06 14:40" &
                         (! is.na(RDO_mgL_5)),NA, RDO_mgL_5),
    Flag_DO_5_obs = ifelse(DateTime >= "2020-07-13 13:00" & DateTime < "2020-07-13 14:50" & 
                             (! is.na(RDO_mgL_5)) ,2, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime >= "2020-07-13 13:00" & DateTime < "2020-07-13 14:50" & 
                             (! is.na(RDOsat_percent_5)) ,2, Flag_DO_5_sat),
    RDOsat_percent_5 = ifelse(DateTime >= "2020-07-13 13:00" & DateTime < "2020-07-13 14:50" & 
                                (! is.na(RDOsat_percent_5)) ,NA, RDOsat_percent_5), 
    RDO_mgL_5 = ifelse(DateTime >= "2020-07-13 13:00" & DateTime < "2020-07-13 14:50" &
                         (! is.na(RDO_mgL_5)),NA, RDO_mgL_5),
    Flag_DO_5_obs = ifelse(DateTime >= "2021-06-11 9:41" & DateTime < "2021-09-14 13:30" & 
                             (! is.na(RDO_mgL_5)) ,5, Flag_DO_5_obs),
    Flag_DO_5_sat = ifelse(DateTime >= "2021-06-11 9:41" & DateTime < "2021-09-14 13:30" & 
                             (! is.na(RDOsat_percent_5)) ,5, Flag_DO_5_sat))


#take out questionable values for the 9m DO from 2018-2021
catdata_flag <- catdata_flag %>%
  mutate(
    Flag_DO_9_obs = ifelse(DateTime >= "2018-11-19 15:00" & DateTime < "2018-11-20 12:20" & 
                             (! is.na(RDO_mgL_9)) ,2, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime >= "2018-11-19 15:00" & DateTime < "2018-11-20 12:20" & 
                             (! is.na(RDOsat_percent_9)) ,2, Flag_DO_9_sat),
    RDOsat_percent_9 = ifelse(DateTime >= "2018-11-19 15:00" & DateTime < "2018-11-20 12:20" & 
                                (! is.na(RDOsat_percent_9)),NA, RDOsat_percent_9), 
    RDO_mgL_9 = ifelse(DateTime >= "2018-11-19 15:00" & DateTime < "2018-11-20 12:20" &
                         (! is.na(RDO_mgL_9)),NA, RDO_mgL_9),
    Flag_DO_9_obs = ifelse(DateTime >= "2020-12-02 15:00" & DateTime < "2020-12-02 18:00" & 
                             (! is.na(RDO_mgL_9)) ,2, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime >= "2020-12-02 15:00" & DateTime < "2020-12-02 18:00" & 
                             (! is.na(RDOsat_percent_9)) ,2, Flag_DO_9_sat),
    RDOsat_percent_9 = ifelse(DateTime >= "2020-12-02 15:00" & DateTime < "2020-12-02 18:00" & 
                                (! is.na(RDOsat_percent_9)) ,NA, RDOsat_percent_9),
    RDO_mgL_9 = ifelse(DateTime >= "2020-12-02 15:00" & DateTime < "2020-12-02 18:00" &
                         (! is.na(RDO_mgL_9)),NA, RDO_mgL_9),
    Flag_DO_9_obs = ifelse(DateTime >= "2021-02-08 17:40" & DateTime < "2021-02-08 18:10" & 
                             (! is.na(RDO_mgL_9)) ,2, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime >= "2021-02-08 17:40" & DateTime < "2021-02-08 18:10" & 
                             (! is.na(RDOsat_percent_9)) ,2, Flag_DO_9_sat),
    RDOsat_percent_9 = ifelse(DateTime >= "2021-02-08 17:40" & DateTime < "2021-02-08 18:10" & 
                                (! is.na(RDOsat_percent_9)) ,NA, RDOsat_percent_9),
    RDO_mgL_9 = ifelse(DateTime >= "2021-02-08 17:40" & DateTime < "2021-02-08 18:10" &
                         (! is.na(RDO_mgL_9)),NA, RDO_mgL_9),
    Flag_DO_9_obs = ifelse(DateTime >= "2021-06-07 12:00" & DateTime < "2021-06-08 12:00" & 
                             (! is.na(RDO_mgL_9)) ,2, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime >= "2021-06-07 12:00" & DateTime < "2021-06-08 12:00" & 
                             (! is.na(RDOsat_percent_9) & RDOsat_percent_9>0) ,2, Flag_DO_9_sat),
    RDOsat_percent_9 = ifelse(DateTime >= "2021-06-07 12:00" & DateTime < "2021-06-08 12:00" & 
                                (! is.na(RDOsat_percent_9)& RDOsat_percent_9>0) ,NA, RDOsat_percent_9),
    RDO_mgL_9 = ifelse(DateTime >= "2021-06-07 12:40" & DateTime < "2021-06-08 12:00" &
                         (! is.na(RDO_mgL_9) & RDO_mgL_9>0),NA, RDO_mgL_9),
    Flag_DO_9_obs = ifelse(DateTime >= "2021-06-11 9:41" & DateTime < "2021-09-14 13:30" & 
                             (! is.na(RDO_mgL_9)) ,5, Flag_DO_9_obs),
    Flag_DO_9_sat = ifelse(DateTime >= "2021-06-11 9:41" & DateTime < "2021-09-14 13:30" & 
                             (! is.na(RDOsat_percent_9)) ,5, Flag_DO_9_sat))


#Flag high conductivity values in 2020 but don't remove them. If want to remove later I will but I am
# not convinved it is a sensor malfunction

catdata_flag <- catdata_flag %>%
  mutate(
    #DateTime=as.character(DateTime),#change DateTime to as.character so they line up when making changes
    Flag_Cond = ifelse(DateTime >"2020-05-01 00:00" & DateTime < "2020-08-31 23:59" &
                         EXOCond_uScm_1>42 & Flag_Cond==0 ,5, Flag_Cond),
    Flag_SpCond = ifelse(DateTime >"2020-05-01 00:00" & DateTime < "2020-08-31 23:59" &
                           EXOSpCond_uScm_1>42 & Flag_SpCond==0 ,5, Flag_SpCond),
    Flag_TDS = ifelse(DateTime >"2020-05-01 00:00" & DateTime < "2020-08-31 23:59" &
                        EXOTDS_mgL_1>30 & Flag_TDS==0 ,5, Flag_TDS))

#Convert the time and put it in current time zone so the times line up when changing NAs.
#+5 is during EDT and +4 is during EST(make sure to check this in December)
#Have to do strpttime or you get some NAs in the DateTime column
#Do this so the graphing is nice
 catdata_flag$DateTime<-as.POSIXct(strptime(catdata_flag$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5")
 


################################################################################################################
# #graphing temperature
# 
# 
# #Surface Temp
# #From 2018-current
 surf <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_surface)) +
   geom_point()
 surf
# #Plotly so can pick out questionable values
# ggplotly(surf)
# 
# #Just the current year
 surf21=catdata_flag%>%
   filter(DateTime>start_time & DateTime<end_time)%>%
   ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_surface)) +
   geom_point()
 surf21
# #Plotly so can pick out questionable values
# ggplotly(surf21)
# 
# # check 1m temp data
# #From 2018-current
 m_1 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_1)) +
  geom_point()
 m_1
# #Plotly so can pick out questionable values
 ggplotly(m_1)
# 
# #Just the current year
 m_1_21=catdata_flag%>%
   filter(DateTime>start_time & DateTime<end_time)%>%
   ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_1)) +
   geom_point()
 m_1_21
# #Plotly so can pick out questionable values
# ggplotly(m_1_21)
# 
# # check 2m temp data
# #Plot 2018-current
m_2 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_2)) +
 geom_point()
m_2
# #Plotly so can pick out questionable values
# ggplotly(m_2)
# 
# #Just the current year
m_2_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_2)) +
  geom_point()
m_2_21
# #Plotly so can pick out questionable values
 ggplotly(m_2_21)
# 
# # check 3m temp data
# #Plot From 2018-current
m_3 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_3)) +
 geom_point()
m_3
# #Plotly so can pick out questionable values
 ggplotly(m_3)
# 
# #Just the current year
m_3_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_3)) +
  geom_point()
m_3_21
# #Plotly so can pick out questionable values
# ggplotly(m_3_21)
# 
# # check 4m temp data
# #Plot from 2018-current
m_4 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_4)) +
 geom_point()
m_4
# #Plotly so can pick out questionable values
# ggplotly(m_4)
# 
# # Just from current year
m_4_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_4)) +
  geom_point()
m_4_21
# # Plotly so can pick out questionable values
# ggplotly(m_4_21)
# 
# # check 5m temp data
# # Plot from 2018-current
m_5 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_5)) +
  geom_point()
m_5
# # Plotly so can pick out questionable values
# ggplotly(m_5)
# 
# # Just current year
m_5_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_5)) +
  geom_point()
m_5_21
# # Plotly so can pick out questionable values
# ggplotly(m_5_21)
# 
# # check 6m temp data
# # Plot from 2018-current
m_6 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_6)) +
 geom_point()
m_6
# # Plotly so can pick out questionable values
# ggplotly(m_6)
# 
# # Just the current year
m_6_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_6)) +
  geom_point()
m_6_21
# # Plotly so can pick out questionable values
# ggplotly(m_6_21)
# 
# # check 7m temp data
# # all the temp 2018-current
m_7 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_7)) +
 geom_point()
m_7
# # Plotly so can pick out questionable values
# ggplotly(m_7)
# 
# #filter for the current year
m_7_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_7)) +
  geom_point()
m_7_21
# # plotly so you can pick out the questionable values
# ggplotly(m_7_21)
# 
# # check 8m temp data
# # Plot 2018-current
m_8 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_8)) +
 geom_point()
m_8
# # Plotly so can pick out questionable values
# ggplotly(m_8)
# 
# # Plot just the current year
m_8_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_8)) +
  geom_point()
m_8_21
# # Plotly so can pick out questionable values
# ggplotly(m_8_21)
# 
# 
# # check 9m temp data
# # Plot 2018-current
m_9 <- ggplot(data = catdata_flag, aes(x = DateTime, y = ThermistorTemp_C_9)) +
 geom_point()
m_9
# # Plotly so can pick out questionable values
# ggplotly(m_9)
# 
# # Just the current year
m_9_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = ThermistorTemp_C_9)) +
  geom_point()
m_9_21
# # Plotly so can pick out questionable values
# ggplotly(m_9_21)
# 
# 
# 
# #graph all the temps in 2021on the same graph so use base R
# #create a new data frame for 2021
t2021=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)
# 
# #this part taken from the daily email script
par(mfrow=c(1,1))
par(oma=c(1,1,1,4))
plot(t2021$DateTime,t2021$ThermistorTemp_C_surface, main="Water Temp", xlab="Time", ylab="degrees C", type='l', col="firebrick4", lwd=1.5, ylim=c(0,35))
points(t2021$DateTime, t2021$ThermistorTemp_C_1, col="firebrick1", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_2, col="DarkOrange1", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_3, col="gold", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_4, col="greenyellow", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_5, col="medium sea green", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_6, col="sea green", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_7, col="DeepSkyBlue4", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_8, col="blue2", type='l', lwd=1.5)
points(t2021$DateTime, t2021$ThermistorTemp_C_9, col="blue4", type='l', lwd=1.5)
par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=T)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right",c("0.1m","1m", "2m", "3m", "4m", "5m", "6m", "7m","8m", "9m"),
       text.col=c("firebrick4", "firebrick1", "DarkOrange1", "gold", "greenyellow", "medium sea green", "sea green",
                  "DeepSkyBlue4", "blue2", "blue4"), cex=1, y.intersp=1, x.intersp=0.001, inset=c(0,0), xpd=T, bty='n')
# 
# #check the DO temp compared to temp string in 2021
# 
plot(t2021$DateTime,t2021$ThermistorTemp_C_1, main="EXO vs. Temp String", xlab="Time", ylab="degrees C", type='l', col="firebrick4", lwd=1.5, ylim=c(0,35))
points(t2021$DateTime, t2021$ThermistorTemp_C_2, col="firebrick1", type='l', lwd=1.5)
points(t2021$DateTime, t2021$EXOTemp_C_1, col="black", type='l', lwd=1.5)

plot(t2021$DateTime,t2021$RDOTemp_C_5, main="RDO 5m vs. Temp String", xlab="Time", ylab="degrees C", type='l', col="black", lwd=1.5, ylim=c(0,15))
points(t2021$DateTime, t2021$ThermistorTemp_C_5, col="medium sea green", type='l', lwd=1.5)

plot(t2021$DateTime,t2021$RDOTemp_C_9, main="RDO 9m vs. Temp String", xlab="Time", ylab="degrees C", type='l', col="black", lwd=1.5, ylim=c(0,10))
points(t2021$DateTime, t2021$ThermistorTemp_C_9, col="blue4", type='l', lwd=1.5)






###########################################################################################################################################################################


#Creating a new flag "6" started in 2019: "very questionable value due to potential fouling. Values adjusted using a linear or square root function to match high-resolution CTD profiles are given in RDO_mgL_5 and RDO_sat_percent_5"
#Creating a new flag "5" which means that the values are very questionable due to fouling but not adjusted (starting in 2019)

#check the DO with graphs
  
#convert DateTime to make the graphs nice
  #Convert the time and put it in current time zone so the times line up when changing NAs.
  #+5 is during EDT and +4 is during EST(make sure to check this in December)
  #Have to do strpttime or you get some NAs in the DateTime column
  #Do this so the graphing is nice
  catdata_flag$DateTime<-as.POSIXct(strptime(catdata_flag$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5")

#Check 1.6m EXO DO data
#Plot 2018-current
EXODO <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXODO_mgL_1)) +
  geom_point()
EXODO
# 
# # Plotly so can pick out questionable values
# ggplotly(EXODO)
# 
# #Plot Just the current year
 EXODO_21=catdata_flag%>%
   filter(DateTime>start_time & DateTime<end_time)%>%
   ggplot(.,aes(x = DateTime, y = EXODO_mgL_1)) +
   geom_point()
 EXODO_21
# # Plotly so can pick out questionable values
# ggplotly(EXODO_21)
# 
# # Check the 5m RDO
# # plot of 5m DO from 2018 to present
RDO5 <- ggplot(data = catdata_flag, aes(x = DateTime, y = RDO_mgL_5)) +
  geom_point()
RDO5
# # Plotly so can pick out questionable values
# ggplotly(RDO5)
# 
# # Plot the current year
RDO5_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = RDO_mgL_5)) +
  geom_point()
RDO5_21
# # Plotly so can pick out questionable values
# ggplotly(RDO5_21)
# 
# 
# # Plot the 9m Do
# # From 2018-current 
RDO9=ggplot(catdata_flag, aes(x = DateTime, y = RDO_mgL_9)) +
  geom_point()
RDO9
# # Plotly so can pick out questionable values
# ggplotly(RDO9)
# 
# #Just the current year
 RDO9_21=catdata_flag%>%
   #filter(DateTime>start_time & DateTime<end_time)%>%
   filter(DateTime>"2021-05-31 00:00" & DateTime<"2021-06-09 00:00")
   ggplot(.,aes(x = DateTime, y = RDO_mgL_9)) +
   geom_point()
 RDO9_21
# # Plotly so can pick out questionable values
 ggplotly(RDO9_21)
# 


###########################################################################################################################################################################
# chl and phyco visual qaqc-plot to see if everything looks right

# Chla
# Plot for 2018-current
chl_ugl <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOChla_ugL_1)) +
  geom_point()
chl_ugl
# # Plotly so can pick out questionable values
# ggplotly(chl_ugl)
# 
# # Plot just the current year
chl_ugl_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOChla_ugL_1)) +
  geom_point()
chl_ugl_21
# 
# 
# # plot the daily mean
# # calculate the daily mean
chl_mean <- catdata_flag %>%
 select(DateTime, EXOChla_ugL_1) %>%
 mutate(day = date(DateTime)) %>%
 group_by(day) %>%
 mutate(daily_mean = mean(EXOChla_ugL_1, na.rm = TRUE)) %>%
 distinct(day, .keep_all = TRUE)
# 
# # plot the daily mean
chl_mean_plot <- ggplot(data = chl_mean, aes(x = day, y = daily_mean)) +
 geom_point()
chl_mean_plot
# # Plotly so can pick out questionable values
# ggplotly(chl_mean)
# 
# # Plot the chla and the daily mean on the same graph
plot(catdata_flag$DateTime, catdata_flag$EXOChla_ugL_1)
points(chl_mean$DateTime, chl_mean$daily_mean, type="l", col="green")
# 
# 
# # Chla-RFU
# # Plot 2018-current
chl_rfu <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOChla_RFU_1)) +
  geom_point()
chl_rfu
# # Plotly so can pick out questionable values
# ggplotly(chl_rfu)
# 
# 
# # Just the current year
chl_rfu_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOChla_RFU_1)) +
  geom_point()
chl_rfu_21
# 
# # Phyco-RFU
# # Plot 2018-current
phyco_rfu <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOBGAPC_RFU_1)) +
  geom_point()
phyco_rfu
# # Plotly so can pick out questionable values
# ggplotly(phyco_rfu)
# 
# # Just the current year
phyco_rfu_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOBGAPC_RFU_1)) +
  geom_point()
# # Plotly so can pick out questionable values
# ggplotly(phyco_rfu_21)
# 
# #################################################################################################
# 
# # fDOM-RFU
# # Plot 2018-current
fDOM_rfu <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOfDOM_RFU_1)) +
  geom_point()
fDOM_rfu
# # Plotly so can pick out questionable values
# ggplotly(fDOM_rfu)
# 
# # Just the current year
fDOM_rfu_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOfDOM_RFU_1)) +
  geom_point()
 fDOM_rfu_21
# # Plotly so can pick out questionable values
# ggplotly(fDOM_rfu_21)
# 
# # fDOM-QSU
# # Plot 2018-current
fDOM_qsu <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOfDOM_QSU_1)) +
  geom_point()
fDOM_qsu
# # Plotly so can pick out questionable values
# ggplotly(fDOM_qsu)
# 
# # Just the current year
fDOM_qsu_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOfDOM_QSU_1)) +
  geom_point()
 fDOM_qsu_21
# # Plotly so can pick out questionable values
# ggplotly(fDOM_qsu_21)
# 


###########################################################################################################################################################################
# conductivity visual QAQC

#convert DateTime to make the graphs nice
#Convert the time and put it in current time zone so the times line up when changing NAs.
#+5 is during EDT and +4 is during EST(make sure to check this in December)
#Have to do strpttime or you get some NAs in the DateTime column
#Do this so the graphing is nice
 #catdata_flag$DateTime<-as.POSIXct(strptime(catdata_flag$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5")
# 
# #Plot from 2018-current
Cond <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOCond_uScm_1)) +
  geom_point()
Cond
# # Plotly so can pick out questionable values
# ggplotly(Cond)
 
 #Look at the Flag 5 values
 co=catdata_flag%>%filter(Flag_Cond<5)
 
 plot(catdata_flag$DateTime,catdata_flag$EXOCond_uScm_1, col="red")
 points(co$DateTime,co$EXOCond_uScm_1, col="black")
# 
# 
Con=catdata_flag%>%
  select(c(DateTime,EXOTemp_C_1,EXOCond_uScm_1,EXOSpCond_uScm_1,EXO_wiper,EXO_pressure))
# #Just the current year
# 
Cond_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOCond_uScm_1)) +
  geom_point()
Cond_21
# # Plotly so can pick out questionable values
# ggplotly(Cond_21)
# 
# #Specific Conductivity
# 
# #Plot from 2018-current
SpCond <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOSpCond_uScm_1)) +
  geom_point()
SpCond
# # Plotly so can pick out questionable values
# ggplotly(SpCond)
# 
# #Just the current year
SpCond_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOSpCond_uScm_1)) +
  geom_point()
SpCond_21
# # Plotly so can pick out questionable values
# ggplotly(SpCond_21)
# 
# #Total Dissolved Solids
TDS <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXOTDS_mgL_1)) +
  geom_point()
TDS
# # Plotly so can pick out questionable values
# ggplotly(TDS)
# 
# #Just the current year
TDS_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXOTDS_mgL_1)) +
  geom_point()
TDS_21
# # Plotly so can pick out questionable values
# ggplotly(TDS_21)
# 
# #Depth
Depth <- ggplot(data = catdata_flag, aes(x = DateTime, y = EXO_depth)) +
  geom_point()
Depth
# # Plotly so can pick out questionable values
# ggplotly(Depth)
# 
# #Just the current year
Depth_21=catdata_flag%>%
  filter(DateTime>start_time & DateTime<end_time)%>%
  ggplot(.,aes(x = DateTime, y = EXO_depth)) +
  geom_point()
Depth_21
# # Plotly so can pick out questionable values
# ggplotly(Depth_21)

###########################################################################################################################################################################
# write final csv ----
# some final checking of flags
# for (i in 1:nrow(catdata_all)) {
#   if(is.na(catdata_all$Flag_fDOM[i])){
#     catdata_all$Flag_fDOM[i] <- 0
#   }
# }

#check the flag column
Flags=catdata_flag%>%
  select(starts_with("Flag"))

for(b in 1:nrow(Flags)){
  print(colnames(Flags[b]))
  print(table(Flags[,b], useNA = "always"))
}

#Order by date and time

catdata_flag <- catdata_flag[order(catdata_flag$DateTime),]


# for(b in 1:nrow())
# str(catdata_flag)

#rearrange the cols
catdata_flag <- catdata_flag %>%
  select(Reservoir, Site, DateTime, ThermistorTemp_C_surface:ThermistorTemp_C_9,
         RDO_mgL_5, RDOsat_percent_5, RDO_mgL_5_adjusted, RDOsat_percent_5_adjusted,
         RDOTemp_C_5, RDO_mgL_9, RDOsat_percent_9, RDO_mgL_9_adjusted, RDOsat_percent_9_adjusted, RDOTemp_C_9,
         EXOTemp_C_1, EXOCond_uScm_1, EXOSpCond_uScm_1, EXOTDS_mgL_1, EXODOsat_percent_1,
         EXODO_mgL_1, EXOChla_RFU_1, EXOChla_ugL_1, EXOBGAPC_RFU_1, EXOBGAPC_ugL_1,
         EXOfDOM_RFU_1, EXOfDOM_QSU_1, EXO_pressure_psi, EXO_depth_m, EXO_battery_V, EXO_cablepower_V,
         EXO_wiper_V, Lvl_psi_9, LvlTemp_C_9, RECORD, CR6_Batt_V, CR6Panel_Temp_C,
         Flag_All, Flag_Temp_Surf:Flag_Temp_9,Flag_DO_5_obs, Flag_DO_5_sat, Flag_DO_5_temp,
         Flag_DO_9_obs, Flag_DO_9_sat, Flag_DO_9_temp,Flag_EXOTemp, Flag_Cond, Flag_SpCond,Flag_TDS,
         Flag_DO_1_sat, Flag_DO_1_obs, Flag_Chla_RFU,Flag_Chla_ugL, Flag_Phyco_RFU,Flag_Phyco_ugL,
         Flag_fDOM_RFU,Flag_fDOM_QSU, Flag_LvlPres)

# convert datetimes to characters so that they are properly formatted in the output file
catdata_flag$DateTime <- as.character(catdata_flag$DateTime)
  
  
write.csv(catdata_flag, paste0(folder, '/FCR_Catwalk_2018_2021.csv'), row.names = FALSE)


