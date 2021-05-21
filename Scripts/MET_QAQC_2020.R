###QA/QC script for FCR meteorological data Environmental Data Initiative publishing
###Script written by Bethany Bookout & Cayelan Carey & Adrienne Breef-Pilz
###Last modified 20 Jan 2021 to publish 2015-2020 met data to EDI
###Contact info: Cayelan Carey, Virginia Tech, cayelan@vt.edu

# Modified by A Hounshell on 21 May 2021 to download/clean 2021 data

###1) Install and load packages needed
pacman::p_load(lubridate,tidyverse)

#install.packages("lubridate")
#install.packages("tidyverse")
#library(lubridate)
#library(tidyverse)
#rm(list=ls()) #let's start with a blank slate

###2) Download the "raw" meteorological FCR datasets from GitHub and aggregate into 1 file: 
#a. Past Met data, manual downloads
#download current met data from GitHub
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/FCRmet.csv", "FCRmet.csv")


#original raw files from 2015-2019
#RawMet_1516=read.csv('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2015_2016.csv',header = T) #2015-2016 data
#RawMet_17=read.csv('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2017.csv',header = T) #2017 data
#RawMet_18=read.csv('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2018.csv',header = T) #2018 data
#RawMet_19=read.csv('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2019.csv',header = T) #2019 data
#download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/FCRmet_legacy_2020.csv", "FCRmet_20.csv")
#RawMet_20=read.csv("FCRmet_20.csv", skip = 4, header = F) 

#if(length(names(RawMet_20))>17){ #removes NR01TK_Avg column, which was downloaded on some but not all days
#  RawMet_20$V17<-NULL #remove extra column
#}

#names(RawMet_20) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
#                   "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
#                   "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
#                   "ShortwaveRadiationDown_Average_W_m2", "InfaredRadiationUp_Average_W_m2",
#                   "InfaredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")

#RawMet_20$DateTime<-as.POSIXct(strptime(RawMet_20$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5") #get dates aligned

# Try downloading 2020 data from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2020.csv")
#download.file(inUrl1,infile1,method="curl")

#met_edi <- read.csv("./Data/Met_final_2015_2020.csv", header=T) %>%
#  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

#combine all the years to create a past data frame
#Met_past = rbind(RawMet_1516,RawMet_17,RawMet_18,RawMet_19) #merges 2018 with data

#pesky X column
#Met_past$X<-NULL

#names(Met_past) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
#               "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
#               "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
#               "ShortwaveRadiationDown_Average_W_m2", "InfaredRadiationUp_Average_W_m2",
#               "InfaredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")
#put datetime in a useable form
#Met_past$DateTime<-as.POSIXct(strptime(Met_past$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5")


#Note for publishing the 2015-2018 dataset, I skipped this step, but for future archiving,
# you would need to add in the most recent data that are archived. To do so:

####Aggregate legacy file of 2020 data from Github ####
Met_now=read.csv("FCRmet.csv", skip = 4, header = F) 
#loads in data from SCC_data repository for latest push
if(length(names(Met_now))>17){ #removes NR01TK_Avg column, which was downloaded on some but not all days
  Met_now$V17<-NULL #remove extra column
}

#Because using the finalized version of past data must QAQC the current the metdata to add flags and then add the past file
#renames and reformats dataset for easy bind
names(Met_now) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
                   "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
                   "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
                   "ShortwaveRadiationDown_Average_W_m2", "InfaredRadiationUp_Average_W_m2",
                   "InfaredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")
#met_timechange=max(which(Met_now$TIMESTAMP=="2019-04-15 10:19:00")) #shows time point when met station was switched from GMT -4 to GMT -5
Met_now$DateTime<-as.POSIXct(strptime(Met_now$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5") #get dates aligned
#Met_now$TIMESTAMP[c(1:met_timechange-1)]<-with_tz(force_tz(Met_now$TIMESTAMP[c(1:met_timechange-1)],"Etc/GMT+4"), "Etc/GMT+5") #pre time change data gets assigned proper timezone then corrected to GMT -5 to match the rest of the data set

#Met_now=Met_now[-c(which(is.na(Met_now$TIMESTAMP))),] #checks for NA from failed parse, and deletes from dataset; no NAs present
#Met_now$PAR_Tot_Tot=as.numeric(Met_now$PAR_Tot_Tot) #matching str to past data
#str(Met_now); str(Met_past) #checks structure for matching of historical dataset with
# current data


#Set end date
#Met_now=Met_now[year(Met_now$DateTime)<2021,]

####3) Aggregate data set for QA/QC ####
#add the past and now together
Met=Met_now

 #reset data so you don't have to load from scratch 
#Met = Met_past #if you are *only* archiving past data - note! you won't do this if you are adding in new data after 2018
#Met$TIMESTAMP=ymd_hms(Met$TIMESTAMP, tz="Etc/GMT+5") #formats timestamp as double check

#order data by timestamp
Met=Met[order(Met$DateTime),]
Met$DOY=yday(Met$DateTime)


#check record for gaps
#daily record gaps by day of year
for(i in 2:nrow(Met)){ #this identifies if there are any data gaps in the long-term record, and where they are by record number
  if(Met$DOY[i]-Met$DOY[i-1]>1){
    print(c(Met$DateTime[i-1],Met$DateTime[i]))
  }
}
#sub-daily record gaps by record number
for(i in 2:length(Met$Record)){ #this identifies if there are any data gaps in the long-term record, and where they are by record number
  if(abs(Met$Record[i]-Met$Record[i-1])>1){
    print(c(Met$DateTime[i-1],Met$DateTime[i]))
  }
}


#EDI Column names
names(Met) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
               "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
               "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
               "ShortwaveRadiationDown_Average_W_m2", "InfaredRadiationUp_Average_W_m2",
               "InfaredRadiationDown_Average_W_m2", "Albedo_Average_W_m2", "DOY") #finalized column names
Met$Reservoir="FCR" #add reservoir name for EDI archiving
Met$Site=50 #add site column for EDI archiving
Met_raw=Met #Met=Met_raw; reset your data, compare QAQC

#plot(Met$DateTime, Met$BP_Average_kPa, type= "l")

####4) Load in maintenance txt file #### 
# the maintenance file tracks when sensors were repaired or offline due to maintenance
RemoveMet=read.table("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/MET_MaintenanceLog.txt", sep = ",", header = T)
#str(RemoveMet)
RemoveMet$TIMESTAMP_start=ymd_hms(RemoveMet$TIMESTAMP_start, tz="Etc/GMT+5")#setting time zone
RemoveMet$TIMESTAMP_end=ymd_hms(RemoveMet$TIMESTAMP_end, tz="Etc/GMT+5") #setting time zone
RemoveMet$notes=as.character(RemoveMet$notes)

####5) Create data flags for publishing ####
#get rid of NaNs
#create flag + notes columns for data columns c(5:17)
#set flag 2
for(i in 5:17) { #for loop to create new columns in data frame
  Met[,paste0("Flag_",colnames(Met[i]))] <- 0 #creates flag column + name of variable
  Met[,paste0("Note_",colnames(Met[i]))] <- NA #creates note column + names of variable
  Met[which(is.nan(Met[,i])),i] <- NA
  Met[c(which(is.na(Met[,i]))),paste0("Flag_",colnames(Met[i]))] <-2 #puts in flag 2
  Met[c(which(is.na(Met[,i]))),paste0("Note_",colnames(Met[i]))] <- "Sample not collected" #note for flag 2
}

#change the rain totals to per a minute for the time that was recorded in 5 minute intervals
#Met$Flag_Rain_Total_mm=ifelse((Met$DateTime>"2015-07-01 00:00:00"&Met$DateTime<"2015-07-13 12:20:00"&Met$Rain_Total_mm>0),4,Met$Flag_Rain_Total_mm)
#Met$Note_Rain_Total_mm=ifelse((Met$DateTime>"2015-07-01 00:00:00"&Met$DateTime<"2015-07-13 12:20:00"&Met$Rain_Total_mm>0),"Change_to_mm_per_min",Met$Note_Rain_Total_mm)
#Met$Rain_Total_mm=ifelse((Met$DateTime>"2015-07-01 00:00:00"&Met$DateTime<"2015-07-13 12:20:00"&Met$Rain_Total_mm>0), (Met$Rain_Total_mm/5) , Met$Rain_Total_mm)
#Air temperature data cleaning
#separate data by date; before and after air temp filter installed 2019-02-21 11:50:00
#Met_prefilter=Met[Met$DateTime<"2019-02-21 11:50:00",]
#Met_postfilter=Met[Met$DateTime>="2019-02-21 11:50:00",]
#create linear model between the panel temp and the air temperature sensor for 2015,
# (lm_Panel2015) and then apply correction to air temperature dataset

#Check the air filter to the panel figure to see what it looks like
#plot(Met$DateTime, Met$AirTemp_Average_C, type='l')
#points(Met$DateTime, Met$CR3000Panel_temp_C, col="red", type='l')
#legend("topright", c("Air Temp", "Panel Temp"), text.col=c("black", "red"), cex=0.75)

#No pre and post filter but run all through the QAQC 
#MetAir_2015=Met[Met$DateTime<"2016-01-01 00:00:00",c(1,4,8)]
#lm_Panel2015=lm(MetAir_2015$AirTemp_Average_C ~ MetAir_2015$CR3000Panel_temp_C)
#summary(lm_Panel2015)#gives data on linear model parameters

#make the panelTemp from character into numeric 
#Met$CR3000Panel_temp_C=as.numeric(Met$CR3000Panel_temp_C)

#if Air - Panel > 3 sd(lm_Panel2015) then replace with PanelTemp predicted by lm equation rather than raw value
#Met$Flag_AirTemp_Average_C=ifelse((Met$AirTemp_Average_C - (1.6278+(0.9008*Met$CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)), 4, Met$Flag_AirTemp_Average_C)
#Met$Note_AirTemp_Average_C=ifelse((Met$AirTemp_Average_C - (1.6278+(0.9008*Met$CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)),"Substituted_value_calculated_from_Panel_Temp_and_linear_model", Met$Note_AirTemp_Average_C)
#Met$AirTemp_Average_C=ifelse((Met$AirTemp_Average_C - (1.6278+(0.9008*Met$CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)),(1.6278+(0.9008*Met$CR3000Panel_temp_C)), Met$AirTemp_Average_C)

#merge back air temp correction data
#Met<-rbind(Met_prefilter,Met_postfilter)

#put all airtemp through the panel QAQC check
#plot(Met_raw$DateTime, Met_raw$AirTemp_Average_C, type= 'l')
#points(Met$DateTime, Met$AirTemp_Average_C, type = "l", col="red")

#Air temp maximum set
Met$Flag_AirTemp_Average_C=ifelse(Met$AirTemp_Average_C>40.6,4,Met$Flag_AirTemp_Average_C)
Met$Note_AirTemp_Average_C=ifelse(Met$AirTemp_Average_C>40.6,"Outlier_set_to_NA",Met$Note_AirTemp_Average_C)
Met$AirTemp_Average_C=ifelse(Met$AirTemp_Average_C>40.6,NA,Met$AirTemp_Average_C)

#check the infared radiation-if there are any low values
#plot(Met_raw$DateTime, Met_raw$InfaredRadiationDown_Average_W_m2,main = "Raw Infared Radiation 2020", ylab = "Average_W_m2", type = "l")

#Infared radiation cleaning
#fix infrared radiation voltage reading after airtemp correction
#only need to do this for data from 2015 to  2016-07-25 10:12:00

#Met$Flag_InfaredRadiationUp_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationUp_Average_W_m2<100,4,Met$Flag_InfaredRadiationUp_Average_W_m2)
#Met$Note_InfaredRadiationUp_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationUp_Average_W_m2<100,"Value_corrected_from_Voltage_with_InfRadUp_equation_as_described_in_metadata",Met$Note_InfaredRadiationUp_Average_W_m2)
#Met$InfaredRadiationUp_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationUp_Average_W_m2<100,Met$InfaredRadiationUp_Average_W_m2+5.67*10^-8*(Met$AirTemp_Average_C+273.15)^4,Met$InfaredRadiationUp_Average_W_m2)

#Met$Flag_InfaredRadiationDown_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationDown_Average_W_m2<250,4,Met$Flag_InfaredRadiationDown_Average_W_m2)
#Met$Note_InfaredRadiationDown_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationDown_Average_W_m2<250,"Value_corrected_from_Voltage_with_InfRadDn_equation_as_described_in_metadata",Met$Note_InfaredRadiationDown_Average_W_m2)
#Met$InfaredRadiationDown_Average_W_m2=ifelse(Met$DateTime<"2016-07-25 10:12:00" &Met$InfaredRadiationDown_Average_W_m2<250,5.67*10^-8*(Met$AirTemp_Average_C+273.15)^4,Met$InfaredRadiationDown_Average_W_m2)
  

#Mean correction for InfRadDown (needs to be after voltage correction)
#Using 2018 data, taking the mean and sd of values on DOY to correct to
#Met$DOY=yday(Met$DateTime)
#Met_infrad=Met[year(Met$DateTime)<2018,]
#Met_infrad$infradavg=ave(Met_infrad$InfaredRadiationDown_Average_W_m2, Met_infrad$DOY) #creating column with mean of infraddown by day of year
#Met_infrad$infradsd=ave(Met_infrad$InfaredRadiationDown_Average_W_m2, Met_infrad$DOY, FUN = sd) #creating column with sd of infraddown by day of year
#Met_infrad=unique(Met_infrad[,c(18,47,48)])

#Met=merge(Met, Met_infrad, by = "DOY") #putting in columns for infrared mean and sd by DOY into main data set
#Met=Met[order(Met$DateTime),] #ordering table after merging and removing unnecessary columns

#Met$Flag_InfaredRadiationDown_Average_W_m2=ifelse((Met$InfaredRadiationDown_Average_W_m2-Met$infradavg)<(-3*Met$infradsd),4,Met$Flag_InfaredRadiationDown_Average_W_m2)
#Met$Note_InfaredRadiationDown_Average_W_m2=ifelse((Met$InfaredRadiationDown_Average_W_m2-Met$infradavg)<(-3*Met$infradsd),"Value_corrected_from_mean_InfRadDn_before_fouling_as_described_in_metadata",Met$Note_InfaredRadiationDown_Average_W_m2)
#Met$InfaredRadiationDown_Average_W_m2=ifelse((Met$InfaredRadiationDown_Average_W_m2-Met$infradavg)<(-3*Met$infradsd),Met$infradavg,Met$InfaredRadiationDown_Average_W_m2)

#Met=Met[,-c(1,47,48)]

#take out impossible outliers and Infinite values before the other outliers are removed
#set flag 3 (see metadata: this corrects for impossible outliers)
#for(i in 5:17) { #for loop to create new columns in data frame
#  Met[c(which(is.infinite(Met[,i]))),paste0("Flag_",colnames(Met[i]))] <-3 #puts in flag 3
#  Met[c(which(is.infinite(Met[,i]))),paste0("Note_",colnames(Met[i]))] <- "Infinite_value_set_to_NA" #note for flag 3
#  Met[c(which(is.infinite(Met[,i]))),i] <- NA #set infinite vals to NA
  
#  if(i!=8) { #flag 3 for negative values for everything except air temp
#    Met[c(which((Met[,i]<0))),paste0("Flag_",colnames(Met[i]))] <- 3
#    Met[c(which((Met[,i]<0))),paste0("Note_",colnames(Met[i]))] <- "Negative_value_set_to_0"
#    Met[c(which((Met[,i]<0))),i] <- 0 #replaces value with 0
#  }
#  if(i==9) { #flag for RH over 100
#    Met[c(which((Met[,i]>100))),paste0("Flag_",colnames(Met[i]))] <- 3
#    Met[c(which((Met[,i]>100))),paste0("Note_",colnames(Met[i]))] <- "Value_set_to_100"
#    Met[c(which((Met[,i]>100))),i] <- 100 #replaces value with 100
#  }
#}


#full data set QAQC
#Inf outliers, must go after corrections
Met$Flag_InfaredRadiationUp_Average_W_m2=ifelse(Met$InfaredRadiationUp_Average_W_m2<150,4,Met$Flag_InfaredRadiationUp_Average_W_m2)
Met$Note_InfaredRadiationUp_Average_W_m2=ifelse(Met$InfaredRadiationUp_Average_W_m2<150,"Outlier_set_to_NA",Met$Note_InfaredRadiationUp_Average_W_m2)
Met$InfaredRadiationUp_Average_W_m2=ifelse(Met$InfaredRadiationUp_Average_W_m2<150,NA,Met$InfaredRadiationUp_Average_W_m2)

#Remove barometric pressure outliers
Met$Flag_BP_Average_kPa=ifelse(Met$BP_Average_kPa<98.5, 4,Met$Flag_BP_Average_kPa)
Met$Note_BP_Average_kPa=ifelse(Met$BP_Average_kPa<98.5,"Outlier_set_to_NA",Met$Note_BP_Average_kPa)
Met$BP_Average_kPa=ifelse(Met$BP_Average_kPa<98.5,NA,Met$BP_Average_kPa)

#Remove total PAR (PAR_Tot) outliers
Met$Flag_PAR_Total_mmol_m2=ifelse(Met$PAR_Total_mmol_m2>200, 4, Met$Flag_PAR_Total_mmol_m2)
Met$Note_PAR_Total_mmol_m2=ifelse(Met$PAR_Total_mmol_m2>200, "Outlier_set_to_NA", Met$Note_PAR_Total_mmol_m2)
Met$PAR_Total_mmol_m2=ifelse(Met$PAR_Total_mmol_m2>200, NA, Met$PAR_Total_mmol_m2)

#flagging IR down values in the summer of 2020 that are lower than 400
#Met$Flag_InfaredRadiationDown_Average_W_m2=ifelse((Met$DateTime>"2020-06-10 00:00:00"&Met$DateTime<"2020-07-15 23:59:00"&Met$InfaredRadiationDown_Average_W_m2<400),5,Met$Flag_InfaredRadiationDown_Average_W_m2)
#Met$Note_InfaredRadiationDown_Average_W_m2=ifelse((Met$DateTime>"2020-06-10 00:00:00"&Met$DateTime<"2020-07-15 23:59:00"&Met$InfaredRadiationDown_Average_W_m2<400),"Questionable_value",Met$Note_InfaredRadiationDown_Average_W_m2)
#Met$InfaredRadiationDown_Average_W_m2=ifelse((Met$DateTime>"2020-06-10 00:00:00"&Met$DateTime<"2020-07-15 23:59:00"&Met$InfaredRadiationDown_Average_W_m2<400),(Met$InfaredRadiationDown_Average_W_m2),Met$InfaredRadiationDown_Average_W_m2)


#PAR TOT and AVG 2019 QAQC
#Take out bad data due to sensor failure. Aug 4-5; Aug 17 - Dec 13
#PAR_failure=c("2019-08-04":"2019-08-05", "2019-08-15":"2019-12-13")
#Met$Flag_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), 4, Met$Flag_PAR_Total_mmol_m2)
#Met$Note_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Total_mmol_m2)
#Met$PAR_Total_mmol_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), NA, Met$PAR_Total_mmol_m2)

#Met$Flag_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), 4, Met$Flag_PAR_Average_umol_s_m2)
#Met$Note_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Average_umol_s_m2)
#Met$PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2019-08-04 00:00:00"&Met$DateTime<"2019-08-05 00:00:00"), NA, Met$PAR_Average_umol_s_m2)

#Met$Flag_PAR_Total_mmol_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", 4, Met$Flag_PAR_Total_mmol_m2)
#Met$Note_PAR_Total_mmol_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Total_mmol_m2)
#Met$PAR_Total_mmol_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", NA, Met$PAR_Total_mmol_m2)

#Met$Flag_PAR_Average_umol_s_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", 4, Met$Flag_PAR_Average_umol_s_m2)
#Met$Note_PAR_Average_umol_s_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Average_umol_s_m2)
#Met$PAR_Average_umol_s_m2=ifelse(Met$DateTime>"2019-08-15 00:00:00"&Met$DateTime<"2020-01-01 00:00:00", NA, Met$PAR_Average_umol_s_m2)

#PAR TOT and AVG 2020 QAQC
#Take out bad data due to sensor failure. Jul 3-21; Nov 1-24; Dec 1-31
#PAR_failure=c("2020-07-03 11:56":"2020-07-21 13:38", "2020-11-01 12:37":"2020-11-24 12:32", "2020-12-01 09:54":"2020-12-31 23:59")
#Met$Flag_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), 4, Met$Flag_PAR_Total_mmol_m2)
#Met$Note_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Total_mmol_m2)
#Met$PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), NA, Met$PAR_Total_mmol_m2)

#Met$Flag_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), 4, Met$Flag_PAR_Average_umol_s_m2)
#Met$Note_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Average_umol_s_m2)
#Met$PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-07-03 11:56:00"&Met$DateTime<"2020-07-21 13:38:00"), NA, Met$PAR_Average_umol_s_m2)

#Met$Flag_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), 4, Met$Flag_PAR_Total_mmol_m2)
#Met$Note_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Total_mmol_m2)
#Met$PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), NA, Met$PAR_Total_mmol_m2)

#Met$Flag_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), 4, Met$Flag_PAR_Average_umol_s_m2)
#Met$Note_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Average_umol_s_m2)
#Met$PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-11-01 12:37:00"&Met$DateTime<"2020-11-24 12:32:00"), NA, Met$PAR_Average_umol_s_m2)

Met$Flag_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), 4, Met$Flag_PAR_Total_mmol_m2)
Met$Note_PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Total_mmol_m2)
Met$PAR_Total_mmol_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), NA, Met$PAR_Total_mmol_m2)

Met$Flag_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), 4, Met$Flag_PAR_Average_umol_s_m2)
Met$Note_PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), "Outlier_set_to_NA_due_to_probable_sensor_failure", Met$Note_PAR_Average_umol_s_m2)
Met$PAR_Average_umol_s_m2=ifelse((Met$DateTime>"2020-12-01 09:54:00"), NA, Met$PAR_Average_umol_s_m2)

#Remove shortwave radiation outliers
#first shortwave upwelling
Met$Flag_ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1500, 4, Met$Flag_ShortwaveRadiationUp_Average_W_m2)
Met$Note_ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1500, "Outlier_set_to_NA", Met$Note_ShortwaveRadiationUp_Average_W_m2)
Met$ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1500, NA, Met$ShortwaveRadiationUp_Average_W_m2)


#add a flag for suspect values from 1499 to 1300
#above 1500 already set to NA
Met$Flag_ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1300, 5, Met$Flag_ShortwaveRadiationUp_Average_W_m2)
Met$Note_ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1300, "Questionable_value", Met$Note_ShortwaveRadiationUp_Average_W_m2)
Met$ShortwaveRadiationUp_Average_W_m2=ifelse(Met$ShortwaveRadiationUp_Average_W_m2>1300, (Met$ShortwaveRadiationUp_Average_W_m2), Met$ShortwaveRadiationUp_Average_W_m2)

#Eliminate the High Shortwave Radiation Up in 2018 before cleaning
#Met$Flag_ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), 4, Met$Flag_ShortwaveRadiationUp_Average_W_m2)
#Met$Note_ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), "Outlier_set_to_NA", Met$Note_ShortwaveRadiationUp_Average_W_m2)
#Met$ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), NA, Met$ShortwaveRadiationUp_Average_W_m2)

#and then shortwave downwelling (what goes up must come down)
Met$Flag_ShortwaveRadiationDown_Average_W_m2=ifelse(Met$ShortwaveRadiationDown_Average_W_m2>300, 4, Met$Flag_ShortwaveRadiationDown_Average_W_m2)
Met$Note_ShortwaveRadiationDown_Average_W_m2=ifelse(Met$ShortwaveRadiationDown_Average_W_m2>300, "Outlier_set_to_NA", Met$Note_ShortwaveRadiationDown_Average_W_m2)
Met$ShortwaveRadiationDown_Average_W_m2=ifelse(Met$ShortwaveRadiationDown_Average_W_m2>300, NA, Met$ShortwaveRadiationDown_Average_W_m2)

#Eliminate the High Shortwave Radiation Up in 2018 before cleaning
#Met$Flag_ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), 4, Met$Flag_ShortwaveRadiationUp_Average_W_m2)
#Met$Note_ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), "Outlier_set_to_NA", Met$Note_ShortwaveRadiationUp_Average_W_m2)
#Met$ShortwaveRadiationUp_Average_W_m2=ifelse((Met$DateTime>"2018-09-05 12:00:00"&Met$DateTime<"2018-09-07 20:00:00"), NA, Met$ShortwaveRadiationUp_Average_W_m2)


#shortwave downwelling random point in 2016
#Met$Flag_ShortwaveRadiationDown_Average_W_m2=ifelse((Met$DateTime=="2016-07-18 13:38:00"), 4, Met$Flag_ShortwaveRadiationDown_Average_W_m2)
#Met$Note_ShortwaveRadiationDown_Average_W_m2=ifelse((Met$DateTime=="2016-07-18 13:38:00"), "Outlier_set_to_NA", Met$Note_ShortwaveRadiationDown_Average_W_m2)
#Met$ShortwaveRadiationDown_Average_W_m2=ifelse((Met$DateTime=="2016-07-18 13:38:00"), NA, Met$ShortwaveRadiationDown_Average_W_m2)


#Remove albedo outliers
#over 1000
Met$Flag_Albedo_Average_W_m2=ifelse(Met$Albedo_Average_W_m2>1000, 4, Met$Flag_Albedo_Average_W_m2)
Met$Note_Albedo_Average_W_m2=ifelse(Met$Albedo_Average_W_m2>1000, "Outlier_set_to_NA", Met$Note_Albedo_Average_W_m2)
Met$Albedo_Average_W_m2=ifelse(Met$Albedo_Average_W_m2>1000, NA, Met$Albedo_Average_W_m2)

#set to NA when shortwave radiation up is equal to NA
Met$Flag_Albedo_Average_W_m2=ifelse(is.na(Met$ShortwaveRadiationUp_Average_W_m2)|is.na(Met$ShortwaveRadiationDown_Average_W_m2), 4, Met$Flag_Albedo_Average_W_m2)
Met$Note_Albedo_Average_W_m2=ifelse(is.na(Met$ShortwaveRadiationUp_Average_W_m2)|is.na(Met$ShortwaveRadiationDown_Average_W_m2), "Set_to_NA_because_Shortwave_equals_NA", Met$Note_Albedo_Average_W_m2)
Met$Albedo_Average_W_m2=ifelse(is.na(Met$ShortwaveRadiationUp_Average_W_m2)|is.na(Met$ShortwaveRadiationDown_Average_W_m2), NA, Met$Albedo_Average_W_m2)

#create loop putting in maintenance flags 1 + 4 (these are flags for values removed due
# to maintenance and also flags potentially questionable values)
for(j in 1:nrow(RemoveMet)){
  print(j) # #if statement to only write in flag 4 if there are no other flags
  if(RemoveMet$flag[j]==4){
    
    Met[c(which(Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3] & (Met[,paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]==0))), paste0("Note_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$notes[j]#same as above, but for notes
    Met[c(which(Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3] & (Met[,paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]==0))), paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$flag[j]#when met timestamp is between remove timestamp
    #print(j)#and met column derived from remove column
    #matching time frame, inserting flag
  }
  #if flag == 1, set parameter to NA, overwrites any other flag
  
  if(RemoveMet$flag[j]==1){
    #print(j)
    Met[c(which((Met[,1]>=RemoveMet[j,2]) & (Met[,1]<=RemoveMet[j,3]))),paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))] = RemoveMet$flag[j] #when met timestamp is between remove timestamp
    #and met column derived from remove column
    #matching time frame, inserting flag
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], paste0("Note_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$notes[j]#same as above, but for notes
    
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], RemoveMet$colnumber[j]] = NA
  } #replaces value of var with NA
}


#Reorder so the flags are all next to each other
Met=Met%>%
  select(c("Reservoir","Site","DateTime","Record","CR3000_Batt_V","CR3000Panel_temp_C","PAR_Average_umol_s_m2","PAR_Total_mmol_m2","BP_Average_kPa",                          
           "AirTemp_Average_C","RH_percent","Rain_Total_mm","WindSpeed_Average_m_s","WindDir_degrees","ShortwaveRadiationUp_Average_W_m2",       
           "ShortwaveRadiationDown_Average_W_m2","InfaredRadiationUp_Average_W_m2","InfaredRadiationDown_Average_W_m2","Albedo_Average_W_m2",
           "Flag_PAR_Average_umol_s_m2","Flag_PAR_Total_mmol_m2","Flag_BP_Average_kPa","Flag_AirTemp_Average_C","Flag_Rain_Total_mm","Flag_RH_percent",
           "Flag_WindSpeed_Average_m_s","Flag_WindDir_degrees","Flag_ShortwaveRadiationUp_Average_W_m2","Flag_ShortwaveRadiationDown_Average_W_m2",
           "Flag_InfaredRadiationUp_Average_W_m2","Flag_InfaredRadiationDown_Average_W_m2","Flag_Albedo_Average_W_m2","Note_PAR_Average_umol_s_m2",              
           "Note_PAR_Total_mmol_m2","Note_BP_Average_kPa","Note_AirTemp_Average_C","Note_RH_percent","Note_Rain_Total_mm","Note_WindSpeed_Average_m_s",              
           "Note_WindDir_degrees","Note_ShortwaveRadiationUp_Average_W_m2","Note_ShortwaveRadiationDown_Average_W_m2","Note_InfaredRadiationUp_Average_W_m2",       
           "Note_InfaredRadiationDown_Average_W_m2","Note_Albedo_Average_W_m2"))


#prints table of flag frequency
#can't get it to work

#for(b in 20:32) {
  #print(colnames(Met[b]))
  #print(table(Met[,paste0("Flag_*",colnames(Met[b]))]))}



####6) Make plots to view data #####
#plots to check for any wonkiness
x11(); par(mfrow=c(2,2))
plot(Met$DateTime, Met$CR3000_Batt_V, type = 'l')
plot(Met$DateTime, Met$CR3000Panel_temp_C, type = 'l')
#PAR
plot(Met_raw$DateTime, Met_raw$PAR_Average_umol_s_m2, col="red", type='l')
plot(Met$DateTime, Met$PAR_Average_umol_s_m2, type = 'l')
plot(Met_raw$DateTime, Met_raw$PAR_Total_mmol_m2, col="red", type='l')
plot(Met$DateTime, Met$PAR_Total_mmol_m2, type = 'l')

#BP
plot(Met_raw$DateTime, Met_raw$BP_Average_kPa, col="red", type='l')
points(Met$DateTime, Met$BP_Average_kPa, type = 'l')
#Air Temp
plot(Met_raw$DateTime, Met_raw$AirTemp_Average_C, col="red", type='l')
points(Met$DateTime, Met$AirTemp_Average_C, type = 'l')
#RH
plot(Met_raw$DateTime, Met_raw$RH_percent, col="red", type='l')
points(Met$DateTime, Met$RH_percent, type = 'l')
#Rain
plot(Met_raw$DateTime, Met_raw$Rain_Total_mm, col="red", type='h')
points(Met$DateTime, Met$Rain_Total_mm, type = 'h')
#Wind
plot(Met$DateTime, Met$WindSpeed_Average_m_s, type = 'l')
hist(Met$WindDir_degrees)
#SW Radiation
plot(Met_raw$DateTime, Met_raw$ShortwaveRadiationUp_Average_W_m2, col="red", type='l')
points(Met$DateTime, Met$ShortwaveRadiationUp_Average_W_m2, type = 'l')
plot(Met_raw$DateTime, Met_raw$ShortwaveRadiationDown_Average_W_m2, col="red", type='l')
points(Met$DateTime, Met$ShortwaveRadiationDown_Average_W_m2, type = 'l')
#Albedo
plot(Met_raw$DateTime, Met_raw$Albedo_Average_W_m2, col="red", type='l')
plot(Met$DateTime, Met$Albedo_Average_W_m2, type = 'l')
#InfRad
plot(Met_raw$DateTime, Met_raw$InfaredRadiationUp_Average_W_m2, col="red", type='l')
points(Met$DateTime, Met$InfaredRadiationUp_Average_W_m2, type = 'l')
plot(Met_raw$DateTime, Met_raw$InfaredRadiationDown_Average_W_m2, col="red", type='l')
points(Met$DateTime, Met$InfaredRadiationDown_Average_W_m2, type = 'l')


#Met unique values for notes
#for (u in 33:45) {
#  print(colnames(Met[u]))
#  print(unique(Met[,u]))
#}

###Prep RemoveMet for final file version
#names(RemoveMet) = c("Station", "DateTime_start","DateTime_end", "Parameter", "ColumnNumber", "Flag", "Notes") #finalized column names
#RemoveMet$Reservoir= "FCR"#add reservoir name for EDI archiving
#RemoveMet$Site=50 #add site column for EDI archiving

#RemoveMet=RemoveMet[,c(8:9,1:7)]

###7) Write file with final cleaned dataset! ###
Met_final=Met%>%
  select(c("Site","Reservoir","DateTime","Record","CR3000_Batt_V","CR3000Panel_temp_C","PAR_Average_umol_s_m2","PAR_Total_mmol_m2","BP_Average_kPa",                          
           "AirTemp_Average_C","RH_percent","Rain_Total_mm","WindSpeed_Average_m_s","WindDir_degrees","ShortwaveRadiationUp_Average_W_m2",       
           "ShortwaveRadiationDown_Average_W_m2","InfaredRadiationUp_Average_W_m2","InfaredRadiationDown_Average_W_m2","Albedo_Average_W_m2",
           "Flag_PAR_Average_umol_s_m2","Note_PAR_Average_umol_s_m2","Flag_PAR_Total_mmol_m2","Note_PAR_Total_mmol_m2","Flag_BP_Average_kPa",
           "Note_BP_Average_kPa","Flag_AirTemp_Average_C","Note_AirTemp_Average_C","Flag_RH_percent","Note_RH_percent","Flag_Rain_Total_mm",
           "Note_Rain_Total_mm","Flag_WindSpeed_Average_m_s","Note_WindSpeed_Average_m_s",
           "Flag_WindDir_degrees","Note_WindDir_degrees","Flag_ShortwaveRadiationUp_Average_W_m2","Note_ShortwaveRadiationUp_Average_W_m2",
           "Flag_ShortwaveRadiationDown_Average_W_m2","Note_ShortwaveRadiationDown_Average_W_m2",
           "Flag_InfaredRadiationUp_Average_W_m2","Note_InfaredRadiationUp_Average_W_m2",
           "Flag_InfaredRadiationDown_Average_W_m2","Note_InfaredRadiationDown_Average_W_m2","Flag_Albedo_Average_W_m2","Note_Albedo_Average_W_m2"))
#setwd('./Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEML_MetData')
write.csv(Met_final, "./Data/Met_GitHub_2021.csv", row.names=F, quote=F)
#write.csv(RemoveMet, "Met_Maintenance_2015_2020.csv", row.names=F, quote = F)

