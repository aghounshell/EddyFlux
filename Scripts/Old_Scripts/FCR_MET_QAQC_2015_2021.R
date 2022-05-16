###QA/QC script for FCR meteorological data Environmental Data Initiative publishing
###Script written by Bethany Bookout & Cayelan Carey & Adrienne Breef-Pilz
###Last modified 20 Jan 2021 to publish 2015-2020 met data to EDI
###Contact info: Cayelan Carey, Virginia Tech, cayelan@vt.edu

### Modified by A. Hounshell, 06 May 2022, to download and correct 2022 data

###############################################################################

###1) Install and load packages needed

pacman::p_load("RCurl","tidyverse","lubridate", "plotly", "magrittr", "suncalc")

rm(list=ls()) #let's start with a blank slate

# Set working directory
wd <- getwd()
setwd(wd)

#folder <- "./Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEML_FCRMetData/2021/"

###2) Download the "raw" meteorological FCR datasets from GitHub and aggregate into 1 file: 
#a. Past Met data, manual downloads
#download current met data from GitHub
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/FCRmet.csv", paste0(wd, "/Data/FCRmet_2022.csv"))

#download maintenance file
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/MET_MaintenanceLog.txt", paste0(wd, "/Data/FCR_Met_Maintenance_2015_2021.txt"))

#original raw files from 2015-2020
#download.file('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2015_2016.csv',paste0(folder, "misc_data_files/RawMetData_2015_2016.csv")) #2015-2016 data
#download.file('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2017.csv',paste0(folder, "misc_data_files/RawMetData_2017.csv")) #2017 data
#download.file('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2018.csv',paste0(folder, "misc_data_files/RawMetData_2018.csv")) #2018 data
#download.file('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2019.csv',paste0(folder, "misc_data_files/RawMetData_2019.csv")) #2019 data
#download.file('https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataAlreadyUploadedToEDI/CollatedDataForEDI/MetData/RawMetData_2020.csv',paste0(folder, "misc_data_files/RawMetData_2020.csv")) #2019 data


#merge previous metdata

#mydir = "Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEML_FCRMetData/2021/misc_data_files"
#myfiles = list.files(path=mydir, pattern="Raw*", full.names=TRUE)#list the files from BVR platform


#create dataframe for the for loop
#out.file<-""

#combine all of the files into one data sheet, have to come back and fix this loop
#for(k in 1:length(myfiles)){
#  files<-read.csv(myfiles[k],skip= 0, header=T) #get header minus wonky Campbell rows
#  if(length(names(files))==17){
#  names(files) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
#                   "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
#                  "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
#                  "ShortwaveRadiationDown_Average_W_m2", "InfraredRadiationUp_Average_W_m2",
#                   "InfraredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")#rename headers
#  }
#  if(length(names(files))>17){ #removes NR01TK_Avg column, which was downloaded on some but not all days
#    files$NR01TK_Avg<-NULL #remove column
#    names(files) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
#                        "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
#                        "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
#                        "ShortwaveRadiationDown_Average_W_m2", "InfraredRadiationUp_Average_W_m2",
#                        "InfraredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")#rename headers
#  }
#  out.file=rbind(out.file, files)
#}

#out.file=out.file%>%filter(Record!="")

#put datetime in a useable form
#out.file$DateTime<-as.POSIXct(strptime(out.file$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5")


####Add current file from Github ####

Met_now=read.csv("./Data/FCRmet_2022.csv", skip = 4, header = F) 

#loads in data from SCC_data repository for latest push
if(length(names(Met_now))>17){ #removes NR01TK_Avg column, which was downloaded on some but not all days
  Met_now$V17<-NULL #remove extra column
}

#Because using the finalized version of past data must QAQC the current the metdata to add flags and then add the past file
#renames and reformats dataset for easy bind
names(Met_now) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
                   "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
                   "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
                   "ShortwaveRadiationDown_Average_W_m2", "InfraredRadiationUp_Average_W_m2",
                   "InfraredRadiationDown_Average_W_m2", "Albedo_Average_W_m2")

####3) Aggregate data set for QA/QC ####
#add the past and now together
#Met=rbind(out.file,Met_now)

#change the columns from as.character to as.numeric after the merge
Met_now[, c(2:17)] <- sapply(Met_now[, c(2:17)], as.numeric)



#met_na=Met%>%filter(is.na(DateTime))

#Change DateTime when it was changed from EDT to EST
#Met_now$DateTime<-as.POSIXct(strptime(Met_now$DateTime, "%Y-%m-%d %H:%M:%S"), tz = "Etc/GMT+5")
#met_timechange=max(which(Met_now$DateTime=="2019-04-15 10:19:00")) #shows time point when met station was switched from GMT -4(EDT) to GMT -5(EST)
#Met$DateTime<-as.POSIXct(strptime(Met$DateTime, "%Y-%m-%d %H:%M"), tz = "Etc/GMT+5") #get dates aligned
#Met_now$DateTime[c(1:met_timechange-1)]<-with_tz(force_tz(Met_now$DateTime[c(1:met_timechange-1)],"Etc/GMT+4"), "Etc/GMT+5") #pre time change data gets assigned proper timezone then corrected to GMT -5 to match the rest of the data set

#Set end date to years working with 
Met=Met_now%>%filter(DateTime>="2022-01-01 00:00:00")

######################################################################
#Check for record gaps and day gpas

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
# #sub-daily record gaps by record number
 for(i in 2:length(Met$Record)){ #this identifies if there are any data gaps in the long-term record, and where they are by record number
   if(abs(Met$Record[i]-Met$Record[i-1])>1){
     print(c(Met$DateTime[i-1],Met$DateTime[i]))
   }
 }


#EDI Column names
names(Met) = c("DateTime","Record", "CR3000_Batt_V", "CR3000Panel_temp_C", 
               "PAR_Average_umol_s_m2", "PAR_Total_mmol_m2", "BP_Average_kPa", "AirTemp_Average_C", 
               "RH_percent", "Rain_Total_mm", "WindSpeed_Average_m_s", "WindDir_degrees", "ShortwaveRadiationUp_Average_W_m2",
               "ShortwaveRadiationDown_Average_W_m2", "InfraredRadiationUp_Average_W_m2",
               "InfraredRadiationDown_Average_W_m2", "Albedo_Average_W_m2", "DOY") #finalized column names
Met$Reservoir="FCR" #add reservoir name for EDI archiving
Met$Site=50 #add site column for EDI archiving




Met_raw=Met #Met=Met_raw; reset your data, compare QAQC

#Met=Met_raw



####4) Load in maintenance txt file #### 
# the maintenance file tracks when sensors were repaired or offline due to maintenance
RemoveMet=read.csv("./Data/FCR_Met_Maintenance_2015_2021.txt", header = T)
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

#create loop putting in maintenance flags 1 + 4 (these are flags for values removed due
# to maintenance and also flags potentially questionable values)
for(j in 1:nrow(RemoveMet)){
  #print(j) # #if statement to only write in flag 4 if there are no other flags
  if(RemoveMet$flag[j]==4){
    print(j) # #if statement to only write in flag 4 if there are no other flags
    Met[c(which(Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3] & (Met[,paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]==0))), paste0("Note_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$notes[j]#same as above, but for notes
    Met[c(which(Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3] & (Met[,paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]==0))), paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$flag[j]#when met timestamp is between remove timestamp
    #print(j)#and met column derived from remove column
    #matching time frame, inserting flag
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], RemoveMet$colnumber[j]] = NA
  }
  #if flag == 1, set parameter to NA, overwrites any other flag
  
  if(RemoveMet$flag[j]==1){
    print(j)
    Met[c(which((Met[,1]>=RemoveMet[j,2]) & (Met[,1]<=RemoveMet[j,3]))),paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))] = RemoveMet$flag[j] #when met timestamp is between remove timestamp
    #and met column derived from remove column
    #matching time frame, inserting flag
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], paste0("Note_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$notes[j]#same as above, but for notes
    
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], RemoveMet$colnumber[j]] = NA
  } #replaces value of var with NA
  
  if(RemoveMet$flag[j]==5){
    print(j)
    Met[c(which((Met[,1]>=RemoveMet[j,2]) & (Met[,1]<=RemoveMet[j,3]))),paste0("Flag_",colnames(Met[RemoveMet$colnumber[j]]))] = RemoveMet$flag[j] #when met timestamp is between remove timestamp
    #and met column derived from remove column
    #matching time frame, inserting flag
    Met[Met[,1]>=RemoveMet[j,2] & Met[,1]<=RemoveMet[j,3], paste0("Note_",colnames(Met[RemoveMet$colnumber[j]]))]=RemoveMet$notes[j]#same as above, but for notes
  }
}


# # #Filter for just flags
# Flags=Met%>%
#   select(starts_with("Flag"))
# for(b in 1:nrow(Flags)){
#   print(colnames(Flags[b]))
#   print(table(Flags[,b],useNA="always"))
# }


#change the rain totals to per a minute for the time that was recorded in 5 minute intervals
Met=Met%>%
  mutate(
    Rain_Total_mm=as.numeric(Rain_Total_mm))
#  Flag_Rain_Total_mm=ifelse(DateTime>"2015-07-01 00:00:00"& DateTime<"2015-07-13 12:20:00"&Rain_Total_mm>0,4,Flag_Rain_Total_mm),
#  Note_Rain_Total_mm=ifelse(DateTime>"2015-07-01 00:00:00"& DateTime<"2015-07-13 12:20:00"&Rain_Total_mm>0,"Change_to_mm_per_min",Note_Rain_Total_mm),
#  Rain_Total_mm=ifelse(DateTime>"2015-07-01 00:00:00"& DateTime<"2015-07-13 12:20:00"&Rain_Total_mm>0, (Rain_Total_mm/5) ,Rain_Total_mm))

#make columns

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
#Replace missing values with PanelTamp predicted by lm equation


Met=Met%>%
  mutate(
    AirTemp_Average_C=ifelse(Flag_AirTemp_Average_C==2,(1.6278+(0.9008*CR3000Panel_temp_C)), AirTemp_Average_C),
    Note_AirTemp_Average_C=ifelse(Flag_AirTemp_Average_C==2,"Substituted_value_calculated_from_Panel_Temp_and_linear_model", Note_AirTemp_Average_C),
    Flag_AirTemp_Average_C=ifelse(Flag_AirTemp_Average_C==2, 4, Flag_AirTemp_Average_C))


Met=Met%>%
  mutate(
    Note_AirTemp_Average_C=ifelse((AirTemp_Average_C - (1.6278+(0.9008*CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)) & !is.na(AirTemp_Average_C),"Substituted_value_calculated_from_Panel_Temp_and_linear_model", Note_AirTemp_Average_C),
    Flag_AirTemp_Average_C=ifelse((AirTemp_Average_C - (1.6278+(0.9008*CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)) & !is.na(AirTemp_Average_C), 4, Flag_AirTemp_Average_C),
  AirTemp_Average_C=ifelse((AirTemp_Average_C - (1.6278+(0.9008*CR3000Panel_temp_C)))>(3*sd(lm_Panel2015$residuals)) & !is.na(AirTemp_Average_C),(1.6278+(0.9008*CR3000Panel_temp_C)), AirTemp_Average_C))
  
  
  

#merge back air temp correction data
#Met<-rbind(Met_prefilter,Met_postfilter)

#put all airtemp through the panel QAQC check
#plot(Met_raw$DateTime, Met_raw$AirTemp_Average_C, type= 'l')
#points(Met$DateTime, Met$AirTemp_Average_C, type = "l", col="red")

#Air temp maximum set
Met=Met%>%
  mutate(
  Flag_AirTemp_Average_C=ifelse(AirTemp_Average_C>40.6 & !is.na(AirTemp_Average_C),4,Flag_AirTemp_Average_C),
  Note_AirTemp_Average_C=ifelse(AirTemp_Average_C>40.6 & is.na(AirTemp_Average_C),"Outlier_set_to_NA",Note_AirTemp_Average_C),
  AirTemp_Average_C=ifelse(AirTemp_Average_C>40.6 & is.na(AirTemp_Average_C),NA,AirTemp_Average_C))

#check the infared radiation-if there are any low values
#plot(Met_raw$DateTime, Met_raw$InfraredRadiationDown_Average_W_m2,main = "Raw Infared Radiation 2020", ylab = "Average_W_m2", type = "l")

#Infared radiation cleaning
#fix infrared radiation voltage reading after airtemp correction
#only need to do this for data from 2015 to  2016-07-25 10:12:00

Met=Met%>%
  mutate(  
    Note_InfraredRadiationDown_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationDown_Average_W_m2<250 & !is.na(InfraredRadiationDown_Average_W_m2)
                                                ,"Value_corrected_with_InfRadDn_equation_as_described_in_metadata",Note_InfraredRadiationDown_Average_W_m2),
  Flag_InfraredRadiationDown_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationDown_Average_W_m2<250 & !is.na(InfraredRadiationDown_Average_W_m2)
                                                ,4,Flag_InfraredRadiationDown_Average_W_m2),
    InfraredRadiationDown_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationDown_Average_W_m2<250 & !is.na(InfraredRadiationDown_Average_W_m2)
                                             ,5.67*10^-8*(AirTemp_Average_C+273.15)^4,InfraredRadiationDown_Average_W_m2),
  Note_InfraredRadiationUp_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationUp_Average_W_m2<100 & !is.na(InfraredRadiationUp_Average_W_m2)
                                              ,"Value_corrected_with_InfRadUp_equation_as_described_in_metadata",Note_InfraredRadiationUp_Average_W_m2),
 Flag_InfraredRadiationUp_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationUp_Average_W_m2<100 & is.na(InfraredRadiationUp_Average_W_m2)
                                                ,4,Flag_InfraredRadiationUp_Average_W_m2),
 InfraredRadiationUp_Average_W_m2=ifelse(DateTime<"2016-07-25 10:12:00" & InfraredRadiationUp_Average_W_m2<100 & !is.na(InfraredRadiationUp_Average_W_m2)
                                         ,InfraredRadiationUp_Average_W_m2+5.67*10^-8*(AirTemp_Average_C+273.15)^4,InfraredRadiationUp_Average_W_m2))

#Mean correction for InfRadDown (needs to be after voltage correction)
#Using 2018 data, taking the mean and sd of values on DOY to correct to
Met$DOY=yday(Met$DateTime)


Met_infrad=Met[year(Met$DateTime)<2018,]
Met_infrad$infradavg=ave(Met_infrad$InfraredRadiationDown_Average_W_m2, Met_infrad$DOY) #creating column with mean of infraddown by day of year
Met_infrad$infradsd=ave(Met_infrad$InfraredRadiationDown_Average_W_m2, Met_infrad$DOY, FUN = sd) #creating column with sd of infraddown by day of year
Met_infrad=unique(Met_infrad[,c(18,47,48)])

Met=merge(Met, Met_infrad, by = "DOY") #putting in columns for infrared mean and sd by DOY into main data set
Met=Met[order(Met$DateTime),] #ordering table after merging and removing unnecessary columns




#If the IR Down is greater than 3SD by DOY and replace with the equation from the manual(5.67*10-8(AirTemp_Average_C+273.15)^4)
Met=Met%>%
  mutate(
  Flag_InfraredRadiationDown_Average_W_m2=ifelse((abs(InfraredRadiationDown_Average_W_m2-infradavg)>(2*infradsd)) & !is.na(InfraredRadiationDown_Average_W_m2),4,Flag_InfraredRadiationDown_Average_W_m2),
  Note_InfraredRadiationDown_Average_W_m2=ifelse((abs(InfraredRadiationDown_Average_W_m2-infradavg)>(2*infradsd)) & !is.na(InfraredRadiationDown_Average_W_m2),"Greater_than_2SD_Value_corrected_with_InfRadDn_equation_as_described_in_metadata",Note_InfraredRadiationDown_Average_W_m2),
  InfraredRadiationDown_Average_W_m2=ifelse((abs(InfraredRadiationDown_Average_W_m2-infradavg)>(2*infradsd)) & !is.na(InfraredRadiationDown_Average_W_m2),5.67*10^-8*(AirTemp_Average_C+273.15)^4,InfraredRadiationDown_Average_W_m2))
  
Met=Met[,-c(1,47,48)]


#Inf outliers, must go after corrections
Met=Met%>%
  mutate(
    Flag_InfraredRadiationUp_Average_W_m2=ifelse(InfraredRadiationUp_Average_W_m2<150 & !is.na(InfraredRadiationUp_Average_W_m2) ,4,Flag_InfraredRadiationUp_Average_W_m2),
    Note_InfraredRadiationUp_Average_W_m2=ifelse(InfraredRadiationUp_Average_W_m2<150 & !is.na(InfraredRadiationUp_Average_W_m2),"Outlier_set_to_NA",Note_InfraredRadiationUp_Average_W_m2),
    InfraredRadiationUp_Average_W_m2=ifelse(InfraredRadiationUp_Average_W_m2<150 & !is.na(InfraredRadiationUp_Average_W_m2),NA,InfraredRadiationUp_Average_W_m2))

Met=Met%>%
  mutate(
    Flag_InfraredRadiationDown_Average_W_m2=ifelse(InfraredRadiationDown_Average_W_m2>540 & !is.na(InfraredRadiationDown_Average_W_m2)|InfraredRadiationDown_Average_W_m2<200 & !is.na(InfraredRadiationDown_Average_W_m2),4,Flag_InfraredRadiationDown_Average_W_m2),
    Note_InfraredRadiationDown_Average_W_m2=ifelse(InfraredRadiationDown_Average_W_m2>540 & !is.na(InfraredRadiationDown_Average_W_m2)|InfraredRadiationDown_Average_W_m2<200 & !is.na(InfraredRadiationDown_Average_W_m2),"outlier_set_to_NA_due_to_probable_sensor_failure",Note_InfraredRadiationDown_Average_W_m2),
    InfraredRadiationDown_Average_W_m2=ifelse(InfraredRadiationDown_Average_W_m2>540  & !is.na(InfraredRadiationDown_Average_W_m2) |InfraredRadiationDown_Average_W_m2<200 & !is.na(InfraredRadiationDown_Average_W_m2),NA,InfraredRadiationDown_Average_W_m2)
  )

#Replace missing values with estimations from the equation
Met=Met%>%
  mutate( 
    InfraredRadiationDown_Average_W_m2=ifelse(Flag_InfraredRadiationDown_Average_W_m2==2,5.67*10^-8*(AirTemp_Average_C+273.15)^4,InfraredRadiationDown_Average_W_m2),
    Note_InfraredRadiationDown_Average_W_m2=ifelse(Flag_InfraredRadiationDown_Average_W_m2==2,"Value_corrected_with_InfRadDn_equation_as_described_in_metadata",Note_InfraredRadiationDown_Average_W_m2),
    Flag_InfraredRadiationDown_Average_W_m2=ifelse(Flag_InfraredRadiationDown_Average_W_m2==2,4,Flag_InfraredRadiationDown_Average_W_m2),
    Note_InfraredRadiationUp_Average_W_m2=ifelse(Flag_InfraredRadiationUp_Average_W_m2==2,"Value_corrected_with_InfRadUp_equation_as_described_in_metadata",Note_InfraredRadiationUp_Average_W_m2),
    InfraredRadiationUp_Average_W_m2=ifelse(Flag_InfraredRadiationUp_Average_W_m2==2,InfraredRadiationUp_Average_W_m2+5.67*10^-8*(AirTemp_Average_C+273.15)^4,InfraredRadiationUp_Average_W_m2),
    Flag_InfraredRadiationUp_Average_W_m2=ifelse(Flag_InfraredRadiationUp_Average_W_m2==2,4,Flag_InfraredRadiationUp_Average_W_m2))

#take out impossible outliers and Infinite values before the other outliers are removed
#set flag 3 (see metadata: this corrects for impossible outliers)
for(i in 5:17) { #for loop to create new columns in data frame
  Met[c(which(is.infinite(Met[,i]))),paste0("Flag_",colnames(Met[i]))] <-3 #puts in flag 3
  Met[c(which(is.infinite(Met[,i]))),paste0("Note_",colnames(Met[i]))] <- "Infinite_value_set_to_NA" #note for flag 3
  Met[c(which(is.infinite(Met[,i]))),i] <- NA #set infinite vals to NA
  
  if(i!=8) { #flag 3 for negative values for everything except air temp
    Met[c(which((Met[,i]<0))),paste0("Flag_",colnames(Met[i]))] <- 3
    Met[c(which((Met[,i]<0))),paste0("Note_",colnames(Met[i]))] <- "Negative_value_set_to_0"
    Met[c(which((Met[,i]<0))),i] <- 0 #replaces value with 0
  }
  if(i==9) { #flag for RH over 100
    Met[c(which((Met[,i]>100))),paste0("Flag_",colnames(Met[i]))] <- 3
    Met[c(which((Met[,i]>100))),paste0("Note_",colnames(Met[i]))] <- "Value_set_to_100"
    Met[c(which((Met[,i]>100))),i] <- 100 #replaces value with 100
  }
}

#full data set QAQC


#Remove barometric pressure outliers
Met=Met%>%
  mutate(
  Flag_BP_Average_kPa=ifelse(BP_Average_kPa<98.5 & !is.na(BP_Average_kPa), 4,Flag_BP_Average_kPa),
  Note_BP_Average_kPa=ifelse(BP_Average_kPa<98.5 & !is.na(BP_Average_kPa),"Outlier_set_to_NA",Note_BP_Average_kPa),
  BP_Average_kPa=ifelse(BP_Average_kPa<98.5 & !is.na(BP_Average_kPa),NA,BP_Average_kPa))


#remove high PAR values at night
#get sunrise and sunset times
suntimes=getSunlightTimes(date = seq.Date(Sys.Date()-2500, Sys.Date(), by = 1),
                          keep = c("sunrise",  "sunset"),
                          lat = 37.30, lon = -79.83, tz = "Etc/GMT+5")

#create date column
Met$date <- as.Date(Met$DateTime)

#create subset to join


#now merge the datasets to get daylight time
Met <- left_join(Met, suntimes, by = "date") %>%
  mutate(daylight_intvl = interval(sunrise, sunset)) %>%
  mutate(during_day = DateTime %within% daylight_intvl)

#Remove PAR Tot
Met=Met%>%
  mutate(
    Flag_PAR_Total_mmol_m2=ifelse(during_day==FALSE & PAR_Total_mmol_m2>1 & !is.na(PAR_Total_mmol_m2) , 4, Flag_PAR_Total_mmol_m2),
    Note_PAR_Total_mmol_m2=ifelse(during_day==FALSE & PAR_Total_mmol_m2>1 & !is.na(PAR_Total_mmol_m2), "Outlier_set_to_NA", Note_PAR_Total_mmol_m2),
    PAR_Total_mmol_m2=ifelse(during_day==FALSE & PAR_Total_mmol_m2>1 & !is.na(PAR_Total_mmol_m2), NA, PAR_Total_mmol_m2))

Met=Met%>%
  mutate(
    Flag_PAR_Average_umol_s_m2=ifelse(during_day==FALSE & PAR_Average_umol_s_m2> 5 & !is.na(PAR_Average_umol_s_m2), 4, Flag_PAR_Average_umol_s_m2),
    Note_PAR_Average_umol_s_m2=ifelse(during_day==FALSE & PAR_Average_umol_s_m2>5 & !is.na(PAR_Average_umol_s_m2), "Outlier_set_to_NA", Note_PAR_Average_umol_s_m2),
    PAR_Average_umol_s_m2=ifelse(during_day==FALSE & PAR_Average_umol_s_m2>5 & !is.na(PAR_Average_umol_s_m2), NA, PAR_Average_umol_s_m2))

#Remove total PAR (PAR_Tot) outliers
Met=Met%>%
  mutate(
  Flag_PAR_Total_mmol_m2=ifelse(PAR_Total_mmol_m2>200& !is.na(PAR_Total_mmol_m2), 4, Flag_PAR_Total_mmol_m2),
  Flag_PAR_Average_umol_s_m2=ifelse(PAR_Average_umol_s_m2>3000 & !is.na(PAR_Average_umol_s_m2), 4, Flag_PAR_Average_umol_s_m2),
  Note_PAR_Total_mmol_m2=ifelse(PAR_Total_mmol_m2>200 & !is.na(PAR_Total_mmol_m2), "Outlier_set_to_NA", Note_PAR_Total_mmol_m2),
  Note_PAR_Average_umol_s_m2=ifelse(PAR_Average_umol_s_m2>3000 & !is.na(PAR_Average_umol_s_m2), "Outlier_set_to_NA", Note_PAR_Average_umol_s_m2),
  PAR_Total_mmol_m2=ifelse(PAR_Total_mmol_m2>200 & !is.na(PAR_Total_mmol_m2), NA, PAR_Total_mmol_m2),
  PAR_Average_umol_s_m2=ifelse(PAR_Average_umol_s_m2>3000 & !is.na(PAR_Average_umol_s_m2), NA, PAR_Average_umol_s_m2))

#Flag high total values for PAR Avg as over 3000
Met=Met%>%
  mutate(
    Flag_PAR_Average_umol_s_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Average_umol_s_m2>2500& !is.na(PAR_Average_umol_s_m2), 4, Flag_PAR_Average_umol_s_m2),
    Note_PAR_Average_umol_s_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Average_umol_s_m2>2500& !is.na(PAR_Average_umol_s_m2), "Outlier_set_to_NA", Note_PAR_Average_umol_s_m2),
    PAR_Average_umol_s_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Average_umol_s_m2>2500& !is.na(PAR_Average_umol_s_m2), NA, PAR_Average_umol_s_m2),
    Flag_PAR_Total_mmol_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Total_mmol_m2>150 & !is.na(PAR_Total_mmol_m2), 4, Flag_PAR_Total_mmol_m2),
    Note_PAR_Total_mmol_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Total_mmol_m2>150 & !is.na(PAR_Total_mmol_m2), "Outlier_set_to_NA", Note_PAR_Total_mmol_m2),
    PAR_Total_mmol_m2=ifelse(DateTime>"2021-05-01 00:00:00"&DateTime<"2021-07-01 23:59:00"& PAR_Total_mmol_m2>150 & !is.na(PAR_Total_mmol_m2), NA, PAR_Total_mmol_m2)
    )

#Remove shortwave radiation outliers
#first shortwave upwelling
Met=Met%>%
  mutate(
  Flag_ShortwaveRadiationUp_Average_W_m2=ifelse(ShortwaveRadiationUp_Average_W_m2>1500 & !is.na(ShortwaveRadiationUp_Average_W_m2), 4, Flag_ShortwaveRadiationUp_Average_W_m2),
  Note_ShortwaveRadiationUp_Average_W_m2=ifelse(ShortwaveRadiationUp_Average_W_m2>1500 & !is.na(ShortwaveRadiationUp_Average_W_m2), "Outlier_set_to_NA", Note_ShortwaveRadiationUp_Average_W_m2),
  ShortwaveRadiationUp_Average_W_m2=ifelse(ShortwaveRadiationUp_Average_W_m2>1500 & !is.na(ShortwaveRadiationUp_Average_W_m2), NA, ShortwaveRadiationUp_Average_W_m2))


#add a flag for suspect values from 1499 to 1300
#above 1500 already set to NA
Met=Met%>%
  mutate(
  Flag_ShortwaveRadiationUp_Average_W_m2=ifelse(ShortwaveRadiationUp_Average_W_m2>1300 & !is.na(ShortwaveRadiationUp_Average_W_m2), 5, Flag_ShortwaveRadiationUp_Average_W_m2),
  Note_ShortwaveRadiationUp_Average_W_m2=ifelse(ShortwaveRadiationUp_Average_W_m2>1300 & !is.na(ShortwaveRadiationUp_Average_W_m2), "Questionable_value", Note_ShortwaveRadiationUp_Average_W_m2))


#Eliminate the High Shortwave Radiation Up in 2018 before cleaning-put in maintenance log


#and then shortwave downwelling (what goes up must come down)
Met=Met%>%
  mutate(
  Flag_ShortwaveRadiationDown_Average_W_m2=ifelse(ShortwaveRadiationDown_Average_W_m2>300 & !is.na(ShortwaveRadiationDown_Average_W_m2), 4, Flag_ShortwaveRadiationDown_Average_W_m2),
  Note_ShortwaveRadiationDown_Average_W_m2=ifelse(ShortwaveRadiationDown_Average_W_m2>300 & !is.na(ShortwaveRadiationDown_Average_W_m2), "Outlier_set_to_NA", Note_ShortwaveRadiationDown_Average_W_m2),
  ShortwaveRadiationDown_Average_W_m2=ifelse(ShortwaveRadiationDown_Average_W_m2>300 & !is.na(ShortwaveRadiationDown_Average_W_m2), NA, ShortwaveRadiationDown_Average_W_m2))

#Eliminate the High Shortwave Radiation Up in 2018 before cleaning-put in maintenance log



#shortwave downwelling random point in 2016-put in maintenance log



#Remove albedo outliers
#over 1000
Met=Met%>%
  mutate(
  Flag_Albedo_Average_W_m2=ifelse(Albedo_Average_W_m2>1000 & !is.na(Albedo_Average_W_m2), 4, Flag_Albedo_Average_W_m2),
  Note_Albedo_Average_W_m2=ifelse(Albedo_Average_W_m2>1000 & !is.na(Albedo_Average_W_m2), "Outlier_set_to_NA", Note_Albedo_Average_W_m2),
  Albedo_Average_W_m2=ifelse(Albedo_Average_W_m2>1000 & !is.na(Albedo_Average_W_m2), NA, Albedo_Average_W_m2))

#set to NA when shortwave radiation up is equal to NA
Met=Met%>%
  mutate(
  Flag_Albedo_Average_W_m2=ifelse(is.na(ShortwaveRadiationUp_Average_W_m2)|is.na(ShortwaveRadiationDown_Average_W_m2), 4, Flag_Albedo_Average_W_m2),
  Note_Albedo_Average_W_m2=ifelse(is.na(ShortwaveRadiationUp_Average_W_m2)|is.na(ShortwaveRadiationDown_Average_W_m2), "Set_to_NA_because_Shortwave_equals_NA", Note_Albedo_Average_W_m2),
  Albedo_Average_W_m2=ifelse(is.na(ShortwaveRadiationUp_Average_W_m2)|is.na(ShortwaveRadiationDown_Average_W_m2), NA, Albedo_Average_W_m2))


#Reorder so the flags are all next to each other
Met=Met%>%
  select(-c(date))%>%
  select(c("Reservoir","Site","DateTime","Record","CR3000_Batt_V","CR3000Panel_temp_C","PAR_Average_umol_s_m2","PAR_Total_mmol_m2","BP_Average_kPa",                          
           "AirTemp_Average_C","RH_percent","Rain_Total_mm","WindSpeed_Average_m_s","WindDir_degrees","ShortwaveRadiationUp_Average_W_m2",       
           "ShortwaveRadiationDown_Average_W_m2","InfraredRadiationUp_Average_W_m2","InfraredRadiationDown_Average_W_m2","Albedo_Average_W_m2",
           "Flag_PAR_Average_umol_s_m2","Flag_PAR_Total_mmol_m2","Flag_BP_Average_kPa","Flag_AirTemp_Average_C","Flag_Rain_Total_mm","Flag_RH_percent",
           "Flag_WindSpeed_Average_m_s","Flag_WindDir_degrees","Flag_ShortwaveRadiationUp_Average_W_m2","Flag_ShortwaveRadiationDown_Average_W_m2",
           "Flag_InfraredRadiationUp_Average_W_m2","Flag_InfraredRadiationDown_Average_W_m2","Flag_Albedo_Average_W_m2","Note_PAR_Average_umol_s_m2",              
           "Note_PAR_Total_mmol_m2","Note_BP_Average_kPa","Note_AirTemp_Average_C","Note_RH_percent","Note_Rain_Total_mm","Note_WindSpeed_Average_m_s",              
           "Note_WindDir_degrees","Note_ShortwaveRadiationUp_Average_W_m2","Note_ShortwaveRadiationDown_Average_W_m2","Note_InfraredRadiationUp_Average_W_m2",       
           "Note_InfraredRadiationDown_Average_W_m2","Note_Albedo_Average_W_m2"))


###7) Write file with final cleaned dataset! ###
Met_final=Met%>%
  select(c("Reservoir","Site","DateTime","Record","CR3000_Batt_V","CR3000Panel_temp_C","PAR_Average_umol_s_m2","PAR_Total_mmol_m2","BP_Average_kPa",                          
           "AirTemp_Average_C","RH_percent","Rain_Total_mm","WindSpeed_Average_m_s","WindDir_degrees","ShortwaveRadiationUp_Average_W_m2",       
           "ShortwaveRadiationDown_Average_W_m2","InfraredRadiationUp_Average_W_m2","InfraredRadiationDown_Average_W_m2","Albedo_Average_W_m2",
           "Flag_PAR_Average_umol_s_m2","Note_PAR_Average_umol_s_m2","Flag_PAR_Total_mmol_m2","Note_PAR_Total_mmol_m2","Flag_BP_Average_kPa",
           "Note_BP_Average_kPa","Flag_AirTemp_Average_C","Note_AirTemp_Average_C","Flag_RH_percent","Note_RH_percent","Flag_Rain_Total_mm",
           "Note_Rain_Total_mm","Flag_WindSpeed_Average_m_s","Note_WindSpeed_Average_m_s",
           "Flag_WindDir_degrees","Note_WindDir_degrees","Flag_ShortwaveRadiationUp_Average_W_m2","Note_ShortwaveRadiationUp_Average_W_m2",
           "Flag_ShortwaveRadiationDown_Average_W_m2","Note_ShortwaveRadiationDown_Average_W_m2",
           "Flag_InfraredRadiationUp_Average_W_m2","Note_InfraredRadiationUp_Average_W_m2",
           "Flag_InfraredRadiationDown_Average_W_m2","Note_InfraredRadiationDown_Average_W_m2","Flag_Albedo_Average_W_m2","Note_Albedo_Average_W_m2"))

### Check with plots
ggplot(Met_final,mapping=aes(x=DateTime,y=PAR_Average_umol_s_m2))+
  geom_line()

## Save final Met data
write.csv(Met_final, "./Data/FCR_Met_final_2022.csv", row.names=F, quote=F)
