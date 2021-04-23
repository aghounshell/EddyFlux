### Script to check raw data for CO2/CH4 troubleshooting
### A Hounshell, 23 Apr 2021

setwd("C:/Users/ahoun/OneDrive/Desktop/EddyFlux")

# Load in packages
pacman::p_load(tidyverse,ggplot2,ggpubr)

# Load in data
files <- list.files("./RawData/Mar_2021", full.names = TRUE)
files
length(files)

df <- read.delim(files[1],header=TRUE,sep="\t")[-1,]
head(df)
tail(df)

for( i in files) {
  # Reading in file ----
  message("Starting ", i)
  df <- read.delim(i,header=TRUE,sep="\t")[-1,]
  if(i == files[1]) {
    df2 <- df
  } else {
    df2 <- rbind(df2, df)
  }
  # Finished loop ----
}

# Strange - files starting on 2021-04-06 have two extra columns (air_t_mean;	co2_signal_strength_7500_mean
# ) at the end of the file? Load in two separate data columns to see
files <- list.files("./RawData/Apr_2021", full.names = TRUE)
files
length(files)

df3 <- read.delim(files[1],header=TRUE,sep="\t")[-1,]
head(df3)
tail(df3)

for( i in files) {
  # Reading in file ----
  message("Starting ", i)
  df3 <- read.delim(i,header=TRUE,sep="\t")[-1,]
  if(i == files[1]) {
    df4 <- df3
  } else {
    df4 <- rbind(df4, df3)
  }
  # Finished loop ----
}

# But then goes back to 'normal' on 8 Apr 2021
files <- list.files("./RawData/Apr_2021_a", full.names = TRUE)
files
length(files)

df5 <- read.delim(files[1],header=TRUE,sep="\t")[-1,]
head(df5)
tail(df5)

for( i in files) {
  # Reading in file ----
  message("Starting ", i)
  df5 <- read.delim(i,header=TRUE,sep="\t")[-1,]
  if(i == files[1]) {
    df6 <- df5
  } else {
    df6 <- rbind(df6, df5)
  }
  # Finished loop ----
}

# df2 = pre 05 Apr 2021; df4 = 05 Apr 2021 - 07 Apr 2021; df6 = 08 Apr 2021 - 16 Apr 2021
# Concatenate full data frame (minus the two 'extra' files)

df4_a <- df4 %>% 
  select(-c("air_t_mean","co2_signal_strength_7500_mean"))

full <- rbind(df2,df4_a,df6)

