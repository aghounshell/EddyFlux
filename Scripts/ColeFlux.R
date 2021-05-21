### Calculating diffusive fluxes from dissolved GHGs using the Cole method
### Following McClure et al. 2018
### 21 May 2021, A Hounshell

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse,ncdf4,ggplot2,ggpubr,LakeMetabolizer,zoo,scales)

# First load in wind data from Met station at FCR ----
# For now, will need to download from GitHub


# Download 2020 Met data from EDI
# Downloaded from EDI: 21 May 2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431" 
#infile1 <- paste0(getwd(),"/Data/Met_final_2015_2020.csv")
#download.file(inUrl1,infile1,method="curl")

met_edi <- read.csv("./Data/Met_final_2015_2020.csv", header=T) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2019-12-31"))

# Load met data from 2021 (from GitHub, cleaned w/ script: MET_QAQC_2020.R)
met_21 <- read.csv("./Data/Met_GitHub_2021.csv", header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST")))

# Combine into one data frame
met_all <- rbind(met_edi,met_21)

# Select data every 30 minutes from Jan 2020 to end of met data
met_30 <- as.data.frame(seq(as.POSIXct("2020-01-01", tz = "EST"), as.POSIXct("2021-05-19", tz = "EST"), by = "30 min"))

met_30 <- met_30 %>% 
  rename(DateTime = `seq(as.POSIXct("2020-01-01", tz = "EST"), as.POSIXct("2021-05-19", tz = "EST"), by = "30 min")`)

met_30 <- left_join(met_30,met_all,by="DateTime")

