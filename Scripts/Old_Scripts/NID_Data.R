### Looking at NID data for dam stats
### A Hounshell, 28 June 2021

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(tidyverse)

# Load in data
nid <- read.csv("./Data/NID2019_U.csv", header=T)

# Select reservoirs with surface area < 1 km2
nid_small <- nid %>% 
  mutate(SURFACE_AREA_km2 = SURFACE_AREA*0.00404686) %>% 
  filter(SURFACE_AREA_km2 <= 1)

# Select old reservoirs
nid_small_old <- nid_small %>% 
  filter(YEAR_COMPLETED <= 1921)
