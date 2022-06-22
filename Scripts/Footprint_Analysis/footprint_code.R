###########################################################################

##### TO ESTIMATE FOOTPRINT
## From: Brenda on 8 June 2022

###########################################################################

## Clear workspace
rm(list = ls())

## Install necessary packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage", force = TRUE)

pacman::p_load(EBImage, spatialfil, raster, rasterVis, tidyverse, bigleaf,
               lubridate, rgdal)

## Set working directory
wd <- getwd()
setwd(wd)

#library(EBImage)
#library(spatialfil)
#library(raster)
#library(rasterVis)

#######################################

# change dates so that they match eddy dat

ts <- seq.POSIXt(as.POSIXct("2020-05-01 00:00:00",'%Y-%m-%d %H:%M:%S', tz = 'Etc/GMT+5'), 
                 as.POSIXct("2022-04-30 23:30:00",'%Y-%m-%d %H:%M:%S', tz = 'Etc/GMT+5'), by = "30 min")
ts2 <- data.frame(datetime = ts)

# ro is the name of the dataframe that has 
# the gapfilled eddy covariance data and that also includes
# L, sigma_v, Ustar, wind_dir, u, z_d_L which are part of the eddy output
# and then we need to estimate z0 (see below)

# Load in gap-filled eddy covariance data
ro <- read.csv("./Data/20220506_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-05-01 01:00:00") & DateTime < as.POSIXct("2022-05-01 01:00:00"))

###################################

# finding kinematic viscosity of air
# visc is the kinematic viscosity of the air (m2 s-1), given by (Massman 1999b) from bigleaf package

ro$visc <- bigleaf::kinematic.viscosity(Tair = ro$Tair, pressure = ro$pressure)

# z0 was estimated according to Zilitinkevich (1969) which combines
# charnock et al 1955 (ustar > 0.1), and Roll et al 1948 (ustar < 0.1)

z0_estimate <- function(c1, visc, ustar, g, c2) {
  z0 = c1*visc/ustar + ustar^2/(c2*g)
}

# estimate z0 in meters, c1 and c2 where obtained from Zilitinkevich 1969 and summarized in Foken

ro$z0 <- z0_estimate(c1 = 0.1, visc = ro$visc, ustar = ro$Ustar, g = 9.81, c2 = 20.8)


#######################################

ro <- as.data.frame(ro)

df2 <- ro %>% dplyr::select(DateTime, L, sigma_v, Ustar, wind_dir, u, z0, z_d_L)

df2$year <- year(ro$DateTime)
df2$month <- month(ro$DateTime)
df2$day <- day(ro$DateTime)
df2$HH_UTC <- hour(ro$DateTime)
df2$MM <- minute(ro$DateTime)

# remove -9999

df2[df2==-9999] <- NA


# modelling h

# f is coriolis = 0.881*10^-4 s for FCR's latitude

blh <- function(ustar, L) {
  f = 0.881*10^-4
  h = (L/3.8)*(-1 + (1 + 2.28*(ustar/(f*L)))^0.5) 
  return(h)
}

df2$h <- blh(ustar = df2$Ustar, L = df2$L)

df2$h <- ifelse(is.nan(df2$h), 1500, df2$h)

df2$h <- ifelse(df2$h <= 10, NA, df2$h)

df2$zm <- 2.9

df2$z_d_L <- ifelse(df2$z_d_L < -15.5, NA, df2$z_d_L)

df2 <- na.omit(df2)

# then get the center in easting and northing where the EC tower is
center = c(603066.34, 4129121.07) # center in easting (x), northing (y)

n_df <- nrow(df2)


# better to test with a small dataset because it takes a long time otherwise
df2_may2021 <- df2 %>% filter(DateTime <= '2021-05-15')



to_footprint_full <- df2 %>% dplyr::select(yyyy = year, 
                                    mm = month,
                                    day,
                                    HH_UTC, 
                                    MM,
                                    zm, 
                                    z0,
                                    L,
                                    umean = u,
                                    wind_dir,
                                    u_star = Ustar,
                                    sigma_v) %>% 
  mutate(d = 0)

to_footprint_may <- df2_may2021 %>% dplyr::select(yyyy = year, 
                                           mm = month,
                                           day,
                                           HH_UTC, 
                                           MM,
                                           zm, 
                                           z0,
                                           L,
                                           umean = u,
                                           wind_dir,
                                           u_star = Ustar,
                                           sigma_v) %>% 
  mutate(d = 0)


# write just in case
write_csv(to_footprint_full, './Data/full_footprint.csv')
write_csv(to_footprint_may, './Data/may_footprint.csv')


############################################################################

# reading the raster with water/vegetation
library(raster)
library(rasterVis)

fig_new <- raster("./Scripts/Footprint_Analysis/class_fcr.tif")
plot(fig_new)

#fig_new <- raster('~/Falling creek reservoir data/fcr_map_gee_modified.tif')

list(fig_new)
sr <- "+proj=utm +zone=17S +datum=NAD83"

projected_fig <- projectRaster(fig_new, crs = sr)

fig2 <- projected_fig

plot(fig2)

fig2[fig2 < 0] <- NA
fig2[fig2 > 100] <- NA

hist_fcr <- hist(fig2,
                 main = "Distribution of raster cell values",
                 xlab = "col", ylab = "Number of Pixels")

# creating reclassified matrix

reclass_df <- c(0.5, 1.5, 0,
                1.5, 3,1)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

fig4 <- reclassify(fig2, reclass_m)

plot(fig_new)


##########################################################################

# footprint climatology --> cumulative footprint for the whole study period
# call footprint function
library(EBImage)


# choose where the script to estimate the footprint is saved in your computer
source("./Scripts/Footprint_Analysis/calc_footprint_FFP_climatology.R")

FFPclim <- calc_footprint_FFP_climatology(zm = df2_may2021$zm, 
                                          z0 = df2_may2021$z0, 
                                          umean = df2_may2021$u, 
                                          h = df2_may2021$h, 
                                          ol = df2_may2021$L, 
                                          sigmav = df2_may2021$sigma_v, 
                                          ustar = df2_may2021$Ustar, 
                                          wind_dir = df2_may2021$wind_dir,
                                          domain = c(-500, 500, -500, 500))

##########################################################################

# contour lines

FFPclim$xr[[1]] <- c(FFPclim$xr[[1]]) + center[1]
FFPclim$yr[[1]] <- c(FFPclim$yr[[1]]) + center[2]
FFPclim$xr[[2]] <- c(FFPclim$xr[[2]]) + center[1]
FFPclim$yr[[2]] <- c(FFPclim$yr[[2]]) + center[2]
FFPclim$xr[[3]] <- c(FFPclim$xr[[3]]) + center[1]
FFPclim$yr[[3]] <- c(FFPclim$yr[[3]]) + center[2]
FFPclim$xr[[4]] <- c(FFPclim$xr[[4]]) + center[1]
FFPclim$yr[[4]] <- c(FFPclim$yr[[4]]) + center[2]
FFPclim$xr[[5]] <- c(FFPclim$xr[[5]]) + center[1]
FFPclim$yr[[5]] <- c(FFPclim$yr[[5]]) + center[2]
FFPclim$xr[[6]] <- c(FFPclim$xr[[6]]) + center[1]
FFPclim$yr[[6]] <- c(FFPclim$yr[[6]]) + center[2]
FFPclim$xr[[7]] <- c(FFPclim$xr[[7]]) + center[1]
FFPclim$yr[[7]] <- c(FFPclim$yr[[7]]) + center[2]
FFPclim$xr[[8]] <- c(FFPclim$xr[[8]]) + center[1]
FFPclim$yr[[8]] <- c(FFPclim$yr[[8]]) + center[2]
# create dataframe with latitude and longitude using the output from above  

mydfclim1 <- data.frame(
  lon = FFPclim$xr[[1]],
  lat = FFPclim$yr[[1]],
  ID = 1
)

mydfclim2 <- data.frame(
  lon = FFPclim$xr[[2]],
  lat = FFPclim$yr[[2]],
  ID = 2
)

mydfclim3 <- data.frame(
  lon = FFPclim$xr[[3]],
  lat = FFPclim$yr[[3]],
  ID = 3
)

mydfclim4 <- data.frame(
  lon = FFPclim$xr[[4]],
  lat = FFPclim$yr[[4]],
  ID = 4
)

mydfclim5 <- data.frame(
  lon = FFPclim$xr[[5]],
  lat = FFPclim$yr[[5]],
  ID = 5
)

mydfclim6 <- data.frame(
  lon = FFPclim$xr[[6]],
  lat = FFPclim$yr[[6]],
  ID = 6
)

mydfclim7 <- data.frame(
  lon = FFPclim$xr[[7]],
  lat = FFPclim$yr[[7]],
  ID = 7
)

mydfclim8 <- data.frame(
  lon = FFPclim$xr[[8]],
  lat = FFPclim$yr[[8]],
  ID = 8
)


r1 <- rbind(mydfclim1, mydfclim2)
r2 <- rbind(r1, mydfclim3)
r3 <- rbind(r2, mydfclim4, mydfclim5)
datalines <- rbind(r3, mydfclim6, mydfclim7, mydfclim8)

library(sp)

coordinates(datalines) <- ~lon+lat

x <- lapply(split(datalines, datalines$ID), function(x) Lines(list(Line(coordinates(x))), x$ID[1L]))

dat <- SpatialLines(x, CRS("+proj=utm +zone=17S +datum=NAD83"))



library(mapview)
mapview(dat, color = "red", legend = FALSE, 
        map.types = "Esri.WorldImagery", layer.name = 'Falling Creek reservoir')



# converting to FFP dataframe to polygon and converting the figure with classification to 
# spatial polygon to extract the values from classification figure
Piclim1 = Polygon(mydfclim1[,1:2])
Psiclim1 = SpatialPolygons(list(Polygons(list(Piclim1), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim2 = Polygon(mydfclim2[,1:2])
Psiclim2 = SpatialPolygons(list(Polygons(list(Piclim2), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim3 = Polygon(mydfclim3[,1:2])
Psiclim3 = SpatialPolygons(list(Polygons(list(Piclim3), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim4 = Polygon(mydfclim4[,1:2])
Psiclim4 = SpatialPolygons(list(Polygons(list(Piclim4), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim5 = Polygon(mydfclim5[,1:2])
Psiclim5 = SpatialPolygons(list(Polygons(list(Piclim5), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim6 = Polygon(mydfclim6[,1:2])
Psiclim6 = SpatialPolygons(list(Polygons(list(Piclim6), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim7 = Polygon(mydfclim7[,1:2])
Psiclim7 = SpatialPolygons(list(Polygons(list(Piclim7), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))
Piclim8 = Polygon(mydfclim8[,1:2])
Psiclim8 = SpatialPolygons(list(Polygons(list(Piclim8), ID = "a")), 
                           proj4string=CRS("+proj=utm +zone=17S +datum=NAD83"))





levelplot(fig4, margin=FALSE, col.regions=terrain.colors(3)) + 
  layer(sp.polygons(Psiclim1, col = 'black', lwd =2)) +
  layer(sp.polygons(Psiclim2, col = 'red', lwd =2)) +
  layer(sp.polygons(Psiclim3, col = 'blue', lwd =2)) +
  layer(sp.polygons(Psiclim4, col = 'green', lwd =2)) +
  layer(sp.polygons(Psiclim5, col = 'magenta', lwd =2)) +
  layer(sp.polygons(Psiclim6, col = 'turquoise', lwd =2)) +
  layer(sp.polygons(Psiclim7, col = 'red', lty =2, lwd =2)) +
  layer(sp.polygons(Psiclim8, col = 'black', lty =2, lwd =2))

levelplot(fig4, margin=FALSE, col.regions=terrain.colors(3))

vi7 <- raster::extract(fig4, Psiclim7) # extracting all the values from the classification fig.
vi8 <- raster::extract(fig4, Psiclim8) # extracting all the values from the classification fig.

# for me: 1 means vegetation, 0 means water
# estimate how much water and vegetation are in the 70% footprint isoline
output7  <-  data.frame(vi7)
total7 <- nrow(output7)
non_water7 <- sum(output7)
reservoir7 <- (total7 - non_water7)
non_water7 <- non_water7/total7*100
reservoir7 <- reservoir7/total7*100

# estimate how much water and vegetation are in the 80% footprint isoline

output8  <-  data.frame(vi8)
total8 <- nrow(output8)
non_water8 <- sum(output8)
reservoir8 <- (total8 - non_water8)
non_water8 <- non_water8/total8*100
reservoir8 <- reservoir8/total8*100

# end of footprint estimation


