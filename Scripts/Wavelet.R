### Script to conduct wavelet analysis on Eddy Flux Data
# Following Wavelets_for_fDOM_final_plots_lr_DWH.R
# A Hounshell, 16 July 2021

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr,lubridate)

# Load in Eddy Flux data (cleaned using FCR_Process_BD)
eddy_flux <- read_csv("./Data/20210615_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")) %>% 
  filter(DateTime >= as.POSIXct("2020-04-05 20:00:00") & DateTime < as.POSIXct("2021-04-05 20:00:00"))

# Aggregate to hourly (most common time step used for all analyses)
fcr_hourly <- eddy_flux %>% 
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

# Format for wavelet analysis: start with gap-filled data every hour
co2_data <- fcr_hourly %>% 
  select(DateTime,NEE) %>% 
  mutate(Time = 1:length(fcr_hourly$DateTime)) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

headers<-names(co2_data)
all<-headers[2]
temp<-matrix(-99,length(co2_data$Time),length(all),byrow=F)
head(temp)

temp <- scale(co2_data$NEE)

snt_co2_data<-as.data.frame(temp)
snt_co2_data<-setNames(snt_co2_data, all)#new dataframe with standard normal transformed data
head(snt_co2_data)

#### Setting up plotting function ####

#make new graphical function to fix period vs. scale issue
cols1<-c('blue3', 'blue', "dodgerblue3", "cyan", "green", "greenyellow", "yellow","orange","red", "red3") 

#run this so the COI countor appears for entire plot
options("max.contour.segments" = 250000)

wavelet.plot.new<-function (wave.list, wavelet.levels = quantile(wave.list$Power, 
                                                                 probs = seq(from = 0, to = 1, by = 0.1)), add.coi = TRUE, 
                            add.sig = TRUE, x.lab = gettext("Time"), period.lab = gettext("Period (days)"), 
                            #crn.lab = gettext("RWI"), 
                            key.cols = cols1, key.lab = parse(text = paste0("\"", gettext("Power"), 
                                                                            "\"^2")), add.spline = FALSE, f = 0.5, nyrs = NULL, crn.col = "black", 
                            crn.lwd = 1, crn.ylim = range(wave.list$y) * 1.1, 
                            side.by.side = FALSE) 
{
  y <- wave.list$y
  x <- wave.list$x
  x <- x/24  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi/24 # this fixes COI to correspond w/ editing period values so they = days 
  coi[coi == 0] <- 1e-12
  Power <- wave.list$Power 
  siglvl <- wave.list$siglvl
  if (any(diff(x) <= 0) || any(diff(period) <= 0)) {
    stop("'wave.list$x' and 'wave.list$period' must be strictly ascending")
  }
  if (period[1] <= 0) {
    stop("'wave.list$period' must be positive")
  }
  Signif <- t(matrix(Signif, dim(wave)[2], dim(wave)[1]))
  Signif <- Power/Signif
  period2 <- log2(period)
  ytick <- unique(trunc(period2))
  ytickv <- round ( (2^(ytick)) , digits = 2) #added round command to clean up y-axis so there wasn't 7 decimal points
  coi2 <- log2(coi)
  coi2[coi2 < 0] <- 0
  coi2.yy <- c(coi2, rep(max(period2, na.rm = TRUE), length(coi2)))
  coi2.yy[is.na(coi2.yy)] <- coi[2]
  yr.vec.xx <- c(x, rev(x))
  par.orig <- par(c("mar", "las", "mfrow"))
  on.exit(par(par.orig))
  nlevels <- length(wavelet.levels)
  seq.level <- seq_len(nlevels - 1)
  key.labs <- formatC(wavelet.levels, digits = 3, format = "f") #digits was 4, 3 gets last number, 2 gets second to last
  asp <- NA
  xaxs <- "i"
  yaxs <- "i"
  las <- 1
  xlim <- range(x, finite = TRUE)  
  ylim <- range(period2, finite = TRUE)
  z <- Power
  if (side.by.side) {
    layout(matrix(c(3, 2, 1), nrow = 1, byrow = TRUE), widths = c(1, 
                                                                  1, 0.2))
    #mar <- c(3, 1, 3, 3)
    mar <- c(3,3,3,3)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0), las = las)
    # plot.new()
    # plot.window(ylim = c(1, nlevels), xlim = c(0, 1), xaxs = xaxs, 
    #             yaxs = yaxs, asp = asp)
    # rect(0, seq.level, 1, 2:nlevels, col = key.cols)
    # axis(4, at = seq_along(wavelet.levels), labels = key.labs)
    # title(key.lab, cex.main = 1)
    mar <- c(3, 3, 3, 3)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0)) #***  was:tcl = 0.5, mgp = c(1.5, 0.25, 0)
    plot.new()  #controls bottom one 
    plot.window(xlim, ylim = c(0,5000), "", xaxs = xaxs, yaxs = yaxs, 
                asp = asp, las = las)
    .filled.contour(as.double(x), as.double(period2), z, 
                    as.double(wavelet.levels), key.cols)
    if (add.sig) {
      contour(x, period2, Signif, levels = 1, labels = siglvl, 
              drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, 
              add = TRUE, lwd = 2, col = "black", max.contour.segments > 25000) # still need to run max.contour.segments on line 81
    }
    if (add.coi) {
      polygon(yr.vec.xx, coi2.yy, density = c(10, 20), 
              angle = c(-45, 45), col = "black")
    }
    axis(1)
    #axis(3)
    axis(2, at = ytick, labels = ytickv)  #add in cex.lab to chagne font size ie. cex.lab = 1.2
    #axis(4, at = ytick, labels = ytickv)#ditto to abbove #
    title(xlab = x.lab, ylab = period.lab)
    box()
    mar <- c(3, 3, 3, 3)
    par(mar = mar, las = 0)
    plot(x, y, type = "l", xlim, ylim, xaxs = xaxs, yaxs = yaxs, 
         asp = asp, xlab = "", ylab = "", axes = FALSE, col = crn.col, 
         lwd = crn.lwd, ylim = crn.ylim, cex.lab = 1.3) # to try and increase font size 
    if (add.spline) {
      spl <- y
      tmp <- na.omit(spl)
      if (is.null(nyrs)) {
        nyrs2 <- length(tmp) * 0.33
      }
      else {
        nyrs2 <- nyrs
      }
      tmp <- ffcsaps(y = tmp, x = seq_along(tmp), nyrs = nyrs2, 
                     f = f)
      spl[!is.na(spl)] <- tmp
      lines(x, spl, col = "red", lwd = 2)
    }
    axis(1)
    #axis(3)
    axis(2)
    #axis(4)
    title(xlab = x.lab, ylab = crn.lab)
    box()
  }
  else {                                                              
    layout(matrix(c(2, 1), ncol = 1, byrow = TRUE), heights = c(1, 0.3)) # YES # Removed 3, from matrrix(c()), and 1, from hiehgts = c. This gets just wavelet and power grid on plot
    mar <- c(4,4,1,4)  #changed third number to 1 from 0.1. These values affect placement of plot (4414, 3315 worked best so far, all fit besides farright power values, past xxx5 doesnt help  )
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0), las = las) #chaning this isn't affecting power grid
    plot.new()
    plot.window(xlim = c(1, nlevels), ylim = c(0, 1), xaxs = xaxs, 
                yaxs = yaxs, asp = asp)
    rect(seq.level, 0, 2:nlevels, 1, col = key.cols) #change to shape power rectangle, was 2:nlevels
    axis(1, at = seq_along(wavelet.levels), labels = key.labs)
    title(sub = key.lab, cex.sub = 1, line = 1.5)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0))
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
                asp = asp, las = las)
    .filled.contour(as.double(x), as.double(period2), z, 
                    as.double(wavelet.levels), key.cols)
    if (add.sig) {
      contour(x, period2, Signif, levels = 1, labels = siglvl, 
              drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, 
              add = TRUE, lwd = 2, col = "black")
    }
    if (add.coi) {
      polygon(yr.vec.xx, coi2.yy, density = c(10, 20), 
              angle = c(-45, 45), col = "black")
    }
    axis(1)
    axis(2, at = ytick, labels = ytickv)
    #axis(3, labels = NA)
    #axis(4, at = ytick, labels = NA)
    title(xlab = x.lab, ylab = period.lab)
    # box()
    # mar <- c(0.1, 3, 3, 3)
    # par(mar = mar, las = 0)
    # plot(x, y, type = "l", xlim, xaxs = xaxs, yaxs = yaxs,
    #      asp = asp, xlab = "", ylab = "", axes = FALSE, col = crn.col, 
    #      lwd = crn.lwd, ylim = crn.ylim)
    if (add.spline) {
      spl <- y
      tmp <- na.omit(spl)
      if (is.null(nyrs)) {
        nyrs2 <- length(tmp) * 0.33
      }
      else {
        nyrs2 <- nyrs
      }
      tmp <- ffcsaps(y = tmp, x = seq_along(tmp), nyrs = nyrs2, 
                     f = f)
      spl[!is.na(spl)] <- tmp
      lines(x, spl, col = "red", lwd = 2)
    }
    axis(1, labels = NA)
    axis(2, labels = NA)
    #axis(3)
    #axis(4)
    #mtext(crn.lab, side = 4, line = 1.5, cex = 0.75)
    box()
  }
  invisible()
}

#### Getting Wavelet co2_output and plotting ####
Time<-co2_data$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_data<-cbind(Time, snt_co2_data)
head(snt_co2_data)
snt_co2_data <- snt_co2_data %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_data))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_data)

co2_output<-morlet(snt_co2_data$NEE, snt_co2_data$Time, dj=(1/12), siglvl = 0.95, p2= 12) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

co2_output$period <- round((co2_output$period/24), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
co2_dummy <- c(1:145)
for (j in 1:ncol(co2_output$Power)) {
  co2_dummy[j] <- mean(co2_output$Power[,j])
}

co2_powerplot <- as.data.frame(cbind(co2_dummy, co2_output$period))
co2_powerplot <- rename(co2_powerplot, c( "mean_power" = "co2_dummy"))
co2_powerplot <- rename(co2_powerplot, c( "period" = "V2"))

# Month mean power plot 
co2_powerplot <- co2_powerplot %>% 
  filter(period < 182)

###################
# Methane!
# Format for wavelet analysis: start with gap-filled averaged to hourly
ch4_data <- fcr_hourly %>% 
  select(DateTime,CH4) %>% 
  mutate(Time = 1:length(fcr_hourly$DateTime)) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

headers<-names(ch4_data)
all<-headers[2]
temp<-matrix(-99,length(ch4_data$Time),length(all),byrow=F)
head(temp)

temp <- scale(ch4_data$CH4)

snt_ch4_data<-as.data.frame(temp)
snt_ch4_data<-setNames(snt_ch4_data, all)#new dataframe with standard normal transformed data
head(snt_ch4_data)

#### Setting up plotting function ####

#make new graphical function to fix period vs. scale issue
cols1<-c('blue3', 'blue', "dodgerblue3", "cyan", "green", "greenyellow", "yellow","orange","red", "red3") 

#run this so the COI countor appears for entire plot
options("max.contour.segments" = 250000)

wavelet.plot.new<-function (wave.list, wavelet.levels = quantile(wave.list$Power, 
                                                                 probs = seq(from = 0, to = 1, by = 0.1)), add.coi = TRUE, 
                            add.sig = TRUE, x.lab = gettext("Time"), period.lab = gettext("Period (days)"), 
                            #crn.lab = gettext("RWI"), 
                            key.cols = cols1, key.lab = parse(text = paste0("\"", gettext("Power"), 
                                                                            "\"^2")), add.spline = FALSE, f = 0.5, nyrs = NULL, crn.col = "black", 
                            crn.lwd = 1, crn.ylim = range(wave.list$y) * 1.1, 
                            side.by.side = FALSE) 
{
  y <- wave.list$y
  x <- wave.list$x
  x <- x/24  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi/24 # this fixes COI to correspond w/ editing period values so they = days 
  coi[coi == 0] <- 1e-12
  Power <- wave.list$Power 
  siglvl <- wave.list$siglvl
  if (any(diff(x) <= 0) || any(diff(period) <= 0)) {
    stop("'wave.list$x' and 'wave.list$period' must be strictly ascending")
  }
  if (period[1] <= 0) {
    stop("'wave.list$period' must be positive")
  }
  Signif <- t(matrix(Signif, dim(wave)[2], dim(wave)[1]))
  Signif <- Power/Signif
  period2 <- log2(period)
  ytick <- unique(trunc(period2))
  ytickv <- round ( (2^(ytick)) , digits = 2) #added round command to clean up y-axis so there wasn't 7 decimal points
  coi2 <- log2(coi)
  coi2[coi2 < 0] <- 0
  coi2.yy <- c(coi2, rep(max(period2, na.rm = TRUE), length(coi2)))
  coi2.yy[is.na(coi2.yy)] <- coi[2]
  yr.vec.xx <- c(x, rev(x))
  par.orig <- par(c("mar", "las", "mfrow"))
  on.exit(par(par.orig))
  nlevels <- length(wavelet.levels)
  seq.level <- seq_len(nlevels - 1)
  key.labs <- formatC(wavelet.levels, digits = 3, format = "f") #digits was 4, 3 gets last number, 2 gets second to last
  asp <- NA
  xaxs <- "i"
  yaxs <- "i"
  las <- 1
  xlim <- range(x, finite = TRUE)  
  ylim <- range(period2, finite = TRUE)
  z <- Power
  if (side.by.side) {
    layout(matrix(c(3, 2, 1), nrow = 1, byrow = TRUE), widths = c(1, 
                                                                  1, 0.2))
    #mar <- c(3, 1, 3, 3)
    mar <- c(3,3,3,3)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0), las = las)
    # plot.new()
    # plot.window(ylim = c(1, nlevels), xlim = c(0, 1), xaxs = xaxs, 
    #             yaxs = yaxs, asp = asp)
    # rect(0, seq.level, 1, 2:nlevels, col = key.cols)
    # axis(4, at = seq_along(wavelet.levels), labels = key.labs)
    # title(key.lab, cex.main = 1)
    mar <- c(3, 3, 3, 3)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0)) #***  was:tcl = 0.5, mgp = c(1.5, 0.25, 0)
    plot.new()  #controls bottom one 
    plot.window(xlim, ylim = c(0,5000), "", xaxs = xaxs, yaxs = yaxs, 
                asp = asp, las = las)
    .filled.contour(as.double(x), as.double(period2), z, 
                    as.double(wavelet.levels), key.cols)
    if (add.sig) {
      contour(x, period2, Signif, levels = 1, labels = siglvl, 
              drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, 
              add = TRUE, lwd = 2, col = "black", max.contour.segments > 25000) # still need to run max.contour.segments on line 81
    }
    if (add.coi) {
      polygon(yr.vec.xx, coi2.yy, density = c(10, 20), 
              angle = c(-45, 45), col = "black")
    }
    axis(1)
    #axis(3)
    axis(2, at = ytick, labels = ytickv)  #add in cex.lab to chagne font size ie. cex.lab = 1.2
    #axis(4, at = ytick, labels = ytickv)#ditto to abbove #
    title(xlab = x.lab, ylab = period.lab)
    box()
    mar <- c(3, 3, 3, 3)
    par(mar = mar, las = 0)
    plot(x, y, type = "l", xlim, ylim, xaxs = xaxs, yaxs = yaxs, 
         asp = asp, xlab = "", ylab = "", axes = FALSE, col = crn.col, 
         lwd = crn.lwd, ylim = crn.ylim, cex.lab = 1.3) # to try and increase font size 
    if (add.spline) {
      spl <- y
      tmp <- na.omit(spl)
      if (is.null(nyrs)) {
        nyrs2 <- length(tmp) * 0.33
      }
      else {
        nyrs2 <- nyrs
      }
      tmp <- ffcsaps(y = tmp, x = seq_along(tmp), nyrs = nyrs2, 
                     f = f)
      spl[!is.na(spl)] <- tmp
      lines(x, spl, col = "red", lwd = 2)
    }
    axis(1)
    #axis(3)
    axis(2)
    #axis(4)
    title(xlab = x.lab, ylab = crn.lab)
    box()
  }
  else {                                                              
    layout(matrix(c(2, 1), ncol = 1, byrow = TRUE), heights = c(1, 0.3)) # YES # Removed 3, from matrrix(c()), and 1, from hiehgts = c. This gets just wavelet and power grid on plot
    mar <- c(4,4,1,4)  #changed third number to 1 from 0.1. These values affect placement of plot (4414, 3315 worked best so far, all fit besides farright power values, past xxx5 doesnt help  )
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0), las = las) #chaning this isn't affecting power grid
    plot.new()
    plot.window(xlim = c(1, nlevels), ylim = c(0, 1), xaxs = xaxs, 
                yaxs = yaxs, asp = asp)
    rect(seq.level, 0, 2:nlevels, 1, col = key.cols) #change to shape power rectangle, was 2:nlevels
    axis(1, at = seq_along(wavelet.levels), labels = key.labs)
    title(sub = key.lab, cex.sub = 1, line = 1.5)
    par(mar = mar, tcl = 0.5, mgp = c(1.5, 0.25, 0))
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
                asp = asp, las = las)
    .filled.contour(as.double(x), as.double(period2), z, 
                    as.double(wavelet.levels), key.cols)
    if (add.sig) {
      contour(x, period2, Signif, levels = 1, labels = siglvl, 
              drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, 
              add = TRUE, lwd = 2, col = "black")
    }
    if (add.coi) {
      polygon(yr.vec.xx, coi2.yy, density = c(10, 20), 
              angle = c(-45, 45), col = "black")
    }
    axis(1)
    axis(2, at = ytick, labels = ytickv)
    #axis(3, labels = NA)
    #axis(4, at = ytick, labels = NA)
    title(xlab = x.lab, ylab = period.lab)
    # box()
    # mar <- c(0.1, 3, 3, 3)
    # par(mar = mar, las = 0)
    # plot(x, y, type = "l", xlim, xaxs = xaxs, yaxs = yaxs,
    #      asp = asp, xlab = "", ylab = "", axes = FALSE, col = crn.col, 
    #      lwd = crn.lwd, ylim = crn.ylim)
    if (add.spline) {
      spl <- y
      tmp <- na.omit(spl)
      if (is.null(nyrs)) {
        nyrs2 <- length(tmp) * 0.33
      }
      else {
        nyrs2 <- nyrs
      }
      tmp <- ffcsaps(y = tmp, x = seq_along(tmp), nyrs = nyrs2, 
                     f = f)
      spl[!is.na(spl)] <- tmp
      lines(x, spl, col = "red", lwd = 2)
    }
    axis(1, labels = NA)
    axis(2, labels = NA)
    #axis(3)
    #axis(4)
    #mtext(crn.lab, side = 4, line = 1.5, cex = 0.75)
    box()
  }
  invisible()
}

#### Getting Wavelet ch4_output and plotting ####
Time<-ch4_data$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_data<-cbind(Time, snt_ch4_data)
head(snt_ch4_data)
snt_ch4_data <- snt_ch4_data %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_data))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_data)

ch4_output<-morlet(snt_ch4_data$CH4, snt_ch4_data$Time, dj=(1/12), siglvl = 0.95, p2= 12) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

ch4_output$period <- round((ch4_output$period/24), digits = 4)  #This works with correcting period to days for y axis 

#### Making mean global power plots and calculations ####
#calculating mean power per period 
ch4_dummy <- c(1:145)
for (j in 1:ncol(ch4_output$Power)) {
  ch4_dummy[j] <- mean(ch4_output$Power[,j])
}

ch4_powerplot <- as.data.frame(cbind(ch4_dummy, ch4_output$period))
ch4_powerplot <- rename(ch4_powerplot, c( "mean_power" = "ch4_dummy"))
ch4_powerplot <- rename(ch4_powerplot, c( "period" = "V2"))

# Month mean power plot 
ch4_powerplot <- ch4_powerplot %>% 
  filter(period < 182)

# Format graphs for MS
# Wavelet analysis
jpeg('./Fig_output/wavelet_co2_final.jpg',width=1050,height=650)
wavelet.plot.new(co2_output) #export in 977 x 701 for all numbers to show on color grid; 1000 x 650
dev.off()

jpeg('./Fig_output/wavelet_ch4_final.jpg',width=1050,height=650)
wavelet.plot.new(ch4_output)
dev.off()

jpeg('./Fig_output/Wavelet_both.jpg',width=1100,height = 1300)
ggarrange(wavelet.plot.new(co2_output),wavelet.plot.new(ch4_output),nrow=2,ncol=1,labels=c("A.","B."),
          font.label=list(face="plain",size=30))
dev.off()

# Global power spectra
rect1 <- data.frame (xmin = 0, xmax = 2, ymin=-Inf, ymax=Inf) #making rectangle to show range in daily plot

co2_monthlyPower <- ggplot(data = co2_powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 7,linetype="dashed")+
  annotate(geom="text",x = 4,y = 30,label = "7 d.")+
  geom_vline(xintercept = 14,linetype="dashed")+
  annotate(geom="text",x = 19,y = 30,label = "14 d.")+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 35,y = 30,label = "30 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  xlim(0,60)+
  ylim(0,30)+
  theme_classic(base_size = 15)

#Day mean power plot 
co2_powerplotday <- co2_powerplot %>% 
  filter(period <= 2)

co2_dailyPower <- ggplot(data = co2_powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.15,y = 7,label = "1 d.")+
  geom_vline(xintercept = 0.5,linetype="dashed")+
  annotate(geom="text",x = 0.7,y = 7,label = "12 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  theme_classic(base_size = 15)

ch4_monthlyPower <- ggplot(data = ch4_powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 7,linetype="dashed")+
  annotate(geom="text",x = 4,y = 60,label = "7 d.")+
  geom_vline(xintercept = 14,linetype="dashed")+
  annotate(geom="text",x = 19,y = 60,label = "14 d.")+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 35,y = 60,label = "30 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  xlim(0,60)+
  ylim(0,60)+
  theme_classic(base_size = 15)

#Day mean power plot 
ch4_powerplotday <- ch4_powerplot %>% 
  filter(period <= 2)

ch4_dailyPower <- ggplot(data = ch4_powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.15,y = 5,label = "1 d.")+
  geom_vline(xintercept = 0.5,linetype="dashed")+
  annotate(geom="text",x = 0.7,y = 5,label = "12 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  theme_classic(base_size = 15)

ggarrange(co2_monthlyPower,co2_dailyPower,ch4_monthlyPower,ch4_dailyPower,ncol=2,nrow=2,labels=c("A.","B.","C.","D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/GlobalPower_All.jpg",width=9,height=7,units="in",dpi=320)
