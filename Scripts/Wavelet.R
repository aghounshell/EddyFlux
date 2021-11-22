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

# Format for wavelet analysis: Following M. Johnson's comment, conduct wavelet on half-hourly data
co2_data <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f) %>% 
  mutate(Time = 1:length(eddy_flux$DateTime)) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

headers<-names(co2_data)
all<-headers[2]
temp<-matrix(-99,length(co2_data$Time),length(all),byrow=F)
head(temp)

temp <- scale(co2_data$NEE_uStar_f)

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
  x <- x/48  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi/48 # this fixes COI to correspond w/ editing period values so they = days 
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

co2_output<-morlet(snt_co2_data$NEE_uStar_f, snt_co2_data$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

co2_output$period <- round((co2_output$period/48), digits = 4)  #This works with correcting period to days for y axis 

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
ch4_data <- eddy_flux %>% 
  select(DateTime,ch4_flux_uStar_f) %>% 
  mutate(Time = 1:length(eddy_flux$DateTime)) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

headers<-names(ch4_data)
all<-headers[2]
temp<-matrix(-99,length(ch4_data$Time),length(all),byrow=F)
head(temp)

temp <- scale(ch4_data$ch4_flux_uStar_f)

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
  x <- x/48  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi/48 # this fixes COI to correspond w/ editing period values so they = days 
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

ch4_output<-morlet(snt_ch4_data$ch4_flux_uStar_f, snt_ch4_data$Time, dj=(1/12), siglvl = 0.95, p2= 13) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

ch4_output$period <- round((ch4_output$period/48), digits = 4)  #This works with correcting period to days for y axis 

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
jpeg('./Fig_output/wavelet_co2_halfhour.jpg',width=1050,height=650)
wavelet.plot.new(co2_output) #export in 977 x 701 for all numbers to show on color grid; 1000 x 650
dev.off()

jpeg('./Fig_output/wavelet_ch4_halfhour.jpg',width=1050,height=650)
wavelet.plot.new(ch4_output)
dev.off()

jpeg('./Fig_output/Wavelet_both_halfhour.jpg',width=1100,height = 1300)
ggarrange(wavelet.plot.new(co2_output),wavelet.plot.new(ch4_output),nrow=2,ncol=1,labels=c("A.","B."),
          font.label=list(face="plain",size=30))
dev.off()

# Global power spectra
rect1 <- data.frame (xmin = 0, xmax = 2, ymin=-Inf, ymax=Inf) #making rectangle to show range in daily plot

co2_monthlyPower <- ggplot(data = co2_powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 7,linetype="dashed")+
  annotate(geom="text",x = 4,y = 45,label = "7 d.")+
  geom_vline(xintercept = 14,linetype="dashed")+
  annotate(geom="text",x = 19,y = 45,label = "14 d.")+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 35,y = 45,label = "30 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  xlim(0,60)+
  ylim(0,45)+
  theme_classic(base_size = 15)

#Day mean power plot 
co2_powerplotday <- co2_powerplot %>% 
  filter(period <= 2)

co2_dailyPower <- ggplot(data = co2_powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.15,y = 10,label = "1 d.")+
  geom_vline(xintercept = 0.5,linetype="dashed")+
  annotate(geom="text",x = 0.7,y = 10,label = "12 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  theme_classic(base_size = 15)

ch4_monthlyPower <- ggplot(data = ch4_powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 7,linetype="dashed")+
  annotate(geom="text",x = 4,y = 100,label = "7 d.")+
  geom_vline(xintercept = 14,linetype="dashed")+
  annotate(geom="text",x = 19,y = 100,label = "14 d.")+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 35,y = 100,label = "30 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  xlim(0,60)+
  ylim(0,100)+
  theme_classic(base_size = 15)

#Day mean power plot 
ch4_powerplotday <- ch4_powerplot %>% 
  filter(period <= 2)

ch4_dailyPower <- ggplot(data = ch4_powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.15,y = 7,label = "1 d.")+
  geom_vline(xintercept = 0.5,linetype="dashed")+
  annotate(geom="text",x = 0.7,y = 7,label = "12 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  theme_classic(base_size = 15)

ggarrange(co2_monthlyPower,co2_dailyPower,ch4_monthlyPower,ch4_dailyPower,ncol=2,nrow=2,labels=c("A.","B.","C.","D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/GlobalPower_All_halfhourly.jpg",width=9,height=7,units="in",dpi=320)

### Check for days when daily or diel time points are significant ----
#taking calcs from wavelet.plot.new function
signif_test <- t(matrix(co2_output$Signif, dim(co2_output$wave)[2], dim(co2_output$wave)[1]))
signif_testa <- co2_output$Power/signif_test            #puts signif in 52549 x 16 matrix 

a <- signif_testa 
b <- co2_output$Power

# Check time periods
time_period <- as.data.frame(round((2^(log2(co2_output$period))),digits=2))

### YES!! This works to assign 1 or 0 to true or false values at one scale 
c <- as.data.frame(b >= signif_test)  #this tells us if the power value > signif value 

#c7 gets T/F values for if signif at daily timescale 
c_day <- c %>%  
  select(V42:V43,V54:V57)
c_day <- as.vector(c_day)

#if else gives us numerics for logistic regression. 1 = T, 0 = F
coutput_day <- ifelse(c_day$V54 == TRUE & c_day$V55 == TRUE & c_day$V56 == TRUE & c_day$V57 == TRUE, 1, 0)

coutput_diel <- ifelse(c_day$V42 == TRUE & c_day$V43 == TRUE, 1,0)

coutput <- cbind(co2_data,coutput_day,coutput_diel)

coutput <- coutput %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour == 12)

coutput_sum <- coutput %>% 
  summarise(day_sum = sum(coutput_day),
            diel_sum = sum(coutput_diel))
# Day = 105 days
# Diel = 57 days - from 5-23 to 11-08

# Do same thing for CH4 data
#taking calcs from wavelet.plot.new function
signif_test_ch4 <- t(matrix(ch4_output$Signif, dim(ch4_output$wave)[2], dim(ch4_output$wave)[1]))
signif_testa_ch4 <- ch4_output$Power/signif_test_ch4            #puts signif in 52549 x 16 matrix 

a_ch4 <- signif_testa_ch4 
b_ch4 <- ch4_output$Power

### YES!! This works to assign 1 or 0 to true or false values at one scale 
c_ch4 <- as.data.frame(b_ch4 >= signif_test_ch4)  #this tells us if the power value > signif value 

#c7 gets T/F values for if signif at daily timescale 
c_day_ch4 <- c_ch4 %>%  
  select(V42:V43,V54:V57)
c_day_ch4 <- as.vector(c_day_ch4)

#if else gives us numerics for logistic regression. 1 = T, 0 = F
coutput_day_ch4 <- ifelse(c_day_ch4$V54 == TRUE & c_day_ch4$V55 == TRUE & c_day_ch4$V56 == TRUE & c_day_ch4$V57 == TRUE, 1, 0)

coutput_diel_ch4 <- ifelse(c_day_ch4$V42 == TRUE & c_day_ch4$V43 == TRUE, 1,0)

coutput_ch4 <- cbind(ch4_data,coutput_day_ch4,coutput_diel_ch4)

# Select only days with noon
coutput_ch4 <- coutput_ch4 %>% 
  mutate(hour = hour(DateTime)) %>% 
  filter(hour == 12)

coutput_sum_ch4 <- coutput_ch4 %>% 
  summarise(day_sum = sum(coutput_day_ch4),
            diel_sum = sum(coutput_diel_ch4))
# Day = 89 days
# Diel = 64 days

# Plot significant timescales by date
day_sig <- ggplot()+
  annotate("text",x=as.POSIXct("2021-04-05"),y=1.05,label="Sig.",size=5)+
  annotate("text",x=as.POSIXct("2021-04-01"),y=0.2,label="Not Sig.",size=5)+
  geom_point(coutput,mapping=aes(x=DateTime,y=coutput_day,color="co2"),alpha=0.5)+
  geom_point(coutput_ch4,mapping=aes(x=DateTime,y=coutput_day_ch4+0.1,color="ch4"),alpha=0.5)+
  ylab("Daily Significance")+
  xlab("")+
  scale_color_manual(breaks=c("co2","ch4"), labels=c((expression(~CO[2])),expression(~CH[4])),
                     values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title=element_blank())

diel_sig <- ggplot()+
  annotate("text",x=as.POSIXct("2021-04-05"),y=1.05,label="Sig.",size=5)+
  annotate("text",x=as.POSIXct("2021-04-01"),y=0.2,label="Not Sig.",size=5)+
  geom_point(coutput,mapping=aes(x=DateTime,y=coutput_diel,color="co2"),alpha=0.5)+
  geom_point(coutput_ch4,mapping=aes(x=DateTime,y=coutput_diel_ch4+0.1,color="ch4"),alpha=0.5)+
  ylab("Diel Significance")+
  xlab("")+
  scale_color_manual(breaks=c("co2","ch4"), labels=c((expression(~CO[2])),expression(~CH[4])),
                     values=c("#E63946","#4c8bfe"))+
  theme_classic(base_size = 15)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title=element_blank())

ggarrange(day_sig,diel_sig,nrow=2,ncol=1,common.legend=TRUE,
          labels=c("A.","B."),font.label = list(face="plain",size=15))

ggsave("./Fig_Output/SI_Timescale_Sig_halfhourly.jpg",width=8,height=6,units="in",dpi=320)

##### Plot diel comparisons for summer stratified period (June, July) ######
# Separate into day (0900 to 1500) and night (2100 to 0300)
daynight <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,NEE_uStar_orig,ch4_flux_uStar_f,ch4_flux_uStar_orig,u) %>% 
  mutate(Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 9 & Hour <= 14, "Day",
                       ifelse(Hour >= 21 | Hour <= 2, "Night", NA))) %>% 
  filter(diel == "Day" | diel == "Night") %>% 
  mutate(DateTime = format(as.POSIXct(DateTime,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d'))

# Separate into max/min emissions: max = 0600 and min = 1800 (following Shao et al. 2015)
maxmin <- eddy_flux %>% 
  select(DateTime,NEE_uStar_f,NEE_uStar_orig,ch4_flux_uStar_f,ch4_flux_uStar_orig,u) %>% 
  mutate(Hour = hour(DateTime)) %>% 
  mutate(diel = ifelse(Hour >= 5 & Hour <= 7, "Min",
                       ifelse(Hour >= 17 & Hour <= 19, "Max", NA))) %>% 
  filter(diel == "Min" | diel == "Max") %>% 
  mutate(DateTime = format(as.POSIXct(DateTime,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d'))

# Calculate median for each day/night time period
# Add seasonal designations
daynight_med <- daynight %>% 
  group_by(DateTime,diel) %>% 
  summarise_all(median,na.rm=TRUE) %>% 
  mutate(Month = month(DateTime)) %>% 
  mutate(season = ifelse(DateTime < as.POSIXct("2020-06-01"), "Spring",
                         ifelse(DateTime >= as.POSIXct("2020-06-01") & DateTime < as.POSIXct("2020-09-01"), "Summer",
                                ifelse(DateTime >= as.POSIXct("2020-09-01") & DateTime < as.POSIXct("2020-11-01"), "Fall",
                                       ifelse(DateTime >= as.POSIXct("2020-11-01") & DateTime < as.POSIXct("2021-03-01"), "Winter",
                                              ifelse(DateTime >= as.POSIXct("2021-03-01"), "Spring", NA))))))

maxmin_med <- maxmin %>% 
  group_by(DateTime,diel) %>% 
  summarise_all(median,na.rm=TRUE) %>% 
  mutate(Month = month(DateTime)) %>% 
  mutate(season = ifelse(DateTime < as.POSIXct("2020-06-01"), "Spring",
                         ifelse(DateTime >= as.POSIXct("2020-06-01") & DateTime < as.POSIXct("2020-09-01"), "Summer",
                                ifelse(DateTime >= as.POSIXct("2020-09-01") & DateTime < as.POSIXct("2020-11-01"), "Fall",
                                       ifelse(DateTime >= as.POSIXct("2020-11-01") & DateTime < as.POSIXct("2021-03-01"), "Winter",
                                              ifelse(DateTime >= as.POSIXct("2021-03-01"), "Spring", NA))))))

# Plot Max-Min
maxmin_med %>% 
  ggplot(mapping=aes(x=season,y=NEE_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Min','Max'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15) 

maxmin_med %>% 
  ggplot(mapping=aes(x=diel,y=ch4_flux_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Min','Max'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 


# Plot
ggplot(daynight,mapping=aes(x=DateTime,y=NEE_uStar_f,color=diel))+
  geom_line()+
  geom_point()+
  theme_classic(base_size = 15)

# Conduct statistical analysis - Wilcoxon signed rank tests
# Need to first pivot_wide for paired analyses
daynight_wide_co2_f <- daynight_med %>% 
  select(DateTime,diel,NEE_uStar_f) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_f)

wilcox.test(daynight_wide_co2_f$Day,daynight_wide_co2_f$Night,paired=TRUE)

daynight_wide_ch4_f <- daynight_med %>% 
  select(DateTime,diel,ch4_flux_uStar_f) %>% 
  pivot_wider(names_from = diel, values_from = ch4_flux_uStar_f)

wilcox.test(daynight_wide_ch4_f$Day,daynight_wide_ch4_f$Night,paired=TRUE)

# Plot?!
# CO2 measured
daynight_wide_co2_orig <- daynight_med %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

wilcox.test(daynight_wide_co2_orig$Day,daynight_wide_co2_orig$Night,paired=TRUE)

co2_diel <- daynight_med %>% 
  ggplot(mapping=aes(x=diel,y=NEE_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label ="p = 0.03*
n = 225",
            x = "Day", hjust = 0, y = 25, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 

daynight_wide_ch4_orig <- daynight_med %>% 
  select(DateTime,diel,ch4_flux_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = ch4_flux_uStar_orig) %>% 
  drop_na

wilcox.test(daynight_wide_ch4_orig$Day,daynight_wide_ch4_orig$Night,paired=TRUE)

ch4_diel <- daynight_med %>% 
  ggplot(mapping=aes(x=diel,y=ch4_flux_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  annotate("text",label = "p = 0.28
n = 165",
           x = "Day", hjust = 0, y = 0.07, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 

# Include wind speed?
daynight_wide_u <- daynight_med %>% 
  select(DateTime,diel,u) %>% 
  pivot_wider(names_from = diel, values_from = u) %>% 
  drop_na()

wilcox.test(daynight_wide_u$Day,daynight_wide_u$Night,paired=TRUE)

wind <- daynight_med %>% 
  ggplot(mapping=aes(x=diel,y=u,color=diel))+
  annotate("text",label = "p < 0.005*
n = 282",
           x = "Day", hjust = 0, y = 4, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("Wind speed (m s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.position = "none") 

ggarrange(co2_diel,ch4_diel,wind,nrow=1,ncol=3,
          labels=c("A.","B.","C."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/summer_diel_wind.jpg",width = 9, height=4, units="in",dpi=320)

# Plot by month
co2_month <- daynight_med %>% 
  ggplot(mapping=aes(x=as.factor(Month),y=NEE_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  #annotate("text",label = "*",x = "Spring", y = 30, size = 5)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CO"[2]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank()) 

ch4_month <- daynight_med %>% 
  ggplot(mapping=aes(x=as.factor(Month),y=ch4_flux_uStar_orig,color=diel))+
  geom_hline(yintercept = 0, linetype = "dashed", color="darkgrey", size = 0.8)+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("CH"[4]*" (",mu,"mol C m"^-2*" s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank()) 

u_month <- daynight_med %>% 
  ggplot(mapping=aes(x=as.factor(Month),y=u,color=diel))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_point(position=position_jitterdodge(),alpha=0.1)+
  scale_color_manual(breaks=c('Day','Night'),values=c("#F4BB01","#183662"))+
  ylab(expression(paste("Wind speed (m s"^-1*")")))+
  xlab("")+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank()) 

ggarrange(co2_month,ch4_month,u_month,nrow=3,ncol=1,common.legend = TRUE,
          labels=c("A.","B.","C."), font.label = list(face="plain",size=15))

ggsave("./Fig_Output/month_diel.jpg",width = 9, height=10, units="in",dpi=320)

temp <- daynight_med %>% 
  filter(Month == 1) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 2) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 3) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 4) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 5) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 6) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 7) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 8) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 9) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 10) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 11) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

temp <- daynight_med %>% 
  filter(Month == 12) %>% 
  select(DateTime,diel,NEE_uStar_orig) %>% 
  pivot_wider(names_from = diel, values_from = NEE_uStar_orig) %>% 
  drop_na()

