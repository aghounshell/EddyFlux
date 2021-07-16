### Script to conduct wavelet analysis on Eddy Flux Data
# Following Wavelets_for_fDOM_final_plots_lr_DWH.R
# A Hounshell, 16 July 2021

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load in libraries
pacman::p_load(zoo,dplR,dplyr,tidyverse,ggplot2,ggpubr)

# Load in Eddy Flux data (cleaned using FCR_Process_BD)
eddy_flux <- read_csv("./Data/20210615_EC_processed.csv") %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST"))

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

# Format for wavelet analysis: start with gap-filled data every 30 minutes?
co2_data <- fcr_hourly %>% 
  select(DateTime,NEE) %>% 
  mutate(Time = 1:length(fcr_hourly$DateTime)) %>% 
  filter(DateTime <= as.POSIXct("2021-04-05 16:00:00"))

headers<-names(co2_data)
all<-headers[2]
temp<-matrix(-99,length(co2_data$Time),length(all),byrow=F)
head(temp)

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(co2_data)) {
  
  temp[,j] <- (co2_data$NEE - mean(co2_data$NEE, na.rm = TRUE)) / sd(co2_data$NEE, na.rm = TRUE)
  
}

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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

#### Getting Wavelet output and plotting ####

Time<-co2_data$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_data<-cbind(Time, snt_co2_data)
head(snt_co2_data)
snt_co2_data <- snt_co2_data %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_data))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_data)

output<-morlet(snt_co2_data$NEE, snt_co2_data$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

jpeg('./Fig_Output/Wavelet_co2.jpg')
wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 
dev.off()

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
View(output$Power)
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}

View(dummy)

powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot <- powerplot %>% 
  filter(period < 182)

rect1 <- data.frame (xmin = 0, xmax = 7, ymin=-Inf, ymax=Inf) #making rectangle to show range in daily plot

monthlyPower <- ggplot(data = powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 30,y = 400,label = "30 d.")+
  geom_vline(xintercept = 60,linetype="dashed")+
  annotate(geom="text",x = 65,y = 400,label = "60 d.")+
  geom_vline(xintercept = 120,linetype="dashed")+
  annotate(geom="text",x = 120,y = 400,label = "120 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  mytheme_AS
monthlyPower

#Day mean power plot 
powerplotday <- powerplot %>% 
  filter(period < 7)

dailyPower <- ggplot(data = powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 25,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 25,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
dailyPower

ggarrange(ggarrange(co2_wave,ncol=1,labels = c("A."),font.label = list(face="plain",size=15)),
          ggarrange(monthlyPower,dailyPower,ncol=2,labels = c("B.","C."),font.label=list(face="plain",size=15)),
          nrow=2)

ggsave("./Fig_Output/Co2_wave.jpg",width=10,height=8,units="in",dpi=320)


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

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(ch4_data)) {
  
  temp[,j] <- (ch4_data$CH4 - mean(ch4_data$CH4, na.rm = TRUE)) / sd(ch4_data$CH4, na.rm = TRUE)
  
}

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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



#### Getting Wavelet output and plotting ####

Time<-ch4_data$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_data<-cbind(Time, snt_ch4_data)
head(snt_ch4_data)
snt_ch4_data <- snt_ch4_data %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_data))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_data)

output<-morlet(snt_ch4_data$CH4, snt_ch4_data$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
View(output$Power)
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}

View(dummy)

powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot <- powerplot %>% 
  filter(period < 182)

rect1 <- data.frame (xmin = 0, xmax = 7, ymin=-Inf, ymax=Inf) #making rectangle to show range in daily plot

monthlyPower <- ggplot(data = powerplot, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 30,linetype="dashed")+
  annotate(geom="text",x = 30,y = 700,label = "30 d.")+
  geom_vline(xintercept = 60,linetype="dashed")+
  annotate(geom="text",x = 65,y = 700,label = "60 d.")+
  geom_vline(xintercept = 120,linetype="dashed")+
  annotate(geom="text",x = 120,y = 700,label = "120 d.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  mytheme_AS
monthlyPower

#Day mean power plot 
powerplotday <- powerplot %>% 
  filter(period < 8)

dailyPower <- ggplot(data = powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_vline(xintercept = 7,linetype="dashed")+
  annotate(geom="text",x = 7,y = 60,label = "7 d.")+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 60,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 60,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
dailyPower

ggarrange(ggarrange(co2_wave,ncol=1,labels = c("A."),font.label = list(face="plain",size=15)),
          ggarrange(monthlyPower,dailyPower,ncol=2,labels = c("B.","C."),font.label=list(face="plain",size=15)),
          nrow=2)

ggsave("./Fig_Output/ch4_wave.jpg",width=10,height=8,units="in",dpi=320)

##########################
# Curious to see if we can elucidate diel cycles using the summer vs. winter time periods
# Start with summer CO2
summer_hourly <- fcr_hourly %>% 
  filter(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 06:00:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2020-06-29 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-29 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2020-06-30 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-06-30 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2020-07-01 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-01 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2020-07-02 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-02 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2020-07-03 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-03 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2020-07-04 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-04 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2020-07-05 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2020-07-05 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2020-07-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2020-07-06 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

# Format for wavelet analysis: start with gap-filled data every 30 minutes?
co2_summer <- summer_hourly %>% 
  select(DateTime,NEE) %>% 
  mutate(Time = 1:length(summer_hourly$DateTime))

headers<-names(co2_summer)
all<-headers[2]
temp<-matrix(-99,length(co2_summer$Time),length(all),byrow=F)
head(temp)

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(co2_summer)) {
  
  temp[,j] <- (co2_summer$NEE - mean(co2_summer$NEE, na.rm = TRUE)) / sd(co2_summer$NEE, na.rm = TRUE)
  
}

snt_co2_summer<-as.data.frame(temp)
snt_co2_summer<-setNames(snt_co2_summer, all)#new summerframe with standard normal transformed summer
head(snt_co2_summer)

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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

#### Getting Wavelet output and plotting ####

Time<-co2_summer$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_summer<-cbind(Time, snt_co2_summer)
head(snt_co2_summer)
snt_co2_summer <- snt_co2_summer %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_summer))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_summer)

output<-morlet(snt_co2_summer$NEE, snt_co2_summer$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
View(output$Power)
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}

View(dummy)

powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot_co2_summer <- powerplot %>% 
  filter(period < 2)

co2_power_summer <- ggplot(data = powerplot_co2, mapping = aes(x = period, y = mean_power))+
  annotate(geom="text",x=1.75,y=1.5,label="CO2 summer")+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 20,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 20,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
co2_power_summer

# Winter CO2
# Start with winter
winter_hourly <- fcr_hourly %>% 
  filter(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 06:00:00", tz = "EST")) %>% 
  mutate(day = ifelse(DateTime >= as.POSIXct("2021-01-04 07:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-04 23:30:00", tz="EST"), 1, 
                      ifelse(DateTime >= as.POSIXct("2021-01-05 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-05 23:30:00", tz = "EST"), 2,
                             ifelse(DateTime >= as.POSIXct("2021-01-06 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-06 23:30:00", tz = "EST"), 3,
                                    ifelse(DateTime >= as.POSIXct("2021-01-07 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-07 23:30:00", tz = "EST"), 4,
                                           ifelse(DateTime >= as.POSIXct("2021-01-08 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-08 23:30:00", tz = "EST"), 5,
                                                  ifelse(DateTime >= as.POSIXct("2021-01-09 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-09 23:30:00", tz = "EST"), 6,
                                                         ifelse(DateTime >= as.POSIXct("2021-01-10 00:00:00", tz= "EST") & DateTime <= as.POSIXct("2021-01-10 23:30:00", tz = "EST"), 7,
                                                                ifelse(DateTime >= as.POSIXct("2021-01-11 00:00:00", tz = "EST") & DateTime <= as.POSIXct("2021-01-11 23:30:00", tz = "EST"), 8,
                                                                       NA))))))))) %>% 
  mutate(hour = hour(DateTime))

# Format for wavelet analysis: start with gap-filled data every 30 minutes?
co2_winter <- winter_hourly %>% 
  select(DateTime,NEE) %>% 
  mutate(Time = 1:length(winter_hourly$DateTime))

headers<-names(co2_winter)
all<-headers[2]
temp<-matrix(-99,length(co2_winter$Time),length(all),byrow=F)
head(temp)

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(co2_winter)) {
  
  temp[,j] <- (co2_winter$NEE - mean(co2_winter$NEE, na.rm = TRUE)) / sd(co2_winter$NEE, na.rm = TRUE)
  
}

snt_co2_winter<-as.data.frame(temp)
snt_co2_winter<-setNames(snt_co2_winter, all)#new winterframe with standard normal transformed winter
head(snt_co2_winter)

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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

#### Getting Wavelet output and plotting ####

Time<-co2_winter$Time
head(Time)
#Time <- Time[1:21586]
snt_co2_winter<-cbind(Time, snt_co2_winter)
head(snt_co2_winter)
snt_co2_winter <- snt_co2_winter %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_co2_winter))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_co2_winter)

output<-morlet(snt_co2_winter$NEE, snt_co2_winter$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
View(output$Power)
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}

View(dummy)

powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot_co2_winter <- powerplot %>% 
  filter(period < 2)

co2_power_winter <- ggplot(data = powerplot_co2_winter, mapping = aes(x = period, y = mean_power))+
  annotate(geom="text",x=1.75,y=1.5,label="CO2 Winter")+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 20,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 20,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
co2_power_winter

# CH4 summer
# Format for wavelet analysis: start with gap-filled data every 30 minutes?
ch4_summer <- summer_hourly %>% 
  select(DateTime,CH4) %>% 
  mutate(Time = 1:length(summer_hourly$DateTime))

headers<-names(ch4_summer)
all<-headers[2]
temp<-matrix(-99,length(ch4_summer$Time),length(all),byrow=F)
head(temp)

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(ch4_summer)) {
  
  temp[,j] <- (ch4_summer$CH4 - mean(ch4_summer$CH4, na.rm = TRUE)) / sd(ch4_summer$CH4, na.rm = TRUE)
  
}

snt_ch4_summer<-as.data.frame(temp)
snt_ch4_summer<-setNames(snt_ch4_summer, all)#new summerframe with standard normal transformed summer
head(snt_ch4_summer)

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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
              add = TRUE, lwd = 2, col = "black", max.contour.segments > 25000) # still CH4d to run max.contour.segments on line 81
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

#### Getting Wavelet output and plotting ####

Time<-ch4_summer$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_summer<-cbind(Time, snt_ch4_summer)
head(snt_ch4_summer)
snt_ch4_summer <- snt_ch4_summer %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_summer))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_summer)

output<-morlet(snt_ch4_summer$CH4, snt_ch4_summer$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
View(output$Power)
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}

View(dummy)

powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot_ch4_summer <- powerplot %>% 
  filter(period < 2)

ch4_power_summer <- ggplot(data = powerplot_ch4_summer, mapping = aes(x = period, y = mean_power))+
  annotate(geom="text",x=1.75,y=1.5,label="CH4 summer")+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 20,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 20,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
ch4_power_summer

# CH4 winter
# Format for wavelet analysis: start with gap-filled data every 30 minutes?
ch4_winter <- winter_hourly %>% 
  select(DateTime,CH4) %>% 
  mutate(Time = 1:length(winter_hourly$DateTime))

headers<-names(ch4_winter)
all<-headers[2]
temp<-matrix(-99,length(ch4_winter$Time),length(all),byrow=F)
head(temp)

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(ch4_winter)) {
  
  temp[,j] <- (ch4_winter$CH4 - mean(ch4_winter$CH4, na.rm = TRUE)) / sd(ch4_winter$CH4, na.rm = TRUE)
  
}

snt_ch4_winter<-as.data.frame(temp)
snt_ch4_winter<-setNames(snt_ch4_winter, all)#new winterframe with standard normal transformed winter
head(snt_ch4_winter)

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
  x <- ((x*10)/(60*24))  ## added here; makes x axis = # of days since 2018-10-02 
  wave <- wave.list$wave
  period <- wave.list$period
  Signif <- wave.list$Signif
  coi <- wave.list$coi * (10/(60*24)) # this fixes COI to correspond w/ editing period values so they = days 
  coi<- coi              ### edit here to get COI, removed *14 
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
              add = TRUE, lwd = 2, col = "black", max.contour.segments > 25000) # still CH4d to run max.contour.segments on line 81
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

#### Getting Wavelet output and plotting ####

Time<-ch4_winter$Time
head(Time)
#Time <- Time[1:21586]
snt_ch4_winter<-cbind(Time, snt_ch4_winter)
head(snt_ch4_winter)
snt_ch4_winter <- snt_ch4_winter %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_ch4_winter))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_ch4_winter)

output<-morlet(snt_ch4_winter$CH4, snt_ch4_winter$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis 

wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 

#### Making mean global power plots and calculations ####

#creating theme for plots 
mytheme_AS <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    legend.key = element_blank(),legend.background = element_blank(),
                    legend.title = element_text(size = 12),
                    legend.text=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.title=element_text(size=12,
                                            #face="bold"
                    ),
                    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

#calculating mean power per period 
dummy <- c(1:175)
for (j in 1:ncol(output$Power)) {
  dummy[j] <- mean(output$Power[,j])
}


powerplot <- as.data.frame(cbind(dummy, output$period))
head(powerplot)
powerplot <- rename(powerplot, c( "mean_power" = "dummy"))
powerplot <- rename(powerplot, c( "period" = "V2"))
head(powerplot)

# Month mean power plot 
powerplot_ch4_winter <- powerplot %>% 
  filter(period < 2)

ch4_power_winter <- ggplot(data = powerplot_ch4_winter, mapping = aes(x = period, y = mean_power))+
  annotate(geom="text",x=1.75,y=1.5,label="CH4 Winter")+
  geom_vline(xintercept = 1,linetype="dashed")+
  annotate(geom="text",x = 1.1,y = 20,label = "1 d.")+
  geom_vline(xintercept = 0.12,linetype="dashed")+
  annotate(geom="text",x = 0.17,y = 20,label = "4 hr.")+
  geom_line()+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
ch4_power_winter

ggarrange(co2_power_summer,ch4_power_summer,co2_power_winter,ch4_power_winter)

ggsave("./Fig_Output/Wavelet_week.jpg",width=10,height=8,units="in",dpi=320)
