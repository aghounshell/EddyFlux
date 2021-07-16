###Final Wavelet Analysis 
#DWH modified script originally from CCC 
#Fall 2019 - Summer 2020 


###Reading in packages 
library(zoo)
library(dplR)#Use this instead of Rwave for the morlet analysis
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

#### Getting Data formated ####

temp<-read.csv("./Data/exo_for_wavelet.csv", header=TRUE)
head(temp)
min(temp$EXO_fdom_qsu, na.rm = TRUE)
max(temp$EXO_fdom_qsu, na.rm = TRUE)
temp_na <- temp %>% 
  filter(is.na(fdom))

head(temp)
temp <- subset(temp, select = -c(EXO_fdom_qsu, EXOTemp_C_1))
head(temp)
temp <- rename(temp, "fdom_withNA" = 5)
head(temp)  

fdom<-na.approx(temp$fdom_withNA) #replaces NA w/ interpolation

temp <-cbind(temp, fdom)
head(temp)
temp <- subset(temp, select = -c(fdom_withNA))
head(temp)


#name temp as phytos to fit earlier code
phytos<-temp
head(phytos)


headers<-names(phytos)
allspecies<-headers[5]
temp<-matrix(-99,length(phytos$Time),length(allspecies),byrow=F)
head(temp)


#old loop using genABLE package 
# for(j in 1:length(allspecies)){
#   
#   temp[,j]<-ztransform(eval(parse(text=paste0("phytos$",allspecies[j]))), phytos)#ztransform the data prior to the analysis
#   
# }

#new loop for ztransform just using manual calculation for ztransform 
for (j in 1:length(allspecies)) {
  
  temp[,j] <- (phytos$fdom - mean(phytos$fdom, na.rm = TRUE)) / sd(phytos$fdom, na.rm = TRUE)
  
}


snt_phytos<-as.data.frame(temp)
snt_phytos<-setNames(snt_phytos, allspecies)#new dataframe with standard normal transformed data
head(snt_phytos)



#other z transform test, using for AR. works the exact same! 
# test_scale <- scale(phytos$fdom, center = TRUE, scale = TRUE) # initial function from MEL that's used in AR drivers
# test_manualcalc <- (phytos$fdom - mean(phytos$fdom, na.rm = TRUE)) / sd(phytos$fdom, na.rm = TRUE) #manual z-score calc from MEL
# testall <- cbind(snt_phytos, test_scale, test_manualcalc)

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

Time<-phytos$Time
head(Time)
#Time <- Time[1:21586]
snt_phytos<-cbind(Time, snt_phytos)
head(snt_phytos)
snt_phytos <- snt_phytos %>% 
  mutate(Time = round(as.numeric(c(0.00001:nrow(snt_phytos))), digits = 0))  ##changing time to 0:xx instead 1:xx to make x axis represent days since 2 oct 2018
head(snt_phytos)


output<-morlet(snt_phytos$fdom, snt_phytos$Time, dj=(1/12), siglvl = 0.95, p2= 14.5) ###p2 is 2^ whatever value thats set. # NULL sets p2 to 15. #Chose 14.5 so that plot has minimal area above COI

output$period <- round( (output$period * (10/(60*24)) ), digits = 4)  #This works with correcting period to days for y axis  
#View(output$period)


wavelet.plot.new(output) #export in 977 x 701 for all numbers to show on color grid 



##to save wavelet plot
# for(j in 1:length(allspecies)){
#   
#   pdf(paste0("WaveletAnalysis_",allspecies[j],"_24oct19.pdf"))
#   output<-morlet(eval(parse(text=paste0("snt_phytos$",allspecies[j]))),snt_phytos$Time, dj=(1/12), siglvl = 0.95, p2=12) #p2=7.3)#ztransform the data prior to the analysis
#   wavelet.plot.new(output)		
#   dev.off()		
# }



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
  filter(period < 80)

rect1 <- data.frame (xmin = 0, xmax = 5, ymin=-Inf, ymax=Inf) #making rectangle to show range in daily plot

monthlyPower <- ggplot(data = powerplot, mapping = aes(x = period, y = mean_power))+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  geom_rect(data= rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE) +
  mytheme_AS
monthlyPower

#Day mean power plot 
powerplotday <- powerplot %>% 
  filter(period < 5)

dailyPower <- ggplot(data = powerplotday, mapping = aes(x = period, y = mean_power))+
  geom_point(color = "black", size = 1)+
  labs(x = "Period (days)",
       y = (expression(paste("Mean ",Power^2,))))+
  mytheme_AS
dailyPower

#combining two plots and saving 
fig5meanpower <- ggarrange(monthlyPower, dailyPower,labels = c("b", "c"), ncol = 2, nrow = 1)

fig5meanpower

ggsave(filename = "./Images/Wavelet_meanpowerplots_dailyhighlighted.tiff", 
       fig5meanpower, device = "tiff", width = 230, height = 100, units = "mm")



#### Logistic regressions ####
#Used to determine how many days in our study period exhibited significant daily time scale of variability and look at what environemntal driver conditions were on those days 

#view output and variables working with 
View(output)
View(output$x)
View(output$Signif)
View(output$period)
View(output$wave)


#taking calcs from wavelet.plot.new function
View(output$Signif)
signif_test <- t(matrix(output$Signif, dim(output$wave)[2], dim(output$wave)[1]))
signif_testa <- output$Power/signif_test            #puts signif in 52549 x 16 matrix 
View(output$Power)

a <- signif_testa 
b <- output$Power

#c <- array( c( a, b) , dim = c(52549, 16, 2 ) )

#hold <- c(1:52549)

#hold <- data.frame(matrix(ncol = 16, nrow = 52549))

### mean of signifcant value per row, DOESNT WORK
# outputa <- vector("double", nrow(signif_testa))
# for(j in 1:nrow(signif_testa)){
#   outputa[j] <- mean(signif_testa[j,])
# }
# outputa
# logoutput <- log2(outputa)
# View(logoutput)


### YES!! This works to assign 1 or 0 to true or false values at one scale 
c <- as.data.frame(b >= signif_test)  #this tells us if the power value > signif value 

#c7 gets T/F values for if signif at daily timescale 
c7 <- c %>%  
  select(V84)
c7 <- as.vector(c7)

#sig7 pulls out signif values at daily timescale; not using now 
# sig7 <- as.data.frame(signif_test)
# sig7 <- sig7 %>% 
#   select(V84)


#if else gives us numerics for logistic regression. 1 = T, 0 = F
coutput <- ifelse(c7 == TRUE, 1, 0)

View(coutput)

#binding logsitic info to rest of data 
bind <- cbind(phytos, coutput) 
bind <- rename(bind, "logistic" = 6)
head(bind)
class(bind$logistic)

#summarizing logistics to daily scale 
logis <- bind %>% 
  mutate(Date = as.Date(DateTime)) %>% 
  group_by(Date) %>% 
  summarise(daily_fdom = mean(fdom),
            daily_logisticvalue = mean(logistic))

head(logis)

#rounding logistic values to 0 or 1
logis$roundlogistic_05 <- round(logis$daily_logisticvalue/1)*1  # rounding to 0/1 based on 0.5
logis$roundlogistic_anysig <- ifelse(logis$daily_logisticvalue > 0, 1, 0)

View(logis) #This gaves us 1/0 for each day for if the daily time scale was significant 
summary(logis)
sig_05 <- logis %>% 
  filter(roundlogistic_05 == 1) ##41 significant days at 0.5 cut off 
sig_any <- logis %>% 
  filter(roundlogistic_anysig == 1) ##62 significant days at > 0 cut off

sig_any_summer <- sig_any %>% 
  filter(Date > "2019-05-31",
         Date < "2019-09-01") # 50 of 62 days were during the summer 


#bringing in drivers to do logistic regression
daily_drivers <- read_csv("./Data/ar_daily_data_joined.csv")
head(daily_drivers)

drivers_filt <- daily_drivers %>% 
  select(Date, daily_EXO_wtr_temp, daily_SRup_mean, daily_EXO_chla_ugL, daily_EXO_do_persat, daily_rain_mm, 
         WVWA_Flow_cms_daily_mean, WRT_days_daily)
drivers_filt<- slice(drivers_filt, -1) #remvoing so dates will line up 
head(drivers_filt)
tail(drivers_filt)


logistic_filt <- logis %>% 
  select(Date, roundlogistic_05, roundlogistic_anysig, daily_fdom)
head(logistic_filt)
tail(logistic_filt)

binded_logistic <- left_join(logistic_filt, drivers_filt, by = "Date")
head(binded_logistic)
tail(binded_logistic)

### Models and plots w/ 0.5 as logistic cut off ###
{
#temp
tempLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_wtr_temp, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("Temp")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
tempLR
model_tempLR <- glm(roundlogistic_05 ~ daily_EXO_wtr_temp, family = "binomial" , data = binded_logistic)
summary(model_tempLR)  # p < 0.001 ***
coef(summary(model_tempLR))  # p < 0.001 ***

#SW
swLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_SRup_mean, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("Shortwave")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
swLR
model_swLR <- glm(roundlogistic_05 ~ daily_SRup_mean, family = "binomial" , data = binded_logistic)
summary(model_swLR) # p < 0.001 ***
tidy(model_swLR)

#chla
chlaLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_chla_ugL, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("Chl-a")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
chlaLR
model_chlaLR <- glm(roundlogistic_05 ~ daily_EXO_chla_ugL, family = "binomial" , data = binded_logistic)
summary(model_chlaLR) # p < 0.001
tidy(model_chlaLR)

#DO %
doLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_do_persat, y = roundlogistic_05))+
  geom_point(alpha = .15)+
    ggtitle("DO % sat")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
doLR
model_doLR <- glm(roundlogistic_05 ~ daily_EXO_do_persat, family = "binomial" , data = binded_logistic)
summary(model_doLR)

#rain
rainLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_rain_mm, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("Precipitation")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
rainLR
model_rainLR <- glm(roundlogistic_05 ~ daily_rain_mm, family = "binomial" , data = binded_logistic)
summary(model_rainLR)

#flow
flowLR <- ggplot(data = binded_logistic, mapping = aes(x = WVWA_Flow_cms_daily_mean, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("Inflow")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
flowLR
model_flowLR <- glm(roundlogistic_05 ~ WVWA_Flow_cms_daily_mean, family = "binomial" , data = binded_logistic)
summary(model_flowLR)

#wrt 
wrtLR <- ggplot(data = binded_logistic, mapping = aes(x = WRT_days_daily, y = roundlogistic_05))+
  geom_point(alpha = .15)+
  ggtitle("WRT")+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) 
wrtLR
model_wrtLR <- glm(roundlogistic_05 ~ WRT_days_daily, family = "binomial" , data = binded_logistic)
summary(model_wrtLR)

model_allLR <- glm(roundlogistic_05 ~ daily_EXO_wtr_temp + daily_SRup_mean + daily_EXO_chla_ugL + 
                     daily_EXO_do_persat + daily_rain_mm   + WRT_days_daily, 
                   family = "binomial" , data = binded_logistic)
summary(model_allLR)


driversLR <- ggarrange(tempLR, swLR, chlaLR, doLR, rainLR, flowLR, wrtLR,
                    #labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                    ncol = 4, nrow = 2)
driversLR

#ggsave(filename = "./Final_Data_Scripts/Images/LogisticRegression_drivers_roundedlogic_05.tiff", 
#        driversLR, device = "tiff", width = 230, height = 140, units = "mm")

}


### Models and plots w/ logistic > 0 as signif ###
{
  #temp
  tempLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_wtr_temp, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("Temp")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  tempLR
  model_tempLR <- glm(roundlogistic_anysig ~ daily_EXO_wtr_temp, family = "binomial" , data = binded_logistic)
  summary(model_tempLR)  # p < 0.001 ***

  #SW
  swLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_SRup_mean, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("Shortwave")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  swLR
  model_swLR <- glm(roundlogistic_anysig ~ daily_SRup_mean, family = "binomial" , data = binded_logistic)
  summary(model_swLR) # p < 0.001 ***

  #chla
  chlaLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_chla_ugL, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("Chl-a")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  chlaLR
  model_chlaLR <- glm(roundlogistic_anysig ~ daily_EXO_chla_ugL, family = "binomial" , data = binded_logistic)
  summary(model_chlaLR) # p < 0.001

  #DO %
  doLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_EXO_do_persat, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("DO % sat")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  doLR
  model_doLR <- glm(roundlogistic_anysig ~ daily_EXO_do_persat, family = "binomial" , data = binded_logistic)
  summary(model_doLR) # p < 0.001
  
  #rain
  rainLR <- ggplot(data = binded_logistic, mapping = aes(x = daily_rain_mm, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("Precipitation")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  rainLR
  model_rainLR <- glm(roundlogistic_anysig ~ daily_rain_mm, family = "binomial" , data = binded_logistic)
  summary(model_rainLR) # not sig, p = 0.9
  
  #flow
  flowLR <- ggplot(data = binded_logistic, mapping = aes(x = WVWA_Flow_cms_daily_mean, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("Inflow")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  flowLR
  model_flowLR <- glm(roundlogistic_anysig ~ WVWA_Flow_cms_daily_mean, family = "binomial" , data = binded_logistic)
  summary(model_flowLR) # p < 0.001 for flow but not intercept 
  
  #wrt 
  wrtLR <- ggplot(data = binded_logistic, mapping = aes(x = WRT_days_daily, y = roundlogistic_anysig))+
    geom_point(alpha = .15)+
    ggtitle("WRT")+
    geom_smooth(method = "glm", method.args = list(family = "binomial")) 
  wrtLR
  model_wrtLR <- glm(roundlogistic_anysig ~ WRT_days_daily, family = "binomial" , data = binded_logistic)
  summary(model_wrtLR) # p < 0.001
  
  model_allLR <- glm(roundlogistic_anysig ~ daily_EXO_wtr_temp + daily_SRup_mean + daily_EXO_chla_ugL + 
                       daily_EXO_do_persat + daily_rain_mm   + WRT_days_daily, 
                     family = "binomial" , data = binded_logistic)
  summary(model_allLR) # intercept and temp: p < 0.001, SW, Chla, WRT: p < 0.05
  
  
  driversLR <- ggarrange(tempLR, swLR, chlaLR, doLR, rainLR, flowLR, wrtLR,
                         #labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                         ncol = 4, nrow = 2)
  driversLR
  
#ggsave(filename = "./Final_Data_Scripts/Images/LogisticRegression_drivers_roundedlogic_anysig.tiff", 
#          driversLR, device = "tiff", width = 230, height = 140, units = "mm")
  
}



#### Breakpoints for LR: Looking at breakpoints for when/how daily scales emerge ####

count(binded_logistic, vars = roundlogistic_anysig) # 62 significant days for all sig 
count(binded_logistic, vars = roundlogistic_05) # 41 significant days for sig > 0.5

#temp > 25 C
temp25 <- binded_logistic %>% 
  select(Date, roundlogistic_anysig, daily_EXO_wtr_temp) %>% 
  filter(daily_EXO_wtr_temp > 25)
  #filter(daily_EXO_wtr_temp > 25 & roundlogistic_anysig == 1)
head(temp25)
count(temp25, vars = roundlogistic_anysig) # 43 of 60 days above 25C were significant 
round(((43/62)*100), digits = 0) # 43 of 62 significant days occured above 25C, 43/62 = 69% 

temp25a <- binded_logistic %>% 
  select(Date, roundlogistic_05, daily_EXO_wtr_temp) %>% 
  filter(daily_EXO_wtr_temp > 25)
#filter(daily_EXO_wtr_temp > 25 & roundlogistic_anysig == 1)
head(temp25a)
count(temp25a, vars = roundlogistic_05) # 31 of 60 days above 25C were significant 
round(((31/41)*100), digits = 0) # 31 of 41 significant days occured above 25C, 31/41 = 76%


#do > 100
do100 <- binded_logistic %>% 
  select(Date, roundlogistic_anysig, daily_EXO_do_persat) %>% 
  filter(daily_EXO_do_persat > 100)
head(do100)
count(do100, vars = roundlogistic_anysig) # 62 of 278 days were significant 
round(((62/62)*100), digits = 0) # 62 of 62 significant days , 62/62 = 100% 

do100a <- binded_logistic %>% 
  select(Date, roundlogistic_05, daily_EXO_do_persat) %>% 
  filter(daily_EXO_do_persat > 100)
head(do100a)
count(do100a, vars = roundlogistic_05) # 41 of 278 days were significant 
round(((41/41)*100), digits = 0) # 41 of 41 significant days, 41/41 = 100%


#Sw > 200
mean(binded_logistic$daily_SRup_mean)
sw200 <- binded_logistic %>% 
  select(Date, roundlogistic_anysig, daily_SRup_mean) %>% 
  filter(daily_SRup_mean > 200)
head(sw200)
count(sw200, vars = roundlogistic_anysig) # 49 of 149 days were significant 
round(((49/62)*100), digits = 0) # 49 of 62 significant days, 49/62 = 79% 
#mean SW = 164, 55 of 62 days = 89% occured above mean SW 

sw200a <- binded_logistic %>% 
  select(Date, roundlogistic_05, daily_SRup_mean) %>% 
  filter(daily_SRup_mean > 200)
head(sw200a)
count(sw200a, vars = roundlogistic_05) # 32 of 149  were significant 
round(((32/41)*100), digits = 0) # 41 of 41 significant days, 32/41 = 79%
#mean SW = 164, 37 of 41 days = 90% occured above mean SW 







