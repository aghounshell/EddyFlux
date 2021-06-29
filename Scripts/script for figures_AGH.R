library(readr)
library(dplyr)
library(lubridate)
library(zoo)
library(ggplot2)
library(ggpubr)
library(openair)
library(REddyProc)

pacman::p_load(OneR)

# Set working directory
wd <- getwd()
setwd(wd)

# read data 

fcr_gf <- read_csv("./Data/20210615_EC_processed.csv")

# set DateTime to EST

fcr_gf$DateTime <- as_datetime(fcr_gf$DateTime, tz = 'EST')

# if solar radiation more than 12, it's day, otherwise, it's night

fcr_gf$Day_time <- ifelse(fcr_gf$Rg_f >= 12, 1, 0)

#############################################################
#COUNTING NA

nas_percent <- fcr_gf %>% select(DateTime, NEE, ch4_flux, 
																 NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
						ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)


nas_count <- fcr_gf %>% select(DateTime, NEE, ch4_flux,
															 NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	summarise(co2_available = n() - sum(is.na(NEE_uStar_orig)),
						ch4_available = n()-sum(is.na(ch4_flux_uStar_orig)))

fcr_gf %>% select(DateTime, NEE, ch4_flux, 
									NEE_uStar_orig, ch4_flux_uStar_orig) %>% 
	group_by(year(DateTime), month(DateTime)) %>% 
	summarise(co2_available = 100- sum(is.na(NEE_uStar_orig))/n()*100,
						ch4_available = 100-sum(is.na(ch4_flux_uStar_orig))/n()*100)


#############################################################

#############################################################
# Check where fluxes are coming from: using NON gap-filled data
ggplot(fcr_gf,mapping=aes(x=wind_dir,y=NEE_uStar_orig,color=Ustar))+
  geom_point()

ggplot(fcr_gf,mapping=aes(x=wind_dir,y=ch4_flux_uStar_orig,color=Ustar))+
  geom_point()

fcr_gf %>% ggplot(aes(wind_dir, NEE_uStar_orig,color=Ustar)) + 
  geom_hline(yintercept = 40, col = 'goldenrod2', lwd = 2) +
  geom_hline(yintercept = 50, col = 'green', lwd = 1.4) +
  geom_hline(yintercept = 100, col = 'blue', lwd = 1.4) +
  geom_hline(yintercept = 120, col = 'gray2', lwd = 1.4) +
  geom_hline(yintercept = 150, col = 'red',lwd = 1.4) +
  geom_point() +
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 360, 45)) +
  theme_bw() + 
  coord_polar()

# plotting 30min values
#air temperature

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, Tair_f)) + 
	geom_line(aes(DateTime, Tair_f)) +
	theme_bw()

# wind speed

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, u)) + 
	geom_line(aes(DateTime, u)) +
	theme_bw()

# NEE

fcr_gf %>% ggplot() + 
	geom_point(aes(DateTime, NEE_uStar_f), color = 'red') +
	geom_line(aes(DateTime, NEE_uStar_orig), color = 'gray1') + 
	theme_bw()


# CH4

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, ch4_flux_uStar_f), col = 'red') + 
	geom_line(aes(DateTime, ch4_flux_uStar_orig), color = 'gray1') +
	theme_bw()


# air pressure

fcr_gf %>% ggplot() +
	geom_point(aes(DateTime, airP)) + 
	geom_line(aes(DateTime, airP)) +
	theme_bw()


#####

# Group by hour, plot hourly data
fcr_hourly <- fcr_gf %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST")) %>% 
  group_by(DateTime) %>% 
	summarise(NEE = mean(NEE_uStar_f, na.rm = TRUE),
						NEE05 = mean(NEE_U05_f, na.rm = TRUE),
						NEE50 = mean(NEE_U50_f, na.rm = TRUE),
						NEE95 = mean(NEE_U95_f, na.rm = TRUE),
						NEE_sd = sqrt(sum(NEE_uStar_fsd^2, na.rm = TRUE)),
						CH4 = mean(ch4_flux_uStar_f, na.rm = TRUE),
						CH405 = mean(ch4_flux_U05_f, na.rm = TRUE),
						CH450 = mean(ch4_flux_U50_f, na.rm = TRUE),
						CH495 = mean(ch4_flux_U95_f, na.rm = TRUE),
						CH4_sd = sqrt(sum(ch4_flux_uStar_fsd^2, na.rm = TRUE)),
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

# Plot hourly CO2 data over the year
ggplot(fcr_hourly,mapping=aes(x=DateTime))

# plot co2 hourly data
fcr_gf  %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																	"Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																	'Dec', 'Jan', 'Feb', 'Mar')), 
										 hour = hour(DateTime)) %>% 
	summarise(NEE_m = mean(NEE_uStar_f, na.rm = TRUE),
						NEEsd = sd(NEE_uStar_f, na.rm = TRUE)) %>% 
	ggplot(aes(hour, NEE_m, col = factor(year))) + 
	geom_point(size = 3) + geom_line(size = 1) +
	facet_wrap(~month) + ylab(expression(CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(" ") +
	geom_errorbar(aes(ymin = NEE_m - NEEsd, ymax = NEE_m + NEEsd), 
								width = 0.5, size = 0.8) +
	geom_hline(yintercept = 0, lty = 2, col = 'red') +
	theme_bw() + labs(col='Year') + scale_color_brewer(palette = "Dark2")

# plot ch4 hourly data

fcr_gf %>% 
	dplyr::group_by(year = year(DateTime), 
									month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																				"Jul", "Aug", "Sep", "Oct", 'Nov', 
																																				'Dec', 'Jan', 'Feb', 'Mar')), 
									hour = hour(DateTime)) %>% 
	dplyr::summarise(ch4_m = mean(ch4_flux_uStar_f, na.rm = TRUE), 
									 ch4_sd = sd(ch4_flux_uStar_f)) %>% ungroup() %>% 
	ggplot(aes(hour, ch4_m, col = factor(year))) + 
	geom_point(size = 3) + geom_line(size = 1) +
	facet_wrap(~month) + ylab(expression(CH[4]~(mu~mol~m^-2~s^-1))) +
	xlab(" ") +
	geom_errorbar(aes(ymin = ch4_m - ch4_sd, ymax = ch4_m + ch4_sd), 
								width = 0.5, size = 0.8) +
	geom_hline(yintercept = 0, lty = 2, col = 'red') +
	theme_bw() + labs(col='Year') + scale_color_brewer(palette = "Dark2")

# plot radiation hourly

fcr_gf %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																 "Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																 'Dec', 'Jan', 'Feb', 'Mar')), 
										hour = hour(DateTime)) %>% 
	summarise(Rg = mean(Rg, na.rm = TRUE)) %>% 
	ggplot(aes(hour, Rg, col = factor(year))) + geom_point() + 
	facet_wrap(~month) + theme_bw() + ylab('Incoming radiation') + xlab("") +
	scale_color_brewer(palette = "Dark2")


# air temperature 

fcr_gf %>% group_by(year = year(DateTime), month = factor(month.abb[month(DateTime)], levels = c("Apr", "May", "Jun",
																																																 "Jul", "Aug", "Sep", "Oct", 'Nov', 
																																																 'Dec', 'Jan', 'Feb', 'Mar')), 
										hour = hour(DateTime)) %>% 
	summarise(airt = mean(Tair_f, na.rm = TRUE)) %>% 
	ggplot(aes(hour, airt, col = factor(year))) + geom_point() + 
	facet_wrap(~month) + theme_bw() + ylab('Air temperature') + xlab("") +
	scale_color_brewer(palette = "Dark2")

# group by hour 

fcr_hourly$datetime <- with(fcr_hourly, ymd_h(paste(Year, Month, Day, Hour, sep= ' ')))

# plot data averaged hourly

fcr_hourly %>%
	ggplot() +
	geom_point(aes(datetime, CH4)) + 
	geom_line(aes(datetime, CH4)) +
	ylab(expression(CH[4]~(mu~mol~m^-2~s^-1))) +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), lty =2) +
	scale_x_datetime(date_breaks = '2 week', date_labels = '%b-%d') +
	xlab(" ") +
	theme_bw()


fcr_hourly %>%
	ggplot() +
	geom_point(aes(datetime, NEE)) + 
	geom_line(aes(datetime, NEE)) +
	ylab(expression(CO[2]~(mu~mol~m^-2~s^-1))) +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), lty =2) +
	scale_x_datetime(date_breaks = '2 weeks', date_labels = '%b-%d') +
	xlab(" ") +
	theme_bw()

#####################################################################################

# Group by day -> plot

# data in umolm2s
fcr_daily <- fcr_hourly %>% dplyr::group_by(Year, Month, Day) %>% 
	dplyr::summarise(NEE_mean = mean(NEE, na.rm = TRUE),
									 NEE_U5 = mean(NEE05, na.rm = TRUE),
									 NEE_U50 = mean(NEE50, na.rm = TRUE),
									 NEE_U95 = mean(NEE95, na.rm = TRUE),
									 NEE_sd = sd(NEE, na.rm = TRUE),
									 CH4_mean = mean(CH4, na.rm = TRUE),
									 CH4_U5 = mean(CH405, na.rm = TRUE),
									 CH4_U50 = mean(CH450, na.rm = TRUE),
									 CH4_U95 = mean(CH495, na.rm = TRUE),
									 CH4_sd = sd(CH4, na.rm = TRUE),
									 umean = mean(umean, na.rm = TRUE),
									 umax = max(umax, na.rm = TRUE),
									 umin= min(umin, na.rm = TRUE),
									 Tmin = min(Tmin, na.rm = TRUE),
									 Tmax = max(Tmax, na.rm = TRUE),
									 Tmean = mean(Tmean, na.rm = TRUE),
									 pressure = mean(pressure, na.rm = TRUE),
									 minpress = min(minpress, na.rm = TRUE),
									 maxpress = max(maxpress, na.rm = TRUE),
									 H = mean(H, na.rm = TRUE),
									 LE = mean(LE, na.rm = TRUE),
									 VPD = mean(VPD, na.rm = TRUE),
									 RH = mean(RH, na.rm = TRUE),
									 PAR = mean(PAR_tot, na.rm = TRUE),
									 PPT_tot = sum(precip_sum, na.rm =TRUE),
									 LW_in= mean(LW_in, na.rm = TRUE),
									 LW_out = mean(LW_out, na.rm = TRUE),
									 Rn = mean(Rn, na.rm = TRUE),
									 Rg = mean(Rg, na.rm = TRUE),
									 SW_out = mean(SW_out, na.rm = TRUE))

fcr_daily$Date <- as.POSIXct(paste(fcr_daily$Year, fcr_daily$Month, fcr_daily$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')


# nee value in gCm2day, CH4 values = mgCm2day

fcr_dailysum <- fcr_hourly %>% dplyr::group_by(Year, Month, Day) %>% 
	dplyr::summarise(
		NEE_sum = sum(NEE*3600*12/1000000, na.rm = TRUE),
		NEE_sum05 = sum(NEE05*3600*12/1000000, na.rm = TRUE),
		NEE_sum50 = sum(NEE50*3600*12/1000000, na.rm = TRUE),
		NEE_sum95 = sum(NEE95*3600*12/1000000, na.rm = TRUE),
		CH4_sum = sum(CH4*3600*12/1000, na.rm = TRUE),
		CH4_sum05 = sum(CH405*3600*12/1000, na.rm = TRUE),
		CH4_sum50 = sum(CH450*3600*12/1000, na.rm = TRUE),
		CH4_sum95 = sum(CH495*3600*12/1000, na.rm = TRUE),
		Tmin = min(Tmin, na.rm = TRUE),
		Tmax = max(Tmax, na.rm = TRUE),
		Tmean = mean(Tmean, na.rm = TRUE),
		H = sum(H, na.rm = TRUE),
		LE = sum(LE, na.rm = TRUE),
		pressure = mean(pressure, na.rm = TRUE),
		VPD = mean(VPD, na.rm = TRUE),
		RH = mean(RH, na.rm = TRUE),
		umean = mean(umean, na.rm = TRUE),
		umax = max(umax, na.rm = TRUE),
		PAR = sum(PAR_tot, na.rm = TRUE),
		PPT_tot = sum(precip_sum, na.rm = TRUE),
		LW_in = mean(LW_in, na.rm = TRUE),
		LW_out = mean(LW_out, na.rm = TRUE),
		Rg = mean(Rg, na.rm = TRUE),
		SW_out = mean(SW_out, na.rm = TRUE),
		Rn = mean(Rn , na.rm = TRUE))


fcr_dailysum$Date <- as.POSIXct(paste(fcr_dailysum$Year, fcr_dailysum$Month, fcr_dailysum$Day, sep = '-'), "%Y-%m-%d", tz = 'EST')

# daily means
dailynee <- fcr_daily %>% 
	ggplot() +
	geom_line(aes(Date, NEE_mean)) +
	xlab("") +
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	# annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#          alpha = .2, fill = 'gray70') +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '7 day', date_labels = '%d %b') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 12, colour = 'black'))


dailych4 <- fcr_daily %>% 
	ggplot() +
	geom_line(aes(Date, CH4_mean*1000)) +
	xlab("") +
	ylab(expression(~CH[4]~(nmol~m^-2~s^-1)))+
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '7 day', date_labels = '%d %b') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 12, colour = 'black'))


dailynee
dailych4

ggarrange(dailynee, dailych4, nrow = 2, ncol = 1, align = "v")

########################################3#
# plot meteorological variables


tfig <- fcr_daily %>% ggplot(aes(Date, Tmean)) + 
	geom_line(col = 'red', size = 1) +
	geom_ribbon(aes(ymin = Tmin, ymax = Tmax), fill = 'darkred', alpha = 0.2) +
	ylab(expression(T[a]~(degree~C))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw()  +theme(axis.text = element_text(size = 12, colour = 'black'),
										 axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
										 axis.title = element_text(size = 16, colour = 'black'))


pptfig <- fcr_daily %>% ggplot(aes(Date, PPT_tot)) + 
	geom_bar(stat= 'identity', fill = 'darkblue') +
	ylab(expression(Ppt~(mm~day^-1))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

parfig <- fcr_hourly %>% ggplot(aes(datetime, PAR_avg)) + 
	geom_line(col = 'goldenrod3') +
	ylab(expression(PAR~(mu~mol~m^-2~s^-1))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

rnfig <- fcr_hourly %>% ggplot(aes(datetime, Rn)) + 
	geom_line(col = 'goldenrod3') +
	ylab(expression(Rn~(W~m^-2))) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


pressfig <- fcr_daily %>% ggplot(aes(Date, pressure)) + 
	geom_line(col = 'black', size = 1) +
	ylab(expression(P~(KPa))) +
	geom_ribbon(aes(ymin = minpress, ymax = maxpress), fill = 'black', alpha = 0.2) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ufig <- fcr_daily %>% ggplot(aes(Date, umean)) + 
	geom_line(col = 'black', size = 1) +
	ylab(expression(u~(m~s^-1))) +
	geom_ribbon(aes(ymin = umin, ymax = umax), fill = 'darkblue', alpha = 0.2) +
	xlab("") +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
	theme_bw() + 
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


ggarrange(tfig, pptfig, 
					nrow = 2, align = "v")

ggarrange(ufig, pressfig, 
					nrow = 2, align = "v")



# daily totals 

fcr_dailysum %>% 
	ggplot(aes(Date, NEE_sum)) +
	geom_point() + geom_line() +
	#  geom_vline(xintercept = as.POSIXct('2020-06-29 -05'), lty = 2, col = 'red') +
	theme_bw() + ylab(expression(CO[2]~flux~(g~C~m^-2~d^-1))) +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-11-20 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1, lty =2) + 
	geom_hline(yintercept = 0, lty =2) +
	scale_x_datetime(date_breaks = '14 days', date_labels = '%b/%d') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


fcr_dailysum %>% 
	ggplot(aes(Date, CH4_sum)) +
	geom_point() + geom_line() +
	ylab(expression(CH[4]~(mg~C~m^-2~d^-1))) +
	#  geom_vline(xintercept = as.POSIXct('2020-06-29 -05'), lty = 2, col = 'red') +
	geom_hline(yintercept = 0, lty = 2) +
	theme_bw() +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#          alpha = .2, fill = 'gray70') +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-11-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1, lty =2) + 
	scale_x_datetime(date_breaks = '14 days', date_labels = '%b/%d') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))


###################################################################
# CUMULATIVE VALUES
###################################################################


fcr_dailysum <- fcr_dailysum %>% ungroup() %>% 
	mutate(cumnee = cumsum(NEE_sum), 
				 cumch4 = cumsum(CH4_sum))


#

fcr_daily$cumnee <- fcr_dailysum$cumnee
fcr_daily$cumch4 <- fcr_dailysum$cumch4

coeff = 40

dailynee <- fcr_daily %>% 
	ggplot() +
	geom_line(aes(Date, NEE_mean*3600*24*12/1000000)) +
	geom_line(aes(Date, cumnee/coeff), col = 'darkgreen', size =1) +
	xlab("") +
	# annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#         alpha = .2, fill = 'gray70') +
	#annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#        alpha = .2, fill = 'gray70') +
	#geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%d %b') +
	scale_y_continuous(name = expression(CO[2]~(g~C~m^-2~d^-1)), 
										 sec.axis = sec_axis(trans= ~.*coeff, 
										 										name = expression(Cumulative~flux~(g~C~m^-2~d^-1)))) +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 13, colour = 'black'))

coeff2 <- 30

dailych4 <- fcr_daily %>% 
	ggplot() +
	geom_line(aes(Date, CH4_mean*3600*24*12/1000)) +
	geom_line(aes(Date, cumch4/coeff2), col = 'goldenrod3', size = 1) +
	xlab("") +
	#  annotate("rect", xmin = as.POSIXct(as.Date('2020-06-29 -5')), xmax = as.POSIXct(as.Date('2020-09-11 -5')), ymin = -Inf, ymax = Inf,
	#           alpha = .2, fill = 'gray70') +
	# annotate("rect", xmin = as.POSIXct(as.Date('2020-09-25 -5')), xmax = as.POSIXct(as.Date('2020-12-14 -5')), ymin = -Inf, ymax = Inf,
	#         alpha = .2, fill = 'gray70') +
	scale_y_continuous(name = expression(CH[4]~(mg~C~m^-2~d^-1)), 
										 sec.axis = sec_axis(trans= ~.*coeff2, 
										 										name = expression(Cumulative~flux~(mg~C~m^-2~d^-1)))) +
	#  geom_vline(xintercept = as.POSIXct('2020-11-01 18:40:00 -5'), col = 'black', size = 1) + 
	geom_hline(yintercept = 0, lty = 2) +
	scale_x_datetime(date_breaks = '1 month', date_labels = '%d %b') +
	theme_bw() +
	theme(axis.text = element_text(size = 12, colour = 'black'),
				axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 13, colour = 'black'))


dailynee
dailych4

ggarrange(dailynee, dailych4, nrow = 2, ncol = 1, align = 'v')

########################################################

# montly mean

fcr_monthly <- fcr_daily %>% group_by(Year, Month) %>% 
	summarise(NEE = mean(NEE_mean, na.rm = TRUE),
						NEE05 = mean(NEE_U5, na.rm = TRUE),
						NEE50 = mean(NEE_U50, na.rm = TRUE),
						NEE95 = mean(NEE_U95, na.rm = TRUE),
						NEE_sd = sd(NEE_mean, na.rm = TRUE),
						CH4 = mean(CH4_mean, na.rm = TRUE),
						CH405 = mean(CH4_U5, na.rm = TRUE),
						CH450 = mean(CH4_U50, na.rm = TRUE),
						CH495 = mean(CH4_U95, na.rm = TRUE),
						CH4_sd = sd(CH4_mean, na.rm = TRUE),
						Tair = mean(Tmean, na.rm = TRUE),
						Tmin = min(Tmin, na.rm =TRUE),
						Tmax = max(Tmax, na.rm =TRUE),
						pressure = mean(pressure, na.rm =TRUE),
						maxpress = max(maxpress, na.rm = TRUE),
						minpress = min(minpress, na.rm = TRUE),
						H = mean(H, na.rm = TRUE),
						LE = mean(LE, na.rm = TRUE),
						VPD = mean(VPD, na.rm = TRUE),
						RH = mean(RH, na.rm = TRUE),
						u = mean(umean, na.rm = TRUE),
						umax = max(umax, na.rm = TRUE),
						umin = min(umin, na.rm = TRUE),
						PAR = mean(PAR, na.rm = TRUE),
						PPT = sum(PPT_tot, na.rm = TRUE),
						LW_in = mean(LW_in, na.rm = TRUE),
						LW_out = mean(LW_out, na.rm = TRUE),
						Rg = mean(Rg, na.rm = TRUE),
						SW_out = mean(SW_out, na.rm = TRUE),
						Rn = mean(Rn, na.rm = TRUE))



# monthly sum NEE in gCm2month, CH4 in mgCm2month

fcr_monthlysum <- fcr_dailysum %>% 
	dplyr::group_by(Year, Month) %>% 
	dplyr::summarise(NEE = sum(NEE_sum, na.rm = TRUE),
									 NEEU5 = sum(NEE_sum05, na.rm = TRUE),
									 NEEU50 = sum(NEE_sum50, na.rm = TRUE),
									 NEEU95 = sum(NEE_sum95, na.rm = TRUE),
									 NEE_sd = sd(NEE_sum, na.rm = TRUE),
									 CH4 = sum(CH4_sum, na.rm = TRUE),
									 CH4U5 = sum(CH4_sum05, na.rm = TRUE),
									 CH4U50 = sum(CH4_sum50, na.rm = TRUE),
									 CH4U95 = sum(CH4_sum95, na.rm = TRUE),
									 CH4_sd = sd(CH4_sum, na.rm = TRUE),
									 Tair = mean(Tmean, na.rm = TRUE),
									 Tmax = max(Tmax, na.rm = TRUE),
									 Tmin = min(Tmin, na.rm = TRUE),
									 H = sum(H, na.rm = TRUE),
									 LE = sum(LE, na.rm = TRUE),
									 Rn = sum(Rn, na.rm = TRUE),
									 VPD = mean(VPD, na.rm = TRUE),
									 RH = mean(RH, na.rm = TRUE),
									 u = mean(umean, na.rm = TRUE),
									 umax= max(umax, na.rm = TRUE),
									 PAR = sum(PAR, na.rm = TRUE),
									 PPT = sum(PPT_tot, na.rm = TRUE),
									 LW_in = mean(LW_in, na.rm = TRUE),
									 LW_out = mean(LW_out, na.rm = TRUE),
									 Rg = mean(Rg, na.rm = TRUE),
									 SW_out = mean(SW_out, na.rm = TRUE))


# monthly bar plot with totals NEE and CH4

fcr_monthlysum$yearmon <- with(fcr_monthlysum, sprintf("%d-%02d", Year, Month))

neemon <- fcr_monthlysum %>% ggplot() +
	geom_bar(aes(yearmon, NEE), stat = 'identity', fill = 'darkgreen', alpha = 0.8, width = 0.3) +
	geom_errorbar(aes(yearmon, ymin = NEE - NEE_sd, ymax = NEE + NEE_sd), width = 0.2) +
	ylab(expression(~CO[2]~(g~C~m^-2~month^-1))) +
	xlab("") +
	theme_bw() +
	theme(axis.text = element_text(size = 14, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ch4mon <- fcr_monthlysum %>% ggplot() +
	geom_bar(aes(yearmon, CH4), stat = 'identity', fill = 'goldenrod3', alpha = 0.6, width = 0.3) +
	geom_errorbar(aes(yearmon, ymin = CH4 - CH4_sd, ymax = CH4 + CH4_sd), width = 0.2) +
	xlab('') + ylab(expression(~CH[4]~(mg~C~m^-2~month^-1))) +
	theme_bw() +
	theme(axis.text = element_text(size = 14, colour = 'black'),
				axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
				axis.title = element_text(size = 16, colour = 'black'))

ggarrange(neemon, ch4mon, align = 'v')


# get some numbers:

# hourly means

# CO2 in mmolm2h
mean(fcr_hourly$NEE*3600/1000)

# in mgCO2m2h
mean(fcr_hourly$NEE*3600*44/1000)

# CH4 in umolm2h
mean(fcr_hourly$CH4*3600)

# CH4 in ugCH4m2h
mean(fcr_hourly$CH4*3600*16)


# daily means

# mean of the sum - CO2 in gCm2d
mean(fcr_dailysum$NEE_sum)

# mean of the mean
mean(fcr_daily$NEE_mean*3600*24*12/1000000)

# mean of the sum - CH4 in mgCm2d
mean(fcr_dailysum$CH4_sum)

# mean of the mean mgCm2d
mean(fcr_daily$CH4_mean*3600*24*12/1000)

# monthly means
mean(fcr_monthly$NEE*3600*24*12*30/1000000) #gCm2month
mean(fcr_monthly$CH4*3600*24*12*30/1000) #mgCm2month

mean(fcr_monthlysum$NEE) #gCm2month
mean(fcr_monthlysum$CH4) #mgCm2month


# total values (more than 1 year of data)

fcr_monthlysum %>% select(yearmon, NEE, NEEU5, NEEU50, NEEU95, CH4, CH4U5, CH4U50, CH4U95)


annualco2 <- sum(fcr_monthlysum$NEE) # gCm2
sdco2 <- sqrt(sum(fcr_monthlysum$NEE_sd^2, na.rm = TRUE))
annualch4 <- sum(fcr_monthlysum$CH4) # mgCm2
sdch4 <- sqrt(sum(fcr_monthlysum$CH4_sd^2, na.rm = TRUE))

annualco2
sdco2
annualch4
sdch4



#########

# Looking for some relationships
# temp
fcr_dailysum %>% ggplot(aes(Tmean, NEE_sum)) + 
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~Temperature~(degree~C))) +
	geom_point() + theme_bw() +
	geom_smooth(span = 0.9)

templm <- lm(log(fcr_dailysum$NEE_sum)~fcr_dailysum$Tmean)
summary(templm)

# LW

fcr_daily %>% ggplot(aes(LW_in, NEE_mean)) + 
	geom_point() + theme_bw() +
	ylab(expression(~CO[2]~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~LW[incoming]~(W~m^-2))) +
	geom_smooth(span = 0.9) 

lwco2lm <- lm(log(fcr_daily$NEE_mean)~fcr_daily$LW_in)
summary(lwco2lm)



fcr_dailysum %>% ggplot(aes(Rn, CH4_sum)) + 
	geom_point() + theme_bw() +
	ylab(expression(~CH4~(mu~mol~m^-2~s^-1))) +
	xlab(expression(~Rn~(W~m^-2~s^-1))) +
	geom_smooth(method= 'lm')

rnlm <- lm(fcr_daily$CH4_mean~fcr_daily$Rn)
summary(rnlm)
