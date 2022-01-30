### Script to clean up files produced by LiCor Eddy Pro
### 8 October 2021, A. Hounshell

### Updated with new data on 21 Jan 2022

#####################################################

# Clear workspace
rm(list = ls())

# Download/load libraries
pacman::p_load(lubridate,tidyverse)

# Set working directory
wd <- getwd()
setwd(wd)

# Read compiled file: From Eddy Pro using basic processing
# Original file from Brenda on 11 May 2021
ec <- read_csv("./Data/FCR_Eddy_up_to_2022-02-11.csv")

ec <- ec %>% 
  select(date,time,DOY,Tau,qc_Tau,H,qc_H,LE,qc_LE,co2_flux,qc_co2_flux,h2o_flux,
         qc_h2o_flux,ch4_flux,qc_ch4_flux,`co2_v-adv`,`h2o_v-adv`,`ch4_v-adv`,
         co2_molar_density,co2_mole_fraction,co2_mixing_ratio,co2_time_lag,
         co2_def_timelag,h2o_molar_density,h2o_mole_fraction,h2o_mixing_ratio,
         h2o_time_lag,h2o_def_timelag,ch4_molar_density,ch4_mole_fraction,
         ch4_mixing_ratio,ch4_time_lag,ch4_def_timelag,sonic_temperature,
         air_temperature,air_pressure,air_density,air_heat_capacity,
         air_molar_volume,ET,water_vapor_density,e,es,specific_humidity,RH,VPD,
         Tdew,wind_speed,max_wind_speed,wind_dir,`u*`,TKE,L,`(z-d)/L`,bowen_ratio,
         `T*`,x_peak,x_offset,`x_10%`,`x_30%`,`x_50%`,`x_70%`,`x_90%`,un_Tau,
         Tau_scf,un_H,H_scf,un_LE,LE_scf,un_co2_flux,co2_scf,un_h2o_flux,
         h2o_scf,un_ch4_flux,ch4_scf,u_var,v_var,w_var,rssi_77_mean) %>% 
  rename(Tau_kgms2 = Tau,
         H_wm2 = H,
         LE_wm2 = LE,
         co2_flux_umolm2s = co2_flux,
         h2o_flux_umolm2s = h2o_flux,
         ch4_flux_umolm2s = ch4_flux,
         co2_v_adv_umolm2s = `co2_v-adv`,
         h2o_v_adv_umolm2s = `h2o_v-adv`,
         ch4_v_adv_umolm2s = `ch4_v-adv`,
         co2_molar_density_mmolm3 = co2_molar_density,
         co2_mole_fraction_umolmol = co2_mole_fraction,
         co2_mixing_ratio_umolmol = co2_mixing_ratio,
         co2_time_lag_s = co2_time_lag,
         h2o_molar_density_mmolm3 = h2o_molar_density,
         h2o_mole_fraction_umolmol = h2o_mole_fraction,
         h2o_mixing_ratio_umolmol = h2o_mixing_ratio,
         h2o_time_lag_s = h2o_time_lag,
         ch4_molar_density_mmolm3 = ch4_molar_density,
         ch4_mole_fraction_umolmol = ch4_mole_fraction,
         ch4_mixing_ratio_umolmol = ch4_mixing_ratio,
         ch4_time_lag_s = ch4_time_lag,
         sonic_temperature_k = sonic_temperature,
         air_temperature_k = air_temperature,
         air_pressure_pa = air_pressure,
         air_density_kgm3 = air_density,
         air_heat_capacity_jkkg = air_heat_capacity,
         air_molar_volume_m3mol = air_molar_volume,
         ET_mmhr = ET,
         water_vapor_density_kgm3 = water_vapor_density,
         e_pa = e,
         es_pa = es,
         specific_humidity_kgkg = specific_humidity,
         VPD_pa = VPD,
         Tdew_k = Tdew,
         wind_speed_ms = wind_speed,
         max_wind_speed_ms = max_wind_speed,
         u_star_ms = `u*`,
         TKE_m2s2 = TKE,
         L_m = L,
         MO_stability = `(z-d)/L`,
         scale_T_k = `T*`,
         x_peak_m = x_peak,
         x_offset_m = x_offset,
         x_10_m = `x_10%`,
         x_30_m = `x_30%`,
         x_50_m = `x_50%`,
         x_70_m = `x_70%`,
         x_90_m = `x_90%`,
         un_Tau_kgms2 = un_Tau,
         un_H_wm2 = un_H,
         un_LE_wm2 = un_LE,
         un_co2_flux_umolm2s = un_co2_flux,
         un_h2o_flux_umolm2s = un_h2o_flux,
         un_ch4_flux_umolm2s = un_ch4_flux,
         u_var_ms = u_var,
         v_var_ms = v_var,
         w_var_ms = w_var)

# Remove -9999 and replace with NAs
ec[ec ==-9999] <- NA

## Add flag for missing data: 3 = missing data
# For: qc_tau, qc_H, qc_LE, qc_co2_flux, qc_h2o_flux, qc_ch4_flux
ec <- ec %>% 
  mutate(qc_Tau = ifelse(is.na(Tau_kgms2), 3, qc_Tau),
         qc_H = ifelse(is.na(H_wm2), 3, qc_H),
         qc_LE = ifelse(is.na(LE_wm2), 3, qc_LE),
         qc_co2_flux = ifelse(is.na(co2_flux_umolm2s), 3, qc_co2_flux),
         qc_h2o_flux = ifelse(is.na(h2o_flux_umolm2s), 3, qc_h2o_flux),
         qc_ch4_flux = ifelse(is.na(ch4_flux_umolm2s), 3, qc_ch4_flux))


# Output data
write_csv(ec, "./Data/20220130_EddyPro_Cleaned.csv")

write_csv(ec, "./EDI/20220130_EddyPro_Cleaned.csv")
