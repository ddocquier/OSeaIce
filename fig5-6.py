#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 5: Plot Arctic sea-ice area (10^6 km^2), computed via compute_area.py
    Fig. 6: Plot Arctic sea-ice volume (10^3 km^3), computed via compute_volume.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    21/10/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Options
save_fig = True
save_var = False

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/'
dir_t613 = dir_input + 'post-proc/t613/'
dir_D000 = dir_input + 'post-proc/D000/'
dir_D012 = dir_input + 'post-proc/D012/'
dir_D013 = dir_input + 'post-proc/D013/'
dir_D014 = dir_input + 'post-proc/D014/'
dir_D015 = dir_input + 'post-proc/D015/'
dir_D016 = dir_input + 'post-proc/D016/'
dir_D017 = dir_input + 'post-proc/D017/'
dir_D018 = dir_input + 'post-proc/D018/'
dir_D019 = dir_input + 'post-proc/D019/'
dir_D020 = dir_input + 'post-proc/D020/'
dir_D021 = dir_input + 'post-proc/D021/'
dir_D022 = dir_input + 'post-proc/D022/'
dir_D023 = dir_input + 'post-proc/D023/'
dir_D024 = dir_input + 'post-proc/D024/'
dir_D025 = dir_input + 'post-proc/D025/'
dir_D027 = dir_input + 'post-proc/D027/'
dir_D028 = dir_input + 'post-proc/D028/'
dir_D029 = dir_input + 'post-proc/D029/'
dir_D030 = dir_input + 'post-proc/D030/'
dir_piomas = dir_input + 'obs/'
dir_fig = dir_input + 'OSeaIce_Paper/'

# Function to retrieve monthly mean and compute seasonal cycle
def compute_monthly(nm,nmy,var):
    nyears = int(nm / nmy)
    var_mon = np.zeros((nmy,nyears))
    var_ym = np.zeros(nmy)
    sd_var = np.zeros(nmy)
    for i in np.arange(nmy):
        for k in np.arange(nyears):
            if np.isnan(var[i+k*nmy]) == False:
                var_mon[i,k] = var[i+k*nmy]
                var_ym[i] = var_ym[i] + var[i+k*nmy]
        var_ym[i] = var_ym[i] / nyears
    for i in np.arange(nmy):
        for k in np.arange(nyears):
            if np.isnan(var[i+k*nmy]) == False:
                sd_var[i] = sd_var[i] + (var[i+k*nmy] - var_ym[i])**2.
        sd_var[i] = np.sqrt(sd_var[i] / nyears)
    return var_mon,var_ym,sd_var

# Function to compute trend + SD of trend
def compute_trend(nyears,var):
    ndecades = nyears / 10
    year = np.arange(nyears) + 1
    idx = np.isfinite(year) & np.isfinite(var)
    a,b = np.polyfit(year[idx],var[idx],1)
    lreg = a * year + b
    trend_decade = a * nyears / ndecades
    trend_percent = 100 * trend_decade / lreg[0]
    n_freedom = nm-2
    s_yx = np.sum((var[idx]-lreg[idx])**2)/n_freedom
    SS_xx = np.sum((year-np.mean(year))**2)
    s_a = np.sqrt(s_yx/SS_xx)
    sd_trend = s_a * nyears / ndecades
    t_a = a / s_a
    alpha = 0.05
    t_crit = stats.t.ppf(1-alpha/2,n_freedom)
    sig_a = 0
    if np.abs(t_a) > t_crit:
        sig_a = 1
    return trend_decade,sd_trend,trend_percent,sig_a

# Function to compute significativity of difference between two means
def compute_sig(n,var1,var2):
    mean_var1 = np.nanmean(var1)
    mean_var2 = np.nanmean(var2)
    diff = mean_var1 - mean_var2
    reldiff = diff / mean_var2
    sum_var1 = 0.
    sum_var2 = 0.
    for i in np.arange(n):
        if np.isnan(var1[i]) == False:
            sum_var1 = sum_var1 + (var1[i] - mean_var1)**2.
        if np.isnan(var2[i]) == False:
            sum_var2 = sum_var2 + (var2[i] - mean_var2)**2.
    sd = np.sqrt((sum_var1 + sum_var2) / (2. * n - 2.))
    es = sd * np.sqrt(2. / n)
    t = diff / es
    alpha = 0.05
    t_crit = stats.t.ppf(1-alpha/2,(2.*n-2.))
    sig = 0
    if np.abs(t) > t_crit:
        sig = 1
    return sig

# Load sea-ice area EC-Earth D000
filename = dir_D000 + 'SIarea_D000.npy'
area_D000,area_barents_D000,area_greenland_D000,area_labrador_D000,area_laptev_D000,area_chukchi_D000,area_beaufort_D000,area_central_D000 = np.load(filename)

# Load sea-ice area EC-Earth D012
filename = dir_D012 + 'SIarea_D012.npy'
area_D012a,area_barents_D012a,area_greenland_D012a,area_labrador_D012a,area_laptev_D012a,area_chukchi_D012a,area_beaufort_D012a,area_central_D012a = np.load(filename)

# Load sea-ice area EC-Earth D012 from 2180-2229
filename = dir_D012 + 'SIarea_D012_2180.npy'
area_D012b,area_barents_D012b,area_greenland_D012b,area_labrador_D012b,area_laptev_D012b,area_chukchi_D012b,area_beaufort_D012b,area_central_D012b = np.load(filename)
area_D012 = np.concatenate((area_D012a,area_D012b),axis=0)

# Load sea-ice area EC-Earth D013
filename = dir_D013 + 'SIarea_D013.npy'
area_D013,area_barents_D013,area_greenland_D013,area_labrador_D013,area_laptev_D013,area_chukchi_D013,area_beaufort_D013,area_central_D013 = np.load(filename)

# Load sea-ice area EC-Earth D014
filename = dir_D014 + 'SIarea_D014.npy'
area_D014,area_barents_D014,area_greenland_D014,area_labrador_D014,area_laptev_D014,area_chukchi_D014,area_beaufort_D014,area_central_D014 = np.load(filename)

# Load sea-ice area EC-Earth D015
filename = dir_D015 + 'SIarea_D015.npy'
area_D015,area_barents_D015,area_greenland_D015,area_labrador_D015,area_laptev_D015,area_chukchi_D015,area_beaufort_D015,area_central_D015 = np.load(filename)

# Load sea-ice area EC-Earth D016
filename = dir_D016 + 'SIarea_D016.npy'
area_D016,area_barents_D016,area_greenland_D016,area_labrador_D016,area_laptev_D016,area_chukchi_D016,area_beaufort_D016,area_central_D016 = np.load(filename)

# Load sea-ice area EC-Earth D017
filename = dir_D017 + 'SIarea_D017.npy'
area_D017,area_barents_D017,area_greenland_D017,area_labrador_D017,area_laptev_D017,area_chukchi_D017,area_beaufort_D017,area_central_D017 = np.load(filename)

# Load sea-ice area EC-Earth D018
filename = dir_D018 + 'SIarea_D018.npy'
area_D018,area_barents_D018,area_greenland_D018,area_labrador_D018,area_laptev_D018,area_chukchi_D018,area_beaufort_D018,area_central_D018 = np.load(filename)

# Load sea-ice area EC-Earth D019
filename = dir_D019 + 'SIarea_D019.npy'
area_D019,area_barents_D019,area_greenland_D019,area_labrador_D019,area_laptev_D019,area_chukchi_D019,area_beaufort_D019,area_central_D019 = np.load(filename)

# Load sea-ice area EC-Earth D020
filename = dir_D020 + 'SIarea_D020.npy'
area_D020,area_barents_D020,area_greenland_D020,area_labrador_D020,area_laptev_D020,area_chukchi_D020,area_beaufort_D020,area_central_D020 = np.load(filename)

# Load sea-ice area EC-Earth D021
filename = dir_D021 + 'SIarea_D021.npy'
area_D021a,area_barents_D021a,area_greenland_D021a,area_labrador_D021a,area_laptev_D021a,area_chukchi_D021a,area_beaufort_D021a,area_central_D021a = np.load(filename)

# Load sea-ice area EC-Earth D021 from 2180-2229
filename = dir_D021 + 'SIarea_D021_2180.npy'
area_D021b,area_barents_D021b,area_greenland_D021b,area_labrador_D021b,area_laptev_D021b,area_chukchi_D021b,area_beaufort_D021b,area_central_D021b = np.load(filename)
area_D021 = np.concatenate((area_D021a,area_D021b),axis=0)

# Load sea-ice area EC-Earth D022
filename = dir_D022 + 'SIarea_D022.npy'
area_D022,area_barents_D022,area_greenland_D022,area_labrador_D022,area_laptev_D022,area_chukchi_D022,area_beaufort_D022,area_central_D022 = np.load(filename)

# Load sea-ice area EC-Earth D023
filename = dir_D023 + 'SIarea_D023.npy'
area_D023,area_barents_D023,area_greenland_D023,area_labrador_D023,area_laptev_D023,area_chukchi_D023,area_beaufort_D023,area_central_D023 = np.load(filename)

# Load sea-ice area EC-Earth D024
filename = dir_D024 + 'SIarea_D024.npy'
area_D024,area_barents_D024,area_greenland_D024,area_labrador_D024,area_laptev_D024,area_chukchi_D024,area_beaufort_D024,area_central_D024 = np.load(filename)

# Load sea-ice area EC-Earth D025
filename = dir_D025 + 'SIarea_D025.npy'
area_D025,area_barents_D025,area_greenland_D025,area_labrador_D025,area_laptev_D025,area_chukchi_D025,area_beaufort_D025,area_central_D025 = np.load(filename)

# Load sea-ice area EC-Earth D027
filename = dir_D027 + 'SIarea_D027.npy'
area_D027,area_barents_D027,area_greenland_D027,area_labrador_D027,area_laptev_D027,area_chukchi_D027,area_beaufort_D027,area_central_D027 = np.load(filename)

# Load sea-ice area EC-Earth D028
filename = dir_D028 + 'SIarea_D028.npy'
area_D028,area_barents_D028,area_greenland_D028,area_labrador_D028,area_laptev_D028,area_chukchi_D028,area_beaufort_D028,area_central_D028 = np.load(filename)

# Load sea-ice area EC-Earth D029
filename = dir_D029 + 'SIarea_D029.npy'
area_D029,area_barents_D029,area_greenland_D029,area_labrador_D029,area_laptev_D029,area_chukchi_D029,area_beaufort_D029,area_central_D029 = np.load(filename)

# Load sea-ice area EC-Earth D030
filename = dir_D030 + 'SIarea_D030.npy'
area_D030,area_barents_D030,area_greenland_D030,area_labrador_D030,area_laptev_D030,area_chukchi_D030,area_beaufort_D030,area_central_D030 = np.load(filename)

# Load sea-ice area EC-Earth t613
filename = dir_t613 + 'SIarea_t613_1950-2014.npy'
area,area_barents,area_greenland,area_labrador,area_laptev,area_chukchi,area_beaufort,area_central = np.load(filename)

# Load sea-ice area OSI SAF (OSI-450) 1979-2015
filename = dir_t613 + 'SIarea_OSI-450_1979-2015.npy'
area_obs,area_barents_obs,area_greenland_obs,area_labrador_obs,area_laptev_obs,area_chukchi_obs,area_beaufort_obs,area_central_obs = np.load(filename)

# Load sea-ice volume EC-Earth D000
filename = dir_D000 + 'SIvolume_D000.npy'
volume_D000,volume_barents_D000,volume_greenland_D000,volume_labrador_D000,volume_laptev_D000,volume_chukchi_D000,volume_beaufort_D000,volume_central_D000 = np.load(filename)

# Load sea-ice volume EC-Earth D012
filename = dir_D012 + 'SIvolume_D012.npy'
volume_D012a,volume_barents_D012a,volume_greenland_D012a,volume_labrador_D012a,volume_laptev_D012a,volume_chukchi_D012a,volume_beaufort_D012a,volume_central_D012a = np.load(filename)

# Load sea-ice volume EC-Earth D012 from 2180-2229
filename = dir_D012 + 'SIvolume_D012_2180.npy'
volume_D012b,volume_barents_D012b,volume_greenland_D012b,volume_labrador_D012b,volume_laptev_D012b,volume_chukchi_D012b,volume_beaufort_D012b,volume_central_D012b = np.load(filename)
volume_D012 = np.concatenate((volume_D012a,volume_D012b),axis=0)

# Load sea-ice volume EC-Earth D013
filename = dir_D013 + 'SIvolume_D013.npy'
volume_D013,volume_barents_D013,volume_greenland_D013,volume_labrador_D013,volume_laptev_D013,volume_chukchi_D013,volume_beaufort_D013,volume_central_D013 = np.load(filename)

# Load sea-ice volume EC-Earth D014
filename = dir_D014 + 'SIvolume_D014.npy'
volume_D014,volume_barents_D014,volume_greenland_D014,volume_labrador_D014,volume_laptev_D014,volume_chukchi_D014,volume_beaufort_D014,volume_central_D014 = np.load(filename)

# Load sea-ice volume EC-Earth D015
filename = dir_D015 + 'SIvolume_D015.npy'
volume_D015,volume_barents_D015,volume_greenland_D015,volume_labrador_D015,volume_laptev_D015,volume_chukchi_D015,volume_beaufort_D015,volume_central_D015 = np.load(filename)

# Load sea-ice volume EC-Earth D016
filename = dir_D016 + 'SIvolume_D016.npy'
volume_D016,volume_barents_D016,volume_greenland_D016,volume_labrador_D016,volume_laptev_D016,volume_chukchi_D016,volume_beaufort_D016,volume_central_D016 = np.load(filename)

# Load sea-ice volume EC-Earth D017
filename = dir_D017 + 'SIvolume_D017.npy'
volume_D017,volume_barents_D017,volume_greenland_D017,volume_labrador_D017,volume_laptev_D017,volume_chukchi_D017,volume_beaufort_D017,volume_central_D017 = np.load(filename)

# Load sea-ice volume EC-Earth D018
filename = dir_D018 + 'SIvolume_D018.npy'
volume_D018,volume_barents_D018,volume_greenland_D018,volume_labrador_D018,volume_laptev_D018,volume_chukchi_D018,volume_beaufort_D018,volume_central_D018 = np.load(filename)

# Load sea-ice volume EC-Earth D019
filename = dir_D019 + 'SIvolume_D019.npy'
volume_D019,volume_barents_D019,volume_greenland_D019,volume_labrador_D019,volume_laptev_D019,volume_chukchi_D019,volume_beaufort_D019,volume_central_D019 = np.load(filename)

# Load sea-ice volume EC-Earth D020
filename = dir_D020 + 'SIvolume_D020.npy'
volume_D020,volume_barents_D020,volume_greenland_D020,volume_labrador_D020,volume_laptev_D020,volume_chukchi_D020,volume_beaufort_D020,volume_central_D020 = np.load(filename)

# Load sea-ice volume EC-Earth D021
filename = dir_D021 + 'SIvolume_D021.npy'
volume_D021a,volume_barents_D021a,volume_greenland_D021a,volume_labrador_D021a,volume_laptev_D021a,volume_chukchi_D021a,volume_beaufort_D021a,volume_central_D021a = np.load(filename)

# Load sea-ice volume EC-Earth D021 from 2180-2229
filename = dir_D021 + 'SIvolume_D021_2180.npy'
volume_D021b,volume_barents_D021b,volume_greenland_D021b,volume_labrador_D021b,volume_laptev_D021b,volume_chukchi_D021b,volume_beaufort_D021b,volume_central_D021b = np.load(filename)
volume_D021 = np.concatenate((volume_D021a,volume_D021b),axis=0)

# Load sea-ice volume EC-Earth D022
filename = dir_D022 + 'SIvolume_D022.npy'
volume_D022,volume_barents_D022,volume_greenland_D022,volume_labrador_D022,volume_laptev_D022,volume_chukchi_D022,volume_beaufort_D022,volume_central_D022 = np.load(filename)

# Load sea-ice volume EC-Earth D023
filename = dir_D023 + 'SIvolume_D023.npy'
volume_D023,volume_barents_D023,volume_greenland_D023,volume_labrador_D023,volume_laptev_D023,volume_chukchi_D023,volume_beaufort_D023,volume_central_D023 = np.load(filename)

# Load sea-ice volume EC-Earth D024
filename = dir_D024 + 'SIvolume_D024.npy'
volume_D024,volume_barents_D024,volume_greenland_D024,volume_labrador_D024,volume_laptev_D024,volume_chukchi_D024,volume_beaufort_D024,volume_central_D024 = np.load(filename)

# Load sea-ice volume EC-Earth D025
filename = dir_D025 + 'SIvolume_D025.npy'
volume_D025,volume_barents_D025,volume_greenland_D025,volume_labrador_D025,volume_laptev_D025,volume_chukchi_D025,volume_beaufort_D025,volume_central_D025 = np.load(filename)

# Load sea-ice volume EC-Earth D027
filename = dir_D027 + 'SIvolume_D027.npy'
volume_D027,volume_barents_D027,volume_greenland_D027,volume_labrador_D027,volume_laptev_D027,volume_chukchi_D027,volume_beaufort_D027,volume_central_D027 = np.load(filename)

# Load sea-ice volume EC-Earth D028
filename = dir_D028 + 'SIvolume_D028.npy'
volume_D028,volume_barents_D028,volume_greenland_D028,volume_labrador_D028,volume_laptev_D028,volume_chukchi_D028,volume_beaufort_D028,volume_central_D028 = np.load(filename)

# Load sea-ice volume EC-Earth D029
filename = dir_D029 + 'SIvolume_D029.npy'
volume_D029,volume_barents_D029,volume_greenland_D029,volume_labrador_D029,volume_laptev_D029,volume_chukchi_D029,volume_beaufort_D029,volume_central_D029 = np.load(filename)

# Load sea-ice volume EC-Earth D030
filename = dir_D030 + 'SIvolume_D030.npy'
volume_D030,volume_barents_D030,volume_greenland_D030,volume_labrador_D030,volume_laptev_D030,volume_chukchi_D030,volume_beaufort_D030,volume_central_D030 = np.load(filename)

# Load sea-ice volume EC-Earth t613
filename = dir_t613 + 'SIvolume_t613_1950-2014.npy'
volume,volume_barents,volume_greenland,volume_labrador,volume_laptev,volume_chukchi,volume_beaufort,volume_central = np.load(filename)

# Time parameters
nmy = 12
nm = np.size(area)
nm_D000 = np.size(area_D000)
nm_D012 = np.size(area_D012)
nm_D013 = np.size(area_D013)
nm_obs = np.size(area_obs)
nyears = int(nm / nmy)
nyears_D000 = int(nm_D000 / nmy)
nyears_D012 = int(nm_D012 / nmy)
nyears_D013 = int(nm_D013 / nmy)
nyears_obs = int(nm_obs / nmy)

# Load sea-ice volume PIOMAS 1979-2019 (2 last months missing at the time of writing)
year,vol1,vol2,vol3,vol4,vol5,vol6,vol7,vol8,vol9,vol10,vol11,vol12 = np.loadtxt(dir_piomas+'PIOMAS.2sst.monthly.Current.v2.1.txt',unpack=True)
nyears_piomas = np.size(year)
nm_piomas = nyears_piomas*nmy
volume_piomas = np.zeros(nm_piomas)
k = 0
for j in np.arange(nyears_piomas):
    volume_piomas[k] = vol1[j]
    volume_piomas[k+1] = vol2[j]
    volume_piomas[k+2] = vol3[j]
    volume_piomas[k+3] = vol4[j]
    volume_piomas[k+4] = vol5[j]
    volume_piomas[k+5] = vol6[j]
    volume_piomas[k+6] = vol7[j]
    volume_piomas[k+7] = vol8[j]
    volume_piomas[k+8] = vol9[j]
    volume_piomas[k+9] = vol10[j]
    volume_piomas[k+10] = vol11[j]
    volume_piomas[k+11] = vol12[j]
    k = k + nmy
volume_piomas[-1] = np.nan
volume_piomas[-2] = np.nan

# Retrieve monthly mean SIA EC-Earth
area_mon,area_ym,sd_area = compute_monthly(nm,nmy,area) # t613 1950-2014
nm_reduced = int((2014-1979+1)*12)
begin_reduced = int((1979-1950)*12)
notused,area_ym2,sd_area2 = compute_monthly(nm_reduced,nmy,area[begin_reduced::]) # t613 1979-2014
area_mon_D000,notused,notused = compute_monthly(nm_D000,nmy,area_D000) # 2014-2213
begin_reduced_D000 = int((2130-2014)*12)
area_mon_D000b,area_ym_D000,sd_area_D000 = compute_monthly(nm_D013,nmy,area_D000[begin_reduced_D000:begin_reduced_D000+nm_D013]) # 2130-2179
area_mon_D012,notused,notused = compute_monthly(nm_D012,nmy,area_D012)
area_mon_D012b,area_ym_D012,sd_area_D012 = compute_monthly(nm_D013,nmy,area_D012a)
area_mon_D013,area_ym_D013,sd_area_D013 = compute_monthly(nm_D013,nmy,area_D013)
area_mon_D014,area_ym_D014,sd_area_D014 = compute_monthly(nm_D013,nmy,area_D014)
area_mon_D015,area_ym_D015,sd_area_D015 = compute_monthly(nm_D013,nmy,area_D015)
area_mon_D016,area_ym_D016,sd_area_D016 = compute_monthly(nm_D013,nmy,area_D016)
area_mon_D017,area_ym_D017,sd_area_D017 = compute_monthly(nm_D013,nmy,area_D017)
area_mon_D018,area_ym_D018,sd_area_D018 = compute_monthly(nm_D013,nmy,area_D018)
area_mon_D019,area_ym_D019,sd_area_D019 = compute_monthly(nm_D013,nmy,area_D019)
area_mon_D020,area_ym_D020,sd_area_D020 = compute_monthly(nm_D013,nmy,area_D020)
area_mon_D021,notused,notused = compute_monthly(nm_D012,nmy,area_D021)
area_mon_D021b,area_ym_D021,sd_area_D021 = compute_monthly(nm_D013,nmy,area_D021a)
area_mon_D022,area_ym_D022,sd_area_D022 = compute_monthly(nm_D013,nmy,area_D022)
area_mon_D023,area_ym_D023,sd_area_D023 = compute_monthly(nm_D013,nmy,area_D023)
area_mon_D024,area_ym_D024,sd_area_D024 = compute_monthly(nm_D013,nmy,area_D024)
area_mon_D025,area_ym_D025,sd_area_D025 = compute_monthly(nm_D013,nmy,area_D025)
area_mon_D027,area_ym_D027,sd_area_D027 = compute_monthly(nm_D013,nmy,area_D027)
area_mon_D028,area_ym_D028,sd_area_D028 = compute_monthly(nm_D013,nmy,area_D028)
area_mon_D029,area_ym_D029,sd_area_D029 = compute_monthly(nm_D013,nmy,area_D029)
area_mon_D030,area_ym_D030,sd_area_D030 = compute_monthly(nm_D013,nmy,area_D030)

# Retrieve monthly mean SIV EC-Earth
volume_mon,volume_ym,sd_volume = compute_monthly(nm,nmy,volume) # 1950-2014
notused,volume_ym2,sd_volume2 = compute_monthly(nm_reduced,nmy,volume[begin_reduced::]) # 1979-2014
volume_mon_D000,notused,notused = compute_monthly(nm_D000,nmy,volume_D000) # 2014-2213
volume_mon_D000b,volume_ym_D000,sd_volume_D000 = compute_monthly(nm_D013,nmy,volume_D000[begin_reduced_D000:begin_reduced_D000+nm_D013]) # 2130-2179
volume_mon_D012,notused,notused = compute_monthly(nm_D012,nmy,volume_D012)
volume_mon_D012b,volume_ym_D012,sd_volume_D012 = compute_monthly(nm_D013,nmy,volume_D012a)
volume_mon_D013,volume_ym_D013,sd_volume_D013 = compute_monthly(nm_D013,nmy,volume_D013)
volume_mon_D014,volume_ym_D014,sd_volume_D014 = compute_monthly(nm_D013,nmy,volume_D014)
volume_mon_D015,volume_ym_D015,sd_volume_D015 = compute_monthly(nm_D013,nmy,volume_D015)
volume_mon_D016,volume_ym_D016,sd_volume_D016 = compute_monthly(nm_D013,nmy,volume_D016)
volume_mon_D017,volume_ym_D017,sd_volume_D017 = compute_monthly(nm_D013,nmy,volume_D017)
volume_mon_D018,volume_ym_D018,sd_volume_D018 = compute_monthly(nm_D013,nmy,volume_D018)
volume_mon_D019,volume_ym_D019,sd_volume_D019 = compute_monthly(nm_D013,nmy,volume_D019)
volume_mon_D020,volume_ym_D020,sd_volume_D020 = compute_monthly(nm_D013,nmy,volume_D020)
volume_mon_D021,notused,notused = compute_monthly(nm_D012,nmy,volume_D021)
volume_mon_D021b,volume_ym_D021,sd_volume_D021 = compute_monthly(nm_D013,nmy,volume_D021a)
volume_mon_D022,volume_ym_D022,sd_volume_D022 = compute_monthly(nm_D013,nmy,volume_D022)
volume_mon_D023,volume_ym_D023,sd_volume_D023 = compute_monthly(nm_D013,nmy,volume_D023)
volume_mon_D024,volume_ym_D024,sd_volume_D024 = compute_monthly(nm_D013,nmy,volume_D024)
volume_mon_D025,volume_ym_D025,sd_volume_D025 = compute_monthly(nm_D013,nmy,volume_D025)
volume_mon_D027,volume_ym_D027,sd_volume_D027 = compute_monthly(nm_D013,nmy,volume_D027)
volume_mon_D028,volume_ym_D028,sd_volume_D028 = compute_monthly(nm_D013,nmy,volume_D028)
volume_mon_D029,volume_ym_D029,sd_volume_D029 = compute_monthly(nm_D013,nmy,volume_D029)
volume_mon_D030,volume_ym_D030,sd_volume_D030 = compute_monthly(nm_D013,nmy,volume_D030)

# Retrieve monthly mean SIA Barents/Kara Sea
area_mon_barents_D000,notused,notused = compute_monthly(nm_D000,nmy,area_barents_D000)
area_mon_barents_D012,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D012a)
area_mon_barents_D013,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D013)
area_mon_barents_D014,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D014)
area_mon_barents_D015,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D015)
area_mon_barents_D016,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D016)
area_mon_barents_D017,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D017)
area_mon_barents_D018,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D018)
area_mon_barents_D019,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D019)
area_mon_barents_D020,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D020)
area_mon_barents_D021,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D021a)
area_mon_barents_D022,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D022)
area_mon_barents_D023,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D023)
area_mon_barents_D024,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D024)
area_mon_barents_D025,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D025)
area_mon_barents_D027,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D027)
area_mon_barents_D028,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D028)
area_mon_barents_D029,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D029)
area_mon_barents_D030,notused,notused = compute_monthly(nm_D013,nmy,area_barents_D030)

# Retrieve monthly mean SIA OSI SAF
area_obs_mon,area_obs_ym,sd_area_obs = compute_monthly(nm_obs,nmy,area_obs) # 1979-2015
nm_obs_reduced = int((2014-1979+1)*12)
notused,area_obs_ym2,sd_area_obs2 = compute_monthly(nm_obs_reduced,nmy,area_obs[0:nm_obs_reduced]) # 1979-2014

# Retrieve monthly mean SIV PIOMAS
volume_piomas_mon,volume_piomas_ym,sd_volume_piomas = compute_monthly(nm_piomas,nmy,volume_piomas) # 1979-2019 (missing 11 and 12/2019)
notused,volume_piomas_ym2,sd_volume_piomas2 = compute_monthly(nm_obs_reduced,nmy,volume_piomas[0:nm_obs_reduced]) # 1979-2014

# Save March/September SIA and SIV
if save_var == True:
    filename = dir_fig + 'SIA_March.npy'
    np.save(filename,[area_mon_D000b[2,:],area_mon_D012b[2,:],area_mon_D013[2,:],area_mon_D014[2,:],area_mon_D015[2,:],area_mon_D016[2,:],area_mon_D017[2,:],area_mon_D018[2,:],area_mon_D019[2,:],area_mon_D020[2,:],area_mon_D021b[2,:],area_mon_D022[2,:],area_mon_D023[2,:],area_mon_D024[2,:],area_mon_D025[2,:],area_mon_D027[2,:],area_mon_D028[2,:],area_mon_D029[2,:],area_mon_D030[2,:]])
    filename = dir_fig + 'SIV_March.npy'
    np.save(filename,[volume_mon_D000b[2,:],volume_mon_D012b[2,:],volume_mon_D013[2,:],volume_mon_D014[2,:],volume_mon_D015[2,:],volume_mon_D016[2,:],volume_mon_D017[2,:],volume_mon_D018[2,:],volume_mon_D019[2,:],volume_mon_D020[2,:],volume_mon_D021b[2,:],volume_mon_D022[2,:],volume_mon_D023[2,:],volume_mon_D024[2,:],volume_mon_D025[2,:],volume_mon_D027[2,:],volume_mon_D028[2,:],volume_mon_D029[2,:],volume_mon_D030[2,:]])
    filename = dir_fig + 'SIA_September.npy'
    np.save(filename,[area_mon_D000b[8,:],area_mon_D012b[8,:],area_mon_D013[8,:],area_mon_D014[8,:],area_mon_D015[8,:],area_mon_D016[8,:],area_mon_D017[8,:],area_mon_D018[8,:],area_mon_D019[8,:],area_mon_D020[8,:],area_mon_D021b[8,:],area_mon_D022[8,:],area_mon_D023[8,:],area_mon_D024[8,:],area_mon_D025[8,:],area_mon_D027[8,:],area_mon_D028[8,:],area_mon_D029[8,:],area_mon_D030[8,:]])
    filename = dir_fig + 'SIV_September.npy'
    np.save(filename,[volume_mon_D000b[8,:],volume_mon_D012b[8,:],volume_mon_D013[8,:],volume_mon_D014[8,:],volume_mon_D015[8,:],volume_mon_D016[8,:],volume_mon_D017[8,:],volume_mon_D018[8,:],volume_mon_D019[8,:],volume_mon_D020[8,:],volume_mon_D021b[8,:],volume_mon_D022[8,:],volume_mon_D023[8,:],volume_mon_D024[8,:],volume_mon_D025[8,:],volume_mon_D027[8,:],volume_mon_D028[8,:],volume_mon_D029[8,:],volume_mon_D030[8,:]])

###########
# Table 2 #
###########

# Print mean numbers - March Arctic SIA
month = 2
begin_reduced = int(2130-2014)
print('March Arctic SIA (10^6 km^2)')
print('CTRL =',np.round(area_ym_D000[month],1),'+/-',np.round(sd_area_D000[month],1))
print('Obs =',np.round(area_obs_ym[month],1),'+/-',np.round(sd_area_obs[month],1))
print('ATL1+1K =',np.round(area_ym_D013[month],1),'+/-',np.round(sd_area_D013[month],1),compute_sig(nyears_D013,area_mon_D013[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+3K =',np.round(area_ym_D012[month],1),'+/-',np.round(sd_area_D012[month],1),compute_sig(nyears_D013,area_mon_D012[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+5K =',np.round(area_ym_D014[month],1),'+/-',np.round(sd_area_D014[month],1),compute_sig(nyears_D013,area_mon_D014[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+1K =',np.round(area_ym_D016[month],1),'+/-',np.round(sd_area_D016[month],1),compute_sig(nyears_D013,area_mon_D016[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+3K =',np.round(area_ym_D015[month],1),'+/-',np.round(sd_area_D015[month],1),compute_sig(nyears_D013,area_mon_D015[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+5K =',np.round(area_ym_D017[month],1),'+/-',np.round(sd_area_D017[month],1),compute_sig(nyears_D013,area_mon_D017[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+1K =',np.round(area_ym_D019[month],1),'+/-',np.round(sd_area_D019[month],1),compute_sig(nyears_D013,area_mon_D019[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+3K =',np.round(area_ym_D018[month],1),'+/-',np.round(sd_area_D018[month],1),compute_sig(nyears_D013,area_mon_D018[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+5K =',np.round(area_ym_D020[month],1),'+/-',np.round(sd_area_D020[month],1),compute_sig(nyears_D013,area_mon_D020[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+1K =',np.round(area_ym_D027[month],1),'+/-',np.round(sd_area_D027[month],1),compute_sig(nyears_D013,area_mon_D027[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+3K =',np.round(area_ym_D021[month],1),'+/-',np.round(sd_area_D021[month],1),compute_sig(nyears_D013,area_mon_D021[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+5K =',np.round(area_ym_D028[month],1),'+/-',np.round(sd_area_D028[month],1),compute_sig(nyears_D013,area_mon_D028[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+1K =',np.round(area_ym_D029[month],1),'+/-',np.round(sd_area_D029[month],1),compute_sig(nyears_D013,area_mon_D029[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K =',np.round(area_ym_D022[month],1),'+/-',np.round(sd_area_D022[month],1),compute_sig(nyears_D013,area_mon_D022[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+5K =',np.round(area_ym_D030[month],1),'+/-',np.round(sd_area_D030[month],1),compute_sig(nyears_D013,area_mon_D030[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+1K =',np.round(area_ym_D024[month],1),'+/-',np.round(sd_area_D024[month],1),compute_sig(nyears_D013,area_mon_D024[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+3K =',np.round(area_ym_D023[month],1),'+/-',np.round(sd_area_D023[month],1),compute_sig(nyears_D013,area_mon_D023[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+5K =',np.round(area_ym_D025[month],1),'+/-',np.round(sd_area_D025[month],1),compute_sig(nyears_D013,area_mon_D025[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K - PAC3+3K =',compute_sig(nyears_D013,area_mon_D022[month,:],area_mon_D023[month,:]))
print('ATL2+3K - ATL2+1K =',compute_sig(nyears_D013,area_mon_D015[month,:],area_mon_D016[month,:]))
print('ATL3+3K - ATL3+1K =',compute_sig(nyears_D013,area_mon_D018[month,:],area_mon_D019[month,:]))
print('PAC3+3K - PAC3+5K =',compute_sig(nyears_D013,area_mon_D023[month,:],area_mon_D025[month,:]))
print('------------')

# Print mean numbers - September Arctic SIA
month = 8
print('September Arctic SIA (10^6 km^2)')
print('CTRL =',np.round(area_ym_D000[month],1),'+/-',np.round(sd_area_D000[month],1))
print('Obs =',np.round(area_obs_ym[month],1),'+/-',np.round(sd_area_obs[month],1))
print('ATL1+1K =',np.round(area_ym_D013[month],1),'+/-',np.round(sd_area_D013[month],1),compute_sig(nyears_D013,area_mon_D013[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+3K =',np.round(area_ym_D012[month],1),'+/-',np.round(sd_area_D012[month],1),compute_sig(nyears_D013,area_mon_D012[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+5K =',np.round(area_ym_D014[month],1),'+/-',np.round(sd_area_D014[month],1),compute_sig(nyears_D013,area_mon_D014[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+1K =',np.round(area_ym_D016[month],1),'+/-',np.round(sd_area_D016[month],1),compute_sig(nyears_D013,area_mon_D016[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+3K =',np.round(area_ym_D015[month],1),'+/-',np.round(sd_area_D015[month],1),compute_sig(nyears_D013,area_mon_D015[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+5K =',np.round(area_ym_D017[month],1),'+/-',np.round(sd_area_D017[month],1),compute_sig(nyears_D013,area_mon_D017[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+1K =',np.round(area_ym_D019[month],1),'+/-',np.round(sd_area_D019[month],1),compute_sig(nyears_D013,area_mon_D019[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+3K =',np.round(area_ym_D018[month],1),'+/-',np.round(sd_area_D018[month],1),compute_sig(nyears_D013,area_mon_D018[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+5K =',np.round(area_ym_D020[month],1),'+/-',np.round(sd_area_D020[month],1),compute_sig(nyears_D013,area_mon_D020[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+1K =',np.round(area_ym_D027[month],1),'+/-',np.round(sd_area_D027[month],1),compute_sig(nyears_D013,area_mon_D027[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+3K =',np.round(area_ym_D021[month],1),'+/-',np.round(sd_area_D021[month],1),compute_sig(nyears_D013,area_mon_D021[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+5K =',np.round(area_ym_D028[month],1),'+/-',np.round(sd_area_D028[month],1),compute_sig(nyears_D013,area_mon_D028[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+1K =',np.round(area_ym_D029[month],1),'+/-',np.round(sd_area_D029[month],1),compute_sig(nyears_D013,area_mon_D029[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K =',np.round(area_ym_D022[month],1),'+/-',np.round(sd_area_D022[month],1),compute_sig(nyears_D013,area_mon_D022[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+5K =',np.round(area_ym_D030[month],1),'+/-',np.round(sd_area_D030[month],1),compute_sig(nyears_D013,area_mon_D030[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+1K =',np.round(area_ym_D024[month],1),'+/-',np.round(sd_area_D024[month],1),compute_sig(nyears_D013,area_mon_D024[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+3K =',np.round(area_ym_D023[month],1),'+/-',np.round(sd_area_D023[month],1),compute_sig(nyears_D013,area_mon_D023[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+5K =',np.round(area_ym_D025[month],1),'+/-',np.round(sd_area_D025[month],1),compute_sig(nyears_D013,area_mon_D025[month,:],area_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K - PAC3+3K =',compute_sig(nyears_D013,area_mon_D022[month,:],area_mon_D023[month,:]))
print('ATL2+3K - ATL2+1K =',compute_sig(nyears_D013,area_mon_D015[month,:],area_mon_D016[month,:]))
print('ATL3+3K - ATL3+1K =',compute_sig(nyears_D013,area_mon_D018[month,:],area_mon_D019[month,:]))
print('PAC3+3K - PAC3+5K =',compute_sig(nyears_D013,area_mon_D023[month,:],area_mon_D025[month,:]))
print('------------')

# Print mean numbers - March Arctic SIV
month = 2
print('March Arctic SIV (10^3 km^3)')
print('CTRL =',np.round(volume_ym_D000[month],1),'+/-',np.round(sd_volume_D000[month],1))
print('PIOMAS =',np.round(volume_piomas_ym[month],1),'+/-',np.round(sd_volume_piomas[month],1))
print('ATL1+1K =',np.round(volume_ym_D013[month],1),'+/-',np.round(sd_volume_D013[month],1),compute_sig(nyears_D013,volume_mon_D013[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+3K =',np.round(volume_ym_D012[month],1),'+/-',np.round(sd_volume_D012[month],1),compute_sig(nyears_D013,volume_mon_D012[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+5K =',np.round(volume_ym_D014[month],1),'+/-',np.round(sd_volume_D014[month],1),compute_sig(nyears_D013,volume_mon_D014[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+1K =',np.round(volume_ym_D016[month],1),'+/-',np.round(sd_volume_D016[month],1),compute_sig(nyears_D013,volume_mon_D016[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+3K =',np.round(volume_ym_D015[month],1),'+/-',np.round(sd_volume_D015[month],1),compute_sig(nyears_D013,volume_mon_D015[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+5K =',np.round(volume_ym_D017[month],1),'+/-',np.round(sd_volume_D017[month],1),compute_sig(nyears_D013,volume_mon_D017[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+1K =',np.round(volume_ym_D019[month],1),'+/-',np.round(sd_volume_D019[month],1),compute_sig(nyears_D013,volume_mon_D019[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+3K =',np.round(volume_ym_D018[month],1),'+/-',np.round(sd_volume_D018[month],1),compute_sig(nyears_D013,volume_mon_D018[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+5K =',np.round(volume_ym_D020[month],1),'+/-',np.round(sd_volume_D020[month],1),compute_sig(nyears_D013,volume_mon_D020[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+1K =',np.round(volume_ym_D027[month],1),'+/-',np.round(sd_volume_D027[month],1),compute_sig(nyears_D013,volume_mon_D027[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+3K =',np.round(volume_ym_D021[month],1),'+/-',np.round(sd_volume_D021[month],1),compute_sig(nyears_D013,volume_mon_D021[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+5K =',np.round(volume_ym_D028[month],1),'+/-',np.round(sd_volume_D028[month],1),compute_sig(nyears_D013,volume_mon_D028[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+1K =',np.round(volume_ym_D029[month],1),'+/-',np.round(sd_volume_D029[month],1),compute_sig(nyears_D013,volume_mon_D029[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K =',np.round(volume_ym_D022[month],1),'+/-',np.round(sd_volume_D022[month],1),compute_sig(nyears_D013,volume_mon_D022[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+5K =',np.round(volume_ym_D030[month],1),'+/-',np.round(sd_volume_D030[month],1),compute_sig(nyears_D013,volume_mon_D030[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+1K =',np.round(volume_ym_D024[month],1),'+/-',np.round(sd_volume_D024[month],1),compute_sig(nyears_D013,volume_mon_D024[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+3K =',np.round(volume_ym_D023[month],1),'+/-',np.round(sd_volume_D023[month],1),compute_sig(nyears_D013,volume_mon_D023[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+5K =',np.round(volume_ym_D025[month],1),'+/-',np.round(sd_volume_D025[month],1),compute_sig(nyears_D013,volume_mon_D025[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K - PAC3+3K =',compute_sig(nyears_D013,volume_mon_D022[month,:],volume_mon_D023[month,:]))
print('ATL2+3K - ATL2+1K =',compute_sig(nyears_D013,volume_mon_D015[month,:],volume_mon_D016[month,:]))
print('ATL3+3K - ATL3+1K =',compute_sig(nyears_D013,volume_mon_D018[month,:],volume_mon_D019[month,:]))
print('PAC3+3K - PAC3+5K =',compute_sig(nyears_D013,volume_mon_D023[month,:],volume_mon_D025[month,:]))
print('------------')

# Print mean numbers - September Arctic SIV
month = 8
print('September Arctic SIV (10^3 km^3)')
print('CTRL =',np.round(volume_ym_D000[month],1),'+/-',np.round(sd_volume_D000[month],1))
print('PIOMAS =',np.round(volume_piomas_ym[month],1),'+/-',np.round(sd_volume_piomas[month],1))
print('ATL1+1K =',np.round(volume_ym_D013[month],1),'+/-',np.round(sd_volume_D013[month],1),compute_sig(nyears_D013,volume_mon_D013[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+3K =',np.round(volume_ym_D012[month],1),'+/-',np.round(sd_volume_D012[month],1),compute_sig(nyears_D013,volume_mon_D012[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL1+5K =',np.round(volume_ym_D014[month],1),'+/-',np.round(sd_volume_D014[month],1),compute_sig(nyears_D013,volume_mon_D014[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+1K =',np.round(volume_ym_D016[month],1),'+/-',np.round(sd_volume_D016[month],1),compute_sig(nyears_D013,volume_mon_D016[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+3K =',np.round(volume_ym_D015[month],1),'+/-',np.round(sd_volume_D015[month],1),compute_sig(nyears_D013,volume_mon_D015[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL2+5K =',np.round(volume_ym_D017[month],1),'+/-',np.round(sd_volume_D017[month],1),compute_sig(nyears_D013,volume_mon_D017[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+1K =',np.round(volume_ym_D019[month],1),'+/-',np.round(sd_volume_D019[month],1),compute_sig(nyears_D013,volume_mon_D019[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+3K =',np.round(volume_ym_D018[month],1),'+/-',np.round(sd_volume_D018[month],1),compute_sig(nyears_D013,volume_mon_D018[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('ATL3+5K =',np.round(volume_ym_D020[month],1),'+/-',np.round(sd_volume_D020[month],1),compute_sig(nyears_D013,volume_mon_D020[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+1K =',np.round(volume_ym_D027[month],1),'+/-',np.round(sd_volume_D027[month],1),compute_sig(nyears_D013,volume_mon_D027[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+3K =',np.round(volume_ym_D021[month],1),'+/-',np.round(sd_volume_D021[month],1),compute_sig(nyears_D013,volume_mon_D021[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC1+5K =',np.round(volume_ym_D028[month],1),'+/-',np.round(sd_volume_D028[month],1),compute_sig(nyears_D013,volume_mon_D028[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+1K =',np.round(volume_ym_D029[month],1),'+/-',np.round(sd_volume_D029[month],1),compute_sig(nyears_D013,volume_mon_D029[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K =',np.round(volume_ym_D022[month],1),'+/-',np.round(sd_volume_D022[month],1),compute_sig(nyears_D013,volume_mon_D022[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+5K =',np.round(volume_ym_D030[month],1),'+/-',np.round(sd_volume_D030[month],1),compute_sig(nyears_D013,volume_mon_D030[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+1K =',np.round(volume_ym_D024[month],1),'+/-',np.round(sd_volume_D024[month],1),compute_sig(nyears_D013,volume_mon_D024[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+3K =',np.round(volume_ym_D023[month],1),'+/-',np.round(sd_volume_D023[month],1),compute_sig(nyears_D013,volume_mon_D023[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC3+5K =',np.round(volume_ym_D025[month],1),'+/-',np.round(sd_volume_D025[month],1),compute_sig(nyears_D013,volume_mon_D025[month,:],volume_mon_D000[month,begin_reduced:begin_reduced+nyears_D013]))
print('PAC2+3K - PAC3+3K =',compute_sig(nyears_D013,volume_mon_D022[month,:],volume_mon_D023[month,:]))
print('ATL2+3K - ATL2+1K =',compute_sig(nyears_D013,volume_mon_D015[month,:],volume_mon_D016[month,:]))
print('ATL3+3K - ATL3+1K =',compute_sig(nyears_D013,volume_mon_D018[month,:],volume_mon_D019[month,:]))
print('PAC3+3K - PAC3+5K =',compute_sig(nyears_D013,volume_mon_D023[month,:],volume_mon_D025[month,:]))
print('------------')


# Compute trend in March sea-ice area
month = 2
trend_area,sd_area,trend_area_percent,sig_area = compute_trend(nyears,area_mon[month,:])
trend_area_D000,sd_area_D000,trend_area_percent_D000,sig_area_D000 = compute_trend(nyears_D000,area_mon_D000[month,:])
trend_area_obs,sd_area_obs,trend_area_percent_obs,sig_area_obs = compute_trend(nyears_obs,area_obs_mon[month,:])
print('Trend in March sea-ice area (10^6 km^2/decade)')
print('hist =',trend_area,'(',sig_area,')')
print('CTRL =',trend_area_D000,'(',sig_area_D000,')')
print('OSI SAF =',trend_area_obs,'(',sig_area_obs,')')

# Compute trend in September sea-ice area
month = 8
trend_area,sd_area,trend_area_percent,sig_area = compute_trend(nyears,area_mon[month,:])
trend_area_D000,sd_area_D000,trend_area_percent_D000,sig_area_D000 = compute_trend(nyears_D000,area_mon_D000[month,:])
trend_area_obs,sd_area_obs,trend_area_percent_obs,sig_area_obs = compute_trend(nyears_obs,area_obs_mon[month,:])
print('Trend in September sea-ice area (10^6 km^2/decade)')
print('hist =',trend_area,'(',sig_area,')')
print('CTRL =',trend_area_D000,'(',sig_area_D000,')')
print('OSI SAF =',trend_area_obs,'(',sig_area_obs,')')

# Compute trend in March sea-ice volume
month = 2
trend_volume,sd_volume,trend_volume_percent,sig_volume = compute_trend(nyears,volume_mon[month,:])
trend_volume_D000,sd_volume_D000,trend_volume_percent_D000,sig_volume_D000 = compute_trend(nyears_D000,volume_mon_D000[month,:])
trend_volume_piomas,sd_volume_piomas,trend_volume_percent_piomas,sig_volume_piomas = compute_trend(nyears_piomas,volume_piomas_mon[month,:])
print('Trend in March sea-ice volume (10^3 km^3/decade)')
print('hist =',trend_volume,'(',sig_volume,')')
print('CTRL =',trend_volume_D000,'(',sig_volume_D000,')')
print('PIOMAS =',trend_volume_piomas,'(',sig_volume_piomas,')')

# Compute trend in September sea-ice volume
month = 8
trend_volume,sd_volume,trend_volume_percent,sig_volume = compute_trend(nyears,volume_mon[month,:])
trend_volume_D000,sd_volume_D000,trend_volume_percent_D000,sig_volume_D000 = compute_trend(nyears_D000,volume_mon_D000[month,:])
trend_volume_piomas,sd_volume_piomas,trend_volume_percent_piomas,sig_volume_piomas = compute_trend(nyears_piomas,volume_piomas_mon[month,:])
print('Trend in September sea-ice volume (10^3 km^3/decade)')
print('hist =',trend_volume,'(',sig_volume,')')
print('CTRL =',trend_volume_D000,'(',sig_volume_D000,')')
print('PIOMAS =',trend_volume_piomas,'(',sig_volume_piomas,')')

# Labels
#name_xticks = ['1960','1980','2000','2020','2040','2060','2080','2100','2120','2140','2160','2180','2200','2220']
name_xticks = ['','','','','','','','','','','','','','','-30','-20','-10','0','10','20','30','40','50','60','70','80','90','100']
name_xticks2 = ['J','F','M','A','M','J','J','A','S','O','N','D']


# Fig. 5 - Time series of total Arctic sea-ice area
xrangea = np.arange(11, 282, 10)
fig,ax = plt.subplots(2,1,figsize=(18,12))
fig.subplots_adjust(left=0.08,bottom=0.09,right=0.95,top=0.95,wspace=None,hspace=0.3)

# March area
#ax[0].plot(np.arange(nyears) + 1,area_mon[2,:],'-',color='red',label='CMIP6 historical')
ax[0].plot(np.arange(nyears_D012) + 1 + 180,area_mon_D012[2,:],'-',color='red',label='ATL1+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D015[2,:],'-',color='gray',label='ATL2+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D018[2,:],'-',color='blue',label='ATL3+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D000) + 1 + 64,area_mon_D000[2,:],'-',color='black',label='CTRL',linewidth=3)
ax[0].plot(np.arange(nyears_D012) + 1 + 180,area_mon_D021[2,:],'--',color='green',label='PAC1+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D022[2,:],'--',color='brown',label='PAC2+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D023[2,:],'--',color='orange',label='PAC3+3$^\circ$C',linewidth=2)
#ax[0].plot(np.arange(nyears_obs) + 1 + 29,area_obs_mon[2,:],'-',color='black',label='Observations',linewidth=2)
ax[0].legend(loc='upper left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[0].set_ylabel('March sea-ice area (10$^6$ km$^2$)',fontsize=24)
ax[0].set_xticks(xrangea)
ax[0].set_xticklabels(name_xticks)
ax[0].set_yticks(np.arange(10, 16.1, 1))
ax[0].tick_params(axis='both',labelsize=20)
ax[0].axis([150, 281, 10, 16])
ax[0].grid(linestyle='--')
ax[0].set_title('a',loc='left',fontsize=25,fontweight='bold')
    
# September area
#ax[1].plot(np.arange(nyears) + 1,area_mon[8,:],'-',color='red',label='CMIP6 historical')
ax[1].plot(np.arange(nyears_D012) + 1 + 180,area_mon_D012[8,:],'-',color='red',label='ATL1+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D015[8,:],'-',color='gray',label='ATL2+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D018[8,:],'-',color='blue',label='ATL3+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D000) + 1 + 64,area_mon_D000[8,:],'-',color='black',label='CTRL',linewidth=3)
ax[1].plot(np.arange(nyears_D012) + 1 + 180,area_mon_D021[8,:],'--',color='green',label='PAC1+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D022[8,:],'--',color='brown',label='PAC2+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,area_mon_D023[8,:],'--',color='orange',label='PAC3+3$^\circ$C',linewidth=2)
#ax[1].plot(np.arange(nyears_obs) + 1 + 29,area_obs_mon[8,:],'-',color='black',label='Observations',linewidth=2)
ax[1].legend(loc='upper left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[1].set_xlabel('Year',fontsize=24)
ax[1].set_ylabel('September sea-ice area (10$^6$ km$^2$)',fontsize=24)
ax[1].set_xticks(xrangea)
ax[1].set_xticklabels(name_xticks)
ax[1].set_yticks(np.arange(0, 8.1, 1))
ax[1].tick_params(axis='both',labelsize=20)
ax[1].axis([150, 281, 0, 8])
ax[1].grid(linestyle='--')
ax[1].set_title('b',loc='left',fontsize=25,fontweight='bold')

if save_fig == True:
    fig.savefig(dir_fig + 'fig5.png')
    fig.savefig(dir_fig + 'fig5.eps',dpi=300)


# Supp. Fig. 5b - Time series of total Arctic sea-ice volume
fig,ax = plt.subplots(2,1,figsize=(18,12))
fig.subplots_adjust(left=0.08,bottom=0.09,right=0.95,top=0.95,wspace=None,hspace=0.2)

# March volume
#ax[0].plot(np.arange(nyears) + 1,volume_mon[2,:],'-',color='red',label='CMIP6 historical')
ax[0].plot(np.arange(nyears_D012) + 1 + 180,volume_mon_D012[2,:],'--',color='red',label='ATL1+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D015[2,:],'--',color='blue',label='ATL2+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D018[2,:],'--',color='green',label='ATL3+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D000) + 1 + 64,volume_mon_D000[2,:],'-',color='blue',label='CTRL',linewidth=2)
ax[0].plot(np.arange(nyears_D012) + 1 + 180,volume_mon_D021[2,:],'--',color='cyan',label='PAC1+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D022[2,:],'--',color='orange',label='PAC2+3$^\circ$C',linewidth=2)
ax[0].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D023[2,:],'--',color='gray',label='PAC3+3$^\circ$C',linewidth=2)
#ax[0].plot(np.arange(nyears_piomas) + 1 + 29,volume_piomas_mon[2,:],'-',color='black',label='Reanalysis',linewidth=2)
ax[0].legend(loc='upper left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[0].set_ylabel('March sea-ice volume (10$^3$ km$^3$)',fontsize=22)
ax[0].set_xticks(xrangea)
ax[0].set_xticklabels(name_xticks)
ax[0].set_yticks(np.arange(10, 35.1, 5))
ax[0].tick_params(axis='both',labelsize=20)
ax[0].axis([150, 281, 10, 35])
ax[0].grid(linestyle='--')
ax[0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# September volume
#ax[1].plot(np.arange(nyears) + 1,volume_mon[8,:],'-',color='red',label='CMIP6 historical')
ax[1].plot(np.arange(nyears_D012) + 1 + 180,volume_mon_D012[8,:],'--',color='red',label='ATL1+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D015[8,:],'--',color='blue',label='ATL2+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D018[8,:],'--',color='green',label='ATL3+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D000) + 1 + 64,volume_mon_D000[8,:],'-',color='blue',label='CTRL',linewidth=2)
ax[1].plot(np.arange(nyears_D012) + 1 + 180,volume_mon_D021[8,:],'--',color='cyan',label='PAC1+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D022[8,:],'--',color='orange',label='PAC2+3$^\circ$C',linewidth=2)
ax[1].plot(np.arange(nyears_D013) + 1 + 180,volume_mon_D023[8,:],'--',color='gray',label='PAC3+3$^\circ$C',linewidth=2)
#ax[1].plot(np.arange(nyears_piomas) + 1 + 29,volume_piomas_mon[8,:],'-',color='black',label='Reanalysis',linewidth=2)
ax[1].legend(loc='upper left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[1].set_xlabel('Year',fontsize=20)
ax[1].set_ylabel('September sea-ice volume (10$^3$ km$^3$)',fontsize=22)
ax[1].set_xticks(xrangea)
ax[1].set_xticklabels(name_xticks)
ax[1].set_yticks(np.arange(0, 25.1, 5))
ax[1].tick_params(axis='both',labelsize=20)
ax[1].axis([150, 281, 0, 25])
ax[1].grid(linestyle='--')
ax[1].set_title('b',loc='left',fontsize=25,fontweight='bold')

#if save_fig == True:
#    fig.savefig(dir_fig + 'fig5b.png')
    

# Fig. 6 - Seasonal cycles of total Arctic sea-ice area and volume
fig,ax = plt.subplots(2,2,figsize=(24,18))
fig.subplots_adjust(left=0.08,bottom=0.1,right=0.95,top=0.95,wspace=0.15,hspace=0.3)    

# SIA Atlantic SST experiments
ax[0,0].set_title('Atlantic SST+3$^\circ$C experiments',fontsize=32)
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='black',linewidth=4,label='CTRL')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D012,'.-',color='red',linewidth=3,label='ATL1+3$^\circ$C ('+str(np.round(area_ym_D012[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D012[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D015,'.-',color='gray',linewidth=3,label='ATL2+3$^\circ$C ('+str(np.round(area_ym_D015[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D015[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D018,'.-',color='blue',linewidth=3,label='ATL3+3$^\circ$C ('+str(np.round(area_ym_D018[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D018[8]-area_ym_D000[8],2))+')')
ax[0,0].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,0].set_ylabel('Arctic sea-ice area (10$^6$ km$^2$)',fontsize=30)
ax[0,0].set_xticks(np.arange(nmy)+1)
ax[0,0].set_xticklabels(name_xticks2)
ax[0,0].set_yticks(np.arange(0, 14.1, 2))
ax[0,0].tick_params(axis='both',labelsize=26)
ax[0,0].axis([0, 13, 0, 14])
ax[0,0].grid(linestyle='--')
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')

# SIA Pacific SST experiments
ax[0,1].set_title('Pacific SST+3$^\circ$C experiments',fontsize=32)
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='black',linewidth=4,label='CTRL')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D021,'.-',color='green',linewidth=3,label='PAC1+3$^\circ$C ('+str(np.round(area_ym_D021[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D021[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D022,'.-',color='brown',linewidth=3,label='PAC2+3$^\circ$C ('+str(np.round(area_ym_D022[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D022[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D023,'.-',color='orange',linewidth=3,label='PAC3+3$^\circ$C ('+str(np.round(area_ym_D023[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D023[8]-area_ym_D000[8],2))+')')
ax[0,1].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,1].set_xticks(np.arange(nmy)+1)
ax[0,1].set_xticklabels(name_xticks2)
ax[0,1].set_yticks(np.arange(0, 14.1, 2))
ax[0,1].tick_params(axis='both',labelsize=26)
ax[0,1].axis([0, 13, 0, 14])
ax[0,1].grid(linestyle='--')
ax[0,1].set_title('b',loc='left',fontsize=32,fontweight='bold')

# SIV Atlantic SST experiments
ax[1,0].set_title('Atlantic SST+3$^\circ$C experiments',fontsize=32)
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='black',linewidth=4,label='CTRL')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D012,'.-',color='red',linewidth=3,label='ATL1+3$^\circ$C ('+str(np.round(volume_ym_D012[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D012[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D015,'.-',color='gray',linewidth=3,label='ATL2+3$^\circ$C ('+str(np.round(volume_ym_D015[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D015[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D018,'.-',color='blue',linewidth=3,label='ATL3+3$^\circ$C ('+str(np.round(volume_ym_D018[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D018[8]-volume_ym_D000[8],2))+')')
ax[1,0].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[1,0].set_ylabel('Arctic sea-ice volume (10$^3$ km$^3$)',fontsize=30)
ax[1,0].set_xlabel('Month',fontsize=30)
ax[1,0].set_xticks(np.arange(nmy)+1)
ax[1,0].set_xticklabels(name_xticks2)
ax[1,0].set_yticks(np.arange(0, 30.1, 5))
ax[1,0].tick_params(axis='both',labelsize=26)
ax[1,0].axis([0, 13, 0, 30])
ax[1,0].grid(linestyle='--')
ax[1,0].set_title('c',loc='left',fontsize=32,fontweight='bold')

# SIV Pacific SST experiments
ax[1,1].set_title('Pacific SST+3$^\circ$C experiments',fontsize=32)
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='black',linewidth=4,label='CTRL')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D021,'.-',color='green',linewidth=3,label='PAC1+3$^\circ$C ('+str(np.round(volume_ym_D021[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D021[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D022,'.-',color='brown',linewidth=3,label='PAC2+3$^\circ$C ('+str(np.round(volume_ym_D022[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D022[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D023,'.-',color='orange',linewidth=3,label='PAC3+3$^\circ$C ('+str(np.round(volume_ym_D023[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D023[8]-volume_ym_D000[8],2))+')')
ax[1,1].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[1,1].set_xlabel('Month',fontsize=28)
ax[1,1].set_xticks(np.arange(nmy)+1)
ax[1,1].set_xticklabels(name_xticks2)
ax[1,1].set_yticks(np.arange(0, 30.1, 5))
ax[1,1].tick_params(axis='both',labelsize=26)
ax[1,1].axis([0, 13, 0, 30])
ax[1,1].grid(linestyle='--')
ax[1,1].set_title('d',loc='left',fontsize=32,fontweight='bold')

# Save figure
if save_fig == True:
    fig.savefig(dir_fig + 'fig6.png')
    fig.savefig(dir_fig + 'fig6.eps',dpi=300)


# Fig. Supp. 6a (SST+1K) - Seasonal cycles of total Arctic sea-ice area and volume
fig,ax = plt.subplots(2,2,figsize=(24,18))
fig.subplots_adjust(left=0.08,bottom=0.1,right=0.95,top=0.95,wspace=0.15,hspace=0.3)    

# SIA Atlantic SST experiments
ax[0,0].set_title('Atlantic SST+1$^\circ$C experiments',fontsize=32)
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D013,'.--',color='purple',linewidth=2,label='ATL1+1$^\circ$C ('+str(np.round(area_ym_D013[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D013[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D016,'.--',color='red',linewidth=2,label='ATL2+1$^\circ$C ('+str(np.round(area_ym_D016[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D016[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D019,'.--',color='lightcoral',linewidth=2,label='ATL3+1$^\circ$C ('+str(np.round(area_ym_D019[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D019[8]-area_ym_D000[8],2))+')')
ax[0,0].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,0].set_ylabel('Arctic sea-ice area (10$^6$ km$^2$)',fontsize=30)
ax[0,0].set_xticks(np.arange(nmy)+1)
ax[0,0].set_xticklabels(name_xticks2)
ax[0,0].set_yticks(np.arange(0, 14.1, 2))
ax[0,0].tick_params(axis='both',labelsize=26)
ax[0,0].axis([0, 13, 0, 14])
ax[0,0].grid(linestyle='--')
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')

# SIA Pacific SST experiments
ax[0,1].set_title('Pacific SST+1$^\circ$C experiments',fontsize=32)
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D027,'.--',color='purple',linewidth=2,label='PAC1+1$^\circ$C ('+str(np.round(area_ym_D027[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D027[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D029,'.--',color='red',linewidth=2,label='PAC2+1$^\circ$C ('+str(np.round(area_ym_D029[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D029[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D024,'.--',color='lightcoral',linewidth=2,label='PAC3+1$^\circ$C ('+str(np.round(area_ym_D024[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D024[8]-area_ym_D000[8],2))+')')
ax[0,1].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,1].set_xticks(np.arange(nmy)+1)
ax[0,1].set_xticklabels(name_xticks2)
ax[0,1].set_yticks(np.arange(0, 14.1, 2))
ax[0,1].tick_params(axis='both',labelsize=26)
ax[0,1].axis([0, 13, 0, 14])
ax[0,1].grid(linestyle='--')
ax[0,1].set_title('b',loc='left',fontsize=32,fontweight='bold')

# SIV Atlantic SST experiments
ax[1,0].set_title('Atlantic SST+1$^\circ$C experiments',fontsize=32)
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D013,'.--',color='purple',linewidth=2,label='ATL1+1$^\circ$C ('+str(np.round(volume_ym_D013[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D013[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D016,'.--',color='red',linewidth=2,label='ATL2+1$^\circ$C ('+str(np.round(volume_ym_D016[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D016[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D019,'.--',color='lightcoral',linewidth=2,label='ATL3+1$^\circ$C ('+str(np.round(volume_ym_D019[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D019[8]-volume_ym_D000[8],2))+')')
ax[1,0].legend(loc='upper right',shadow=True,frameon=False,fontsize=22)
ax[1,0].set_ylabel('Arctic sea-ice volume (10$^3$ km$^3$)',fontsize=30)
ax[1,0].set_xlabel('Month',fontsize=30)
ax[1,0].set_xticks(np.arange(nmy)+1)
ax[1,0].set_xticklabels(name_xticks2)
ax[1,0].set_yticks(np.arange(0, 30.1, 5))
ax[1,0].tick_params(axis='both',labelsize=26)
ax[1,0].axis([0, 13, 0, 30])
ax[1,0].grid(linestyle='--')
ax[1,0].set_title('c',loc='left',fontsize=32,fontweight='bold')

# SIV Pacific SST experiments
ax[1,1].set_title('Pacific SST+1$^\circ$C experiments',fontsize=32)
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D027,'.--',color='purple',linewidth=2,label='PAC1+1$^\circ$C ('+str(np.round(volume_ym_D027[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D027[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D029,'.--',color='red',linewidth=2,label='PAC2+1$^\circ$C ('+str(np.round(volume_ym_D029[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D029[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D024,'.--',color='lightcoral',linewidth=2,label='PAC3+1$^\circ$C ('+str(np.round(volume_ym_D024[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D024[8]-volume_ym_D000[8],2))+')')
ax[1,1].legend(loc='upper right',shadow=True,frameon=False,fontsize=22)
ax[1,1].set_xlabel('Month',fontsize=28)
ax[1,1].set_xticks(np.arange(nmy)+1)
ax[1,1].set_xticklabels(name_xticks2)
ax[1,1].set_yticks(np.arange(0, 30.1, 5))
ax[1,1].tick_params(axis='both',labelsize=26)
ax[1,1].axis([0, 13, 0, 30])
ax[1,1].grid(linestyle='--')
ax[1,1].set_title('d',loc='left',fontsize=32,fontweight='bold')

# Save figure
if save_fig == True:
    fig.savefig(dir_fig + 'fig6a.png')


# Fig. Supp. 6b (SST+5K) - Seasonal cycles of total Arctic sea-ice area and volume
fig,ax = plt.subplots(2,2,figsize=(24,18))
fig.subplots_adjust(left=0.08,bottom=0.1,right=0.95,top=0.95,wspace=0.15,hspace=0.3)    

# SIA Atlantic SST experiments
ax[0,0].set_title('Atlantic SST+5$^\circ$C experiments',fontsize=32)
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D014,'.--',color='purple',linewidth=2,label='ATL1+5$^\circ$C ('+str(np.round(area_ym_D014[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D014[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D017,'.--',color='red',linewidth=2,label='ATL2+5$^\circ$C ('+str(np.round(area_ym_D017[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D017[8]-area_ym_D000[8],2))+')')
ax[0,0].plot(np.arange(nmy) + 1,area_ym_D020,'.--',color='lightcoral',linewidth=2,label='ATL3+5$^\circ$C ('+str(np.round(area_ym_D020[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D020[8]-area_ym_D000[8],2))+')')
ax[0,0].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,0].set_ylabel('Arctic sea-ice area (10$^6$ km$^2$)',fontsize=30)
ax[0,0].set_xticks(np.arange(nmy)+1)
ax[0,0].set_xticklabels(name_xticks2)
ax[0,0].set_yticks(np.arange(0, 14.1, 2))
ax[0,0].tick_params(axis='both',labelsize=26)
ax[0,0].axis([0, 13, 0, 14])
ax[0,0].grid(linestyle='--')
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')

# SIA Pacific SST experiments
ax[0,1].set_title('Pacific SST+5$^\circ$C experiments',fontsize=32)
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D028,'.--',color='purple',linewidth=2,label='PAC1+5$^\circ$C ('+str(np.round(area_ym_D028[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D028[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D030,'.--',color='red',linewidth=2,label='PAC2+5$^\circ$C ('+str(np.round(area_ym_D030[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D030[8]-area_ym_D000[8],2))+')')
ax[0,1].plot(np.arange(nmy) + 1,area_ym_D025,'.--',color='lightcoral',linewidth=2,label='PAC3+5$^\circ$C ('+str(np.round(area_ym_D025[2]-area_ym_D000[2],2))+'; '+str(np.round(area_ym_D025[8]-area_ym_D000[8],2))+')')
ax[0,1].legend(loc='lower left',shadow=True,frameon=False,fontsize=22)
ax[0,1].set_xticks(np.arange(nmy)+1)
ax[0,1].set_xticklabels(name_xticks2)
ax[0,1].set_yticks(np.arange(0, 14.1, 2))
ax[0,1].tick_params(axis='both',labelsize=26)
ax[0,1].axis([0, 13, 0, 14])
ax[0,1].grid(linestyle='--')
ax[0,1].set_title('b',loc='left',fontsize=32,fontweight='bold')

# SIV Atlantic SST experiments
ax[1,0].set_title('Atlantic SST+5$^\circ$C experiments',fontsize=32)
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D014,'.--',color='purple',linewidth=2,label='ATL1+5$^\circ$C ('+str(np.round(volume_ym_D014[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D014[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D017,'.--',color='red',linewidth=2,label='ATL2+5$^\circ$C ('+str(np.round(volume_ym_D017[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D017[8]-volume_ym_D000[8],2))+')')
ax[1,0].plot(np.arange(nmy) + 1,volume_ym_D020,'.--',color='lightcoral',linewidth=2,label='ATL3+5$^\circ$C ('+str(np.round(volume_ym_D020[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D020[8]-volume_ym_D000[8],2))+')')
ax[1,0].legend(loc='upper right',shadow=True,frameon=False,fontsize=22)
ax[1,0].set_ylabel('Arctic sea-ice volume (10$^3$ km$^3$)',fontsize=30)
ax[1,0].set_xlabel('Month',fontsize=30)
ax[1,0].set_xticks(np.arange(nmy)+1)
ax[1,0].set_xticklabels(name_xticks2)
ax[1,0].set_yticks(np.arange(0, 30.1, 5))
ax[1,0].tick_params(axis='both',labelsize=26)
ax[1,0].axis([0, 13, 0, 30])
ax[1,0].grid(linestyle='--')
ax[1,0].set_title('c',loc='left',fontsize=32,fontweight='bold')

# SIV Pacific SST experiments
ax[1,1].set_title('Pacific SST+5$^\circ$C experiments',fontsize=32)
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D000,'o-',color='blue',linewidth=2,label='CTRL')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D028,'.--',color='purple',linewidth=2,label='PAC1+5$^\circ$C ('+str(np.round(volume_ym_D028[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D028[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D030,'.--',color='red',linewidth=2,label='PAC2+5$^\circ$C ('+str(np.round(volume_ym_D030[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D030[8]-volume_ym_D000[8],2))+')')
ax[1,1].plot(np.arange(nmy) + 1,volume_ym_D025,'.--',color='lightcoral',linewidth=2,label='PAC3+5$^\circ$C ('+str(np.round(volume_ym_D025[2]-volume_ym_D000[2],2))+'; '+str(np.round(volume_ym_D025[8]-volume_ym_D000[8],2))+')')
ax[1,1].legend(loc='upper right',shadow=True,frameon=False,fontsize=22)
ax[1,1].set_xlabel('Month',fontsize=28)
ax[1,1].set_xticks(np.arange(nmy)+1)
ax[1,1].set_xticklabels(name_xticks2)
ax[1,1].set_yticks(np.arange(0, 30.1, 5))
ax[1,1].tick_params(axis='both',labelsize=26)
ax[1,1].axis([0, 13, 0, 30])
ax[1,1].grid(linestyle='--')
ax[1,1].set_title('d',loc='left',fontsize=32,fontweight='bold')

# Save figure
#if save_fig == True:
#    fig.savefig(dir_fig + 'fig6b.png')
