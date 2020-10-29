#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 9: OHT-sea ice scatter plots
    Fig. 16: SHF-OHT scatter plots
    Sea-ice files saved with fig5-6.py
    OHT files saved with fig2.py
    SHF_atl and SHF_pac files saved with plot_nemo_surface_fluxes.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    29/10/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats.stats import pearsonr

# Option
save_fig = True

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_fig = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Function to compute slope + SD of slope
def compute_slope(nm,x,y):
    idx = np.isfinite(x) & np.isfinite(y)
    a,b = np.polyfit(x[idx],y[idx],1)
    lreg = a*x[idx]+b
    n_freedom = nm-2
    s_yx = np.sum((y[idx]-lreg)**2)/n_freedom
    SS_xx = np.sum((x[idx]-np.mean(x[idx]))**2)
    sd_a = np.sqrt(s_yx/SS_xx)
    t_a = a / sd_a
    alpha = 0.05
    t_crit = stats.t.ppf(1-alpha/2,n_freedom)
    sig_a = 0
    if np.abs(t_a) > t_crit:
        sig_a = 1
    return a,b,sd_a,sig_a

# Load annual mean total OHT
oht_annmean_total_D000,oht_annmean_total_D012,oht_annmean_total_D013,oht_annmean_total_D014,\
    oht_annmean_total_D015,oht_annmean_total_D016,oht_annmean_total_D017,oht_annmean_total_D018,\
    oht_annmean_total_D019,oht_annmean_total_D020,oht_annmean_total_D021,oht_annmean_total_D022,\
    oht_annmean_total_D023,oht_annmean_total_D024,oht_annmean_total_D025,oht_annmean_total_D027,\
    oht_annmean_total_D028,oht_annmean_total_D029,oht_annmean_total_D030 = np.load(dir_fig\
    + 'OHT_annmean.npy')

# Load annual mean BSO OHT
oht_annmean_bso_D000,oht_annmean_bso_D012,oht_annmean_bso_D013,oht_annmean_bso_D014,\
    oht_annmean_bso_D015,oht_annmean_bso_D016,oht_annmean_bso_D017,oht_annmean_bso_D018,\
    oht_annmean_bso_D019,oht_annmean_bso_D020,oht_annmean_bso_D021,oht_annmean_bso_D022,\
    oht_annmean_bso_D023,oht_annmean_bso_D024,oht_annmean_bso_D025,oht_annmean_bso_D027,\
    oht_annmean_bso_D028,oht_annmean_bso_D029,oht_annmean_bso_D030 = np.load(dir_fig\
    + 'OHT_annmean_BSO.npy')

# Load annual mean Bering Strait OHT
oht_annmean_bering_D000,oht_annmean_bering_D012,oht_annmean_bering_D013,oht_annmean_bering_D014,\
    oht_annmean_bering_D015,oht_annmean_bering_D016,oht_annmean_bering_D017,oht_annmean_bering_D018,\
    oht_annmean_bering_D019,oht_annmean_bering_D020,oht_annmean_bering_D021,oht_annmean_bering_D022,\
    oht_annmean_bering_D023,oht_annmean_bering_D024,oht_annmean_bering_D025,oht_annmean_bering_D027,\
    oht_annmean_bering_D028,oht_annmean_bering_D029,oht_annmean_bering_D030 = np.load(dir_fig\
    + 'OHT_annmean_Bering.npy')

# Load annual mean Fram Strait OHT
oht_annmean_fram_D000,oht_annmean_fram_D012,oht_annmean_fram_D013,oht_annmean_fram_D014,\
    oht_annmean_fram_D015,oht_annmean_fram_D016,oht_annmean_fram_D017,oht_annmean_fram_D018,\
    oht_annmean_fram_D019,oht_annmean_fram_D020,oht_annmean_fram_D021,oht_annmean_fram_D022,\
    oht_annmean_fram_D023,oht_annmean_fram_D024,oht_annmean_fram_D025,oht_annmean_fram_D027,\
    oht_annmean_fram_D028,oht_annmean_fram_D029,oht_annmean_fram_D030 = np.load(dir_fig\
    + 'OHT_annmean_Fram.npy')

# Load annual mean Davis Strait OHT
oht_annmean_davis_D000,oht_annmean_davis_D012,oht_annmean_davis_D013,oht_annmean_davis_D014,\
    oht_annmean_davis_D015,oht_annmean_davis_D016,oht_annmean_davis_D017,oht_annmean_davis_D018,\
    oht_annmean_davis_D019,oht_annmean_davis_D020,oht_annmean_davis_D021,oht_annmean_davis_D022,\
    oht_annmean_davis_D023,oht_annmean_davis_D024,oht_annmean_davis_D025,oht_annmean_davis_D027,\
    oht_annmean_davis_D028,oht_annmean_davis_D029,oht_annmean_davis_D030 = np.load(dir_fig\
    + 'OHT_annmean_Davis.npy')

# Compute annual mean Atlantic OHT (BSO + Fram)
oht_annmean_atl_D000 = oht_annmean_bso_D000 + oht_annmean_fram_D000
oht_annmean_atl_D012 = oht_annmean_bso_D012 + oht_annmean_fram_D012
oht_annmean_atl_D013 = oht_annmean_bso_D013 + oht_annmean_fram_D013
oht_annmean_atl_D014 = oht_annmean_bso_D014 + oht_annmean_fram_D014
oht_annmean_atl_D015 = oht_annmean_bso_D015 + oht_annmean_fram_D015
oht_annmean_atl_D016 = oht_annmean_bso_D016 + oht_annmean_fram_D016
oht_annmean_atl_D017 = oht_annmean_bso_D017 + oht_annmean_fram_D017
oht_annmean_atl_D018 = oht_annmean_bso_D018 + oht_annmean_fram_D018
oht_annmean_atl_D019 = oht_annmean_bso_D019 + oht_annmean_fram_D019
oht_annmean_atl_D020 = oht_annmean_bso_D020 + oht_annmean_fram_D020
oht_annmean_atl_D021 = oht_annmean_bso_D021 + oht_annmean_fram_D021
oht_annmean_atl_D022 = oht_annmean_bso_D022 + oht_annmean_fram_D022
oht_annmean_atl_D023 = oht_annmean_bso_D023 + oht_annmean_fram_D023
oht_annmean_atl_D024 = oht_annmean_bso_D024 + oht_annmean_fram_D024
oht_annmean_atl_D025 = oht_annmean_bso_D025 + oht_annmean_fram_D025
oht_annmean_atl_D027 = oht_annmean_bso_D027 + oht_annmean_fram_D027
oht_annmean_atl_D028 = oht_annmean_bso_D028 + oht_annmean_fram_D028
oht_annmean_atl_D029 = oht_annmean_bso_D029 + oht_annmean_fram_D029
oht_annmean_atl_D030 = oht_annmean_bso_D030 + oht_annmean_fram_D030

# Load March Arctic SIA
area_march_D000,area_march_D012,area_march_D013,area_march_D014,area_march_D015,area_march_D016,\
    area_march_D017,area_march_D018,area_march_D019,area_march_D020,area_march_D021,area_march_D022,\
    area_march_D023,area_march_D024,area_march_D025,area_march_D027,area_march_D028,area_march_D029,\
    area_march_D030 = np.load(dir_fig + 'SIA_March.npy',allow_pickle=True)

# Load September Arctic SIA
area_sept_D000,area_sept_D012,area_sept_D013,area_sept_D014,area_sept_D015,area_sept_D016,\
    area_sept_D017,area_sept_D018,area_sept_D019,area_sept_D020,area_sept_D021,area_sept_D022,\
    area_sept_D023,area_sept_D024,area_sept_D025,area_sept_D027,area_sept_D028,area_sept_D029,\
    area_sept_D030 = np.load(dir_fig + 'SIA_September.npy',allow_pickle=True)

# Load March Arctic SIV
volume_march_D000,volume_march_D012,volume_march_D013,volume_march_D014,volume_march_D015,volume_march_D016,\
    volume_march_D017,volume_march_D018,volume_march_D019,volume_march_D020,volume_march_D021,volume_march_D022,\
    volume_march_D023,volume_march_D024,volume_march_D025,volume_march_D027,volume_march_D028,volume_march_D029,\
    volume_march_D030 = np.load(dir_fig + 'SIV_March.npy',allow_pickle=True)

# Load September Arctic SIV
volume_sept_D000,volume_sept_D012,volume_sept_D013,volume_sept_D014,volume_sept_D015,volume_sept_D016,\
    volume_sept_D017,volume_sept_D018,volume_sept_D019,volume_sept_D020,volume_sept_D021,volume_sept_D022,\
    volume_sept_D023,volume_sept_D024,volume_sept_D025,volume_sept_D027,volume_sept_D028,volume_sept_D029,\
    volume_sept_D030 = np.load(dir_fig + 'SIV_September.npy',allow_pickle=True)

# Compute mean total OHT
oht_mean_total_D000 = np.nanmean(oht_annmean_total_D000)
print('Mean total OHT CTRL (TW):',oht_mean_total_D000)
oht_mean_total_D012 = np.nanmean(oht_annmean_total_D012)
oht_mean_total_D013 = np.nanmean(oht_annmean_total_D013)
oht_mean_total_D014 = np.nanmean(oht_annmean_total_D014)
oht_mean_total_D015 = np.nanmean(oht_annmean_total_D015)
oht_mean_total_D016 = np.nanmean(oht_annmean_total_D016)
oht_mean_total_D017 = np.nanmean(oht_annmean_total_D017)
oht_mean_total_D018 = np.nanmean(oht_annmean_total_D018)
oht_mean_total_D019 = np.nanmean(oht_annmean_total_D019)
oht_mean_total_D020 = np.nanmean(oht_annmean_total_D020)
oht_mean_total_D021 = np.nanmean(oht_annmean_total_D021)
oht_mean_total_D022 = np.nanmean(oht_annmean_total_D022)
oht_mean_total_D023 = np.nanmean(oht_annmean_total_D023)
oht_mean_total_D024 = np.nanmean(oht_annmean_total_D024)
oht_mean_total_D025 = np.nanmean(oht_annmean_total_D025)
oht_mean_total_D027 = np.nanmean(oht_annmean_total_D027)
oht_mean_total_D028 = np.nanmean(oht_annmean_total_D028)
oht_mean_total_D029 = np.nanmean(oht_annmean_total_D029)
oht_mean_total_D030 = np.nanmean(oht_annmean_total_D030)
oht_mean_total_array = np.array([oht_mean_total_D012-oht_mean_total_D000,oht_mean_total_D013-oht_mean_total_D000,\
                        oht_mean_total_D014-oht_mean_total_D000,oht_mean_total_D015-oht_mean_total_D000,\
                        oht_mean_total_D016-oht_mean_total_D000,oht_mean_total_D017-oht_mean_total_D000,\
                        oht_mean_total_D018-oht_mean_total_D000,oht_mean_total_D019-oht_mean_total_D000,\
                        oht_mean_total_D020-oht_mean_total_D000,oht_mean_total_D021-oht_mean_total_D000,\
                        oht_mean_total_D022-oht_mean_total_D000,oht_mean_total_D023-oht_mean_total_D000,\
                        oht_mean_total_D024-oht_mean_total_D000,oht_mean_total_D025-oht_mean_total_D000,\
                        oht_mean_total_D027-oht_mean_total_D000,oht_mean_total_D028-oht_mean_total_D000,\
                        oht_mean_total_D029-oht_mean_total_D000,oht_mean_total_D030-oht_mean_total_D000])

# Compute mean Atlantic OHT (BSO + Fram)
oht_mean_atl_D000 = np.nanmean(oht_annmean_atl_D000)
print('Mean Atlantic OHT CTRL (TW):',oht_mean_atl_D000)
oht_mean_atl_D012 = np.nanmean(oht_annmean_atl_D012)
oht_mean_atl_D013 = np.nanmean(oht_annmean_atl_D013)
oht_mean_atl_D014 = np.nanmean(oht_annmean_atl_D014)
oht_mean_atl_D015 = np.nanmean(oht_annmean_atl_D015)
oht_mean_atl_D016 = np.nanmean(oht_annmean_atl_D016)
oht_mean_atl_D017 = np.nanmean(oht_annmean_atl_D017)
oht_mean_atl_D018 = np.nanmean(oht_annmean_atl_D018)
oht_mean_atl_D019 = np.nanmean(oht_annmean_atl_D019)
oht_mean_atl_D020 = np.nanmean(oht_annmean_atl_D020)
oht_mean_atl_D021 = np.nanmean(oht_annmean_atl_D021)
oht_mean_atl_D022 = np.nanmean(oht_annmean_atl_D022)
oht_mean_atl_D023 = np.nanmean(oht_annmean_atl_D023)
oht_mean_atl_D024 = np.nanmean(oht_annmean_atl_D024)
oht_mean_atl_D025 = np.nanmean(oht_annmean_atl_D025)
oht_mean_atl_D027 = np.nanmean(oht_annmean_atl_D027)
oht_mean_atl_D028 = np.nanmean(oht_annmean_atl_D028)
oht_mean_atl_D029 = np.nanmean(oht_annmean_atl_D029)
oht_mean_atl_D030 = np.nanmean(oht_annmean_atl_D030)
oht_mean_atl_array = np.array([oht_mean_atl_D012-oht_mean_atl_D000,oht_mean_atl_D013-oht_mean_atl_D000,\
                        oht_mean_atl_D014-oht_mean_atl_D000,oht_mean_atl_D015-oht_mean_atl_D000,\
                        oht_mean_atl_D016-oht_mean_atl_D000,oht_mean_atl_D017-oht_mean_atl_D000,\
                        oht_mean_atl_D018-oht_mean_atl_D000,oht_mean_atl_D019-oht_mean_atl_D000,\
                        oht_mean_atl_D020-oht_mean_atl_D000,oht_mean_atl_D021-oht_mean_atl_D000,\
                        oht_mean_atl_D022-oht_mean_atl_D000,oht_mean_atl_D023-oht_mean_atl_D000,\
                        oht_mean_atl_D024-oht_mean_atl_D000,oht_mean_atl_D025-oht_mean_atl_D000,\
                        oht_mean_atl_D027-oht_mean_atl_D000,oht_mean_atl_D028-oht_mean_atl_D000,\
                        oht_mean_atl_D029-oht_mean_atl_D000,oht_mean_atl_D030-oht_mean_atl_D000])

# Compute mean Pacific OHT (Bering)
oht_mean_bering_D000 = np.nanmean(oht_annmean_bering_D000)
print('Mean Pacific OHT CTRL (TW):',oht_mean_bering_D000)
oht_mean_bering_D012 = np.nanmean(oht_annmean_bering_D012)
oht_mean_bering_D013 = np.nanmean(oht_annmean_bering_D013)
oht_mean_bering_D014 = np.nanmean(oht_annmean_bering_D014)
oht_mean_bering_D015 = np.nanmean(oht_annmean_bering_D015)
oht_mean_bering_D016 = np.nanmean(oht_annmean_bering_D016)
oht_mean_bering_D017 = np.nanmean(oht_annmean_bering_D017)
oht_mean_bering_D018 = np.nanmean(oht_annmean_bering_D018)
oht_mean_bering_D019 = np.nanmean(oht_annmean_bering_D019)
oht_mean_bering_D020 = np.nanmean(oht_annmean_bering_D020)
oht_mean_bering_D021 = np.nanmean(oht_annmean_bering_D021)
oht_mean_bering_D022 = np.nanmean(oht_annmean_bering_D022)
oht_mean_bering_D023 = np.nanmean(oht_annmean_bering_D023)
oht_mean_bering_D024 = np.nanmean(oht_annmean_bering_D024)
oht_mean_bering_D025 = np.nanmean(oht_annmean_bering_D025)
oht_mean_bering_D027 = np.nanmean(oht_annmean_bering_D027)
oht_mean_bering_D028 = np.nanmean(oht_annmean_bering_D028)
oht_mean_bering_D029 = np.nanmean(oht_annmean_bering_D029)
oht_mean_bering_D030 = np.nanmean(oht_annmean_bering_D030)
oht_mean_bering_array = np.array([oht_mean_bering_D012-oht_mean_bering_D000,oht_mean_bering_D013-oht_mean_bering_D000,\
                        oht_mean_bering_D014-oht_mean_bering_D000,oht_mean_bering_D015-oht_mean_bering_D000,\
                        oht_mean_bering_D016-oht_mean_bering_D000,oht_mean_bering_D017-oht_mean_bering_D000,\
                        oht_mean_bering_D018-oht_mean_bering_D000,oht_mean_bering_D019-oht_mean_bering_D000,\
                        oht_mean_bering_D020-oht_mean_bering_D000,oht_mean_bering_D021-oht_mean_bering_D000,\
                        oht_mean_bering_D022-oht_mean_bering_D000,oht_mean_bering_D023-oht_mean_bering_D000,\
                        oht_mean_bering_D024-oht_mean_bering_D000,oht_mean_bering_D025-oht_mean_bering_D000,\
                        oht_mean_bering_D027-oht_mean_bering_D000,oht_mean_bering_D028-oht_mean_bering_D000,\
                        oht_mean_bering_D029-oht_mean_bering_D000,oht_mean_bering_D030-oht_mean_bering_D000])

# Compute mean March Arctic SIA
area_mean_march_D000 = np.nanmean(area_march_D000)
print('Mean March SIA CTRL (1e6km2):',area_mean_march_D000)
area_mean_march_D012 = np.nanmean(area_march_D012)
area_mean_march_D013 = np.nanmean(area_march_D013)
area_mean_march_D014 = np.nanmean(area_march_D014)
area_mean_march_D015 = np.nanmean(area_march_D015)
area_mean_march_D016 = np.nanmean(area_march_D016)
area_mean_march_D017 = np.nanmean(area_march_D017)
area_mean_march_D018 = np.nanmean(area_march_D018)
area_mean_march_D019 = np.nanmean(area_march_D019)
area_mean_march_D020 = np.nanmean(area_march_D020)
area_mean_march_D021 = np.nanmean(area_march_D021)
area_mean_march_D022 = np.nanmean(area_march_D022)
area_mean_march_D023 = np.nanmean(area_march_D023)
area_mean_march_D024 = np.nanmean(area_march_D024)
area_mean_march_D025 = np.nanmean(area_march_D025)
area_mean_march_D027 = np.nanmean(area_march_D027)
area_mean_march_D028 = np.nanmean(area_march_D028)
area_mean_march_D029 = np.nanmean(area_march_D029)
area_mean_march_D030 = np.nanmean(area_march_D030)
area_mean_march_array = np.array([area_mean_march_D012-area_mean_march_D000,area_mean_march_D013-area_mean_march_D000,\
                        area_mean_march_D014-area_mean_march_D000,area_mean_march_D015-area_mean_march_D000,\
                        area_mean_march_D016-area_mean_march_D000,area_mean_march_D017-area_mean_march_D000,\
                        area_mean_march_D018-area_mean_march_D000,area_mean_march_D019-area_mean_march_D000,\
                        area_mean_march_D020-area_mean_march_D000,area_mean_march_D021-area_mean_march_D000,\
                        area_mean_march_D022-area_mean_march_D000,area_mean_march_D023-area_mean_march_D000,\
                        area_mean_march_D024-area_mean_march_D000,area_mean_march_D025-area_mean_march_D000,\
                        area_mean_march_D027-area_mean_march_D000,area_mean_march_D028-area_mean_march_D000,\
                        area_mean_march_D029-area_mean_march_D000,area_mean_march_D030-area_mean_march_D000])

# Compute mean September Arctic SIA
area_mean_sept_D000 = np.nanmean(area_sept_D000)
print('Mean September SIA CTRL (1e6 km2):',area_mean_sept_D000)
area_mean_sept_D012 = np.nanmean(area_sept_D012)
area_mean_sept_D013 = np.nanmean(area_sept_D013)
area_mean_sept_D014 = np.nanmean(area_sept_D014)
area_mean_sept_D015 = np.nanmean(area_sept_D015)
area_mean_sept_D016 = np.nanmean(area_sept_D016)
area_mean_sept_D017 = np.nanmean(area_sept_D017)
area_mean_sept_D018 = np.nanmean(area_sept_D018)
area_mean_sept_D019 = np.nanmean(area_sept_D019)
area_mean_sept_D020 = np.nanmean(area_sept_D020)
area_mean_sept_D021 = np.nanmean(area_sept_D021)
area_mean_sept_D022 = np.nanmean(area_sept_D022)
area_mean_sept_D023 = np.nanmean(area_sept_D023)
area_mean_sept_D024 = np.nanmean(area_sept_D024)
area_mean_sept_D025 = np.nanmean(area_sept_D025)
area_mean_sept_D027 = np.nanmean(area_sept_D027)
area_mean_sept_D028 = np.nanmean(area_sept_D028)
area_mean_sept_D029 = np.nanmean(area_sept_D029)
area_mean_sept_D030 = np.nanmean(area_sept_D030)
area_mean_sept_array = np.array([area_mean_sept_D012-area_mean_sept_D000,area_mean_sept_D013-area_mean_sept_D000,\
                        area_mean_sept_D014-area_mean_sept_D000,area_mean_sept_D015-area_mean_sept_D000,\
                        area_mean_sept_D016-area_mean_sept_D000,area_mean_sept_D017-area_mean_sept_D000,\
                        area_mean_sept_D018-area_mean_sept_D000,area_mean_sept_D019-area_mean_sept_D000,\
                        area_mean_sept_D020-area_mean_sept_D000,area_mean_sept_D021-area_mean_sept_D000,\
                        area_mean_sept_D022-area_mean_sept_D000,area_mean_sept_D023-area_mean_sept_D000,\
                        area_mean_sept_D024-area_mean_sept_D000,area_mean_sept_D025-area_mean_sept_D000,\
                        area_mean_sept_D027-area_mean_sept_D000,area_mean_sept_D028-area_mean_sept_D000,\
                        area_mean_sept_D029-area_mean_sept_D000,area_mean_sept_D030-area_mean_sept_D000])

# Compute mean March Arctic SIV
volume_mean_march_D000 = np.nanmean(volume_march_D000)
print('Mean March SIV CTRL (1e3 km3):',volume_mean_march_D000)
volume_mean_march_D012 = np.nanmean(volume_march_D012)
volume_mean_march_D013 = np.nanmean(volume_march_D013)
volume_mean_march_D014 = np.nanmean(volume_march_D014)
volume_mean_march_D015 = np.nanmean(volume_march_D015)
volume_mean_march_D016 = np.nanmean(volume_march_D016)
volume_mean_march_D017 = np.nanmean(volume_march_D017)
volume_mean_march_D018 = np.nanmean(volume_march_D018)
volume_mean_march_D019 = np.nanmean(volume_march_D019)
volume_mean_march_D020 = np.nanmean(volume_march_D020)
volume_mean_march_D021 = np.nanmean(volume_march_D021)
volume_mean_march_D022 = np.nanmean(volume_march_D022)
volume_mean_march_D023 = np.nanmean(volume_march_D023)
volume_mean_march_D024 = np.nanmean(volume_march_D024)
volume_mean_march_D025 = np.nanmean(volume_march_D025)
volume_mean_march_D027 = np.nanmean(volume_march_D027)
volume_mean_march_D028 = np.nanmean(volume_march_D028)
volume_mean_march_D029 = np.nanmean(volume_march_D029)
volume_mean_march_D030 = np.nanmean(volume_march_D030)
volume_mean_march_array = np.array([volume_mean_march_D012-volume_mean_march_D000,volume_mean_march_D013-volume_mean_march_D000,\
                        volume_mean_march_D014-volume_mean_march_D000,volume_mean_march_D015-volume_mean_march_D000,\
                        volume_mean_march_D016-volume_mean_march_D000,volume_mean_march_D017-volume_mean_march_D000,\
                        volume_mean_march_D018-volume_mean_march_D000,volume_mean_march_D019-volume_mean_march_D000,\
                        volume_mean_march_D020-volume_mean_march_D000,volume_mean_march_D021-volume_mean_march_D000,\
                        volume_mean_march_D022-volume_mean_march_D000,volume_mean_march_D023-volume_mean_march_D000,\
                        volume_mean_march_D024-volume_mean_march_D000,volume_mean_march_D025-volume_mean_march_D000,\
                        volume_mean_march_D027-volume_mean_march_D000,volume_mean_march_D028-volume_mean_march_D000,\
                        volume_mean_march_D029-volume_mean_march_D000,volume_mean_march_D030-volume_mean_march_D000])

# Compute mean September Arctic SIV
volume_mean_sept_D000 = np.nanmean(volume_sept_D000)
print('Mean September SIV CTRL (1e3 km3):',volume_mean_sept_D000)
volume_mean_sept_D012 = np.nanmean(volume_sept_D012)
volume_mean_sept_D013 = np.nanmean(volume_sept_D013)
volume_mean_sept_D014 = np.nanmean(volume_sept_D014)
volume_mean_sept_D015 = np.nanmean(volume_sept_D015)
volume_mean_sept_D016 = np.nanmean(volume_sept_D016)
volume_mean_sept_D017 = np.nanmean(volume_sept_D017)
volume_mean_sept_D018 = np.nanmean(volume_sept_D018)
volume_mean_sept_D019 = np.nanmean(volume_sept_D019)
volume_mean_sept_D020 = np.nanmean(volume_sept_D020)
volume_mean_sept_D021 = np.nanmean(volume_sept_D021)
volume_mean_sept_D022 = np.nanmean(volume_sept_D022)
volume_mean_sept_D023 = np.nanmean(volume_sept_D023)
volume_mean_sept_D024 = np.nanmean(volume_sept_D024)
volume_mean_sept_D025 = np.nanmean(volume_sept_D025)
volume_mean_sept_D027 = np.nanmean(volume_sept_D027)
volume_mean_sept_D028 = np.nanmean(volume_sept_D028)
volume_mean_sept_D029 = np.nanmean(volume_sept_D029)
volume_mean_sept_D030 = np.nanmean(volume_sept_D030)
volume_mean_sept_array = np.array([volume_mean_sept_D012-volume_mean_sept_D000,volume_mean_sept_D013-volume_mean_sept_D000,\
                        volume_mean_sept_D014-volume_mean_sept_D000,volume_mean_sept_D015-volume_mean_sept_D000,\
                        volume_mean_sept_D016-volume_mean_sept_D000,volume_mean_sept_D017-volume_mean_sept_D000,\
                        volume_mean_sept_D018-volume_mean_sept_D000,volume_mean_sept_D019-volume_mean_sept_D000,\
                        volume_mean_sept_D020-volume_mean_sept_D000,volume_mean_sept_D021-volume_mean_sept_D000,\
                        volume_mean_sept_D022-volume_mean_sept_D000,volume_mean_sept_D023-volume_mean_sept_D000,\
                        volume_mean_sept_D024-volume_mean_sept_D000,volume_mean_sept_D025-volume_mean_sept_D000,\
                        volume_mean_sept_D027-volume_mean_sept_D000,volume_mean_sept_D028-volume_mean_sept_D000,\
                        volume_mean_sept_D029-volume_mean_sept_D000,volume_mean_sept_D030-volume_mean_sept_D000])

# Load net surface heat flux Atlantic
qtoce_atl_D000,qtoce_atl_D012,qtoce_atl_D013,qtoce_atl_D014,qtoce_atl_D015,qtoce_atl_D016,qtoce_atl_D017,qtoce_atl_D018,\
    qtoce_atl_D019,qtoce_atl_D020,qtoce_atl_D021,qtoce_atl_D022,qtoce_atl_D023,qtoce_atl_D024,qtoce_atl_D025,qtoce_atl_D027,\
    qtoce_atl_D028,qtoce_atl_D029,qtoce_atl_D030 = np.load(dir_fig + 'TotalSurfaceFluxes_atl.npy',allow_pickle=True)
qtoce_atl_array = np.array([qtoce_atl_D012-qtoce_atl_D000,qtoce_atl_D013-qtoce_atl_D000,qtoce_atl_D014-qtoce_atl_D000,qtoce_atl_D015-qtoce_atl_D000,\
                        qtoce_atl_D016-qtoce_atl_D000,qtoce_atl_D017-qtoce_atl_D000,qtoce_atl_D018-qtoce_atl_D000,qtoce_atl_D019-qtoce_atl_D000,\
                        qtoce_atl_D020-qtoce_atl_D000,qtoce_atl_D021-qtoce_atl_D000,qtoce_atl_D022-qtoce_atl_D000,qtoce_atl_D023-qtoce_atl_D000,\
                        qtoce_atl_D024-qtoce_atl_D000,qtoce_atl_D025-qtoce_atl_D000,qtoce_atl_D027-qtoce_atl_D000,qtoce_atl_D028-qtoce_atl_D000,\
                        qtoce_atl_D029-qtoce_atl_D000,qtoce_atl_D030-qtoce_atl_D000])

# Load net surface heat flux Pacific
qtoce_pac_D000,qtoce_pac_D012,qtoce_pac_D013,qtoce_pac_D014,qtoce_pac_D015,qtoce_pac_D016,qtoce_pac_D017,qtoce_pac_D018,\
    qtoce_pac_D019,qtoce_pac_D020,qtoce_pac_D021,qtoce_pac_D022,qtoce_pac_D023,qtoce_pac_D024,qtoce_pac_D025,qtoce_pac_D027,\
    qtoce_pac_D028,qtoce_pac_D029,qtoce_pac_D030 = np.load(dir_fig + 'TotalSurfaceFluxes_pac.npy',allow_pickle=True)
qtoce_pac_array = np.array([qtoce_pac_D012-qtoce_pac_D000,qtoce_pac_D013-qtoce_pac_D000,qtoce_pac_D014-qtoce_pac_D000,qtoce_pac_D015-qtoce_pac_D000,\
                        qtoce_pac_D016-qtoce_pac_D000,qtoce_pac_D017-qtoce_pac_D000,qtoce_pac_D018-qtoce_pac_D000,qtoce_pac_D019-qtoce_pac_D000,\
                        qtoce_pac_D020-qtoce_pac_D000,qtoce_pac_D021-qtoce_pac_D000,qtoce_pac_D022-qtoce_pac_D000,qtoce_pac_D023-qtoce_pac_D000,\
                        qtoce_pac_D024-qtoce_pac_D000,qtoce_pac_D025-qtoce_pac_D000,qtoce_pac_D027-qtoce_pac_D000,qtoce_pac_D028-qtoce_pac_D000,\
                        qtoce_pac_D029-qtoce_pac_D000,qtoce_pac_D030-qtoce_pac_D000])

# Compute total OHT-SIA correlations
R_arctic_march,p_arctic_march = pearsonr(oht_mean_total_array,area_mean_march_array)
R_arctic_atl_march,p_arctic_atl_march = pearsonr(oht_mean_total_array[0:9],area_mean_march_array[0:9])
R_arctic_pac_march,p_arctic_pac_march = pearsonr(oht_mean_total_array[9::],area_mean_march_array[9::])
R_arctic_sept,p_arctic_sept = pearsonr(oht_mean_total_array,area_mean_sept_array)
R_arctic_atl_sept,p_arctic_atl_sept = pearsonr(oht_mean_total_array[0:9],area_mean_sept_array[0:9])
R_arctic_pac_sept,p_arctic_pac_sept = pearsonr(oht_mean_total_array[9::],area_mean_sept_array[9::])

# Compute total OHT-SIV correlations
R_arctic_march_siv,p_arctic_march_siv = pearsonr(oht_mean_total_array,volume_mean_march_array)
R_arctic_atl_march_siv,p_arctic_atl_march_siv = pearsonr(oht_mean_total_array[0:9],volume_mean_march_array[0:9])
R_arctic_pac_march_siv,p_arctic_pac_march_siv = pearsonr(oht_mean_total_array[9::],volume_mean_march_array[9::])
R_arctic_sept_siv,p_arctic_sept_siv = pearsonr(oht_mean_total_array,volume_mean_sept_array)
R_arctic_atl_sept_siv,p_arctic_atl_sept_siv = pearsonr(oht_mean_total_array[0:9],volume_mean_sept_array[0:9])
R_arctic_pac_sept_siv,p_arctic_pac_sept_siv = pearsonr(oht_mean_total_array[9::],volume_mean_sept_array[9::])

# Compute total OHT-SIA slopes
a_arctic_march,b_arctic_march,sd_arctic_march,sig_arctic_march = compute_slope(18,oht_mean_total_array,area_mean_march_array)
a_arctic_atl_march,b_arctic_atl_march,sd_arctic_atl_march,sig_arctic_atl_march = compute_slope(9,oht_mean_total_array[0:9],area_mean_march_array[0:9])
a_arctic_pac_march,b_arctic_pac_march,sd_arctic_pac_march,sig_arctic_pac_march = compute_slope(9,oht_mean_total_array[9::],area_mean_march_array[9::])
a_arctic_sept,b_arctic_sept,sd_arctic_sept,sig_arctic_sept = compute_slope(18,oht_mean_total_array,area_mean_sept_array)
a_arctic_atl_sept,b_arctic_atl_sept,sd_arctic_atl_sept,sig_arctic_atl_sept = compute_slope(9,oht_mean_total_array[0:9],area_mean_sept_array[0:9])
a_arctic_pac_sept,b_arctic_pac_sept,sd_arctic_pac_sept,sig_arctic_pac_sept = compute_slope(9,oht_mean_total_array[9::],area_mean_sept_array[9::])

# Compute total OHT-SIV slopes
a_arctic_march_siv,b_arctic_march_siv,sd_arctic_march_siv,sig_arctic_march_siv = compute_slope(18,oht_mean_total_array,volume_mean_march_array)
a_arctic_atl_march_siv,b_arctic_atl_march_siv,sd_arctic_atl_march_siv,sig_arctic_atl_march_siv = compute_slope(9,oht_mean_total_array[0:9],volume_mean_march_array[0:9])
a_arctic_pac_march_siv,b_arctic_pac_march_siv,sd_arctic_pac_march_siv,sig_arctic_pac_march_siv = compute_slope(9,oht_mean_total_array[9::],volume_mean_march_array[9::])
a_arctic_sept_siv,b_arctic_sept_siv,sd_arctic_sept_siv,sig_arctic_sept_siv = compute_slope(18,oht_mean_total_array,volume_mean_sept_array)
a_arctic_atl_sept_siv,b_arctic_atl_sept_siv,sd_arctic_atl_sept_siv,sig_arctic_atl_sept_siv = compute_slope(9,oht_mean_total_array[0:9],volume_mean_sept_array[0:9])
a_arctic_pac_sept_siv,b_arctic_pac_sept_siv,sd_arctic_pac_sept_siv,sig_arctic_pac_sept_siv = compute_slope(9,oht_mean_total_array[9::],volume_mean_sept_array[9::])

# Compute Atlantic OHT-SIA correlations
R_atl_march,p_atl_march = pearsonr(oht_mean_atl_array,area_mean_march_array)
R_atl_atl_march,p_atl_atl_march = pearsonr(oht_mean_atl_array[0:9],area_mean_march_array[0:9])
R_atl_pac_march,p_atl_pac_march = pearsonr(oht_mean_atl_array[9::],area_mean_march_array[9::])
R_atl_sept,p_atl_sept = pearsonr(oht_mean_atl_array,area_mean_sept_array)
R_atl_atl_sept,p_atl_atl_sept = pearsonr(oht_mean_atl_array[0:9],area_mean_sept_array[0:9])
R_atl_pac_sept,p_atl_pac_sept = pearsonr(oht_mean_atl_array[9::],area_mean_sept_array[9::])

# Compute Atlantic OHT-SIV correlations
R_atl_march_siv,p_atl_march_siv = pearsonr(oht_mean_atl_array,volume_mean_march_array)
R_atl_atl_march_siv,p_atl_atl_march_siv = pearsonr(oht_mean_atl_array[0:9],volume_mean_march_array[0:9])
R_atl_pac_march_siv,p_atl_pac_march_siv = pearsonr(oht_mean_atl_array[9::],volume_mean_march_array[9::])
R_atl_sept_siv,p_atl_sept_siv = pearsonr(oht_mean_atl_array,volume_mean_sept_array)
R_atl_atl_sept_siv,p_atl_atl_sept_siv = pearsonr(oht_mean_atl_array[0:9],volume_mean_sept_array[0:9])
R_atl_pac_sept_siv,p_atl_pac_sept_siv = pearsonr(oht_mean_atl_array[9::],volume_mean_sept_array[9::])

# Compute Atlantic OHT-SIA slopes
a_atl_march,b_atl_march,sd_atl_march,sig_atl_march = compute_slope(18,oht_mean_atl_array,area_mean_march_array)
a_atl_atl_march,b_atl_atl_march,sd_atl_atl_march,sig_atl_atl_march = compute_slope(9,oht_mean_atl_array[0:9],area_mean_march_array[0:9])
a_atl_pac_march,b_atl_pac_march,sd_atl_pac_march,sig_atl_pac_march = compute_slope(9,oht_mean_atl_array[9::],area_mean_march_array[9::])
a_atl_sept,b_atl_sept,sd_atl_sept,sig_atl_sept = compute_slope(18,oht_mean_atl_array,area_mean_sept_array)
a_atl_atl_sept,b_atl_atl_sept,sd_atl_atl_sept,sig_atl_atl_sept = compute_slope(9,oht_mean_atl_array[0:9],area_mean_sept_array[0:9])
a_atl_pac_sept,b_atl_pac_sept,sd_atl_pac_sept,sig_atl_pac_sept = compute_slope(9,oht_mean_atl_array[9::],area_mean_sept_array[9::])

# Compute Atlantic OHT-SIV slopes
a_atl_march_siv,b_atl_march_siv,sd_atl_march_siv,sig_atl_march_siv = compute_slope(18,oht_mean_atl_array,volume_mean_march_array)
a_atl_atl_march_siv,b_atl_atl_march_siv,sd_atl_atl_march_siv,sig_atl_atl_march_siv = compute_slope(9,oht_mean_atl_array[0:9],volume_mean_march_array[0:9])
a_atl_pac_march_siv,b_atl_pac_march_siv,sd_atl_pac_march_siv,sig_atl_pac_march_siv = compute_slope(9,oht_mean_atl_array[9::],volume_mean_march_array[9::])
a_atl_sept_siv,b_atl_sept_siv,sd_atl_sept_siv,sig_atl_sept_siv = compute_slope(18,oht_mean_atl_array,volume_mean_sept_array)
a_atl_atl_sept_siv,b_atl_atl_sept_siv,sd_atl_atl_sept_siv,sig_atl_atl_sept_siv = compute_slope(9,oht_mean_atl_array[0:9],volume_mean_sept_array[0:9])
a_atl_pac_sept_siv,b_atl_pac_sept_siv,sd_atl_pac_sept_siv,sig_atl_pac_sept_siv = compute_slope(9,oht_mean_atl_array[9::],volume_mean_sept_array[9::])

# Compute Bering OHT-SIA correlations
R_bering_march,p_bering_march = pearsonr(oht_mean_bering_array,area_mean_march_array)
R_bering_atl_march,p_bering_atl_march = pearsonr(oht_mean_bering_array[0:9],area_mean_march_array[0:9])
R_bering_pac_march,p_bering_pac_march = pearsonr(oht_mean_bering_array[9::],area_mean_march_array[9::])
R_bering_sept,p_bering_sept = pearsonr(oht_mean_bering_array,area_mean_sept_array)
R_bering_atl_sept,p_bering_atl_sept = pearsonr(oht_mean_bering_array[0:9],area_mean_sept_array[0:9])
R_bering_pac_sept,p_bering_pac_sept = pearsonr(oht_mean_bering_array[9::],area_mean_sept_array[9::])

# Compute Bering OHT-SIV correlations
R_bering_march_siv,p_bering_march_siv = pearsonr(oht_mean_bering_array,volume_mean_march_array)
R_bering_atl_march_siv,p_bering_atl_march_siv = pearsonr(oht_mean_bering_array[0:9],volume_mean_march_array[0:9])
R_bering_pac_march_siv,p_bering_pac_march_siv = pearsonr(oht_mean_bering_array[9::],volume_mean_march_array[9::])
R_bering_sept_siv,p_bering_sept_siv = pearsonr(oht_mean_bering_array,volume_mean_sept_array)
R_bering_atl_sept_siv,p_bering_atl_sept_siv = pearsonr(oht_mean_bering_array[0:9],volume_mean_sept_array[0:9])
R_bering_pac_sept_siv,p_bering_pac_sept_siv = pearsonr(oht_mean_bering_array[9::],volume_mean_sept_array[9::])

# Compute Bering OHT-SIA slopes
a_bering_march,b_bering_march,sd_bering_march,sig_bering_march = compute_slope(18,oht_mean_bering_array,area_mean_march_array)
a_bering_atl_march,b_bering_atl_march,sd_bering_atl_march,sig_bering_atl_march = compute_slope(9,oht_mean_bering_array[0:9],area_mean_march_array[0:9])
a_bering_pac_march,b_bering_pac_march,sd_bering_pac_march,sig_bering_pac_march = compute_slope(9,oht_mean_bering_array[9::],area_mean_march_array[9::])
a_bering_sept,b_bering_sept,sd_bering_sept,sig_bering_sept = compute_slope(18,oht_mean_bering_array,area_mean_sept_array)
a_bering_atl_sept,b_bering_atl_sept,sd_bering_atl_sept,sig_bering_atl_sept = compute_slope(9,oht_mean_bering_array[0:9],area_mean_sept_array[0:9])
a_bering_pac_sept,b_bering_pac_sept,sd_bering_pac_sept,sig_bering_pac_sept = compute_slope(9,oht_mean_bering_array[9::],area_mean_sept_array[9::])

# Compute Bering OHT-SIV slopes
a_bering_march_siv,b_bering_march_siv,sd_bering_march_siv,sig_bering_march_siv = compute_slope(18,oht_mean_bering_array,volume_mean_march_array)
a_bering_atl_march_siv,b_bering_atl_march_siv,sd_bering_atl_march_siv,sig_bering_atl_march_siv = compute_slope(9,oht_mean_bering_array[0:9],volume_mean_march_array[0:9])
a_bering_pac_march_siv,b_bering_pac_march_siv,sd_bering_pac_march_siv,sig_bering_pac_march_siv = compute_slope(9,oht_mean_bering_array[9::],volume_mean_march_array[9::])
a_bering_sept_siv,b_bering_sept_siv,sd_bering_sept_siv,sig_bering_sept_siv = compute_slope(18,oht_mean_bering_array,volume_mean_sept_array)
a_bering_atl_sept_siv,b_bering_atl_sept_siv,sd_bering_atl_sept_siv,sig_bering_atl_sept_siv = compute_slope(9,oht_mean_bering_array[0:9],volume_mean_sept_array[0:9])
a_bering_pac_sept_siv,b_bering_pac_sept_siv,sd_bering_pac_sept_siv,sig_bering_pac_sept_siv = compute_slope(9,oht_mean_bering_array[9::],volume_mean_sept_array[9::])

# Compute SHF_Atl-OHT_Atl correlations
R_shf_atl,p_shf_atl = pearsonr(qtoce_atl_array,oht_mean_atl_array)
R_shf_atl_atl,p_shf_atl_atl = pearsonr(qtoce_atl_array[0:9],oht_mean_atl_array[0:9])
R_shf_atl_pac,p_shf_atl_pac = pearsonr(qtoce_atl_array[9::],oht_mean_atl_array[9::])

# Compute SHF_Pac-OHT_Bering correlations
R_shf_pac,p_shf_pac = pearsonr(qtoce_pac_array,oht_mean_bering_array)
R_shf_pac_atl,p_shf_pac_atl = pearsonr(qtoce_pac_array[0:9],oht_mean_bering_array[0:9])
R_shf_pac_pac,p_shf_pac_pac = pearsonr(qtoce_pac_array[9::],oht_mean_bering_array[9::])

# Compute SHF_Atl-OHT_Atl slopes
a_shf_atl_atl,b_shfl_atl_atl,sd_shf_atl_atl,sig_shf_atl_atl = compute_slope(9,qtoce_atl_array[0:9],oht_mean_atl_array[0:9])
a_shf_atl_pac,b_shfl_atl_pac,sd_shf_atl_pac,sig_shf_atl_pac = compute_slope(9,qtoce_atl_array[9::],oht_mean_atl_array[9::])
print('SHF_Atl-OHT_Atl in ATL experiments:',a_shf_atl_atl)
print('SHF_Atl-OHT_Atl in PAC experiments:',a_shf_atl_pac)
print('PAC1+3K:',qtoce_atl_array[9],oht_mean_atl_array[9])
print('PAC2+3K:',qtoce_atl_array[10],oht_mean_atl_array[10])


# Fig. 9 - OHT-SIA/SIV Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(2,2,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.18,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Total OHT - March Arctic SIA
ax[0,0].plot(oht_mean_total_D013-oht_mean_total_D000,area_mean_march_D013-area_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D012-oht_mean_total_D000,area_mean_march_D012-area_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D014-oht_mean_total_D000,area_mean_march_D014-area_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D016-oht_mean_total_D000,area_mean_march_D016-area_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D015-oht_mean_total_D000,area_mean_march_D015-area_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D017-oht_mean_total_D000,area_mean_march_D017-area_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D019-oht_mean_total_D000,area_mean_march_D019-area_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D018-oht_mean_total_D000,area_mean_march_D018-area_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D020-oht_mean_total_D000,area_mean_march_D020-area_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D027-oht_mean_total_D000,area_mean_march_D027-area_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D021-oht_mean_total_D000,area_mean_march_D021-area_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D028-oht_mean_total_D000,area_mean_march_D028-area_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D029-oht_mean_total_D000,area_mean_march_D029-area_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D022-oht_mean_total_D000,area_mean_march_D022-area_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D030-oht_mean_total_D000,area_mean_march_D030-area_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D024-oht_mean_total_D000,area_mean_march_D024-area_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D023-oht_mean_total_D000,area_mean_march_D023-area_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D025-oht_mean_total_D000,area_mean_march_D025-area_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,0].set_ylabel('$\Delta$SIA$_{Arctic,March}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,0].set_xticks(np.arange(0, 140.1, 20))
ax[0,0].set_yticks(np.arange(-5, 0.3, 1))
ax[0,0].axis([0, 140, -5, 0.3])
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].annotate('$a$ = '+str(int(np.round(a_arctic_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Total OHT - September Arctic SIA
ax[0,1].plot(oht_mean_total_D013-oht_mean_total_D000,area_mean_sept_D013-area_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D012-oht_mean_total_D000,area_mean_sept_D012-area_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D014-oht_mean_total_D000,area_mean_sept_D014-area_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D016-oht_mean_total_D000,area_mean_sept_D016-area_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D015-oht_mean_total_D000,area_mean_sept_D015-area_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D017-oht_mean_total_D000,area_mean_sept_D017-area_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D019-oht_mean_total_D000,area_mean_sept_D019-area_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D018-oht_mean_total_D000,area_mean_sept_D018-area_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D020-oht_mean_total_D000,area_mean_sept_D020-area_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D027-oht_mean_total_D000,area_mean_sept_D027-area_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D021-oht_mean_total_D000,area_mean_sept_D021-area_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D028-oht_mean_total_D000,area_mean_sept_D028-area_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D029-oht_mean_total_D000,area_mean_sept_D029-area_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D022-oht_mean_total_D000,area_mean_sept_D022-area_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D030-oht_mean_total_D000,area_mean_sept_D030-area_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D024-oht_mean_total_D000,area_mean_sept_D024-area_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D023-oht_mean_total_D000,area_mean_sept_D023-area_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D025-oht_mean_total_D000,area_mean_sept_D025-area_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,1].set_ylabel('$\Delta$SIA$_{Arctic,September}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,1].set_xticks(np.arange(0, 140.1, 20))
ax[0,1].set_yticks(np.arange(-5, 0.3, 1))
ax[0,1].axis([0, 140, -5, 0.3])
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Total OHT - March Arctic SIV
ax[1,0].plot(oht_mean_total_D013-oht_mean_total_D000,volume_mean_march_D013-volume_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D012-oht_mean_total_D000,volume_mean_march_D012-volume_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D014-oht_mean_total_D000,volume_mean_march_D014-volume_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D016-oht_mean_total_D000,volume_mean_march_D016-volume_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D015-oht_mean_total_D000,volume_mean_march_D015-volume_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D017-oht_mean_total_D000,volume_mean_march_D017-volume_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D019-oht_mean_total_D000,volume_mean_march_D019-volume_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D018-oht_mean_total_D000,volume_mean_march_D018-volume_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D020-oht_mean_total_D000,volume_mean_march_D020-volume_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D027-oht_mean_total_D000,volume_mean_march_D027-volume_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D021-oht_mean_total_D000,volume_mean_march_D021-volume_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D028-oht_mean_total_D000,volume_mean_march_D028-volume_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D029-oht_mean_total_D000,volume_mean_march_D029-volume_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D022-oht_mean_total_D000,volume_mean_march_D022-volume_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D030-oht_mean_total_D000,volume_mean_march_D030-volume_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D024-oht_mean_total_D000,volume_mean_march_D024-volume_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D023-oht_mean_total_D000,volume_mean_march_D023-volume_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D025-oht_mean_total_D000,volume_mean_march_D025-volume_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,0].set_xlabel('$\Delta$OHT$_{Arctic}$  (TW)',fontsize=20)
ax[1,0].set_ylabel('$\Delta$SIV$_{Arctic,March}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,0].set_xticks(np.arange(0, 140.1, 20))
ax[1,0].set_yticks(np.arange(-15, 1, 3))
ax[1,0].axis([0, 140, -15, 1])
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].annotate('$a$ = '+str(int(np.round(a_arctic_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Total OHT - September Arctic SIV
ax[1,1].plot(oht_mean_total_D013-oht_mean_total_D000,volume_mean_sept_D013-volume_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D012-oht_mean_total_D000,volume_mean_sept_D012-volume_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D014-oht_mean_total_D000,volume_mean_sept_D014-volume_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D016-oht_mean_total_D000,volume_mean_sept_D016-volume_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D015-oht_mean_total_D000,volume_mean_sept_D015-volume_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D017-oht_mean_total_D000,volume_mean_sept_D017-volume_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D019-oht_mean_total_D000,volume_mean_sept_D019-volume_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D018-oht_mean_total_D000,volume_mean_sept_D018-volume_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D020-oht_mean_total_D000,volume_mean_sept_D020-volume_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D027-oht_mean_total_D000,volume_mean_sept_D027-volume_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D021-oht_mean_total_D000,volume_mean_sept_D021-volume_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D028-oht_mean_total_D000,volume_mean_sept_D028-volume_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D029-oht_mean_total_D000,volume_mean_sept_D029-volume_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D022-oht_mean_total_D000,volume_mean_sept_D022-volume_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D030-oht_mean_total_D000,volume_mean_sept_D030-volume_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D024-oht_mean_total_D000,volume_mean_sept_D024-volume_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D023-oht_mean_total_D000,volume_mean_sept_D023-volume_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D025-oht_mean_total_D000,volume_mean_sept_D025-volume_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,1].set_xlabel('$\Delta$OHT$_{Arctic}$  (TW)',fontsize=20)
ax[1,1].set_ylabel('$\Delta$SIV$_{Arctic,September}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,1].set_xticks(np.arange(0, 140.1, 20))
ax[1,1].set_yticks(np.arange(-15, 1, 3))
ax[1,1].axis([0, 140, -15, 1])
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')
ax[1,1].legend(shadow=True,frameon=False,fontsize=16,bbox_to_anchor=(1,-0.2),ncol=6)

if save_fig == True:
    fig.savefig(dir_fig + 'fig9.png')
    fig.savefig(dir_fig + 'fig9.eps',dpi=300)


# Fig. 9b - Atlantic OHT-SIA/SIV Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(2,2,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.18,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Atlantic OHT - March Arctic SIA
ax[0,0].plot(oht_mean_atl_D013-oht_mean_atl_D000,area_mean_march_D013-area_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D012-oht_mean_atl_D000,area_mean_march_D012-area_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D014-oht_mean_atl_D000,area_mean_march_D014-area_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,0].plot(oht_mean_atl_D016-oht_mean_atl_D000,area_mean_march_D016-area_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D015-oht_mean_atl_D000,area_mean_march_D015-area_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D017-oht_mean_atl_D000,area_mean_march_D017-area_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,0].plot(oht_mean_atl_D019-oht_mean_atl_D000,area_mean_march_D019-area_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D018-oht_mean_atl_D000,area_mean_march_D018-area_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D020-oht_mean_atl_D000,area_mean_march_D020-area_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,0].plot(oht_mean_atl_D027-oht_mean_atl_D000,area_mean_march_D027-area_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D021-oht_mean_atl_D000,area_mean_march_D021-area_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D028-oht_mean_atl_D000,area_mean_march_D028-area_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,0].plot(oht_mean_atl_D029-oht_mean_atl_D000,area_mean_march_D029-area_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D022-oht_mean_atl_D000,area_mean_march_D022-area_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D030-oht_mean_atl_D000,area_mean_march_D030-area_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,0].plot(oht_mean_atl_D024-oht_mean_atl_D000,area_mean_march_D024-area_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,0].plot(oht_mean_atl_D023-oht_mean_atl_D000,area_mean_march_D023-area_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,0].plot(oht_mean_atl_D025-oht_mean_atl_D000,area_mean_march_D025-area_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,0].set_ylabel('$\Delta$SIA$_{Arctic,March}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,0].set_xticks(np.arange(0, 140.1, 20))
ax[0,0].set_yticks(np.arange(-5, 0.3, 1))
ax[0,0].axis([0, 140, -5, 0.3])
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].annotate('$a$ = '+str(int(np.round(a_atl_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{atl}$ = '+str(int(np.round(a_atl_atl_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{pac}$ = '+str(int(np.round(a_atl_pac_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Atlantic OHT - September Arctic SIA
ax[0,1].plot(oht_mean_atl_D013-oht_mean_atl_D000,area_mean_sept_D013-area_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D012-oht_mean_atl_D000,area_mean_sept_D012-area_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D014-oht_mean_atl_D000,area_mean_sept_D014-area_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,1].plot(oht_mean_atl_D016-oht_mean_atl_D000,area_mean_sept_D016-area_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D015-oht_mean_atl_D000,area_mean_sept_D015-area_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D017-oht_mean_atl_D000,area_mean_sept_D017-area_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,1].plot(oht_mean_atl_D019-oht_mean_atl_D000,area_mean_sept_D019-area_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D018-oht_mean_atl_D000,area_mean_sept_D018-area_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D020-oht_mean_atl_D000,area_mean_sept_D020-area_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,1].plot(oht_mean_atl_D027-oht_mean_atl_D000,area_mean_sept_D027-area_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D021-oht_mean_atl_D000,area_mean_sept_D021-area_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D028-oht_mean_atl_D000,area_mean_sept_D028-area_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,1].plot(oht_mean_atl_D029-oht_mean_atl_D000,area_mean_sept_D029-area_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D022-oht_mean_atl_D000,area_mean_sept_D022-area_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D030-oht_mean_atl_D000,area_mean_sept_D030-area_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,1].plot(oht_mean_atl_D024-oht_mean_atl_D000,area_mean_sept_D024-area_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,1].plot(oht_mean_atl_D023-oht_mean_atl_D000,area_mean_sept_D023-area_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,1].plot(oht_mean_atl_D025-oht_mean_atl_D000,area_mean_sept_D025-area_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,1].set_ylabel('$\Delta$SIA$_{Arctic,September}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,1].set_xticks(np.arange(0, 140.1, 20))
ax[0,1].set_yticks(np.arange(-5, 0.3, 1))
ax[0,1].axis([0, 140, -5, 0.3])
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].annotate('$a$ = '+str(int(np.round(a_atl_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{atl}$ = '+str(int(np.round(a_atl_atl_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{pac}$ = '+str(int(np.round(a_atl_pac_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Atlantic OHT - March Arctic SIV
ax[1,0].plot(oht_mean_atl_D013-oht_mean_atl_D000,volume_mean_march_D013-volume_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D012-oht_mean_atl_D000,volume_mean_march_D012-volume_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D014-oht_mean_atl_D000,volume_mean_march_D014-volume_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,0].plot(oht_mean_atl_D016-oht_mean_atl_D000,volume_mean_march_D016-volume_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D015-oht_mean_atl_D000,volume_mean_march_D015-volume_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D017-oht_mean_atl_D000,volume_mean_march_D017-volume_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,0].plot(oht_mean_atl_D019-oht_mean_atl_D000,volume_mean_march_D019-volume_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D018-oht_mean_atl_D000,volume_mean_march_D018-volume_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D020-oht_mean_atl_D000,volume_mean_march_D020-volume_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,0].plot(oht_mean_atl_D027-oht_mean_atl_D000,volume_mean_march_D027-volume_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D021-oht_mean_atl_D000,volume_mean_march_D021-volume_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D028-oht_mean_atl_D000,volume_mean_march_D028-volume_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,0].plot(oht_mean_atl_D029-oht_mean_atl_D000,volume_mean_march_D029-volume_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D022-oht_mean_atl_D000,volume_mean_march_D022-volume_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D030-oht_mean_atl_D000,volume_mean_march_D030-volume_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,0].plot(oht_mean_atl_D024-oht_mean_atl_D000,volume_mean_march_D024-volume_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,0].plot(oht_mean_atl_D023-oht_mean_atl_D000,volume_mean_march_D023-volume_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,0].plot(oht_mean_atl_D025-oht_mean_atl_D000,volume_mean_march_D025-volume_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,0].set_xlabel('$\Delta$OHT$_{Atlantic}$  (TW)',fontsize=20)
ax[1,0].set_ylabel('$\Delta$SIV$_{Arctic,March}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,0].set_xticks(np.arange(0, 140.1, 20))
ax[1,0].set_yticks(np.arange(-15, 1, 3))
ax[1,0].axis([0, 140, -15, 1])
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].annotate('$a$ = '+str(int(np.round(a_atl_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{atl}$ = '+str(int(np.round(a_atl_atl_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{pac}$ = '+str(int(np.round(a_atl_pac_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Atlantic OHT - September Arctic SIV
ax[1,1].plot(oht_mean_atl_D013-oht_mean_atl_D000,volume_mean_sept_D013-volume_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D012-oht_mean_atl_D000,volume_mean_sept_D012-volume_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D014-oht_mean_atl_D000,volume_mean_sept_D014-volume_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,1].plot(oht_mean_atl_D016-oht_mean_atl_D000,volume_mean_sept_D016-volume_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D015-oht_mean_atl_D000,volume_mean_sept_D015-volume_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D017-oht_mean_atl_D000,volume_mean_sept_D017-volume_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,1].plot(oht_mean_atl_D019-oht_mean_atl_D000,volume_mean_sept_D019-volume_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D018-oht_mean_atl_D000,volume_mean_sept_D018-volume_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D020-oht_mean_atl_D000,volume_mean_sept_D020-volume_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,1].plot(oht_mean_atl_D027-oht_mean_atl_D000,volume_mean_sept_D027-volume_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D021-oht_mean_atl_D000,volume_mean_sept_D021-volume_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D028-oht_mean_atl_D000,volume_mean_sept_D028-volume_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,1].plot(oht_mean_atl_D029-oht_mean_atl_D000,volume_mean_sept_D029-volume_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D022-oht_mean_atl_D000,volume_mean_sept_D022-volume_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D030-oht_mean_atl_D000,volume_mean_sept_D030-volume_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,1].plot(oht_mean_atl_D024-oht_mean_atl_D000,volume_mean_sept_D024-volume_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,1].plot(oht_mean_atl_D023-oht_mean_atl_D000,volume_mean_sept_D023-volume_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,1].plot(oht_mean_atl_D025-oht_mean_atl_D000,volume_mean_sept_D025-volume_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,1].set_xlabel('$\Delta$OHT$_{Atlantic}$  (TW)',fontsize=20)
ax[1,1].set_ylabel('$\Delta$SIV$_{Arctic,September}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,1].set_xticks(np.arange(0, 140.1, 20))
ax[1,1].set_yticks(np.arange(-15, 1, 3))
ax[1,1].axis([0, 140, -15, 1])
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].annotate('$a$ = '+str(int(np.round(a_atl_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{atl}$ = '+str(int(np.round(a_atl_atl_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{pac}$ = '+str(int(np.round(a_atl_pac_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')
ax[1,1].legend(shadow=True,frameon=False,fontsize=16,bbox_to_anchor=(1,-0.2),ncol=6)

#if save_fig == True:
#    fig.savefig(dir_fig + 'fig9b.png')


# Fig. 9c - Bering OHT-SIA/SIV Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(2,2,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.18,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Bering OHT - March Arctic SIA
ax[0,0].plot(oht_mean_bering_D013-oht_mean_bering_D000,area_mean_march_D013-area_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D012-oht_mean_bering_D000,area_mean_march_D012-area_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D014-oht_mean_bering_D000,area_mean_march_D014-area_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,0].plot(oht_mean_bering_D016-oht_mean_bering_D000,area_mean_march_D016-area_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D015-oht_mean_bering_D000,area_mean_march_D015-area_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D017-oht_mean_bering_D000,area_mean_march_D017-area_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,0].plot(oht_mean_bering_D019-oht_mean_bering_D000,area_mean_march_D019-area_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D018-oht_mean_bering_D000,area_mean_march_D018-area_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D020-oht_mean_bering_D000,area_mean_march_D020-area_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,0].plot(oht_mean_bering_D027-oht_mean_bering_D000,area_mean_march_D027-area_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D021-oht_mean_bering_D000,area_mean_march_D021-area_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D028-oht_mean_bering_D000,area_mean_march_D028-area_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,0].plot(oht_mean_bering_D029-oht_mean_bering_D000,area_mean_march_D029-area_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D022-oht_mean_bering_D000,area_mean_march_D022-area_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D030-oht_mean_bering_D000,area_mean_march_D030-area_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,0].plot(oht_mean_bering_D024-oht_mean_bering_D000,area_mean_march_D024-area_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,0].plot(oht_mean_bering_D023-oht_mean_bering_D000,area_mean_march_D023-area_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,0].plot(oht_mean_bering_D025-oht_mean_bering_D000,area_mean_march_D025-area_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,0].set_ylabel('$\Delta$SIA$_{Arctic,March}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,0].set_xticks(np.arange(0, 30.1, 5))
ax[0,0].set_yticks(np.arange(-5, 0.3, 1))
ax[0,0].axis([0, 30, -5, 0.3])
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].annotate('$a$ = '+str(int(np.round(a_bering_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{atl}$ = '+str(int(np.round(a_bering_atl_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{pac}$ = '+str(int(np.round(a_bering_pac_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Bering OHT - September Arctic SIA
ax[0,1].plot(oht_mean_bering_D013-oht_mean_bering_D000,area_mean_sept_D013-area_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D012-oht_mean_bering_D000,area_mean_sept_D012-area_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D014-oht_mean_bering_D000,area_mean_sept_D014-area_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0,1].plot(oht_mean_bering_D016-oht_mean_bering_D000,area_mean_sept_D016-area_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D015-oht_mean_bering_D000,area_mean_sept_D015-area_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D017-oht_mean_bering_D000,area_mean_sept_D017-area_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0,1].plot(oht_mean_bering_D019-oht_mean_bering_D000,area_mean_sept_D019-area_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D018-oht_mean_bering_D000,area_mean_sept_D018-area_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D020-oht_mean_bering_D000,area_mean_sept_D020-area_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0,1].plot(oht_mean_bering_D027-oht_mean_bering_D000,area_mean_sept_D027-area_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D021-oht_mean_bering_D000,area_mean_sept_D021-area_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D028-oht_mean_bering_D000,area_mean_sept_D028-area_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0,1].plot(oht_mean_bering_D029-oht_mean_bering_D000,area_mean_sept_D029-area_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D022-oht_mean_bering_D000,area_mean_sept_D022-area_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D030-oht_mean_bering_D000,area_mean_sept_D030-area_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0,1].plot(oht_mean_bering_D024-oht_mean_bering_D000,area_mean_sept_D024-area_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0,1].plot(oht_mean_bering_D023-oht_mean_bering_D000,area_mean_sept_D023-area_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0,1].plot(oht_mean_bering_D025-oht_mean_bering_D000,area_mean_sept_D025-area_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0,1].set_ylabel('$\Delta$SIA$_{Arctic,September}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,1].set_xticks(np.arange(0, 30.1, 5))
ax[0,1].set_yticks(np.arange(-5, 0.3, 1))
ax[0,1].axis([0, 30, -5, 0.3])
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].annotate('$a$ = '+str(int(np.round(a_bering_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{atl}$ = '+str(int(np.round(a_bering_atl_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{pac}$ = '+str(int(np.round(a_bering_pac_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Bering OHT - March Arctic SIV
ax[1,0].plot(oht_mean_bering_D013-oht_mean_bering_D000,volume_mean_march_D013-volume_mean_march_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D012-oht_mean_bering_D000,volume_mean_march_D012-volume_mean_march_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D014-oht_mean_bering_D000,volume_mean_march_D014-volume_mean_march_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,0].plot(oht_mean_bering_D016-oht_mean_bering_D000,volume_mean_march_D016-volume_mean_march_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D015-oht_mean_bering_D000,volume_mean_march_D015-volume_mean_march_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D017-oht_mean_bering_D000,volume_mean_march_D017-volume_mean_march_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,0].plot(oht_mean_bering_D019-oht_mean_bering_D000,volume_mean_march_D019-volume_mean_march_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D018-oht_mean_bering_D000,volume_mean_march_D018-volume_mean_march_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D020-oht_mean_bering_D000,volume_mean_march_D020-volume_mean_march_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,0].plot(oht_mean_bering_D027-oht_mean_bering_D000,volume_mean_march_D027-volume_mean_march_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D021-oht_mean_bering_D000,volume_mean_march_D021-volume_mean_march_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D028-oht_mean_bering_D000,volume_mean_march_D028-volume_mean_march_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,0].plot(oht_mean_bering_D029-oht_mean_bering_D000,volume_mean_march_D029-volume_mean_march_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D022-oht_mean_bering_D000,volume_mean_march_D022-volume_mean_march_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D030-oht_mean_bering_D000,volume_mean_march_D030-volume_mean_march_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,0].plot(oht_mean_bering_D024-oht_mean_bering_D000,volume_mean_march_D024-volume_mean_march_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,0].plot(oht_mean_bering_D023-oht_mean_bering_D000,volume_mean_march_D023-volume_mean_march_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,0].plot(oht_mean_bering_D025-oht_mean_bering_D000,volume_mean_march_D025-volume_mean_march_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,0].set_xlabel('$\Delta$OHT$_{Bering}$  (TW)',fontsize=20)
ax[1,0].set_ylabel('$\Delta$SIV$_{Arctic,March}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,0].set_xticks(np.arange(0, 30.1, 5))
ax[1,0].set_yticks(np.arange(-15, 1, 3))
ax[1,0].axis([0, 30, -15, 1])
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].annotate('$a$ = '+str(int(np.round(a_bering_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{atl}$ = '+str(int(np.round(a_bering_atl_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{pac}$ = '+str(int(np.round(a_bering_pac_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Bering OHT - September Arctic SIV
ax[1,1].plot(oht_mean_bering_D013-oht_mean_bering_D000,volume_mean_sept_D013-volume_mean_sept_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D012-oht_mean_bering_D000,volume_mean_sept_D012-volume_mean_sept_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D014-oht_mean_bering_D000,volume_mean_sept_D014-volume_mean_sept_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1,1].plot(oht_mean_bering_D016-oht_mean_bering_D000,volume_mean_sept_D016-volume_mean_sept_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D015-oht_mean_bering_D000,volume_mean_sept_D015-volume_mean_sept_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D017-oht_mean_bering_D000,volume_mean_sept_D017-volume_mean_sept_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1,1].plot(oht_mean_bering_D019-oht_mean_bering_D000,volume_mean_sept_D019-volume_mean_sept_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D018-oht_mean_bering_D000,volume_mean_sept_D018-volume_mean_sept_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D020-oht_mean_bering_D000,volume_mean_sept_D020-volume_mean_sept_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1,1].plot(oht_mean_bering_D027-oht_mean_bering_D000,volume_mean_sept_D027-volume_mean_sept_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D021-oht_mean_bering_D000,volume_mean_sept_D021-volume_mean_sept_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D028-oht_mean_bering_D000,volume_mean_sept_D028-volume_mean_sept_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1,1].plot(oht_mean_bering_D029-oht_mean_bering_D000,volume_mean_sept_D029-volume_mean_sept_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D022-oht_mean_bering_D000,volume_mean_sept_D022-volume_mean_sept_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D030-oht_mean_bering_D000,volume_mean_sept_D030-volume_mean_sept_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1,1].plot(oht_mean_bering_D024-oht_mean_bering_D000,volume_mean_sept_D024-volume_mean_sept_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1,1].plot(oht_mean_bering_D023-oht_mean_bering_D000,volume_mean_sept_D023-volume_mean_sept_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1,1].plot(oht_mean_bering_D025-oht_mean_bering_D000,volume_mean_sept_D025-volume_mean_sept_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1,1].set_xlabel('$\Delta$OHT$_{Bering}$  (TW)',fontsize=20)
ax[1,1].set_ylabel('$\Delta$SIV$_{Arctic,September}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,1].set_xticks(np.arange(0, 30.1, 5))
ax[1,1].set_yticks(np.arange(-15, 1, 3))
ax[1,1].axis([0, 30, -15, 1])
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].annotate('$a$ = '+str(int(np.round(a_bering_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{atl}$ = '+str(int(np.round(a_bering_atl_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{pac}$ = '+str(int(np.round(a_bering_pac_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')
ax[1,1].legend(shadow=True,frameon=False,fontsize=16,bbox_to_anchor=(1,-0.2),ncol=6)

#if save_fig == True:
#    fig.savefig(dir_fig + 'fig9c')


# Fig. 16 - Net surface heat flux in the Atlantic - OHT Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(1,2,figsize=(15,7))
fig.subplots_adjust(left=0.1,bottom=0.3,right=0.95,top=0.9,wspace=0.3,hspace=None)

# Net surface heat flux in the Atlantic - Atlantic OHT
ax[0].plot(qtoce_atl_D013-qtoce_atl_D000,oht_mean_atl_D013-oht_mean_atl_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[0].plot(qtoce_atl_D012-qtoce_atl_D000,oht_mean_atl_D012-oht_mean_atl_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[0].plot(qtoce_atl_D014-qtoce_atl_D000,oht_mean_atl_D014-oht_mean_atl_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[0].plot(qtoce_atl_D016-qtoce_atl_D000,oht_mean_atl_D016-oht_mean_atl_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[0].plot(qtoce_atl_D015-qtoce_atl_D000,oht_mean_atl_D015-oht_mean_atl_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[0].plot(qtoce_atl_D017-qtoce_atl_D000,oht_mean_atl_D017-oht_mean_atl_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[0].plot(qtoce_atl_D019-qtoce_atl_D000,oht_mean_atl_D019-oht_mean_atl_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[0].plot(qtoce_atl_D018-qtoce_atl_D000,oht_mean_atl_D018-oht_mean_atl_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[0].plot(qtoce_atl_D020-qtoce_atl_D000,oht_mean_atl_D020-oht_mean_atl_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[0].plot(qtoce_atl_D027-qtoce_atl_D000,oht_mean_atl_D027-oht_mean_atl_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[0].plot(qtoce_atl_D021-qtoce_atl_D000,oht_mean_atl_D021-oht_mean_atl_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[0].plot(qtoce_atl_D028-qtoce_atl_D000,oht_mean_atl_D028-oht_mean_atl_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[0].plot(qtoce_atl_D029-qtoce_atl_D000,oht_mean_atl_D029-oht_mean_atl_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[0].plot(qtoce_atl_D022-qtoce_atl_D000,oht_mean_atl_D022-oht_mean_atl_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[0].plot(qtoce_atl_D030-qtoce_atl_D000,oht_mean_atl_D030-oht_mean_atl_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[0].plot(qtoce_atl_D024-qtoce_atl_D000,oht_mean_atl_D024-oht_mean_atl_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[0].plot(qtoce_atl_D023-qtoce_atl_D000,oht_mean_atl_D023-oht_mean_atl_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[0].plot(qtoce_atl_D025-qtoce_atl_D000,oht_mean_atl_D025-oht_mean_atl_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[0].set_xlabel('$\Delta$SHF$_{Atlantic}$ (TW)',fontsize=20)
ax[0].set_ylabel('$\Delta$OHT$_{Atlantic}$ (TW)',fontsize=20)
ax[0].set_yticks(np.arange(0, 140.1, 20))
ax[0].set_xticks(np.arange(-150, 100.1, 50))
ax[0].axis([-150, 100, 0, 140])
ax[0].tick_params(axis='both',labelsize=16)
ax[0].grid(linestyle='--')
ax[0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Net surface heat flux in the Pacific - Bering Strait OHT
ax[1].plot(qtoce_pac_D013-qtoce_pac_D000,oht_mean_bering_D013-oht_mean_bering_D000,'o',color='lightcoral',markersize=10,label='ATL1+1$^\circ$C')
ax[1].plot(qtoce_pac_D012-qtoce_pac_D000,oht_mean_bering_D012-oht_mean_bering_D000,'ro',markersize=10,label='ATL1+3$^\circ$C')
ax[1].plot(qtoce_pac_D014-qtoce_pac_D000,oht_mean_bering_D014-oht_mean_bering_D000,'o',color='purple',markersize=10,label='ATL1+5$^\circ$C')
ax[1].plot(qtoce_pac_D016-qtoce_pac_D000,oht_mean_bering_D016-oht_mean_bering_D000,'o',color='lightblue',markersize=10,label='ATL2+1$^\circ$C')
ax[1].plot(qtoce_pac_D015-qtoce_pac_D000,oht_mean_bering_D015-oht_mean_bering_D000,'bo',markersize=10,label='ATL2+3$^\circ$C')
ax[1].plot(qtoce_pac_D017-qtoce_pac_D000,oht_mean_bering_D017-oht_mean_bering_D000,'o',color='darkblue',markersize=10,label='ATL2+5$^\circ$C')
ax[1].plot(qtoce_pac_D019-qtoce_pac_D000,oht_mean_bering_D019-oht_mean_bering_D000,'o',color='lightgreen',markersize=10,label='ATL3+1$^\circ$C')
ax[1].plot(qtoce_pac_D018-qtoce_pac_D000,oht_mean_bering_D018-oht_mean_bering_D000,'go',markersize=10,label='ATL3+3$^\circ$C')
ax[1].plot(qtoce_pac_D020-qtoce_pac_D000,oht_mean_bering_D020-oht_mean_bering_D000,'o',color='darkgreen',markersize=10,label='ATL3+5$^\circ$C')
ax[1].plot(qtoce_pac_D027-qtoce_pac_D000,oht_mean_bering_D027-oht_mean_bering_D000,'x',color='lightcoral',markersize=10,label='PAC1+1$^\circ$C')
ax[1].plot(qtoce_pac_D021-qtoce_pac_D000,oht_mean_bering_D021-oht_mean_bering_D000,'rx',markersize=10,label='PAC1+3$^\circ$C')
ax[1].plot(qtoce_pac_D028-qtoce_pac_D000,oht_mean_bering_D028-oht_mean_bering_D000,'x',color='purple',markersize=10,label='PAC1+5$^\circ$C')
ax[1].plot(qtoce_pac_D029-qtoce_pac_D000,oht_mean_bering_D029-oht_mean_bering_D000,'x',color='lightblue',markersize=10,label='PAC1+1$^\circ$C')
ax[1].plot(qtoce_pac_D022-qtoce_pac_D000,oht_mean_bering_D022-oht_mean_bering_D000,'bx',markersize=10,label='PAC2+3$^\circ$C')
ax[1].plot(qtoce_pac_D030-qtoce_pac_D000,oht_mean_bering_D030-oht_mean_bering_D000,'x',color='darkblue',markersize=10,label='PAC2+5$^\circ$C')
ax[1].plot(qtoce_pac_D024-qtoce_pac_D000,oht_mean_bering_D024-oht_mean_bering_D000,'x',color='lightgreen',markersize=10,label='PAC3+1$^\circ$C')
ax[1].plot(qtoce_pac_D023-qtoce_pac_D000,oht_mean_bering_D023-oht_mean_bering_D000,'gx',markersize=10,label='PAC3+3$^\circ$C')
ax[1].plot(qtoce_pac_D025-qtoce_pac_D000,oht_mean_bering_D025-oht_mean_bering_D000,'x',color='darkgreen',markersize=10,label='PAC3+5$^\circ$C')
ax[1].set_xlabel('$\Delta$SHF$_{Pacific}$ (TW)',fontsize=20)
ax[1].set_ylabel('$\Delta$OHT$_{Bering}$ (TW)',fontsize=20)
ax[1].set_yticks(np.arange(0, 20.1, 4))
ax[1].set_xticks(np.arange(-200, 50.1, 50))
ax[1].axis([-200, 50, 0, 20])
ax[1].tick_params(axis='both',labelsize=16)
ax[1].grid(linestyle='--')
ax[1].set_title('b',loc='left',fontsize=25,fontweight='bold')
ax[1].legend(shadow=True,frameon=False,fontsize=16,bbox_to_anchor=(1,-0.2),ncol=6)

if save_fig == True:
    fig.savefig(dir_fig + 'fig16.png')
    fig.savefig(dir_fig + 'fig16.eps',dpi=300)
