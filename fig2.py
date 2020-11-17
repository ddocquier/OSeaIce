#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 2: Plot OHT at Arctic straits
    OHT computed for each transect via compute_oht_transect.py, which follows compute_oht.py 
    Table 1: Mean OHT
PROGRAMMER
    D. Docquier
LAST UPDATE
    17/11/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Conversion factor to use if OHT was computed using cst values from CDFTOOLS
rho_lien = 1027.
cp_lien = 3985.
rho_cdftools = 1000.
cp_cdftools = 4000.
conv_fac = (rho_lien * cp_lien) / (rho_cdftools * cp_cdftools)

# Options
save_fig = True
save_var = True

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

# Working directories
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_D000 = 'D000/OHT_transects/' # CTRL
dir_D012 = 'D012/OHT_transects/' # ATL1+3K
dir_D013 = 'D013/OHT_transects/' # ATL1+1K
dir_D014 = 'D014/OHT_transects/' # ATL1+5K
dir_D015 = 'D015/OHT_transects/' # ATL2+3K
dir_D016 = 'D016/OHT_transects/' # ATL2+1K
dir_D017 = 'D017/OHT_transects/' # ATL2+5K
dir_D018 = 'D018/OHT_transects/' # ATL3+3K
dir_D019 = 'D019/OHT_transects/' # ATL3+1K
dir_D020 = 'D020/OHT_transects/' # ATL3+5K
dir_D021 = 'D021/OHT_transects/' # PAC1+3K
dir_D022 = 'D022/OHT_transects/' # PAC2+3K
dir_D023 = 'D023/OHT_transects/' # PAC3+3K
dir_D024 = 'D024/OHT_transects/' # PAC3+1K
dir_D025 = 'D025/OHT_transects/' # PAC3+5K
dir_D026 = 'D026/OHT_transects/' # PAC3+3Kb
dir_D027 = 'D027/OHT_transects/' # PAC1+1K
dir_D028 = 'D028/OHT_transects/' # PAC1+5K
dir_D029 = 'D029/OHT_transects/' # PAC2+1K
dir_D030 = 'D030/OHT_transects/' # PAC2+5K
dir_fig = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Load OHT D000
oht_full_barents_D000 = np.load(dir_output + dir_D000 + 'oht_full_barents_D000.npy')
oht_bering_D000 = np.load(dir_output + dir_D000 + 'oht_bering_D000.npy')
oht_fram_D000 = np.load(dir_output + dir_D000 + 'oht_fram_D000.npy')
oht_davis_D000 = np.load(dir_output + dir_D000 + 'oht_davis_D000.npy')
oht_total_D000 = oht_full_barents_D000 + oht_bering_D000 + oht_fram_D000 + oht_davis_D000

# Load OHT D012
oht_full_barents_D012 = np.load(dir_output + dir_D012 + 'oht_full_barents_D012.npy')
oht_bering_D012 = np.load(dir_output + dir_D012 + 'oht_bering_D012.npy')
oht_fram_D012 = np.load(dir_output + dir_D012 + 'oht_fram_D012.npy')
oht_davis_D012 = np.load(dir_output + dir_D012 + 'oht_davis_D012.npy')
oht_total_D012 = oht_full_barents_D012 + oht_bering_D012 + oht_fram_D012 + oht_davis_D012

# Load OHT D013
oht_full_barents_D013 = np.load(dir_output + dir_D013 + 'oht_full_barents_D013.npy')
oht_bering_D013 = np.load(dir_output + dir_D013 + 'oht_bering_D013.npy')
oht_fram_D013 = np.load(dir_output + dir_D013 + 'oht_fram_D013.npy')
oht_davis_D013 = np.load(dir_output + dir_D013 + 'oht_davis_D013.npy')
oht_total_D013 = oht_full_barents_D013 + oht_bering_D013 + oht_fram_D013 + oht_davis_D013

# Load OHT D014
oht_full_barents_D014 = np.load(dir_output + dir_D014 + 'oht_full_barents_D014.npy')
oht_bering_D014 = np.load(dir_output + dir_D014 + 'oht_bering_D014.npy')
oht_fram_D014 = np.load(dir_output + dir_D014 + 'oht_fram_D014.npy')
oht_davis_D014 = np.load(dir_output + dir_D014 + 'oht_davis_D014.npy')
oht_total_D014 = oht_full_barents_D014 + oht_bering_D014 + oht_fram_D014 + oht_davis_D014

# Load OHT D015
oht_full_barents_D015 = np.load(dir_output + dir_D015 + 'oht_full_barents_D015.npy')
oht_bering_D015 = np.load(dir_output + dir_D015 + 'oht_bering_D015.npy')
oht_fram_D015 = np.load(dir_output + dir_D015 + 'oht_fram_D015.npy')
oht_davis_D015 = np.load(dir_output + dir_D015 + 'oht_davis_D015.npy')
oht_total_D015 = oht_full_barents_D015 + oht_bering_D015 + oht_fram_D015 + oht_davis_D015

# Load OHT D016
oht_full_barents_D016 = np.load(dir_output + dir_D016 + 'oht_full_barents_D016.npy')
oht_bering_D016 = np.load(dir_output + dir_D016 + 'oht_bering_D016.npy')
oht_fram_D016 = np.load(dir_output + dir_D016 + 'oht_fram_D016.npy')
oht_davis_D016 = np.load(dir_output + dir_D016 + 'oht_davis_D016.npy')
oht_total_D016 = oht_full_barents_D016 + oht_bering_D016 + oht_fram_D016 + oht_davis_D016

# Load OHT D017
oht_full_barents_D017 = np.load(dir_output + dir_D017 + 'oht_full_barents_D017.npy')
oht_bering_D017 = np.load(dir_output + dir_D017 + 'oht_bering_D017.npy')
oht_fram_D017 = np.load(dir_output + dir_D017 + 'oht_fram_D017.npy')
oht_davis_D017 = np.load(dir_output + dir_D017 + 'oht_davis_D017.npy')
oht_total_D017 = oht_full_barents_D017 + oht_bering_D017 + oht_fram_D017 + oht_davis_D017

# Load OHT D018
oht_full_barents_D018 = np.load(dir_output + dir_D018 + 'oht_full_barents_D018.npy')
oht_bering_D018 = np.load(dir_output + dir_D018 + 'oht_bering_D018.npy')
oht_fram_D018 = np.load(dir_output + dir_D018 + 'oht_fram_D018.npy')
oht_davis_D018 = np.load(dir_output + dir_D018 + 'oht_davis_D018.npy')
oht_total_D018 = oht_full_barents_D018 + oht_bering_D018 + oht_fram_D018 + oht_davis_D018

# Load OHT D019
oht_full_barents_D019 = np.load(dir_output + dir_D019 + 'oht_full_barents_D019.npy')
oht_bering_D019 = np.load(dir_output + dir_D019 + 'oht_bering_D019.npy')
oht_fram_D019 = np.load(dir_output + dir_D019 + 'oht_fram_D019.npy')
oht_davis_D019 = np.load(dir_output + dir_D019 + 'oht_davis_D019.npy')
oht_total_D019 = oht_full_barents_D019 + oht_bering_D019 + oht_fram_D019 + oht_davis_D019

# Load OHT D020
oht_full_barents_D020 = np.load(dir_output + dir_D020 + 'oht_full_barents_D020.npy')
oht_bering_D020 = np.load(dir_output + dir_D020 + 'oht_bering_D020.npy')
oht_fram_D020 = np.load(dir_output + dir_D020 + 'oht_fram_D020.npy')
oht_davis_D020 = np.load(dir_output + dir_D020 + 'oht_davis_D020.npy')
oht_total_D020 = oht_full_barents_D020 + oht_bering_D020 + oht_fram_D020 + oht_davis_D020

# Load OHT D021
oht_full_barents_D021 = np.load(dir_output + dir_D021 + 'oht_full_barents_D021.npy')
oht_bering_D021 = np.load(dir_output + dir_D021 + 'oht_bering_D021.npy')
oht_bering_D021b = np.load(dir_output + dir_D021 + 'oht_bering_D021_othercst.npy')
oht_fram_D021 = np.load(dir_output + dir_D021 + 'oht_fram_D021.npy')
oht_davis_D021 = np.load(dir_output + dir_D021 + 'oht_davis_D021.npy')
oht_total_D021 = oht_full_barents_D021 + oht_bering_D021 + oht_fram_D021 + oht_davis_D021

# Load OHT D022
oht_full_barents_D022 = np.load(dir_output + dir_D022 + 'oht_full_barents_D022.npy')
oht_bering_D022 = np.load(dir_output + dir_D022 + 'oht_bering_D022.npy')
oht_fram_D022 = np.load(dir_output + dir_D022 + 'oht_fram_D022.npy')
oht_davis_D022 = np.load(dir_output + dir_D022 + 'oht_davis_D022.npy')
oht_total_D022 = oht_full_barents_D022 + oht_bering_D022 + oht_fram_D022 + oht_davis_D022

# Load OHT D023
oht_full_barents_D023 = np.load(dir_output + dir_D023 + 'oht_full_barents_D023.npy')
oht_bering_D023 = np.load(dir_output + dir_D023 + 'oht_bering_D023.npy')
oht_fram_D023 = np.load(dir_output + dir_D023 + 'oht_fram_D023.npy')
oht_davis_D023 = np.load(dir_output + dir_D023 + 'oht_davis_D023.npy')
oht_total_D023 = oht_full_barents_D023 + oht_bering_D023 + oht_fram_D023 + oht_davis_D023

# Load OHT D024
oht_full_barents_D024 = np.load(dir_output + dir_D024 + 'oht_full_barents_D024.npy')
oht_bering_D024 = np.load(dir_output + dir_D024 + 'oht_bering_D024.npy')
oht_fram_D024 = np.load(dir_output + dir_D024 + 'oht_fram_D024.npy')
oht_davis_D024 = np.load(dir_output + dir_D024 + 'oht_davis_D024.npy')
oht_total_D024 = oht_full_barents_D024 + oht_bering_D024 + oht_fram_D024 + oht_davis_D024

# Load OHT D025
oht_full_barents_D025 = np.load(dir_output + dir_D025 + 'oht_full_barents_D025.npy')
oht_bering_D025 = np.load(dir_output + dir_D025 + 'oht_bering_D025.npy')
oht_fram_D025 = np.load(dir_output + dir_D025 + 'oht_fram_D025.npy')
oht_davis_D025 = np.load(dir_output + dir_D025 + 'oht_davis_D025.npy')
oht_total_D025 = oht_full_barents_D025 + oht_bering_D025 + oht_fram_D025 + oht_davis_D025

# Load OHT D026
oht_full_barents_D026 = np.load(dir_output + dir_D026 + 'oht_full_barents_D026.npy')
oht_bering_D026 = np.load(dir_output + dir_D026 + 'oht_bering_D026.npy')
oht_fram_D026 = np.load(dir_output + dir_D026 + 'oht_fram_D026.npy')
oht_davis_D026 = np.load(dir_output + dir_D026 + 'oht_davis_D026.npy')
oht_total_D026 = oht_full_barents_D026 + oht_bering_D026 + oht_fram_D026 + oht_davis_D026

# Load OHT D027
oht_full_barents_D027 = np.load(dir_output + dir_D027 + 'oht_full_barents_D027.npy')
oht_bering_D027 = np.load(dir_output + dir_D027 + 'oht_bering_D027.npy')
oht_fram_D027 = np.load(dir_output + dir_D027 + 'oht_fram_D027.npy')
oht_davis_D027 = np.load(dir_output + dir_D027 + 'oht_davis_D027.npy')
oht_total_D027 = oht_full_barents_D027 + oht_bering_D027 + oht_fram_D027 + oht_davis_D027

# Load OHT D028
oht_full_barents_D028 = np.load(dir_output + dir_D028 + 'oht_full_barents_D028.npy')
oht_bering_D028 = np.load(dir_output + dir_D028 + 'oht_bering_D028.npy')
oht_fram_D028 = np.load(dir_output + dir_D028 + 'oht_fram_D028.npy')
oht_davis_D028 = np.load(dir_output + dir_D028 + 'oht_davis_D028.npy')
oht_total_D028 = oht_full_barents_D028 + oht_bering_D028 + oht_fram_D028 + oht_davis_D028

# Load OHT D029
oht_full_barents_D029 = np.load(dir_output + dir_D029 + 'oht_full_barents_D029.npy')
oht_bering_D029 = np.load(dir_output + dir_D029 + 'oht_bering_D029.npy')
oht_fram_D029 = np.load(dir_output + dir_D029 + 'oht_fram_D029.npy')
oht_davis_D029 = np.load(dir_output + dir_D029 + 'oht_davis_D029.npy')
oht_total_D029 = oht_full_barents_D029 + oht_bering_D029 + oht_fram_D029 + oht_davis_D029

# Load OHT D030
oht_full_barents_D030 = np.load(dir_output + dir_D030 + 'oht_full_barents_D030.npy')
oht_bering_D030 = np.load(dir_output + dir_D030 + 'oht_bering_D030.npy')
oht_fram_D030 = np.load(dir_output + dir_D030 + 'oht_fram_D030.npy')
oht_davis_D030 = np.load(dir_output + dir_D030 + 'oht_davis_D030.npy')
oht_total_D030 = oht_full_barents_D030 + oht_bering_D030 + oht_fram_D030 + oht_davis_D030

# Parameters
nmy = 12
nm = np.size(oht_full_barents_D000)
nyears = int(nm / nmy)

# Compute annual mean OHT D000
oht_annmean_full_barents_D000 = np.nanmean(oht_full_barents_D000*conv_fac,axis=1)
oht_annmean_bering_D000 = np.nanmean(oht_bering_D000*conv_fac,axis=1)
oht_annmean_fram_D000 = np.nanmean(oht_fram_D000*conv_fac,axis=1)
oht_annmean_davis_D000 = np.nanmean(oht_davis_D000*conv_fac,axis=1)
oht_annmean_total_D000 = np.nanmean(oht_total_D000*conv_fac,axis=1)

# Compute annual mean OHT D012
oht_annmean_full_barents_D012 = np.nanmean(oht_full_barents_D012*conv_fac,axis=1)
oht_annmean_bering_D012 = np.nanmean(oht_bering_D012*conv_fac,axis=1)
oht_annmean_fram_D012 = np.nanmean(oht_fram_D012*conv_fac,axis=1)
oht_annmean_davis_D012 = np.nanmean(oht_davis_D012*conv_fac,axis=1)
oht_annmean_total_D012 = np.nanmean(oht_total_D012*conv_fac,axis=1)

# Compute annual mean OHT D013
oht_annmean_full_barents_D013 = np.nanmean(oht_full_barents_D013*conv_fac,axis=1)
oht_annmean_bering_D013 = np.nanmean(oht_bering_D013*conv_fac,axis=1)
oht_annmean_fram_D013 = np.nanmean(oht_fram_D013*conv_fac,axis=1)
oht_annmean_davis_D013 = np.nanmean(oht_davis_D013*conv_fac,axis=1)
oht_annmean_total_D013 = np.nanmean(oht_total_D013*conv_fac,axis=1)

# Compute annual mean OHT D014
oht_annmean_full_barents_D014 = np.nanmean(oht_full_barents_D014*conv_fac,axis=1)
oht_annmean_bering_D014 = np.nanmean(oht_bering_D014*conv_fac,axis=1)
oht_annmean_fram_D014 = np.nanmean(oht_fram_D014*conv_fac,axis=1)
oht_annmean_davis_D014 = np.nanmean(oht_davis_D014*conv_fac,axis=1)
oht_annmean_total_D014 = np.nanmean(oht_total_D014*conv_fac,axis=1)

# Compute annual mean OHT D015
oht_annmean_full_barents_D015 = np.nanmean(oht_full_barents_D015*conv_fac,axis=1)
oht_annmean_bering_D015 = np.nanmean(oht_bering_D015*conv_fac,axis=1)
oht_annmean_fram_D015 = np.nanmean(oht_fram_D015*conv_fac,axis=1)
oht_annmean_davis_D015 = np.nanmean(oht_davis_D015*conv_fac,axis=1)
oht_annmean_total_D015 = np.nanmean(oht_total_D015*conv_fac,axis=1)

# Compute annual mean OHT D016
oht_annmean_full_barents_D016 = np.nanmean(oht_full_barents_D016*conv_fac,axis=1)
oht_annmean_bering_D016 = np.nanmean(oht_bering_D016*conv_fac,axis=1)
oht_annmean_fram_D016 = np.nanmean(oht_fram_D016*conv_fac,axis=1)
oht_annmean_davis_D016 = np.nanmean(oht_davis_D016*conv_fac,axis=1)
oht_annmean_total_D016 = np.nanmean(oht_total_D016*conv_fac,axis=1)

# Compute annual mean OHT D017
oht_annmean_full_barents_D017 = np.nanmean(oht_full_barents_D017*conv_fac,axis=1)
oht_annmean_bering_D017 = np.nanmean(oht_bering_D017*conv_fac,axis=1)
oht_annmean_fram_D017 = np.nanmean(oht_fram_D017*conv_fac,axis=1)
oht_annmean_davis_D017 = np.nanmean(oht_davis_D017*conv_fac,axis=1)
oht_annmean_total_D017 = np.nanmean(oht_total_D017*conv_fac,axis=1)

# Compute annual mean OHT D018
oht_annmean_full_barents_D018 = np.nanmean(oht_full_barents_D018*conv_fac,axis=1)
oht_annmean_bering_D018 = np.nanmean(oht_bering_D018*conv_fac,axis=1)
oht_annmean_fram_D018 = np.nanmean(oht_fram_D018*conv_fac,axis=1)
oht_annmean_davis_D018 = np.nanmean(oht_davis_D018*conv_fac,axis=1)
oht_annmean_total_D018 = np.nanmean(oht_total_D018*conv_fac,axis=1)

# Compute annual mean OHT D019
oht_annmean_full_barents_D019 = np.nanmean(oht_full_barents_D019*conv_fac,axis=1)
oht_annmean_bering_D019 = np.nanmean(oht_bering_D019*conv_fac,axis=1)
oht_annmean_fram_D019 = np.nanmean(oht_fram_D019*conv_fac,axis=1)
oht_annmean_davis_D019 = np.nanmean(oht_davis_D019*conv_fac,axis=1)
oht_annmean_total_D019 = np.nanmean(oht_total_D019*conv_fac,axis=1)

# Compute annual mean OHT D020
oht_annmean_full_barents_D020 = np.nanmean(oht_full_barents_D020*conv_fac,axis=1)
oht_annmean_bering_D020 = np.nanmean(oht_bering_D020*conv_fac,axis=1)
oht_annmean_fram_D020 = np.nanmean(oht_fram_D020*conv_fac,axis=1)
oht_annmean_davis_D020 = np.nanmean(oht_davis_D020*conv_fac,axis=1)
oht_annmean_total_D020 = np.nanmean(oht_total_D020*conv_fac,axis=1)

# Compute annual mean OHT D021
oht_annmean_full_barents_D021 = np.nanmean(oht_full_barents_D021*conv_fac,axis=1)
oht_annmean_bering_D021 = np.nanmean(oht_bering_D021*conv_fac,axis=1)
oht_annmean_fram_D021 = np.nanmean(oht_fram_D021*conv_fac,axis=1)
oht_annmean_davis_D021 = np.nanmean(oht_davis_D021*conv_fac,axis=1)
oht_annmean_total_D021 = np.nanmean(oht_total_D021*conv_fac,axis=1)

# Compute annual mean OHT D022
oht_annmean_full_barents_D022 = np.nanmean(oht_full_barents_D022*conv_fac,axis=1)
oht_annmean_bering_D022 = np.nanmean(oht_bering_D022*conv_fac,axis=1)
oht_annmean_fram_D022 = np.nanmean(oht_fram_D022*conv_fac,axis=1)
oht_annmean_davis_D022 = np.nanmean(oht_davis_D022*conv_fac,axis=1)
oht_annmean_total_D022 = np.nanmean(oht_total_D022*conv_fac,axis=1)

# Compute annual mean OHT D023
oht_annmean_full_barents_D023 = np.nanmean(oht_full_barents_D023*conv_fac,axis=1)
oht_annmean_bering_D023 = np.nanmean(oht_bering_D023*conv_fac,axis=1)
oht_annmean_fram_D023 = np.nanmean(oht_fram_D023*conv_fac,axis=1)
oht_annmean_davis_D023 = np.nanmean(oht_davis_D023*conv_fac,axis=1)
oht_annmean_total_D023 = np.nanmean(oht_total_D023*conv_fac,axis=1)

# Compute annual mean OHT D024
oht_annmean_full_barents_D024 = np.nanmean(oht_full_barents_D024*conv_fac,axis=1)
oht_annmean_bering_D024 = np.nanmean(oht_bering_D024*conv_fac,axis=1)
oht_annmean_fram_D024 = np.nanmean(oht_fram_D024*conv_fac,axis=1)
oht_annmean_davis_D024 = np.nanmean(oht_davis_D024*conv_fac,axis=1)
oht_annmean_total_D024 = np.nanmean(oht_total_D024*conv_fac,axis=1)

# Compute annual mean OHT D025
oht_annmean_full_barents_D025 = np.nanmean(oht_full_barents_D025*conv_fac,axis=1)
oht_annmean_bering_D025 = np.nanmean(oht_bering_D025*conv_fac,axis=1)
oht_annmean_fram_D025 = np.nanmean(oht_fram_D025*conv_fac,axis=1)
oht_annmean_davis_D025 = np.nanmean(oht_davis_D025*conv_fac,axis=1)
oht_annmean_total_D025 = np.nanmean(oht_total_D025*conv_fac,axis=1)

# Compute annual mean OHT D026
oht_annmean_full_barents_D026 = np.nanmean(oht_full_barents_D026*conv_fac,axis=1)
oht_annmean_bering_D026 = np.nanmean(oht_bering_D026*conv_fac,axis=1)
oht_annmean_fram_D026 = np.nanmean(oht_fram_D026*conv_fac,axis=1)
oht_annmean_davis_D026 = np.nanmean(oht_davis_D026*conv_fac,axis=1)
oht_annmean_total_D026 = np.nanmean(oht_total_D026*conv_fac,axis=1)

# Compute annual mean OHT D027
oht_annmean_full_barents_D027 = np.nanmean(oht_full_barents_D027*conv_fac,axis=1)
oht_annmean_bering_D027 = np.nanmean(oht_bering_D027*conv_fac,axis=1)
oht_annmean_fram_D027 = np.nanmean(oht_fram_D027*conv_fac,axis=1)
oht_annmean_davis_D027 = np.nanmean(oht_davis_D027*conv_fac,axis=1)
oht_annmean_total_D027 = np.nanmean(oht_total_D027*conv_fac,axis=1)

# Compute annual mean OHT D028
oht_annmean_full_barents_D028 = np.nanmean(oht_full_barents_D028*conv_fac,axis=1)
oht_annmean_bering_D028 = np.nanmean(oht_bering_D028*conv_fac,axis=1)
oht_annmean_fram_D028 = np.nanmean(oht_fram_D028*conv_fac,axis=1)
oht_annmean_davis_D028 = np.nanmean(oht_davis_D028*conv_fac,axis=1)
oht_annmean_total_D028 = np.nanmean(oht_total_D028*conv_fac,axis=1)

# Compute annual mean OHT D029
oht_annmean_full_barents_D029 = np.nanmean(oht_full_barents_D029*conv_fac,axis=1)
oht_annmean_bering_D029 = np.nanmean(oht_bering_D029*conv_fac,axis=1)
oht_annmean_fram_D029 = np.nanmean(oht_fram_D029*conv_fac,axis=1)
oht_annmean_davis_D029 = np.nanmean(oht_davis_D029*conv_fac,axis=1)
oht_annmean_total_D029 = np.nanmean(oht_total_D029*conv_fac,axis=1)

# Compute annual mean OHT D030
oht_annmean_full_barents_D030 = np.nanmean(oht_full_barents_D030*conv_fac,axis=1)
oht_annmean_bering_D030 = np.nanmean(oht_bering_D030*conv_fac,axis=1)
oht_annmean_fram_D030 = np.nanmean(oht_fram_D030*conv_fac,axis=1)
oht_annmean_davis_D030 = np.nanmean(oht_davis_D030*conv_fac,axis=1)
oht_annmean_total_D030 = np.nanmean(oht_total_D030*conv_fac,axis=1)

# Table 1 - Print mean numbers - Total Arctic
print('Total Arctic')
print('CTRL = ',np.round(np.nanmean(oht_annmean_total_D000),1),'+/-',np.round(np.nanstd(oht_annmean_total_D000),1),' TW')
print('ATL1+1K = ',np.round(np.nanmean(oht_annmean_total_D013),1),'+/-',np.round(np.nanstd(oht_annmean_total_D013),1),' TW',compute_sig(nyears,oht_annmean_total_D013,oht_annmean_total_D000))
print('ATL1+3K = ',np.round(np.nanmean(oht_annmean_total_D012),1),'+/-',np.round(np.nanstd(oht_annmean_total_D012),1),' TW',compute_sig(nyears,oht_annmean_total_D012,oht_annmean_total_D000))
print('ATL1+5K = ',np.round(np.nanmean(oht_annmean_total_D014),1),'+/-',np.round(np.nanstd(oht_annmean_total_D014),1),' TW',compute_sig(nyears,oht_annmean_total_D014,oht_annmean_total_D000))
print('ATL2+1K = ',np.round(np.nanmean(oht_annmean_total_D016),1),'+/-',np.round(np.nanstd(oht_annmean_total_D016),1),' TW',compute_sig(nyears,oht_annmean_total_D016,oht_annmean_total_D000))
print('ATL2+3K = ',np.round(np.nanmean(oht_annmean_total_D015),1),'+/-',np.round(np.nanstd(oht_annmean_total_D015),1),' TW',compute_sig(nyears,oht_annmean_total_D015,oht_annmean_total_D000))
print('ATL2+5K = ',np.round(np.nanmean(oht_annmean_total_D017),1),'+/-',np.round(np.nanstd(oht_annmean_total_D017),1),' TW',compute_sig(nyears,oht_annmean_total_D017,oht_annmean_total_D000))
print('ATL3+1K = ',np.round(np.nanmean(oht_annmean_total_D019),1),'+/-',np.round(np.nanstd(oht_annmean_total_D019),1),' TW',compute_sig(nyears,oht_annmean_total_D019,oht_annmean_total_D000))
print('ATL3+3K = ',np.round(np.nanmean(oht_annmean_total_D018),1),'+/-',np.round(np.nanstd(oht_annmean_total_D018),1),' TW',compute_sig(nyears,oht_annmean_total_D018,oht_annmean_total_D000))
print('ATL3+5K = ',np.round(np.nanmean(oht_annmean_total_D020),1),'+/-',np.round(np.nanstd(oht_annmean_total_D020),1),' TW',compute_sig(nyears,oht_annmean_total_D020,oht_annmean_total_D000))
print('PAC1+1K = ',np.round(np.nanmean(oht_annmean_total_D027),1),'+/-',np.round(np.nanstd(oht_annmean_total_D027),1),' TW',compute_sig(nyears,oht_annmean_total_D027,oht_annmean_total_D000))
print('PAC1+3K = ',np.round(np.nanmean(oht_annmean_total_D021),1),'+/-',np.round(np.nanstd(oht_annmean_total_D021),1),' TW',compute_sig(nyears,oht_annmean_total_D021,oht_annmean_total_D000))
print('PAC1+5K = ',np.round(np.nanmean(oht_annmean_total_D028),1),'+/-',np.round(np.nanstd(oht_annmean_total_D028),1),' TW',compute_sig(nyears,oht_annmean_total_D028,oht_annmean_total_D000))
print('PAC2+1K = ',np.round(np.nanmean(oht_annmean_total_D029),1),'+/-',np.round(np.nanstd(oht_annmean_total_D029),1),' TW',compute_sig(nyears,oht_annmean_total_D029,oht_annmean_total_D000))
print('PAC2+3K = ',np.round(np.nanmean(oht_annmean_total_D022),1),'+/-',np.round(np.nanstd(oht_annmean_total_D022),1),' TW',compute_sig(nyears,oht_annmean_total_D022,oht_annmean_total_D000))
print('PAC2+5K = ',np.round(np.nanmean(oht_annmean_total_D030),1),'+/-',np.round(np.nanstd(oht_annmean_total_D030),1),' TW',compute_sig(nyears,oht_annmean_total_D030,oht_annmean_total_D000))
print('PAC3+1K = ',np.round(np.nanmean(oht_annmean_total_D024),1),'+/-',np.round(np.nanstd(oht_annmean_total_D024),1),' TW',compute_sig(nyears,oht_annmean_total_D024,oht_annmean_total_D000))
print('PAC3+3K = ',np.round(np.nanmean(oht_annmean_total_D023),1),'+/-',np.round(np.nanstd(oht_annmean_total_D023),1),' TW',compute_sig(nyears,oht_annmean_total_D023,oht_annmean_total_D000))
print('PAC3+3Kb = ',np.round(np.nanmean(oht_annmean_total_D026),1),'+/-',np.round(np.nanstd(oht_annmean_total_D026),1),' TW',compute_sig(nyears,oht_annmean_total_D026,oht_annmean_total_D000))
print('PAC3+5K = ',np.round(np.nanmean(oht_annmean_total_D025),1),'+/-',np.round(np.nanstd(oht_annmean_total_D025),1),' TW',compute_sig(nyears,oht_annmean_total_D025,oht_annmean_total_D000))
print('PAC3+3K - PAC2+3K = ',np.round(np.nanmean(oht_annmean_total_D023)-np.nanmean(oht_annmean_total_D022),1),' TW',compute_sig(nyears,oht_annmean_total_D023,oht_annmean_total_D022))
print('PAC3+3K - PAC3+5K = ',np.round(np.nanmean(oht_annmean_total_D023)-np.nanmean(oht_annmean_total_D025),1),' TW',compute_sig(nyears,oht_annmean_total_D023,oht_annmean_total_D025))
print('ATL3+3K - ATL3+1K = ',np.round(np.nanmean(oht_annmean_total_D018)-np.nanmean(oht_annmean_total_D019),1),' TW',compute_sig(nyears,oht_annmean_total_D018,oht_annmean_total_D019))
print('------------')
print('BSO')
print('CTRL = ',np.round(np.nanmean(oht_annmean_full_barents_D000),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D000),1),' TW')
print('ATL1+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D013),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D013),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D013,oht_annmean_full_barents_D000))
print('ATL1+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D012),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D012),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D012,oht_annmean_full_barents_D000))
print('ATL1+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D014),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D014),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D014,oht_annmean_full_barents_D000))
print('ATL2+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D016),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D016),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D016,oht_annmean_full_barents_D000))
print('ATL2+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D015),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D015),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D015,oht_annmean_full_barents_D000))
print('ATL2+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D017),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D017),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D017,oht_annmean_full_barents_D000))
print('ATL3+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D019),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D019),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D019,oht_annmean_full_barents_D000))
print('ATL3+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D018),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D018),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D018,oht_annmean_full_barents_D000))
print('ATL3+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D020),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D020),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D020,oht_annmean_full_barents_D000))
print('PAC1+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D027),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D027),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D027,oht_annmean_full_barents_D000))
print('PAC1+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D021),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D021),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D021,oht_annmean_full_barents_D000))
print('PAC1+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D028),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D028),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D028,oht_annmean_full_barents_D000))
print('PAC2+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D029),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D029),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D029,oht_annmean_full_barents_D000))
print('PAC2+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D022),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D022),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D022,oht_annmean_full_barents_D000))
print('PAC2+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D030),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D030),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D030,oht_annmean_full_barents_D000))
print('PAC3+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D024),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D024),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D024,oht_annmean_full_barents_D000))
print('PAC3+3K = ',np.round(np.nanmean(oht_annmean_full_barents_D023),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D023),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D023,oht_annmean_full_barents_D000))
print('PAC3+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D025),1),'+/-',np.round(np.nanstd(oht_annmean_full_barents_D025),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D025,oht_annmean_full_barents_D000))
print('PAC3+3K - PAC3+5K = ',np.round(np.nanmean(oht_annmean_full_barents_D023)-np.nanmean(oht_annmean_full_barents_D025),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D023,oht_annmean_full_barents_D025))
print('ATL3+3K - ATL3+1K = ',np.round(np.nanmean(oht_annmean_full_barents_D018)-np.nanmean(oht_annmean_full_barents_D019),1),' TW',compute_sig(nyears,oht_annmean_full_barents_D018,oht_annmean_full_barents_D019))
print('------------')
print('Bering Strait')
print('CTRL = ',np.round(np.nanmean(oht_annmean_bering_D000),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D000),1),' TW')
print('ATL1+1K = ',np.round(np.nanmean(oht_annmean_bering_D013),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D013),1),' TW',compute_sig(nyears,oht_annmean_bering_D013,oht_annmean_bering_D000))
print('ATL1+3K = ',np.round(np.nanmean(oht_annmean_bering_D012),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D012),1),' TW',compute_sig(nyears,oht_annmean_bering_D012,oht_annmean_bering_D000))
print('ATL1+5K = ',np.round(np.nanmean(oht_annmean_bering_D014),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D014),1),' TW',compute_sig(nyears,oht_annmean_bering_D014,oht_annmean_bering_D000))
print('ATL2+1K = ',np.round(np.nanmean(oht_annmean_bering_D016),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D016),1),' TW',compute_sig(nyears,oht_annmean_bering_D016,oht_annmean_bering_D000))
print('ATL2+3K = ',np.round(np.nanmean(oht_annmean_bering_D015),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D015),1),' TW',compute_sig(nyears,oht_annmean_bering_D015,oht_annmean_bering_D000))
print('ATL2+5K = ',np.round(np.nanmean(oht_annmean_bering_D017),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D017),1),' TW',compute_sig(nyears,oht_annmean_bering_D017,oht_annmean_bering_D000))
print('ATL3+1K = ',np.round(np.nanmean(oht_annmean_bering_D019),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D019),1),' TW',compute_sig(nyears,oht_annmean_bering_D019,oht_annmean_bering_D000))
print('ATL3+3K = ',np.round(np.nanmean(oht_annmean_bering_D018),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D018),1),' TW',compute_sig(nyears,oht_annmean_bering_D018,oht_annmean_bering_D000))
print('ATL3+5K = ',np.round(np.nanmean(oht_annmean_bering_D020),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D020),1),' TW',compute_sig(nyears,oht_annmean_bering_D020,oht_annmean_bering_D000))
print('PAC1+1K = ',np.round(np.nanmean(oht_annmean_bering_D027),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D027),1),' TW',compute_sig(nyears,oht_annmean_bering_D027,oht_annmean_bering_D000))
print('PAC1+3K = ',np.round(np.nanmean(oht_annmean_bering_D021),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D021),1),' TW',compute_sig(nyears,oht_annmean_bering_D021,oht_annmean_bering_D000))
print('PAC1+5K = ',np.round(np.nanmean(oht_annmean_bering_D028),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D028),1),' TW',compute_sig(nyears,oht_annmean_bering_D028,oht_annmean_bering_D000))
print('PAC2+1K = ',np.round(np.nanmean(oht_annmean_bering_D029),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D029),1),' TW',compute_sig(nyears,oht_annmean_bering_D029,oht_annmean_bering_D000))
print('PAC2+3K = ',np.round(np.nanmean(oht_annmean_bering_D022),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D022),1),' TW',compute_sig(nyears,oht_annmean_bering_D022,oht_annmean_bering_D000))
print('PAC2+5K = ',np.round(np.nanmean(oht_annmean_bering_D030),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D030),1),' TW',compute_sig(nyears,oht_annmean_bering_D030,oht_annmean_bering_D000))
print('PAC3+1K = ',np.round(np.nanmean(oht_annmean_bering_D024),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D024),1),' TW',compute_sig(nyears,oht_annmean_bering_D024,oht_annmean_bering_D000))
print('PAC3+3K = ',np.round(np.nanmean(oht_annmean_bering_D023),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D023),1),' TW',compute_sig(nyears,oht_annmean_bering_D023,oht_annmean_bering_D000))
print('PAC3+5K = ',np.round(np.nanmean(oht_annmean_bering_D025),1),'+/-',np.round(np.nanstd(oht_annmean_bering_D025),1),' TW',compute_sig(nyears,oht_annmean_bering_D025,oht_annmean_bering_D000))
print('------------')
print('Fram Strait')
print('CTRL = ',np.round(np.nanmean(oht_annmean_fram_D000),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D000),1),' TW')
print('ATL1+1K = ',np.round(np.nanmean(oht_annmean_fram_D013),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D013),1),' TW',compute_sig(nyears,oht_annmean_fram_D013,oht_annmean_fram_D000))
print('ATL1+3K = ',np.round(np.nanmean(oht_annmean_fram_D012),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D012),1),' TW',compute_sig(nyears,oht_annmean_fram_D012,oht_annmean_fram_D000))
print('ATL1+5K = ',np.round(np.nanmean(oht_annmean_fram_D014),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D014),1),' TW',compute_sig(nyears,oht_annmean_fram_D014,oht_annmean_fram_D000))
print('ATL2+1K = ',np.round(np.nanmean(oht_annmean_fram_D016),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D016),1),' TW',compute_sig(nyears,oht_annmean_fram_D016,oht_annmean_fram_D000))
print('ATL2+3K = ',np.round(np.nanmean(oht_annmean_fram_D015),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D015),1),' TW',compute_sig(nyears,oht_annmean_fram_D015,oht_annmean_fram_D000))
print('ATL2+5K = ',np.round(np.nanmean(oht_annmean_fram_D017),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D017),1),' TW',compute_sig(nyears,oht_annmean_fram_D017,oht_annmean_fram_D000))
print('ATL3+1K = ',np.round(np.nanmean(oht_annmean_fram_D019),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D019),1),' TW',compute_sig(nyears,oht_annmean_fram_D019,oht_annmean_fram_D000))
print('ATL3+3K = ',np.round(np.nanmean(oht_annmean_fram_D018),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D018),1),' TW',compute_sig(nyears,oht_annmean_fram_D018,oht_annmean_fram_D000))
print('ATL3+5K = ',np.round(np.nanmean(oht_annmean_fram_D020),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D020),1),' TW',compute_sig(nyears,oht_annmean_fram_D020,oht_annmean_fram_D000))
print('PAC1+1K = ',np.round(np.nanmean(oht_annmean_fram_D027),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D027),1),' TW',compute_sig(nyears,oht_annmean_fram_D027,oht_annmean_fram_D000))
print('PAC1+3K = ',np.round(np.nanmean(oht_annmean_fram_D021),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D021),1),' TW',compute_sig(nyears,oht_annmean_fram_D021,oht_annmean_fram_D000))
print('PAC1+5K = ',np.round(np.nanmean(oht_annmean_fram_D028),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D028),1),' TW',compute_sig(nyears,oht_annmean_fram_D028,oht_annmean_fram_D000))
print('PAC2+1K = ',np.round(np.nanmean(oht_annmean_fram_D029),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D029),1),' TW',compute_sig(nyears,oht_annmean_fram_D029,oht_annmean_fram_D000))
print('PAC2+3K = ',np.round(np.nanmean(oht_annmean_fram_D022),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D022),1),' TW',compute_sig(nyears,oht_annmean_fram_D022,oht_annmean_fram_D000))
print('PAC2+5K = ',np.round(np.nanmean(oht_annmean_fram_D030),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D030),1),' TW',compute_sig(nyears,oht_annmean_fram_D030,oht_annmean_fram_D000))
print('PAC3+1K = ',np.round(np.nanmean(oht_annmean_fram_D024),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D024),1),' TW',compute_sig(nyears,oht_annmean_fram_D024,oht_annmean_fram_D000))
print('PAC3+3K = ',np.round(np.nanmean(oht_annmean_fram_D023),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D023),1),' TW',compute_sig(nyears,oht_annmean_fram_D023,oht_annmean_fram_D000))
print('PAC3+5K = ',np.round(np.nanmean(oht_annmean_fram_D025),1),'+/-',np.round(np.nanstd(oht_annmean_fram_D025),1),' TW',compute_sig(nyears,oht_annmean_fram_D025,oht_annmean_fram_D000))
print('------------')
print('Davis Strait')
print('CTRL = ',np.round(np.nanmean(oht_annmean_davis_D000),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D000),1),' TW')
print('ATL1+1K = ',np.round(np.nanmean(oht_annmean_davis_D013),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D013),1),' TW',compute_sig(nyears,oht_annmean_davis_D013,oht_annmean_davis_D000))
print('ATL1+3K = ',np.round(np.nanmean(oht_annmean_davis_D012),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D012),1),' TW',compute_sig(nyears,oht_annmean_davis_D012,oht_annmean_davis_D000))
print('ATL1+5K = ',np.round(np.nanmean(oht_annmean_davis_D014),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D014),1),' TW',compute_sig(nyears,oht_annmean_davis_D014,oht_annmean_davis_D000))
print('ATL2+1K = ',np.round(np.nanmean(oht_annmean_davis_D016),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D016),1),' TW',compute_sig(nyears,oht_annmean_davis_D016,oht_annmean_davis_D000))
print('ATL2+3K = ',np.round(np.nanmean(oht_annmean_davis_D015),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D015),1),' TW',compute_sig(nyears,oht_annmean_davis_D015,oht_annmean_davis_D000))
print('ATL2+5K = ',np.round(np.nanmean(oht_annmean_davis_D017),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D017),1),' TW',compute_sig(nyears,oht_annmean_davis_D017,oht_annmean_davis_D000))
print('ATL3+1K = ',np.round(np.nanmean(oht_annmean_davis_D019),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D019),1),' TW',compute_sig(nyears,oht_annmean_davis_D019,oht_annmean_davis_D000))
print('ATL3+3K = ',np.round(np.nanmean(oht_annmean_davis_D018),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D018),1),' TW',compute_sig(nyears,oht_annmean_davis_D018,oht_annmean_davis_D000))
print('ATL3+5K = ',np.round(np.nanmean(oht_annmean_davis_D020),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D020),1),' TW',compute_sig(nyears,oht_annmean_davis_D020,oht_annmean_davis_D000))
print('PAC1+1K = ',np.round(np.nanmean(oht_annmean_davis_D027),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D027),1),' TW',compute_sig(nyears,oht_annmean_davis_D027,oht_annmean_davis_D000))
print('PAC1+3K = ',np.round(np.nanmean(oht_annmean_davis_D021),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D021),1),' TW',compute_sig(nyears,oht_annmean_davis_D021,oht_annmean_davis_D000))
print('PAC1+5K = ',np.round(np.nanmean(oht_annmean_davis_D028),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D028),1),' TW',compute_sig(nyears,oht_annmean_davis_D028,oht_annmean_davis_D000))
print('PAC2+1K = ',np.round(np.nanmean(oht_annmean_davis_D029),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D029),1),' TW',compute_sig(nyears,oht_annmean_davis_D029,oht_annmean_davis_D000))
print('PAC2+3K = ',np.round(np.nanmean(oht_annmean_davis_D022),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D022),1),' TW',compute_sig(nyears,oht_annmean_davis_D022,oht_annmean_davis_D000))
print('PAC2+5K = ',np.round(np.nanmean(oht_annmean_davis_D030),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D030),1),' TW',compute_sig(nyears,oht_annmean_davis_D030,oht_annmean_davis_D000))
print('PAC3+1K = ',np.round(np.nanmean(oht_annmean_davis_D024),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D024),1),' TW',compute_sig(nyears,oht_annmean_davis_D024,oht_annmean_davis_D000))
print('PAC3+3K = ',np.round(np.nanmean(oht_annmean_davis_D023),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D023),1),' TW',compute_sig(nyears,oht_annmean_davis_D023,oht_annmean_davis_D000))
print('PAC3+5K = ',np.round(np.nanmean(oht_annmean_davis_D025),1),'+/-',np.round(np.nanstd(oht_annmean_davis_D025),1),' TW',compute_sig(nyears,oht_annmean_davis_D025,oht_annmean_davis_D000))
print('-----------')
print('PAC3+3K - PAC2+3K = ',np.round(np.nanmean(oht_annmean_total_D023)-np.nanmean(oht_annmean_total_D022),1),' TW',compute_sig(nyears,oht_annmean_total_D023,oht_annmean_total_D022))
print('PAC3+3Kb - PAC2+3K = ',np.round(np.nanmean(oht_annmean_total_D026)-np.nanmean(oht_annmean_total_D022),1),' TW',compute_sig(nyears,oht_annmean_total_D026,oht_annmean_total_D022))
print('PAC3+3K - PAC1+3K = ',np.round(np.nanmean(oht_annmean_total_D023)-np.nanmean(oht_annmean_total_D021),1),' TW',compute_sig(nyears,oht_annmean_total_D023,oht_annmean_total_D021))
print('ATL3+3K - ATL2+3K = ',np.round(np.nanmean(oht_annmean_total_D018)-np.nanmean(oht_annmean_total_D015),1),' TW',compute_sig(nyears,oht_annmean_total_D018,oht_annmean_total_D015))
print('ATL2+3K - ATL1+3K = ',np.round(np.nanmean(oht_annmean_total_D015)-np.nanmean(oht_annmean_total_D012),1),' TW',compute_sig(nyears,oht_annmean_total_D015,oht_annmean_total_D012))

# Save annual mean OHT
if save_var == True:
    filename = dir_fig + 'OHT_annmean.npy'
    np.save(filename,[oht_annmean_total_D000,oht_annmean_total_D012,oht_annmean_total_D013,oht_annmean_total_D014,oht_annmean_total_D015,oht_annmean_total_D016,oht_annmean_total_D017,oht_annmean_total_D018,oht_annmean_total_D019,oht_annmean_total_D020,oht_annmean_total_D021,oht_annmean_total_D022,oht_annmean_total_D023,oht_annmean_total_D024,oht_annmean_total_D025,oht_annmean_total_D027,oht_annmean_total_D028,oht_annmean_total_D029,oht_annmean_total_D030])
    filename = dir_fig + 'OHT_annmean_BSO.npy'
    np.save(filename,[oht_annmean_full_barents_D000,oht_annmean_full_barents_D012,oht_annmean_full_barents_D013,oht_annmean_full_barents_D014,oht_annmean_full_barents_D015,oht_annmean_full_barents_D016,oht_annmean_full_barents_D017,oht_annmean_full_barents_D018,oht_annmean_full_barents_D019,oht_annmean_full_barents_D020,oht_annmean_full_barents_D021,oht_annmean_full_barents_D022,oht_annmean_full_barents_D023,oht_annmean_full_barents_D024,oht_annmean_full_barents_D025,oht_annmean_full_barents_D027,oht_annmean_full_barents_D028,oht_annmean_full_barents_D029,oht_annmean_full_barents_D030])
    filename = dir_fig + 'OHT_annmean_Bering.npy'
    np.save(filename,[oht_annmean_bering_D000,oht_annmean_bering_D012,oht_annmean_bering_D013,oht_annmean_bering_D014,oht_annmean_bering_D015,oht_annmean_bering_D016,oht_annmean_bering_D017,oht_annmean_bering_D018,oht_annmean_bering_D019,oht_annmean_bering_D020,oht_annmean_bering_D021,oht_annmean_bering_D022,oht_annmean_bering_D023,oht_annmean_bering_D024,oht_annmean_bering_D025,oht_annmean_bering_D027,oht_annmean_bering_D028,oht_annmean_bering_D029,oht_annmean_bering_D030])
    filename = dir_fig + 'OHT_annmean_Fram.npy'
    np.save(filename,[oht_annmean_fram_D000,oht_annmean_fram_D012,oht_annmean_fram_D013,oht_annmean_fram_D014,oht_annmean_fram_D015,oht_annmean_fram_D016,oht_annmean_fram_D017,oht_annmean_fram_D018,oht_annmean_fram_D019,oht_annmean_fram_D020,oht_annmean_fram_D021,oht_annmean_fram_D022,oht_annmean_fram_D023,oht_annmean_fram_D024,oht_annmean_fram_D025,oht_annmean_fram_D027,oht_annmean_fram_D028,oht_annmean_fram_D029,oht_annmean_fram_D030])
    filename = dir_fig + 'OHT_annmean_Davis.npy'
    np.save(filename,[oht_annmean_davis_D000,oht_annmean_davis_D012,oht_annmean_davis_D013,oht_annmean_davis_D014,oht_annmean_davis_D015,oht_annmean_davis_D016,oht_annmean_davis_D017,oht_annmean_davis_D018,oht_annmean_davis_D019,oht_annmean_davis_D020,oht_annmean_davis_D021,oht_annmean_davis_D022,oht_annmean_davis_D023,oht_annmean_davis_D024,oht_annmean_davis_D025,oht_annmean_davis_D027,oht_annmean_davis_D028,oht_annmean_davis_D029,oht_annmean_davis_D030])

# Labels
name_xticks = ['0','10','20','30','40','50']
 

# Fig. 2 - Time series of annual mean total OHT to the Arctic
fig,ax = plt.subplots(2,1,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.09,right=0.95,top=0.95,wspace=None,hspace=0.2)

# Atlantic SST experiments
ax[0].set_title('Atlantic SST+3$^\circ$C experiments',fontsize=26)
ax[0].plot(np.arange(nyears) + 1,oht_annmean_total_D000,'k-',linewidth=3,label='CTRL')
ax[0].plot(np.arange(nyears) + 1,oht_annmean_total_D012,'r-',linewidth=2,label='ATL1+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D012-oht_annmean_total_D000),1))+')')
ax[0].plot(np.arange(nyears) + 1,oht_annmean_total_D015,'-',color='gray',linewidth=2,label='ATL2+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D015-oht_annmean_total_D000),1))+')')
ax[0].plot(np.arange(nyears) + 1,oht_annmean_total_D018,'b-',linewidth=2,label='ATL3+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D018-oht_annmean_total_D000),1))+')')
ax[0].legend(loc='lower left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[0].set_ylabel('Total Arctic OHT (TW)',fontsize=24)
ax[0].set_title('a',loc='left',fontsize=24,fontweight='bold')
ax[0].set_xticks(np.arange(1, 53, 10))
ax[0].set_xticklabels(name_xticks)
ax[0].set_yticks(np.arange(0, 271, 30))
ax[0].tick_params(axis='both',labelsize=16)
ax[0].axis([0, 52, 0, 280])
ax[0].grid()

# Pacific SST experiments
ax[1].set_title('Pacific SST+3$^\circ$C experiments',fontsize=26)
ax[1].plot(np.arange(nyears) + 1,oht_annmean_total_D000,'k-',linewidth=3,label='CTRL')
ax[1].plot(np.arange(nyears) + 1,oht_annmean_total_D021,'g-',linewidth=2,label='PAC1+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D021-oht_annmean_total_D000),1))+')')
ax[1].plot(np.arange(nyears) + 1,oht_annmean_total_D022,'-',color='brown',linewidth=2,label='PAC2+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D022-oht_annmean_total_D000),1))+')')
ax[1].plot(np.arange(nyears) + 1,oht_annmean_total_D023,'-',color='orange',linewidth=2,label='PAC3+3$^\circ$C ('+str(np.round(np.nanmean(oht_annmean_total_D023-oht_annmean_total_D000),1))+')')
ax[1].legend(loc='lower left',shadow=True,frameon=False,fontsize=18,ncol=2)
ax[1].set_xlabel('Year',fontsize=24)
ax[1].set_ylabel('Total Arctic OHT (TW)',fontsize=24)
ax[1].set_title('b',loc='left',fontsize=24,fontweight='bold')
ax[1].set_xticks(np.arange(1, 53, 10))
ax[1].set_xticklabels(name_xticks)
ax[1].set_yticks(np.arange(0, 271, 30))
ax[1].tick_params(axis='both',labelsize=16)
ax[1].axis([0, 52, 0, 280])
ax[1].grid()

if save_fig == True:
    fig.savefig(dir_fig + 'fig2.png')
    fig.savefig(dir_fig + 'fig2.eps',dpi=300)
