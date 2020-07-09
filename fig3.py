#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 3: Plot OHT changes du to velocity and temperature anomalies at Arctic straits
    OHT velocity and temperature anomalies computed for each transect via decompose_oht_transect.py, which follows decompose_oht.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    09/06/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt

# Option
save_fig = False

# Function to retrieve monthly mean
def compute_monthly(nm,nmy,var):
    nyears = int(nm / nmy)
    var_mon = np.zeros((nmy,nyears))
    for i in np.arange(nmy):
        for k in np.arange(nyears):
            if np.isnan(var[i+k*nmy]) == False:
                var_mon[i,k] = var[i+k*nmy]
    return var_mon

# Working directories
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_D000 = 'D000/OHT_transects/'
dir_D012 = 'D012/OHT_transects/'
dir_D013 = 'D013/OHT_transects/'
dir_D014 = 'D014/OHT_transects/'
dir_D015 = 'D015/OHT_transects/'
dir_D016 = 'D016/OHT_transects/'
dir_D017 = 'D017/OHT_transects/'
dir_D018 = 'D018/OHT_transects/'
dir_D019 = 'D019/OHT_transects/'
dir_D020 = 'D020/OHT_transects/'
dir_D021 = 'D021/OHT_transects/'
dir_D022 = 'D022/OHT_transects/'
dir_D023 = 'D023/OHT_transects/'
dir_D024 = 'D024/OHT_transects/'
dir_D025 = 'D025/OHT_transects/'
dir_D027 = 'D027/OHT_transects/'
dir_D028 = 'D028/OHT_transects/'
dir_D029 = 'D029/OHT_transects/'
dir_D030 = 'D030/OHT_transects/'
dir_fig = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Load OHT D012
vel_barents_D012,temp_barents_D012,cov_barents_D012 = np.load(dir_output + dir_D012 + 'decoht_barents_D012.npy')
vel_bering_D012,temp_bering_D012,cov_bering_D012 = np.load(dir_output + dir_D012 + 'decoht_bering_D012.npy')
vel_fram_D012,temp_fram_D012,cov_fram_D012 = np.load(dir_output + dir_D012 + 'decoht_fram_D012.npy')
vel_davis_D012,temp_davis_D012,cov_davis_D012 = np.load(dir_output + dir_D012 + 'decoht_davis_D012.npy')
vel_total_D012 = vel_barents_D012 + vel_bering_D012 + vel_fram_D012 + vel_davis_D012
temp_total_D012 = temp_barents_D012 + temp_bering_D012 + temp_fram_D012 + temp_davis_D012
cov_total_D012 = cov_barents_D012 + cov_bering_D012 + cov_fram_D012 + cov_davis_D012

# Load OHT D013
vel_barents_D013,temp_barents_D013,cov_barents_D013 = np.load(dir_output + dir_D013 + 'decoht_barents_D013.npy')
vel_bering_D013,temp_bering_D013,cov_bering_D013 = np.load(dir_output + dir_D013 + 'decoht_bering_D013.npy')
vel_fram_D013,temp_fram_D013,cov_fram_D013 = np.load(dir_output + dir_D013 + 'decoht_fram_D013.npy')
vel_davis_D013,temp_davis_D013,cov_davis_D013 = np.load(dir_output + dir_D013 + 'decoht_davis_D013.npy')
vel_total_D013 = vel_barents_D013 + vel_bering_D013 + vel_fram_D013 + vel_davis_D013
temp_total_D013 = temp_barents_D013 + temp_bering_D013 + temp_fram_D013 + temp_davis_D013
cov_total_D013 = cov_barents_D013 + cov_bering_D013 + cov_fram_D013 + cov_davis_D013

# Load OHT D014
vel_barents_D014,temp_barents_D014,cov_barents_D014 = np.load(dir_output + dir_D014 + 'decoht_barents_D014.npy')
vel_bering_D014,temp_bering_D014,cov_bering_D014 = np.load(dir_output + dir_D014 + 'decoht_bering_D014.npy')
vel_fram_D014,temp_fram_D014,cov_fram_D014 = np.load(dir_output + dir_D014 + 'decoht_fram_D014.npy')
vel_davis_D014,temp_davis_D014,cov_davis_D014 = np.load(dir_output + dir_D014 + 'decoht_davis_D014.npy')
vel_total_D014 = vel_barents_D014 + vel_bering_D014 + vel_fram_D014 + vel_davis_D014
temp_total_D014 = temp_barents_D014 + temp_bering_D014 + temp_fram_D014 + temp_davis_D014
cov_total_D014 = cov_barents_D014 + cov_bering_D014 + cov_fram_D014 + cov_davis_D014

# Load OHT D015
vel_barents_D015,temp_barents_D015,cov_barents_D015 = np.load(dir_output + dir_D015 + 'decoht_barents_D015.npy')
vel_bering_D015,temp_bering_D015,cov_bering_D015 = np.load(dir_output + dir_D015 + 'decoht_bering_D015.npy')
vel_fram_D015,temp_fram_D015,cov_fram_D015 = np.load(dir_output + dir_D015 + 'decoht_fram_D015.npy')
vel_davis_D015,temp_davis_D015,cov_davis_D015 = np.load(dir_output + dir_D015 + 'decoht_davis_D015.npy')
vel_total_D015 = vel_barents_D015 + vel_bering_D015 + vel_fram_D015 + vel_davis_D015
temp_total_D015 = temp_barents_D015 + temp_bering_D015 + temp_fram_D015 + temp_davis_D015
cov_total_D015 = cov_barents_D015 + cov_bering_D015 + cov_fram_D015 + cov_davis_D015

# Load OHT D016
vel_barents_D016,temp_barents_D016,cov_barents_D016 = np.load(dir_output + dir_D016 + 'decoht_barents_D016.npy')
vel_bering_D016,temp_bering_D016,cov_bering_D016 = np.load(dir_output + dir_D016 + 'decoht_bering_D016.npy')
vel_fram_D016,temp_fram_D016,cov_fram_D016 = np.load(dir_output + dir_D016 + 'decoht_fram_D016.npy')
vel_davis_D016,temp_davis_D016,cov_davis_D016 = np.load(dir_output + dir_D016 + 'decoht_davis_D016.npy')
vel_total_D016 = vel_barents_D016 + vel_bering_D016 + vel_fram_D016 + vel_davis_D016
temp_total_D016 = temp_barents_D016 + temp_bering_D016 + temp_fram_D016 + temp_davis_D016
cov_total_D016 = cov_barents_D016 + cov_bering_D016 + cov_fram_D016 + cov_davis_D016

# Load OHT D017
vel_barents_D017,temp_barents_D017,cov_barents_D017 = np.load(dir_output + dir_D017 + 'decoht_barents_D017.npy')
vel_bering_D017,temp_bering_D017,cov_bering_D017 = np.load(dir_output + dir_D017 + 'decoht_bering_D017.npy')
vel_fram_D017,temp_fram_D017,cov_fram_D017 = np.load(dir_output + dir_D017 + 'decoht_fram_D017.npy')
vel_davis_D017,temp_davis_D017,cov_davis_D017 = np.load(dir_output + dir_D017 + 'decoht_davis_D017.npy')
vel_total_D017 = vel_barents_D017 + vel_bering_D017 + vel_fram_D017 + vel_davis_D017
temp_total_D017 = temp_barents_D017 + temp_bering_D017 + temp_fram_D017 + temp_davis_D017
cov_total_D017 = cov_barents_D017 + cov_bering_D017 + cov_fram_D017 + cov_davis_D017

# Load OHT D018
vel_barents_D018,temp_barents_D018,cov_barents_D018 = np.load(dir_output + dir_D018 + 'decoht_barents_D018.npy')
vel_bering_D018,temp_bering_D018,cov_bering_D018 = np.load(dir_output + dir_D018 + 'decoht_bering_D018.npy')
vel_fram_D018,temp_fram_D018,cov_fram_D018 = np.load(dir_output + dir_D018 + 'decoht_fram_D018.npy')
vel_davis_D018,temp_davis_D018,cov_davis_D018 = np.load(dir_output + dir_D018 + 'decoht_davis_D018.npy')
vel_total_D018 = vel_barents_D018 + vel_bering_D018 + vel_fram_D018 + vel_davis_D018
temp_total_D018 = temp_barents_D018 + temp_bering_D018 + temp_fram_D018 + temp_davis_D018
cov_total_D018 = cov_barents_D018 + cov_bering_D018 + cov_fram_D018 + cov_davis_D018

# Load OHT D019
vel_barents_D019,temp_barents_D019,cov_barents_D019 = np.load(dir_output + dir_D019 + 'decoht_barents_D019.npy')
vel_bering_D019,temp_bering_D019,cov_bering_D019 = np.load(dir_output + dir_D019 + 'decoht_bering_D019.npy')
vel_fram_D019,temp_fram_D019,cov_fram_D019 = np.load(dir_output + dir_D019 + 'decoht_fram_D019.npy')
vel_davis_D019,temp_davis_D019,cov_davis_D019 = np.load(dir_output + dir_D019 + 'decoht_davis_D019.npy')
vel_total_D019 = vel_barents_D019 + vel_bering_D019 + vel_fram_D019 + vel_davis_D019
temp_total_D019 = temp_barents_D019 + temp_bering_D019 + temp_fram_D019 + temp_davis_D019
cov_total_D019 = cov_barents_D019 + cov_bering_D019 + cov_fram_D019 + cov_davis_D019

# Load OHT D020
vel_barents_D020,temp_barents_D020,cov_barents_D020 = np.load(dir_output + dir_D020 + 'decoht_barents_D020.npy')
vel_bering_D020,temp_bering_D020,cov_bering_D020 = np.load(dir_output + dir_D020 + 'decoht_bering_D020.npy')
vel_fram_D020,temp_fram_D020,cov_fram_D020 = np.load(dir_output + dir_D020 + 'decoht_fram_D020.npy')
vel_davis_D020,temp_davis_D020,cov_davis_D020 = np.load(dir_output + dir_D020 + 'decoht_davis_D020.npy')
vel_total_D020 = vel_barents_D020 + vel_bering_D020 + vel_fram_D020 + vel_davis_D020
temp_total_D020 = temp_barents_D020 + temp_bering_D020 + temp_fram_D020 + temp_davis_D020
cov_total_D020 = cov_barents_D020 + cov_bering_D020 + cov_fram_D020 + cov_davis_D020

# Load OHT D021
vel_barents_D021,temp_barents_D021,cov_barents_D021 = np.load(dir_output + dir_D021 + 'decoht_barents_D021.npy')
vel_bering_D021,temp_bering_D021,cov_bering_D021 = np.load(dir_output + dir_D021 + 'decoht_bering_D021.npy')
vel_fram_D021,temp_fram_D021,cov_fram_D021 = np.load(dir_output + dir_D021 + 'decoht_fram_D021.npy')
vel_davis_D021,temp_davis_D021,cov_davis_D021 = np.load(dir_output + dir_D021 + 'decoht_davis_D021.npy')
vel_total_D021 = vel_barents_D021 + vel_bering_D021 + vel_fram_D021 + vel_davis_D021
temp_total_D021 = temp_barents_D021 + temp_bering_D021 + temp_fram_D021 + temp_davis_D021
cov_total_D021 = cov_barents_D021 + cov_bering_D021 + cov_fram_D021 + cov_davis_D021

# Load OHT D022
vel_barents_D022,temp_barents_D022,cov_barents_D022 = np.load(dir_output + dir_D022 + 'decoht_barents_D022.npy')
vel_bering_D022,temp_bering_D022,cov_bering_D022 = np.load(dir_output + dir_D022 + 'decoht_bering_D022.npy')
vel_fram_D022,temp_fram_D022,cov_fram_D022 = np.load(dir_output + dir_D022 + 'decoht_fram_D022.npy')
vel_davis_D022,temp_davis_D022,cov_davis_D022 = np.load(dir_output + dir_D022 + 'decoht_davis_D022.npy')
vel_total_D022 = vel_barents_D022 + vel_bering_D022 + vel_fram_D022 + vel_davis_D022
temp_total_D022 = temp_barents_D022 + temp_bering_D022 + temp_fram_D022 + temp_davis_D022
cov_total_D022 = cov_barents_D022 + cov_bering_D022 + cov_fram_D022 + cov_davis_D022

# Load OHT D023
vel_barents_D023,temp_barents_D023,cov_barents_D023 = np.load(dir_output + dir_D023 + 'decoht_barents_D023.npy')
vel_bering_D023,temp_bering_D023,cov_bering_D023 = np.load(dir_output + dir_D023 + 'decoht_bering_D023.npy')
vel_fram_D023,temp_fram_D023,cov_fram_D023 = np.load(dir_output + dir_D023 + 'decoht_fram_D023.npy')
vel_davis_D023,temp_davis_D023,cov_davis_D023 = np.load(dir_output + dir_D023 + 'decoht_davis_D023.npy')
vel_total_D023 = vel_barents_D023 + vel_bering_D023 + vel_fram_D023 + vel_davis_D023
temp_total_D023 = temp_barents_D023 + temp_bering_D023 + temp_fram_D023 + temp_davis_D023
cov_total_D023 = cov_barents_D023 + cov_bering_D023 + cov_fram_D023 + cov_davis_D023

# Load OHT D024
vel_barents_D024,temp_barents_D024,cov_barents_D024 = np.load(dir_output + dir_D024 + 'decoht_barents_D024.npy')
vel_bering_D024,temp_bering_D024,cov_bering_D024 = np.load(dir_output + dir_D024 + 'decoht_bering_D024.npy')
vel_fram_D024,temp_fram_D024,cov_fram_D024 = np.load(dir_output + dir_D024 + 'decoht_fram_D024.npy')
vel_davis_D024,temp_davis_D024,cov_davis_D024 = np.load(dir_output + dir_D024 + 'decoht_davis_D024.npy')
vel_total_D024 = vel_barents_D024 + vel_bering_D024 + vel_fram_D024 + vel_davis_D024
temp_total_D024 = temp_barents_D024 + temp_bering_D024 + temp_fram_D024 + temp_davis_D024
cov_total_D024 = cov_barents_D024 + cov_bering_D024 + cov_fram_D024 + cov_davis_D024

# Load OHT D025
vel_barents_D025,temp_barents_D025,cov_barents_D025 = np.load(dir_output + dir_D025 + 'decoht_barents_D025.npy')
vel_bering_D025,temp_bering_D025,cov_bering_D025 = np.load(dir_output + dir_D025 + 'decoht_bering_D025.npy')
vel_fram_D025,temp_fram_D025,cov_fram_D025 = np.load(dir_output + dir_D025 + 'decoht_fram_D025.npy')
vel_davis_D025,temp_davis_D025,cov_davis_D025 = np.load(dir_output + dir_D025 + 'decoht_davis_D025.npy')
vel_total_D025 = vel_barents_D025 + vel_bering_D025 + vel_fram_D025 + vel_davis_D025
temp_total_D025 = temp_barents_D025 + temp_bering_D025 + temp_fram_D025 + temp_davis_D025
cov_total_D025 = cov_barents_D025 + cov_bering_D025 + cov_fram_D025 + cov_davis_D025

# Load OHT D027
vel_barents_D027,temp_barents_D027,cov_barents_D027 = np.load(dir_output + dir_D027 + 'decoht_barents_D027.npy')
vel_bering_D027,temp_bering_D027,cov_bering_D027 = np.load(dir_output + dir_D027 + 'decoht_bering_D027.npy')
vel_fram_D027,temp_fram_D027,cov_fram_D027 = np.load(dir_output + dir_D027 + 'decoht_fram_D027.npy')
vel_davis_D027,temp_davis_D027,cov_davis_D027 = np.load(dir_output + dir_D027 + 'decoht_davis_D027.npy')
vel_total_D027 = vel_barents_D027 + vel_bering_D027 + vel_fram_D027 + vel_davis_D027
temp_total_D027 = temp_barents_D027 + temp_bering_D027 + temp_fram_D027 + temp_davis_D027
cov_total_D027 = cov_barents_D027 + cov_bering_D027 + cov_fram_D027 + cov_davis_D027

# Load OHT D028
vel_barents_D028,temp_barents_D028,cov_barents_D028 = np.load(dir_output + dir_D028 + 'decoht_barents_D028.npy')
vel_bering_D028,temp_bering_D028,cov_bering_D028 = np.load(dir_output + dir_D028 + 'decoht_bering_D028.npy')
vel_fram_D028,temp_fram_D028,cov_fram_D028 = np.load(dir_output + dir_D028 + 'decoht_fram_D028.npy')
vel_davis_D028,temp_davis_D028,cov_davis_D028 = np.load(dir_output + dir_D028 + 'decoht_davis_D028.npy')
vel_total_D028 = vel_barents_D028 + vel_bering_D028 + vel_fram_D028 + vel_davis_D028
temp_total_D028 = temp_barents_D028 + temp_bering_D028 + temp_fram_D028 + temp_davis_D028
cov_total_D028 = cov_barents_D028 + cov_bering_D028 + cov_fram_D028 + cov_davis_D028

# Load OHT D029
vel_barents_D029,temp_barents_D029,cov_barents_D029 = np.load(dir_output + dir_D029 + 'decoht_barents_D029.npy')
vel_bering_D029,temp_bering_D029,cov_bering_D029 = np.load(dir_output + dir_D029 + 'decoht_bering_D029.npy')
vel_fram_D029,temp_fram_D029,cov_fram_D029 = np.load(dir_output + dir_D029 + 'decoht_fram_D029.npy')
vel_davis_D029,temp_davis_D029,cov_davis_D029 = np.load(dir_output + dir_D029 + 'decoht_davis_D029.npy')
vel_total_D029 = vel_barents_D029 + vel_bering_D029 + vel_fram_D029 + vel_davis_D029
temp_total_D029 = temp_barents_D029 + temp_bering_D029 + temp_fram_D029 + temp_davis_D029
cov_total_D029 = cov_barents_D029 + cov_bering_D029 + cov_fram_D029 + cov_davis_D029

# Load OHT D030
vel_barents_D030,temp_barents_D030,cov_barents_D030 = np.load(dir_output + dir_D030 + 'decoht_barents_D030.npy')
vel_bering_D030,temp_bering_D030,cov_bering_D030 = np.load(dir_output + dir_D030 + 'decoht_bering_D030.npy')
vel_fram_D030,temp_fram_D030,cov_fram_D030 = np.load(dir_output + dir_D030 + 'decoht_fram_D030.npy')
vel_davis_D030,temp_davis_D030,cov_davis_D030 = np.load(dir_output + dir_D030 + 'decoht_davis_D030.npy')
vel_total_D030 = vel_barents_D030 + vel_bering_D030 + vel_fram_D030 + vel_davis_D030
temp_total_D030 = temp_barents_D030 + temp_bering_D030 + temp_fram_D030 + temp_davis_D030
cov_total_D030 = cov_barents_D030 + cov_bering_D030 + cov_fram_D030 + cov_davis_D030

# Parameters
nmy = 12
nm = np.size(vel_barents_D012)
nyears = int(nm / nmy)

# Compute annual mean OHT D012
vel_annmean_barents_D012 = np.nanmean(vel_barents_D012,axis=1)
vel_annmean_bering_D012 = np.nanmean(vel_bering_D012,axis=1)
vel_annmean_fram_D012 = np.nanmean(vel_fram_D012,axis=1)
vel_annmean_davis_D012 = np.nanmean(vel_davis_D012,axis=1)
vel_annmean_total_D012 = np.nanmean(vel_total_D012,axis=1)
temp_annmean_barents_D012 = np.nanmean(temp_barents_D012,axis=1)
temp_annmean_bering_D012 = np.nanmean(temp_bering_D012,axis=1)
temp_annmean_fram_D012 = np.nanmean(temp_fram_D012,axis=1)
temp_annmean_davis_D012 = np.nanmean(temp_davis_D012,axis=1)
temp_annmean_total_D012 = np.nanmean(temp_total_D012,axis=1)
cov_annmean_barents_D012 = np.nanmean(cov_barents_D012,axis=1)
cov_annmean_bering_D012 = np.nanmean(cov_bering_D012,axis=1)
cov_annmean_fram_D012 = np.nanmean(cov_fram_D012,axis=1)
cov_annmean_davis_D012 = np.nanmean(cov_davis_D012,axis=1)
cov_annmean_total_D012 = np.nanmean(cov_total_D012,axis=1)
oht_annmean_barents_D012 = vel_annmean_barents_D012 + temp_annmean_barents_D012 + cov_annmean_barents_D012
oht_annmean_bering_D012 = vel_annmean_bering_D012 + temp_annmean_bering_D012 + cov_annmean_bering_D012
oht_annmean_fram_D012 = vel_annmean_fram_D012 + temp_annmean_fram_D012 + cov_annmean_fram_D012
oht_annmean_davis_D012 = vel_annmean_davis_D012 + temp_annmean_davis_D012 + cov_annmean_davis_D012
oht_annmean_total_D012 = vel_annmean_total_D012 + temp_annmean_total_D012 + cov_annmean_total_D012

# Compute annual mean OHT D013
vel_annmean_barents_D013 = np.nanmean(vel_barents_D013,axis=1)
vel_annmean_bering_D013 = np.nanmean(vel_bering_D013,axis=1)
vel_annmean_fram_D013 = np.nanmean(vel_fram_D013,axis=1)
vel_annmean_davis_D013 = np.nanmean(vel_davis_D013,axis=1)
vel_annmean_total_D013 = np.nanmean(vel_total_D013,axis=1)
temp_annmean_barents_D013 = np.nanmean(temp_barents_D013,axis=1)
temp_annmean_bering_D013 = np.nanmean(temp_bering_D013,axis=1)
temp_annmean_fram_D013 = np.nanmean(temp_fram_D013,axis=1)
temp_annmean_davis_D013 = np.nanmean(temp_davis_D013,axis=1)
temp_annmean_total_D013 = np.nanmean(temp_total_D013,axis=1)
cov_annmean_barents_D013 = np.nanmean(cov_barents_D013,axis=1)
cov_annmean_bering_D013 = np.nanmean(cov_bering_D013,axis=1)
cov_annmean_fram_D013 = np.nanmean(cov_fram_D013,axis=1)
cov_annmean_davis_D013 = np.nanmean(cov_davis_D013,axis=1)
cov_annmean_total_D013 = np.nanmean(cov_total_D013,axis=1)
oht_annmean_barents_D013 = vel_annmean_barents_D013 + temp_annmean_barents_D013 + cov_annmean_barents_D013
oht_annmean_bering_D013 = vel_annmean_bering_D013 + temp_annmean_bering_D013 + cov_annmean_bering_D013
oht_annmean_fram_D013 = vel_annmean_fram_D013 + temp_annmean_fram_D013 + cov_annmean_fram_D013
oht_annmean_davis_D013 = vel_annmean_davis_D013 + temp_annmean_davis_D013 + cov_annmean_davis_D013
oht_annmean_total_D013 = vel_annmean_total_D013 + temp_annmean_total_D013 + cov_annmean_total_D013

# Compute annual mean OHT D014
vel_annmean_barents_D014 = np.nanmean(vel_barents_D014,axis=1)
vel_annmean_bering_D014 = np.nanmean(vel_bering_D014,axis=1)
vel_annmean_fram_D014 = np.nanmean(vel_fram_D014,axis=1)
vel_annmean_davis_D014 = np.nanmean(vel_davis_D014,axis=1)
vel_annmean_total_D014 = np.nanmean(vel_total_D014,axis=1)
temp_annmean_barents_D014 = np.nanmean(temp_barents_D014,axis=1)
temp_annmean_bering_D014 = np.nanmean(temp_bering_D014,axis=1)
temp_annmean_fram_D014 = np.nanmean(temp_fram_D014,axis=1)
temp_annmean_davis_D014 = np.nanmean(temp_davis_D014,axis=1)
temp_annmean_total_D014 = np.nanmean(temp_total_D014,axis=1)
cov_annmean_barents_D014 = np.nanmean(cov_barents_D014,axis=1)
cov_annmean_bering_D014 = np.nanmean(cov_bering_D014,axis=1)
cov_annmean_fram_D014 = np.nanmean(cov_fram_D014,axis=1)
cov_annmean_davis_D014 = np.nanmean(cov_davis_D014,axis=1)
cov_annmean_total_D014 = np.nanmean(cov_total_D014,axis=1)
oht_annmean_barents_D014 = vel_annmean_barents_D014 + temp_annmean_barents_D014 + cov_annmean_barents_D014
oht_annmean_bering_D014 = vel_annmean_bering_D014 + temp_annmean_bering_D014 + cov_annmean_bering_D014
oht_annmean_fram_D014 = vel_annmean_fram_D014 + temp_annmean_fram_D014 + cov_annmean_fram_D014
oht_annmean_davis_D014 = vel_annmean_davis_D014 + temp_annmean_davis_D014 + cov_annmean_davis_D014
oht_annmean_total_D014 = vel_annmean_total_D014 + temp_annmean_total_D014 + cov_annmean_total_D014

# Compute annual mean OHT D015
vel_annmean_barents_D015 = np.nanmean(vel_barents_D015,axis=1)
vel_annmean_bering_D015 = np.nanmean(vel_bering_D015,axis=1)
vel_annmean_fram_D015 = np.nanmean(vel_fram_D015,axis=1)
vel_annmean_davis_D015 = np.nanmean(vel_davis_D015,axis=1)
vel_annmean_total_D015 = np.nanmean(vel_total_D015,axis=1)
temp_annmean_barents_D015 = np.nanmean(temp_barents_D015,axis=1)
temp_annmean_bering_D015 = np.nanmean(temp_bering_D015,axis=1)
temp_annmean_fram_D015 = np.nanmean(temp_fram_D015,axis=1)
temp_annmean_davis_D015 = np.nanmean(temp_davis_D015,axis=1)
temp_annmean_total_D015 = np.nanmean(temp_total_D015,axis=1)
cov_annmean_barents_D015 = np.nanmean(cov_barents_D015,axis=1)
cov_annmean_bering_D015 = np.nanmean(cov_bering_D015,axis=1)
cov_annmean_fram_D015 = np.nanmean(cov_fram_D015,axis=1)
cov_annmean_davis_D015 = np.nanmean(cov_davis_D015,axis=1)
cov_annmean_total_D015 = np.nanmean(cov_total_D015,axis=1)
oht_annmean_barents_D015 = vel_annmean_barents_D015 + temp_annmean_barents_D015 + cov_annmean_barents_D015
oht_annmean_bering_D015 = vel_annmean_bering_D015 + temp_annmean_bering_D015 + cov_annmean_bering_D015
oht_annmean_fram_D015 = vel_annmean_fram_D015 + temp_annmean_fram_D015 + cov_annmean_fram_D015
oht_annmean_davis_D015 = vel_annmean_davis_D015 + temp_annmean_davis_D015 + cov_annmean_davis_D015
oht_annmean_total_D015 = vel_annmean_total_D015 + temp_annmean_total_D015 + cov_annmean_total_D015

# Compute annual mean OHT D016
vel_annmean_barents_D016 = np.nanmean(vel_barents_D016,axis=1)
vel_annmean_bering_D016 = np.nanmean(vel_bering_D016,axis=1)
vel_annmean_fram_D016 = np.nanmean(vel_fram_D016,axis=1)
vel_annmean_davis_D016 = np.nanmean(vel_davis_D016,axis=1)
vel_annmean_total_D016 = np.nanmean(vel_total_D016,axis=1)
temp_annmean_barents_D016 = np.nanmean(temp_barents_D016,axis=1)
temp_annmean_bering_D016 = np.nanmean(temp_bering_D016,axis=1)
temp_annmean_fram_D016 = np.nanmean(temp_fram_D016,axis=1)
temp_annmean_davis_D016 = np.nanmean(temp_davis_D016,axis=1)
temp_annmean_total_D016 = np.nanmean(temp_total_D016,axis=1)
cov_annmean_barents_D016 = np.nanmean(cov_barents_D016,axis=1)
cov_annmean_bering_D016 = np.nanmean(cov_bering_D016,axis=1)
cov_annmean_fram_D016 = np.nanmean(cov_fram_D016,axis=1)
cov_annmean_davis_D016 = np.nanmean(cov_davis_D016,axis=1)
cov_annmean_total_D016 = np.nanmean(cov_total_D016,axis=1)
oht_annmean_barents_D016 = vel_annmean_barents_D016 + temp_annmean_barents_D016 + cov_annmean_barents_D016
oht_annmean_bering_D016 = vel_annmean_bering_D016 + temp_annmean_bering_D016 + cov_annmean_bering_D016
oht_annmean_fram_D016 = vel_annmean_fram_D016 + temp_annmean_fram_D016 + cov_annmean_fram_D016
oht_annmean_davis_D016 = vel_annmean_davis_D016 + temp_annmean_davis_D016 + cov_annmean_davis_D016
oht_annmean_total_D016 = vel_annmean_total_D016 + temp_annmean_total_D016 + cov_annmean_total_D016

# Compute annual mean OHT D017
vel_annmean_barents_D017 = np.nanmean(vel_barents_D017,axis=1)
vel_annmean_bering_D017 = np.nanmean(vel_bering_D017,axis=1)
vel_annmean_fram_D017 = np.nanmean(vel_fram_D017,axis=1)
vel_annmean_davis_D017 = np.nanmean(vel_davis_D017,axis=1)
vel_annmean_total_D017 = np.nanmean(vel_total_D017,axis=1)
temp_annmean_barents_D017 = np.nanmean(temp_barents_D017,axis=1)
temp_annmean_bering_D017 = np.nanmean(temp_bering_D017,axis=1)
temp_annmean_fram_D017 = np.nanmean(temp_fram_D017,axis=1)
temp_annmean_davis_D017 = np.nanmean(temp_davis_D017,axis=1)
temp_annmean_total_D017 = np.nanmean(temp_total_D017,axis=1)
cov_annmean_barents_D017 = np.nanmean(cov_barents_D017,axis=1)
cov_annmean_bering_D017 = np.nanmean(cov_bering_D017,axis=1)
cov_annmean_fram_D017 = np.nanmean(cov_fram_D017,axis=1)
cov_annmean_davis_D017 = np.nanmean(cov_davis_D017,axis=1)
cov_annmean_total_D017 = np.nanmean(cov_total_D017,axis=1)
oht_annmean_barents_D017 = vel_annmean_barents_D017 + temp_annmean_barents_D017 + cov_annmean_barents_D017
oht_annmean_bering_D017 = vel_annmean_bering_D017 + temp_annmean_bering_D017 + cov_annmean_bering_D017
oht_annmean_fram_D017 = vel_annmean_fram_D017 + temp_annmean_fram_D017 + cov_annmean_fram_D017
oht_annmean_davis_D017 = vel_annmean_davis_D017 + temp_annmean_davis_D017 + cov_annmean_davis_D017
oht_annmean_total_D017 = vel_annmean_total_D017 + temp_annmean_total_D017 + cov_annmean_total_D017

# Compute annual mean OHT D018
vel_annmean_barents_D018 = np.nanmean(vel_barents_D018,axis=1)
vel_annmean_bering_D018 = np.nanmean(vel_bering_D018,axis=1)
vel_annmean_fram_D018 = np.nanmean(vel_fram_D018,axis=1)
vel_annmean_davis_D018 = np.nanmean(vel_davis_D018,axis=1)
vel_annmean_total_D018 = np.nanmean(vel_total_D018,axis=1)
temp_annmean_barents_D018 = np.nanmean(temp_barents_D018,axis=1)
temp_annmean_bering_D018 = np.nanmean(temp_bering_D018,axis=1)
temp_annmean_fram_D018 = np.nanmean(temp_fram_D018,axis=1)
temp_annmean_davis_D018 = np.nanmean(temp_davis_D018,axis=1)
temp_annmean_total_D018 = np.nanmean(temp_total_D018,axis=1)
cov_annmean_barents_D018 = np.nanmean(cov_barents_D018,axis=1)
cov_annmean_bering_D018 = np.nanmean(cov_bering_D018,axis=1)
cov_annmean_fram_D018 = np.nanmean(cov_fram_D018,axis=1)
cov_annmean_davis_D018 = np.nanmean(cov_davis_D018,axis=1)
cov_annmean_total_D018 = np.nanmean(cov_total_D018,axis=1)
oht_annmean_barents_D018 = vel_annmean_barents_D018 + temp_annmean_barents_D018 + cov_annmean_barents_D018
oht_annmean_bering_D018 = vel_annmean_bering_D018 + temp_annmean_bering_D018 + cov_annmean_bering_D018
oht_annmean_fram_D018 = vel_annmean_fram_D018 + temp_annmean_fram_D018 + cov_annmean_fram_D018
oht_annmean_davis_D018 = vel_annmean_davis_D018 + temp_annmean_davis_D018 + cov_annmean_davis_D018
oht_annmean_total_D018 = vel_annmean_total_D018 + temp_annmean_total_D018 + cov_annmean_total_D018

# Compute annual mean OHT D019
vel_annmean_barents_D019 = np.nanmean(vel_barents_D019,axis=1)
vel_annmean_bering_D019 = np.nanmean(vel_bering_D019,axis=1)
vel_annmean_fram_D019 = np.nanmean(vel_fram_D019,axis=1)
vel_annmean_davis_D019 = np.nanmean(vel_davis_D019,axis=1)
vel_annmean_total_D019 = np.nanmean(vel_total_D019,axis=1)
temp_annmean_barents_D019 = np.nanmean(temp_barents_D019,axis=1)
temp_annmean_bering_D019 = np.nanmean(temp_bering_D019,axis=1)
temp_annmean_fram_D019 = np.nanmean(temp_fram_D019,axis=1)
temp_annmean_davis_D019 = np.nanmean(temp_davis_D019,axis=1)
temp_annmean_total_D019 = np.nanmean(temp_total_D019,axis=1)
cov_annmean_barents_D019 = np.nanmean(cov_barents_D019,axis=1)
cov_annmean_bering_D019 = np.nanmean(cov_bering_D019,axis=1)
cov_annmean_fram_D019 = np.nanmean(cov_fram_D019,axis=1)
cov_annmean_davis_D019 = np.nanmean(cov_davis_D019,axis=1)
cov_annmean_total_D019 = np.nanmean(cov_total_D019,axis=1)
oht_annmean_barents_D019 = vel_annmean_barents_D019 + temp_annmean_barents_D019 + cov_annmean_barents_D019
oht_annmean_bering_D019 = vel_annmean_bering_D019 + temp_annmean_bering_D019 + cov_annmean_bering_D019
oht_annmean_fram_D019 = vel_annmean_fram_D019 + temp_annmean_fram_D019 + cov_annmean_fram_D019
oht_annmean_davis_D019 = vel_annmean_davis_D019 + temp_annmean_davis_D019 + cov_annmean_davis_D019
oht_annmean_total_D019 = vel_annmean_total_D019 + temp_annmean_total_D019 + cov_annmean_total_D019

# Compute annual mean OHT D020
vel_annmean_barents_D020 = np.nanmean(vel_barents_D020,axis=1)
vel_annmean_bering_D020 = np.nanmean(vel_bering_D020,axis=1)
vel_annmean_fram_D020 = np.nanmean(vel_fram_D020,axis=1)
vel_annmean_davis_D020 = np.nanmean(vel_davis_D020,axis=1)
vel_annmean_total_D020 = np.nanmean(vel_total_D020,axis=1)
temp_annmean_barents_D020 = np.nanmean(temp_barents_D020,axis=1)
temp_annmean_bering_D020 = np.nanmean(temp_bering_D020,axis=1)
temp_annmean_fram_D020 = np.nanmean(temp_fram_D020,axis=1)
temp_annmean_davis_D020 = np.nanmean(temp_davis_D020,axis=1)
temp_annmean_total_D020 = np.nanmean(temp_total_D020,axis=1)
cov_annmean_barents_D020 = np.nanmean(cov_barents_D020,axis=1)
cov_annmean_bering_D020 = np.nanmean(cov_bering_D020,axis=1)
cov_annmean_fram_D020 = np.nanmean(cov_fram_D020,axis=1)
cov_annmean_davis_D020 = np.nanmean(cov_davis_D020,axis=1)
cov_annmean_total_D020 = np.nanmean(cov_total_D020,axis=1)
oht_annmean_barents_D020 = vel_annmean_barents_D020 + temp_annmean_barents_D020 + cov_annmean_barents_D020
oht_annmean_bering_D020 = vel_annmean_bering_D020 + temp_annmean_bering_D020 + cov_annmean_bering_D020
oht_annmean_fram_D020 = vel_annmean_fram_D020 + temp_annmean_fram_D020 + cov_annmean_fram_D020
oht_annmean_davis_D020 = vel_annmean_davis_D020 + temp_annmean_davis_D020 + cov_annmean_davis_D020
oht_annmean_total_D020 = vel_annmean_total_D020 + temp_annmean_total_D020 + cov_annmean_total_D020

# Compute annual mean OHT D021
vel_annmean_barents_D021 = np.nanmean(vel_barents_D021,axis=1)
vel_annmean_bering_D021 = np.nanmean(vel_bering_D021,axis=1)
vel_annmean_fram_D021 = np.nanmean(vel_fram_D021,axis=1)
vel_annmean_davis_D021 = np.nanmean(vel_davis_D021,axis=1)
vel_annmean_total_D021 = np.nanmean(vel_total_D021,axis=1)
temp_annmean_barents_D021 = np.nanmean(temp_barents_D021,axis=1)
temp_annmean_bering_D021 = np.nanmean(temp_bering_D021,axis=1)
temp_annmean_fram_D021 = np.nanmean(temp_fram_D021,axis=1)
temp_annmean_davis_D021 = np.nanmean(temp_davis_D021,axis=1)
temp_annmean_total_D021 = np.nanmean(temp_total_D021,axis=1)
cov_annmean_barents_D021 = np.nanmean(cov_barents_D021,axis=1)
cov_annmean_bering_D021 = np.nanmean(cov_bering_D021,axis=1)
cov_annmean_fram_D021 = np.nanmean(cov_fram_D021,axis=1)
cov_annmean_davis_D021 = np.nanmean(cov_davis_D021,axis=1)
cov_annmean_total_D021 = np.nanmean(cov_total_D021,axis=1)
oht_annmean_barents_D021 = vel_annmean_barents_D021 + temp_annmean_barents_D021 + cov_annmean_barents_D021
oht_annmean_bering_D021 = vel_annmean_bering_D021 + temp_annmean_bering_D021 + cov_annmean_bering_D021
oht_annmean_fram_D021 = vel_annmean_fram_D021 + temp_annmean_fram_D021 + cov_annmean_fram_D021
oht_annmean_davis_D021 = vel_annmean_davis_D021 + temp_annmean_davis_D021 + cov_annmean_davis_D021
oht_annmean_total_D021 = vel_annmean_total_D021 + temp_annmean_total_D021 + cov_annmean_total_D021

# Compute annual mean OHT D022
vel_annmean_barents_D022 = np.nanmean(vel_barents_D022,axis=1)
vel_annmean_bering_D022 = np.nanmean(vel_bering_D022,axis=1)
vel_annmean_fram_D022 = np.nanmean(vel_fram_D022,axis=1)
vel_annmean_davis_D022 = np.nanmean(vel_davis_D022,axis=1)
vel_annmean_total_D022 = np.nanmean(vel_total_D022,axis=1)
temp_annmean_barents_D022 = np.nanmean(temp_barents_D022,axis=1)
temp_annmean_bering_D022 = np.nanmean(temp_bering_D022,axis=1)
temp_annmean_fram_D022 = np.nanmean(temp_fram_D022,axis=1)
temp_annmean_davis_D022 = np.nanmean(temp_davis_D022,axis=1)
temp_annmean_total_D022 = np.nanmean(temp_total_D022,axis=1)
cov_annmean_barents_D022 = np.nanmean(cov_barents_D022,axis=1)
cov_annmean_bering_D022 = np.nanmean(cov_bering_D022,axis=1)
cov_annmean_fram_D022 = np.nanmean(cov_fram_D022,axis=1)
cov_annmean_davis_D022 = np.nanmean(cov_davis_D022,axis=1)
cov_annmean_total_D022 = np.nanmean(cov_total_D022,axis=1)
oht_annmean_barents_D022 = vel_annmean_barents_D022 + temp_annmean_barents_D022 + cov_annmean_barents_D022
oht_annmean_bering_D022 = vel_annmean_bering_D022 + temp_annmean_bering_D022 + cov_annmean_bering_D022
oht_annmean_fram_D022 = vel_annmean_fram_D022 + temp_annmean_fram_D022 + cov_annmean_fram_D022
oht_annmean_davis_D022 = vel_annmean_davis_D022 + temp_annmean_davis_D022 + cov_annmean_davis_D022
oht_annmean_total_D022 = vel_annmean_total_D022 + temp_annmean_total_D022 + cov_annmean_total_D022

# Compute annual mean OHT D023
vel_annmean_barents_D023 = np.nanmean(vel_barents_D023,axis=1)
vel_annmean_bering_D023 = np.nanmean(vel_bering_D023,axis=1)
vel_annmean_fram_D023 = np.nanmean(vel_fram_D023,axis=1)
vel_annmean_davis_D023 = np.nanmean(vel_davis_D023,axis=1)
vel_annmean_total_D023 = np.nanmean(vel_total_D023,axis=1)
temp_annmean_barents_D023 = np.nanmean(temp_barents_D023,axis=1)
temp_annmean_bering_D023 = np.nanmean(temp_bering_D023,axis=1)
temp_annmean_fram_D023 = np.nanmean(temp_fram_D023,axis=1)
temp_annmean_davis_D023 = np.nanmean(temp_davis_D023,axis=1)
temp_annmean_total_D023 = np.nanmean(temp_total_D023,axis=1)
cov_annmean_barents_D023 = np.nanmean(cov_barents_D023,axis=1)
cov_annmean_bering_D023 = np.nanmean(cov_bering_D023,axis=1)
cov_annmean_fram_D023 = np.nanmean(cov_fram_D023,axis=1)
cov_annmean_davis_D023 = np.nanmean(cov_davis_D023,axis=1)
cov_annmean_total_D023 = np.nanmean(cov_total_D023,axis=1)
oht_annmean_barents_D023 = vel_annmean_barents_D023 + temp_annmean_barents_D023 + cov_annmean_barents_D023
oht_annmean_bering_D023 = vel_annmean_bering_D023 + temp_annmean_bering_D023 + cov_annmean_bering_D023
oht_annmean_fram_D023 = vel_annmean_fram_D023 + temp_annmean_fram_D023 + cov_annmean_fram_D023
oht_annmean_davis_D023 = vel_annmean_davis_D023 + temp_annmean_davis_D023 + cov_annmean_davis_D023
oht_annmean_total_D023 = vel_annmean_total_D023 + temp_annmean_total_D023 + cov_annmean_total_D023

# Compute annual mean OHT D024
vel_annmean_barents_D024 = np.nanmean(vel_barents_D024,axis=1)
vel_annmean_bering_D024 = np.nanmean(vel_bering_D024,axis=1)
vel_annmean_fram_D024 = np.nanmean(vel_fram_D024,axis=1)
vel_annmean_davis_D024 = np.nanmean(vel_davis_D024,axis=1)
vel_annmean_total_D024 = np.nanmean(vel_total_D024,axis=1)
temp_annmean_barents_D024 = np.nanmean(temp_barents_D024,axis=1)
temp_annmean_bering_D024 = np.nanmean(temp_bering_D024,axis=1)
temp_annmean_fram_D024 = np.nanmean(temp_fram_D024,axis=1)
temp_annmean_davis_D024 = np.nanmean(temp_davis_D024,axis=1)
temp_annmean_total_D024 = np.nanmean(temp_total_D024,axis=1)
cov_annmean_barents_D024 = np.nanmean(cov_barents_D024,axis=1)
cov_annmean_bering_D024 = np.nanmean(cov_bering_D024,axis=1)
cov_annmean_fram_D024 = np.nanmean(cov_fram_D024,axis=1)
cov_annmean_davis_D024 = np.nanmean(cov_davis_D024,axis=1)
cov_annmean_total_D024 = np.nanmean(cov_total_D024,axis=1)
oht_annmean_barents_D024 = vel_annmean_barents_D024 + temp_annmean_barents_D024 + cov_annmean_barents_D024
oht_annmean_bering_D024 = vel_annmean_bering_D024 + temp_annmean_bering_D024 + cov_annmean_bering_D024
oht_annmean_fram_D024 = vel_annmean_fram_D024 + temp_annmean_fram_D024 + cov_annmean_fram_D024
oht_annmean_davis_D024 = vel_annmean_davis_D024 + temp_annmean_davis_D024 + cov_annmean_davis_D024
oht_annmean_total_D024 = vel_annmean_total_D024 + temp_annmean_total_D024 + cov_annmean_total_D024

# Compute annual mean OHT D025
vel_annmean_barents_D025 = np.nanmean(vel_barents_D025,axis=1)
vel_annmean_bering_D025 = np.nanmean(vel_bering_D025,axis=1)
vel_annmean_fram_D025 = np.nanmean(vel_fram_D025,axis=1)
vel_annmean_davis_D025 = np.nanmean(vel_davis_D025,axis=1)
vel_annmean_total_D025 = np.nanmean(vel_total_D025,axis=1)
temp_annmean_barents_D025 = np.nanmean(temp_barents_D025,axis=1)
temp_annmean_bering_D025 = np.nanmean(temp_bering_D025,axis=1)
temp_annmean_fram_D025 = np.nanmean(temp_fram_D025,axis=1)
temp_annmean_davis_D025 = np.nanmean(temp_davis_D025,axis=1)
temp_annmean_total_D025 = np.nanmean(temp_total_D025,axis=1)
cov_annmean_barents_D025 = np.nanmean(cov_barents_D025,axis=1)
cov_annmean_bering_D025 = np.nanmean(cov_bering_D025,axis=1)
cov_annmean_fram_D025 = np.nanmean(cov_fram_D025,axis=1)
cov_annmean_davis_D025 = np.nanmean(cov_davis_D025,axis=1)
cov_annmean_total_D025 = np.nanmean(cov_total_D025,axis=1)
oht_annmean_barents_D025 = vel_annmean_barents_D025 + temp_annmean_barents_D025 + cov_annmean_barents_D025
oht_annmean_bering_D025 = vel_annmean_bering_D025 + temp_annmean_bering_D025 + cov_annmean_bering_D025
oht_annmean_fram_D025 = vel_annmean_fram_D025 + temp_annmean_fram_D025 + cov_annmean_fram_D025
oht_annmean_davis_D025 = vel_annmean_davis_D025 + temp_annmean_davis_D025 + cov_annmean_davis_D025
oht_annmean_total_D025 = vel_annmean_total_D025 + temp_annmean_total_D025 + cov_annmean_total_D025

# Compute annual mean OHT D027
vel_annmean_barents_D027 = np.nanmean(vel_barents_D027,axis=1)
vel_annmean_bering_D027 = np.nanmean(vel_bering_D027,axis=1)
vel_annmean_fram_D027 = np.nanmean(vel_fram_D027,axis=1)
vel_annmean_davis_D027 = np.nanmean(vel_davis_D027,axis=1)
vel_annmean_total_D027 = np.nanmean(vel_total_D027,axis=1)
temp_annmean_barents_D027 = np.nanmean(temp_barents_D027,axis=1)
temp_annmean_bering_D027 = np.nanmean(temp_bering_D027,axis=1)
temp_annmean_fram_D027 = np.nanmean(temp_fram_D027,axis=1)
temp_annmean_davis_D027 = np.nanmean(temp_davis_D027,axis=1)
temp_annmean_total_D027 = np.nanmean(temp_total_D027,axis=1)
cov_annmean_barents_D027 = np.nanmean(cov_barents_D027,axis=1)
cov_annmean_bering_D027 = np.nanmean(cov_bering_D027,axis=1)
cov_annmean_fram_D027 = np.nanmean(cov_fram_D027,axis=1)
cov_annmean_davis_D027 = np.nanmean(cov_davis_D027,axis=1)
cov_annmean_total_D027 = np.nanmean(cov_total_D027,axis=1)
oht_annmean_barents_D027 = vel_annmean_barents_D027 + temp_annmean_barents_D027 + cov_annmean_barents_D027
oht_annmean_bering_D027 = vel_annmean_bering_D027 + temp_annmean_bering_D027 + cov_annmean_bering_D027
oht_annmean_fram_D027 = vel_annmean_fram_D027 + temp_annmean_fram_D027 + cov_annmean_fram_D027
oht_annmean_davis_D027 = vel_annmean_davis_D027 + temp_annmean_davis_D027 + cov_annmean_davis_D027
oht_annmean_total_D027 = vel_annmean_total_D027 + temp_annmean_total_D027 + cov_annmean_total_D027

# Compute annual mean OHT D028
vel_annmean_barents_D028 = np.nanmean(vel_barents_D028,axis=1)
vel_annmean_bering_D028 = np.nanmean(vel_bering_D028,axis=1)
vel_annmean_fram_D028 = np.nanmean(vel_fram_D028,axis=1)
vel_annmean_davis_D028 = np.nanmean(vel_davis_D028,axis=1)
vel_annmean_total_D028 = np.nanmean(vel_total_D028,axis=1)
temp_annmean_barents_D028 = np.nanmean(temp_barents_D028,axis=1)
temp_annmean_bering_D028 = np.nanmean(temp_bering_D028,axis=1)
temp_annmean_fram_D028 = np.nanmean(temp_fram_D028,axis=1)
temp_annmean_davis_D028 = np.nanmean(temp_davis_D028,axis=1)
temp_annmean_total_D028 = np.nanmean(temp_total_D028,axis=1)
cov_annmean_barents_D028 = np.nanmean(cov_barents_D028,axis=1)
cov_annmean_bering_D028 = np.nanmean(cov_bering_D028,axis=1)
cov_annmean_fram_D028 = np.nanmean(cov_fram_D028,axis=1)
cov_annmean_davis_D028 = np.nanmean(cov_davis_D028,axis=1)
cov_annmean_total_D028 = np.nanmean(cov_total_D028,axis=1)
oht_annmean_barents_D028 = vel_annmean_barents_D028 + temp_annmean_barents_D028 + cov_annmean_barents_D028
oht_annmean_bering_D028 = vel_annmean_bering_D028 + temp_annmean_bering_D028 + cov_annmean_bering_D028
oht_annmean_fram_D028 = vel_annmean_fram_D028 + temp_annmean_fram_D028 + cov_annmean_fram_D028
oht_annmean_davis_D028 = vel_annmean_davis_D028 + temp_annmean_davis_D028 + cov_annmean_davis_D028
oht_annmean_total_D028 = vel_annmean_total_D028 + temp_annmean_total_D028 + cov_annmean_total_D028

# Compute annual mean OHT D029
vel_annmean_barents_D029 = np.nanmean(vel_barents_D029,axis=1)
vel_annmean_bering_D029 = np.nanmean(vel_bering_D029,axis=1)
vel_annmean_fram_D029 = np.nanmean(vel_fram_D029,axis=1)
vel_annmean_davis_D029 = np.nanmean(vel_davis_D029,axis=1)
vel_annmean_total_D029 = np.nanmean(vel_total_D029,axis=1)
temp_annmean_barents_D029 = np.nanmean(temp_barents_D029,axis=1)
temp_annmean_bering_D029 = np.nanmean(temp_bering_D029,axis=1)
temp_annmean_fram_D029 = np.nanmean(temp_fram_D029,axis=1)
temp_annmean_davis_D029 = np.nanmean(temp_davis_D029,axis=1)
temp_annmean_total_D029 = np.nanmean(temp_total_D029,axis=1)
cov_annmean_barents_D029 = np.nanmean(cov_barents_D029,axis=1)
cov_annmean_bering_D029 = np.nanmean(cov_bering_D029,axis=1)
cov_annmean_fram_D029 = np.nanmean(cov_fram_D029,axis=1)
cov_annmean_davis_D029 = np.nanmean(cov_davis_D029,axis=1)
cov_annmean_total_D029 = np.nanmean(cov_total_D029,axis=1)
oht_annmean_barents_D029 = vel_annmean_barents_D029 + temp_annmean_barents_D029 + cov_annmean_barents_D029
oht_annmean_bering_D029 = vel_annmean_bering_D029 + temp_annmean_bering_D029 + cov_annmean_bering_D029
oht_annmean_fram_D029 = vel_annmean_fram_D029 + temp_annmean_fram_D029 + cov_annmean_fram_D029
oht_annmean_davis_D029 = vel_annmean_davis_D029 + temp_annmean_davis_D029 + cov_annmean_davis_D029
oht_annmean_total_D029 = vel_annmean_total_D029 + temp_annmean_total_D029 + cov_annmean_total_D029

# Compute annual mean OHT D030
vel_annmean_barents_D030 = np.nanmean(vel_barents_D030,axis=1)
vel_annmean_bering_D030 = np.nanmean(vel_bering_D030,axis=1)
vel_annmean_fram_D030 = np.nanmean(vel_fram_D030,axis=1)
vel_annmean_davis_D030 = np.nanmean(vel_davis_D030,axis=1)
vel_annmean_total_D030 = np.nanmean(vel_total_D030,axis=1)
temp_annmean_barents_D030 = np.nanmean(temp_barents_D030,axis=1)
temp_annmean_bering_D030 = np.nanmean(temp_bering_D030,axis=1)
temp_annmean_fram_D030 = np.nanmean(temp_fram_D030,axis=1)
temp_annmean_davis_D030 = np.nanmean(temp_davis_D030,axis=1)
temp_annmean_total_D030 = np.nanmean(temp_total_D030,axis=1)
cov_annmean_barents_D030 = np.nanmean(cov_barents_D030,axis=1)
cov_annmean_bering_D030 = np.nanmean(cov_bering_D030,axis=1)
cov_annmean_fram_D030 = np.nanmean(cov_fram_D030,axis=1)
cov_annmean_davis_D030 = np.nanmean(cov_davis_D030,axis=1)
cov_annmean_total_D030 = np.nanmean(cov_total_D030,axis=1)
oht_annmean_barents_D030 = vel_annmean_barents_D030 + temp_annmean_barents_D030 + cov_annmean_barents_D030
oht_annmean_bering_D030 = vel_annmean_bering_D030 + temp_annmean_bering_D030 + cov_annmean_bering_D030
oht_annmean_fram_D030 = vel_annmean_fram_D030 + temp_annmean_fram_D030 + cov_annmean_fram_D030
oht_annmean_davis_D030 = vel_annmean_davis_D030 + temp_annmean_davis_D030 + cov_annmean_davis_D030
oht_annmean_total_D030 = vel_annmean_total_D030 + temp_annmean_total_D030 + cov_annmean_total_D030

# Labels
name_xticks = ['2130','2140','2150','2160','2170','2180']
    

# Fig. 3 - Bar plots of changes in OHT (decomposition)
fig,ax = plt.subplots(2,2,figsize=(18,12))
fig.subplots_adjust(left=0.1,bottom=0.2,right=0.95,top=0.95,wspace=None,hspace=0.7)
index = np.arange(18)
bar_width = 1
name_xticks = ['ATL1+1K','ATL1+3K','ATL1+5K','ATL2+1K','ATL2+3K','ATL2+5K','ATL3+1K','ATL3+3K','ATL3+5K','PAC1+1K','PAC1+3K','PAC1+5K','PAC2+1K','PAC2+3K','PAC2+5K','PAC3+1K','PAC3+3K','PAC3+5K']

# Total OHT changes
ax[0,0].set_title('Total OHT changes $v_sT_s - v_cT_c$',fontsize=26)
ax[0,0].set_title('a',loc='left',fontsize=26,fontweight='bold')
ax[0,0].bar(index[0],np.nanmean(oht_annmean_total_D013),bar_width,color='lightcoral')
ax[0,0].bar(index[1],np.nanmean(oht_annmean_total_D012),bar_width,color='red')
ax[0,0].bar(index[2],np.nanmean(oht_annmean_total_D014),bar_width,color='purple')
ax[0,0].bar(index[3],np.nanmean(oht_annmean_total_D016),bar_width,color='lightblue')
ax[0,0].bar(index[4],np.nanmean(oht_annmean_total_D015),bar_width,color='blue')
ax[0,0].bar(index[5],np.nanmean(oht_annmean_total_D017),bar_width,color='darkblue')
ax[0,0].bar(index[6],np.nanmean(oht_annmean_total_D019),bar_width,color='lightgreen')
ax[0,0].bar(index[7],np.nanmean(oht_annmean_total_D018),bar_width,color='green')
ax[0,0].bar(index[8],np.nanmean(oht_annmean_total_D020),bar_width,color='darkgreen')
ax[0,0].bar(index[9],np.nanmean(oht_annmean_total_D027),bar_width,color='lightcyan')
ax[0,0].bar(index[10],np.nanmean(oht_annmean_total_D021),bar_width,color='cyan')
ax[0,0].bar(index[11],np.nanmean(oht_annmean_total_D028),bar_width,color='darkcyan')
ax[0,0].bar(index[12],np.nanmean(oht_annmean_total_D029),bar_width,color='orange')
ax[0,0].bar(index[13],np.nanmean(oht_annmean_total_D022),bar_width,color='goldenrod')
ax[0,0].bar(index[14],np.nanmean(oht_annmean_total_D030),bar_width,color='brown')
ax[0,0].bar(index[15],np.nanmean(oht_annmean_total_D024),bar_width,color='lightgray')
ax[0,0].bar(index[16],np.nanmean(oht_annmean_total_D023),bar_width,color='gray')
ax[0,0].bar(index[17],np.nanmean(oht_annmean_total_D025),bar_width,color='black')
ax[0,0].set_ylabel('Arctic OHT PERT-CTRL \n (TW)',fontsize=24)
ax[0,0].set_yticks(np.arange(-20, 141, 20))
ax[0,0].set_xticks(index)
ax[0,0].set_xticklabels(name_xticks,rotation=90)
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].axhline(c='k')
ax[0,0].set_xlabel('Experiment',fontsize=24)

# OHT changes due to velocity
ax[0,1].set_title('OHT changes due to velocity $v^{\prime}T_c$',fontsize=24)
ax[0,1].set_title('b',loc='left',fontsize=26,fontweight='bold')
ax[0,1].bar(index[0],np.nanmean(vel_annmean_total_D013),bar_width,color='lightcoral')
ax[0,1].bar(index[1],np.nanmean(vel_annmean_total_D012),bar_width,color='red')
ax[0,1].bar(index[2],np.nanmean(vel_annmean_total_D014),bar_width,color='purple')
ax[0,1].bar(index[3],np.nanmean(vel_annmean_total_D016),bar_width,color='lightblue')
ax[0,1].bar(index[4],np.nanmean(vel_annmean_total_D015),bar_width,color='blue')
ax[0,1].bar(index[5],np.nanmean(vel_annmean_total_D017),bar_width,color='darkblue')
ax[0,1].bar(index[6],np.nanmean(vel_annmean_total_D019),bar_width,color='lightgreen')
ax[0,1].bar(index[7],np.nanmean(vel_annmean_total_D018),bar_width,color='green')
ax[0,1].bar(index[8],np.nanmean(vel_annmean_total_D020),bar_width,color='darkgreen')
ax[0,1].bar(index[9],np.nanmean(vel_annmean_total_D027),bar_width,color='lightcyan')
ax[0,1].bar(index[10],np.nanmean(vel_annmean_total_D021),bar_width,color='cyan')
ax[0,1].bar(index[11],np.nanmean(vel_annmean_total_D028),bar_width,color='darkcyan')
ax[0,1].bar(index[12],np.nanmean(vel_annmean_total_D029),bar_width,color='orange')
ax[0,1].bar(index[13],np.nanmean(vel_annmean_total_D022),bar_width,color='goldenrod')
ax[0,1].bar(index[14],np.nanmean(vel_annmean_total_D030),bar_width,color='brown')
ax[0,1].bar(index[15],np.nanmean(vel_annmean_total_D024),bar_width,color='lightgray')
ax[0,1].bar(index[16],np.nanmean(vel_annmean_total_D023),bar_width,color='gray')
ax[0,1].bar(index[17],np.nanmean(vel_annmean_total_D025),bar_width,color='black')
ax[0,1].set_yticks(np.arange(-20, 141, 20))
ax[0,1].set_xticks(index)
ax[0,1].set_xticklabels(name_xticks,rotation=90)
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].axhline(c='k')
ax[0,1].set_xlabel('Experiment',fontsize=24)

# OHT changes due to temperature
ax[1,0].set_title('OHT changes due to temperature $v_cT^{\prime}$',fontsize=24)
ax[1,0].set_title('c',loc='left',fontsize=26,fontweight='bold')
ax[1,0].bar(index[0],np.nanmean(temp_annmean_total_D013),bar_width,color='lightcoral')
ax[1,0].bar(index[1],np.nanmean(temp_annmean_total_D012),bar_width,color='red')
ax[1,0].bar(index[2],np.nanmean(temp_annmean_total_D014),bar_width,color='purple')
ax[1,0].bar(index[3],np.nanmean(temp_annmean_total_D016),bar_width,color='lightblue')
ax[1,0].bar(index[4],np.nanmean(temp_annmean_total_D015),bar_width,color='blue')
ax[1,0].bar(index[5],np.nanmean(temp_annmean_total_D017),bar_width,color='darkblue')
ax[1,0].bar(index[6],np.nanmean(temp_annmean_total_D019),bar_width,color='lightgreen')
ax[1,0].bar(index[7],np.nanmean(temp_annmean_total_D018),bar_width,color='green')
ax[1,0].bar(index[8],np.nanmean(temp_annmean_total_D020),bar_width,color='darkgreen')
ax[1,0].bar(index[9],np.nanmean(temp_annmean_total_D027),bar_width,color='lightcyan')
ax[1,0].bar(index[10],np.nanmean(temp_annmean_total_D021),bar_width,color='cyan')
ax[1,0].bar(index[11],np.nanmean(temp_annmean_total_D028),bar_width,color='darkcyan')
ax[1,0].bar(index[12],np.nanmean(temp_annmean_total_D029),bar_width,color='orange')
ax[1,0].bar(index[13],np.nanmean(temp_annmean_total_D022),bar_width,color='goldenrod')
ax[1,0].bar(index[14],np.nanmean(temp_annmean_total_D030),bar_width,color='brown')
ax[1,0].bar(index[15],np.nanmean(temp_annmean_total_D024),bar_width,color='lightgray')
ax[1,0].bar(index[16],np.nanmean(temp_annmean_total_D023),bar_width,color='gray')
ax[1,0].bar(index[17],np.nanmean(temp_annmean_total_D025),bar_width,color='black')
ax[1,0].set_ylabel('Arctic OHT PERT-CTRL \n (TW)',fontsize=24)
ax[1,0].set_yticks(np.arange(-20, 141, 20))
ax[1,0].set_xlabel('Experiment',fontsize=24)
ax[1,0].set_xticks(index)
ax[1,0].set_xticklabels(name_xticks,rotation=90)
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].axhline(c='k')

# OHT changes due to covariance
ax[1,1].set_title('OHT changes due to covariance $v^{\prime}T^{\prime}$',fontsize=24)
ax[1,1].set_title('d',loc='left',fontsize=26,fontweight='bold')
ax[1,1].bar(index[0],np.nanmean(cov_annmean_total_D013),bar_width,color='lightcoral')
ax[1,1].bar(index[1],np.nanmean(cov_annmean_total_D012),bar_width,color='red')
ax[1,1].bar(index[2],np.nanmean(cov_annmean_total_D014),bar_width,color='purple')
ax[1,1].bar(index[3],np.nanmean(cov_annmean_total_D016),bar_width,color='lightblue')
ax[1,1].bar(index[4],np.nanmean(cov_annmean_total_D015),bar_width,color='blue')
ax[1,1].bar(index[5],np.nanmean(cov_annmean_total_D017),bar_width,color='darkblue')
ax[1,1].bar(index[6],np.nanmean(cov_annmean_total_D019),bar_width,color='lightgreen')
ax[1,1].bar(index[7],np.nanmean(cov_annmean_total_D018),bar_width,color='green')
ax[1,1].bar(index[8],np.nanmean(cov_annmean_total_D020),bar_width,color='darkgreen')
ax[1,1].bar(index[9],np.nanmean(cov_annmean_total_D027),bar_width,color='lightcyan')
ax[1,1].bar(index[10],np.nanmean(cov_annmean_total_D021),bar_width,color='cyan')
ax[1,1].bar(index[11],np.nanmean(cov_annmean_total_D028),bar_width,color='darkcyan')
ax[1,1].bar(index[12],np.nanmean(cov_annmean_total_D029),bar_width,color='orange')
ax[1,1].bar(index[13],np.nanmean(cov_annmean_total_D022),bar_width,color='goldenrod')
ax[1,1].bar(index[14],np.nanmean(cov_annmean_total_D030),bar_width,color='brown')
ax[1,1].bar(index[15],np.nanmean(cov_annmean_total_D024),bar_width,color='lightgray')
ax[1,1].bar(index[16],np.nanmean(cov_annmean_total_D023),bar_width,color='gray')
ax[1,1].bar(index[17],np.nanmean(cov_annmean_total_D025),bar_width,color='black')
ax[1,1].set_yticks(np.arange(-20, 141, 20))
ax[1,1].set_xlabel('Experiment',fontsize=24)
ax[1,1].set_xticks(index)
ax[1,1].set_xticklabels(name_xticks,rotation=90)
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].axhline(c='k')

if save_fig == True:
    fig.savefig(dir_fig + 'fig3.png')
