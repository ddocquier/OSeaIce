#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Save ocean temperature (thetao) for Atlantic basin (including Atlantic part of the Arctic)
PROGRAMMER
    D. Docquier
LAST UPDATE
    29/04/2020
'''

# Standard libraries
import numpy as np
from netCDF4 import Dataset

# Options
exp = 'D013'
save_var = True
start_year = 2130
end_year = 2179

# Time parameters
start_folder = int(start_year - 2130 + 281)
nyears = int(end_year-start_year+1)
nmy = int(12) # number of months in a year

# Working directory
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/run/'+str(exp)+'/output/nemo/'
dir_mask = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/run/D000/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'+str(exp)+'/'

# Load dimensions
filename = dir_input + '281/' + str(exp)+'_1m_21300101_21301231_grid_T.nc'
fh = Dataset(filename, mode='r')
test = fh.variables['thetao'][:]
notused,nz,ny,nx = test.shape
fh.close()

# Load ocean mask
filename = dir_mask + 'subbasins.nc'
fh = Dataset(filename, mode='r')
atlmsk = fh.variables['atlmsk'][:]
fh.close()

# Compute 2D mean Atlantic ocean temperature (mean over all longitudes) averaged over 50 years
thetao_2D = np.zeros((nyears,nz,ny))
for year in np.arange(nyears):
    print(start_year+year)
    filename = dir_input + str(start_folder+year) + '/' + str(exp)+'_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_T.nc'
    fh = Dataset(filename, mode='r')
    thetao = fh.variables['thetao'][:]
    thetao[thetao>60.] = np.nan
    thetao = np.nanmean(thetao,axis=0) # annual mean
    fh.close()
    thetao_atl = thetao
    for z in np.arange(nz):
        thetao_atl[z,:,:][atlmsk == 0] = np.nan # NaN if not in the Atlantic
    thetao_2D[year,:,:] = np.nanmean(thetao_atl,axis=2) # average over all longitudes of the Atlantic
thetao_2D_mean = np.nanmean(thetao_2D,axis=0) # time mean over whole period

# Save thetao
if save_var == True:
    filename = dir_output + 'thetao_' + str(exp) + '.npy'
    np.save(filename,thetao_2D_mean)
