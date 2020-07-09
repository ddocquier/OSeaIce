#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Compute time-mean OHF in each grid cell
PROGRAMMER
    D. Docquier
LAST UPDATE
    29/04/2020
'''

# Options
exp = 'D012'
save_var = True
start_year = 2130
end_year = 2179

# Standard libraries
import numpy as np
from netCDF4 import Dataset
import time
start_time = time.time()

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/OHT_transects/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'
dir_output = dir_input

# Load modeled OHT (computed with compute_oht.py)
filename = dir_input + 'oht_' + str(exp) + '.npy'
oht_u,oht_v = np.load(filename)

# Load grid sizes in X (gridx), Y (gridy)
filename = dir_grid + 'coordinates.nc'
fh = Dataset(filename, mode='r')
gridx = fh.variables['e1u'][:]
gridy = fh.variables['e2u'][:]
fh.close()

# Load grid size in Z (gridz) and depth
filename = dir_grid + 'mesh_zgr.nc'
fh = Dataset(filename, mode='r')
gridz = fh.variables['e3u'][:]
gridz = gridz[0,:,:,:]
mbathy = fh.variables['mbathy'][:]
mbathy = mbathy[0,:,:]
gdepu = fh.variables['gdepu'][:]
gdepu = gdepu[0,:,:,:]
depth=np.zeros((np.size(mbathy,0),np.size(mbathy,1)))
for i in np.arange(np.size(mbathy,0)):
    for j in np.arange(np.size(mbathy,1)):
        depth[i,j] = gdepu[mbathy[i,j],i,j]
nz,ny,nx = gridz.shape
fh.close()

# Dimensions
nyears = np.size(oht_u,0)
nmy = np.size(oht_u,1)
ny = np.size(oht_u,2)
nx = np.size(oht_u,3)

# Compute time-mean OHF (W/m^2) in each grid cell
oht_u_mean = np.nanmean(oht_u,axis=(0,1))
oht_u_mean = oht_u_mean / (gridy * depth)
oht_v_mean = np.nanmean(oht_v,axis=(0,1))
oht_v_mean = oht_v_mean / (gridx * depth)

# Check dimensions of OHF
print(np.size(oht_u_mean))
print(np.size(oht_u_mean,0))
print(np.size(oht_u_mean,1))

# Save variables
if save_var == True:
    filename = dir_output + 'oht_mean_' + str(exp) + '.npy'
    np.save(filename,[oht_u_mean,oht_v_mean])

# Computing time
print("--- %s seconds ---" % (time.time() - start_time))
