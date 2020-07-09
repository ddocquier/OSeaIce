#!/usr/bin/env python

'''
GOAL
    Compute OHT across Arctic transects (Barents Sea Opening, Bering Strait, Fram Strait, Davis Strait)
    By default, it is better to use method = 2 in the options below
PROGRAMMER
    D. Docquier
LAST UPDATE
    27/04/2020
'''

# Options
exp = 'D012'
start_year = 2130
transect = 2 # 1: Barents Sea Opening (AW: 20E, 71.5-73.5N; Arthun et al. [2012]); 2: Bering Strait (65.7N, 168-170.5W); 3: Fram Strait (79N, 20W-11E); 4: Barents Sea Opening (northern Norway to Bear Island: 20E, 70-74.5N); 5: Barents Sea Opening (NCC: 20E, 70-71N); 6: Barents Sea Opening full section until Svalbard (20E, 70-77N); 7: Nares Strait (78N, 76-73W); 8: Barrow Strait (90W, 74-75N); 9: Davis Strait (69N, 67-51W)
method = 2 # 1: Only u velocity; 2: u and v velocities (in case the next grid point of the transect is not on the same y line)
save_var = True

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import time
start_time = time.time()

# Working directory
dir_in = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/OHT_transects/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'

# Function to find nearest neighbour to a certain target
def find_nearest(point,target):
    a = np.abs(point-target)
    index = np.argmin(a)
    return index

# Load modeled OHT (computed with compute_oht.py)
filename = dir_in + 'oht_' + str(exp) + '.npy'
#filename = dir_in + 'oht_' + str(exp) + '_othercst.npy'
oht_u,oht_v = np.load(filename)

# Dimensions
nyears = np.size(oht_u,0)
nmy = np.size(oht_u,1)
ny = np.size(oht_u,2)
nx = np.size(oht_u,3)

# Load latitude and longitude of model
filename = dir_grid + 'mesh_zgr.nc'
fh = Dataset(filename, mode='r')
lon_u = fh.variables['nav_lon'][:]
lat_u = fh.variables['nav_lat'][:]
fh.close()

# Find grid points closest to transect
mask_lonlat = np.zeros((ny,nx))
if transect == 1 or transect == 4 or transect == 5 or transect == 6 or transect == 8: # Barents Sea Opening or Barrow Strait
    if transect == 8:
        lon_target = -90.
    else:
        lon_target = 20.
    if transect == 1:
        lat_min = 71.5
        lat_max = 73.5
    elif transect == 4:
        lat_min = 70.
        lat_max = 74.5
    elif transect == 5:
        lat_min = 70.
        lat_max = 71.
    elif transect == 6:
        lat_min = 70.
        lat_max = 77.
    elif transect == 8:
        lat_min = 74.
        lat_max = 75.
    mask_lat = np.zeros((ny,nx))
    mask_lat[(lat_u >= lat_min) * (lat_u <= lat_max)] = 1
    for y in np.arange(ny):
        x_lon = find_nearest(lon_u[y,:] * mask_lat[y,:],lon_target)
        mask_lonlat[y,x_lon] = mask_lat[y,x_lon]
    lon_threshold = 1. 
    mask_lonlat[(np.abs(lon_u - lon_target) > lon_threshold)] = 0
else:
    if transect == 2: # Bering Strait
        lat_target = 65.7
        lon_min = -170.5
        lon_max = -168.
    elif transect == 3: # Fram Strait
        lat_target = 79.
        lon_min = -20.
        lon_max = 11.
    elif transect == 7: # Nares Strait
        lat_target = 78.
        lon_min = -76.
        lon_max = -73.
    elif transect == 9: # Davis Strait
        lat_target = 69.
        lon_min = -67.
        lon_max = -51.
    mask_lon = np.zeros((ny,nx))
    mask_lon[(lon_u >= lon_min) * (lon_u <= lon_max)] = 1
    for x in np.arange(nx):
        y_lat = find_nearest(lat_u[:,x] * mask_lon[:,x],lat_target)
        mask_lonlat[y_lat,x] = mask_lon[y_lat,x]
    lat_threshold = 1.
    mask_lonlat[(np.abs(lat_u - lat_target) > lat_threshold)] = 0

# Compute OHT across transect (in TW)
oht_transect = np.zeros((nyears,nmy))
for year in np.arange(nyears):
    print(start_year+year)
    for t in np.arange(nmy):
        if transect == 1 or transect == 4 or transect == 5 or transect == 6 or transect == 8:
            oht_transect[year,t] = np.nansum(oht_u[year,t,:,:] * mask_lonlat)
        else:
            oht_transect[year,t] = np.nansum(oht_v[year,t,:,:] * mask_lonlat)
        if method == 2:
            if transect == 1 or transect == 4 or transect == 5 or transect == 6 or transect == 8:
                for y in np.arange(ny-1):
                    for x in np.arange(nx-1):
                        if mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x+1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + oht_v[year,t,y,x+1]
                        elif mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x-1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + oht_v[year,t,y,x]
            else:
                for x in np.arange(nx-1):
                    for y in np.arange(ny-2):
                        if mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x+1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + oht_u[year,t,y+1,x]
                        elif mask_lonlat[y,x] == 1 and mask_lonlat[y-1,x+1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + oht_u[year,t,y,x]
                        elif mask_lonlat[y,x] == 1 and mask_lonlat[y+2,x+1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + (oht_u[year,t,y+1,x] + oht_u[year,t,y+2,x]) / 2.
                        elif mask_lonlat[y,x] == 1 and mask_lonlat[y-2,x+1] == 1:
                            oht_transect[year,t] = oht_transect[year,t] + (oht_u[year,t,y,x] + oht_u[year,t,y-1,x]) / 2.
        oht_transect[year,t] = oht_transect[year,t] / 1.e12
    print(np.nanmean(oht_transect[year,:]),'TW')

# Save variables
if save_var == True:
    if method == 1:
        if transect == 1:
            filename = dir_in + 'oht_barents_aw_' + str(exp) + '.npy'
        elif transect == 2:
            filename = dir_in + 'oht_bering_' + str(exp) + '.npy'
        elif transect == 3:
            filename = dir_in + 'oht_fram_' + str(exp) + '.npy'
    elif method == 2:
        if transect == 1:
            filename = dir_in + 'oht_barents_aw_' + str(exp) + '.npy'
        elif transect == 2:
            filename = dir_in + 'oht_bering_' + str(exp) + '.npy'
            #filename = dir_in + 'oht_bering_' + str(exp) + '_othercst.npy'
        elif transect == 3:
            filename = dir_in + 'oht_fram_' + str(exp) + '.npy'
        elif transect == 4:
            filename = dir_in + 'oht_barents_bear_' + str(exp) + '.npy'
        elif transect == 5:
            filename = dir_in + 'oht_barents_ncc_' + str(exp) + '.npy'
        elif transect == 6:
            filename = dir_in + 'oht_full_barents_' + str(exp) + '.npy'
        elif transect == 7:
            filename = dir_in + 'oht_nares_' + str(exp) + '.npy'
        elif transect == 8:
            filename = dir_in + 'oht_barrow_' + str(exp) + '.npy'
        elif transect == 9:
            filename = dir_in + 'oht_davis_' + str(exp) + '.npy'
    np.save(filename,oht_transect)

# Coordinates of points along the transect
print(lat_u[mask_lonlat==1])
print(lon_u[mask_lonlat==1])
print(np.where(mask_lonlat==1))
indices_transect = np.where(mask_lonlat==1)

# Computing time
print("--- %s seconds ---" % (time.time() - start_time))
