#!/usr/bin/env python

'''
GOAL
    Compute OHT velocity and temperature anomalies across Arctic transects (Barents Sea Opening, Bering Strait, Fram Strait, Davis Strait)
PROGRAMMER
    D. Docquier
LAST UPDATE
    23/04/2020
'''

# Options
exp = 'D012'
start_year = 2130
transect = 1 # 1: Barents Sea Opening (20E, 70-77N); 2: Bering Strait (65.7N, 168-170.5W); 3: Fram Strait (79N, 20W-11E); 4: Davis Strait (69N, 67-51W)
save_var = True

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import time
start_time = time.time()

# Working directory
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/OHT_transects/'

# Function to find nearest neighbour to a certain target
def find_nearest(point,target):
    a = np.abs(point-target)
    index = np.argmin(a)
    return index

# Load modeled OHT velocity and temperature anomalies (computed with decompose_oht.py)
filename = dir_output + 'decoht_' + str(exp) + '.npy'
u_ano,v_ano,t_u_ano,t_v_ano,cov_u_ano,cov_v_ano = np.load(filename)

# Dimensions
nyears = np.size(u_ano,0)
nmy = np.size(u_ano,1)
ny = np.size(u_ano,2)
nx = np.size(u_ano,3)

# Load latitude and longitude of model
filename = dir_grid + 'mesh_zgr.nc'
fh = Dataset(filename, mode='r')
lon_u = fh.variables['nav_lon'][:]
lat_u = fh.variables['nav_lat'][:]
fh.close()

# Find grid points closest to transect
mask_lonlat = np.zeros((ny,nx))
if transect == 1: # Barents Sea Opening
    lon_target = 20.
    lat_min = 70.
    lat_max = 77.
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
    elif transect == 4: # Davis Strait
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

# Coordinates of points along the transect
print(lat_u[mask_lonlat==1])
print(lon_u[mask_lonlat==1])
print(np.where(mask_lonlat==1))
indices_transect = np.where(mask_lonlat==1)

# Compute OHT anomaly across transect (in TW)
vel_ano_transect = np.zeros((nyears,nmy))
temp_ano_transect = np.zeros((nyears,nmy))
cov_ano_transect = np.zeros((nyears,nmy))
for year in np.arange(nyears):
    print(start_year+year)
    for t in np.arange(nmy):
        if transect == 1:
            vel_ano_transect[year,t] = np.nansum(u_ano[year,t,:,:] * mask_lonlat)
            temp_ano_transect[year,t] = np.nansum(t_u_ano[year,t,:,:] * mask_lonlat)
            cov_ano_transect[year,t] = np.nansum(cov_u_ano[year,t,:,:] * mask_lonlat)
        else:
            vel_ano_transect[year,t] = np.nansum(v_ano[year,t,:,:] * mask_lonlat)
            temp_ano_transect[year,t] = np.nansum(t_v_ano[year,t,:,:] * mask_lonlat)
            cov_ano_transect[year,t] = np.nansum(cov_v_ano[year,t,:,:] * mask_lonlat)
        if transect == 1:
            for y in np.arange(ny-1):
                for x in np.arange(nx-1):
                    if mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x+1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + v_ano[year,t,y,x+1]
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + t_v_ano[year,t,y,x+1]
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + cov_v_ano[year,t,y,x+1]
                    elif mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x-1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + v_ano[year,t,y,x]
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + t_v_ano[year,t,y,x]
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + cov_v_ano[year,t,y,x]
        else:
            for x in np.arange(nx-1):
                for y in np.arange(ny-2):
                    if mask_lonlat[y,x] == 1 and mask_lonlat[y+1,x+1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + u_ano[year,t,y+1,x]
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + t_u_ano[year,t,y+1,x]
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + cov_u_ano[year,t,y+1,x]
                    elif mask_lonlat[y,x] == 1 and mask_lonlat[y-1,x+1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + u_ano[year,t,y,x]
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + t_u_ano[year,t,y,x]
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + cov_u_ano[year,t,y,x]
                    elif mask_lonlat[y,x] == 1 and mask_lonlat[y+2,x+1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + (u_ano[year,t,y+1,x] + u_ano[year,t,y+2,x]) / 2.
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + (t_u_ano[year,t,y+1,x] + t_u_ano[year,t,y+2,x]) / 2.
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + (cov_u_ano[year,t,y+1,x] + cov_u_ano[year,t,y+2,x]) / 2.
                    elif mask_lonlat[y,x] == 1 and mask_lonlat[y-2,x+1] == 1:
                        vel_ano_transect[year,t] = vel_ano_transect[year,t] + (u_ano[year,t,y,x] + u_ano[year,t,y-1,x]) / 2.
                        temp_ano_transect[year,t] = temp_ano_transect[year,t] + (t_u_ano[year,t,y,x] + t_u_ano[year,t,y-1,x]) / 2.
                        cov_ano_transect[year,t] = cov_ano_transect[year,t] + (cov_u_ano[year,t,y,x] + cov_u_ano[year,t,y-1,x]) / 2.
        vel_ano_transect[year,t] = vel_ano_transect[year,t] / 1.e12
        temp_ano_transect[year,t] = temp_ano_transect[year,t] / 1.e12
        cov_ano_transect[year,t] = cov_ano_transect[year,t] / 1.e12
    print(np.nanmean(vel_ano_transect[year,:]),'TW')
    print(np.nanmean(temp_ano_transect[year,:]),'TW')
    print(np.nanmean(cov_ano_transect[year,:]),'TW')

# Save variables
if save_var == True:
    if transect == 1:
        filename = dir_output + 'decoht_barents_' + str(exp) + '.npy'
    elif transect == 2:
        filename = dir_output + 'decoht_bering_' + str(exp) + '.npy'
    elif transect == 3:
        filename = dir_output + 'decoht_fram_' + str(exp) + '.npy'
    elif transect == 4:
        filename = dir_output + 'decoht_davis_' + str(exp) + '.npy'
    np.save(filename,[vel_ano_transect,temp_ano_transect,cov_ano_transect])

# Computing time
print("--- %s seconds ---" % (time.time() - start_time))
