#!/usr/bin/env python


'''
GOAL
    Compute monthly mean Arctic sea-ice area (10^6 km^2) for total Arctic and specific regions (based on Koenigk et al., 2016)
    For specific Arctic regions, see JASMIN: /gws/nopw/j04/primavera1/tools/WP2/Topic7/ICE_AREA_REG_sep06.sh
    EC-Earth (T255-ORCA1)
PROGRAMMER
    D. Docquier
LAST UPDATE
    27/04/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np

# Options
exp = 'D012'
save_var = True # True: save variables in a .npy file; False: don't save variables

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/'
dir_output = dir_input
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'

# Function to compute total Arctic sea-ice area (10^6 km^2) - NEMO grid
def compute_area(nm,lat,lat_threshold,siconc,grid_area,mask):
    area = np.zeros(nm)
    for i in np.arange(nm):
        area[i] = np.nansum((lat >= lat_threshold) * (siconc[i,:,:] / 100.) * (mask == 0) * grid_area)
    area = area / 1.e6
    area[area <= 0.] = np.nan
    return area

# Function to compute Arctic sea-ice area in specific regions (10^6 km^2) - NEMO grid
def compute_area_region(nm,lat,lon,lat1,lat2,lon1,lon2,siconc,grid_area,mask):
    area = np.zeros(nm)
    for i in np.arange(nm):
        area[i] = np.nansum((lat >= lat1) * (lat <= lat2) * (lon >= lon1) * (lon <= lon2) * (siconc[i,:,:] / 100.) * (mask == 0) * grid_area)
    area = area / 1.e6
    area[area <= 0.] = np.nan      
    return area

# Function to compute Arctic sea-ice area in Chukchi Sea (10^6 km^2) - NEMO grid
def compute_area_region2(nm,lat,lon,lat1,lat2,lon1,lon2,lon3,lon4,siconc,grid_area,mask):
    area = np.zeros(nm)
    for i in np.arange(nm):
        area1 = np.nansum((lat >= lat1) * (lat <= lat2) * (lon >= lon1) * (lon <= lon2) * (siconc[i,:,:] / 100.) * (mask == 0) * grid_area)
        area2 = np.nansum((lat >= lat1) * (lat <= lat2) * (lon >= lon3) * (lon <= lon4) * (siconc[i,:,:] / 100.) * (mask == 0) * grid_area)
        area[i] = area1 + area2
    area = area / 1.e6
    area[area <= 0.] = np.nan      
    return area

# Select siconc files
file_siconc = dir_output + 'siconc_' + str(exp) + '_2130.nc'
    
# Load siconc EC-Earth
fh = Dataset(file_siconc, mode='r')
siconc = fh.variables['siconc'][:]
siconc = siconc * 100.
siconc[siconc < 0.] = np.nan
siconc[siconc > 101.] = np.nan
nm,ny,nx = siconc.shape

# Load latitude (-90 to +90) and longitude (-180 to +180) EC-Earth
lat = fh.variables['nav_lat_grid_T'][:]
lon = fh.variables['nav_lon_grid_T'][:]

# Load grid-cell area (km^2) EC-Earth
grid_area = fh.variables['cell_area'][:]
grid_area = grid_area / 1.e6
fh.close()

# Load ocean mask EC-Earth
file_mask = dir_grid + 'masks.nc'
fh = Dataset(file_mask,mode='r')
mask_ocean = fh.variables['O1t0.msk'][:]
fh.close()

# Compute total Arctic sea-ice area (10^6 km^2)
lat_total = 40.
area_total = compute_area(nm,lat,lat_total,siconc,grid_area,mask_ocean)
print('Total Arctic done')

# Compute sea-ice area of Barents / Kara Seas (10^6 km^2)
lat_bar1 = 70.
lat_bar2 = 81.
lon_bar1 = 15.
lon_bar2 = 100.
area_barents =  compute_area_region(nm,lat,lon,lat_bar1,lat_bar2,lon_bar1,lon_bar2,siconc,grid_area,mask_ocean)
print('Barents done')

# Compute sea-ice area of GIN Seas (10^6 km^2)
lat_gre1 = 55.
lat_gre2 = 81.
lon_gre1 = -40.
lon_gre2 = 15.
area_greenland =  compute_area_region(nm,lat,lon,lat_gre1,lat_gre2,lon_gre1,lon_gre2,siconc,grid_area,mask_ocean)
print('GIN done')

# Compute sea-ice area of Labrador Sea / Baffin Bay (10^6 km^2)
lat_lab1 = 50.
lat_lab2 = 81.
lon_lab1 = -90.
lon_lab2 = -40.
area_labrador =  compute_area_region(nm,lat,lon,lat_lab1,lat_lab2,lon_lab1,lon_lab2,siconc,grid_area,mask_ocean)
print('Labrador done')

# Compute sea-ice area of Laptev / East Siberian Seas (10^6 km^2)
lat_lap1 = 70.
lat_lap2 = 81.
lon_lap1 = 100.
lon_lap2 = 170.
area_laptev =  compute_area_region(nm,lat,lon,lat_lap1,lat_lap2,lon_lap1,lon_lap2,siconc,grid_area,mask_ocean)
print('Laptev done')

# Compute sea-ice area of Chukchi / Bering Seas (10^6 km^2)
lat_chu1 = 50.
lat_chu2 = 81.
lon_chu1 = -180.
lon_chu2 = -160.
lon_chu3 = 170.
lon_chu4 = 180.
area_chukchi =  compute_area_region2(nm,lat,lon,lat_chu1,lat_chu2,lon_chu1,lon_chu2,lon_chu3,lon_chu4,siconc,grid_area,mask_ocean)
print('Chukchi done')

# Compute sea-ice area of Beaufort Sea (10^6 km^2)
lat_bea1 = 70.
lat_bea2 = 81.
lon_bea1 = -160.
lon_bea2 = -90.
area_beaufort =  compute_area_region(nm,lat,lon,lat_bea1,lat_bea2,lon_bea1,lon_bea2,siconc,grid_area,mask_ocean)
print('Beaufort done')

# Compute sea-ice area of Central Arctic (10^6 km^2)
lat_cen = 81.
area_central =  compute_area(nm,lat,lat_cen,siconc,grid_area,mask_ocean)
print('Central Arctic done')

# Save variables
if save_var == True:
    filename = dir_output + 'SIarea_' + str(exp) + '.npy'
    np.save(filename,[area_total,area_barents,area_greenland,area_labrador,area_laptev,area_chukchi,area_beaufort,area_central])
