#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Compute mean t2m (Arctic)
PROGRAMMER
    D. Docquier
LAST UPDATE
    13/05/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_D000 = dir_input + 'D000/IFS/'
dir_D012 = dir_input + 'D012/IFS/'
dir_D013 = dir_input + 'D013/IFS/'
dir_D014 = dir_input + 'D014/IFS/'
dir_D015 = dir_input + 'D015/IFS/'
dir_D016 = dir_input + 'D016/IFS/'
dir_D017 = dir_input + 'D017/IFS/'
dir_D018 = dir_input + 'D018/IFS/'
dir_D019 = dir_input + 'D019/IFS/'
dir_D020 = dir_input + 'D020/IFS/'
dir_D021 = dir_input + 'D021/IFS/'
dir_D022 = dir_input + 'D022/IFS/'
dir_D023 = dir_input + 'D023/IFS/'
dir_D024 = dir_input + 'D024/IFS/'
dir_D025 = dir_input + 'D025/IFS/'
dir_D027 = dir_input + 'D027/IFS/'
dir_D028 = dir_input + 'D028/IFS/'
dir_D029 = dir_input + 'D029/IFS/'
dir_D030 = dir_input + 'D030/IFS/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Parameters
start_year = 2130
nyears = 50

# Load lat,lon
filename = dir_D000 + 'ICMGGD000+2130.nc'
fh = Dataset(filename, mode='r')
test = fh.variables['var167'][:]
lat_init = fh.variables['lat'][:]
lon_init = fh.variables['lon'][:]
lon,lat = np.meshgrid(lon_init,lat_init)
lon[lon>180.] = lon[lon>180.] - 360.
nmy,ny,nx = test.shape
fh.close()

# Load grid-cell area
filename = dir_D000 + 'gridarea_ifs.nc'
fh = Dataset(filename, mode='r')
gridarea = fh.variables['cell_area'][:]
fh.close()

# Initialization of tas
tas_D000_mon = np.zeros((nyears,nmy,ny,nx))
tas_D012_mon = np.zeros((nyears,nmy,ny,nx))
tas_D013_mon = np.zeros((nyears,nmy,ny,nx))
tas_D014_mon = np.zeros((nyears,nmy,ny,nx))
tas_D015_mon = np.zeros((nyears,nmy,ny,nx))
tas_D016_mon = np.zeros((nyears,nmy,ny,nx))
tas_D017_mon = np.zeros((nyears,nmy,ny,nx))
tas_D018_mon = np.zeros((nyears,nmy,ny,nx))
tas_D019_mon = np.zeros((nyears,nmy,ny,nx))
tas_D020_mon = np.zeros((nyears,nmy,ny,nx))
tas_D021_mon = np.zeros((nyears,nmy,ny,nx))
tas_D022_mon = np.zeros((nyears,nmy,ny,nx))
tas_D023_mon = np.zeros((nyears,nmy,ny,nx))
tas_D024_mon = np.zeros((nyears,nmy,ny,nx))
tas_D025_mon = np.zeros((nyears,nmy,ny,nx))
tas_D027_mon = np.zeros((nyears,nmy,ny,nx))
tas_D028_mon = np.zeros((nyears,nmy,ny,nx))
tas_D029_mon = np.zeros((nyears,nmy,ny,nx))
tas_D030_mon = np.zeros((nyears,nmy,ny,nx))

# Load tas
for year in np.arange(nyears):
    print(start_year+year)

    # Retrieve tas D000
    filename = dir_D000 + 'ICMGGD000+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D000_init = fh.variables['var167'][:]
    tas_D000_init = tas_D000_init - 273.15
    fh.close()
    tas_D000_mon[year,:,:,:] = tas_D000_init

    # Retrieve tas D012
    filename = dir_D012 + 'ICMGGD012+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D012_init = fh.variables['var167'][:]
    tas_D012_init = tas_D012_init - 273.15
    fh.close()
    tas_D012_mon[year,:,:,:] = tas_D012_init

    # Retrieve tas D013
    filename = dir_D013 + 'ICMGGD013+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D013_init = fh.variables['var167'][:]
    tas_D013_init = tas_D013_init - 273.15
    fh.close()
    tas_D013_mon[year,:,:,:] = tas_D013_init

    # Retrieve tas D014
    filename = dir_D014 + 'ICMGGD014+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D014_init = fh.variables['var167'][:]
    tas_D014_init = tas_D014_init - 273.15
    fh.close()
    tas_D014_mon[year,:,:,:] = tas_D014_init

    # Retrieve tas D015
    filename = dir_D015 + 'ICMGGD015+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D015_init = fh.variables['var167'][:]
    tas_D015_init = tas_D015_init - 273.15
    fh.close()
    tas_D015_mon[year,:,:,:] = tas_D015_init

     # Retrieve tas D016
    filename = dir_D016 + 'ICMGGD016+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D016_init = fh.variables['var167'][:]
    tas_D016_init = tas_D016_init - 273.15
    fh.close()
    tas_D016_mon[year,:,:,:] = tas_D016_init

    # Retrieve tas D017
    filename = dir_D017 + 'ICMGGD017+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D017_init = fh.variables['var167'][:]
    tas_D017_init = tas_D017_init - 273.15
    fh.close()
    tas_D017_mon[year,:,:,:] = tas_D017_init

    # Retrieve tas D018
    filename = dir_D018 + 'ICMGGD018+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D018_init = fh.variables['var167'][:]
    tas_D018_init = tas_D018_init - 273.15
    fh.close()
    tas_D018_mon[year,:,:,:] = tas_D018_init

    # Retrieve tas D019
    filename = dir_D019 + 'ICMGGD019+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D019_init = fh.variables['var167'][:]
    tas_D019_init = tas_D019_init - 273.15
    fh.close()
    tas_D019_mon[year,:,:,:] = tas_D019_init

    # Retrieve tas D020
    filename = dir_D020 + 'ICMGGD020+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D020_init = fh.variables['var167'][:]
    tas_D020_init = tas_D020_init - 273.15
    fh.close()
    tas_D020_mon[year,:,:,:] = tas_D020_init

    # Retrieve tas D021
    filename = dir_D021 + 'ICMGGD021+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D021_init = fh.variables['var167'][:]
    tas_D021_init = tas_D021_init - 273.15
    fh.close()
    tas_D021_mon[year,:,:,:] = tas_D021_init

    # Retrieve tas D022
    filename = dir_D022 + 'ICMGGD022+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D022_init = fh.variables['var167'][:]
    tas_D022_init = tas_D022_init - 273.15
    fh.close()
    tas_D022_mon[year,:,:,:] = tas_D022_init

    # Retrieve tas D023
    filename = dir_D023 + 'ICMGGD023+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D023_init = fh.variables['var167'][:]
    tas_D023_init = tas_D023_init - 273.15
    fh.close()
    tas_D023_mon[year,:,:,:] = tas_D023_init

    # Retrieve tas D024
    filename = dir_D024 + 'ICMGGD024+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D024_init = fh.variables['var167'][:]
    tas_D024_init = tas_D024_init - 273.15
    fh.close()
    tas_D024_mon[year,:,:,:] = tas_D024_init

    # Retrieve tas D025
    filename = dir_D025 + 'ICMGGD025+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D025_init = fh.variables['var167'][:]
    tas_D025_init = tas_D025_init - 273.15
    fh.close()
    tas_D025_mon[year,:,:,:] = tas_D025_init

    # Retrieve tas D027
    filename = dir_D027 + 'ICMGGD027+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D027_init = fh.variables['var167'][:]
    tas_D027_init = tas_D027_init - 273.15
    fh.close()
    tas_D027_mon[year,:,:,:] = tas_D027_init

    # Retrieve tas D028
    filename = dir_D028 + 'ICMGGD028+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D028_init = fh.variables['var167'][:]
    tas_D028_init = tas_D028_init - 273.15
    fh.close()
    tas_D028_mon[year,:,:,:] = tas_D028_init
    
    # Retrieve tas D029
    filename = dir_D029 + 'ICMGGD029+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D029_init = fh.variables['var167'][:]
    tas_D029_init = tas_D029_init - 273.15
    fh.close()
    tas_D029_mon[year,:,:,:] = tas_D029_init

    # Retrieve tas D030
    filename = dir_D030 + 'ICMGGD030+'+str(start_year+year)+'.nc'
    fh = Dataset(filename, mode='r')
    tas_D030_init = fh.variables['var167'][:]
    tas_D030_init = tas_D030_init - 273.15
    fh.close()
    tas_D030_mon[year,:,:,:] = tas_D030_init

# Compute annual mean
tas_D000_annmean = np.nanmean(tas_D000_mon,axis=1)
tas_D012_annmean = np.nanmean(tas_D012_mon,axis=1)
tas_D013_annmean = np.nanmean(tas_D013_mon,axis=1)
tas_D014_annmean = np.nanmean(tas_D014_mon,axis=1)
tas_D015_annmean = np.nanmean(tas_D015_mon,axis=1)
tas_D016_annmean = np.nanmean(tas_D016_mon,axis=1)
tas_D017_annmean = np.nanmean(tas_D017_mon,axis=1)
tas_D018_annmean = np.nanmean(tas_D018_mon,axis=1)
tas_D019_annmean = np.nanmean(tas_D029_mon,axis=1)
tas_D020_annmean = np.nanmean(tas_D020_mon,axis=1)
tas_D021_annmean = np.nanmean(tas_D021_mon,axis=1)
tas_D022_annmean = np.nanmean(tas_D022_mon,axis=1)
tas_D023_annmean = np.nanmean(tas_D023_mon,axis=1)
tas_D024_annmean = np.nanmean(tas_D024_mon,axis=1)
tas_D025_annmean = np.nanmean(tas_D025_mon,axis=1)
tas_D027_annmean = np.nanmean(tas_D027_mon,axis=1)
tas_D028_annmean = np.nanmean(tas_D028_mon,axis=1)
tas_D029_annmean = np.nanmean(tas_D029_mon,axis=1)
tas_D030_annmean = np.nanmean(tas_D030_mon,axis=1)

# Save variables
filename = dir_output + 'tas_annmean.npy'
np.save(filename,[tas_D000_annmean,tas_D012_annmean,tas_D013_annmean,tas_D014_annmean,tas_D015_annmean,tas_D016_annmean,tas_D017_annmean,tas_D018_annmean,tas_D019_annmean,tas_D020_annmean,tas_D021_annmean,tas_D022_annmean,tas_D023_annmean,tas_D024_annmean,tas_D025_annmean,tas_D027_annmean,tas_D028_annmean,tas_D029_annmean,tas_D030_annmean])

