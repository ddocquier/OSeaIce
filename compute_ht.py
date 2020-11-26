#!/usr/bin/env python

'''
GOAL
    Compute meridional total heat transport
PROGRAMMER
    D. Docquier
LAST UPDATE
    19/11/2020
'''

# Options
exp = 'D023'
save_var = True

# Standard libraries
from netCDF4 import Dataset
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/'
dir_output = dir_input

# Function to compute total heat transport in PW - adapted from B. Rose, http://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture13%20--%20Heat%20transport.html
def inferred_heat_transport(energy_in,lat_deg):
    lat_rad = np.deg2rad(lat_deg)
    radius = 6.371e6 # Earth radius (m)
    ht = 2. * np.math.pi * radius**2 * integrate.cumtrapz(np.cos(lat_rad)*energy_in,x=lat_rad,initial=0.)
    ht = ht / 1.e15 # Convert from W to PW
    return ht

# Load top net solar (shortwave) radiation (J/m^2, positive downwards)
filename = dir_input + 'tsr_'+str(exp)+'_2130.nc'
fh = Dataset(filename, mode='r')
tsr = fh.variables['var178'][:]
lat = fh.variables['lat'][:]
fh.close()
nm,ny,nx = tsr.shape

# Load top net thermal (longwave) radiation (J/m^2, positive downwards)
filename = dir_input + 'ttr_'+str(exp)+'_2130.nc'
fh = Dataset(filename, mode='r')
ttr = fh.variables['var179'][:]
fh.close()

# Load restoring flux (W/m^2)
if exp != 'D000':
    filename = dir_input + 'hfcorr_ifs_'+str(exp)+'_2130-2179.nc'
    fh = Dataset(filename, mode='r')
    hfcorr = fh.variables['hfcorr'][:]
    hfcorr[hfcorr>1.e10] = 0.
    hfcorr[hfcorr<-1.e10] = 0.
    fh.close()

# Load surface ocean heat flux (W/m^2)
filename = dir_input + 'qtoce_ifs_'+str(exp)+'_2130.nc'
fh = Dataset(filename, mode='r')
ohfl = fh.variables['qt_oce'][:]
ohfl[ohfl>1.e10] = 0.
ohfl[ohfl<-1.e10] = 0.
fh.close()

# Convert radiations from J/m^2 to W/m^2
if exp == 'D000': # Output every 6h
    tsr = tsr / (6*3600)
    ttr = ttr / (6*3600)
else: # Output every 12h
    tsr = tsr / (12*3600)
    ttr = ttr / (12*3600)

# Zonal means
tsr_mean = np.nanmean(tsr,axis=2)
ttr_mean = np.nanmean(ttr,axis=2)
ohfl_mean = np.nanmean(ohfl,axis=2)
if exp != 'D000':
    hfcorr_mean = np.nanmean(hfcorr,axis=2)
    ohfl_mean = ohfl_mean + hfcorr_mean

# Compute top net downward radiation (sum of solar and thermal radiation, if both positive downwards)
Rt = tsr_mean + ttr_mean

# Compute surface net heat flux into the atmosphere
if exp != 'D000':
    Fatmin = Rt - (ohfl_mean - hfcorr_mean)
else:
    Fatmin = Rt - ohfl_mean

# Compute total, atmospheric and ocean heat transports from South Pole
lat = np.flip(lat)
Rt = np.flip(Rt,1)
Fatmin = np.flip(Fatmin,1)
ohfl_mean = np.flip(ohfl_mean,1)
ht_total = np.zeros((nm,ny))
aht = np.zeros((nm,ny))
oht = np.zeros((nm,ny))
for t in np.arange(nm):
    ht_total[t,:] = inferred_heat_transport(Rt[t,:],lat)
    aht[t,:] = inferred_heat_transport(Fatmin[t,:],lat)
    oht[t,:] = inferred_heat_transport(ohfl_mean[t,:],lat)

# Save variables
if save_var == True:
    filename = dir_output + 'ht_' + str(exp) + '.npy'
    np.save(filename,[ht_total,aht,oht,lat])
