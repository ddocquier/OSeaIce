#!/usr/bin/env python

'''
GOAL
    Compute meridional total heat transport
PROGRAMMER
    D. Docquier
LAST UPDATE
    09/04/2020
'''

# Options
exp = 'D021'
save_var = True
save_fig = False

# Standard libraries
from netCDF4 import Dataset
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Working directories
dir_input = '/nobackup/rossby20/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/'
dir_output = dir_input

# Function to compute total heat transport in PW - adapted from B. Rose, http://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture13%20--%20Heat%20transport.html
def inferred_heat_transport(energy_in,lat_deg):
    lat_rad = np.deg2rad(lat_deg)
    radius = 6.371e6 # Earth radius (m)
    ht = 2. * np.math.pi * radius**2 * integrate.cumtrapz(np.cos(lat_rad)*energy_in,x=lat_rad,initial=0.)
    ht = ht / 1.e15 # Convert from W to PW
    return ht

# Load top net solar (shortwave) radiation (positive downwards)
#filename = dir_input + 'tsr_'+str(exp)+'_2130.nc'
filename = dir_input + 'tsr_'+str(exp)+'_2180.nc'
fh = Dataset(filename, mode='r')
tsr = fh.variables['var178'][:]
tsr[tsr>1.e10] = np.nan
lat = fh.variables['lat'][:]
fh.close()
nm,ny,nx = tsr.shape

# Load top net thermal (longwave) radiation (positive downwards)
#filename = dir_input + 'ttr_'+str(exp)+'_2130.nc'
filename = dir_input + 'ttr_'+str(exp)+'_2180.nc'
fh = Dataset(filename, mode='r')
ttr = fh.variables['var179'][:]
ttr[ttr>1.e10] = np.nan
fh.close()

# Convert radiations from J/m^2 to W/m^2
if exp == 'D000': # Output every 6h
    tsr = tsr / (6*3600)
    ttr = ttr / (6*3600)
else: # Output every 12h
    tsr = tsr / (12*3600)
    ttr = ttr / (12*3600)

# Zonal mean of solar and thermal radiations
tsr_mean = np.nanmean(tsr,axis=2)
ttr_mean = np.nanmean(ttr,axis=2)

# Compute top net downward radiation (sum of solar and thermal radiation, if both positive downwards)
Rt = tsr_mean + ttr_mean

# Compute total heat transport from South Pole
lat = np.flip(lat)
Rt = np.flip(Rt)
ht_total = np.zeros((nm,ny))
for t in np.arange(nm):
    ht_total[t,:] = inferred_heat_transport(Rt[t,:],lat)

# Save variables
if save_var == True:
    #filename = dir_output + 'ht_' + str(exp) + '.npy'
    filename = dir_output + 'ht_' + str(exp) + '_2180.npy'
    np.save(filename,[ht_total,lat])

# Figure
ticks = [-90,-60,-30,0,30,60,90]
fig,ax = plt.subplots()
ax.plot(lat,np.nanmean(tsr_mean,axis=0),label='Net downward solar')
ax.plot(lat,np.nanmean(-ttr_mean,axis=0),label='Upward thermal')
ax.plot(lat,np.nanmean(tsr_mean+ttr_mean,axis=0),label='Net downward ($R_t$)')
ax.set_ylabel('TOA radiation (W m$^{-2})$')
ax.set_xlabel('Latitude ($^{\circ}$)')
ax.set_xlim(-90,90)
ax.set_ylim(-150,350)
ax.set_xticks(ticks);
ax.legend()
ax.grid()
if save_fig == True:
    fig.savefig(dir_output + 'Fluxes_TOA_'+str(exp)+'.png')
