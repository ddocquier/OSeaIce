#!/usr/bin/env python

'''
GOAL
    Save surface fluxes from SST experiments (IFS outputs)
PROGRAMMER
    D. Docquier
LAST UPDATE
    14/10/2020
'''

# Options
exp = 'D013'
save_var = True
plot_fig = False

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/'
dir_output = dir_input

# Load surface net solar (shortwave) radiation (positive downwards)
filename = dir_input + 'ssr_'+ str(exp) + '_2130.nc'
fh = Dataset(filename, mode='r')
swrad = fh.variables['var176'][:]
swrad[swrad>1.e10] = np.nan
swrad[swrad<-1.e10] = np.nan
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
fh.close()
nm,ny,nx = swrad.shape

# Load surface net thermal (longwave) radiation (positive downwards)
filename = dir_input + 'str_' + str(exp) + '_2130.nc'
fh = Dataset(filename, mode='r')
lwrad = fh.variables['var177'][:]
lwrad[lwrad>1.e10] = np.nan
lwrad[lwrad<-1.e10] = np.nan
fh.close()

# Load surface sensible heat flux (positive downwards)
filename = dir_input + 'sshf_' + str(exp) + '_2130.nc'
fh = Dataset(filename, mode='r')
shflux = fh.variables['var146'][:]
shflux[shflux>1.e10] = np.nan
shflux[shflux<-1.e10] = np.nan
fh.close()

# Load surface latent heat flux (positive downwards)
filename = dir_input + 'slhf_' + str(exp) + '_2130.nc'
fh = Dataset(filename, mode='r')
lhflux = fh.variables['var147'][:]
lhflux[lhflux>1.e10] = np.nan
lhflux[lhflux<-1.e10] = np.nan
fh.close()

# Convert radiations from J/m^2 to W/m^2 (1 W = 1 J/s)
if exp == 'D000': # Output every 6h
    swrad = swrad / (6*3600)
    lwrad = lwrad / (6*3600)
    shflux = shflux / (6*3600)
    lhflux = lhflux / (6*3600)
else: # Output every 12h
    swrad = swrad / (12*3600)
    lwrad = lwrad / (12*3600)
    shflux = shflux / (12*3600)
    lhflux = lhflux / (12*3600)

# Zonal mean of solar and thermal radiations
swrad_mean = np.nanmean(swrad,axis=2)
lwrad_mean = np.nanmean(lwrad,axis=2)
shflux_mean = np.nanmean(shflux,axis=2)
lhflux_mean = np.nanmean(lhflux,axis=2)

# Compute surface net downward radiation
Rts_mean = swrad_mean + lwrad_mean + shflux_mean + lhflux_mean
Rts = np.nanmean(swrad + lwrad + shflux + lhflux,axis=0)

# Save variables
if save_var == True:
    filename = dir_output + 'SurfaceFluxes_' + str(exp) + '.npy'
    np.save(filename,[swrad_mean,lwrad_mean,shflux_mean,lhflux_mean,Rts_mean,Rts,lat,lon])

# Figure
if plot_fig == True:
    ticks = [-90,-60,-30,0,30,60,90]
    fig,ax = plt.subplots()
    ax.plot(lat,np.nanmean(swrad_mean,axis=0),label='Shortwave radiation')
    ax.plot(lat,np.nanmean(lwrad_mean,axis=0),label='Longwave radiation')
    ax.plot(lat,np.nanmean(shflux_mean,axis=0),label='Sensible heat flux')
    ax.plot(lat,np.nanmean(lhflux_mean,axis=0),label='Latent heat flux')
    ax.plot(lat,np.nanmean(Rts_mean,axis=0),label='Net heat flux')
    ax.set_ylabel('Surface heat flux (W m$^{-2})$')
    ax.set_xlabel('Latitude ($^{\circ}$)')
    ax.set_xlim(-90,90)
    ax.set_ylim(-150,300)
    ax.set_xticks(ticks);
    ax.legend()
    ax.grid()
    fig.savefig(dir_output + 'Fluxes_Surface_' + str(exp) + '.png')
