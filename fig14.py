#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 14: Plot total HT, computed via compute_ht.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    19/11/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Option
save_fig = True

# Working directory
dir_D000 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D000/'
dir_D001 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D001/'
dir_D012 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D012/'
dir_D015 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D015/'
dir_D018 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D018/'
dir_D021 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D021/'
dir_D022 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D022/'
dir_D023 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D023/'
dir_obs = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/obs/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Load total HT D000
filename = dir_D000 + 'ht_D000.npy'
ht_D000,aht_D000,oht_D000b,notused = np.load(filename,allow_pickle=True)

# Load total HT D012
filename = dir_D012 + 'ht_D012.npy'
ht_D012,aht_D012,oht_D012b,lat_ifs = np.load(filename,allow_pickle=True)
lat_ifs = np.array(lat_ifs)
nm,ny = ht_D012.shape

# Load total HT D015
filename = dir_D015 + 'ht_D015.npy'
ht_D015,aht_D015,oht_D015b,notused = np.load(filename,allow_pickle=True)

# Load total HT D018
filename = dir_D018 + 'ht_D018.npy'
ht_D018,aht_D018,oht_D018b,notused = np.load(filename,allow_pickle=True)

# Load total HT D021
filename = dir_D021 + 'ht_D021.npy'
ht_D021,aht_D021,oht_D021b,notused = np.load(filename,allow_pickle=True)

# Load total HT D022
filename = dir_D022 + 'ht_D022.npy'
ht_D022,aht_D022,oht_D022b,notused = np.load(filename,allow_pickle=True)

# Load total HT D023
filename = dir_D023 + 'ht_D023.npy'
ht_D023,aht_D023,oht_D023b,notused = np.load(filename,allow_pickle=True)

# Load OHT D000
filename = dir_D001 + 'OHT_D001_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D000 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Load OHT D012
filename = dir_D012 + 'OHT_D012_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D012 = fh.variables['sopht_vt'][:,:,260]
lat = fh.variables['nav_lat'][:,260]
fh.close()

# Load OHT D015
filename = dir_D015 + 'OHT_D015_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D015 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Load OHT D018
filename = dir_D018 + 'OHT_D018_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D018 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Load OHT D021
filename = dir_D021 + 'OHT_D021_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D021 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Load OHT D022
filename = dir_D022 + 'OHT_D022_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D022 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Load OHT D023
filename = dir_D023 + 'OHT_D023_2130-2179.nc'
fh = Dataset(filename, mode='r')
oht_D023 = fh.variables['sopht_vt'][:,:,260]
fh.close()

# Compute mean total HT 50 years
start = 0
ht_mean_D000 = np.nanmean(ht_D000[start::,:],axis=0)
ht_mean_D012 = np.nanmean(ht_D012[start::,:],axis=0)
ht_mean_D015 = np.nanmean(ht_D015[start::,:],axis=0)
ht_mean_D018 = np.nanmean(ht_D018[start::,:],axis=0)
ht_mean_D021 = np.nanmean(ht_D021[start::,:],axis=0)
ht_mean_D022 = np.nanmean(ht_D022[start::,:],axis=0)
ht_mean_D023 = np.nanmean(ht_D023[start::,:],axis=0)

# Compute mean OHT 50 years
oht_mean_D000 = np.nanmean(oht_D000[start::,:],axis=0)
oht_mean_D012 = np.nanmean(oht_D012[start::,:],axis=0)
oht_mean_D015 = np.nanmean(oht_D015[start::,:],axis=0)
oht_mean_D018 = np.nanmean(oht_D018[start::,:],axis=0)
oht_mean_D021 = np.nanmean(oht_D021[start::,:],axis=0)
oht_mean_D022 = np.nanmean(oht_D022[start::,:],axis=0)
oht_mean_D023 = np.nanmean(oht_D023[start::,:],axis=0)

# Compute mean AHT 50 years
aht_mean_D000 = np.nanmean(aht_D000[start::,:],axis=0)
aht_mean_D012 = np.nanmean(aht_D012[start::,:],axis=0)
aht_mean_D015 = np.nanmean(aht_D015[start::,:],axis=0)
aht_mean_D018 = np.nanmean(aht_D018[start::,:],axis=0)
aht_mean_D021 = np.nanmean(aht_D021[start::,:],axis=0)
aht_mean_D022 = np.nanmean(aht_D022[start::,:],axis=0)
aht_mean_D023 = np.nanmean(aht_D023[start::,:],axis=0)


# Latitudinal transect of HT (50 years)
fig,ax = plt.subplots(3,2,figsize=(18,18))
fig.subplots_adjust(left=0.1,bottom=0.08,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# CTRL Atlantic SST experiments
ax[0,0].set_title('CTRL',fontsize=24)
ax[0,0].plot(lat_ifs,ht_mean_D000,'k-',linewidth=3,label='Total HT CTRL')
ax[0,0].plot(lat_ifs,aht_mean_D000,'k:',linewidth=2,label='AHT CTRL')
ax[0,0].plot(lat,oht_mean_D000,'k--',linewidth=2,label='OHT CTRL')
ax[0,0].legend(loc='upper left',shadow=True,frameon=False,fontsize=20)
ax[0,0].set_ylabel('Northward heat transport (PW)',fontsize=22)
ax[0,0].set_xticks(np.arange(-90, 91, 30))
ax[0,0].set_yticks(np.arange(-6, 6.1, 2))
ax[0,0].tick_params(axis='both',labelsize=20)
ax[0,0].axis([-90, 90, -6, 6])
ax[0,0].grid(linestyle='--')
ax[0,0].set_title('a',loc='left',fontsize=24,fontweight='bold')

# Remove axes that are not necessary
fig.delaxes(ax[0,1])

# AHT Atlantic SST experiments
ax[1,0].set_title('Atlantic SST+3$^\circ$C experiments',fontsize=24)
ax[1,0].plot(lat_ifs,aht_mean_D012-aht_mean_D000,'-',color='red',linewidth=2,label='ATL1+3$^\circ$C')
ax[1,0].plot(lat_ifs,aht_mean_D015-aht_mean_D000,'-',color='gray',linewidth=2,label='ATL2+3$^\circ$C')
ax[1,0].plot(lat_ifs,aht_mean_D018-aht_mean_D000,'-',color='blue',linewidth=2,label='ATL3+3$^\circ$C')
ax[1,0].legend(loc='upper left',shadow=True,frameon=False,fontsize=20)
ax[1,0].set_ylabel('dAHT (PW)',fontsize=22)
ax[1,0].set_xticks(np.arange(-90, 91, 30))
ax[1,0].set_yticks(np.arange(-0.4, 0.41, 0.2))
ax[1,0].tick_params(axis='both',labelsize=20)
ax[1,0].axis([-90, 90, -0.4, 0.4])
ax[1,0].grid(linestyle='--')
ax[1,0].set_title('b',loc='left',fontsize=24,fontweight='bold')

# AHT Pacific SST experiments
ax[1,1].set_title('Pacific SST+3$^\circ$C experiments',fontsize=24)
ax[1,1].plot(lat_ifs,aht_mean_D021-aht_mean_D000,'-',color='green',linewidth=2,label='PAC1+3$^\circ$C')
ax[1,1].plot(lat_ifs,aht_mean_D022-aht_mean_D000,'-',color='brown',linewidth=2,label='PAC2+3$^\circ$C')
ax[1,1].plot(lat_ifs,aht_mean_D023-aht_mean_D000,'-',color='orange',linewidth=2,label='PAC3+3$^\circ$C')
ax[1,1].legend(loc='upper left',shadow=True,frameon=False,fontsize=20)
ax[1,1].set_ylabel('dAHT (PW)',fontsize=22)
ax[1,1].set_xticks(np.arange(-90, 91, 30))
ax[1,1].set_yticks(np.arange(-0.4, 0.41, 0.2))
ax[1,1].tick_params(axis='both',labelsize=20)
ax[1,1].axis([-90, 90, -0.4, 0.4])
ax[1,1].grid(linestyle='--')
ax[1,1].set_title('c',loc='left',fontsize=24,fontweight='bold')

# OHT Atlantic SST experiments
ax[2,0].set_title('Atlantic SST+3$^\circ$C experiments',fontsize=24)
ax[2,0].plot(lat,oht_mean_D012-oht_mean_D000,'-',color='red',linewidth=2,label='ATL1+3$^\circ$C')
ax[2,0].plot(lat,oht_mean_D015-oht_mean_D000,'-',color='gray',linewidth=2,label='ATL2+3$^\circ$C')
ax[2,0].plot(lat,oht_mean_D018-oht_mean_D000,'-',color='blue',linewidth=2,label='ATL3+3$^\circ$C')
ax[2,0].legend(loc='upper left',shadow=True,frameon=False,fontsize=20)
ax[2,0].set_xlabel('Latitude ($^{\circ}$)',fontsize=22)
ax[2,0].set_ylabel('dOHT (PW)',fontsize=22)
ax[2,0].set_xticks(np.arange(-90, 91, 30))
ax[2,0].set_yticks(np.arange(-0.4, 0.41, 0.2))
ax[2,0].tick_params(axis='both',labelsize=20)
ax[2,0].axis([-90, 90, -0.4, 0.4])
ax[2,0].grid(linestyle='--')
ax[2,0].set_title('d',loc='left',fontsize=24,fontweight='bold')

# OHT Pacific SST experiments
ax[2,1].set_title('Pacific SST+3$^\circ$C experiments',fontsize=24)
ax[2,1].plot(lat,oht_mean_D021-oht_mean_D000,'-',color='green',linewidth=2,label='PAC1+3$^\circ$C')
ax[2,1].plot(lat,oht_mean_D022-oht_mean_D000,'-',color='brown',linewidth=2,label='PAC2+3$^\circ$C')
ax[2,1].plot(lat,oht_mean_D023-oht_mean_D000,'-',color='orange',linewidth=2,label='PAC3+3$^\circ$C')
ax[2,1].legend(loc='upper left',shadow=True,frameon=False,fontsize=18)
ax[2,1].set_xlabel('Latitude ($^{\circ}$)',fontsize=22)
ax[2,1].set_ylabel('dOHT (PW)',fontsize=20)
ax[2,1].set_xticks(np.arange(-90, 91, 30))
ax[2,1].set_yticks(np.arange(-0.4, 0.41, 0.2))
ax[2,1].tick_params(axis='both',labelsize=20)
ax[2,1].axis([-90, 90, -0.4, 0.4])
ax[2,1].grid(linestyle='--')
ax[2,1].set_title('e',loc='left',fontsize=24,fontweight='bold')

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig14.png')
    fig.savefig(dir_output + 'fig14.eps',dpi=300)
