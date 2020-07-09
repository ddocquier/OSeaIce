#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig.1: Map of SST
PROGRAMMER
    D. Docquier
LAST UPDATE
    28/04/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Working directories
dir_D000 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D000/'
dir_D012 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D012/'
dir_D013 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D013/'
dir_D014 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D014/'
dir_D015 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D015/'
dir_D016 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D016/'
dir_D017 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D017/'
dir_D018 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D018/'
dir_D019 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D019/'
dir_D020 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D020/'
dir_D021 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D021/'
dir_D022 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D022/'
dir_D023 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D023/'
dir_D024 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D024/'
dir_D025 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D025/'
dir_D027 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D027/'
dir_D028 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D028/'
dir_D029 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D029/'
dir_D030 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D030/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'
dir_mask = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/surface_restoring/'

# Parameters
save_fig = False

# Load SST D000
filename = dir_D000 + 'tos_D000_ym.nc'
fh = Dataset(filename, mode='r')
tos_D000 = fh.variables['tos'][:]
lat = fh.variables['nav_lat'][:]
lon = fh.variables['nav_lon'][:]
nm,ny,nx = tos_D000.shape
fh.close()
tos_D000_annmean = np.nanmean(tos_D000,axis=0)

# Load SST D012
filename = dir_D012 + 'tos_D012_ym.nc'
fh = Dataset(filename, mode='r')
tos_D012 = fh.variables['tos'][:]
fh.close()
tos_D012_annmean = np.nanmean(tos_D012,axis=0)

# Load SST D013
filename = dir_D013 + 'tos_D013_ym.nc'
fh = Dataset(filename, mode='r')
tos_D013 = fh.variables['tos'][:]
fh.close()
tos_D013_annmean = np.nanmean(tos_D013,axis=0)

# Load SST D014
filename = dir_D014 + 'tos_D014_ym.nc'
fh = Dataset(filename, mode='r')
tos_D014 = fh.variables['tos'][:]
fh.close()
tos_D014_annmean = np.nanmean(tos_D014,axis=0)

# Load SST D015
filename = dir_D015 + 'tos_D015_ym.nc'
fh = Dataset(filename, mode='r')
tos_D015 = fh.variables['tos'][:]
fh.close()
tos_D015_annmean = np.nanmean(tos_D015,axis=0)

# Load SST D016
filename = dir_D016 + 'tos_D016_ym.nc'
fh = Dataset(filename, mode='r')
tos_D016 = fh.variables['tos'][:]
fh.close()
tos_D016_annmean = np.nanmean(tos_D016,axis=0)

# Load SST D017
filename = dir_D017 + 'tos_D017_ym.nc'
fh = Dataset(filename, mode='r')
tos_D017 = fh.variables['tos'][:]
fh.close()
tos_D017_annmean = np.nanmean(tos_D017,axis=0)

# Load SST D018
filename = dir_D018 + 'tos_D018_ym.nc'
fh = Dataset(filename, mode='r')
tos_D018 = fh.variables['tos'][:]
fh.close()
tos_D018_annmean = np.nanmean(tos_D018,axis=0)

# Load SST D019
filename = dir_D019 + 'tos_D019_ym.nc'
fh = Dataset(filename, mode='r')
tos_D019 = fh.variables['tos'][:]
fh.close()
tos_D019_annmean = np.nanmean(tos_D019,axis=0)

# Load SST D020
filename = dir_D020 + 'tos_D020_ym.nc'
fh = Dataset(filename, mode='r')
tos_D020 = fh.variables['tos'][:]
fh.close()
tos_D020_annmean = np.nanmean(tos_D020,axis=0)

# Load SST D021
filename = dir_D021 + 'tos_D021_ym.nc'
fh = Dataset(filename, mode='r')
tos_D021 = fh.variables['tos'][:]
fh.close()
tos_D021_annmean = np.nanmean(tos_D021,axis=0)

# Load SST D022
filename = dir_D022 + 'tos_D022_ym.nc'
fh = Dataset(filename, mode='r')
tos_D022 = fh.variables['tos'][:]
fh.close()
tos_D022_annmean = np.nanmean(tos_D022,axis=0)

# Load SST D023
filename = dir_D023 + 'tos_D023_ym.nc'
fh = Dataset(filename, mode='r')
tos_D023 = fh.variables['tos'][:]
fh.close()
tos_D023_annmean = np.nanmean(tos_D023,axis=0)

# Load SST D024
filename = dir_D024 + 'tos_D024_ym.nc'
fh = Dataset(filename, mode='r')
tos_D024 = fh.variables['tos'][:]
fh.close()
tos_D024_annmean = np.nanmean(tos_D024,axis=0)

# Load SST D025
filename = dir_D025 + 'tos_D025_ym.nc'
fh = Dataset(filename, mode='r')
tos_D025 = fh.variables['tos'][:]
fh.close()
tos_D025_annmean = np.nanmean(tos_D025,axis=0)

# Load SST D027
filename = dir_D027 + 'tos_D027_ym.nc'
fh = Dataset(filename, mode='r')
tos_D027 = fh.variables['tos'][:]
fh.close()
tos_D027_annmean = np.nanmean(tos_D027,axis=0)

# Load SST D028
filename = dir_D028 + 'tos_D028_ym.nc'
fh = Dataset(filename, mode='r')
tos_D028 = fh.variables['tos'][:]
fh.close()
tos_D028_annmean = np.nanmean(tos_D028,axis=0)

# Load SST D029
filename = dir_D029 + 'tos_D029_ym.nc'
fh = Dataset(filename, mode='r')
tos_D029 = fh.variables['tos'][:]
fh.close()
tos_D029_annmean = np.nanmean(tos_D029,axis=0)

# Load SST D030
filename = dir_D030 + 'tos_D030_ym.nc'
fh = Dataset(filename, mode='r')
tos_D030 = fh.variables['tos'][:]
fh.close()
tos_D030_annmean = np.nanmean(tos_D030,axis=0)

# Load mask ATL1
filename = dir_mask + 'mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_atl1 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask ATL2
filename = dir_mask + 'northernNA/mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_atl2 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask ATL3
filename = dir_mask + 'BSO/mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_atl3 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC1
filename = dir_mask + 'Pacific/mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_pac1 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC2
filename = dir_mask + 'easternPacific/mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_pac2 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC3
filename = dir_mask + 'Bering/mask_restore.nc'
fh = Dataset(filename, mode='r')
mask_pac3 = fh.variables['mask_ssr'][:]
fh.close()

# Mean SST difference in SST domain
print('Mean SST difference in SST domain')
print('ATL1+1K:',np.nanmean(tos_D013_annmean[mask_atl1==1]-tos_D000_annmean[mask_atl1==1]))
print('ATL1+3K:',np.nanmean(tos_D012_annmean[mask_atl1==1]-tos_D000_annmean[mask_atl1==1]))
print('ATL1+5K:',np.nanmean(tos_D014_annmean[mask_atl1==1]-tos_D000_annmean[mask_atl1==1]))
print('ATL2+1K:',np.nanmean(tos_D016_annmean[mask_atl2==1]-tos_D000_annmean[mask_atl2==1]))
print('ATL2+3K:',np.nanmean(tos_D015_annmean[mask_atl2==1]-tos_D000_annmean[mask_atl2==1]))
print('ATL2+5K:',np.nanmean(tos_D017_annmean[mask_atl2==1]-tos_D000_annmean[mask_atl2==1]))
print('ATL3+1K:',np.nanmean(tos_D019_annmean[mask_atl3==1]-tos_D000_annmean[mask_atl3==1]))
print('ATL3+3K:',np.nanmean(tos_D018_annmean[mask_atl3==1]-tos_D000_annmean[mask_atl3==1]))
print('ATL3+5K:',np.nanmean(tos_D020_annmean[mask_atl3==1]-tos_D000_annmean[mask_atl3==1]))
print('PAC1+1K:',np.nanmean(tos_D027_annmean[mask_pac1==1]-tos_D000_annmean[mask_pac1==1]))
print('PAC1+3K:',np.nanmean(tos_D021_annmean[mask_pac1==1]-tos_D000_annmean[mask_pac1==1]))
print('PAC1+5K:',np.nanmean(tos_D028_annmean[mask_pac1==1]-tos_D000_annmean[mask_pac1==1]))
print('PAC2+1K:',np.nanmean(tos_D029_annmean[mask_pac2==1]-tos_D000_annmean[mask_pac2==1]))
print('PAC2+3K:',np.nanmean(tos_D022_annmean[mask_pac2==1]-tos_D000_annmean[mask_pac2==1]))
print('PAC2+5K:',np.nanmean(tos_D030_annmean[mask_pac2==1]-tos_D000_annmean[mask_pac2==1]))
print('PAC3+1K:',np.nanmean(tos_D024_annmean[mask_pac3==1]-tos_D000_annmean[mask_pac3==1]))
print('PAC3+3K:',np.nanmean(tos_D023_annmean[mask_pac3==1]-tos_D000_annmean[mask_pac3==1]))
print('PAC3+5K:',np.nanmean(tos_D025_annmean[mask_pac3==1]-tos_D000_annmean[mask_pac3==1]))
                         
# Map projection
boundlat = 30.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)

# Palettes and plot parameters
palette_var = plt.cm.Reds._resample(20)
min_var = 0.
max_var = 30.
palette_diff = plt.cm.seismic._resample(40)
min_diff = -3.
max_diff = 3.


# Fig. 1 - SST+3K maps - 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,tos_D000_annmean,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,5,10,15,20,25,30],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST ($^\circ$C)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
cs=map.pcolor(x,y,tos_D012_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,tos_D015_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,tos_D018_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,tos_D021_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,tos_D022_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,tos_D023_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,mask_pac3,range(1,2,1),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST PERT - CTRL ($^\circ$C)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig1.png')


# Supp. Fig. 1b -  SST+1K maps - 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,tos_D000_annmean,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,5,10,15,20,25,30],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST ($^\circ$C)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D013
cs=map.pcolor(x,y,tos_D013_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+1$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D016
cs=map.pcolor(x,y,tos_D016_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+1$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D019
cs=map.pcolor(x,y,tos_D019_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+1$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D027
cs=map.pcolor(x,y,tos_D027_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+1$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D029
cs=map.pcolor(x,y,tos_D029_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+1$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D024
cs=map.pcolor(x,y,tos_D024_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,mask_pac3,range(1,2,1),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+1$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST PERT - CTRL ($^\circ$C)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig1b.png')


# Supp. Fig. 1c - SST+5K maps - 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,tos_D000_annmean,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,5,10,15,20,25,30],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST ($^\circ$C)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D014
cs=map.pcolor(x,y,tos_D014_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+5$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D017
cs=map.pcolor(x,y,tos_D017_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+5$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D020
cs=map.pcolor(x,y,tos_D020_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+5$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D028
cs=map.pcolor(x,y,tos_D028_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+5$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D030
cs=map.pcolor(x,y,tos_D030_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+5$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D025
cs=map.pcolor(x,y,tos_D025_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,mask_pac3,range(1,2,1),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+5$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST PERT - CTRL ($^\circ$C)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig1c.png')
