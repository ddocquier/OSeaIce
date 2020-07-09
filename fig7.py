#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 7: Maps of sea-ice concentration (siconc)
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
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/'
dir_D000 = dir_input + 'post-proc/D000/'
dir_D012 = dir_input + 'post-proc/D012/'
dir_D013 = dir_input + 'post-proc/D013/'
dir_D014 = dir_input + 'post-proc/D014/'
dir_D015 = dir_input + 'post-proc/D015/'
dir_D016 = dir_input + 'post-proc/D016/'
dir_D017 = dir_input + 'post-proc/D017/'
dir_D018 = dir_input + 'post-proc/D018/'
dir_D019 = dir_input + 'post-proc/D019/'
dir_D020 = dir_input + 'post-proc/D020/'
dir_D021 = dir_input + 'post-proc/D021/'
dir_D022 = dir_input + 'post-proc/D022/'
dir_D023 = dir_input + 'post-proc/D023/'
dir_D024 = dir_input + 'post-proc/D024/'
dir_D025 = dir_input + 'post-proc/D025/'
dir_D027 = dir_input + 'post-proc/D027/'
dir_D028 = dir_input + 'post-proc/D028/'
dir_D029 = dir_input + 'post-proc/D029/'
dir_D030 = dir_input + 'post-proc/D030/'
dir_output = dir_input + 'OSeaIce_Paper/'
dir_mask = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/surface_restoring/'

# Parameters
save_fig = False
month = 3
month2 = 9

# Load siconc D000
filename = dir_D000 + 'siconc_D000_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D000 = fh.variables['siconc'][:]
siconc_D000 = siconc_D000 * 100.
lat = fh.variables['nav_lat_grid_T'][:]
lon = fh.variables['nav_lon_grid_T'][:]
nm,ny,nx = siconc_D000.shape
fh.close()

# Load siconc D012
filename = dir_D012 + 'siconc_D012_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D012 = fh.variables['siconc'][:]
siconc_D012 = siconc_D012 * 100.
fh.close()

# Load siconc D013
filename = dir_D013 + 'siconc_D013_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D013 = fh.variables['siconc'][:]
siconc_D013 = siconc_D013 * 100.
fh.close()

# Load siconc D014
filename = dir_D014 + 'siconc_D014_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D014 = fh.variables['siconc'][:]
siconc_D014 = siconc_D014 * 100.
fh.close()

# Load siconc D015
filename = dir_D015 + 'siconc_D015_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D015 = fh.variables['siconc'][:]
siconc_D015 = siconc_D015 * 100.
fh.close()

# Load siconc D016
filename = dir_D016 + 'siconc_D016_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D016 = fh.variables['siconc'][:]
siconc_D016 = siconc_D016 * 100.
fh.close()

# Load siconc D017
filename = dir_D017 + 'siconc_D017_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D017 = fh.variables['siconc'][:]
siconc_D017 = siconc_D017 * 100.
fh.close()

# Load siconc D018
filename = dir_D018 + 'siconc_D018_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D018 = fh.variables['siconc'][:]
siconc_D018 = siconc_D018 * 100.
fh.close()

# Load siconc D019
filename = dir_D019 + 'siconc_D019_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D019 = fh.variables['siconc'][:]
siconc_D019 = siconc_D019 * 100.
fh.close()

# Load siconc D020
filename = dir_D020 + 'siconc_D020_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D020 = fh.variables['siconc'][:]
siconc_D020 = siconc_D020 * 100.
fh.close()

# Load siconc D021
filename = dir_D021 + 'siconc_D021_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D021 = fh.variables['siconc'][:]
siconc_D021 = siconc_D021 * 100.
fh.close()

# Load siconc D022
filename = dir_D022 + 'siconc_D022_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D022 = fh.variables['siconc'][:]
siconc_D022 = siconc_D022 * 100.
fh.close()

# Load siconc D023
filename = dir_D023 + 'siconc_D023_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D023 = fh.variables['siconc'][:]
siconc_D023 = siconc_D023 * 100.
fh.close()

# Load siconc D024
filename = dir_D024 + 'siconc_D024_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D024 = fh.variables['siconc'][:]
siconc_D024 = siconc_D024 * 100.
fh.close()

# Load siconc D025
filename = dir_D025 + 'siconc_D025_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D025 = fh.variables['siconc'][:]
siconc_D025 = siconc_D025 * 100.
fh.close()

# Load siconc D027
filename = dir_D027 + 'siconc_D027_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D027 = fh.variables['siconc'][:]
siconc_D027 = siconc_D027 * 100.
fh.close()

# Load siconc D028
filename = dir_D028 + 'siconc_D028_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D028 = fh.variables['siconc'][:]
siconc_D028 = siconc_D028 * 100.
fh.close()

# Load siconc D029
filename = dir_D029 + 'siconc_D029_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D029 = fh.variables['siconc'][:]
siconc_D029 = siconc_D029 * 100.
fh.close()

# Load siconc D030
filename = dir_D030 + 'siconc_D030_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D030 = fh.variables['siconc'][:]
siconc_D030 = siconc_D030 * 100.
fh.close()

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
                    
# Map projection
boundlat = 40.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)

# Palettes and plot parameters
palette_var = plt.cm.Blues_r._resample(20)
min_var = 0.
max_var = 100.
palette_diff = plt.cm.seismic_r._resample(40)
min_diff = -50.
max_diff = 50.


# Fig. 7 - Maps of difference in March siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D012
cs = map.pcolor(x,y,np.squeeze(siconc_D012[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D015
cs = map.pcolor(x,y,np.squeeze(siconc_D015[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D018
cs=map.pcolor(x,y,np.squeeze(siconc_D018[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D021
cs = map.pcolor(x,y,np.squeeze(siconc_D021[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D022
cs = map.pcolor(x,y,np.squeeze(siconc_D022[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D023
cs = map.pcolor(x,y,np.squeeze(siconc_D023[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7.png')


# Supp. Fig. 7b (SST+1K) - Maps of difference in March siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D013
cs = map.pcolor(x,y,np.squeeze(siconc_D013[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+1$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D016
cs = map.pcolor(x,y,np.squeeze(siconc_D016[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+1$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D019
cs=map.pcolor(x,y,np.squeeze(siconc_D019[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+1$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D027
cs = map.pcolor(x,y,np.squeeze(siconc_D027[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+1$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D029
cs = map.pcolor(x,y,np.squeeze(siconc_D029[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+1$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D024
cs = map.pcolor(x,y,np.squeeze(siconc_D024[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7b.png')


# Supp. Fig. 7c (SST+3K) - Maps of difference in March siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D014
cs = map.pcolor(x,y,np.squeeze(siconc_D014[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+5$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D017
cs = map.pcolor(x,y,np.squeeze(siconc_D017[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+5$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D020
cs=map.pcolor(x,y,np.squeeze(siconc_D020[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+5$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D028
cs = map.pcolor(x,y,np.squeeze(siconc_D028[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+5$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D030
cs = map.pcolor(x,y,np.squeeze(siconc_D030[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+5$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D025
cs = map.pcolor(x,y,np.squeeze(siconc_D025[month-1,:,:]-siconc_D000[month-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7c.png')


# Supp. Fig. 7d - Maps of difference in September siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month2-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D012
cs = map.pcolor(x,y,np.squeeze(siconc_D012[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D015
cs = map.pcolor(x,y,np.squeeze(siconc_D015[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D018
cs=map.pcolor(x,y,np.squeeze(siconc_D015[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D021
cs = map.pcolor(x,y,np.squeeze(siconc_D021[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D022
cs = map.pcolor(x,y,np.squeeze(siconc_D022[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D023
cs = map.pcolor(x,y,np.squeeze(siconc_D023[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7d.png')


# Supp. Fig. 7e (SST+1K) - Maps of difference in September siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month2-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D013
cs = map.pcolor(x,y,np.squeeze(siconc_D013[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+1$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D016
cs = map.pcolor(x,y,np.squeeze(siconc_D016[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+1$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D019
cs=map.pcolor(x,y,np.squeeze(siconc_D019[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+1$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D027
cs = map.pcolor(x,y,np.squeeze(siconc_D027[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+1$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D029
cs = map.pcolor(x,y,np.squeeze(siconc_D029[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+1$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D024
cs = map.pcolor(x,y,np.squeeze(siconc_D024[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7e.png')


# Supp. Fig. 7f (SST+3K) - Maps of difference in September siconc between the SST restoring experiments and CTRL
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,siconc_D000[month2-1,:,:],vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,25,50,75,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC ($\%$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# March D014
cs = map.pcolor(x,y,np.squeeze(siconc_D014[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+5$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# March D017
cs = map.pcolor(x,y,np.squeeze(siconc_D017[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+5$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# March D020
cs=map.pcolor(x,y,np.squeeze(siconc_D020[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+5$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# March D028
cs = map.pcolor(x,y,np.squeeze(siconc_D028[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+5$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# March D030
cs = map.pcolor(x,y,np.squeeze(siconc_D030[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+5$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# March D025
cs = map.pcolor(x,y,np.squeeze(siconc_D025[month2-1,:,:]-siconc_D000[month2-1,:,:]),vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SIC PERT - CTRL ($\%$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig7f.png')
