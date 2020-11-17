#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig.1: Map of SST
PROGRAMMER
    D. Docquier
LAST UPDATE
    17/11/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec

# Working directories
dir_D000 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D000/' # CTRL
dir_D012 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D012/' # ATL1+3K
dir_D015 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D015/' # ATL3+3K
dir_D018 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D018/' # ATL5+3K
dir_D021 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D021/' # PAC1+3K
dir_D022 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D022/' # PAC2+3K
dir_D023 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D023/' # PAC3+3K
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'
dir_mask = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/surface_restoring/'

# Parameters
save_fig = True

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

# Load SST D015
filename = dir_D015 + 'tos_D015_ym.nc'
fh = Dataset(filename, mode='r')
tos_D015 = fh.variables['tos'][:]
fh.close()
tos_D015_annmean = np.nanmean(tos_D015,axis=0)

# Load SST D018
filename = dir_D018 + 'tos_D018_ym.nc'
fh = Dataset(filename, mode='r')
tos_D018 = fh.variables['tos'][:]
fh.close()
tos_D018_annmean = np.nanmean(tos_D018,axis=0)

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
boundlat = 30.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)
boundlat2 = 40.
map2 = Basemap(projection='nplaea',boundinglat=boundlat2,lon_0=l0,resolution='c')
x2,y2 = map2(lon,lat)

# Palettes and plot parameters
palette_var = plt.cm.Reds._resample(20)
min_var = 0.
max_var = 30.
palette_diff = plt.cm.seismic._resample(40)
min_diff = -3.
max_diff = 3.

# Coordinates arrows
lat_arrow1 = 73.75
lon_arrow1 = 15.
x_arrow1,y_arrow1 = map2(lon_arrow1,lat_arrow1) # BSO
x2_arrow1,y2_arrow1 = map2(lon_arrow1+20,lat_arrow1)
lat_arrow2 = 77.
lon_arrow2 = -4.5
x_arrow2,y_arrow2 = map2(lon_arrow2,lat_arrow2) # Fram Strait
x2_arrow2,y2_arrow2 = map2(lon_arrow2,lat_arrow2+5)
lat_arrow3 = 66.
lon_arrow3 = -59.
x_arrow3,y_arrow3 = map2(lon_arrow3,lat_arrow3) # Davis Strait
x2_arrow3,y2_arrow3 = map2(lon_arrow3,lat_arrow3+5)
lat_arrow4 = 63.
lon_arrow4 = -169.25
x_arrow4,y_arrow4 = map2(lon_arrow4,lat_arrow4) # Bering Strait
x2_arrow4,y2_arrow4 = map2(lon_arrow4,lat_arrow4+5)

# Coordinates texts
lat_text1 = 69.
lon_text1 = 5.
x_text1,y_text1 = map2(lon_text1,lat_text1) # BSO
lat_text2 = 72.
lon_text2 = -10.
x_text2,y_text2 = map2(lon_text2,lat_text2) # Fram Strait
lat_text3 = 61.
lon_text3 = -58.
x_text3,y_text3 = map2(lon_text3,lat_text3) # Davis Strait
lat_text4 = 60.
lon_text4 = -165.
x_text4,y_text4 = map2(lon_text4,lat_text4) # Bering Strait


# Fig. 1 - SST+3K maps - 50-year average
fig,ax=plt.subplots(4,3,figsize=(18,24))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
gs = gridspec.GridSpec(4,4)
axbig = fig.add_subplot(gs[0:2,0:2])
cs=map2.pcolor(x2,y2,tos_D000_annmean,vmin=min_var,vmax=max_var,cmap=palette_var,ax=axbig)
map2.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=axbig)
map2.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=axbig)
map2.drawcoastlines(ax=axbig)
map2.fillcontinents(color='grey',lake_color='w',ax=axbig)
axbig.set_title('CTRL',fontsize=32)
axbig.set_title('a',loc='left',fontsize=32,fontweight='bold')
axbig.yaxis.set_label_coords(-0.05,0.9)
plt.arrow(x_arrow1,y_arrow1,x2_arrow1-x_arrow1,y2_arrow1-y_arrow1,fc='k',ec='k',linewidth=10,head_width=50000,head_length=50000)
plt.text(x_text1,y_text1,'1',fontsize=24)
plt.arrow(x_arrow2,y_arrow2,x2_arrow2-x_arrow2,y2_arrow2-y_arrow2,fc='k',ec='k',linewidth=10,head_width=30000,head_length=30000)
plt.text(x_text2,y_text2,'2',fontsize=24)
plt.arrow(x_arrow3,y_arrow3,x2_arrow3-x_arrow3,y2_arrow3-y_arrow3,fc='k',ec='k',linewidth=10,head_width=20000,head_length=20000)
plt.text(x_text3,y_text3,'3',fontsize=24)
plt.arrow(x_arrow4,y_arrow4,x2_arrow4-x_arrow4,y2_arrow4-y_arrow4,fc='k',ec='k',linewidth=10,head_width=20000,head_length=20000)
plt.text(x_text4,y_text4,'4',fontsize=24)

# Add color bar absolute value
cb_ax = fig.add_axes([0.52, 0.6, 0.02, 0.28])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,5,10,15,20,25,30],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST ($^\circ$C)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,0])
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])
fig.delaxes(ax[1,0])
fig.delaxes(ax[1,1])
fig.delaxes(ax[1,2])

# D012
cs=map.pcolor(x,y,tos_D012_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,tos_D015_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,tos_D018_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[2,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,tos_D021_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[3,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[3,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[3,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[3,0])
map.drawcoastlines(ax=ax[3,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[3,0])
ax[3,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[3,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[3,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,tos_D022_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[3,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[3,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[3,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[3,1])
map.drawcoastlines(ax=ax[3,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[3,1])
ax[3,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[3,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[3,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,tos_D023_annmean-tos_D000_annmean,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[3,2])
map.contour(x,y,mask_pac3,range(1,2,1),colors='black',ax=ax[3,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[3,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[3,2])
map.drawcoastlines(ax=ax[3,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[3,2])
ax[3,2].set_title('PAC3+3$^{\circ}$C',fontsize=32)
ax[3,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[3,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.62, 0.6, 0.33, 0.015])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('SST PERT - CTRL ($^\circ$C)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig1.png')
    fig.savefig(dir_output + 'fig1.eps',dpi=300)
