#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 4: Map of ocean heat flux in each grid cell
    Ocean heat flux computed via average_ohf.py, which follows compute_oht.py
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

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Options
save_fig = True

# Load OHF D000
filename = dir_input + 'D000/OHT_transects/oht_mean_D000.npy'
[oht_u_D000,oht_v_D000] = np.load(filename)
oht_D000 = np.sqrt(oht_u_D000**2. + oht_v_D000**2.)
oht_D000 = oht_D000 / 1.e3

# Load OHF D012
filename = dir_input + 'D012/OHT_transects/oht_mean_D012.npy'
[oht_u_D012,oht_v_D012] = np.load(filename)
oht_D012 = np.sqrt(oht_u_D012**2. + oht_v_D012**2.)
oht_D012 = oht_D012 / 1.e3

# Load OHF D015
filename = dir_input + 'D015/OHT_transects/oht_mean_D015.npy'
[oht_u_D015,oht_v_D015] = np.load(filename)
oht_D015 = np.sqrt(oht_u_D015**2. + oht_v_D015**2.)
oht_D015 = oht_D015 / 1.e3

# Load OHF D018
filename = dir_input + 'D018/OHT_transects/oht_mean_D018.npy'
[oht_u_D018,oht_v_D018] = np.load(filename)
oht_D018 = np.sqrt(oht_u_D018**2. + oht_v_D018**2.)
oht_D018 = oht_D018 / 1.e3

# Load OHF D021
filename = dir_input + 'D021/OHT_transects/oht_mean_D021.npy'
[oht_u_D021,oht_v_D021] = np.load(filename)
oht_D021 = np.sqrt(oht_u_D021**2. + oht_v_D021**2.)
oht_D021 = oht_D021 / 1.e3

# Load OHF D022
filename = dir_input + 'D022/OHT_transects/oht_mean_D022.npy'
[oht_u_D022,oht_v_D022] = np.load(filename)
oht_D022 = np.sqrt(oht_u_D022**2. + oht_v_D022**2.)
oht_D022 = oht_D022 / 1.e3

# Load OHF D023
filename = dir_input + 'D023/OHT_transects/oht_mean_D023.npy'
[oht_u_D023,oht_v_D023] = np.load(filename)
oht_D023 = np.sqrt(oht_u_D023**2. + oht_v_D023**2.)
oht_D023 = oht_D023 / 1.e3

# Load siconc D000
filename = dir_input + 'D000/siconc_D000_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D000 = fh.variables['siconc'][:]
siconc_D000 = siconc_D000 * 100.
nm,ny,nx = siconc_D000.shape
fh.close()

# Load siconc D012
filename = dir_input + 'D012/siconc_D012_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D012 = fh.variables['siconc'][:]
siconc_D012 = siconc_D012 * 100.
fh.close()

# Load siconc D015
filename = dir_input + 'D015/siconc_D015_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D015 = fh.variables['siconc'][:]
siconc_D015 = siconc_D015 * 100.
fh.close()

# Load siconc D018
filename = dir_input + 'D018/siconc_D018_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D018 = fh.variables['siconc'][:]
siconc_D018 = siconc_D018 * 100.
fh.close()

# Load siconc D021
filename = dir_input + 'D021/siconc_D021_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D021 = fh.variables['siconc'][:]
siconc_D021 = siconc_D021 * 100.
fh.close()

# Load siconc D022
filename = dir_input + 'D022/siconc_D022_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D022 = fh.variables['siconc'][:]
siconc_D022 = siconc_D022 * 100.
fh.close()

# Load siconc D023
filename = dir_input + 'D023/siconc_D023_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D023 = fh.variables['siconc'][:]
siconc_D023 = siconc_D023 * 100.
fh.close()

# Load latitude and longitude of model
filename = dir_grid + 'mesh_zgr.nc'
fh = Dataset(filename, mode='r')
lon = fh.variables['nav_lon'][:]
lat = fh.variables['nav_lat'][:]
fh.close()

# Map projection
boundlat = 40.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)

# Palettes
palette_oht = plt.cm.cubehelix_r._resample(20)
min_oht = 0.
max_oht = 1000.
palette_diff = plt.cm.seismic._resample(40)
min_diff = -400.
max_diff = 400.

    
# Fig. 4 - OHF maps
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,oht_D000,vmin=min_oht,vmax=max_oht,cmap=palette_oht,ax=ax[0,0])
map.contour(x,y,siconc_D000[2,:,:],range(15,16,5),colors='blue',ax=ax[0,0],linewidths=3)
map.contour(x,y,siconc_D000[8,:,:],range(15,16,5),colors='black',ax=ax[0,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,250,500,750,1000],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF (kW m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
cs=map.pcolor(x,y,oht_D012-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,siconc_D012[2,:,:],range(15,16,5),colors='blue',ax=ax[1,0],linewidths=3)
map.contour(x,y,siconc_D012[8,:,:],range(15,16,5),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,oht_D015-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,siconc_D015[2,:,:],range(15,16,5),colors='blue',ax=ax[1,1],linewidths=3)
map.contour(x,y,siconc_D015[8,:,:],range(15,16,5),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,oht_D018-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,siconc_D018[2,:,:],range(15,16,5),colors='blue',ax=ax[1,2],linewidths=3)
map.contour(x,y,siconc_D018[8,:,:],range(15,16,5),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,oht_D021-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,siconc_D021[2,:,:],range(15,16,5),colors='blue',ax=ax[2,0],linewidths=3)
map.contour(x,y,siconc_D021[8,:,:],range(15,16,5),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,oht_D022-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,siconc_D022[2,:,:],range(15,16,5),colors='blue',ax=ax[2,1],linewidths=3)
map.contour(x,y,siconc_D022[8,:,:],range(15,16,5),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,oht_D023-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,siconc_D023[2,:,:],range(15,16,5),colors='blue',ax=ax[2,2],linewidths=3)
map.contour(x,y,siconc_D023[8,:,:],range(15,16,5),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-400,-200,0,200,400],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF PERT - CTRL (kW m$^{-2}$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig4.png')
    fig.savefig(dir_output + 'fig4.eps',dpi=300)
