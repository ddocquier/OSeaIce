#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 4: Map of ocean heat flux in each grid cell
    Ocean heat flux computed via average_ohf.py, which follows compute_oht.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    23/10/2020
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

# Load OHF D013
filename = dir_input + 'D013/OHT_transects/oht_mean_D013.npy'
[oht_u_D013,oht_v_D013] = np.load(filename)
oht_D013 = np.sqrt(oht_u_D013**2. + oht_v_D013**2.)
oht_D013 = oht_D013 / 1.e3

# Load OHF D014
filename = dir_input + 'D014/OHT_transects/oht_mean_D014.npy'
[oht_u_D014,oht_v_D014] = np.load(filename)
oht_D014 = np.sqrt(oht_u_D014**2. + oht_v_D014**2.)
oht_D014 = oht_D014 / 1.e3

# Load OHF D015
filename = dir_input + 'D015/OHT_transects/oht_mean_D015.npy'
[oht_u_D015,oht_v_D015] = np.load(filename)
oht_D015 = np.sqrt(oht_u_D015**2. + oht_v_D015**2.)
oht_D015 = oht_D015 / 1.e3

# Load OHF D016
filename = dir_input + 'D016/OHT_transects/oht_mean_D016.npy'
[oht_u_D016,oht_v_D016] = np.load(filename)
oht_D016 = np.sqrt(oht_u_D016**2. + oht_v_D016**2.)
oht_D016 = oht_D016 / 1.e3

# Load OHF D017
filename = dir_input + 'D017/OHT_transects/oht_mean_D017.npy'
[oht_u_D017,oht_v_D017] = np.load(filename)
oht_D017 = np.sqrt(oht_u_D017**2. + oht_v_D017**2.)
oht_D017 = oht_D017 / 1.e3

# Load OHF D018
filename = dir_input + 'D018/OHT_transects/oht_mean_D018.npy'
[oht_u_D018,oht_v_D018] = np.load(filename)
oht_D018 = np.sqrt(oht_u_D018**2. + oht_v_D018**2.)
oht_D018 = oht_D018 / 1.e3

# Load OHF D019
filename = dir_input + 'D019/OHT_transects/oht_mean_D019.npy'
[oht_u_D019,oht_v_D019] = np.load(filename)
oht_D019 = np.sqrt(oht_u_D019**2. + oht_v_D019**2.)
oht_D019 = oht_D019 / 1.e3

# Load OHF D020
filename = dir_input + 'D020/OHT_transects/oht_mean_D020.npy'
[oht_u_D020,oht_v_D020] = np.load(filename)
oht_D020 = np.sqrt(oht_u_D020**2. + oht_v_D020**2.)
oht_D020 = oht_D020 / 1.e3

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

# Load OHF D024
filename = dir_input + 'D024/OHT_transects/oht_mean_D024.npy'
[oht_u_D024,oht_v_D024] = np.load(filename)
oht_D024 = np.sqrt(oht_u_D024**2. + oht_v_D024**2.)
oht_D024 = oht_D024 / 1.e3

# Load OHF D025
filename = dir_input + 'D025/OHT_transects/oht_mean_D025.npy'
[oht_u_D025,oht_v_D025] = np.load(filename)
oht_D025 = np.sqrt(oht_u_D025**2. + oht_v_D025**2.)
oht_D025 = oht_D025 / 1.e3

# Load OHF D027
filename = dir_input + 'D027/OHT_transects/oht_mean_D027.npy'
[oht_u_D027,oht_v_D027] = np.load(filename)
oht_D027 = np.sqrt(oht_u_D027**2. + oht_v_D027**2.)
oht_D027 = oht_D027 / 1.e3

# Load OHF D028
filename = dir_input + 'D028/OHT_transects/oht_mean_D028.npy'
[oht_u_D028,oht_v_D028] = np.load(filename)
oht_D028 = np.sqrt(oht_u_D028**2. + oht_v_D028**2.)
oht_D028 = oht_D028 / 1.e3

# Load OHF D029
filename = dir_input + 'D029/OHT_transects/oht_mean_D029.npy'
[oht_u_D029,oht_v_D029] = np.load(filename)
oht_D029 = np.sqrt(oht_u_D029**2. + oht_v_D029**2.)
oht_D029 = oht_D029 / 1.e3

# Load OHF D030
filename = dir_input + 'D030/OHT_transects/oht_mean_D030.npy'
[oht_u_D030,oht_v_D030] = np.load(filename)
oht_D030 = np.sqrt(oht_u_D030**2. + oht_v_D030**2.)
oht_D030 = oht_D030 / 1.e3

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

# Load siconc D013
filename = dir_input + 'D013/siconc_D013_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D013 = fh.variables['siconc'][:]
siconc_D013 = siconc_D013 * 100.
fh.close()

# Load siconc D014
filename = dir_input + 'D014/siconc_D014_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D014 = fh.variables['siconc'][:]
siconc_D014 = siconc_D014 * 100.
fh.close()

# Load siconc D015
filename = dir_input + 'D015/siconc_D015_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D015 = fh.variables['siconc'][:]
siconc_D015 = siconc_D015 * 100.
fh.close()

# Load siconc D016
filename = dir_input + 'D016/siconc_D016_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D016 = fh.variables['siconc'][:]
siconc_D016 = siconc_D016 * 100.
fh.close()

# Load siconc D017
filename = dir_input + 'D017/siconc_D017_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D017 = fh.variables['siconc'][:]
siconc_D017 = siconc_D017 * 100.
fh.close()

# Load siconc D018
filename = dir_input + 'D018/siconc_D018_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D018 = fh.variables['siconc'][:]
siconc_D018 = siconc_D018 * 100.
fh.close()

# Load siconc D019
filename = dir_input + 'D019/siconc_D019_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D019 = fh.variables['siconc'][:]
siconc_D019 = siconc_D019 * 100.
fh.close()

# Load siconc D020
filename = dir_input + 'D020/siconc_D020_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D020 = fh.variables['siconc'][:]
siconc_D020 = siconc_D020 * 100.
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

# Load siconc D024
filename = dir_input + 'D024/siconc_D024_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D024 = fh.variables['siconc'][:]
siconc_D024 = siconc_D024 * 100.
fh.close()

# Load siconc D025
filename = dir_input + 'D025/siconc_D025_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D025 = fh.variables['siconc'][:]
siconc_D025 = siconc_D025 * 100.
fh.close()

# Load siconc D027
filename = dir_input + 'D027/siconc_D027_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D027 = fh.variables['siconc'][:]
siconc_D027 = siconc_D027 * 100.
fh.close()

# Load siconc D028
filename = dir_input + 'D028/siconc_D028_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D028 = fh.variables['siconc'][:]
siconc_D028 = siconc_D028 * 100.
fh.close()

# Load siconc D029
filename = dir_input + 'D029/siconc_D029_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D029 = fh.variables['siconc'][:]
siconc_D029 = siconc_D029 * 100.
fh.close()

# Load siconc D030
filename = dir_input + 'D030/siconc_D030_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D030 = fh.variables['siconc'][:]
siconc_D030 = siconc_D030 * 100.
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


# Fig. Supp. 4a (SST+1K) - OHF maps
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,oht_D000,vmin=min_oht,vmax=max_oht,cmap=palette_oht,ax=ax[0,0])
map.contour(x,y,siconc_D000[2,:,:],range(15,16,5),colors='gray',ax=ax[0,0],linewidths=3)
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,200,400,600,800],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF (kW m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D013
cs=map.pcolor(x,y,oht_D013-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,siconc_D013[2,:,:],range(15,16,5),colors='gray',ax=ax[1,0],linewidths=3)
map.contour(x,y,siconc_D013[8,:,:],range(15,16,5),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+1$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D016
cs=map.pcolor(x,y,oht_D016-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,siconc_D016[2,:,:],range(15,16,5),colors='gray',ax=ax[1,1],linewidths=3)
map.contour(x,y,siconc_D016[8,:,:],range(15,16,5),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+1$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D019
cs=map.pcolor(x,y,oht_D019-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,siconc_D019[2,:,:],range(15,16,5),colors='gray',ax=ax[1,2],linewidths=3)
map.contour(x,y,siconc_D019[8,:,:],range(15,16,5),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+1$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D027
cs=map.pcolor(x,y,oht_D027-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,siconc_D027[2,:,:],range(15,16,5),colors='gray',ax=ax[2,0],linewidths=3)
map.contour(x,y,siconc_D027[8,:,:],range(15,16,5),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+1$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D029
cs=map.pcolor(x,y,oht_D029-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,siconc_D029[2,:,:],range(15,16,5),colors='gray',ax=ax[2,1],linewidths=3)
map.contour(x,y,siconc_D029[8,:,:],range(15,16,5),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+1$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D024
cs=map.pcolor(x,y,oht_D024-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,siconc_D024[2,:,:],range(15,16,5),colors='gray',ax=ax[2,2],linewidths=3)
map.contour(x,y,siconc_D024[8,:,:],range(15,16,5),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+1$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-400,-200,0,200,400],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF PERT - CTRL (kW m$^{-2}$)',fontsize=28)

# Save figure
#if save_fig == True:
#    fig.savefig(dir_output + 'fig4a.png')


# Fig. Supp. 4b (SST+5K) - OHF maps
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,oht_D000,vmin=min_oht,vmax=max_oht,cmap=palette_oht,ax=ax[0,0])
map.contour(x,y,siconc_D000[2,:,:],range(15,16,5),colors='gray',ax=ax[0,0],linewidths=3)
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
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[0,200,400,600,800],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF (kW m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D014
cs=map.pcolor(x,y,oht_D014-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,siconc_D014[2,:,:],range(15,16,5),colors='gray',ax=ax[1,0],linewidths=3)
map.contour(x,y,siconc_D014[8,:,:],range(15,16,5),colors='black',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+5$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D017
cs=map.pcolor(x,y,oht_D017-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,siconc_D017[2,:,:],range(15,16,5),colors='gray',ax=ax[1,1],linewidths=3)
map.contour(x,y,siconc_D017[8,:,:],range(15,16,5),colors='black',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+5$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D020
cs=map.pcolor(x,y,oht_D020-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,siconc_D020[2,:,:],range(15,16,5),colors='gray',ax=ax[1,2],linewidths=3)
map.contour(x,y,siconc_D020[8,:,:],range(15,16,5),colors='black',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+5$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D028
cs=map.pcolor(x,y,oht_D028-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,siconc_D028[2,:,:],range(15,16,5),colors='gray',ax=ax[2,0],linewidths=3)
map.contour(x,y,siconc_D028[8,:,:],range(15,16,5),colors='black',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+5$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D030
cs=map.pcolor(x,y,oht_D030-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,siconc_D030[2,:,:],range(15,16,5),colors='gray',ax=ax[2,1],linewidths=3)
map.contour(x,y,siconc_D030[8,:,:],range(15,16,5),colors='black',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+5$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D025
cs=map.pcolor(x,y,oht_D025-oht_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,siconc_D025[2,:,:],range(15,16,5),colors='gray',ax=ax[2,2],linewidths=3)
map.contour(x,y,siconc_D025[8,:,:],range(15,16,5),colors='black',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
ax[2,2].set_title('PAC3+5$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-400,-200,0,200,400],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Horizontal OHF PERT - CTRL (kW m$^{-2}$)',fontsize=28)

# Save figure
#if save_fig == True:
#    fig.savefig(dir_output + 'fig4b.png')
