#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 15: Plot surface fluxes from SST experiments (IFS outputs)
    SHF files saved with save_surface_fluxes.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    17/10/2020
'''

# Options
save_fig = True

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'
dir_mask = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/surface_restoring/ifs/'

# Load D000
filename = dir_input + 'D000/SurfaceFluxes_D000.npy'
swrad_D000,lwrad_D000,shflux_D000,lhflux_D000,Rts_D000,Rts_spatial_D000,lat_init,lon_init = np.load(filename,allow_pickle=True)
lon,lat = np.meshgrid(lon_init,lat_init)
lon[lon>180.] = lon[lon>180.] - 360.

# Load D012
filename = dir_input + 'D012/SurfaceFluxes_D012.npy'
swrad_D012,lwrad_D012,shflux_D012,lhflux_D012,Rts_D012,Rts_spatial_D012,notused,notused = np.load(filename,allow_pickle=True)

# Load D015
filename = dir_input + 'D015/SurfaceFluxes_D015.npy'
swrad_D015,lwrad_D015,shflux_D015,lhflux_D015,Rts_D015,Rts_spatial_D015,notused,notused = np.load(filename,allow_pickle=True)

# Load D018
filename = dir_input + 'D018/SurfaceFluxes_D018.npy'
swrad_D018,lwrad_D018,shflux_D018,lhflux_D018,Rts_D018,Rts_spatial_D018,notused,notused = np.load(filename,allow_pickle=True)

# Load D021
filename = dir_input + 'D021/SurfaceFluxes_D021.npy'
swrad_D021,lwrad_D021,shflux_D021,lhflux_D021,Rts_D021,Rts_spatial_D021,notused,notused = np.load(filename,allow_pickle=True)

# Load D022
filename = dir_input + 'D022/SurfaceFluxes_D022.npy'
swrad_D022,lwrad_D022,shflux_D022,lhflux_D022,Rts_D022,Rts_spatial_D022,notused,notused = np.load(filename,allow_pickle=True)

# Load D023
filename = dir_input + 'D023/SurfaceFluxes_D023.npy'
swrad_D023,lwrad_D023,shflux_D023,lhflux_D023,Rts_D023,Rts_spatial_D023,notused,notused = np.load(filename,allow_pickle=True)

# Load mask ATL1
filename = dir_mask + 'mask_restore_ATL1.nc'
fh = Dataset(filename, mode='r')
mask_atl1 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask ATL2
filename = dir_mask + 'mask_restore_ATL2.nc'
fh = Dataset(filename, mode='r')
mask_atl2 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask ATL3
filename = dir_mask + 'mask_restore_ATL3.nc'
fh = Dataset(filename, mode='r')
mask_atl3 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC1
filename = dir_mask + 'mask_restore_PAC1.nc'
fh = Dataset(filename, mode='r')
mask_pac1 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC2
filename = dir_mask + 'mask_restore_PAC2.nc'
fh = Dataset(filename, mode='r')
mask_pac2 = fh.variables['mask_ssr'][:]
fh.close()

# Load mask PAC3
filename = dir_mask + 'mask_restore_PAC3.nc'
fh = Dataset(filename, mode='r')
mask_pac3 = fh.variables['mask_ssr'][:]
fh.close()

# Ocean mask
filename = dir_mask + 'mask_ocean_ifs.nc'
fh = Dataset(filename, mode='r')
mask_ocean = fh.variables['glomsk'][:]
fh.close()

# Define Atlantic and Pacific domains for SHF
mask_atlantic = np.zeros((np.size(Rts_spatial_D000,0),np.size(Rts_spatial_D000,1)))
mask_pacific = np.zeros((np.size(Rts_spatial_D000,0),np.size(Rts_spatial_D000,1)))
mask_atlantic[(lat>=50.)*(lat<=80.)*(lon>=-40.)*(lon<=20.)*(mask_ocean==1)] = 1
mask_pacific[(lat>=45.)*(lat<=66.)*(lon>=150.)*(lon<=180.)*(mask_ocean==1)] = 1
mask_pacific[(lat>=45.)*(lat<=66.)*(lon>=-180.)*(lon<=-150.)*(mask_ocean==1)] = 1

# Map projection
boundlat = 40.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)

# Palettes and plot parameters
palette_var = plt.cm.seismic._resample(20)
min_var = -100.
max_var = 100.
palette_diff = plt.cm.seismic._resample(20)
min_diff = -50.
max_diff = 50.


# Figure - Zonal means
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.1,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)
ticks = [-90,-60,-30,0,30,60,90]

# D000
ax[0,0].plot(lat_init,np.nanmean(swrad_D000,axis=0),'r-',label='Shortwave radiation')
ax[0,0].plot(lat_init,np.nanmean(lwrad_D000,axis=0),'b-',label='Longwave radiation')
ax[0,0].plot(lat_init,np.nanmean(shflux_D000,axis=0),'m--',label='Sensible heat flux')
ax[0,0].plot(lat_init,np.nanmean(lhflux_D000,axis=0),'c--',label='Latent heat flux')
ax[0,0].plot(lat_init,np.nanmean(Rts_D000,axis=0),'k-',label='Net heat flux')
ax[0,0].set_ylabel('Surface heat flux (W m$^{-2})$',fontsize=20)
ax[0,0].axis([-90,90,-150,300])
ax[0,0].set_xticks(ticks);
ax[0,0].legend(shadow=True,frameon=False,fontsize=20,bbox_to_anchor=(1,1))
ax[0,0].grid(linestyle='--')
ax[0,0].tick_params(axis='both',labelsize=18)
ax[0,0].set_title('CTRL',fontsize=20)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012-D000
ax[1,0].plot(lat_init,np.nanmean(swrad_D012-swrad_D000,axis=0),'r-')
ax[1,0].plot(lat_init,np.nanmean(lwrad_D012-lwrad_D000,axis=0),'b-')
ax[1,0].plot(lat_init,np.nanmean(shflux_D012-shflux_D000,axis=0),'m--')
ax[1,0].plot(lat_init,np.nanmean(lhflux_D012-lhflux_D000,axis=0),'c--')
ax[1,0].plot(lat_init,np.nanmean(Rts_D012-Rts_D000,axis=0),'k-')
ax[1,0].set_ylabel('Change in surface heat flux \n (W m$^{-2})$',fontsize=20)
ax[1,0].axis([-90,90,-10,10])
ax[1,0].set_xticks(ticks);
ax[1,0].grid(linestyle='--')
ax[1,0].tick_params(axis='both',labelsize=18)
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=20)

# D015-D000
ax[1,1].plot(lat_init,np.nanmean(swrad_D015-swrad_D000,axis=0),'r-')
ax[1,1].plot(lat_init,np.nanmean(lwrad_D015-lwrad_D000,axis=0),'b-')
ax[1,1].plot(lat_init,np.nanmean(shflux_D015-shflux_D000,axis=0),'m--')
ax[1,1].plot(lat_init,np.nanmean(lhflux_D015-lhflux_D000,axis=0),'c--')
ax[1,1].plot(lat_init,np.nanmean(Rts_D015-Rts_D000,axis=0),'k-')
ax[1,1].axis([-90,90,-10,10])
ax[1,1].set_xticks(ticks);
ax[1,1].grid(linestyle='--')
ax[1,1].tick_params(axis='both',labelsize=18)
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=20)

# D018-D000
ax[1,2].plot(lat_init,np.nanmean(swrad_D018-swrad_D000,axis=0),'r-')
ax[1,2].plot(lat_init,np.nanmean(lwrad_D018-lwrad_D000,axis=0),'b-')
ax[1,2].plot(lat_init,np.nanmean(shflux_D018-shflux_D000,axis=0),'m--')
ax[1,2].plot(lat_init,np.nanmean(lhflux_D018-lhflux_D000,axis=0),'c--')
ax[1,2].plot(lat_init,np.nanmean(Rts_D018-Rts_D000,axis=0),'k-')
ax[1,2].axis([-90,90,-10,10])
ax[1,2].set_xticks(ticks);
ax[1,2].grid(linestyle='--')
ax[1,2].tick_params(axis='both',labelsize=18)
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=20)

# D021-D000
ax[2,0].plot(lat_init,np.nanmean(swrad_D021-swrad_D000,axis=0),'r-')
ax[2,0].plot(lat_init,np.nanmean(lwrad_D021-lwrad_D000,axis=0),'b-')
ax[2,0].plot(lat_init,np.nanmean(shflux_D021-shflux_D000,axis=0),'m--')
ax[2,0].plot(lat_init,np.nanmean(lhflux_D021-lhflux_D000,axis=0),'c--')
ax[2,0].plot(lat_init,np.nanmean(Rts_D021-Rts_D000,axis=0),'k-')
ax[2,0].set_ylabel('Change in surface heat flux \n (W m$^{-2})$',fontsize=20)
ax[2,0].set_xlabel('Latitude ($^{\circ}$)',fontsize=20)
ax[2,0].axis([-90,90,-10,10])
ax[2,0].set_xticks(ticks);
ax[2,0].grid(linestyle='--')
ax[2,0].tick_params(axis='both',labelsize=18)
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=20)

# D022-D000
ax[2,1].plot(lat_init,np.nanmean(swrad_D022-swrad_D000,axis=0),'r-')
ax[2,1].plot(lat_init,np.nanmean(lwrad_D022-lwrad_D000,axis=0),'b-')
ax[2,1].plot(lat_init,np.nanmean(shflux_D022-shflux_D000,axis=0),'m--')
ax[2,1].plot(lat_init,np.nanmean(lhflux_D022-lhflux_D000,axis=0),'c--')
ax[2,1].plot(lat_init,np.nanmean(Rts_D022-Rts_D000,axis=0),'k-')
ax[2,1].set_xlabel('Latitude ($^{\circ}$)',fontsize=20)
ax[2,1].axis([-90,90,-10,10])
ax[2,1].set_xticks(ticks);
ax[2,1].grid(linestyle='--')
ax[2,1].tick_params(axis='both',labelsize=18)
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=20)

# D015-D000
ax[2,2].plot(lat_init,np.nanmean(swrad_D023-swrad_D000,axis=0),'r-')
ax[2,2].plot(lat_init,np.nanmean(lwrad_D023-lwrad_D000,axis=0),'b-')
ax[2,2].plot(lat_init,np.nanmean(shflux_D023-shflux_D000,axis=0),'m--')
ax[2,2].plot(lat_init,np.nanmean(lhflux_D023-lhflux_D000,axis=0),'c--')
ax[2,2].plot(lat_init,np.nanmean(Rts_D023-Rts_D000,axis=0),'k-')
ax[2,2].set_xlabel('Latitude ($^{\circ}$)',fontsize=20)
ax[2,2].axis([-90,90,-10,10])
ax[2,2].set_xticks(ticks);
ax[2,2].grid(linestyle='--')
ax[2,2].tick_params(axis='both',labelsize=18)
ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=20)

# Save figure
#if save_fig == True:
#    fig.savefig(dir_output + 'LatitudinalTransects_SurfaceFluxes.png')


# Figure - Maps 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,Rts_spatial_D000,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-100,-50,0,50,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Net surface heat flux (W m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
cs=map.pcolor(x,y,Rts_spatial_D012-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.contour(x,y,mask_atl1,range(1,2,1),colors='black',ax=ax[1,0],linewidths=3)
map.contour(x,y,mask_pacific,range(1,2,1),colors='green',ax=ax[1,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,Rts_spatial_D015-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.contour(x,y,mask_atl2,range(1,2,1),colors='black',ax=ax[1,1],linewidths=3)
map.contour(x,y,mask_pacific,range(1,2,1),colors='green',ax=ax[1,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,Rts_spatial_D018-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.contour(x,y,mask_atl3,range(1,2,1),colors='black',ax=ax[1,2],linewidths=3)
map.contour(x,y,mask_pacific,range(1,2,1),colors='green',ax=ax[1,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,Rts_spatial_D021-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.contour(x,y,mask_pac1,range(1,2,1),colors='black',ax=ax[2,0],linewidths=3)
map.contour(x,y,mask_atlantic,range(1,2,1),colors='green',ax=ax[2,0],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,Rts_spatial_D022-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.contour(x,y,mask_pac2,range(1,2,1),colors='black',ax=ax[2,1],linewidths=3)
map.contour(x,y,mask_atlantic,range(1,2,1),colors='green',ax=ax[2,1],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,Rts_spatial_D023-Rts_spatial_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
map.contour(x,y,mask_pac3,range(1,2,1),colors='black',ax=ax[2,2],linewidths=3)
map.contour(x,y,mask_atlantic,range(1,2,1),colors='green',ax=ax[2,2],linewidths=3)
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
map.drawcoastlines(ax=ax[2,2])
ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=32)
ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
ax[2,2].yaxis.set_label_coords(-0.05,0.9)

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-50,-25,0,25,50],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Net SHF PERT - CTRL (W m$^{-2}$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig15.png')
    fig.savefig(dir_output + 'fig15.eps',dpi=300)
