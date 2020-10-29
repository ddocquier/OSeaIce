#!/usr/bin/env python

'''
GOAL
    Plot and save surface fluxes from SST experiments (NEMO outputs)
PROGRAMMER
    D. Docquier
LAST UPDATE
    16/10/2020
'''

# Options
save_fig = False
save_var = True

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Function to compute total net surface heat flux Atlantic
def compute_hflux_atl(latitude,longitude,qt,cell_area,mask):
    lat1 = 50.
    lat2 = 80.
    lon1 = -40.
    lon2 = 20. 
    qt_tot = np.nansum(qt*(latitude>=lat1)*(latitude<=lat2)*(longitude>=lon1)*(longitude<=lon2)*cell_area*(mask==0)) / 1.e12 # TW
    area_tot = np.nansum((latitude>=lat1)*(latitude<=lat2)*(longitude>=lon1)*(longitude<=lon2)*cell_area*(mask==0)) / 1.e12 # 1e6 km2
    return qt_tot,area_tot

# Function to compute total net surface heat flux Pacific
def compute_hflux_pac(latitude,longitude,qt,cell_area,mask):
    lat1 = 45.
    lat2 = 66.
    lon1 = 150.
    lon2 = 180.
    lon3 = -180.
    lon4 = -150.
    qt_tot1 = np.nansum(qt*(latitude>=lat1)*(latitude<=lat2)*(longitude>=lon1)*(longitude<=lon2)*cell_area*(mask==0)) / 1.e12
    qt_tot2 = np.nansum(qt*(latitude>=lat1)*(latitude<=lat2)*(longitude>=lon3)*(longitude<=lon4)*cell_area*(mask==0)) / 1.e12
    qt_tot = qt_tot1 + qt_tot2 # TW
    area_tot1 = np.nansum((latitude>=lat1)*(latitude<=lat2)*(longitude>=lon1)*(longitude<=lon2)*cell_area*(mask==0)) / 1.e12
    area_tot2 = np.nansum((latitude>=lat1)*(latitude<=lat2)*(longitude>=lon3)*(longitude<=lon4)*cell_area*(mask==0)) / 1.e12
    area_tot = area_tot1 + area_tot2 # 1e6 km2
    return qt_tot,area_tot

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/' 

# Load net surface heat flux at ocean surface D000
filename = dir_input + '/D000/qtoce_D000_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D000 = fh.variables['qt_oce'][0,:,:]
grid_area = fh.variables['cell_area'][:]
lat = fh.variables['nav_lat'][:]
lon = fh.variables['nav_lon'][:]
fh.close()

# Load net surface heat flux at ocean surface D012
filename = dir_input + '/D012/qtoce_D012_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D012 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D013
filename = dir_input + '/D013/qtoce_D013_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D013 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D014
filename = dir_input + '/D014/qtoce_D014_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D014 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D015
filename = dir_input + '/D015/qtoce_D015_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D015 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D016
filename = dir_input + '/D016/qtoce_D016_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D016 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D017
filename = dir_input + '/D017/qtoce_D017_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D017 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D018
filename = dir_input + '/D018/qtoce_D018_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D018 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D019
filename = dir_input + '/D019/qtoce_D019_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D019 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D020
filename = dir_input + '/D020/qtoce_D020_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D020 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D021
filename = dir_input + '/D021/qtoce_D021_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D021 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D022
filename = dir_input + '/D022/qtoce_D022_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D022 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D023
filename = dir_input + '/D023/qtoce_D023_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D023 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D024
filename = dir_input + '/D024/qtoce_D024_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D024 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D025
filename = dir_input + '/D025/qtoce_D025_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D025 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D027
filename = dir_input + '/D027/qtoce_D027_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D027 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D028
filename = dir_input + '/D028/qtoce_D028_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D028 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D029
filename = dir_input + '/D029/qtoce_D029_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D029 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load net surface heat flux at ocean surface D030
filename = dir_input + '/D030/qtoce_D030_mean.nc'
fh = Dataset(filename, mode='r')
qtoce_D030 = fh.variables['qt_oce'][0,:,:]
fh.close()

# Load ocean mask EC-Earth
file_mask = dir_grid + 'masks.nc'
fh = Dataset(file_mask,mode='r')
mask_ocean = fh.variables['O1t0.msk'][:]
fh.close()

# Compute total surface heat fluxes in the Atlantic and convert into TW
qtoce_atl_D000,area_atl = compute_hflux_atl(lat,lon,qtoce_D000,grid_area,mask_ocean)
qtoce_atl_D012,notused = compute_hflux_atl(lat,lon,qtoce_D012,grid_area,mask_ocean)
qtoce_atl_D013,notused = compute_hflux_atl(lat,lon,qtoce_D013,grid_area,mask_ocean)
qtoce_atl_D014,notused = compute_hflux_atl(lat,lon,qtoce_D014,grid_area,mask_ocean)
qtoce_atl_D015,notused = compute_hflux_atl(lat,lon,qtoce_D015,grid_area,mask_ocean)
qtoce_atl_D016,notused = compute_hflux_atl(lat,lon,qtoce_D016,grid_area,mask_ocean)
qtoce_atl_D017,notused = compute_hflux_atl(lat,lon,qtoce_D017,grid_area,mask_ocean)
qtoce_atl_D018,notused = compute_hflux_atl(lat,lon,qtoce_D018,grid_area,mask_ocean)
qtoce_atl_D019,notused = compute_hflux_atl(lat,lon,qtoce_D019,grid_area,mask_ocean)
qtoce_atl_D020,notused = compute_hflux_atl(lat,lon,qtoce_D020,grid_area,mask_ocean)
qtoce_atl_D021,notused = compute_hflux_atl(lat,lon,qtoce_D021,grid_area,mask_ocean)
qtoce_atl_D022,notused = compute_hflux_atl(lat,lon,qtoce_D022,grid_area,mask_ocean)
qtoce_atl_D023,notused = compute_hflux_atl(lat,lon,qtoce_D023,grid_area,mask_ocean)
qtoce_atl_D024,notused = compute_hflux_atl(lat,lon,qtoce_D024,grid_area,mask_ocean)
qtoce_atl_D025,notused = compute_hflux_atl(lat,lon,qtoce_D025,grid_area,mask_ocean)
qtoce_atl_D027,notused = compute_hflux_atl(lat,lon,qtoce_D027,grid_area,mask_ocean)
qtoce_atl_D028,notused = compute_hflux_atl(lat,lon,qtoce_D028,grid_area,mask_ocean)
qtoce_atl_D029,notused = compute_hflux_atl(lat,lon,qtoce_D029,grid_area,mask_ocean)
qtoce_atl_D030,notused = compute_hflux_atl(lat,lon,qtoce_D030,grid_area,mask_ocean)
print('Total area Atlantic (1e6 km2):',area_atl)
print('Total net surface heat flux Atlantic CTRL (TW):',qtoce_atl_D000)
print('Total net surface heat flux Atlantic ATL1+3K (TW):',qtoce_atl_D012,qtoce_atl_D012-qtoce_atl_D000)
print('Total net surface heat flux Atlantic PAC1+3K (TW):',qtoce_atl_D021,qtoce_atl_D021-qtoce_atl_D000)
print('Total net surface heat flux Atlantic PAC2+3K (TW):',qtoce_atl_D022,qtoce_atl_D022-qtoce_atl_D000)

# Compute total surface heat fluxes in the Pacific and convert into TW
qtoce_pac_D000,area_pac = compute_hflux_pac(lat,lon,qtoce_D000,grid_area,mask_ocean)
qtoce_pac_D012,notused = compute_hflux_pac(lat,lon,qtoce_D012,grid_area,mask_ocean)
qtoce_pac_D013,notused = compute_hflux_pac(lat,lon,qtoce_D013,grid_area,mask_ocean)
qtoce_pac_D014,notused = compute_hflux_pac(lat,lon,qtoce_D014,grid_area,mask_ocean)
qtoce_pac_D015,notused = compute_hflux_pac(lat,lon,qtoce_D015,grid_area,mask_ocean)
qtoce_pac_D016,notused = compute_hflux_pac(lat,lon,qtoce_D016,grid_area,mask_ocean)
qtoce_pac_D017,notused = compute_hflux_pac(lat,lon,qtoce_D017,grid_area,mask_ocean)
qtoce_pac_D018,notused = compute_hflux_pac(lat,lon,qtoce_D018,grid_area,mask_ocean)
qtoce_pac_D019,notused = compute_hflux_pac(lat,lon,qtoce_D019,grid_area,mask_ocean)
qtoce_pac_D020,notused = compute_hflux_pac(lat,lon,qtoce_D020,grid_area,mask_ocean)
qtoce_pac_D021,notused = compute_hflux_pac(lat,lon,qtoce_D021,grid_area,mask_ocean)
qtoce_pac_D022,notused = compute_hflux_pac(lat,lon,qtoce_D022,grid_area,mask_ocean)
qtoce_pac_D023,notused = compute_hflux_pac(lat,lon,qtoce_D023,grid_area,mask_ocean)
qtoce_pac_D024,notused = compute_hflux_pac(lat,lon,qtoce_D024,grid_area,mask_ocean)
qtoce_pac_D025,notused = compute_hflux_pac(lat,lon,qtoce_D025,grid_area,mask_ocean)
qtoce_pac_D027,notused = compute_hflux_pac(lat,lon,qtoce_D027,grid_area,mask_ocean)
qtoce_pac_D028,notused = compute_hflux_pac(lat,lon,qtoce_D028,grid_area,mask_ocean)
qtoce_pac_D029,notused = compute_hflux_pac(lat,lon,qtoce_D029,grid_area,mask_ocean)
qtoce_pac_D030,notused = compute_hflux_pac(lat,lon,qtoce_D030,grid_area,mask_ocean)
print('Total area Pacific (1e6 km2):',area_pac)
print('Total net surface heat flux Pacific CTRL (TW):',qtoce_pac_D000)
print('Total net surface heat flux Pacific ATL1+3K (TW):',qtoce_pac_D012,qtoce_pac_D012-qtoce_pac_D000)
print('Total net surface heat flux Pacific PAC1+3K (TW):',qtoce_pac_D021,qtoce_pac_D021-qtoce_pac_D000)

# Save net surface heat fluxes
if save_var == True:
    filename = dir_output + 'TotalSurfaceFluxes_atl.npy'
    np.save(filename,[qtoce_atl_D000,qtoce_atl_D012,qtoce_atl_D013,qtoce_atl_D014,qtoce_atl_D015,qtoce_atl_D016,qtoce_atl_D017,qtoce_atl_D018,qtoce_atl_D019,qtoce_atl_D020,qtoce_atl_D021,qtoce_atl_D022,qtoce_atl_D023,qtoce_atl_D024,qtoce_atl_D025,qtoce_atl_D027,qtoce_atl_D028,qtoce_atl_D029,qtoce_atl_D030])
    filename = dir_output + 'TotalSurfaceFluxes_pac.npy'
    np.save(filename,[qtoce_pac_D000,qtoce_pac_D012,qtoce_pac_D013,qtoce_pac_D014,qtoce_pac_D015,qtoce_pac_D016,qtoce_pac_D017,qtoce_pac_D018,qtoce_pac_D019,qtoce_pac_D020,qtoce_pac_D021,qtoce_pac_D022,qtoce_pac_D023,qtoce_pac_D024,qtoce_pac_D025,qtoce_pac_D027,qtoce_pac_D028,qtoce_pac_D029,qtoce_pac_D030])

# Load net surface heat flux at ice surface D000
filename = dir_input + '/D000/qtice_D000_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D000 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D012
filename = dir_input + '/D012/qtice_D012_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D012 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D015
filename = dir_input + '/D015/qtice_D015_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D015 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D018
filename = dir_input + '/D018/qtice_D018_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D018 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D021
filename = dir_input + '/D021/qtice_D021_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D021 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D022
filename = dir_input + '/D022/qtice_D022_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D022 = fh.variables['qt_ice'][0,:,:]
fh.close()

# Load net surface heat flux at ice surface D023
filename = dir_input + '/D023/qtice_D023_mean.nc'
fh = Dataset(filename, mode='r')
qtice_D023 = fh.variables['qt_ice'][0,:,:]
fh.close()

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


# Figure - Maps 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,qtoce_D000,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-100,-50,0,50,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Net ocean surface \n heat flux (W m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
cs=map.pcolor(x,y,qtoce_D012-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,qtoce_D015-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,qtoce_D018-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,qtoce_D021-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,qtoce_D022-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,qtoce_D023-qtoce_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar.set_label('Net OSHF PERT - CTRL (W m$^{-2}$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'Map_Ocean_NetSurfaceFlux.png')


# Figure - Maps 50-year average
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)

# D000
cs=map.pcolor(x,y,qtice_D000,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
map.drawcoastlines(ax=ax[0,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
ax[0,0].set_title('CTRL',fontsize=32)
ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
ax[0,0].yaxis.set_label_coords(-0.05,0.9)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-100,-50,0,50,100],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Net ice surface \n heat flux (W m$^{-2}$)',fontsize=28)

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
cs=map.pcolor(x,y,qtice_D012-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
map.drawcoastlines(ax=ax[1,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
ax[1,0].yaxis.set_label_coords(-0.05,0.9)

# D015
cs=map.pcolor(x,y,qtice_D015-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
map.drawcoastlines(ax=ax[1,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
ax[1,1].yaxis.set_label_coords(-0.05,0.9)

# D018
cs=map.pcolor(x,y,qtice_D018-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
map.drawcoastlines(ax=ax[1,2])
map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
ax[1,2].yaxis.set_label_coords(-0.05,0.9)

# D021
cs=map.pcolor(x,y,qtice_D021-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
map.drawcoastlines(ax=ax[2,0])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
ax[2,0].yaxis.set_label_coords(-0.05,0.9)

# D022
cs=map.pcolor(x,y,qtice_D022-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
map.drawcoastlines(ax=ax[2,1])
map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
ax[2,1].yaxis.set_label_coords(-0.05,0.9)

# D023
cs=map.pcolor(x,y,qtice_D023-qtice_D000,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
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
cbar.set_label('Net ISHF PERT - CTRL (W m$^{-2}$)',fontsize=28)

# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'Map_Ice_NetSurfaceFlux.png')
