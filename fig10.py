#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 10: Plot ocean temperature Arctic (thetao), previously saved with save_thetao.py (Atlantic side) and save_thetao_arctic.py (Pacific side)
PROGRAMMER
    D. Docquier
LAST UPDATE
    17/11/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Option
save_fig = True

# Working directory
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/'
dir_input2 = dir_input + 'run/D000/output/nemo/281/'
dir_D000 = dir_input + 'post-proc/D000/'
dir_D012 = dir_input + 'post-proc/D012/'
dir_D015 = dir_input + 'post-proc/D015/'
dir_D018 = dir_input + 'post-proc/D018/'
dir_D021 = dir_input + 'post-proc/D021/'
dir_D022 = dir_input + 'post-proc/D022/'
dir_D023 = dir_input + 'post-proc/D023/'
dir_output = dir_input + 'OSeaIce_Paper/'

# Load coordinates
filename = dir_input2 + 'D000_1m_21300101_21301231_grid_T.nc'
fh = Dataset(filename, mode='r')
depth = fh.variables['olevel'][:]
lat = fh.variables['nav_lat'][:]
lat_atl = lat[:,265]
lat_pac = lat[:,95]
lon = fh.variables['nav_lon'][:]
fh.close()

# Load thetao Atlantic
thetao_atl_D000 = np.load(dir_D000 + 'thetao_D000.npy')
thetao_atl_D012 = np.load(dir_D012 + 'thetao_D012.npy')
thetao_atl_D015 = np.load(dir_D015 + 'thetao_D015.npy')
thetao_atl_D018 = np.load(dir_D018 + 'thetao_D018.npy')
thetao_atl_D021 = np.load(dir_D021 + 'thetao_D021.npy')
thetao_atl_D022 = np.load(dir_D022 + 'thetao_D022.npy')
thetao_atl_D023 = np.load(dir_D023 + 'thetao_D023.npy')

# Load thetao Pacific
thetao_pac_D000 = np.load(dir_D000 + 'thetao_arc_D000.npy')
thetao_pac_D012 = np.load(dir_D012 + 'thetao_arc_D012.npy')
thetao_pac_D015 = np.load(dir_D015 + 'thetao_arc_D015.npy')
thetao_pac_D018 = np.load(dir_D018 + 'thetao_arc_D018.npy')
thetao_pac_D021 = np.load(dir_D021 + 'thetao_arc_D021.npy')
thetao_pac_D022 = np.load(dir_D022 + 'thetao_arc_D022.npy')
thetao_pac_D023 = np.load(dir_D023 + 'thetao_arc_D023.npy')

# Palettes
palette_var = plt.cm.seismic
min_var = -5.
max_var = 5.
levels_var = np.linspace(min_var,max_var,20)
palette_diff = plt.cm.seismic
min_diff = -3.
max_diff = 3.
levels_diff = np.linspace(min_diff,max_diff,20)


# Fig. 10 - Plots of thetao difference
fig,axes=plt.subplots(6,3,figsize=(18,26))
fig.subplots_adjust(left=0.08,bottom=0.05,right=0.95,top=0.95,wspace=0.3,hspace=0.35)

# D000 Atlantic
axes[0,0].set_title('CTRL',fontsize=30)
axes[0,0].set_title('a',loc='left',fontsize=30,fontweight='bold')
cs = axes[0,0].contourf(lat_atl,depth,thetao_atl_D000,cmap=palette_var,levels=levels_var,extend='both')
axes[0,0].tick_params(axis='both',labelsize=20)
axes[0,0].axis([60, 90, 0, 1000])
axes[0,0].invert_yaxis()
axes[0,0].set_ylabel('Depth (m)',fontsize=25)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.82, 0.015, 0.13])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-5,-2.5,0,2.5,5],extend='both')
cbar.ax.tick_params(labelsize=20)
cbar.set_label('Ocean temperature ($^{\circ}$C)',fontsize=20)

# Remove axes that are not necessary
fig.delaxes(axes[0,1])
fig.delaxes(axes[0,2])

# D012 Atlantic
axes[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=30)
axes[1,0].set_title('b',loc='left',fontsize=30,fontweight='bold')
cs = axes[1,0].contourf(lat_atl,depth,thetao_atl_D012-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[1,0].tick_params(axis='both',labelsize=20)
axes[1,0].axis([60, 90, 0, 1000])
axes[1,0].invert_yaxis()
axes[1,0].set_ylabel('Depth (m)',fontsize=25)

# D015 Atlantic
axes[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=30)
axes[1,1].set_title('c',loc='left',fontsize=30,fontweight='bold')
cs = axes[1,1].contourf(lat_atl,depth,thetao_atl_D015-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[1,1].tick_params(axis='both',labelsize=20)
axes[1,1].axis([60, 90, 0, 1000])
axes[1,1].invert_yaxis()

# D018 Atlantic
axes[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=30)
axes[1,2].set_title('d',loc='left',fontsize=30,fontweight='bold')
cs = axes[1,2].contourf(lat_atl,depth,thetao_atl_D018-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[1,2].tick_params(axis='both',labelsize=20)
axes[1,2].axis([60, 90, 0, 1000])
axes[1,2].invert_yaxis()

# D021 Atlantic
axes[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=30)
axes[2,0].set_title('e',loc='left',fontsize=30,fontweight='bold')
cs = axes[2,0].contourf(lat_atl,depth,thetao_atl_D021-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[2,0].tick_params(axis='both',labelsize=20)
axes[2,0].axis([60, 90, 0, 1000])
axes[2,0].invert_yaxis()
axes[2,0].set_ylabel('Depth (m)',fontsize=25)

# D022 Atlantic
axes[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=30)
axes[2,1].set_title('f',loc='left',fontsize=30,fontweight='bold')
cs = axes[2,1].contourf(lat_atl,depth,thetao_atl_D022-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[2,1].set_xlabel('Latitude ($^{\circ}$)',fontsize=25)
axes[2,1].tick_params(axis='both',labelsize=20)
axes[2,1].axis([60, 90, 0, 1000])
axes[2,1].invert_yaxis()

# D023 Atlantic
axes[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=30)
axes[2,2].set_title('g',loc='left',fontsize=30,fontweight='bold')
cs = axes[2,2].contourf(lat_atl,depth,thetao_atl_D023-thetao_atl_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[2,2].set_xlabel('Latitude ($^{\circ}$)',fontsize=25)
axes[2,2].tick_params(axis='both',labelsize=20)
axes[2,2].axis([60, 90, 0, 1000])
axes[2,2].invert_yaxis()

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.88, 0.4, 0.01])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Temperature PERT-CTRL ($^{\circ}$C)',fontsize=24)


# D000 Pacific
axes[3,0].set_title('CTRL',fontsize=30)
axes[3,0].set_title('h',loc='left',fontsize=30,fontweight='bold')
cs = axes[3,0].contourf(lat_pac,depth,thetao_pac_D000,cmap=palette_var,levels=levels_var,extend='both')
axes[3,0].tick_params(axis='both',labelsize=20)
axes[3,0].axis([60, 90, 0, 1000])
axes[3,0].invert_yaxis()
axes[3,0].set_ylabel('Depth (m)',fontsize=25)

# Add color bar absolute value
cb_ax = fig.add_axes([0.35, 0.35, 0.015, 0.13])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-5,-2.5,0,2.5,5],extend='both')
cbar.ax.tick_params(labelsize=20)
cbar.set_label('Ocean temperature ($^{\circ}$C)',fontsize=20)

# Remove axes that are not necessary
fig.delaxes(axes[3,1])
fig.delaxes(axes[3,2])

# D012 Pacific
axes[4,0].set_title('ATL1+3$^{\circ}$C',fontsize=30)
axes[4,0].set_title('i',loc='left',fontsize=30,fontweight='bold')
cs = axes[4,0].contourf(lat_pac,depth,thetao_pac_D012-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[4,0].tick_params(axis='both',labelsize=20)
axes[4,0].axis([60, 90, 0, 1000])
axes[4,0].invert_yaxis()
axes[4,0].set_ylabel('Depth (m)',fontsize=25)

# D015 Pacific
axes[4,1].set_title('ATL2+3$^{\circ}$C',fontsize=30)
axes[4,1].set_title('j',loc='left',fontsize=30,fontweight='bold')
cs = axes[4,1].contourf(lat_pac,depth,thetao_pac_D015-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[4,1].tick_params(axis='both',labelsize=20)
axes[4,1].axis([60, 90, 0, 1000])
axes[4,1].invert_yaxis()

# D018 Pacific
axes[4,2].set_title('ATL3+3$^{\circ}$C',fontsize=30)
axes[4,2].set_title('k',loc='left',fontsize=30,fontweight='bold')
cs = axes[4,2].contourf(lat_pac,depth,thetao_pac_D018-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[4,2].tick_params(axis='both',labelsize=20)
axes[4,2].axis([60, 90, 0, 1000])
axes[4,2].invert_yaxis()

# D021 Pacific
axes[5,0].set_title('PAC1+3$^{\circ}$C',fontsize=30)
axes[5,0].set_title('l',loc='left',fontsize=30,fontweight='bold')
cs = axes[5,0].contourf(lat_pac,depth,thetao_pac_D021-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[5,0].set_xlabel('Latitude ($^{\circ}$)',fontsize=25)
axes[5,0].tick_params(axis='both',labelsize=20)
axes[5,0].axis([60, 90, 0, 1000])
axes[5,0].invert_yaxis()
axes[5,0].set_ylabel('Depth (m)',fontsize=25)

# D022 Pacific
axes[5,1].set_title('PAC2+3$^{\circ}$C',fontsize=30)
axes[5,1].set_title('m',loc='left',fontsize=30,fontweight='bold')
cs = axes[5,1].contourf(lat_pac,depth,thetao_pac_D022-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[5,1].set_xlabel('Latitude ($^{\circ}$)',fontsize=25)
axes[5,1].tick_params(axis='both',labelsize=20)
axes[5,1].axis([60, 90, 0, 1000])
axes[5,1].invert_yaxis()

# D023 Pacific
axes[5,2].set_title('PAC3+3$^{\circ}$C',fontsize=30)
axes[5,2].set_title('n',loc='left',fontsize=30,fontweight='bold')
cs = axes[5,2].contourf(lat_pac,depth,thetao_pac_D023-thetao_pac_D000,cmap=palette_diff,levels=levels_diff,extend='both')
axes[5,2].set_xlabel('Latitude ($^{\circ}$)',fontsize=25)
axes[5,2].tick_params(axis='both',labelsize=20)
axes[5,2].axis([60, 90, 0, 1000])
axes[5,2].invert_yaxis()

# Add color bar diff
cb_ax = fig.add_axes([0.5, 0.41, 0.4, 0.01])
cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-3,-2,-1,0,1,2,3],extend='both')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Temperature PERT-CTRL ($^{\circ}$C)',fontsize=24)
      
# Save figure
if save_fig == True:
    fig.savefig(dir_output + 'fig10.png')
    fig.savefig(dir_output + 'fig10.eps',dpi=300)
