#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Figs. 12-13: Maps of melt/growth terms (m/year)
    Compute mean quantities for Fig. 11
PROGRAMMER
    D. Docquier
LAST UPDATE
    25/05/2020
'''

# Chosen variable
varname = 'vfxice'

# Parameters
siconc_threshold = 0.
lat_threshold = 80.
plot_fig = False
save_fig = False
save_var = True

# Provide name for color bar based on variable name
if varname == 'vfxice':
    cbarname = 'Net ice growth'
elif varname == 'vfxbog':
    cbarname = 'Basal growth'
elif varname == 'vfxopw':
    cbarname = 'Open-water growth'
elif varname == 'vfxbom':
    cbarname = 'Basal melt'
elif varname == 'vfxsum':
    cbarname = 'Surface melt'
elif varname == 'vfxdyn':
    cbarname = 'Dynamic growth'
elif varname == 'vfxsni':
    cbarname = 'Snow-ice formation'

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/'
dir_D000 = dir_input + 'post-proc/D001/'
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
dir_grid = dir_input + 'grid/'
dir_output = dir_input + 'OSeaIce_Paper/'

# Load var D000
filename = dir_D000 + varname + '_D001_ym.nc'
fh = Dataset(filename, mode='r')
var_D000 = fh.variables[varname][:]
var_D000 = var_D000 * 365.25
lat = fh.variables['nav_lat_grid_T'][:]
lon = fh.variables['nav_lon_grid_T'][:]
nm,ny,nx = var_D000.shape
fh.close()

# Load var D012
filename = dir_D012 + varname + '_D012_ym.nc'
fh = Dataset(filename, mode='r')
var_D012 = fh.variables[varname][:]
var_D012 = var_D012 * 365.25
fh.close()

# Load var D013
filename = dir_D013 + varname + '_D013_ym.nc'
fh = Dataset(filename, mode='r')
var_D013 = fh.variables[varname][:]
var_D013 = var_D013 * 365.25
fh.close()

# Load var D014
filename = dir_D014 + varname + '_D014_ym.nc'
fh = Dataset(filename, mode='r')
var_D014 = fh.variables[varname][:]
var_D014 = var_D014 * 365.25
fh.close()

# Load var D015
filename = dir_D015 + varname + '_D015_ym.nc'
fh = Dataset(filename, mode='r')
var_D015 = fh.variables[varname][:]
var_D015 = var_D015 * 365.25
fh.close()

# Load var D016
filename = dir_D016 + varname + '_D016_ym.nc'
fh = Dataset(filename, mode='r')
var_D016 = fh.variables[varname][:]
var_D016 = var_D016 * 365.25
fh.close()

# Load var D017
filename = dir_D017 + varname + '_D017_ym.nc'
fh = Dataset(filename, mode='r')
var_D017 = fh.variables[varname][:]
var_D017 = var_D017 * 365.25
fh.close()

# Load var D018
filename = dir_D018 + varname + '_D018_ym.nc'
fh = Dataset(filename, mode='r')
var_D018 = fh.variables[varname][:]
var_D018 = var_D018 * 365.25
fh.close()

# Load var D019
filename = dir_D019 + varname + '_D019_ym.nc'
fh = Dataset(filename, mode='r')
var_D019 = fh.variables[varname][:]
var_D019 = var_D019 * 365.25
fh.close()

# Load var D020
filename = dir_D020 + varname + '_D020_ym.nc'
fh = Dataset(filename, mode='r')
var_D020 = fh.variables[varname][:]
var_D020 = var_D020 * 365.25
fh.close()

# Load var D021
filename = dir_D021 + varname + '_D021_ym.nc'
fh = Dataset(filename, mode='r')
var_D021 = fh.variables[varname][:]
var_D021 = var_D021 * 365.25
fh.close()

# Load var D022
filename = dir_D022 + varname + '_D022_ym.nc'
fh = Dataset(filename, mode='r')
var_D022 = fh.variables[varname][:]
var_D022 = var_D022 * 365.25
fh.close()

# Load var D023
filename = dir_D023 + varname + '_D023_ym.nc'
fh = Dataset(filename, mode='r')
var_D023 = fh.variables[varname][:]
var_D023 = var_D023 * 365.25
fh.close()

# Load var D024
filename = dir_D024 + varname + '_D024_ym.nc'
fh = Dataset(filename, mode='r')
var_D024 = fh.variables[varname][:]
var_D024 = var_D024 * 365.25
fh.close()

# Load var D025
filename = dir_D025 + varname + '_D025_ym.nc'
fh = Dataset(filename, mode='r')
var_D025 = fh.variables[varname][:]
var_D025 = var_D025 * 365.25
fh.close()

# Load var D027
filename = dir_D027 + varname + '_D027_ym.nc'
fh = Dataset(filename, mode='r')
var_D027 = fh.variables[varname][:]
var_D027 = var_D027 * 365.25
fh.close()

# Load var D028
filename = dir_D028 + varname + '_D028_ym.nc'
fh = Dataset(filename, mode='r')
var_D028 = fh.variables[varname][:]
var_D028 = var_D028 * 365.25
fh.close()

# Load var D029
filename = dir_D029 + varname + '_D029_ym.nc'
fh = Dataset(filename, mode='r')
var_D029 = fh.variables[varname][:]
var_D029 = var_D029 * 365.25
fh.close()

# Load var D030
filename = dir_D030 + varname + '_D030_ym.nc'
fh = Dataset(filename, mode='r')
var_D030 = fh.variables[varname][:]
var_D030 = var_D030 * 365.25
fh.close()

# Compute mean var (over all months)
var_mean_D000 = np.nanmean(var_D000,axis=0)
var_mean_D012 = np.nanmean(var_D012,axis=0)
var_mean_D013 = np.nanmean(var_D013,axis=0)
var_mean_D014 = np.nanmean(var_D014,axis=0)
var_mean_D015 = np.nanmean(var_D015,axis=0)
var_mean_D016 = np.nanmean(var_D016,axis=0)
var_mean_D017 = np.nanmean(var_D017,axis=0)
var_mean_D018 = np.nanmean(var_D018,axis=0)
var_mean_D019 = np.nanmean(var_D019,axis=0)
var_mean_D020 = np.nanmean(var_D020,axis=0)
var_mean_D021 = np.nanmean(var_D021,axis=0)
var_mean_D022 = np.nanmean(var_D022,axis=0)
var_mean_D023 = np.nanmean(var_D023,axis=0)
var_mean_D024 = np.nanmean(var_D024,axis=0)
var_mean_D025 = np.nanmean(var_D025,axis=0)
var_mean_D027 = np.nanmean(var_D027,axis=0)
var_mean_D028 = np.nanmean(var_D028,axis=0)
var_mean_D029 = np.nanmean(var_D029,axis=0)
var_mean_D030 = np.nanmean(var_D030,axis=0)

# Convert ice production into positive numbers
if varname == 'vfxbog' or varname == 'vfxopw' or varname == 'vfxdyn' or varname == 'vfxsni' or varname == 'vfxice':
    var_mean_D000 = - var_mean_D000
    var_mean_D012 = - var_mean_D012
    var_mean_D013 = - var_mean_D013
    var_mean_D014 = - var_mean_D014
    var_mean_D015 = - var_mean_D015
    var_mean_D016 = - var_mean_D016
    var_mean_D017 = - var_mean_D017
    var_mean_D018 = - var_mean_D018
    var_mean_D019 = - var_mean_D019
    var_mean_D020 = - var_mean_D020
    var_mean_D021 = - var_mean_D021
    var_mean_D022 = - var_mean_D022
    var_mean_D023 = - var_mean_D023
    var_mean_D024 = - var_mean_D024
    var_mean_D025 = - var_mean_D025
    var_mean_D027 = - var_mean_D027
    var_mean_D028 = - var_mean_D028
    var_mean_D029 = - var_mean_D029
    var_mean_D030 = - var_mean_D030

# Load siconc D000
filename = dir_D000 + 'siconc_D001_ym.nc'
fh = Dataset(filename, mode='r')
siconc_D000 = fh.variables['siconc'][:]
siconc_D000 = siconc_D000 * 100.
grid_area = fh.variables['cell_area'][:]
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

# Compute mean siconc (over all months)
siconc_mean_D000 = np.nanmean(siconc_D000,axis=0)
siconc_mean_D012 = np.nanmean(siconc_D012,axis=0)
siconc_mean_D013 = np.nanmean(siconc_D013,axis=0)
siconc_mean_D014 = np.nanmean(siconc_D014,axis=0)
siconc_mean_D015 = np.nanmean(siconc_D015,axis=0)
siconc_mean_D016 = np.nanmean(siconc_D016,axis=0)
siconc_mean_D017 = np.nanmean(siconc_D017,axis=0)
siconc_mean_D018 = np.nanmean(siconc_D018,axis=0)
siconc_mean_D019 = np.nanmean(siconc_D019,axis=0)
siconc_mean_D020 = np.nanmean(siconc_D020,axis=0)
siconc_mean_D021 = np.nanmean(siconc_D021,axis=0)
siconc_mean_D022 = np.nanmean(siconc_D022,axis=0)
siconc_mean_D023 = np.nanmean(siconc_D023,axis=0)
siconc_mean_D024 = np.nanmean(siconc_D024,axis=0)
siconc_mean_D025 = np.nanmean(siconc_D025,axis=0)
siconc_mean_D027 = np.nanmean(siconc_D027,axis=0)
siconc_mean_D028 = np.nanmean(siconc_D028,axis=0)
siconc_mean_D029 = np.nanmean(siconc_D029,axis=0)
siconc_mean_D030 = np.nanmean(siconc_D030,axis=0)

# Compute difference between experiments and CTRL - over all months
var_diff_D012 = var_mean_D012 - var_mean_D000
var_diff_D012[siconc_mean_D012==0.] = 0.
var_diff_D015 = var_mean_D015 - var_mean_D000
var_diff_D015[siconc_mean_D015==0.] = 0.
var_diff_D018 = var_mean_D018 - var_mean_D000
var_diff_D018[siconc_mean_D018==0.] = 0.
var_diff_D021 = var_mean_D021 - var_mean_D000
var_diff_D021[siconc_mean_D021==0.] = 0.
var_diff_D022 = var_mean_D022 - var_mean_D000
var_diff_D022[siconc_mean_D022==0.] = 0.
var_diff_D023 = var_mean_D023 - var_mean_D000
var_diff_D023[siconc_mean_D023==0.] = 0.

# Compute spatial average
var_spatialmean_D000 = np.average(var_mean_D000[(siconc_mean_D000>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D000>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D012 = np.average(var_mean_D012[(siconc_mean_D012>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D012>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D013 = np.average(var_mean_D013[(siconc_mean_D013>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D013>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D014 = np.average(var_mean_D014[(siconc_mean_D014>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D014>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D015 = np.average(var_mean_D015[(siconc_mean_D015>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D015>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D016 = np.average(var_mean_D016[(siconc_mean_D016>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D016>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D017 = np.average(var_mean_D017[(siconc_mean_D017>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D017>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D018 = np.average(var_mean_D018[(siconc_mean_D018>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D018>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D019 = np.average(var_mean_D019[(siconc_mean_D019>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D019>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D020 = np.average(var_mean_D020[(siconc_mean_D020>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D020>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D021 = np.average(var_mean_D021[(siconc_mean_D021>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D021>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D022 = np.average(var_mean_D022[(siconc_mean_D022>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D022>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D023 = np.average(var_mean_D023[(siconc_mean_D023>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D023>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D024 = np.average(var_mean_D024[(siconc_mean_D024>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D024>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D025 = np.average(var_mean_D025[(siconc_mean_D025>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D025>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D027 = np.average(var_mean_D027[(siconc_mean_D027>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D027>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D028 = np.average(var_mean_D028[(siconc_mean_D028>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D028>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D029 = np.average(var_mean_D029[(siconc_mean_D029>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D029>siconc_threshold)*(lat>lat_threshold)])
var_spatialmean_D030 = np.average(var_mean_D030[(siconc_mean_D030>siconc_threshold)*(lat>lat_threshold)],weights=grid_area[(siconc_mean_D030>siconc_threshold)*(lat>lat_threshold)])

# Convert ice loss terms into negative numbers
if varname == 'vfxbom' or varname == 'vfxsum':
    var_spatialmean_D000 = - var_spatialmean_D000
    var_spatialmean_D012 = - var_spatialmean_D012
    var_spatialmean_D013 = - var_spatialmean_D013
    var_spatialmean_D014 = - var_spatialmean_D014
    var_spatialmean_D015 = - var_spatialmean_D015
    var_spatialmean_D016 = - var_spatialmean_D016
    var_spatialmean_D017 = - var_spatialmean_D017
    var_spatialmean_D018 = - var_spatialmean_D018
    var_spatialmean_D019 = - var_spatialmean_D019
    var_spatialmean_D020 = - var_spatialmean_D020
    var_spatialmean_D021 = - var_spatialmean_D021
    var_spatialmean_D022 = - var_spatialmean_D022
    var_spatialmean_D023 = - var_spatialmean_D023
    var_spatialmean_D024 = - var_spatialmean_D024
    var_spatialmean_D025 = - var_spatialmean_D025
    var_spatialmean_D027 = - var_spatialmean_D027
    var_spatialmean_D028 = - var_spatialmean_D028
    var_spatialmean_D029 = - var_spatialmean_D029
    var_spatialmean_D030 = - var_spatialmean_D030
           
# Save spatial averages
if save_var == True:
    if siconc_threshold == 0.:
        if lat_threshold == 80.:
            filename = dir_output + varname + '_siconc0_lat80.npy'
        else:
            filename = dir_output + varname + '_siconc0.npy'
    elif siconc_threshold == 15.:
        filename = dir_output + varname + '_siconc15.npy'
    elif siconc_threshold == 30.:
        filename = dir_output + varname + '_siconc30.npy'
    np.save(filename,[var_spatialmean_D000,var_spatialmean_D012,var_spatialmean_D013,var_spatialmean_D014,var_spatialmean_D015,var_spatialmean_D016,var_spatialmean_D017,var_spatialmean_D018,var_spatialmean_D019,var_spatialmean_D020,var_spatialmean_D021,var_spatialmean_D022,var_spatialmean_D023,var_spatialmean_D024,var_spatialmean_D025,var_spatialmean_D027,var_spatialmean_D028,var_spatialmean_D029,var_spatialmean_D030])
             
# Map projection
boundlat = 50.
l0 = 0.
map = Basemap(projection='nplaea',boundinglat=boundlat,lon_0=l0,resolution='c')
x,y = map(lon,lat)

# Palettes and plot parameters
if varname == 'vfxbom' or varname == 'vfxsum':
    palette_var = plt.cm.seismic._resample(40)
else:
    palette_var = plt.cm.seismic_r._resample(40)
min_var = -4.
max_var = 4.
palette_diff = palette_var
min_diff = -2.
max_diff = 2.


# Maps of difference in mean var between the SST restoring experiments and CTRL
if plot_fig == True:
    fig,ax=plt.subplots(3,3,figsize=(18,18))
    fig.subplots_adjust(left=0.06,bottom=0.05,right=0.95,top=0.95,wspace=0.2,hspace=0.3)
    
    # D000
    cs=map.pcolor(x,y,var_mean_D000,vmin=min_var,vmax=max_var,cmap=palette_var,ax=ax[0,0])
    map.contour(x,y,siconc_mean_D000,range(15,16,5),colors='g',ax=ax[0,0],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[0,0])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[0,0])
    map.drawcoastlines(ax=ax[0,0])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[0,0])
    ax[0,0].set_title('CTRL',fontsize=32)
    ax[0,0].set_title('a',loc='left',fontsize=32,fontweight='bold')
    ax[0,0].yaxis.set_label_coords(-0.05,0.9)
    
    # Delete axes
    fig.delaxes(ax[0,1])
    fig.delaxes(ax[0,2])
    
    # Add color bar absolute value
    cb_ax = fig.add_axes([0.35, 0.7, 0.015, 0.25])
    cbar = fig.colorbar(cs,cax=cb_ax,orientation='vertical',ticks=[-4,-2,0,2,4],extend='both')
    cbar.ax.tick_params(labelsize=24)
    cbar.set_label(cbarname + ' (m year$^{-1}$)',fontsize=28)
    
    # D012
    cs = map.pcolor(x,y,var_diff_D012,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,0])
    map.contour(x,y,siconc_mean_D012,range(15,16,5),colors='g',ax=ax[1,0],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,0])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,0])
    map.drawcoastlines(ax=ax[1,0])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[1,0])
    ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=32)
    ax[1,0].set_title('b',loc='left',fontsize=32,fontweight='bold')
    ax[1,0].yaxis.set_label_coords(-0.05,0.9)
    
    # D015
    cs = map.pcolor(x,y,var_diff_D015,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,1])
    map.contour(x,y,siconc_mean_D015,range(15,16,5),colors='g',ax=ax[1,1],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,1])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,1])
    map.drawcoastlines(ax=ax[1,1])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[1,1])
    ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=32)
    ax[1,1].set_title('c',loc='left',fontsize=32,fontweight='bold')
    ax[1,1].yaxis.set_label_coords(-0.05,0.9)
    
    # D018
    cs=map.pcolor(x,y,var_diff_D018,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[1,2])
    map.contour(x,y,siconc_mean_D018,range(15,16,5),colors='g',ax=ax[1,2],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[1,2])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[1,2])
    map.drawcoastlines(ax=ax[1,2])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[1,2])
    ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=32)
    ax[1,2].set_title('d',loc='left',fontsize=32,fontweight='bold')
    ax[1,2].yaxis.set_label_coords(-0.05,0.9)
    
    # D021
    cs = map.pcolor(x,y,var_diff_D021,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,0])
    map.contour(x,y,siconc_mean_D021,range(15,16,5),colors='g',ax=ax[2,0],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,0])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,0])
    map.drawcoastlines(ax=ax[2,0])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[2,0])
    ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=32)
    ax[2,0].set_title('e',loc='left',fontsize=32,fontweight='bold')
    ax[2,0].yaxis.set_label_coords(-0.05,0.9)
    
    # D022
    cs = map.pcolor(x,y,var_diff_D022,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,1])
    map.contour(x,y,siconc_mean_D022,range(15,16,5),colors='g',ax=ax[2,1],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,1])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,1])
    map.drawcoastlines(ax=ax[2,1])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[2,1])
    ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=32)
    ax[2,1].set_title('f',loc='left',fontsize=32,fontweight='bold')
    ax[2,1].yaxis.set_label_coords(-0.05,0.9)
    
    # D023
    cs = map.pcolor(x,y,var_diff_D023,vmin=min_diff,vmax=max_diff,cmap=palette_diff,ax=ax[2,2])
    map.contour(x,y,siconc_mean_D023,range(15,16,5),colors='g',ax=ax[2,2],linewidths=3)
    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,0],fontsize=24,ax=ax[2,2])
    map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],fontsize=24,ax=ax[2,2])
    map.drawcoastlines(ax=ax[2,2])
    map.fillcontinents(color='grey',lake_color='w',ax=ax[2,2])
    ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=32)
    ax[2,2].set_title('g',loc='left',fontsize=32,fontweight='bold')
    ax[2,2].yaxis.set_label_coords(-0.05,0.9)
    
    # Add color bar diff
    cb_ax = fig.add_axes([0.5, 0.82, 0.4, 0.02])
    cbar = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=[-2,-1,0,1,2],extend='both')
    cbar.ax.tick_params(labelsize=24)
    cbar.set_label(cbarname + ' PERT - CTRL (m year$^{-1}$)',fontsize=28)   
    
    # Save figure
    if save_fig == True:
        fig.savefig(dir_output+varname+'.png')
