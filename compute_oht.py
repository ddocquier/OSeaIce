#!/usr/bin/env python

'''
GOAL
    Compute zonal and meridional OHT at every grid point (integrated over the vertical)
PROGRAMMER
    D. Docquier
LAST UPDATE
    24/04/2020
'''

# Options
exp = 'D012'
save_var = True
start_year = 2130
end_year = 2179

# Time parameters
start_folder = int(start_year - 2130 + 281)
nyears = int(end_year-start_year+1)
nmy = int(12) # number of months in a year

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import time
start_time = time.time()

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/run/' + str(exp) + '/output/nemo/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/grid/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/' + str(exp) + '/OHT_transects/'

# Function to compute moving average and shift to the left (to get T on u-grid)
def movingaverage(variable,window_size):
    window = np.ones(int(window_size))/float(window_size)
    runningmean_u = np.apply_along_axis(lambda m: np.convolve(m,window,'same'),axis=1,arr=variable)
    runningmean_u[:,:-1] = runningmean_u[:,1:]
    runningmean_v = np.apply_along_axis(lambda m: np.convolve(m,window,'same'),axis=0,arr=variable)
    runningmean_v[:-1,:] = runningmean_v[1:,:]
    return runningmean_u,runningmean_v

# Load grid sizes in X (gridx), Y (gridy)
filename = dir_grid + 'coordinates.nc'
fh = Dataset(filename, mode='r')
gridx = fh.variables['e1u'][:]
gridy = fh.variables['e2u'][:]
fh.close()

# Load grid size in Z (gridz)
filename = dir_grid + 'mesh_zgr.nc'
fh = Dataset(filename, mode='r')
gridz = fh.variables['e3u'][:]
gridz = gridz[0,:,:,:]
nz,ny,nx = gridz.shape
fh.close()

# Constant parameters
rho = 1000. # seawater density (1000 kg m^{-3} [CDFTOOLS] or 1027 kg m^{-3} [Lien et al., 2017])
cp = 4000. # specific ocean heat capacity (4000  J kg^{-1} K^{-1} [CDFTOOLS] or 3985 J kg^{-1} K^{-1} [Lien et al., 2017])
tref = 0. # reference ocean temperature, set to 0 degC

# Initialization of zonal and meridional OHT
oht_u = np.zeros((nyears,nmy,ny,nx))
oht_v = np.zeros((nyears,nmy,ny,nx))

# Loop over years
for year in np.arange(nyears):
    print(start_year+year)

    # Retrieve zonal velocity u (m/s)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_U.nc'
    fh = Dataset(filename, mode='r')
    u = fh.variables['uo'][:]
    fh.close()
    
    # Retrieve meridional velocity v (m/s)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_V.nc'
    fh = Dataset(filename, mode='r')
    v = fh.variables['vo'][:]
    fh.close()
    
    # Retrieve potential temperature (degC)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_T.nc'
    fh = Dataset(filename, mode='r')
    temp = fh.variables['thetao'][:]
    fh.close()

    # Compute OHT
    temp_u = np.zeros((nmy,nz,ny,nx))
    ut = np.zeros((nmy,nz,ny,nx))
    temp_v = np.zeros((nmy,nz,ny,nx))
    vt = np.zeros((nmy,nz,ny,nx))
    for t in np.arange(nmy):
        print(t+1)
        for z in np.arange(nz):
            temp_u[t,z,:,:],temp_v[t,z,:,:] = movingaverage(temp[t,z,:,:],2)
            mask_temp_u = np.zeros((ny,nx))
            mask_temp_u[(temp_u[t,z,:,:] >= -3000.) * (temp_u[t,z,:,:] <= 3000.)] = 1
            ut[t,z,:,:] = u[t,z,:,:] * (temp_u[t,z,:,:]-tref) * mask_temp_u
            mask_ut = np.zeros((ny,nx))
            mask_ut[(ut[t,z,:,:] >= -3000.) * (ut[t,z,:,:] <= 3000.)] = 1
            dwkh_u = ut[t,z,:,:] * gridy * gridz[z,:,:] * mask_ut
            oht_u[year,t,:,:] = oht_u[year,t,:,:] + dwkh_u * rho * cp
            mask_temp_v = np.zeros((ny,nx))
            mask_temp_v[(temp_v[t,z,:,:] >= -3000.) * (temp_v[t,z,:,:] <= 3000.)] = 1
            vt[t,z,:,:] = v[t,z,:,:] * (temp_v[t,z,:,:]-tref) * mask_temp_v
            mask_vt = np.zeros((ny,nx))
            mask_vt[(vt[t,z,:,:] >= -3000.) * (vt[t,z,:,:] <= 3000.)] = 1
            dwkh_v = vt[t,z,:,:] * gridx * gridz[z,:,:] * mask_vt
            oht_v[year,t,:,:] = oht_v[year,t,:,:] + dwkh_v * rho * cp

# Save variables
if save_var == True:
    filename = dir_output + 'oht_' + str(exp) + '.npy'
    np.save(filename,[oht_u,oht_v])
    
# Computing time
print("--- %s seconds ---" % (time.time() - start_time))
