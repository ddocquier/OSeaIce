#!/usr/bin/env python

'''
GOAL
    Decompose OHT into velocity and temperature anomalies
PROGRAMMER
    D. Docquier
LAST UPDATE
    23/04/2020
'''

# Options
exp = 'D012'
start_year = 2130
end_year = 2179
save_var = True

# Time parameters
start_folder = int(start_year - 2130 + 281)
nyears = int(end_year-start_year+1)
nmy = int(12) # number of months in a year

# Standard libraries
from netCDF4 import Dataset
import numpy as np
import time
start_time = time.time()

# Working directory
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/run/' + str(exp) + '/output/nemo/'
dir_ctrl = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/run/D000/output/nemo/'
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
rho = 1027. # seawater density (1000 kg m^{-3} [CDFTOOLS] or 1027 kg m^{-3} [Lien et al., 2017])
cp = 3985. # specific ocean heat capacity (4000  J kg^{-1} K^{-1} [CDFTOOLS] or 3985 J kg^{-1} K^{-1} [Lien et al., 2017])
tref = 0. # reference ocean temperature, set to 0 degC
threshold = 3000.

# Initialization of velocity and temperature anomalies in OHT and their covariance
u_ano = np.zeros((nyears,nmy,ny,nx))
v_ano = np.zeros((nyears,nmy,ny,nx))
t_u_ano = np.zeros((nyears,nmy,ny,nx))
t_v_ano = np.zeros((nyears,nmy,ny,nx))
cov_u_ano = np.zeros((nyears,nmy,ny,nx))
cov_v_ano = np.zeros((nyears,nmy,ny,nx))

# Loop over years
for year in np.arange(nyears):
    print(start_year+year)

    # Retrieve perturbed zonal velocity u (m/s)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_U.nc'
    fh = Dataset(filename, mode='r')
    u_s = fh.variables['uo'][:]
    fh.close()

    # Retrieve perturbed meridional velocity v (m/s)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_V.nc'
    fh = Dataset(filename, mode='r')
    v_s = fh.variables['vo'][:]
    fh.close()

    # Retrieve perturbed potential temperature (degC)
    filename = dir_input + str(start_folder+year) + '/' + str(exp) + '_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_T.nc'
    fh = Dataset(filename, mode='r')
    t_s = fh.variables['thetao'][:]
    fh.close()

    # Retrieve control zonal velocity u (m/s)
    filename = dir_ctrl + str(start_folder+year) + '/D000_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_U.nc'
    fh = Dataset(filename, mode='r')
    u_c = fh.variables['uo'][:]
    fh.close()

    # Retrieve control meridional velocity v (m/s)
    filename = dir_ctrl + str(start_folder+year) + '/D000_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_V.nc'
    fh = Dataset(filename, mode='r')
    v_c = fh.variables['vo'][:]
    fh.close()

    # Retrieve control potential temperature (degC)
    filename = dir_ctrl + str(start_folder+year) + '/D000_1m_'+str(start_year+year)+'0101_'+str(start_year+year)+'1231_grid_T.nc'
    fh = Dataset(filename, mode='r')
    t_c = fh.variables['thetao'][:]
    fh.close()
    
    # Compute anomalies
    t_s_u = np.zeros((nmy,nz,ny,nx))
    t_s_v = np.zeros((nmy,nz,ny,nx))
    t_c_u = np.zeros((nmy,nz,ny,nx))
    t_c_v = np.zeros((nmy,nz,ny,nx))
    diff_u = np.zeros((nmy,nz,ny,nx))
    diff_v = np.zeros((nmy,nz,ny,nx))
    diff_t_u = np.zeros((nmy,nz,ny,nx))
    diff_t_v = np.zeros((nmy,nz,ny,nx))
    diff_cov_u = np.zeros((nmy,nz,ny,nx))
    diff_cov_v = np.zeros((nmy,nz,ny,nx))
    for t in np.arange(nmy):
        print(t+1)
        for z in np.arange(nz):
            # Compute temperature on u and v grids
            t_s_u[t,z,:,:],t_s_v[t,z,:,:] = movingaverage(t_s[t,z,:,:],2)
            t_c_u[t,z,:,:],t_c_v[t,z,:,:] = movingaverage(t_c[t,z,:,:],2)

            # Masks to exclude crazy values
            mask_t_s_u = np.zeros((ny,nx))
            mask_t_s_v = np.zeros((ny,nx))
            mask_t_c_u = np.zeros((ny,nx))
            mask_t_c_v = np.zeros((ny,nx))
            mask_t_s_u[(t_s_u[t,z,:,:] >= -threshold) * (t_s_u[t,z,:,:] <= threshold)] = 1
            mask_t_s_v[(t_s_v[t,z,:,:] >= -threshold) * (t_s_v[t,z,:,:] <= threshold)] = 1
            mask_t_c_u[(t_c_u[t,z,:,:] >= -threshold) * (t_c_u[t,z,:,:] <= threshold)] = 1
            mask_t_c_v[(t_c_v[t,z,:,:] >= -threshold) * (t_c_v[t,z,:,:] <= threshold)] = 1

            # u anomalies    
            diff_u[t,z,:,:] = (u_s[t,z,:,:]-u_c[t,z,:,:]) * (t_c_u[t,z,:,:]-tref) * mask_t_c_u
            mask_diff_u = np.zeros((ny,nx))
            mask_diff_u[(diff_u[t,z,:,:] >= -threshold) * (diff_u[t,z,:,:] <= threshold)] = 1
            dwkh_u = diff_u[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_u
            u_ano[year,t,:,:] = u_ano[year,t,:,:] + dwkh_u * rho * cp

            # v anomalies    
            diff_v[t,z,:,:] = (v_s[t,z,:,:]-v_c[t,z,:,:]) * (t_c_v[t,z,:,:]-tref) * mask_t_c_v
            mask_diff_v = np.zeros((ny,nx))
            mask_diff_v[(diff_v[t,z,:,:] >= -threshold) * (diff_v[t,z,:,:] <= threshold)] = 1
            dwkh_v = diff_v[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_v
            v_ano[year,t,:,:] = v_ano[year,t,:,:] + dwkh_v * rho * cp

            # t_u anomalies    
            diff_t_u[t,z,:,:] = u_c[t,z,:,:] * (t_s_u[t,z,:,:]-t_c_u[t,z,:,:]) * mask_t_s_u * mask_t_c_u
            mask_diff_t_u = np.zeros((ny,nx))
            mask_diff_t_u[(diff_t_u[t,z,:,:] >= -threshold) * (diff_t_u[t,z,:,:] <= threshold)] = 1
            dwkh_t_u = diff_t_u[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_t_u
            t_u_ano[year,t,:,:] = t_u_ano[year,t,:,:] + dwkh_t_u * rho * cp

            # t_v anomalies    
            diff_t_v[t,z,:,:] = v_c[t,z,:,:] * (t_s_v[t,z,:,:]-t_c_v[t,z,:,:]) * mask_t_s_v * mask_t_c_v
            mask_diff_t_v = np.zeros((ny,nx))
            mask_diff_t_v[(diff_t_v[t,z,:,:] >= -threshold) * (diff_t_v[t,z,:,:] <= threshold)] = 1
            dwkh_t_v = diff_t_v[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_t_v
            t_v_ano[year,t,:,:] = t_v_ano[year,t,:,:] + dwkh_t_v * rho * cp

            # cov_u anomalies    
            diff_cov_u[t,z,:,:] = (u_s[t,z,:,:]-u_c[t,z,:,:]) * (t_s_u[t,z,:,:]-t_c_u[t,z,:,:]) * mask_t_s_u * mask_t_c_u
            mask_diff_cov_u = np.zeros((ny,nx))
            mask_diff_cov_u[(diff_cov_u[t,z,:,:] >= -threshold) * (diff_cov_u[t,z,:,:] <= threshold)] = 1
            dwkh_cov_u = diff_cov_u[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_cov_u
            cov_u_ano[year,t,:,:] = cov_u_ano[year,t,:,:] + dwkh_cov_u * rho * cp

            # cov_v anomalies    
            diff_cov_v[t,z,:,:] = (v_s[t,z,:,:]-v_c[t,z,:,:]) * (t_s_v[t,z,:,:]-t_c_v[t,z,:,:]) * mask_t_s_v * mask_t_c_v
            mask_diff_cov_v = np.zeros((ny,nx))
            mask_diff_cov_v[(diff_cov_v[t,z,:,:] >= -threshold) * (diff_cov_v[t,z,:,:] <= threshold)] = 1
            dwkh_cov_v = diff_cov_v[t,z,:,:] * gridy * gridz[z,:,:] * mask_diff_cov_v
            cov_v_ano[year,t,:,:] = cov_v_ano[year,t,:,:] + dwkh_cov_v * rho * cp

# Save variables
if save_var == True:
    filename = dir_output + 'decoht_' + str(exp) + '.npy'
    np.save(filename,[u_ano,v_ano,t_u_ano,t_v_ano,cov_u_ano,cov_v_ano])

# Computing time
print("--- %s seconds ---" % (time.time() - start_time))
