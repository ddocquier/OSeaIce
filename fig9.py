'''
GOAL
    Fig. 9: OHT-sea ice scatter plots
    Supp.: tas-sea ice scatter plots
    Sea-ice files saved with fig5-6.py
    OHT files saved with fig2.py
    tas files saved with compute_tas.py 
PROGRAMMER
    D. Docquier
LAST UPDATE
    19/05/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats.stats import pearsonr
from netCDF4 import Dataset

# Option
save_fig = False

# Working directories
dir_fig = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'
dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/D000/IFS/'

# Function to compute slope + SD of slope
def compute_slope(nm,x,y):
    idx = np.isfinite(x) & np.isfinite(y)
    a,b = np.polyfit(x[idx],y[idx],1)
    lreg = a*x[idx]+b
    n_freedom = nm-2
    s_yx = np.sum((y[idx]-lreg)**2)/n_freedom
    SS_xx = np.sum((x[idx]-np.mean(x[idx]))**2)
    sd_a = np.sqrt(s_yx/SS_xx)
    t_a = a / sd_a
    alpha = 0.05
    t_crit = stats.t.ppf(1-alpha/2,n_freedom)
    sig_a = 0
    if np.abs(t_a) > t_crit:
        sig_a = 1
    return a,b,sd_a,sig_a

# Load grid-cell area, lat and lon
filename = dir_grid + 'gridarea_ifs.nc'
fh = Dataset(filename, mode='r')
gridarea = fh.variables['cell_area'][:]
lat_init = fh.variables['lat'][:]
lon_init = fh.variables['lon'][:]
lon,lat = np.meshgrid(lon_init,lat_init)
lon[lon>180.] = lon[lon>180.] - 360.
fh.close()

# Load annual mean total OHT
oht_annmean_total_D000,oht_annmean_total_D012,oht_annmean_total_D013,oht_annmean_total_D014,\
    oht_annmean_total_D015,oht_annmean_total_D016,oht_annmean_total_D017,oht_annmean_total_D018,\
    oht_annmean_total_D019,oht_annmean_total_D020,oht_annmean_total_D021,oht_annmean_total_D022,\
    oht_annmean_total_D023,oht_annmean_total_D024,oht_annmean_total_D025,oht_annmean_total_D027,\
    oht_annmean_total_D028,oht_annmean_total_D029,oht_annmean_total_D030 = np.load(dir_fig\
    + 'OHT_annmean.npy')

# Load annual mean tas
tas_annmean_D000,tas_annmean_D012,tas_annmean_D013,tas_annmean_D014,\
    tas_annmean_D015,tas_annmean_D016,tas_annmean_D017,tas_annmean_D018,\
    tas_annmean_D019,tas_annmean_D020,tas_annmean_D021,tas_annmean_D022,\
    tas_annmean_D023,tas_annmean_D024,tas_annmean_D025,tas_annmean_D027,\
    tas_annmean_D028,tas_annmean_D029,tas_annmean_D030 = np.load(dir_fig + 'tas_annmean.npy')

# Load March Arctic SIA
area_march_D000,area_march_D012,area_march_D013,area_march_D014,area_march_D015,area_march_D016,\
    area_march_D017,area_march_D018,area_march_D019,area_march_D020,area_march_D021,area_march_D022,\
    area_march_D023,area_march_D024,area_march_D025,area_march_D027,area_march_D028,area_march_D029,\
    area_march_D030 = np.load(dir_fig + 'SIA_March.npy',allow_pickle=True)

# Load September Arctic SIA
area_sept_D000,area_sept_D012,area_sept_D013,area_sept_D014,area_sept_D015,area_sept_D016,\
    area_sept_D017,area_sept_D018,area_sept_D019,area_sept_D020,area_sept_D021,area_sept_D022,\
    area_sept_D023,area_sept_D024,area_sept_D025,area_sept_D027,area_sept_D028,area_sept_D029,\
    area_sept_D030 = np.load(dir_fig + 'SIA_September.npy',allow_pickle=True)

# Load March Arctic SIV
volume_march_D000,volume_march_D012,volume_march_D013,volume_march_D014,volume_march_D015,volume_march_D016,\
    volume_march_D017,volume_march_D018,volume_march_D019,volume_march_D020,volume_march_D021,volume_march_D022,\
    volume_march_D023,volume_march_D024,volume_march_D025,volume_march_D027,volume_march_D028,volume_march_D029,\
    volume_march_D030 = np.load(dir_fig + 'SIV_March.npy',allow_pickle=True)

# Load September Arctic SIV
volume_sept_D000,volume_sept_D012,volume_sept_D013,volume_sept_D014,volume_sept_D015,volume_sept_D016,\
    volume_sept_D017,volume_sept_D018,volume_sept_D019,volume_sept_D020,volume_sept_D021,volume_sept_D022,\
    volume_sept_D023,volume_sept_D024,volume_sept_D025,volume_sept_D027,volume_sept_D028,volume_sept_D029,\
    volume_sept_D030 = np.load(dir_fig + 'SIV_September.npy',allow_pickle=True)

# Compute mean total OHT
oht_mean_total_D000 = np.nanmean(oht_annmean_total_D000)
print(oht_mean_total_D000)
oht_mean_total_D012 = np.nanmean(oht_annmean_total_D012)
oht_mean_total_D013 = np.nanmean(oht_annmean_total_D013)
oht_mean_total_D014 = np.nanmean(oht_annmean_total_D014)
oht_mean_total_D015 = np.nanmean(oht_annmean_total_D015)
oht_mean_total_D016 = np.nanmean(oht_annmean_total_D016)
oht_mean_total_D017 = np.nanmean(oht_annmean_total_D017)
oht_mean_total_D018 = np.nanmean(oht_annmean_total_D018)
oht_mean_total_D019 = np.nanmean(oht_annmean_total_D019)
oht_mean_total_D020 = np.nanmean(oht_annmean_total_D020)
oht_mean_total_D021 = np.nanmean(oht_annmean_total_D021)
oht_mean_total_D022 = np.nanmean(oht_annmean_total_D022)
oht_mean_total_D023 = np.nanmean(oht_annmean_total_D023)
oht_mean_total_D024 = np.nanmean(oht_annmean_total_D024)
oht_mean_total_D025 = np.nanmean(oht_annmean_total_D025)
oht_mean_total_D027 = np.nanmean(oht_annmean_total_D027)
oht_mean_total_D028 = np.nanmean(oht_annmean_total_D028)
oht_mean_total_D029 = np.nanmean(oht_annmean_total_D029)
oht_mean_total_D030 = np.nanmean(oht_annmean_total_D030)
oht_mean_total_array = np.array([oht_mean_total_D012-oht_mean_total_D000,oht_mean_total_D013-oht_mean_total_D000,\
                        oht_mean_total_D014-oht_mean_total_D000,oht_mean_total_D015-oht_mean_total_D000,\
                        oht_mean_total_D016-oht_mean_total_D000,oht_mean_total_D017-oht_mean_total_D000,\
                        oht_mean_total_D018-oht_mean_total_D000,oht_mean_total_D019-oht_mean_total_D000,\
                        oht_mean_total_D020-oht_mean_total_D000,oht_mean_total_D021-oht_mean_total_D000,\
                        oht_mean_total_D022-oht_mean_total_D000,oht_mean_total_D023-oht_mean_total_D000,\
                        oht_mean_total_D024-oht_mean_total_D000,oht_mean_total_D025-oht_mean_total_D000,\
                        oht_mean_total_D027-oht_mean_total_D000,oht_mean_total_D028-oht_mean_total_D000,\
                        oht_mean_total_D029-oht_mean_total_D000,oht_mean_total_D030-oht_mean_total_D000])

# Compute mean tas
tas_mean_D000 = np.nanmean(tas_annmean_D000,axis=0)
tas_mean_D012 = np.nanmean(tas_annmean_D012,axis=0)
tas_mean_D013 = np.nanmean(tas_annmean_D013,axis=0)
tas_mean_D014 = np.nanmean(tas_annmean_D014,axis=0)
tas_mean_D015 = np.nanmean(tas_annmean_D015,axis=0)
tas_mean_D016 = np.nanmean(tas_annmean_D016,axis=0)
tas_mean_D017 = np.nanmean(tas_annmean_D017,axis=0)
tas_mean_D018 = np.nanmean(tas_annmean_D018,axis=0)
tas_mean_D019 = np.nanmean(tas_annmean_D019,axis=0)
tas_mean_D020 = np.nanmean(tas_annmean_D020,axis=0)
tas_mean_D021 = np.nanmean(tas_annmean_D021,axis=0)
tas_mean_D022 = np.nanmean(tas_annmean_D022,axis=0)
tas_mean_D023 = np.nanmean(tas_annmean_D023,axis=0)
tas_mean_D024 = np.nanmean(tas_annmean_D024,axis=0)
tas_mean_D025 = np.nanmean(tas_annmean_D025,axis=0)
tas_mean_D027 = np.nanmean(tas_annmean_D027,axis=0)
tas_mean_D028 = np.nanmean(tas_annmean_D028,axis=0)
tas_mean_D029 = np.nanmean(tas_annmean_D029,axis=0)
tas_mean_D030 = np.nanmean(tas_annmean_D030,axis=0)

# Compute mean Arctic tas
tas_arctic_D000 = np.average(tas_mean_D000[lat>70.],weights=gridarea[lat>70.])
print(tas_arctic_D000)
tas_arctic_D012 = np.average(tas_mean_D012[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D013 = np.average(tas_mean_D013[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D014 = np.average(tas_mean_D014[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D015 = np.average(tas_mean_D015[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D016 = np.average(tas_mean_D016[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D017 = np.average(tas_mean_D017[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D018 = np.average(tas_mean_D018[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D019 = np.average(tas_mean_D019[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D020 = np.average(tas_mean_D020[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D021 = np.average(tas_mean_D021[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D022 = np.average(tas_mean_D022[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D023 = np.average(tas_mean_D023[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D024 = np.average(tas_mean_D024[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D025 = np.average(tas_mean_D025[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D027 = np.average(tas_mean_D027[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D028 = np.average(tas_mean_D028[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D029 = np.average(tas_mean_D029[lat>70.],weights=gridarea[lat>70.])
tas_arctic_D030 = np.average(tas_mean_D030[lat>70.],weights=gridarea[lat>70.])
tas_arctic_array = np.array([tas_arctic_D012-tas_arctic_D000,tas_arctic_D013-tas_arctic_D000,\
                        tas_arctic_D014-tas_arctic_D000,tas_arctic_D015-tas_arctic_D000,\
                        tas_arctic_D016-tas_arctic_D000,tas_arctic_D017-tas_arctic_D000,\
                        tas_arctic_D018-tas_arctic_D000,tas_arctic_D019-tas_arctic_D000,\
                        tas_arctic_D020-tas_arctic_D000,tas_arctic_D021-tas_arctic_D000,\
                        tas_arctic_D022-tas_arctic_D000,tas_arctic_D023-tas_arctic_D000,\
                        tas_arctic_D024-tas_arctic_D000,tas_arctic_D025-tas_arctic_D000,\
                        tas_arctic_D027-tas_arctic_D000,tas_arctic_D028-tas_arctic_D000,\
                        tas_arctic_D029-tas_arctic_D000,tas_arctic_D030-tas_arctic_D000])

# Compute mean March Arctic SIA
area_mean_march_D000 = np.nanmean(area_march_D000)
print(area_mean_march_D000)
area_mean_march_D012 = np.nanmean(area_march_D012)
area_mean_march_D013 = np.nanmean(area_march_D013)
area_mean_march_D014 = np.nanmean(area_march_D014)
area_mean_march_D015 = np.nanmean(area_march_D015)
area_mean_march_D016 = np.nanmean(area_march_D016)
area_mean_march_D017 = np.nanmean(area_march_D017)
area_mean_march_D018 = np.nanmean(area_march_D018)
area_mean_march_D019 = np.nanmean(area_march_D019)
area_mean_march_D020 = np.nanmean(area_march_D020)
area_mean_march_D021 = np.nanmean(area_march_D021)
area_mean_march_D022 = np.nanmean(area_march_D022)
area_mean_march_D023 = np.nanmean(area_march_D023)
area_mean_march_D024 = np.nanmean(area_march_D024)
area_mean_march_D025 = np.nanmean(area_march_D025)
area_mean_march_D027 = np.nanmean(area_march_D027)
area_mean_march_D028 = np.nanmean(area_march_D028)
area_mean_march_D029 = np.nanmean(area_march_D029)
area_mean_march_D030 = np.nanmean(area_march_D030)
area_mean_march_array = np.array([area_mean_march_D012-area_mean_march_D000,area_mean_march_D013-area_mean_march_D000,\
                        area_mean_march_D014-area_mean_march_D000,area_mean_march_D015-area_mean_march_D000,\
                        area_mean_march_D016-area_mean_march_D000,area_mean_march_D017-area_mean_march_D000,\
                        area_mean_march_D018-area_mean_march_D000,area_mean_march_D019-area_mean_march_D000,\
                        area_mean_march_D020-area_mean_march_D000,area_mean_march_D021-area_mean_march_D000,\
                        area_mean_march_D022-area_mean_march_D000,area_mean_march_D023-area_mean_march_D000,\
                        area_mean_march_D024-area_mean_march_D000,area_mean_march_D025-area_mean_march_D000,\
                        area_mean_march_D027-area_mean_march_D000,area_mean_march_D028-area_mean_march_D000,\
                        area_mean_march_D029-area_mean_march_D000,area_mean_march_D030-area_mean_march_D000])

# Compute mean September Arctic SIA
area_mean_sept_D000 = np.nanmean(area_sept_D000)
print(area_mean_sept_D000)
area_mean_sept_D012 = np.nanmean(area_sept_D012)
area_mean_sept_D013 = np.nanmean(area_sept_D013)
area_mean_sept_D014 = np.nanmean(area_sept_D014)
area_mean_sept_D015 = np.nanmean(area_sept_D015)
area_mean_sept_D016 = np.nanmean(area_sept_D016)
area_mean_sept_D017 = np.nanmean(area_sept_D017)
area_mean_sept_D018 = np.nanmean(area_sept_D018)
area_mean_sept_D019 = np.nanmean(area_sept_D019)
area_mean_sept_D020 = np.nanmean(area_sept_D020)
area_mean_sept_D021 = np.nanmean(area_sept_D021)
area_mean_sept_D022 = np.nanmean(area_sept_D022)
area_mean_sept_D023 = np.nanmean(area_sept_D023)
area_mean_sept_D024 = np.nanmean(area_sept_D024)
area_mean_sept_D025 = np.nanmean(area_sept_D025)
area_mean_sept_D027 = np.nanmean(area_sept_D027)
area_mean_sept_D028 = np.nanmean(area_sept_D028)
area_mean_sept_D029 = np.nanmean(area_sept_D029)
area_mean_sept_D030 = np.nanmean(area_sept_D030)
area_mean_sept_array = np.array([area_mean_sept_D012-area_mean_sept_D000,area_mean_sept_D013-area_mean_sept_D000,\
                        area_mean_sept_D014-area_mean_sept_D000,area_mean_sept_D015-area_mean_sept_D000,\
                        area_mean_sept_D016-area_mean_sept_D000,area_mean_sept_D017-area_mean_sept_D000,\
                        area_mean_sept_D018-area_mean_sept_D000,area_mean_sept_D019-area_mean_sept_D000,\
                        area_mean_sept_D020-area_mean_sept_D000,area_mean_sept_D021-area_mean_sept_D000,\
                        area_mean_sept_D022-area_mean_sept_D000,area_mean_sept_D023-area_mean_sept_D000,\
                        area_mean_sept_D024-area_mean_sept_D000,area_mean_sept_D025-area_mean_sept_D000,\
                        area_mean_sept_D027-area_mean_sept_D000,area_mean_sept_D028-area_mean_sept_D000,\
                        area_mean_sept_D029-area_mean_sept_D000,area_mean_sept_D030-area_mean_sept_D000])

# Compute mean March Arctic SIV
volume_mean_march_D000 = np.nanmean(volume_march_D000)
print(volume_mean_march_D000)
volume_mean_march_D012 = np.nanmean(volume_march_D012)
volume_mean_march_D013 = np.nanmean(volume_march_D013)
volume_mean_march_D014 = np.nanmean(volume_march_D014)
volume_mean_march_D015 = np.nanmean(volume_march_D015)
volume_mean_march_D016 = np.nanmean(volume_march_D016)
volume_mean_march_D017 = np.nanmean(volume_march_D017)
volume_mean_march_D018 = np.nanmean(volume_march_D018)
volume_mean_march_D019 = np.nanmean(volume_march_D019)
volume_mean_march_D020 = np.nanmean(volume_march_D020)
volume_mean_march_D021 = np.nanmean(volume_march_D021)
volume_mean_march_D022 = np.nanmean(volume_march_D022)
volume_mean_march_D023 = np.nanmean(volume_march_D023)
volume_mean_march_D024 = np.nanmean(volume_march_D024)
volume_mean_march_D025 = np.nanmean(volume_march_D025)
volume_mean_march_D027 = np.nanmean(volume_march_D027)
volume_mean_march_D028 = np.nanmean(volume_march_D028)
volume_mean_march_D029 = np.nanmean(volume_march_D029)
volume_mean_march_D030 = np.nanmean(volume_march_D030)
volume_mean_march_array = np.array([volume_mean_march_D012-volume_mean_march_D000,volume_mean_march_D013-volume_mean_march_D000,\
                        volume_mean_march_D014-volume_mean_march_D000,volume_mean_march_D015-volume_mean_march_D000,\
                        volume_mean_march_D016-volume_mean_march_D000,volume_mean_march_D017-volume_mean_march_D000,\
                        volume_mean_march_D018-volume_mean_march_D000,volume_mean_march_D019-volume_mean_march_D000,\
                        volume_mean_march_D020-volume_mean_march_D000,volume_mean_march_D021-volume_mean_march_D000,\
                        volume_mean_march_D022-volume_mean_march_D000,volume_mean_march_D023-volume_mean_march_D000,\
                        volume_mean_march_D024-volume_mean_march_D000,volume_mean_march_D025-volume_mean_march_D000,\
                        volume_mean_march_D027-volume_mean_march_D000,volume_mean_march_D028-volume_mean_march_D000,\
                        volume_mean_march_D029-volume_mean_march_D000,volume_mean_march_D030-volume_mean_march_D000])

# Compute mean March Arctic SIV
volume_mean_sept_D000 = np.nanmean(volume_sept_D000)
print(volume_mean_sept_D000)
volume_mean_sept_D012 = np.nanmean(volume_sept_D012)
volume_mean_sept_D013 = np.nanmean(volume_sept_D013)
volume_mean_sept_D014 = np.nanmean(volume_sept_D014)
volume_mean_sept_D015 = np.nanmean(volume_sept_D015)
volume_mean_sept_D016 = np.nanmean(volume_sept_D016)
volume_mean_sept_D017 = np.nanmean(volume_sept_D017)
volume_mean_sept_D018 = np.nanmean(volume_sept_D018)
volume_mean_sept_D019 = np.nanmean(volume_sept_D019)
volume_mean_sept_D020 = np.nanmean(volume_sept_D020)
volume_mean_sept_D021 = np.nanmean(volume_sept_D021)
volume_mean_sept_D022 = np.nanmean(volume_sept_D022)
volume_mean_sept_D023 = np.nanmean(volume_sept_D023)
volume_mean_sept_D024 = np.nanmean(volume_sept_D024)
volume_mean_sept_D025 = np.nanmean(volume_sept_D025)
volume_mean_sept_D027 = np.nanmean(volume_sept_D027)
volume_mean_sept_D028 = np.nanmean(volume_sept_D028)
volume_mean_sept_D029 = np.nanmean(volume_sept_D029)
volume_mean_sept_D030 = np.nanmean(volume_sept_D030)
volume_mean_sept_array = np.array([volume_mean_sept_D012-volume_mean_sept_D000,volume_mean_sept_D013-volume_mean_sept_D000,\
                        volume_mean_sept_D014-volume_mean_sept_D000,volume_mean_sept_D015-volume_mean_sept_D000,\
                        volume_mean_sept_D016-volume_mean_sept_D000,volume_mean_sept_D017-volume_mean_sept_D000,\
                        volume_mean_sept_D018-volume_mean_sept_D000,volume_mean_sept_D019-volume_mean_sept_D000,\
                        volume_mean_sept_D020-volume_mean_sept_D000,volume_mean_sept_D021-volume_mean_sept_D000,\
                        volume_mean_sept_D022-volume_mean_sept_D000,volume_mean_sept_D023-volume_mean_sept_D000,\
                        volume_mean_sept_D024-volume_mean_sept_D000,volume_mean_sept_D025-volume_mean_sept_D000,\
                        volume_mean_sept_D027-volume_mean_sept_D000,volume_mean_sept_D028-volume_mean_sept_D000,\
                        volume_mean_sept_D029-volume_mean_sept_D000,volume_mean_sept_D030-volume_mean_sept_D000])

# Compute OHT-SIA correlations
R_arctic_march,p_arctic_march = pearsonr(oht_mean_total_array,area_mean_march_array)
R_arctic_atl_march,p_arctic_atl_march = pearsonr(oht_mean_total_array[0:9],area_mean_march_array[0:9])
R_arctic_pac_march,p_arctic_pac_march = pearsonr(oht_mean_total_array[9::],area_mean_march_array[9::])
R_arctic_sept,p_arctic_sept = pearsonr(oht_mean_total_array,area_mean_sept_array)
R_arctic_atl_sept,p_arctic_atl_sept = pearsonr(oht_mean_total_array[0:9],area_mean_sept_array[0:9])
R_arctic_pac_sept,p_arctic_pac_sept = pearsonr(oht_mean_total_array[9::],area_mean_sept_array[9::])

# Compute OHT-SIV correlations
R_arctic_march_siv,p_arctic_march_siv = pearsonr(oht_mean_total_array,volume_mean_march_array)
R_arctic_atl_march_siv,p_arctic_atl_march_siv = pearsonr(oht_mean_total_array[0:9],volume_mean_march_array[0:9])
R_arctic_pac_march_siv,p_arctic_pac_march_siv = pearsonr(oht_mean_total_array[9::],volume_mean_march_array[9::])
R_arctic_sept_siv,p_arctic_sept_siv = pearsonr(oht_mean_total_array,volume_mean_sept_array)
R_arctic_atl_sept_siv,p_arctic_atl_sept_siv = pearsonr(oht_mean_total_array[0:9],volume_mean_sept_array[0:9])
R_arctic_pac_sept_siv,p_arctic_pac_sept_siv = pearsonr(oht_mean_total_array[9::],volume_mean_sept_array[9::])

# Compute OHT-SIA slopes
a_arctic_march,b_arctic_march,sd_arctic_march,sig_arctic_march = compute_slope(18,oht_mean_total_array,area_mean_march_array)
a_arctic_atl_march,b_arctic_atl_march,sd_arctic_atl_march,sig_arctic_atl_march = compute_slope(9,oht_mean_total_array[0:9],area_mean_march_array[0:9])
a_arctic_pac_march,b_arctic_pac_march,sd_arctic_pac_march,sig_arctic_pac_march = compute_slope(9,oht_mean_total_array[9::],area_mean_march_array[9::])
a_arctic_sept,b_arctic_sept,sd_arctic_sept,sig_arctic_sept = compute_slope(18,oht_mean_total_array,area_mean_sept_array)
a_arctic_atl_sept,b_arctic_atl_sept,sd_arctic_atl_sept,sig_arctic_atl_sept = compute_slope(9,oht_mean_total_array[0:9],area_mean_sept_array[0:9])
a_arctic_pac_sept,b_arctic_pac_sept,sd_arctic_pac_sept,sig_arctic_pac_sept = compute_slope(9,oht_mean_total_array[9::],area_mean_sept_array[9::])

# Compute OHT-SIV slopes
a_arctic_march_siv,b_arctic_march_siv,sd_arctic_march_siv,sig_arctic_march_siv = compute_slope(18,oht_mean_total_array,volume_mean_march_array)
a_arctic_atl_march_siv,b_arctic_atl_march_siv,sd_arctic_atl_march_siv,sig_arctic_atl_march_siv = compute_slope(9,oht_mean_total_array[0:9],volume_mean_march_array[0:9])
a_arctic_pac_march_siv,b_arctic_pac_march_siv,sd_arctic_pac_march_siv,sig_arctic_pac_march_siv = compute_slope(9,oht_mean_total_array[9::],volume_mean_march_array[9::])
a_arctic_sept_siv,b_arctic_sept_siv,sd_arctic_sept_siv,sig_arctic_sept_siv = compute_slope(18,oht_mean_total_array,volume_mean_sept_array)
a_arctic_atl_sept_siv,b_arctic_atl_sept_siv,sd_arctic_atl_sept_siv,sig_arctic_atl_sept_siv = compute_slope(9,oht_mean_total_array[0:9],volume_mean_sept_array[0:9])
a_arctic_pac_sept_siv,b_arctic_pac_sept_siv,sd_arctic_pac_sept_siv,sig_arctic_pac_sept_siv = compute_slope(9,oht_mean_total_array[9::],volume_mean_sept_array[9::])
    
# Compute tas-SIA slopes
a_arctic_march_tas,b_arctic_march_tas,sd_arctic_march_tas,sig_arctic_march_tas = compute_slope(18,tas_arctic_array,area_mean_march_array)
a_arctic_atl_march_tas,b_arctic_atl_march_tas,sd_arctic_atl_march_tas,sig_arctic_atl_march_tas = compute_slope(9,tas_arctic_array[0:9],area_mean_march_array[0:9])
a_arctic_pac_march_tas,b_arctic_pac_march_tas,sd_arctic_pac_march_tas,sig_arctic_pac_march_tas = compute_slope(9,tas_arctic_array[9::],area_mean_march_array[9::])
a_arctic_sept_tas,b_arctic_sept_tas,sd_arctic_sept_tas,sig_arctic_sept_tas = compute_slope(18,tas_arctic_array,area_mean_sept_array)
a_arctic_atl_sept_tas,b_arctic_atl_sept_tas,sd_arctic_atl_sept_tas,sig_arctic_atl_sept_tas = compute_slope(9,tas_arctic_array[0:9],area_mean_sept_array[0:9])
a_arctic_pac_sept_tas,b_arctic_pac_sept_tas,sd_arctic_pac_sept_tas,sig_arctic_pac_sept_tas = compute_slope(9,tas_arctic_array[9::],area_mean_sept_array[9::])

# Compute tas-SIV slopes
a_arctic_march_tas_siv,b_arctic_march_tas_siv,sd_arctic_march_tas_siv,sig_arctic_march_tas_siv = compute_slope(18,tas_arctic_array,volume_mean_march_array)
a_arctic_atl_march_tas_siv,b_arctic_atl_march_tas_siv,sd_arctic_atl_march_tas_siv,sig_arctic_atl_march_tas_siv = compute_slope(9,tas_arctic_array[0:9],volume_mean_march_array[0:9])
a_arctic_pac_march_tas_siv,b_arctic_pac_march_tas_siv,sd_arctic_pac_march_tas_siv,sig_arctic_pac_march_tas_siv = compute_slope(9,tas_arctic_array[9::],volume_mean_march_array[9::])
a_arctic_sept_tas_siv,b_arctic_sept_tas_siv,sd_arctic_sept_tas_siv,sig_arctic_sept_tas_siv = compute_slope(18,tas_arctic_array,volume_mean_sept_array)
a_arctic_atl_sept_tas_siv,b_arctic_atl_sept_tas_siv,sd_arctic_atl_sept_tas_siv,sig_arctic_atl_sept_tas_siv = compute_slope(9,tas_arctic_array[0:9],volume_mean_sept_array[0:9])
a_arctic_pac_sept_tas_siv,b_arctic_pac_sept_tas_siv,sd_arctic_pac_sept_tas_siv,sig_arctic_pac_sept_tas_siv = compute_slope(9,tas_arctic_array[9::],volume_mean_sept_array[9::])


# Fig. 9 - OHT-SIA/SIV Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(2,2,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.18,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Total OHT - March Arctic SIA
ax[0,0].plot(oht_mean_total_D013-oht_mean_total_D000,area_mean_march_D013-area_mean_march_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D012-oht_mean_total_D000,area_mean_march_D012-area_mean_march_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D014-oht_mean_total_D000,area_mean_march_D014-area_mean_march_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D016-oht_mean_total_D000,area_mean_march_D016-area_mean_march_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D015-oht_mean_total_D000,area_mean_march_D015-area_mean_march_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D017-oht_mean_total_D000,area_mean_march_D017-area_mean_march_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D019-oht_mean_total_D000,area_mean_march_D019-area_mean_march_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D018-oht_mean_total_D000,area_mean_march_D018-area_mean_march_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D020-oht_mean_total_D000,area_mean_march_D020-area_mean_march_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D027-oht_mean_total_D000,area_mean_march_D027-area_mean_march_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D021-oht_mean_total_D000,area_mean_march_D021-area_mean_march_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D028-oht_mean_total_D000,area_mean_march_D028-area_mean_march_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D029-oht_mean_total_D000,area_mean_march_D029-area_mean_march_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D022-oht_mean_total_D000,area_mean_march_D022-area_mean_march_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D030-oht_mean_total_D000,area_mean_march_D030-area_mean_march_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[0,0].plot(oht_mean_total_D024-oht_mean_total_D000,area_mean_march_D024-area_mean_march_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[0,0].plot(oht_mean_total_D023-oht_mean_total_D000,area_mean_march_D023-area_mean_march_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[0,0].plot(oht_mean_total_D025-oht_mean_total_D000,area_mean_march_D025-area_mean_march_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[0,0].set_ylabel('$\Delta$SIA$_{Arctic,March}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,0].set_xticks(np.arange(0, 140.1, 20))
ax[0,0].set_yticks(np.arange(-5, 0.3, 1))
ax[0,0].axis([0, 140, -5, 0.3])
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].annotate('$a$ = '+str(int(np.round(a_arctic_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Total OHT - September Arctic SIA
ax[0,1].plot(oht_mean_total_D013-oht_mean_total_D000,area_mean_sept_D013-area_mean_sept_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D012-oht_mean_total_D000,area_mean_sept_D012-area_mean_sept_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D014-oht_mean_total_D000,area_mean_sept_D014-area_mean_sept_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D016-oht_mean_total_D000,area_mean_sept_D016-area_mean_sept_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D015-oht_mean_total_D000,area_mean_sept_D015-area_mean_sept_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D017-oht_mean_total_D000,area_mean_sept_D017-area_mean_sept_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D019-oht_mean_total_D000,area_mean_sept_D019-area_mean_sept_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D018-oht_mean_total_D000,area_mean_sept_D018-area_mean_sept_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D020-oht_mean_total_D000,area_mean_sept_D020-area_mean_sept_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D027-oht_mean_total_D000,area_mean_sept_D027-area_mean_sept_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D021-oht_mean_total_D000,area_mean_sept_D021-area_mean_sept_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D028-oht_mean_total_D000,area_mean_sept_D028-area_mean_sept_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D029-oht_mean_total_D000,area_mean_sept_D029-area_mean_sept_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D022-oht_mean_total_D000,area_mean_sept_D022-area_mean_sept_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D030-oht_mean_total_D000,area_mean_sept_D030-area_mean_sept_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[0,1].plot(oht_mean_total_D024-oht_mean_total_D000,area_mean_sept_D024-area_mean_sept_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[0,1].plot(oht_mean_total_D023-oht_mean_total_D000,area_mean_sept_D023-area_mean_sept_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[0,1].plot(oht_mean_total_D025-oht_mean_total_D000,area_mean_sept_D025-area_mean_sept_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[0,1].set_ylabel('$\Delta$SIA$_{Arctic,September}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,1].set_xticks(np.arange(0, 140.1, 20))
ax[0,1].set_yticks(np.arange(-5, 0.3, 1))
ax[0,1].axis([0, 140, -5, 0.3])
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept*1.e6,0)))+' km$^2$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Total OHT - March Arctic SIV
ax[1,0].plot(oht_mean_total_D013-oht_mean_total_D000,volume_mean_march_D013-volume_mean_march_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D012-oht_mean_total_D000,volume_mean_march_D012-volume_mean_march_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D014-oht_mean_total_D000,volume_mean_march_D014-volume_mean_march_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D016-oht_mean_total_D000,volume_mean_march_D016-volume_mean_march_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D015-oht_mean_total_D000,volume_mean_march_D015-volume_mean_march_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D017-oht_mean_total_D000,volume_mean_march_D017-volume_mean_march_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D019-oht_mean_total_D000,volume_mean_march_D019-volume_mean_march_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D018-oht_mean_total_D000,volume_mean_march_D018-volume_mean_march_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D020-oht_mean_total_D000,volume_mean_march_D020-volume_mean_march_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D027-oht_mean_total_D000,volume_mean_march_D027-volume_mean_march_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D021-oht_mean_total_D000,volume_mean_march_D021-volume_mean_march_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D028-oht_mean_total_D000,volume_mean_march_D028-volume_mean_march_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D029-oht_mean_total_D000,volume_mean_march_D029-volume_mean_march_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D022-oht_mean_total_D000,volume_mean_march_D022-volume_mean_march_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D030-oht_mean_total_D000,volume_mean_march_D030-volume_mean_march_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[1,0].plot(oht_mean_total_D024-oht_mean_total_D000,volume_mean_march_D024-volume_mean_march_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[1,0].plot(oht_mean_total_D023-oht_mean_total_D000,volume_mean_march_D023-volume_mean_march_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[1,0].plot(oht_mean_total_D025-oht_mean_total_D000,volume_mean_march_D025-volume_mean_march_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[1,0].set_xlabel('$\Delta$OHT$_{Arctic}$  (TW)',fontsize=20)
ax[1,0].set_ylabel('$\Delta$SIV$_{Arctic,March}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,0].set_xticks(np.arange(0, 140.1, 20))
ax[1,0].set_yticks(np.arange(-15, 1, 3))
ax[1,0].axis([0, 140, -15, 1])
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].annotate('$a$ = '+str(int(np.round(a_arctic_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Total OHT - September Arctic SIV
ax[1,1].plot(oht_mean_total_D013-oht_mean_total_D000,volume_mean_sept_D013-volume_mean_sept_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D012-oht_mean_total_D000,volume_mean_sept_D012-volume_mean_sept_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D014-oht_mean_total_D000,volume_mean_sept_D014-volume_mean_sept_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D016-oht_mean_total_D000,volume_mean_sept_D016-volume_mean_sept_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D015-oht_mean_total_D000,volume_mean_sept_D015-volume_mean_sept_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D017-oht_mean_total_D000,volume_mean_sept_D017-volume_mean_sept_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D019-oht_mean_total_D000,volume_mean_sept_D019-volume_mean_sept_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D018-oht_mean_total_D000,volume_mean_sept_D018-volume_mean_sept_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D020-oht_mean_total_D000,volume_mean_sept_D020-volume_mean_sept_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D027-oht_mean_total_D000,volume_mean_sept_D027-volume_mean_sept_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D021-oht_mean_total_D000,volume_mean_sept_D021-volume_mean_sept_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D028-oht_mean_total_D000,volume_mean_sept_D028-volume_mean_sept_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D029-oht_mean_total_D000,volume_mean_sept_D029-volume_mean_sept_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D022-oht_mean_total_D000,volume_mean_sept_D022-volume_mean_sept_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D030-oht_mean_total_D000,volume_mean_sept_D030-volume_mean_sept_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[1,1].plot(oht_mean_total_D024-oht_mean_total_D000,volume_mean_sept_D024-volume_mean_sept_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[1,1].plot(oht_mean_total_D023-oht_mean_total_D000,volume_mean_sept_D023-volume_mean_sept_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[1,1].plot(oht_mean_total_D025-oht_mean_total_D000,volume_mean_sept_D025-volume_mean_sept_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[1,1].set_xlabel('$\Delta$OHT$_{Arctic}$  (TW)',fontsize=20)
ax[1,1].set_ylabel('$\Delta$SIV$_{Arctic,September}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,1].set_xticks(np.arange(0, 140.1, 20))
ax[1,1].set_yticks(np.arange(-15, 1, 3))
ax[1,1].axis([0, 140, -15, 1])
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept_siv*1.e3,0)))+' km$^3$ TW$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')
ax[1,1].legend(shadow=True,frameon=False,fontsize=14,bbox_to_anchor=(0.8,-0.2),ncol=6)

# Adjust plot
#plt.subplots_adjust(wspace=0.3)

if save_fig == True:
    fig.savefig(dir_fig + 'fig9.png')


# Supp. Fig. 9b - tas-SIA/SIV Scatter plots (changes compared to CTRL)
fig,ax = plt.subplots(2,2,figsize=(15,12))
fig.subplots_adjust(left=0.1,bottom=0.18,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Arctic tas - March Arctic SIA
ax[0,0].plot(tas_arctic_D013-tas_arctic_D000,area_mean_march_D013-area_mean_march_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[0,0].plot(tas_arctic_D012-tas_arctic_D000,area_mean_march_D012-area_mean_march_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[0,0].plot(tas_arctic_D014-tas_arctic_D000,area_mean_march_D014-area_mean_march_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[0,0].plot(tas_arctic_D016-tas_arctic_D000,area_mean_march_D016-area_mean_march_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[0,0].plot(tas_arctic_D015-tas_arctic_D000,area_mean_march_D015-area_mean_march_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[0,0].plot(tas_arctic_D017-tas_arctic_D000,area_mean_march_D017-area_mean_march_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[0,0].plot(tas_arctic_D019-tas_arctic_D000,area_mean_march_D019-area_mean_march_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[0,0].plot(tas_arctic_D018-tas_arctic_D000,area_mean_march_D018-area_mean_march_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[0,0].plot(tas_arctic_D020-tas_arctic_D000,area_mean_march_D020-area_mean_march_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[0,0].plot(tas_arctic_D027-tas_arctic_D000,area_mean_march_D027-area_mean_march_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[0,0].plot(tas_arctic_D021-tas_arctic_D000,area_mean_march_D021-area_mean_march_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[0,0].plot(tas_arctic_D028-tas_arctic_D000,area_mean_march_D028-area_mean_march_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[0,0].plot(tas_arctic_D029-tas_arctic_D000,area_mean_march_D029-area_mean_march_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[0,0].plot(tas_arctic_D022-tas_arctic_D000,area_mean_march_D022-area_mean_march_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[0,0].plot(tas_arctic_D030-tas_arctic_D000,area_mean_march_D030-area_mean_march_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[0,0].plot(tas_arctic_D024-tas_arctic_D000,area_mean_march_D024-area_mean_march_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[0,0].plot(tas_arctic_D023-tas_arctic_D000,area_mean_march_D023-area_mean_march_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[0,0].plot(tas_arctic_D025-tas_arctic_D000,area_mean_march_D025-area_mean_march_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[0,0].set_ylabel('$\Delta$SIA$_{Arctic,March}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,0].set_xticks(np.arange(0, 6.1, 1))
ax[0,0].set_yticks(np.arange(-5, 0.3, 1))
ax[0,0].axis([-0.2, 6, -5, 0.3])
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].annotate('$a$ = '+str(int(np.round(a_arctic_march_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Arctic tas - September Arctic SIA
ax[0,1].plot(tas_arctic_D013-tas_arctic_D000,area_mean_sept_D013-area_mean_sept_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[0,1].plot(tas_arctic_D012-tas_arctic_D000,area_mean_sept_D012-area_mean_sept_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[0,1].plot(tas_arctic_D014-tas_arctic_D000,area_mean_sept_D014-area_mean_sept_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[0,1].plot(tas_arctic_D016-tas_arctic_D000,area_mean_sept_D016-area_mean_sept_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[0,1].plot(tas_arctic_D015-tas_arctic_D000,area_mean_sept_D015-area_mean_sept_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[0,1].plot(tas_arctic_D017-tas_arctic_D000,area_mean_sept_D017-area_mean_sept_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[0,1].plot(tas_arctic_D019-tas_arctic_D000,area_mean_sept_D019-area_mean_sept_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[0,1].plot(tas_arctic_D018-tas_arctic_D000,area_mean_sept_D018-area_mean_sept_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[0,1].plot(tas_arctic_D020-tas_arctic_D000,area_mean_sept_D020-area_mean_sept_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[0,1].plot(tas_arctic_D027-tas_arctic_D000,area_mean_sept_D027-area_mean_sept_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[0,1].plot(tas_arctic_D021-tas_arctic_D000,area_mean_sept_D021-area_mean_sept_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[0,1].plot(tas_arctic_D028-tas_arctic_D000,area_mean_sept_D028-area_mean_sept_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[0,1].plot(tas_arctic_D029-tas_arctic_D000,area_mean_sept_D029-area_mean_sept_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[0,1].plot(tas_arctic_D022-tas_arctic_D000,area_mean_sept_D022-area_mean_sept_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[0,1].plot(tas_arctic_D030-tas_arctic_D000,area_mean_sept_D030-area_mean_sept_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[0,1].plot(tas_arctic_D024-tas_arctic_D000,area_mean_sept_D024-area_mean_sept_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[0,1].plot(tas_arctic_D023-tas_arctic_D000,area_mean_sept_D023-area_mean_sept_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[0,1].plot(tas_arctic_D025-tas_arctic_D000,area_mean_sept_D025-area_mean_sept_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[0,1].set_ylabel('$\Delta$SIA$_{Arctic,September}$ (10$^6$ km$^2$)',fontsize=20)
ax[0,1].set_xticks(np.arange(0, 6.1, 1))
ax[0,1].set_yticks(np.arange(-5, 0.3, 1))
ax[0,1].axis([-0.2, 6, -5, 0.3])
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept_tas*1.e6,0)))+' km$^2$ $^\circ$C$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Arctic tas - March Arctic SIV
ax[1,0].plot(tas_arctic_D013-tas_arctic_D000,volume_mean_march_D013-volume_mean_march_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[1,0].plot(tas_arctic_D012-tas_arctic_D000,volume_mean_march_D012-volume_mean_march_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[1,0].plot(tas_arctic_D014-tas_arctic_D000,volume_mean_march_D014-volume_mean_march_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[1,0].plot(tas_arctic_D016-tas_arctic_D000,volume_mean_march_D016-volume_mean_march_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[1,0].plot(tas_arctic_D015-tas_arctic_D000,volume_mean_march_D015-volume_mean_march_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[1,0].plot(tas_arctic_D017-tas_arctic_D000,volume_mean_march_D017-volume_mean_march_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[1,0].plot(tas_arctic_D019-tas_arctic_D000,volume_mean_march_D019-volume_mean_march_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[1,0].plot(tas_arctic_D018-tas_arctic_D000,volume_mean_march_D018-volume_mean_march_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[1,0].plot(tas_arctic_D020-tas_arctic_D000,volume_mean_march_D020-volume_mean_march_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[1,0].plot(tas_arctic_D027-tas_arctic_D000,volume_mean_march_D027-volume_mean_march_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[1,0].plot(tas_arctic_D021-tas_arctic_D000,volume_mean_march_D021-volume_mean_march_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[1,0].plot(tas_arctic_D028-tas_arctic_D000,volume_mean_march_D028-volume_mean_march_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[1,0].plot(tas_arctic_D029-tas_arctic_D000,volume_mean_march_D029-volume_mean_march_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[1,0].plot(tas_arctic_D022-tas_arctic_D000,volume_mean_march_D022-volume_mean_march_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[1,0].plot(tas_arctic_D030-tas_arctic_D000,volume_mean_march_D030-volume_mean_march_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[1,0].plot(tas_arctic_D024-tas_arctic_D000,volume_mean_march_D024-volume_mean_march_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[1,0].plot(tas_arctic_D023-tas_arctic_D000,volume_mean_march_D023-volume_mean_march_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[1,0].plot(tas_arctic_D025-tas_arctic_D000,volume_mean_march_D025-volume_mean_march_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[1,0].set_xlabel('$\Delta$T$_{2m,Arctic}$  ($^\circ$C)',fontsize=20)
ax[1,0].set_ylabel('$\Delta$SIV$_{Arctic,March}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,0].set_xticks(np.arange(0, 6.1, 1))
ax[1,0].set_yticks(np.arange(-15, 1, 3))
ax[1,0].axis([-0.2, 6, -15, 1])
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].annotate('$a$ = '+str(int(np.round(a_arctic_march_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_march_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_march_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Arctic tas - September Arctic SIV
ax[1,1].plot(tas_arctic_D013-tas_arctic_D000,volume_mean_sept_D013-volume_mean_sept_D000,'o',color='lightcoral',markersize=8,label='ATL1+1$^\circ$C')
ax[1,1].plot(tas_arctic_D012-tas_arctic_D000,volume_mean_sept_D012-volume_mean_sept_D000,'ro',markersize=8,label='ATL1+3$^\circ$C')
ax[1,1].plot(tas_arctic_D014-tas_arctic_D000,volume_mean_sept_D014-volume_mean_sept_D000,'o',color='purple',markersize=8,label='ATL1+5$^\circ$C')
ax[1,1].plot(tas_arctic_D016-tas_arctic_D000,volume_mean_sept_D016-volume_mean_sept_D000,'o',color='lightblue',markersize=8,label='ATL2+1$^\circ$C')
ax[1,1].plot(tas_arctic_D015-tas_arctic_D000,volume_mean_sept_D015-volume_mean_sept_D000,'bo',markersize=8,label='ATL2+3$^\circ$C')
ax[1,1].plot(tas_arctic_D017-tas_arctic_D000,volume_mean_sept_D017-volume_mean_sept_D000,'o',color='darkblue',markersize=8,label='ATL2+5$^\circ$C')
ax[1,1].plot(tas_arctic_D019-tas_arctic_D000,volume_mean_sept_D019-volume_mean_sept_D000,'o',color='lightgreen',markersize=8,label='ATL3+1$^\circ$C')
ax[1,1].plot(tas_arctic_D018-tas_arctic_D000,volume_mean_sept_D018-volume_mean_sept_D000,'go',markersize=8,label='ATL3+3$^\circ$C')
ax[1,1].plot(tas_arctic_D020-tas_arctic_D000,volume_mean_sept_D020-volume_mean_sept_D000,'o',color='darkgreen',markersize=8,label='ATL3+5$^\circ$C')
ax[1,1].plot(tas_arctic_D027-tas_arctic_D000,volume_mean_sept_D027-volume_mean_sept_D000,'x',color='lightcoral',markersize=8,label='PAC1+1$^\circ$C')
ax[1,1].plot(tas_arctic_D021-tas_arctic_D000,volume_mean_sept_D021-volume_mean_sept_D000,'rx',markersize=8,label='PAC1+3$^\circ$C')
ax[1,1].plot(tas_arctic_D028-tas_arctic_D000,volume_mean_sept_D028-volume_mean_sept_D000,'x',color='purple',markersize=8,label='PAC1+5$^\circ$C')
ax[1,1].plot(tas_arctic_D029-tas_arctic_D000,volume_mean_sept_D029-volume_mean_sept_D000,'x',color='lightblue',markersize=8,label='PAC1+1$^\circ$C')
ax[1,1].plot(tas_arctic_D022-tas_arctic_D000,volume_mean_sept_D022-volume_mean_sept_D000,'bx',markersize=8,label='PAC2+3$^\circ$C')
ax[1,1].plot(tas_arctic_D030-tas_arctic_D000,volume_mean_sept_D030-volume_mean_sept_D000,'x',color='darkblue',markersize=8,label='PAC2+5$^\circ$C')
ax[1,1].plot(tas_arctic_D024-tas_arctic_D000,volume_mean_sept_D024-volume_mean_sept_D000,'x',color='lightgreen',markersize=8,label='PAC3+1$^\circ$C')
ax[1,1].plot(tas_arctic_D023-tas_arctic_D000,volume_mean_sept_D023-volume_mean_sept_D000,'gx',markersize=8,label='PAC3+3$^\circ$C')
ax[1,1].plot(tas_arctic_D025-tas_arctic_D000,volume_mean_sept_D025-volume_mean_sept_D000,'x',color='darkgreen',markersize=8,label='PAC3+5$^\circ$C')
ax[1,1].set_xlabel('$\Delta$T$_{2m,Arctic}$  ($^\circ$C)',fontsize=20)
ax[1,1].set_ylabel('$\Delta$SIV$_{Arctic,September}$ (10$^3$ km$^3$)',fontsize=20)
ax[1,1].set_xticks(np.arange(0, 6.1, 1))
ax[1,1].set_yticks(np.arange(-15, 1, 3))
ax[1,1].axis([-0.2, 6, -15, 1])
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].annotate('$a$ = '+str(int(np.round(a_arctic_sept_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.85),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{atl}$ = '+str(int(np.round(a_arctic_atl_sept_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.78),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].annotate('$a_{pac}$ = '+str(int(np.round(a_arctic_pac_sept_tas_siv*1.e3,0)))+' km$^3$ $^\circ$C$^{-1}$',xy=(0.97,0.71),xycoords='axes fraction',fontsize=15,xytext=(-5,5),textcoords='offset points',ha='right',va='bottom')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')
ax[1,1].legend(shadow=True,frameon=False,fontsize=14,bbox_to_anchor=(0.8,-0.2),ncol=6)

if save_fig == True:
    fig.savefig(dir_fig + 'fig9b.png')
