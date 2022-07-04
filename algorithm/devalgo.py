#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 09:45:56 2022

@author: ssk
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:35:33 2022

@author: ssk
"""

# Calculate regression coefficients based on Alerskans et al., 2020,
# "Construction of a climate data record of sea surface temperature from
# passive microwave measurements"


import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
from scipy.interpolate import interp1d, interp2d
import pytz
from skyfield import api
from skyfield import almanac
from sklearn import linear_model
import matplotlib.pyplot as plt
import scipy

# Citation: see https://scikit-learn.org/stable/about.html#citing-scikit-learn

# =========
# FUNCTIONS
# =========
def relative_angle(azimuth, u, v):
    D2R = np.pi/180
    phi_w = 90 - (np.arctan2(u, v))*D2R
    phi_w[phi_w < 0] = phi_w[phi_w < 0] + 360
    phi_rel = azimuth - phi_w
    phi_rel[phi_rel < 0] = phi_rel[phi_rel < 0] + 360
    return phi_rel


def regression(X,Y):
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)
    intercept = regr.intercept_
    coeffs = regr.coef_
    return intercept, coeffs


def calculate_coeffs_stage_1_ws(data):
    # Construct input data, based on Alerskans 2020, eq.1
    tb_1v   =  data['1vpol'] - 150
    tb_1v2  = (data['1vpol'] - 150)**2
    tb_1h   =  data['1hpol'] - 150
    tb_1h2  = (data['1hpol'] - 150)**2
    tb_6v   =  data['6vpol'] - 150
    tb_6v2  = (data['6vpol'] - 150)**2
    tb_6h   =  data['6hpol'] - 150
    tb_6h2  = (data['6hpol'] - 150)**2
    tb_10v  =  data['10vpol'] - 150
    tb_10v2 = (data['10vpol'] - 150)**2
    tb_10h  =  data['10hpol'] - 150
    tb_10h2 = (data['10hpol'] - 150)**2
    tb_18v  =  data['18vpol'] - 150
    tb_18v2 = (data['18vpol'] - 150)**2
    tb_18h  =  data['18hpol'] - 150
    tb_18h2 = (data['18hpol'] - 150)**2
    tb_36v  =  data['36vpol'] - 150
    tb_36v2 = (data['36vpol'] - 150)**2
    tb_36h  =  data['36hpol'] - 150
    tb_36h2 = (data['36hpol'] - 150)**2
    theta   = data['look_angle'] - 55
    
    X = np.stack((tb_1v,tb_1v2,tb_1h,tb_1h2,tb_6v,tb_6v2,tb_6h,tb_6h2,tb_10v,tb_10v2,tb_10h,tb_10h2,tb_18v,tb_18v2,tb_18h,tb_18h2,tb_36v,tb_36v2,tb_36h,tb_36h2,theta),axis=1)
    Y = data['ws']
    
    # Calculate linear regression coefficients
    intercept, coeffs = regression(X,Y)
    
    # Save coefficients
    coeffs_all = np.append(intercept,coeffs)
    np.save('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws_stage_1.npy',coeffs_all)


def calculate_coeffs_stage_2_ws(data, WSa):
    # Attach WSa to data DataFrame
    data['WSa'] = WSa
    # Wind speed bins
    bins_ws = np.arange(0,20,1, dtype=int)
    
    # Loop through bins
    for bin_ws in bins_ws:
        for bin_ssta in bins_ssta:
            # Find the data that belongs to the wind speed bin and SSTa bin
            data_filt = data.loc[(data['ws'] > bin_ws) & (data['ws'] < bin_ws+1)]
            
            # Construct input data, based in Alerskans 2020, eq.5
            tb_1v   =  data_filt['1vpol'] - 150
            tb_1v2  = (data_filt['1vpol'] - 150)**2
            tb_1h   =  data_filt['1hpol'] - 150
            tb_1h2  = (data_filt['1hpol'] - 150)**2
            tb_6v   =  data_filt['6vpol'] - 150
            tb_6v2  = (data_filt['6vpol'] - 150)**2
            tb_6h   =  data_filt['6hpol'] - 150
            tb_6h2  = (data_filt['6hpol'] - 150)**2
            tb_10v  =  data_filt['10vpol'] - 150
            tb_10v2 = (data_filt['10vpol'] - 150)**2
            tb_10h  =  data_filt['10hpol'] - 150
            tb_10h2 = (data_filt['10hpol'] - 150)**2
            tb_18v  =  data_filt['18vpol'] - 150
            tb_18v2 = (data_filt['18vpol'] - 150)**2
            tb_18h  =  data_filt['18hpol'] - 150
            tb_18h2 = (data_filt['18hpol'] - 150)**2
            tb_36v  =  data_filt['36vpol'] - 150
            tb_36v2 = (data_filt['36vpol'] - 150)**2
            tb_36h  =  data_filt['36hpol'] - 150
            tb_36h2 = (data_filt['36hpol'] - 150)**2
            theta   = data_filt['look_angle'] - 55
                
            X = np.stack((tb_1v,tb_1v2,tb_1h,tb_1h2,tb_6v,tb_6v2,tb_6h,tb_6h2,tb_10v,tb_10v2,tb_10h,tb_10h2,tb_18v,tb_18v2,tb_18h,tb_18h2,tb_36v,tb_36v2,tb_36h,tb_36h2,theta),axis=1)
            Y = data_filt['SSTa']
            
            # the try...continue structure is to account for cases where we don't have data in a particular bin, then the code moves on
            try:
                # Calculate linear regression coefficients
                intercept, coeffs = regression(X,Y)
                
                # Save coefficients
                coeffs_all = np.append(intercept,coeffs)
                np.save('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_2/coeffs_ws'+str(bin_ws)+'.npy',coeffs_all)
            except:
                continue


def calculate_coeffs_stage_1_sst(data):
    # Latitude and time bins
    bins_lat = np.arange(-72,84,2, dtype=int)
    bins_time = ['day', 'night']
    
    # Loop through bins
    for bin_lat in bins_lat:
        for bin_time in bins_time:
            # Find the data that belong to the latitude bin and time bin
            data_filt = data.loc[(data['lat'] > bin_lat) & (data['lat'] < bin_lat+2) & (data['bin_time'] == bin_time)]
            
            # Construct input data, based on Alerskans 2020, eq.6    
            tb_1v   =  data_filt['1vpol'] - 150
            tb_1v2  = (data_filt['1vpol'] - 150)**2
            tb_1h   =  data_filt['1hpol'] - 150
            tb_1h2  = (data_filt['1hpol'] - 150)**2
            tb_6v   =  data_filt['6vpol'] - 150
            tb_6v2  = (data_filt['6vpol'] - 150)**2
            tb_6h   =  data_filt['6hpol'] - 150
            tb_6h2  = (data_filt['6hpol'] - 150)**2
            tb_10v  =  data_filt['10vpol'] - 150
            tb_10v2 = (data_filt['10vpol'] - 150)**2
            tb_10h  =  data_filt['10hpol'] - 150
            tb_10h2 = (data_filt['10hpol'] - 150)**2
            tb_18v  =  data_filt['18vpol'] - 150
            tb_18v2 = (data_filt['18vpol'] - 150)**2
            tb_18h  =  data_filt['18hpol'] - 150
            tb_18h2 = (data_filt['18hpol'] - 150)**2
            tb_36v  =  data_filt['36vpol'] - 150
            tb_36v2 = (data_filt['36vpol'] - 150)**2
            tb_36h  =  data_filt['36hpol'] - 150
            tb_36h2 = (data_filt['36hpol'] - 150)**2
            theta   = data_filt['look_angle'] - 55
            WSr = data_filt['WSr']
            cos_1phiREL = np.cos(data_filt['phi_rel'])
            sin_1phiREL = np.sin(data_filt['phi_rel'])
            cos_2phiREL = np.cos(2*data_filt['phi_rel'])
            sin_2phiREL = np.sin(2*data_filt['phi_rel'])
            
            X = np.stack((tb_1v,tb_1v2,tb_1h,tb_1h2,tb_6v,tb_6v2,tb_6h,tb_6h2,tb_10v,tb_10v2,tb_10h,tb_10h2,tb_18v,tb_18v2,tb_18h,tb_18h2,tb_36v,tb_36v2,tb_36h,tb_36h2,theta,WSr,cos_1phiREL,sin_1phiREL,cos_2phiREL,sin_2phiREL),axis=1)
            Y = data_filt['sst_ref']
            
            # Calculate linear regression coefficients
            intercept, coeffs = regression(X,Y)
            
            # Save coefficients
            coeffs_all = np.append(intercept,coeffs)
            np.save('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_lat)+'_'+bin_time+'.npy',coeffs_all)


def calculate_coeffs_stage_2(data, SSTa):
    # Attach SSTa to data DataFrame
    data['SSTa'] = SSTa
    # Remove all lines with SSTa = NaN
    data = data[data['SSTa'].notna()]
    # Wind speed bins: 0-20 m/s, 2 m/s bins
    bins_ws = np.arange(0,10,2, dtype=int)
    # SSTa bins where we have data: -2 - 34 C, 2 C bins
    bins_ssta = np.range(-2,36,2, dtype=int)
    
    # Loop through bins
    for bin_ws in bins_ws:
        for bin_ssta in bins_ssta:
            # Find the data that belongs to the wind speed bin and SSTa bin
            data_filt = data.loc[(data['ws'] > bin_ws) & (data['ws'] < bin_ws+2) & (data['SSTa'] > bin_ssta) & (data['SSTa'] < bin_ssta+2)]
            
            # Construct input data, based in Alerskans 2020, eq.10
            tb_1v   =  data_filt['1vpol'] - 150
            tb_1v2  = (data_filt['1vpol'] - 150)**2
            tb_1h   =  data_filt['1hpol'] - 150
            tb_1h2  = (data_filt['1hpol'] - 150)**2
            tb_6v   =  data_filt['6vpol'] - 150
            tb_6v2  = (data_filt['6vpol'] - 150)**2
            tb_6h   =  data_filt['6hpol'] - 150
            tb_6h2  = (data_filt['6hpol'] - 150)**2
            tb_10v  =  data_filt['10vpol'] - 150
            tb_10v2 = (data_filt['10vpol'] - 150)**2
            tb_10h  =  data_filt['10hpol'] - 150
            tb_10h2 = (data_filt['10hpol'] - 150)**2
            tb_18v  =  data_filt['18vpol'] - 150
            tb_18v2 = (data_filt['18vpol'] - 150)**2
            tb_18h  =  data_filt['18hpol'] - 150
            tb_18h2 = (data_filt['18hpol'] - 150)**2
            tb_36v  =  data_filt['36vpol'] - 150
            tb_36v2 = (data_filt['36vpol'] - 150)**2
            tb_36h  =  data_filt['36hpol'] - 150
            tb_36h2 = (data_filt['36hpol'] - 150)**2
            theta   = data_filt['look_angle'] - 55
            WSr = data_filt['WSr']
            cos_1phiREL = np.cos(data_filt['phi_rel'])
            sin_1phiREL = np.sin(data_filt['phi_rel'])
            cos_2phiREL = np.cos(2*data_filt['phi_rel'])
            sin_2phiREL = np.sin(2*data_filt['phi_rel'])
            
            X = np.stack((tb_1v,tb_1v2,tb_1h,tb_1h2,tb_6v,tb_6v2,tb_6h,tb_6h2,tb_10v,tb_10v2,tb_10h,tb_10h2,tb_18v,tb_18v2,tb_18h,tb_18h2,tb_36v,tb_36v2,tb_36h,tb_36h2,theta,WSr,cos_1phiREL,sin_1phiREL,cos_2phiREL,sin_2phiREL),axis=1)
            Y = data_filt['SSTa']
            
            # the try...continue structure is to account for cases where we don't have data in a particular bin, then the code moves on
            try:
                # Calculate linear regression coefficients
                intercept, coeffs = regression(X,Y)
                
                # Save coefficients
                coeffs_all = np.append(intercept,coeffs)
                np.save('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_2/coeffs_ws'+str(bin_ws)+'_temp'+str(int(bin_ssta - 273.15))+'.npy',coeffs_all)
            except:
                continue


def create_mdb_regression_ready():
    path_in = '/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/data_matchup/'
    # Load matchup database
    # ---------------------
    # EMIRAD
    data = pd.read_csv(path_in+'emirad_matchup_data_with_sim_mdb.txt')
    data['time'] = pd.to_datetime(data['time'])
    
    # ISAR
    data_isar = xr.open_dataset('/home/ssk/isar_19_matchup_data_mdb.nc', decode_times=False)
    data['sst'] = data_isar.sea_surface_temperature.data
    dates = data_isar.time.data
    date_ref = []
    for date in dates:
        date_ref.append(dt.datetime.fromtimestamp(date+347151600))
    
    # Calculate relative angles
    data['phi_rel'] = np.asarray(relative_angle(data['azimuth'],data['era5_u10'],data['era5_v10']))
    data['ws'] = np.asarray(np.sqrt(data['era5_u10']**2 + data['era5_v10']**2))
    
    # Add empty column to put time flags
    data['bin_time'] = ''
    
    # Needed to use the Skyfield library
    ts = api.load.timescale()
    eph = api.load('de421.bsp')
    
    # Loop through observations and flag them depending on daytime / nighttime
    for i in np.arange(len(data)):
        ilat = data['lat'][i]
        ilon = data['lon'][i]
        itime = data['time'][i]
        itime = itime.replace(tzinfo=pytz.UTC)
        # The following lines follow the example of https://rhodesmill.org/skyfield/almanac.html#sunrise-and-sunset
        bluffton = api.wgs84.latlon(ilat, ilon)
        t0 = ts.utc(itime.year, itime.month, itime.day, 0, 0, 0)
        t1 = ts.utc(itime.year, itime.month, itime.day, 23, 59, 59)
        t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(eph, bluffton))
        if itime >= pd.to_datetime(t[1].utc_iso()) or itime <= pd.to_datetime(t[0].utc_iso()):
            data['bin_time'][i] = 'night' # 0
        else:
            data['bin_time'][i] = 'day' # 1
    
    data.to_csv('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/data_matchup/emirad_matchup_data_REGRESSION_READY.txt', index=False)
    return data


def plot_data(data,SSTr):
    plt.figure()
    plt.plot_date(data['time'], data['sst'], label='ISAR')
    plt.plot_date(data['time'], SSTr, label='EMIRAD')
    plt.xlabel('Date')
    plt.ylabel('Sea surface temperature (K)')
    plt.legend()
    plt.grid(which='both')
    plt.title('observations vs. retrievals')


def print_stats(data,SSTr):
    dif_mean = np.nanmean(SSTr - data['sst'])
    dif_std = np.nanstd(SSTr - data['sst'])
    dif_rmse = np.sqrt(np.nanmean((SSTr - data['sst'])**2))
    
    print('MW retrievals - IR observations')
    print('-------------------------------')
    print('Mean: ' + str(round(dif_mean,2)))
    print('St.d: ' + str(round(dif_std,2)))
    print('RMSE: ' + str(round(dif_rmse,2)))


# =============================================================================
# LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD 
# =============================================================================
path_coef = '/net/isilon/ifs/arch/home/ea/Toshiba_backup/mmd/mmd06c_post_processed/PMW_regression/PMWR_wspd_sst/Parameters/'


# # Delete all coefficients
# folders = ['/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/',
#            '/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_2/']
# for folder in folders:
#     for filename in os.listdir(folder):
#         file_path = os.path.join(folder, filename)
#         try:
#             if os.path.isfile(file_path) or os.path.islink(file_path):
#                 os.unlink(file_path)
#             elif os.path.isdir(file_path):
#                 shutil.rmtree(file_path)
#         except Exception as e:
#             print('Failed to delete %s. Reason: %s' % (file_path, e))

# Create regression ready database
# --------------------------------
data = create_mdb_regression_ready()

# Load regression ready database
# ------------------------------
data = pd.read_csv('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/data_matchup/emirad_matchup_data_REGRESSION_READY.txt')
data['emirad_c_time'] = pd.to_datetime(data['emirad_c_time'])

# =============================================================================
# MULTIPLE LINEAR REGRESSION - MULTIPLE LINEAR REGRESSION - MULTIPLE LINEAR REG
# =============================================================================
# Wind speed retrieval - stage 1
# ------------------------------
# Calculate coefficients
calculate_coeffs_stage_1_ws(data)

# Load coefficients
A = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws_stage_1.npy')

# Calculate WSa
WSa = A[0]                              + \
    A[1]  *  (data['1vpol'] - 150)      + \
    A[2]  * ((data['1vpol'] - 150)**2)  + \
    A[3]  *  (data['1hpol'] - 150)      + \
    A[4]  * ((data['1hpol'] - 150)**2)  + \
    A[5]  *  (data['6vpol'] - 150)      + \
    A[6]  * ((data['6vpol'] - 150)**2)  + \
    A[7]  *  (data['6hpol'] - 150)      + \
    A[8]  * ((data['6hpol'] - 150)**2)  + \
    A[9]  *  (data['10vpol'] - 150)     + \
    A[10] * ((data['10vpol'] - 150)**2) + \
    A[11] *  (data['10hpol'] - 150)     + \
    A[12] * ((data['10hpol'] - 150)**2) + \
    A[13] *  (data['18vpol'] - 150)     + \
    A[14] * ((data['18vpol'] - 150)**2) + \
    A[15] *  (data['18hpol'] - 150)     + \
    A[16] * ((data['18hpol'] - 150)**2) + \
    A[17] *  (data['36vpol'] - 150      + \
    A[18] * ((data['36vpol'] - 150)**2) + \
    A[19] *  (data['36hpol'] - 150)     + \
    A[20] * ((data['36hpol'] - 150)**2) + \
    A[21] *  (data['look_angle'] - 55)


# Wind speed retrieval - stage 2
# ------------------------------
# Calculate coefficients
calculate_coeffs_stage_2_ws(data, WSa)

# Calculate WSr
WSr = []

# Loop through all observations
for i in np.arange(len(data)):
    # Find the wind speed bin that the observation belongs to
    bin_ws = int(np.floor(data['ws'][i]))
    
    # Load the appropriate coefficient file
    B1 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_2/coeffs_ws'+str(bin_ws)+'.npy')

    # Find the nearest latitude bin of our measurement
    if data['ws'][i] >= bin_ws+0.5:
        bin_ws_near = bin_lat+1
    else:
        bin_ws_near = bin_lat-1 

    # Check if there are coefficients for that bin and load
    if os.path.isfile('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_ws+2)+'.npy'):
        B2 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_ws+2)+'.npy')
    else:
        WSr_i = np.nan
        # Append to list
        WSr.append(WSr_i)
        continue

    WSr_i1 = B1[0]                                         + \
             B1[1]  *  (data['1vpol'] - 150)               + \
             B1[2]  * ((data['1vpol'] - 150)**2)           + \
             B1[3]  *  (data['1hpol'] - 150)               + \
             B1[4]  * ((data['1hpol'] - 150)**2)           + \
             B1[5]  *  (data['6vpol'] - 150)               + \
             B1[6]  * ((data['6vpol'] - 150)**2)           + \
             B1[7]  *  (data['6hpol'] - 150)               + \
             B1[8]  * ((data['6hpol'] - 150)**2)           + \
             B1[9]  *  (data['10vpol'] - 150)              + \
             B1[10] * ((data['10vpol'] - 150)**2)          + \
             B1[11] *  (data['10hpol'] - 150)              + \
             B1[12] * ((data['10hpol'] - 150)**2)          + \
             B1[13] *  (data['18vpol'] - 150)              + \
             B1[14] * ((data['18vpol'] - 150)**2)          + \
             B1[15] *  (data['18hpol'] - 150)              + \
             B1[16] * ((data['18hpol'] - 150)**2)          + \
             B1[17] *  (data['36vpol'] - 150               + \
             B1[18] * ((data['36vpol'] - 150)**2)          + \
             B1[19] *  (data['36hpol'] - 150)              + \
             B1[20] * ((data['36hpol'] - 150)**2)          + \
             B1[21] *  (data['look_angle'] - 55)

    WSr_i2 = B2[0]                                         + \
             B2[1]  *  (data['1vpol'] - 150)               + \
             B2[2]  * ((data['1vpol'] - 150)**2)           + \
             B2[3]  *  (data['1hpol'] - 150)               + \
             B2[4]  * ((data['1hpol'] - 150)**2)           + \
             B2[5]  *  (data['6vpol'] - 150)               + \
             B2[6]  * ((data['6vpol'] - 150)**2)           + \
             B2[7]  *  (data['6hpol'] - 150)               + \
             B2[8]  * ((data['6hpol'] - 150)**2)           + \
             B2[9]  *  (data['10vpol'] - 150)              + \
             B2[10] * ((data['10vpol'] - 150)**2)          + \
             B2[11] *  (data['10hpol'] - 150)              + \
             B2[12] * ((data['10hpol'] - 150)**2)          + \
             B2[13] *  (data['18vpol'] - 150)              + \
             B2[14] * ((data['18vpol'] - 150)**2)          + \
             B2[15] *  (data['18hpol'] - 150)              + \
             B2[16] * ((data['18hpol'] - 150)**2)          + \
             B2[17] *  (data['36vpol'] - 150               + \
             B2[18] * ((data['36vpol'] - 150)**2)          + \
             B2[19] *  (data['36hpol'] - 150)              + \
             B2[20] * ((data['36hpol'] - 150)**2)          + \
             B2[21] *  (data['look_angle'] - 55)

    # Apply linear interpolation to get the WSr value
    f_interp = interp1d([bin_ws+0.5,bin_ws_near+1], [WSr_i1,WSr_i2])
    WSr_i = f_interp(data['ws'][i])
    
    # Append to list
    WSr.append(WSr_i)

# Convert list to numpy array
WSr = np.asarray(WSr)
    
# Add retrieved wind speeds to dataframe
data['WSr'] = WSr


# SST retrieval - stage 1
# -----------------------
# Calculate coefficients
calculate_coeffs_stage_1_sst(data)

# Calculate SSTa
SSTa = []

# Loop through all observations
for i in np.arange(len(data)):
    # Find the latitude bin that our measurement belongs to
    bin_lat = int(np.floor(data['lat'][i] / 2) * 2)  
    
    # Find the time bin the measurement belongs to
    bin_time = data['bin_time'][i]
    
    # Load the appropriate coefficient file
    C1 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_lat)+'_'+bin_time+'.npy')
    
    # Find the nearest latitude bin of our measurement
    if data['lat'][i] > bin_lat+1:
        bin_lat_near = bin_lat+2
    else:
        bin_lat_near = bin_lat-2
        
    # Check if there are coefficients for that bin and load
    if os.path.isfile('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_lat+2)+'_'+bin_time+'.npy'):
        C2 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_lat'+str(bin_lat+2)+'_'+bin_time+'.npy')
    else:
        SSTa_i = np.nan
        # Append to list
        SSTa.append(SSTa_i)
        continue

    SSTa_i1 = C1[0]                                     + \
              C1[1]  *  (data_filt['1vpol'] - 150)      + \
              C1[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              C1[3]  *  (data_filt['1hpol'] - 150)      + \
              C1[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              C1[5]  *  (data_filt['6vpol'] - 150)      + \
              C1[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              C1[7]  *  (data_filt['6hpol'] - 150)      + \
              C1[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              C1[9]  *  (data_filt['10vpol'] - 150)     + \
              C1[10] * ((data_filt['10vpol'] - 150)**2) + \
              C1[11] *  (data_filt['10hpol'] - 150 )    + \
              C1[12] * ((data_filt['10hpol'] - 150)**2) + \
              C1[13] *  (data_filt['18vpol'] - 150)     + \
              C1[14] * ((data_filt['18vpol'] - 150)**2) + \
              C1[15] *  (data_filt['18hpol'] - 150)     + \
              C1[16] * ((data_filt['18hpol'] - 150)**2) + \
              C1[17] *  (data_filt['36vpol'] - 150)     + \
              C1[18] * ((data_filt['36vpol'] - 150)**2) + \
              C1[19] *  (data_filt['36hpol'] - 150)     + \
              C1[20] * ((data_filt['36hpol'] - 150)**2) + \
              C1[21] *  (data_filt['look_angle'] - 55)  + \
              C1[22] * data_filt['WSr']                 + \
              C1[23] * np.cos(data_filt['phi_rel'])     + \
              C1[24] * np.sin(data_filt['phi_rel'])     + \
              C1[25] * np.cos(2*data_filt['phi_rel'])   + \
              C1[26] * np.sin(2*data_filt['phi_rel'])   + \
    
    SSTa_i1 = C2[0]                                     + \
              C2[1]  *  (data_filt['1vpol'] - 150)      + \
              C2[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              C2[3]  *  (data_filt['1hpol'] - 150)      + \
              C2[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              C2[5]  *  (data_filt['6vpol'] - 150)      + \
              C2[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              C2[7]  *  (data_filt['6hpol'] - 150)      + \
              C2[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              C2[9]  *  (data_filt['10vpol'] - 150)     + \
              C2[10] * ((data_filt['10vpol'] - 150)**2) + \
              C2[11] *  (data_filt['10hpol'] - 150 )    + \
              C2[12] * ((data_filt['10hpol'] - 150)**2) + \
              C2[13] *  (data_filt['18vpol'] - 150)     + \
              C2[14] * ((data_filt['18vpol'] - 150)**2) + \
              C2[15] *  (data_filt['18hpol'] - 150)     + \
              C2[16] * ((data_filt['18hpol'] - 150)**2) + \
              C2[17] *  (data_filt['36vpol'] - 150)     + \
              C2[18] * ((data_filt['36vpol'] - 150)**2) + \
              C2[19] *  (data_filt['36hpol'] - 150)     + \
              C2[20] * ((data_filt['36hpol'] - 150)**2) + \
              C2[21] *  (data_filt['look_angle'] - 55)  + \
              C2[22] * data_filt['WSr']                 + \
              C2[23] * np.cos(data_filt['phi_rel'])     + \
              C2[24] * np.sin(data_filt['phi_rel'])     + \
              C2[25] * np.cos(2*data_filt['phi_rel'])   + \
              C2[26] * np.sin(2*data_filt['phi_rel'])   + \

    # Apply linear interpolation to get the SSTa value
    f_interp = interp1d([bin_lat+1,bin_lat_near+1], [SSTa_i1,SSTa_i2])
    SSTa_i = f_interp(data['lat'][i])
    
    # Append to list
    SSTa.append(SSTa_i)

# Convert list to numpy array
SSTa = np.asarray(SSTa)


# SST retrieval - stage 2
# -----------------------
# # Calculate coefficients
calculate_coeffs_stage_2(data,SSTa)

# Calculate SSTr
SSTr = []

# Loop through all observations
for i in np.arange(len(data)):
    
    if np.isnan(SSTa[i]):
        SSTr_i = np.nan
        SSTr.append(SSTr_i)
        continue
    
    # Find the wind speed bin that our measurement belongs to
    bin_ws = int(np.floor(data['ws'][i] / 2) * 2)
    
    # Find the SSTa bin the measurement belongs to
    bin_ssta = int(np.floor((SSTa[i] - 273.15) / 2) * 2)
    
    # Load the appropriate coefficient file
    D1 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_2/coeffs_ws'+str(bin_ws)+'_temp'+str(bin_ssta)+'.npy')

    # Find the nearest wind speed bin
    if data['ws'][i] > bin_ws+1:
        bin_ws_near = bin_ws+2
    else:
        bin_ws_near = bin_ws-2 
    
    # Find the nearest SSTa bin
    if SSTa[i] > bin_ssta+1:
        bin_ssta_near = bin_ssta+2
    else:
        bin_ssta_near = bin_ssta-2 
    
    # Check if there are coefficients for that bin and load
    if os.path.isfile('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws+2)+'_temp'+str(bin_ssta)+'.npy'):
        D2 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws+2)+'_temp'+str(bin_ssta)+'.npy')
    else:
        SSTr_i = np.nan
        # Append to list
        SSTr.append(SSTr_i)
        continue
    
    # Check if there are coefficients for that bin and load
    if os.path.isfile('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws)+'_temp'+str(bin_ssta+2)+'.npy'):
        D3 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws)+'_temp'+str(bin_ssta+2)+'.npy')
    else:
        SSTr_i = np.nan
        # Append to list
        SSTr.append(SSTr_i)
        continue
    
    # Check if there are coefficients for that bin and load
    if os.path.isfile('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws+2)+'_temp'+str(bin_ssta+2)+'.npy'):
        D4 = np.load('/data/users/ssk/SHIPS4SST/phase_2/mw_ir_campaign_2021/data/retrieval_coeffs/stage_1/coeffs_ws'+str(bin_ws+2)+'_temp'+str(bin_ssta+2)+'.npy')
    else:
        SSTr_i = np.nan
        # Append to list
        SSTr.append(SSTr_i)
        continue
    
    if TB_data_type == 'obs':   
        SSTr_i1 = D1[0]                                     + \
              D1[1]  *  (data_filt['1vpol'] - 150)      + \
              D1[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              D1[3]  *  (data_filt['1hpol'] - 150)      + \
              D1[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              D1[5]  *  (data_filt['6vpol'] - 150)      + \
              D1[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              D1[7]  *  (data_filt['6hpol'] - 150)      + \
              D1[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              D1[9]  *  (data_filt['10vpol'] - 150)     + \
              D1[10] * ((data_filt['10vpol'] - 150)**2) + \
              D1[11] *  (data_filt['10hpol'] - 150 )    + \
              D1[12] * ((data_filt['10hpol'] - 150)**2) + \
              D1[13] *  (data_filt['18vpol'] - 150)     + \
              D1[14] * ((data_filt['18vpol'] - 150)**2) + \
              D1[15] *  (data_filt['18hpol'] - 150)     + \
              D1[16] * ((data_filt['18hpol'] - 150)**2) + \
              D1[17] *  (data_filt['36vpol'] - 150)     + \
              D1[18] * ((data_filt['36vpol'] - 150)**2) + \
              D1[19] *  (data_filt['36hpol'] - 150)     + \
              D1[20] * ((data_filt['36hpol'] - 150)**2) + \
              D1[21] *  (data_filt['look_angle'] - 55)  + \
              D1[22] * data_filt['WSr']                 + \
              D1[23] * np.cos(data_filt['phi_rel'])     + \
              D1[24] * np.sin(data_filt['phi_rel'])     + \
              D1[25] * np.cos(2*data_filt['phi_rel'])   + \
              D1[26] * np.sin(2*data_filt['phi_rel'])   + \
                  
        SSTr_i2 = D2[0]                                     + \
              D2[1]  *  (data_filt['1vpol'] - 150)      + \
              D2[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              D2[3]  *  (data_filt['1hpol'] - 150)      + \
              D2[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              D2[5]  *  (data_filt['6vpol'] - 150)      + \
              D2[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              D2[7]  *  (data_filt['6hpol'] - 150)      + \
              D2[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              D2[9]  *  (data_filt['10vpol'] - 150)     + \
              D2[10] * ((data_filt['10vpol'] - 150)**2) + \
              D2[11] *  (data_filt['10hpol'] - 150 )    + \
              D2[12] * ((data_filt['10hpol'] - 150)**2) + \
              D2[13] *  (data_filt['18vpol'] - 150)     + \
              D2[14] * ((data_filt['18vpol'] - 150)**2) + \
              D2[15] *  (data_filt['18hpol'] - 150)     + \
              D2[16] * ((data_filt['18hpol'] - 150)**2) + \
              D2[17] *  (data_filt['36vpol'] - 150)     + \
              D2[18] * ((data_filt['36vpol'] - 150)**2) + \
              D2[19] *  (data_filt['36hpol'] - 150)     + \
              D2[20] * ((data_filt['36hpol'] - 150)**2) + \
              D2[21] *  (data_filt['look_angle'] - 55)  + \
              D2[22] * data_filt['WSr']                 + \
              D2[23] * np.cos(data_filt['phi_rel'])     + \
              D2[24] * np.sin(data_filt['phi_rel'])     + \
              D2[25] * np.cos(2*data_filt['phi_rel'])   + \
              D2[26] * np.sin(2*data_filt['phi_rel'])   + \
                  
        SSTr_i3 = D3[0]                                     + \
              D3[1]  *  (data_filt['1vpol'] - 150)      + \
              D3[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              D3[3]  *  (data_filt['1hpol'] - 150)      + \
              D3[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              D3[5]  *  (data_filt['6vpol'] - 150)      + \
              D3[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              D3[7]  *  (data_filt['6hpol'] - 150)      + \
              D3[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              D3[9]  *  (data_filt['10vpol'] - 150)     + \
              D3[10] * ((data_filt['10vpol'] - 150)**2) + \
              D3[11] *  (data_filt['10hpol'] - 150 )    + \
              D3[12] * ((data_filt['10hpol'] - 150)**2) + \
              D3[13] *  (data_filt['18vpol'] - 150)     + \
              D3[14] * ((data_filt['18vpol'] - 150)**2) + \
              D3[15] *  (data_filt['18hpol'] - 150)     + \
              D3[16] * ((data_filt['18hpol'] - 150)**2) + \
              D3[17] *  (data_filt['36vpol'] - 150)     + \
              D3[18] * ((data_filt['36vpol'] - 150)**2) + \
              D3[19] *  (data_filt['36hpol'] - 150)     + \
              D3[20] * ((data_filt['36hpol'] - 150)**2) + \
              D3[21] *  (data_filt['look_angle'] - 55)  + \
              D3[22] * data_filt['WSr']                 + \
              D3[23] * np.cos(data_filt['phi_rel'])     + \
              D3[24] * np.sin(data_filt['phi_rel'])     + \
              D3[25] * np.cos(2*data_filt['phi_rel'])   + \
              D3[26] * np.sin(2*data_filt['phi_rel'])   + \
                  
        SSTr_i4 = D4[0]                                     + \
              D4[1]  *  (data_filt['1vpol'] - 150)      + \
              D4[2]  * ((data_filt['1vpol'] - 150)**2)  + \
              D4[3]  *  (data_filt['1hpol'] - 150)      + \
              D4[4]  * ((data_filt['1hpol'] - 150)**2)  + \
              D4[5]  *  (data_filt['6vpol'] - 150)      + \
              D4[6]  * ((data_filt['6vpol'] - 150)**2)  + \
              D4[7]  *  (data_filt['6hpol'] - 150)      + \
              D4[8]  * ((data_filt['6hpol'] - 150)**2)  + \
              D4[9]  *  (data_filt['10vpol'] - 150)     + \
              D4[10] * ((data_filt['10vpol'] - 150)**2) + \
              D4[11] *  (data_filt['10hpol'] - 150 )    + \
              D4[12] * ((data_filt['10hpol'] - 150)**2) + \
              D4[13] *  (data_filt['18vpol'] - 150)     + \
              D4[14] * ((data_filt['18vpol'] - 150)**2) + \
              D4[15] *  (data_filt['18hpol'] - 150)     + \
              D4[16] * ((data_filt['18hpol'] - 150)**2) + \
              D4[17] *  (data_filt['36vpol'] - 150)     + \
              D4[18] * ((data_filt['36vpol'] - 150)**2) + \
              D4[19] *  (data_filt['36hpol'] - 150)     + \
              D4[20] * ((data_filt['36hpol'] - 150)**2) + \
              D4[21] *  (data_filt['look_angle'] - 55)  + \
              D4[22] * data_filt['WSr']                 + \
              D4[23] * np.cos(data_filt['phi_rel'])     + \
              D4[24] * np.sin(data_filt['phi_rel'])     + \
              D4[25] * np.cos(2*data_filt['phi_rel'])   + \
              D4[26] * np.sin(2*data_filt['phi_rel'])   + \
    
    # Bilinear interpolation
    f_interp2 = interp2d([bin_ws+1,bin_ws_near+1,bin_ws+1,bin_ws_near+1], [bin_ssta+1,bin_ssta+1,bin_ssta_near+1,bin_ssta_near+1], [SSTr_i1,SSTr_i2,SSTr_i3,SSTr_i4])
    SSTr_i = f_interp2(data['ws'][i],SSTa[i])
    
    # Append to list
    SSTr.append(SSTr_i)

# Convert list to numpy array
SSTr = np.asarray(SSTr)


# Plot mw retrievals together with ISAR observations
plot_data(data,SSTr)

# Print statistics
print_stats(data,SSTr)






