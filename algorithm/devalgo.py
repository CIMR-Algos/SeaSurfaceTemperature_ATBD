#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created  on Fri Apr  8 2022 - @author: ssk
Updated  on Fri Apr 22 2022 - @author: ssk
Modified on Tor Jul 28 2022 - @author: ea

"""
# Library imports
import datetime
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sklearn
from sklearn.model_selection import train_test_split 
from sklearn import linear_model
import sys

# Set the random seed
rnseed = 42

# Set project paths
HOME = os.environ['HOME']
PROJECT = HOME + "/Projects/CIMR/DEVALGO/ATBD_SST"
DATA_PATH = PROJECT + "/data"
COEFF_PATH = PROJECT + "/coeffs"

# Settings
Nlim = 100

# Citation: see https://scikit-learn.org/stable/about.html#citing-scikit-learn


# =========
# FUNCTIONS
# =========
def get_vars2extract(data_type):
    if data_type == "mmd_all":
        var_names = ['orbit', 'lat', 'lon', 'solza', 'satza', 'solaz', 'sataz', 'era5_wind_dir', 'ccmp_wind_dir', \
                     'era5_phi_rel', 'era5_ws', 'ccmp_ws', 'era5_sst', 'era5_tcwv', 'era5_clwt', 'tb6vpol', 'tb6hpol', \
                     'tb10vpol', 'tb10hpol', 'tb18vpol', 'tb18hpol', 'tb23vpol', 'tb23Hpol', 'tb36vpol', 'tb36hpol', \
                     'tb89vpol', 'tb89hpol', 'sga', 'sss', 'insitu_sst', 'insitu_time']
    elif data_type == "mmd":
        var_names = ['orbit', 'lat', 'lon', 'satza', 'sataz', 'era5_wind_dir', 'era5_phi_rel', 'era5_ws', 'era5_sst', \
                     'era5_tcwv', 'era5_clwt', 'tb6vpol', 'tb6hpol', 'tb10vpol', 'tb10hpol', 'tb18vpol', 'tb18hpol', \
                     'tb23vpol', 'tb23hpol', 'tb36vpol', 'tb36hpol', 'tb89vpol', 'tb89hpol', 'sga', 'sss', 'insitu_sst', \
                     'insitu_time']
    else:
        print("Could not parse data type: {}".format(data_type))
        print("Exiting...!")
        sys.exit()

    return var_names



def relative_angle(azimuth, u, v):
    D2R = np.pi/180
    #phi_w = 90 - (np.arctan2(u, v))*D2R  # wind to direction
    phiw = 180 + D2R*np.arctan2(u,v)      # wind from direction
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



def ws_algorithm_selection(data,sensor):
    nmatchups = data.shape[0]
    if sensor == 'amsre':
        # Construct input data for AMSR-E configuration
        ninput = 21
        X = np.full((nmatchups,ninput), fill_value=np.nan)

        X[:,0]  =  data['tb6vpol'] - 150
        X[:,1]  = (data['tb6vpol'] - 150)**2
        X[:,2]  =  data['tb6hpol'] - 150
        X[:,3]  = (data['tb6hpol'] - 150)**2
        X[:,4]  =  data['tb10vpol'] - 150
        X[:,5]  = (data['tb10vpol'] - 150)**2
        X[:,6]  =  data['tb10hpol'] - 150
        X[:,7]  = (data['tb10hpol'] - 150)**2
        X[:,8]  =  data['tb18vpol'] - 150
        X[:,9]  = (data['tb18vpol'] - 150)**2
        X[:,10] =  data['tb18hpol'] - 150
        X[:,11] = (data['tb18hpol'] - 150)**2
        X[:,12] =  data['tb23vpol'] - 150
        X[:,13] = (data['tb23vpol'] - 150)**2
        X[:,14] =  data['tb23hpol'] - 150
        X[:,15] = (data['tb23hpol'] - 150)**2
        X[:,16] =  data['tb36vpol'] - 150
        X[:,17] = (data['tb36vpol'] - 150)**2
        X[:,18] =  data['tb36hpol'] - 150
        X[:,19] = (data['tb36hpol'] - 150)**2
        X[:,20] = data['look_angle']
    elif sensor == 'cimr':
        # Construct input data for CIMR configuration
        ninput = 21
        X = np.full((nmatchups,ninput), fill_value=np.nan)

        X[:,0]  =  data['tb1vpol'] - 150
        X[:,1]  = (data['tb1vpol'] - 150)**2
        X[:,2]  =  data['tb1hpol'] - 150
        X[:,3]  = (data['tb1hpol'] - 150)**2
        X[:,4]  =  data['tb6vpol'] - 150
        X[:,5]  = (data['tb6vpol'] - 150)**2
        X[:,6]  =  data['tb6hpol'] - 150
        X[:,7]  = (data['tb6hpol'] - 150)**2
        X[:,8]  =  data['tb10vpol'] - 150
        X[:,9]  = (data['tb10vpol'] - 150)**2
        X[:,10] =  data['tb10hpol'] - 150
        X[:,11] = (data['tb10hpol'] - 150)**2
        X[:,12] =  data['tb18vpol'] - 150
        X[:,13] = (data['tb18vpol'] - 150)**2
        X[:,14] =  data['tb18hpol'] - 150
        X[:,15] = (data['tb18hpol'] - 150)**2
        X[:,16] =  data['tb36vpol'] - 150
        X[:,17] = (data['tb36vpol'] - 150)**2
        X[:,18] =  data['tb36hpol'] - 150
        X[:,19] = (data['tb36hpol'] - 150)**2
        X[:,20] = data['look_angle']
    else:
        print("Sensor {} not understood.\nExiting...!".format(sensor))
        sys.exit()

    return X



def sst_algorithm_selection(data,sensor):
    nmatchups = data.shape[0]
    if sensor == 'amsre':
        # Construct input data for AMSR-E configuration
        ninput = 26
        X = np.full((nmatchups,ninput), fill_value=np.nan)

        X[:,0]  =  data['tb6vpol'] - 150
        X[:,1]  = (data['tb6vpol'] - 150)**2
        X[:,2]  =  data['tb6hpol'] - 150
        X[:,3]  = (data['tb6hpol'] - 150)**2
        X[:,4]  =  data['tb10vpol'] - 150
        X[:,5]  = (data['tb10vpol'] - 150)**2
        X[:,6]  =  data['tb10hpol'] - 150
        X[:,7]  = (data['tb10hpol'] - 150)**2
        X[:,8]  =  data['tb18vpol'] - 150
        X[:,9]  = (data['tb18vpol'] - 150)**2
        X[:,10] =  data['tb18hpol'] - 150
        X[:,11] = (data['tb18hpol'] - 150)**2
        X[:,12] =  data['tb23vpol'] - 150
        X[:,13] = (data['tb23vpol'] - 150)**2
        X[:,14] =  data['tb23hpol'] - 150
        X[:,15] = (data['tb23hpol'] - 150)**2
        X[:,16] =  data['tb36vpol'] - 150
        X[:,17] = (data['tb36vpol'] - 150)**2
        X[:,18] =  data['tb36hpol'] - 150
        X[:,19] = (data['tb36hpol'] - 150)**2
        X[:,20] = data['look_angle']
        X[:,21] = data['WSr']
        X[:,22] = np.cos(data['era5_phi_rel'])
        X[:,23] = np.sin(data['era5_phi_rel'])
        X[:,24] = np.cos(2*data['era5_phi_rel'])
        X[:,25] = np.sin(2*data['era5_phi_rel'])
    elif sensor == 'cimr':
        # Construct input data for AMSR-E configuration
        ninput = 26
        X = np.full((nmatchups,ninput), fill_value=np.nan)

        X[:,0]  =  data['tb1vpol'] - 150
        X[:,1]  = (data['tb1vpol'] - 150)**2
        X[:,2]  =  data['tb1hpol'] - 150
        X[:,3]  = (data['tb1hpol'] - 150)**2
        X[:,4]  =  data['tb6vpol'] - 150
        X[:,5]  = (data['tb6vpol'] - 150)**2
        X[:,6]  =  data['tb6hpol'] - 150
        X[:,7]  = (data['tb6hpol'] - 150)**2
        X[:,8]  =  data['tb10vpol'] - 150
        X[:,9]  = (data['tb10vpol'] - 150)**2
        X[:,10] =  data['tb10hpol'] - 150
        X[:,11] = (data['tb10hpol'] - 150)**2
        X[:,12] =  data['tb18vpol'] - 150
        X[:,13] = (data['tb18vpol'] - 150)**2
        X[:,14] =  data['tb18hpol'] - 150
        X[:,15] = (data['tb18hpol'] - 150)**2
        X[:,16] =  data['tb36vpol'] - 150
        X[:,17] = (data['tb36vpol'] - 150)**2
        X[:,18] =  data['tb36hpol'] - 150
        X[:,19] = (data['tb36hpol'] - 150)**2
        X[:,20] = data['look_angle']
        X[:,21] = data['WSr']
        X[:,22] = np.cos(data['era5_phi_rel'])
        X[:,23] = np.sin(data['era5_phi_rel'])
        X[:,24] = np.cos(2*data['era5_phi_rel'])
        X[:,25] = np.sin(2*data['era5_phi_rel'])
    else:
        print("Sensor {} not understood.\nExiting...!".format(sensor))
        sys.exit()

    return X



def calculate_coeffs_ws(data,sensor):
    # Predictors
    X = ws_algorithm_selection(data,sensor)
    # Predictand
    Y = data['era5_ws'].values
    
    # the try...continue structure is to account for cases where we don't have data in a particular bin, then the code moves on
    try:
        # Calculate linear regression coefficients
        intercept, coeffs = regression(X,Y)
        coeffs_all = np.append(intercept,coeffs)

        return coeffs_all
    except:
        return np.full((ninput+1), fill_value=np.nan)



def calculate_coeffs_sst(data,sensor):
    # Predictors
    X = sst_algorithm_selection(data,sensor)
    # Predictand
    Y = data['insitu_sst'].values
    
    # the try...continue structure is to account for cases where we don't have data in a particular bin, then the code moves on
    try:
        # Calculate linear regression coefficients
        intercept, coeffs = regression(X,Y)
        coeffs_all = np.append(intercept,coeffs)

        return coeffs_all
    except:
        return np.full((ninput+1), fill_value=np.nan)



def retrieve_ws(data,A,sensor):
    # Settings
    nmatchups = data.shape[0]
    # Predictors
    X = ws_algorithm_selection(data,sensor)
    # Add input for the offset/intercept
    X = np.hstack((np.ones((nmatchups,1)),X))
    # Retrieve WS
    WS =  np.dot(A,X.T)

    return WS



def retrieve_sst(data,A,sensor):
    # Settings
    nmatchups = data.shape[0]
    # Predictors
    X = sst_algorithm_selection(data,sensor)
    # Add input for the offset/intercept
    X = np.hstack((np.ones((nmatchups,1)),X))
    # Retrieve SST
    SST =  np.dot(A,X.T)

    return SST



def calculate_coeffs_stage_1_ws(data,sensor,verbose=False):
    print("\nCalculate WS stage 1 coefficients")

    # Get coefficients
    coeffs_all = calculate_coeffs_ws(data,sensor)

    if not np.all(np.isnan(coeffs_all)):
        if verbose:
            print("Save WS stage 1 coefficients")
        coeffs_file = COEFF_PATH + "/ws/coeffs_ws_stage_1.npy"
        np.save(coeffs_file,coeffs_all)



def retrieve_stage_1_ws(data,sensor):
    print("\nRetrieve stage 1 WS")

    # Load coefficients
    coeffs_file = COEFF_PATH + "/ws/coeffs_ws_stage_1.npy"
    A = np.load(coeffs_file)
    
    # Retrieve WS
    WSa =  retrieve_ws(data,A,sensor)

    return WSa



def calculate_coeffs_stage_2_ws(data,ws_bins,sensor,verbose=False):
    print("\nCalculate WS stage 2 coefficients")

    # Settings
    ws_step = ws_bins[1] - ws_bins[0]
    
    # Loop through bins
    for iws in ws_bins:
        # Find the data that belongs to the current wind speed bin 
        mask_sub = ( (data['WSa'].values > iws) & (data['WSa'].values <= iws+ws_step) )

        # Check so that there is enough data
        if np.sum(mask_sub) > Nlim:
            data_sub = data.loc[mask_sub]
                
            # Calculate coefficients
            coeffs_all = calculate_coeffs_ws(data_sub,sensor)
    
            if not np.all(np.isnan(coeffs_all)):
                if verbose:
                    print("Save WS stage 2 coefficients for wind speed bin {}-{} ms-1".format(iws,iws+1))
                coeffs_file = COEFF_PATH + "/ws/coeffs_ws_stage_2_wsbin_"+str(iws)+".npy"
                np.save(coeffs_file,coeffs_all)



def retrieve_stage_2_ws(data,ws_bins,sensor,verbose=False):
    print("\nRetrieve stage 2 WS")

    # Settings
    nmatchups = data.shape[0]
    ws_step = ws_bins[1] - ws_bins[0]

    # Initialize array
    WSr = np.full((nmatchups), fill_value=np.nan, dtype=np.float32)

    # Loop through bins
    for iws in ws_bins:
        if verbose:
            print("Working on wind speed bin {}".format(iws))
        # Find the data that belongs to the current wind speed bin 
        mask_sub = (data['WSa'].values > iws) & (data['WSa'].values <= iws+ws_step)

        # Check so that there is data
        if np.sum(mask_sub) > 0:
            idx_sub = np.argwhere(mask_sub)[:,0]
            data_sub = data.loc[mask_sub]
            data_sub.reset_index(inplace=True,drop=True)
    
            # Load the appropriate coefficient file
            coeffs_file = COEFF_PATH + "/ws/coeffs_ws_stage_2_wsbin_"+str(iws)+".npy"
            if os.path.isfile(coeffs_file):
                B1 = np.load(coeffs_file)
                isnan_B1 = False
            else:
                isnan_B1 = True
                
            # If interpolating, find the nearest wind speed bin of our measurement
            for inear in range(2):
                if inear == 0:
                    # Lower limit
                    iws_near = iws - ws_step
                    mask_int = (data_sub['WSa'].values < iws + ws_step/2)
                    if np.sum(mask_int) == 0:
                        continue
                    idx_int = np.argwhere(mask_int)[:,0]
                    data_int = data_sub.loc[mask_int]
                else:
                    # Upper limit
                    iws_near = iws + ws_step
                    mask_int = (data_sub['WSa'].values >= iws + ws_step/2)
                    if np.sum(mask_int) == 0:
                        continue
                    idx_int = np.argwhere(mask_int)[:,0]
                    data_int = data_sub.loc[mask_int]
    
                # Check if there are coefficients for that bin and load
                coeffs_file_near = COEFF_PATH + "/ws/coeffs_ws_stage_2_wsbin_"+str(iws_near)+".npy"
                if os.path.isfile(coeffs_file_near):
                    B2 = np.load(coeffs_file_near)
                    isnan_B2 = False
                elif os.path.isfile(coeffs_file):
                    isnan_B2 = True
                    B2 = B1.copy()
                else:
                    # No coefficients for current or closest ws bin -> no retrieval
                    if verbose:
                        print("Warning: Coefficients do not exist for current wind speed bin {} or closest bin. WS=NaN...!".format(iws))
                    WSr[idx_sub[idx_int]] = np.nan
                    continue

                # Retrieve WS
                WSr_i1 = retrieve_ws(data_int,B1,sensor)
                WSr_i2 = retrieve_ws(data_int,B2,sensor)
        
                # Define interpolation weights
                w1 = np.abs(data_int['WSa'].values - iws)/ws_step
                w2 = 1 - w1

                # Reset weights as we don't have either B1 or B2 coefficients
                if (isnan_B1 | isnan_B2):
                    if (isnan_B1): w1 = 0.
                    if (isnan_B2): w2 = 0.

                    wsum = w1 + w2
                    w1 = w1 / wsum
                    w2 = w2 / wsum

                # Interpolate between WSr_i1 and WSr_i2
                WSr_int = WSr_i1 * w1 + WSr_i2 * w2
    
                # Assign to the correct WSr elements
                WSr[idx_sub[idx_int]] = WSr_int.copy()

    return WSr



def calculate_coeffs_stage_1_sst(data,orb_bins,lat_bins,sensor,verbose=False):
    print("\nCalculate SST stage 1 coefficients")

    # Settings
    lat_step = lat_bins[1] - lat_bins[0]
    
    # Loop through bins
    for iorb in orb_bins:
        for ilat in lat_bins:
            if verbose:
                print("Working on orbit bin {} and latitude bin {}.".format(iorb,ilat))
            # Find the data that belongs to the current wind speed bin - remove nan values
            mask_sub = ( (data['lat'].values > ilat) & (data['lat'].values <= ilat+lat_step) & (data['orbit'].values == iorb) & \
                         ~np.isnan(data['WSr'].values) )

            # Check so that there is enough data
            if np.sum(mask_sub) > Nlim:
                data_sub = data.loc[mask_sub]
    
                # Calculate coefficients
                coeffs_all = calculate_coeffs_sst(data_sub,sensor)
    
                if not np.all(np.isnan(coeffs_all)):
                    if verbose:
                        print("Save SST stage 1 coefficients for obrit bin {} and latitude bin {}-{} deg".format(iorb,ilat,ilat+1))
                    coeffs_file = COEFF_PATH + "/sst/coeffs_sst_stage_1_orbbin_"+str(iorb)+"_latbin_"+str(ilat)+".npy"
                    np.save(coeffs_file,coeffs_all)



def retrieve_stage_1_sst(data,orb_bins,lat_bins,sensor,verbose=False):
    print("\nRetrieve stage 1 SST")

    # Settings
    nmatchups = data.shape[0]
    orb_step = orb_bins[1] - orb_bins[0]
    lat_step = lat_bins[1] - lat_bins[0]

    # Initialize array
    SSTa = np.full((nmatchups), fill_value=np.nan, dtype=np.float32)

    # Loop through bins
    for iorb in orb_bins:
        for ilat in lat_bins:
            if verbose:
                print("Working on orbit bin {} and latitude bin {}.".format(iorb,ilat))
            # Find the data that belongs to the current orbit and latitude bin 
            mask_sub = ( (data['lat'].values > ilat) & (data['lat'].values <= ilat+lat_step) & (data['orbit'].values == iorb) )

            # Check so that there is data
            if np.sum(mask_sub) > 0:
                idx_sub = np.argwhere( mask_sub )[:,0]
                data_sub = data.loc[ mask_sub ]
                data_sub.reset_index(inplace=True,drop=True)
    
                # Load the appropriate coefficient file
                coeffs_file = COEFF_PATH + "/sst/coeffs_sst_stage_1_orbbin_"+str(iorb)+"_latbin_"+str(ilat)+".npy"
                if os.path.isfile(coeffs_file):
                    C1 = np.load(coeffs_file)
                    isnan_C1 = False
                else:
                    isnan_C1 = True
                    
                # If interpolating (only latitude), find the nearest latitude bin of our measurement
                for inear in range(2):
                    if inear == 0:
                        # Lower limit
                        ilat_near = ilat - lat_step
                        mask_int = (data_sub['lat'].values <  ilat + lat_step/2)
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
                    else:
                        # Upper limit
                        ilat_near = ilat + lat_step
                        mask_int = (data_sub['lat'].values >=  ilat + lat_step/2)
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
        
                    # Check if there are coefficients for that bin and load
                    coeffs_file_near = COEFF_PATH + "/sst/coeffs_sst_stage_1_orbbin_"+str(iorb)+"_latbin_"+str(ilat_near)+".npy"
                    if os.path.isfile(coeffs_file_near):
                        C2 = np.load(coeffs_file_near)
                        isnan_C2 = False
                    elif os.path.isfile(coeffs_file):
                        C2 = C1.copy()
                        isnan_C2 = True
                    else:
                        # No coefficients for current or closest ws bin -> no retrieval
                        if verbose:
                            print("Warning: Coefficients do not exist for current orbit {} and latitude bin {} or closest bin. SST=NaN...!".format(iorb,ilat))
                        SSTa[idx_sub[idx_int]] = np.nan
                        continue

                    # Retrieve SST
                    SSTa_i1 = retrieve_sst(data_int,C1,sensor)
                    SSTa_i2 = retrieve_sst(data_int,C2,sensor)
            
                    # Define interpolation weights
                    w1 = np.abs(data_int['lat'].values - ilat)/lat_step
                    w2 = 1 - w1

                    # Reset weights as we don't have either C1 or C2 coefficients
                    if (isnan_C1 | isnan_C2):
                        if (isnan_C1): w1 = 0.
                        if (isnan_C2): w2 = 0.

                        wsum = w1 + w2
                        w1 = w1 / wsum
                        w2 = w2 / wsum

                    # Interpolate between SSTa_i1 and SSTa_i2
                    SSTa_int = SSTa_i1 * w1 + SSTa_i2 * w2
        
                    # Assign to the correct SSTa elements
                    SSTa[idx_sub[idx_int]] = SSTa_int.copy()

    return SSTa



def calculate_coeffs_stage_2_sst(data,ws_bins,sst_bins,sensor,verbose=False):
    print("\nCalculate SST stage 2 coefficients")

    # Settings
    ws_step = ws_bins[1] - ws_bins[0]
    sst_step = sst_bins[1] - sst_bins[0]
    
    # Loop through bins
    for iws in ws_bins:
        for isst in sst_bins:
            if verbose:
                print("Working on wind speed bin {} and sst bin {}.".format(iws,isst))
            # Find the data that belongs to the current wind speed bin - remove nan values
            mask_sub = ( (data['WSr'].values > iws) & (data['WSr'].values <= iws+ws_step) & \
                         (data['SSTa'].values > (isst+273.15)) & (data['SSTa'].values <= (isst+sst_step+273.15)) & \
                         ~np.isnan(data['WSr'].values) & ~np.isnan(data['SSTa'].values) )

            # Check so that there is enough data
            if np.sum(mask_sub) > Nlim:
                data_sub = data.loc[ mask_sub ]
    
                # Calculate coefficients
                coeffs_all = calculate_coeffs_sst(data_sub,sensor)
    
                if not np.all(np.isnan(coeffs_all)):
                    if verbose:
                        print("Save SST stage 2 coefficients for wind speed bin {}-{} ms-1 and sst bin {}-{} degC".format(iws,iws+1,isst,isst+1))
                    coeffs_file = COEFF_PATH + "/sst/coeffs_sst_stage_2_wsbin_"+str(iws)+"_sstbin_"+str(isst)+".npy"
                    np.save(coeffs_file,coeffs_all)



def retrieve_stage_2_sst(data,ws_bins,sst_bins,sensor,verbose=False):
    print("\nRetrieve stage 2 SST")

    # Settings
    nmatchups = data.shape[0]
    ws_step = ws_bins[1] - ws_bins[0]
    sst_step = sst_bins[1] - sst_bins[0]

    # Number of input to the retrieval algorithm
    nsst_input = 27

    # Initialize array
    SSTr = np.full((nmatchups), fill_value=np.nan, dtype=np.float32)

    # Loop through bins
    for iws in ws_bins:
        for isst in sst_bins:
            if verbose:
                print("Working on wind speed bin {} and sst bin {}.".format(iws,isst))
            # Find the data that belongs to the current orbit and latitude bin 
            mask_sub = ( (data['WSr'].values > iws) & (data['WSr'].values <= iws+ws_step) & \
                         (data['SSTa'].values > (isst+273.15)) & (data['SSTa'].values <= (isst+sst_step+273.15)) )

            # Check so that there is data
            if np.sum(mask_sub) > 0:
                idx_sub = np.argwhere( mask_sub )[:,0]
                data_sub = data.loc[ mask_sub ]
                data_sub.reset_index(inplace=True,drop=True)
    
                # Load the appropriate coefficient file
                coeffs_file = COEFF_PATH + "/sst/coeffs_sst_stage_2_wsbin_"+str(iws)+"_sstbin_"+str(isst)+".npy"
                if os.path.isfile(coeffs_file):
                    D11 = np.load(coeffs_file)
                    isnan_D11 = False
                else:
                    D11 = np.full((nsst_input), fill_value=np.nan)
                    isnan_D11 = True

                # If interpolating (both ws and sst), find the nearest ws and sst bins of our measurement
                for inear in range(4):
                    if inear == 0:
                        # iws-1 & isst-1 limit
                        iws_near = iws - ws_step
                        isst_near = isst - sst_step
                        mask_int = ( (data_sub['WSr'].values < iws + ws_step/2) & (data_sub['SSTa'].values < (isst + sst_step/2 + 273.15)) )
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
                    elif inear == 1:
                        # iws+1 & isst-1 limit
                        iws_near = iws + ws_step
                        isst_near = isst - sst_step
                        mask_int = ( (data_sub['WSr'].values >= iws + ws_step/2) & (data_sub['SSTa'].values < (isst + sst_step/2 + 273.15)) )
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
                    elif inear == 2:
                        # iws-1 & isst+1 limit
                        iws_near = iws - ws_step
                        isst_near = isst + sst_step
                        mask_int = ( (data_sub['WSr'].values < iws + ws_step/2) & (data_sub['SSTa'].values >= (isst + sst_step/2 + 273.15)) )
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
                    elif inear == 3:
                        # iws+1 & isst+1 limit
                        iws_near = iws + ws_step
                        isst_near = isst + sst_step
                        mask_int = ( (data_sub['WSr'].values >= iws + ws_step/2) & (data_sub['SSTa'].values >= (isst + sst_step/2 + 273.15)) )
                        if np.sum(mask_int) == 0:
                            continue
                        idx_int = np.argwhere(mask_int)[:,0]
                        data_int = data_sub.loc[mask_int]
        
                    # Check if there are coefficients for the surrounding bins and load
                    # --iws & isst_near
                    coeffs_file_near1 = COEFF_PATH + "/sst/coeffs_sst_stage_2_wsbin_"+str(iws)+"_sstbin_"+str(isst_near)+".npy"
                    if os.path.isfile(coeffs_file_near1):
                        D21 = np.load(coeffs_file_near1)
                        isnan_D21 = False
                    elif os.path.isfile(coeffs_file):
                        D21 = D11.copy()
                        isnan_D21 = True
                    else:
                        D21 = np.full((nsst_input), fill_value=np.nan)
                        isnan_D21 = True
                    # --iws_near & isst
                    coeffs_file_near2 = COEFF_PATH + "/sst/coeffs_sst_stage_2_wsbin_"+str(iws_near)+"_sstbin_"+str(isst)+".npy"
                    if os.path.isfile(coeffs_file_near2):
                        D12 = np.load(coeffs_file_near2)
                        isnan_D12 = False
                    elif os.path.isfile(coeffs_file):
                        D12 = D11.copy()
                        isnan_D12 = True
                    else:
                        D12 = np.full((nsst_input), fill_value=np.nan)
                        isnan_D12 = True
                    # --iws_near & isst_near
                    coeffs_file_near3 = COEFF_PATH + "/sst/coeffs_sst_stage_2_wsbin_"+str(iws_near)+"_sstbin_"+str(isst_near)+".npy"
                    if os.path.isfile(coeffs_file_near3):
                        D22 = np.load(coeffs_file_near3)
                        isnan_D22 = False
                    elif os.path.isfile(coeffs_file):
                        D22 = D11.copy()
                        isnan_D22 = True
                    else:
                        D22 = np.full((nsst_input), fill_value=np.nan)
                        isnan_D22 = True

                    if ( (not os.path.isfile(coeffs_file)) & (not os.path.isfile(coeffs_file_near1)) & (not os.path.isfile(coeffs_file_near2)) & (not os.path.isfile(coeffs_file_near3)) ):
                        # No coefficients for current or closest bins -> no retrieval
                        if verbose:
                            print("Warning: Coefficients do not exist for current wind speed bin {} and sst bin {} or closest bins. SST=NaN...!".format(iws,isst))
                        SSTr[idx_sub[idx_int]] = np.nan
                        continue
        
                    # Retrieve SST
                    SSTr_i11 = retrieve_sst(data_int,D11,sensor)
                    SSTr_i21 = retrieve_sst(data_int,D21,sensor)
                    SSTr_i12 = retrieve_sst(data_int,D12,sensor)
                    SSTr_i22 = retrieve_sst(data_int,D22,sensor)

                    # Define interpolation weights
                    alpha = np.abs(data_int['SSTa'].values - (isst + 273.15))/sst_step
                    beta  = np.abs(data_int['WSr'].values - iws)/ws_step
                    w11 = (1 - alpha)*(1 - beta)
                    w21 = alpha*(1 - beta)
                    w12 = (1 - alpha)*beta
                    w22 = alpha*beta

                    # Reset weights as we don't have either D11, D21, D12 or D22 coefficients
                    if (isnan_D11 | isnan_D21 | isnan_D12 | isnan_D22):
                        if (isnan_D11): w11 = 0.
                        if (isnan_D21): w21 = 0.
                        if (isnan_D12): w12 = 0.
                        if (isnan_D22): w22 = 0.

                        wsum = w11 + w21 + w12 + w22
                        w11 = w11 / wsum
                        w21 = w21 / wsum
                        w12 = w12 / wsum
                        w22 = w22 / wsum

                    # Interpolate between SSTr_i1, SSTr_i2, SSTr_i3 and SSTr_i4
                    SSTr_int = w11 * SSTr_i11 + w21 * SSTr_i21 + w12 * SSTr_i12 + w22 * SSTr_i22

                    # Assign to the correct SSTr elements
                    SSTr[idx_sub[idx_int]] = SSTr_int.copy()

    return SSTr



def calculate_coeffs_global_sst(data,sensor,verbose=False):
    print("\nCalculate SST global coefficients")

    # Exclude WSr nans
    data_sub = data.dropna(axis=0)

    # Get coefficients
    coeffs_all = calculate_coeffs_sst(data_sub,sensor)

    if not np.all(np.isnan(coeffs_all)):
        if verbose:
            print("Save SST global coefficients")
        coeffs_file = COEFF_PATH + "/sst/coeffs_sst_global.npy"
        np.save(coeffs_file,coeffs_all)



def retrieve_global_sst(data,sensor):
    print("\nRetrieve global SST")

    # Load coefficients
    coeffs_file = COEFF_PATH + "/sst/coeffs_sst_global.npy"
    A = np.load(coeffs_file)
    
    # Retrieve WS
    SSTr =  retrieve_sst(data,A,sensor)

    return SSTr



def print_stats(data,SSTr):
    dif_mean = np.nanmean(SSTr - data['insitu_sst'])
    dif_std = np.nanstd(SSTr - data['insitu_sst'])
    dif_rmse = np.sqrt(np.nanmean((SSTr - data['insitu_sst'])**2))
    
    print('MW retrievals - drifter observations')
    print('-------------------------------')
    print('Mean: ' + str(round(dif_mean,4)))
    print('St.d: ' + str(round(dif_std,4)))
    print('RMSE: ' + str(round(dif_rmse,4)))


def plot_scatter_data(data):
    fig1, ax1 = plt.subplots()
    ax1.scatter(data['insitu_sst']-273.15,data['SSTa']-273.15)
    ax1.plot([-3,40],[-3,40], linestyle='-', color='k', linewidth=1.2)
    ax1.set_title('1st-stage SST')
    ax1.set_xlabel('in situ SST ($^{\circ}$C)')
    ax1.set_ylabel('PMW SST ($^{\circ}$C)')
#    ax1.set_xlim([-2, 34])
#    ax1.set_ylim([-2, 34])

    fig2, ax2 = plt.subplots()
    ax2.scatter(data['insitu_sst']-273.15,data['SSTr']-273.15)
    ax2.plot([-3,40],[-3,40], linestyle='-', color='k', linewidth=1.2)
    ax2.set_title('2nd-stage SST')
    ax2.set_xlabel('in situ SST ($^{\circ}$C)')
    ax2.set_ylabel('PMW SST ($^{\circ}$C)')
#    ax2.set_xlim([-2, 34])
#    ax2.set_ylim([-2, 34])

    fig3, ax3 = plt.subplots()
    ax3.scatter(data['insitu_sst']-273.15,data['SSTr_global']-273.15)
    ax3.plot([-3,40],[-3,40], linestyle='-', color='k', linewidth=1.2)
    ax3.set_title('global SST')
    ax3.set_xlabel('in situ SST ($^{\circ}$C)')
    ax3.set_ylabel('PMW SST ($^{\circ}$C)')
#    ax3.set_xlim([-2, 34])
#    ax3.set_ylim([-2, 34])

    plt.show()






# =============#
# MAIN PROGRAM #
# =============#
if __name__ == "__main__":
    # =====================================================#
    # SETTINGS - SETTINGS - SETTINGS - SETTINGS - SETTINGS #
    # =====================================================#
    derive_stage_1_ws_coeffs = True
    derive_stage_2_ws_coeffs = True
    derive_stage_1_sst_coeffs = True
    derive_stage_2_sst_coeffs = True
    derive_global_sst_coeffs = True
    plot_timeseries = False
    plot_scatter = False
    
    # =============================================================================
    # LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD DATA - LOAD 
    # =============================================================================
    # Load MMD data from regression
    #-------------------------------
    # Settings
    sensor = "amsr2"
    mmd_type = "6b"
    data_type = "mmd"
    year = 2015

    # Wind speed bins
    ws_bins = np.arange(0,20,1, dtype=int)
    # Latitude bins
    lat_bins = np.arange(-72,84,2, dtype=int)
    # Orbit bins
    orb_bins = np.arange(0,2,1, dtype=int)
    # SST bins
    sst_bins = np.arange(-2,36,2, dtype=int)
    
    # Get variable names
    var_names = get_vars2extract(data_type)
    
    # Data
    data_file = DATA_PATH + "MMD" + mmd_type + "_drifter_" + str(year) + ".nc"
    ncid = nc.Dataset(data_file, mode='r', format="NETCDF4_CLASSIC")
    
    # Number of matchups
    nmatchups = ncid.dimensions['matchups'].size
    
    # Get the data
    data = pd.DataFrame(index=np.arange(nmatchups))
    for ivar,var_name in enumerate(var_names):
        data[var_name] = ncid[var_name][:]
    
    # Close the netCDF file
    ncid.close()
    
    
    # Process the data
    data['insitu_datetime'] = pd.to_datetime(data['insitu_time'],unit='s')
    data['look_angle'] = data['satza'] - 55.
    
    
    # Filter the data
    # -Remove nans
    data.dropna(axis=0,inplace=True)
    
    
    # Divide the data into train and test data
    data_train, data_test = train_test_split(data, test_size=0.3, random_state=rnseed)
    data_train.reset_index(inplace=True)
    data_test.reset_index(inplace=True)
    
    
    # =============================================================================
    # MULTIPLE LINEAR REGRESSION - MULTIPLE LINEAR REGRESSION - MULTIPLE LINEAR REG
    # =============================================================================
    
    #============#
    # WIND SPEED #
    #============#
    # ------------------------------
    # Wind speed retrieval - Stage 1
    # ------------------------------
    # Calculate coefficients
    if derive_stage_1_ws_coeffs:
        calculate_coeffs_stage_1_ws(data_train,sensor)
    
    # Retrieve WSa
    WSa = retrieve_stage_1_ws(data_train,sensor)
    data_train['WSa'] = WSa
    WSa = retrieve_stage_1_ws(data_test,sensor)
    data_test['WSa'] = WSa
    
    
    # ------------------------------
    # Wind speed retrieval - Stage 2
    # ------------------------------
    
    # Calculate coefficients
    if derive_stage_2_ws_coeffs:
        calculate_coeffs_stage_2_ws(data_train,ws_bins,sensor)
    
    # Retrieve WSr
    WSr = retrieve_stage_2_ws(data_train,ws_bins,sensor)
    data_train['WSr'] = WSr
    WSr = retrieve_stage_2_ws(data_test,ws_bins,sensor)
    data_test['WSr'] = WSr
    
    
    
    #=========================#
    # SEA SURFACE TEMPERATURE #
    #=========================#
    # -----------------------
    # SST retrieval - Stage 1
    # -----------------------
    
    # Calculate coefficients
    if derive_stage_1_sst_coeffs:
        calculate_coeffs_stage_1_sst(data_train,orb_bins,lat_bins,sensor)
    
    # Retrieve SSTa
    SSTa = retrieve_stage_1_sst(data_train,orb_bins,lat_bins,sensor)
    data_train['SSTa'] = SSTa
    SSTa = retrieve_stage_1_sst(data_test,orb_bins,lat_bins,sensor)
    data_test['SSTa'] = SSTa
    
    
    # -----------------------
    # SST retrieval - Stage 2
    # -----------------------
    
    # Calculate coefficients
    if derive_stage_2_sst_coeffs:
        calculate_coeffs_stage_2_sst(data_train,ws_bins,sst_bins,sensor)
    
    # Retrieve SSTr
    SSTr = retrieve_stage_2_sst(data_train,ws_bins,sst_bins,sensor)
    data_train['SSTr'] = SSTr
    SSTr = retrieve_stage_2_sst(data_test,ws_bins,sst_bins,sensor)
    data_test['SSTr'] = SSTr
    
    
    # -----------------------
    # SST retrieval - Global
    # -----------------------
    # Calculate coefficients
    if derive_global_sst_coeffs:
        calculate_coeffs_global_sst(data_train,sensor)
    
    # Retrieve SSTr
    SSTr = retrieve_global_sst(data_train,sensor)
    data_train['SSTr_global'] = SSTr
    SSTr = retrieve_global_sst(data_test,sensor)
    data_test['SSTr_global'] = SSTr



    #======================================================#
    # ANALYSIS - ANALYSIS - ANALYSIS - ANALYSIS - ANALYSIS #
    #======================================================#
    # Only use the common data
    good_data = ( (~np.isnan(data_train['SSTa'].values)) & (~np.isnan(data_train['SSTa'].values)) & (~np.isnan(data_train['SSTr_global'].values)) )
    data_train = data_train.loc[good_data,:]
    data_train.reset_index(inplace=True,drop=True)
    good_data = ( (~np.isnan(data_test['SSTa'].values)) & (~np.isnan(data_test['SSTa'].values)) & (~np.isnan(data_test['SSTr_global'].values)) )
    data_test = data_test.loc[good_data,:]
    data_test.reset_index(inplace=True,drop=True)

    # Print statistics
    print("\nTRAIN DATASET\n===============================")
    print("1-stage SST retrieval")
    print_stats(data_train,data_train['SSTa'])
    print("\n2-stage SST retrieval")
    print_stats(data_train,data_train['SSTr'])
    print("\nGlobal SST retrieval")
    print_stats(data_train,data_train['SSTr_global'])
    print("\nTEST DATASET\n===============================")
    print("1-stage SST retrieval")
    print_stats(data_test,data_test['SSTa'])
    print("\n2-stage SST retrieval")
    print_stats(data_test,data_test['SSTr'])
    print("\nGlobal SST retrieval")
    print_stats(data_test,data_test['SSTr_global'])


    # Plot retrieved vs in situ SST
    if plot_scatter:
        plot_scatter_data(data_train)
    
    
    
