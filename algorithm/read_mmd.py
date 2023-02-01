import datetime as dt
import h5py
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy
import sys
import time



def get_vars2extract(sensor):
    if sensor == "amsre":
        var_names = ['orbit', 'lat', 'lon', 'solza', 'satza', 'solaz', 'sataz', 'geostationary_reflection_lon', 'geostationary_reflection_lat', \
                     'pixel_quality6V', 'pixel_quality6H', 'pixel_quality10V', 'pixel_quality10H', 'pixel_quality18V', 'pixel_quality18H', \
                     'pixel_quality23V', 'pixel_quality23H', 'pixel_quality36V', 'pixel_quality36H', 'scan_quality', 'land_ocean_flag', 'seaice_fraction',
                     'era5_uwind', 'era5_vwind', 'ccmp_uwind', 'ccmp_vwind', 'era5_sst', 'era5_tcwv', 'era5_clwt', 'tb6vpol', 'tb6hpol', \
                     'tb10vpol', 'tb10hpol', 'tb18vpol', 'tb18hpol', 'tb23vpol', 'tb23hpol', 'tb36vpol', 'tb36hpol', 'tb89vpol', 'tb89hpol', \
                     'sga', 'sss', 'insitu_sst', 'insitu_time']
        mat_names = ['amsre__ascending', 'amsre__latitude', 'amsre__longitude', 'amsre__solar_zenith_angle', 'amsre__satellite_zenith_angle', \
                     'amsre__solar_azimuth_angle', 'amsre__satellite_azimuth_angle', 'amsre__Geostationary_Reflection_Longitude', \
                     'amsre__Geostationary_Reflection_Latitude', 'amsre__pixel_data_quality6V', 'amsre__pixel_data_quality6H', \
                     'amsre__pixel_data_quality10V', 'amsre__pixel_data_quality10H', 'amsre__pixel_data_quality18V', 'amsre__pixel_data_quality18H', \
                     'amsre__pixel_data_quality23V', 'amsre__pixel_data_quality23H', 'amsre__pixel_data_quality36V', 'amsre__pixel_data_quality36H', \
                     'amsre__scan_data_quality', 'amsre__land_ocean_flag_6', 'amsre__nwp__seaice_fraction', 'amsre__era_5__10m_east_wind_component', \
                     'amsre__era_5__10m_north_wind_component', 'amsre__ccmp__10m_east_wind_component', 'amsre__ccmp__10m_north_wind_component', \
                     'amsre__era_5__sea_surface_temperature', 'amsre__era_5__total_column_water_vapour', 'amsre__era_5__cloud_liquid_water_total', \
                     'amsre__brightness_temperature6V', 'amsre__brightness_temperature6H', 'amsre__brightness_temperature10V', 'amsre__brightness_temperature10H', \
                     'amsre__brightness_temperature18V', 'amsre__brightness_temperature18H', 'amsre__brightness_temperature23V', 'amsre__brightness_temperature23H', \
                     'amsre__brightness_temperature36V', 'amsre__brightness_temperature36H', 'amsre__brightness_temperature89V', 'amsre__brightness_temperature89H', \
                     'amsre__Sun_Glint_Angle', 'amsre__cmems__sea_surface_salinity', 'drifter___sst_insitu__sea_surface_temperature', 'drifter___sst_insitu__time']

    elif sensor == "amsr2":
        var_names = ['orbit', 'lat', 'lon', 'solza', 'satza', 'solaz', 'sataz',
                     'pixel_quality6', 'land_ocean_flag6', 'land_ocean_flag10', 'land_ocean_flag23', 'land_ocean_flag36', 'seaice_fraction',
                     'era5_uwind', 'era5_vwind', 'ccmp_uwind', 'ccmp_vwind', 'era5_sst', 'era5_tcwv', 'era5_clwt', 'tb6vpol', 'tb6hpol', \
                     'tb7vpol', 'tb7hpol', 'tb10vpol', 'tb10hpol', 'tb18vpol', 'tb18hpol', 'tb23vpol', 'tb23hpol', 'tb36vpol', 'tb36hpol', \
                     'tb89vpol', 'tb89hpol', 'sga', 'sss', 'insitu_sst', 'insitu_time']

        mat_names = ['amsr2__ascending', 'amsr2__latitude', 'amsr2__longitude', 'amsr2__solar_zenith_angle', 'amsr2__satellite_zenith_angle', \
                     'amsr2__solar_azimuth_angle', 'amsr2__satellite_azimuth_angle', 'amsr2__pixel_data_quality_6', 'amsr2__land_ocean_flag_6', \
                     'amsr2__Land_Ocean_Flag_10', 'amsr2__Land_Ocean_Flag_23', 'amsr2__Land_Ocean_Flag_36', 'amsr2__nwp__seaice_fraction', \
                     'amsr2__era_5__10m_east_wind_component', 'amsr2__era_5__10m_north_wind_component', 'amsr2__ccmp__10m_east_wind_component', \
                     'amsr2__ccmp__10m_north_wind_component', 'amsr2__era_5__sea_surface_temperature', 'amsr2__era_5__total_column_water_vapour', \
                     'amsr2__era_5__cloud_liquid_water_total', 'amsr2__brightness_temperature6V', 'amsr2__brightness_temperature6H', \
                     'amsr2__brightness_temperature7V', 'amsr2__brightness_temperature7H', 'amsr2__brightness_temperature10V', 'amsr2__brightness_temperature10H', \
                     'amsr2__brightness_temperature18V', 'amsr2__brightness_temperature18H', 'amsr2__brightness_temperature23V', 'amsr2__brightness_temperature23H', \
                     'amsr2__brightness_temperature36V', 'amsr2__brightness_temperature36H', 'amsr2__brightness_temperature89V', 'amsr2__brightness_temperature89H', \
                     'amsr2__Sun_Glint_Angle', 'amsr2__cmems__sea_surface_salinity', 'drifter___sst_insitu__sea_surface_temperature', 'drifter___sst_insitu__time']
    else:
        print("Could not parse data type: {}".format(sensor))
        print("Exiting...!")
        sys.exit()

    return var_names, mat_names



def find_rfi(iasc,xlon_ob,xlat_ob,xlon_geo,xlat_geo):
    """
    c.gentemann 7.3.2017
    Based on Gentemann&Hilburn 2014
    This code does not include ascension island, which affects entire scans.
    Output:
    + irfi: RFI mask (=0 NO RFI, ==1 RFI)
    Input:
    + iasc: ascending (=1) or descending (=0) passes
    + xlon_ob = long of observation on earth
    + xlat_ob = lat of observation point on earth
    + xlon_geo = reflection of observation lon
    + xlat_geo = reflection of observation lat
    """

    # Space based
    xlon_reflected = np.array([[235,279], [316,346], [340,360], [0,5], [0,50], [350,360], [0,10], \
                               [340,360], [0,5], [4,25], [289,305], [311,327], [295,315], [4,25]])
    xlat_reflected = np.array([[-22,19], [-16,10], [-15,20], [-14,9], [-30,20], [-10,14], [-10,14], \
                               [-15,4], [-15,4], [-10,13], [-8,12], [-4, 8], [-7, 10], [-10, 13]])
    asc_reflected = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

    # Ground based (without Ascension Island)
    xlon_ground= np.array([[1, 4], [5, 9], [3, 5], [11, 13], [70, 73], [1, 4], [5, 9], [3, 5], [2, 4], [70, 73], [11, 13]])
    xlat_ground= np.array([[58, 62], [63, 67], [52, 55], [34.5, 37], [17.5, 21], [58, 62], \
                           [63, 67], [52, 55], [56, 58], [17.5, 21], [34.5, 37]])
    asc_ground= np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0])

    Nlenr = xlon_reflected.shape[0]
    Nleng = xlon_ground.shape[0]
    Nleno = xlon_ob.shape[0]

    # Initialize rfi array
    irfi = np.zeros(Nleno)

    # Use longitudes in the rage [0, 360]
    xlon_ob_tmp = xlon_ob
    xlon_geo_tmp = xlon_geo
    xlon_ob_tmp[xlon_ob_tmp < 0] += 360.
    xlon_geo_tmp[xlon_geo_tmp < 0] += 360.

    # This is a lot faster (though uglier) than looping through the different set of coordinates!
    for i in range(Nleno):
        # Check if reflected vector at geostationary sat location
        if ( ( (iasc[i] == asc_reflected[0]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[0,0]) & (xlon_geo_tmp[i] <= xlon_reflected[0,1]) & \
               (xlat_geo[i]     >= xlat_reflected[0,0]) & (xlat_geo[i]     <= xlat_reflected[0,1]) ) | \
             ( (iasc[i] == asc_reflected[1]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[1,0]) & (xlon_geo_tmp[i] <= xlon_reflected[1,1]) & \
               (xlat_geo[i]     >= xlat_reflected[1,0]) & (xlat_geo[i]     <= xlat_reflected[1,1]) ) | \
             ( (iasc[i] == asc_reflected[2]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[2,0]) & (xlon_geo_tmp[i] <= xlon_reflected[2,1]) & \
               (xlat_geo[i]     >= xlat_reflected[2,0]) & (xlat_geo[i]     <= xlat_reflected[2,1]) ) | \
             ( (iasc[i] == asc_reflected[3]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[3,0]) & (xlon_geo_tmp[i] <= xlon_reflected[3,1]) & \
               (xlat_geo[i]     >= xlat_reflected[3,0]) & (xlat_geo[i]     <= xlat_reflected[3,1]) ) | \
             ( (iasc[i] == asc_reflected[4]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[4,0]) & (xlon_geo_tmp[i] <= xlon_reflected[4,1]) & \
               (xlat_geo[i]     >= xlat_reflected[4,0]) & (xlat_geo[i]     <= xlat_reflected[4,1]) ) | \
             ( (iasc[i] == asc_reflected[5]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[5,0]) & (xlon_geo_tmp[i] <= xlon_reflected[5,1]) & \
               (xlat_geo[i]     >= xlat_reflected[5,0]) & (xlat_geo[i]     <= xlat_reflected[5,1]) ) | \
             ( (iasc[i] == asc_reflected[6]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[6,0]) & (xlon_geo_tmp[i] <= xlon_reflected[6,1]) & \
               (xlat_geo[i]     >= xlat_reflected[6,0]) & (xlat_geo[i]     <= xlat_reflected[6,1]) ) | \
             ( (iasc[i] == asc_reflected[7]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[7,0]) & (xlon_geo_tmp[i] <= xlon_reflected[7,1]) & \
               (xlat_geo[i]     >= xlat_reflected[7,0]) & (xlat_geo[i]     <= xlat_reflected[7,1]) ) | \
             ( (iasc[i] == asc_reflected[8]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[8,0]) & (xlon_geo_tmp[i] <= xlon_reflected[8,1]) & \
               (xlat_geo[i]     >= xlat_reflected[8,0]) & (xlat_geo[i]     <= xlat_reflected[8,1]) ) | \
             ( (iasc[i] == asc_reflected[9]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[9,0]) & (xlon_geo_tmp[i] <= xlon_reflected[9,1]) & \
               (xlat_geo[i]     >= xlat_reflected[9,0]) & (xlat_geo[i]     <= xlat_reflected[9,1]) ) | \
             ( (iasc[i] == asc_reflected[10]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[10,0]) & (xlon_geo_tmp[i] <= xlon_reflected[10,1]) & \
               (xlat_geo[i]     >= xlat_reflected[10,0]) & (xlat_geo[i]     <= xlat_reflected[10,1]) ) | \
             ( (iasc[i] == asc_reflected[11]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[11,0]) & (xlon_geo_tmp[i] <= xlon_reflected[11,1]) & \
               (xlat_geo[i]     >= xlat_reflected[11,0]) & (xlat_geo[i]     <= xlat_reflected[11,1]) ) | \
             ( (iasc[i] == asc_reflected[12]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[12,0]) & (xlon_geo_tmp[i] <= xlon_reflected[12,1]) & \
               (xlat_geo[i]     >= xlat_reflected[12,0]) & (xlat_geo[i]     <= xlat_reflected[12,1]) ) | \
             ( (iasc[i] == asc_reflected[13]) & \
               (xlon_geo_tmp[i] >= xlon_reflected[13,0]) & (xlon_geo_tmp[i] <= xlon_reflected[13,1]) & \
               (xlat_geo[i]     >= xlat_reflected[13,0]) & (xlat_geo[i]     <= xlat_reflected[13,1]) ) ):
            irfi[i] = 1
        # Check if observation location at ground-based RFI location
        if ( ( (iasc[i] == asc_ground[0]) & \
               (xlon_ob_tmp[i] >= xlon_ground[0,0]) & (xlon_ob_tmp[i] <= xlon_ground[0,1]) & \
               (xlat_ob[i]     >= xlat_ground[0,0]) & (xlat_ob[i]     <= xlat_ground[0,1]) ) | \
             ( (iasc[i] == asc_ground[1]) & \
               (xlon_ob_tmp[i] >= xlon_ground[1,0]) & (xlon_ob_tmp[i] <= xlon_ground[1,1]) & \
               (xlat_ob[i]     >= xlat_ground[1,0]) & (xlat_ob[i]     <= xlat_ground[1,1]) ) | \
             ( (iasc[i] == asc_ground[2]) & \
               (xlon_ob_tmp[i] >= xlon_ground[2,0]) & (xlon_ob_tmp[i] <= xlon_ground[2,1]) & \
               (xlat_ob[i]     >= xlat_ground[2,0]) & (xlat_ob[i]     <= xlat_ground[2,1]) ) | \
             ( (iasc[i] == asc_ground[3]) & \
               (xlon_ob_tmp[i] >= xlon_ground[3,0]) & (xlon_ob_tmp[i] <= xlon_ground[3,1]) & \
               (xlat_ob[i]     >= xlat_ground[3,0]) & (xlat_ob[i]     <= xlat_ground[3,1]) ) | \
             ( (iasc[i] == asc_ground[4]) & \
               (xlon_ob_tmp[i] >= xlon_ground[4,0]) & (xlon_ob_tmp[i] <= xlon_ground[4,1]) & \
               (xlat_ob[i]     >= xlat_ground[4,0]) & (xlat_ob[i]     <= xlat_ground[4,1]) ) | \
             ( (iasc[i] == asc_ground[5]) & \
               (xlon_ob_tmp[i] >= xlon_ground[5,0]) & (xlon_ob_tmp[i] <= xlon_ground[5,1]) & \
               (xlat_ob[i]     >= xlat_ground[5,0]) & (xlat_ob[i]     <= xlat_ground[5,1]) ) | \
             ( (iasc[i] == asc_ground[6]) & \
               (xlon_ob_tmp[i] >= xlon_ground[6,0]) & (xlon_ob_tmp[i] <= xlon_ground[6,1]) & \
               (xlat_ob[i]     >= xlat_ground[6,0]) & (xlat_ob[i]     <= xlat_ground[6,1]) ) | \
             ( (iasc[i] == asc_ground[7]) & \
               (xlon_ob_tmp[i] >= xlon_ground[7,0]) & (xlon_ob_tmp[i] <= xlon_ground[7,1]) & \
               (xlat_ob[i]     >= xlat_ground[7,0]) & (xlat_ob[i]     <= xlat_ground[7,1]) ) | \
             ( (iasc[i] == asc_ground[8]) & \
               (xlon_ob_tmp[i] >= xlon_ground[8,0]) & (xlon_ob_tmp[i] <= xlon_ground[8,1]) & \
               (xlat_ob[i]     >= xlat_ground[8,0]) & (xlat_ob[i]     <= xlat_ground[8,1]) ) | \
             ( (iasc[i] == asc_ground[9]) & \
               (xlon_ob_tmp[i] >= xlon_ground[9,0]) & (xlon_ob_tmp[i] <= xlon_ground[9,1]) & \
               (xlat_ob[i]     >= xlat_ground[9,0]) & (xlat_ob[i]     <= xlat_ground[9,1]) ) | \
             ( (iasc[i] == asc_ground[10]) & \
               (xlon_ob_tmp[i] >= xlon_ground[10,0]) & (xlon_ob_tmp[i] <= xlon_ground[10,1]) & \
               (xlat_ob[i]     >= xlat_ground[10,0]) & (xlat_ob[i]     <= xlat_ground[10,1]) ) ):
            irfi[i] = 1

    return irfi



def find_rfi_ascension(xlon_ob,xlat_ob):
    """
    "Ascension Island affects entire scans.
    The paper Gentemann & Hilburn (2014) masked out asc/dsc ascension island
    RFI but I later realized it was entire scan lines, it wasn't getting into
    the main reflector, it was transmitting and getting into the hot load or
    cold mirror somehow. Better to flag entire scans using cold mirror
    check." - Chelle

    Output:
    + irfi: mask for Ascension Island RFI (=0 NO RFI, ==1 RFI)
    Input:
    + xlon_ob = long of observation on earth
    + xlat_ob = lat of observation point on earth
    """

    # Ascension Island area latitude and longitude
    lon_ascension = np.array([-24, -6])
    lat_ascension = np.array([-18, -2])

    # Get lengths of arrays
    N_ascension=len(lon_ascension)
    Nmatchups=xlon_ob.shape[0]

    # Initialize mask array
    irfi=np.zeros(Nmatchups)

    # Check if amsr2 observation location at ascension island rfi location
    mask = ( (xlon_ob >= lon_ascension[0]) & (xlon_ob <= lon_ascension[1]) & \
             (xlat_ob >= lat_ascension[0]) & (xlat_ob <= lat_ascension[1]) )
    irfi[mask] = 1

    return irfi



def find_diurnal_warming(wspd,sol_zeni,wspd_max):
    """
    Fnd diurnal warming contaminated data (=1 contaminated/erronous data, =0 good data)
    Input:
    + wspd: wind speed
    + sol_zen: solar zenith angle
    + wspd_max: maximum wind speed limit for diurnal warming contamination
    Output:
    + idw: Mask for diurnal warming contaminated data (=False bad data, =True good data)
    """

    # Number of matchups
    Nmatchups = wspd.shape[0]

    # Initialize output array
    idw = np.zeros(Nmatchups)
    idw[:] = True

    # Find the contaminated data
    contaminated_mask = ( (np.abs(sol_zeni) < 90. ) & (wspd < wspd_max) )
    idw[contaminated_mask] = False

    return idw



def nwpis_sigma_filter(SSTnwp,SSTi,time,year0,year1,goodidx):
    """
    Eliminate obviously erroneous data using a 3-sigma filter on the difter and nwp difference [MW 2008]
    Input:
    Output:
    """

    # Dummy index
    dum = 1

    # Get a vector with the years from the time array
    timestamps = pd.to_datetime(time, unit='s')#.year.values
    yyyy = pd.DatetimeIndex(timestamps).year.values
    for iyear in range(year0,year1+1):
        # NWP SST - insitu SST
        nwpisdiff = (SSTnwp) - SSTi
        # Get only the current year
        idx = goodidx[np.argwhere(goodidx[yyyy[goodidx] == iyear])]
        # Calculate std and mean
        Dstd = np.nanstd(nwpisdiff[idx]);
        Dmean = np.nanmean(nwpisdiff[idx]);

        # Get only the good indices
        goodidx2 = idx[np.abs(nwpisdiff[idx]-Dmean) < 3*Dstd];

        # Concatenate the indices arrays
        if dum == 1:
            icheck_nwpis = goodidx2
        else:
            icheck_nwpis = np.concatenate((icheck_nwpis,goodidx2))
        dum = dum + 1;

    return icheck_nwpis



def flagging_bad_data_mmd06c_drifter(MMDall,year0,year1):
    """
    Function for flagging erroneous NWP/in situ/AMSR-E data

    Input:
    + MMDall: Pandas DataFrame with the AMSR-E MMD dataset
    + year0: Start year in the time period we're looking at
    + year1: End year in the time period we're looking at
    Output:
    + good_data_mask: Mask for finding good (=True) and bad (=False) data
    """

    # ---------
    # Settings
    # ---------
    TBmax = 320               # max brightness temperature (TB)
    TBmin = 0                 # min brightness temperature (TB)
    windmax = 50              # max wind speed (gross error)
    lon_lim = 180             # longitude limit
    lat_lim = 90              # latitude limit
    solzen_min = 0            # minimum solar zenith angle
    solzen_max = 180          # maximum solar zenith angle
    TBdiffmax = 0             # maximum difference of 18-36 V-H (consistency check)
    wspd_min = 0              # minimum wind speed (m/s)
    wspd_max = 20             # maximum wind speed (m/s)
    SST_max = 34              # maximum SST (Celsius)
    SST_min = -2              # minimum SST (Celsius)
    rain = 240                # maximum 18V (for rain flagging)
    sun_glint_min = 25        # minimum sun glint angle
    dw_wspdmax = 4            # maximum wspd for diurnal warming flagging
    stdv23V_max = 55          # max stddev of 21x21 pixels for the 23V channel
    stdv23H_max = 35          # max stddev of 21x21 pixels for the 23H channel
    stdv36V_max = 25          # max stddev of 21x21 pixels for the 36V channel
    stdv36H_max = 25          # max stddev of 21x21 pixels for the 36H channel

    # Number of matchups
    Nmatchups = MMDall['lat'].values.shape[0]

    # Index array
    idx_array = np.arange(Nmatchups)

    # Initialize the mask
    good_data_array = np.ones(Nmatchups)


    # Calculate the ERA5 wind speed
    Wspd = np.sqrt(MMDall['era5_uwind'].values**2 + MMDall['era5_vwind'].values**2)


    #--------------------------------
    print('Construct good_data_array')
    #--------------------------------

    # Start execution time taking
    start_time = time.time()

    print(' - Find RFI contaminated matchups')
    icheck_rfi = find_rfi(MMDall['orbit'].values, \
                          MMDall['lon'].values, \
                          MMDall['lat'].values, \
                          MMDall['geostationary_reflection_lon'].values, \
                          MMDall['geostationary_reflection_lat'].values)
    # Get the RFI mask
    rfi_mask = (icheck_rfi == 0)

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))


    # Start execution time taking
    start_time = time.time()

    print(' - Find DW contaminated matchups')
    dw_mask = find_diurnal_warming(Wspd,MMDall['solza'].values,dw_wspdmax)

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))


    # Start execution time taking
    start_time = time.time()

    print('Filter/flag for erroneous matchups')
    # Filter/flag erroneous data and get the indices of the good data (good_data_idx1)
    init_good_mask = ( (MMDall['pixel_quality6V'].values               == 0            ) & \
                       (MMDall['pixel_quality6H'].values               == 0            ) & \
                       (MMDall['pixel_quality10V'].values              == 0            ) & \
                       (MMDall['pixel_quality10H'].values              == 0            ) & \
                       (MMDall['pixel_quality18V'].values              == 0            ) & \
                       (MMDall['pixel_quality18H'].values              == 0            ) & \
                       (MMDall['pixel_quality23V'].values              == 0            ) & \
                       (MMDall['pixel_quality23H'].values              == 0            ) & \
                       (MMDall['pixel_quality36V'].values              == 0            ) & \
                       (MMDall['pixel_quality36H'].values              == 0            ) & \
                       (MMDall['scan_quality'].values                  == 0            ) & \
                       (MMDall['tb6vpol'].values                       >  TBmin        ) & \
                       (MMDall['tb6vpol'].values                       <  TBmax        ) & \
                       (MMDall['tb6hpol'].values                       >  TBmin        ) & \
                       (MMDall['tb6hpol'].values                       <  TBmax        ) & \
                       (MMDall['tb10vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb10vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb10hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb10hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb18vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb18vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb18hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb18hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb23vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb23vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb23hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb23hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb36vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb36vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb36hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb36hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb89vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb89vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb89hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb89hpol'].values                      <  TBmax        ) & \
                       (np.abs(MMDall['lat'].values)                   <= lat_lim      ) & \
                       (np.abs(MMDall['lon'].values)                   <= lon_lim      ) & \
                       (MMDall['solza'].values                         >= solzen_min   ) & \
                       (MMDall['solza'].values                         <= solzen_max   ) & \
                       (MMDall['tb18vpol'].values - MMDall['tb18hpol'].values >= TBdiffmax ) & \
                       (MMDall['tb23vpol'].values - MMDall['tb23hpol'].values >= TBdiffmax ) & \
                       (MMDall['tb36vpol'].values - MMDall['tb36hpol'].values >= TBdiffmax ) & \
                       (Wspd                                           >= wspd_min     ) & \
                       (Wspd                                           <= wspd_max     ) & \
                       (MMDall['insitu_sst'].values - 273.15           >  SST_min      ) & \
                       (MMDall['insitu_sst'].values - 273.15           <= SST_max      ) & \
                       (MMDall['era5_sst'].values - 273.15             >  SST_min      ) & \
                       (MMDall['era5_sst'].values - 273.15             <= SST_max      ) & \
                       (MMDall['land_ocean_flag'].values               == 0            ) & \
                       (MMDall['seaice_fraction'].values               == 0            ) & \
                       (MMDall['sga'].values                           >  sun_glint_min) & \
                       (rfi_mask                                       == True         ) & \
                       (dw_mask                                        == True         ) & \
                       (MMDall['tb18vpol'].values                      <  rain         ) )

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))

    # Get the good indices
    goodidx = idx_array[init_good_mask]


    # Start execution time taking
    start_time = time.time()

    print(' - Perform a 3-sigma filter on NWP SST - in situ SST')
    icheck_nwpis = nwpis_sigma_filter(MMDall['era5_sst'].values,MMDall['insitu_sst'].values, \
                                      MMDall['insitu_time'],year0,year1,goodidx)
    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))

    # Get the common indices between the good data in init_good_mask and the indices from the 3-sigma filter
    good_data_idx = np.intersect1d(icheck_nwpis,goodidx)

    # Final mask - put the good data to 0
    good_data_array[good_data_idx] = 0

    # Get a mask (boolean) where good data = True
    good_data_mask = (good_data_array == 0)


    return good_data_mask



def flagging_bad_data_mmd06b_drifter(MMDall,year0,year1):
    """
    Function for flagging erroneous NWP/in situ/AMSR2 data

    Input:
    + MMDall: Pandas DataFrame with the AMSR2 MMD dataset
    + year0: Start year in the time period we're looking at
    + year1: End year in the time period we're looking at
    Output:
    + good_data_mask: Mask for finding good (=True) and bad (=False) data
    """

    # ---------
    # Settings
    # ---------
    TBmax = 320               # max brightness temperature (TB)
    TBmin = 0                 # min brightness temperature (TB)
    windmax = 50              # max wind speed (gross error)
    lon_lim = 180             # longitude limit
    lat_lim = 90              # latitude limit
    solzen_min = 0            # minimum solar zenith angle
    solzen_max = 180          # maximum solar zenith angle
    TBdiffmax = 0             # maximum difference of 18-36 V-H (consistency check)
    wspd_min = 0              # minimum wind speed (m/s)
    wspd_max = 20             # maximum wind speed (m/s)
    SST_max = 34              # maximum SST (Celsius)
    SST_min = -2              # minimum SST (Celsius)
    rain = 240                # maximum 18V (for rain flagging)
    sun_glint_min = 25        # minimum sun glint angle
    RFImax = 3                # maximum for RFI check on diff between 6-7 GHz
    dw_wspdmax = 6            # maximum wspd for diurnal warming flagging
    stdv23V_max = 55          # max stddev of 21x21 pixels for the 23V channel
    stdv23H_max = 35          # max stddev of 21x21 pixels for the 23H channel
    stdv36V_max = 25          # max stddev of 21x21 pixels for the 36V channel
    stdv36H_max = 25          # max stddev of 21x21 pixels for the 36H channel

    # Number of matchups
    Nmatchups = MMDall['lat'].values.shape[0]

    # Index array
    idx_array = np.arange(Nmatchups)

    # Initialize the mask
    good_data_array = np.ones(Nmatchups)


    # Calculate the ERA5 wind speed
    Wspd = np.sqrt(MMDall['era5_uwind'].values**2 + MMDall['era5_vwind'].values**2)


    #--------------------------------
    print('Construct good_data_array')
    #--------------------------------

    # Start execution time taking
    start_time = time.time()

    print(' - Find RFI contaminated matchups')
    # C-band RFI can be detected using the 6 and 7 GHZ channels.
    # If the difference between the two channels >= RFImax --> flag as RFI
    # contaminated (source: Alsweiss et al. 2016 - Remote Sensing of Sea
    # Surface Temperatrue Using AMSR-2 Measurements)
    icheck_rfi = np.ones(Nmatchups)   # RFI = 1
    icheck_rfi[( (np.abs(MMDall['tb7vpol'] - MMDall['tb6vpol']) < RFImax) & \
                 (np.abs(MMDall['tb7hpol'] - MMDall['tb6hpol']) < RFImax) )] = 0;
    # Need to mask out Ascension Island --> gives a lot of RFI
    irfi = find_rfi_ascension(MMDall['lon'],MMDall['lat'])
    asc_idx = (irfi == 1)
    icheck_rfi[asc_idx] = 1

    # Get the RFI mask - good data = 0
    rfi_mask = (icheck_rfi == 0)

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))


    # Start execution time taking
    start_time = time.time()

    print(' - Find DW contaminated matchups')
    dw_mask = find_diurnal_warming(Wspd,MMDall['solza'].values,dw_wspdmax)

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))


    # Start execution time taking
    start_time = time.time()

    print('Filter/flag for erroneous matchups')
    # Filter/flag erroneous data and get the indices of the good data (good_data_idx1)
    init_good_mask = ( (MMDall['pixel_quality6'].values                == 0            ) & \
                       (MMDall['tb6vpol'].values                       >  TBmin        ) & \
                       (MMDall['tb6vpol'].values                       <  TBmax        ) & \
                       (MMDall['tb6hpol'].values                       >  TBmin        ) & \
                       (MMDall['tb6hpol'].values                       <  TBmax        ) & \
                       (MMDall['tb7vpol'].values                       >  TBmin        ) & \
                       (MMDall['tb7vpol'].values                       <  TBmax        ) & \
                       (MMDall['tb7hpol'].values                       >  TBmin        ) & \
                       (MMDall['tb7hpol'].values                       <  TBmax        ) & \
                       (MMDall['tb10vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb10vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb10hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb10hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb18vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb18vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb18hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb18hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb23vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb23vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb23hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb23hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb36vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb36vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb36hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb36hpol'].values                      <  TBmax        ) & \
                       (MMDall['tb89vpol'].values                      >  TBmin        ) & \
                       (MMDall['tb89vpol'].values                      <  TBmax        ) & \
                       (MMDall['tb89hpol'].values                      >  TBmin        ) & \
                       (MMDall['tb89hpol'].values                      <  TBmax        ) & \
                       (np.abs(MMDall['lat'].values)                   <= lat_lim      ) & \
                       (np.abs(MMDall['lon'].values)                   <= lon_lim      ) & \
                       (MMDall['solza'].values                         >= solzen_min   ) & \
                       (MMDall['solza'].values                         <= solzen_max   ) & \
                       (MMDall['tb18vpol'].values - MMDall['tb18hpol'].values >= TBdiffmax    ) & \
                       (MMDall['tb23vpol'].values - MMDall['tb23hpol'].values >= TBdiffmax    ) & \
                       (MMDall['tb36vpol'].values - MMDall['tb36hpol'].values >= TBdiffmax    ) & \
                       (Wspd                                            >= wspd_min     ) & \
                       (Wspd                                            <= wspd_max     ) & \
                       (MMDall['insitu_sst'].values - 273.15            >  SST_min      ) & \
                       (MMDall['insitu_sst'].values - 273.15            <= SST_max      ) & \
                       (MMDall['era5_sst'].values - 273.15              >  SST_min      ) & \
                       (MMDall['era5_sst'].values - 273.15              <= SST_max      ) & \
                       (MMDall['land_ocean_flag6'].values               == 0            ) & \
                       (MMDall['land_ocean_flag10'].values              == 0            ) & \
                       (MMDall['land_ocean_flag23'].values              == 0            ) & \
                       (MMDall['land_ocean_flag36'].values              == 0            ) & \
                       (MMDall['seaice_fraction'].values                == 0            ) & \
                       (MMDall['sga'].values                            >  sun_glint_min) & \
                       (rfi_mask                                        == True         ) & \
                       (dw_mask                                         == True         ) & \
                       (MMDall['tb18vpol'].values                       <  rain         ) )

    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))

    # Get the good indices
    goodidx = idx_array[init_good_mask]


    # Start execution time taking
    start_time = time.time()

    print(' - Perform a 3-sigma filter on NWP SST - in situ SST')
    icheck_nwpis = nwpis_sigma_filter(MMDall['era5_sst'].values,MMDall['insitu_sst'].values, \
                                      MMDall['insitu_time'],year0,year1,goodidx)
    # End execution time taking
    end_time = time.time()
    print(" --- Execution time: %s seconds ---" % (end_time - start_time))

    # Get the common indices between the good data in init_good_mask and the indices from the 3-sigma filter
    good_data_idx = np.intersect1d(icheck_nwpis,goodidx)

    # Final mask - put the good data to 0
    good_data_array[good_data_idx] = 0

    # Get a mask (boolean) where good data = True
    good_data_mask = (good_data_array == 0)


    return good_data_mask





def calculate_angles(sataz,nwp_u10,nwp_v10):
    """
    Function for calculating angles:
    + wind direction (relative to noth), theta_w
    + satellite azimuth angle, theta_sat
    + wind direction relative to satellite aximuth angle, theta_rel
    Want to have the wind direction relative to north
    arctan2 gives an answer in the unit circle coordinates, which increases
    counterclockwise and have a zero on the x-axis (i.e. pi/2 to the north,
    -pi/2 to the south and +-pi to the west).

    However, we want to know the direction relative to north and not east
    --> need to subtract the atan2(v,u) answer from 90 deg. to get the
    direction the wind is blowing TOWARD relative to north
    """

    # Get satellite azimuth angle in range [0, 360]
    theta_sat = sataz
    theta_sat = np.mod(theta_sat,360)

    # Calculate direction the wind is blowing TOWARDS relative to north (in degrees)
    theta_w = 90. - ( np.arctan2(nwp_v10,nwp_u10) * 180./np.pi )
#    # % Calculate direction the wind is blowing FROM relative to north (in degrees)
#    theta_w = 270. - ( np.arctan2(nwp_v10,nwp_u10) * 180/np.pi )
    # Get theta_w on range [0, 360]
    theta_w = np.mod(theta_w,360)


    # Calculate the relative angle
    theta_rel = theta_sat - theta_w;
    # Get theta_w on range [0, 360]
    theta_rel = np.mod(theta_rel,360)

    return theta_rel, theta_w, theta_sat



def read_mmd(data_file,sensor,data_type,year0,year1):
    """
    Read AMSR-E MMD dataset and perform basic filtering checks
    """
    # Get variable names
    column_names, mat_names = get_vars2extract(sensor)

    # Open mat file
    f = h5py.File(data_file,'r')

    # Read the data
    alldata = np.array([np.array(f['MMDall/'+mat_names[i]]).squeeze() for i,cname in enumerate(column_names)]).T
    f.close()

    # Create a pandas dataframe
    MMD = pd.DataFrame(data=alldata, columns=column_names)

    # Change so that time is seconds since 1970 (epoch) and not 1978
    MMD['insitu_time'] = MMD['insitu_time'].values + 8*60*60*24*365

    # Change in situ SST to Kelvin
    MMD['insitu_sst'] = MMD['insitu_sst'].values + 273.15

    # Perform basic filtering
    if sensor == 'amsre':
        good_data_mask = flagging_bad_data_mmd06c_drifter(MMD,year0,year1)
    elif sensor == 'amsr2':
        good_data_mask = flagging_bad_data_mmd06b_drifter(MMD,year0,year1)

#    print('total=',good_data_mask.shape)
#    print('good=',np.sum(good_data_mask==1))
#    print('bad=',np.sum(good_data_mask==0))

    # Remove the bad data
    MMD = MMD.loc[(good_data_mask == 1),:]


    # Add extra variables
    # - Wind speed
    MMD['era5_wind_speed'] = np.sqrt(MMD['era5_uwind'].values**2 + MMD['era5_vwind'].values**2)
    MMD['ccmp_wind_speed'] = np.sqrt(MMD['ccmp_uwind'].values**2 + MMD['ccmp_vwind'].values**2)
    # - Wind direction and relative angle
    phi_rel, wind_dir, sataz = calculate_angles(MMD['sataz'].values,MMD['era5_uwind'].values,MMD['era5_vwind'])
    MMD['sataz'] = sataz
    MMD['era5_phi_rel'] = phi_rel
    MMD['era5_wind_dir'] = wind_dir
    _, wind_dir, _ = calculate_angles(MMD['sataz'].values,MMD['ccmp_uwind'].values,MMD['ccmp_vwind'])
    MMD['ccmp_wind_dir'] = wind_dir

    return MMD
