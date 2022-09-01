"""
Save MMD from pickle format to netCDF format for the DEVALGO project
"""

import datetime
import netCDF4 as nc
import numpy as np
import os
import pickle
import sys


from read_mmd import read_mmd


if __name__ == "__main__":
    # Paths
    data_dir = "/media/ea/Elements SE/Emy/Projects/PhD/Data/ML_PMW_SST/"

    # Settings
    sensor = 'amsr2'
    mmd_type = '6b'
    data_type = 'mmd_all'
    year = 2015

    print("Read the data from mat-file format")
    data_file = data_dir + "MMD" + mmd_type + "_drifter_" + str(year) + "_era5.mat"
    MMD = read_mmd(data_file,sensor,data_type,year,year)

    # Number of matchups
    nmatchups = MMD.shape[0]

    # Only keep some of the columns 
    cols_to_keep = ['orbit', 'lat', 'lon', 'solza', 'satza', 'solaz', 'sataz', 'era5_wind_dir', \
                    'ccmp_wind_dir', 'era5_phi_rel', 'era5_wind_speed', 'ccmp_wind_speed', 'era5_sst', 'era5_tcwv', 'era5_clwt', \
                    'tb6vpol', 'tb6hpol', 'tb10vpol', 'tb10hpol', 'tb18vpol', 'tb18hpol', 'tb23vpol', 'tb23hpol', 'tb36vpol', \
                    'tb36hpol', 'tb89vpol', 'tb89hpol', 'sga', 'sss', 'insitu_sst', 'insitu_time']
    MMD = MMD[cols_to_keep]

    # Rename columns
    rename_dict = {'era5_wind_speed':'era5_ws', 'ccmp_wind_speed':'ccmp_ws'}
    MMD.rename(columns=rename_dict, inplace=True)
    

    # MMD variables
    mmd_vars = list(MMD.columns)


    print("Create netCDF file")
    # netCDF file
    ncfile = data_dir + "MMD" + mmd_type + "_drifter_" + str(year) + ".nc"

    # Create as empty if it doesn't exist
    ncid = nc.Dataset(ncfile, mode='w', format="NETCDF4_CLASSIC")

    # Create dimensions
    matchups_dim = ncid.createDimension("matchups",nmatchups)


    # Create attributes
    ncid.title = "ESA CCI MMD dataset for year " + str(year)
    ncid.history = "Created " + datetime.datetime.today().strftime("%Y-%d-%m")

    # Create variables
    nc_vars = mmd_vars
    # Standard name
    nc_standard_name = ['satellite_orbit', 'latitude', 'longitude', 'solar_zenith_angle', 'sensor_zenith_angle', 'solar_azimuth_angle', \
                        'sensor_azimuth_angle', 'wind_from_direction', 'wind_from_direction', 'relative_angle', 'wind_speed', 'wind_speed', \
                        'sea_surface_temperature', 'total_column_water_vapor', 'total_column_liquid_water', 'brightness_temperature', \
                        'brightness_temperature', 'brightness_temperature', 'brightness_temperature', 'brightness_temperature', \
                        'brightness_temperature', 'brightness_temperature', 'brightness_temperature', 'brightness_temperature', \
                        'brightness_temperature', 'brightness_temperature', 'brightness_temperature', 'sunglint_angle', 'sea_surface_salinity', \
                        'sea_surface_temperature', 'time']
    # Long name
    nc_long_name = ['satellite orbit', 'latitude', 'longitude', 'solar zenith angle', 'satellite zenith angle', 'solar azimuth angle', \
                    'satellite azimuth angle', 'ERA5 wind direction', 'CCMP wind direction', 'relative angle between wind direction and azimuth angle', \
                    'ERA5 wind speed', 'CCMP wind speed','ERA5 sea surface temperature', 'ERA5 total column water vapor', 'ERA5 total column liquid water', \
                    'brightness temperature of 6V', 'brightness temperature of 6H', 'brightness temperature of 10V', 'brightness temperature of 10H', \
                    'brightness temperature of 18V', 'brightness temperature of 18H', 'brightness temperature of 23V', 'brightness temperature of 23H', \
                    'brightness temperature of 36V', 'brightness temperature of 36H', 'brightness temperature of 89V', 'brightness temperature of 89V', \
                    'sun glint angle', 'CMEMS sea surface salinity', 'in situ sea surface temperature', 'in situ time']

    # Units
    nc_units = ['', 'degree_north', 'degree_east', 'degree', 'degree', 'degree', 'degree', 'degree', 'degree', 'degree', 'm s-1', 'm s-1', 'Kelvin', 'kg m-2', 'g m-3', \
                'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', 'Kelvin', \
                'Kelvin', 'Kelvin', 'degree', 'g kg-1', 'Kelvin', 'seconds']

    # Type
    nc_var_types = [np.int32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, \
                    np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, \
                    np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float64]

    
    # Add variables
    for ivar,var_name in enumerate(nc_vars):
        var = ncid.createVariable(var_name, nc_var_types[ivar], ('matchups'))
        var.units = nc_units[ivar]
        var.standard_name = nc_standard_name[ivar]
        var.long_name = nc_long_name[ivar]
        var[:] = MMD[var_name].to_numpy()
    
