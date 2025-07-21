#Calling libraries
import argparse
import intake
import netCDF4 as nc
import xarray as xr
import numpy as np
import pandas as pd
import copy
import os
import re
import rasterio
import geopandas
import rasterio.plot
import rioxarray
from shapely.geometry import mapping, Polygon
import calendar
import statsmodels.api as sm
import datetime as dt
import scipy.stats as ss
from glob import glob
import xesmf as xe
from pyproj import Transformer, CRS, transform, Proj
from sklearn.neighbors import BallTree


########
#Loads ACCESS-OM2-01 sea ice and ocean data for the Southern Ocean. If ice data is accessed, it corrects the time and coordinate grid to match ocean outputs.
def getACCESSdata_SO(var, start, end, freq, catalog, minlat = -90, maxlat = -45, 
                     exp = '01deg_jra55v140_iaf_cycle4', ice_data = False):
    '''
    Defining function that loads data automatically using in the intake catalog. The
    `getACCESSdata` will achieve the following:  
    - Access data for the experiment and variable of interest at the frequency requested and 
    within the time frame specified  
    - Apply **time corrections** as midnight (00:00:00) is interpreted differently by the
    CICE model and the xarray package.
    - CICE reads *2010-01-01 00:00:00* as the start of 2010-01-01, while xarray interprets 
    it as the start of the following day (2010-01-02). To fix this problem, 12 hours are 
    subtracted from the time dimension (also known as *time coordinate*).  
    - Latitude and longitude will be corrected in the dataset using the `geolon_t` dataset. 
    The coordinate names are replaced by names that are more intuitive.  
    - Minimum and maximum latitudes and longitudes can be specified in the function to 
    access specific areas of the dataset if required.  The **Southern Ocean** is defined as 
    ocean waters south of 45S.

    Inputs:
    var - Short name for the variable of interest
    start - Time from when data has to be returned
    end - Time until when data has to be returned
    freq - Time frequency of the data
    catalog - Intake catalog
    minlat - minimum latitude from which to return data. If not set, defaults to -90 to 
    cover the Southern Ocean.
    maxlat - maximum latitude from which to return data. If not set, defaults to -45 to 
    cover the Southern Ocean.
    exp - Experiment name. Default is 01deg_jra55v140_iaf_cycle4.
    ice_data - Boolean, when True the variable being called is related to sea ice, when 
    False is not. Default is set to False (i.e., it assumes variable is related to the 
    ocean).
        
    Output:
    Data array with corrected time and coordinates within the specified time period and 
    spatial bounding box.
    '''
    
    #If data being accessed is an ice related variable, then apply the following steps
    if ice_data == True:
        #Accessing data
        vararray = (catalog[exp].search(variable = var).
                    to_dask(xarray_open_kwargs = 
                            dict(use_cftime = True, decode_coords = False)).
                    sel(time = slice(start, end)))[var]
        #Accessing corrected coordinate data to update geographical coordinates in the array of interest
        area_t = catalog[exp].search(variable = 'area_t').to_dask()['area_t']
        #Apply time correction so data appears in the middle (12:00) of the day rather than at the beginning of the day (00:00)
        vararray['time'] = vararray.time.to_pandas() - dt.timedelta(hours = 12)
        #Change coordinates so they match ocean dimensions 
        vararray.coords['ni'] = area_t['xt_ocean'].values
        vararray.coords['nj'] = area_t['yt_ocean'].values
        #Rename coordinate variables so they match ocean data
        vararray = vararray.rename(({'ni':'xt_ocean', 'nj':'yt_ocean'}))
        #Drop coordinates that are no longer needed
        if len(vararray.coords) > 3:
            vararray = vararray.drop([i for i in vararray.coords if i not in ['time', 'xt_ocean', 'yt_ocean']])
    else:
        #Accessing data
        vararray = (catalog[exp].search(variable = var).
                    to_dask(xarray_open_kwargs =
                            dict(use_cftime = True)).sel(time = slice(start, end)))[var]
    #Subsetting data to area of interest
    if vararray.name in ['u', 'v']:
        vararray = vararray.sel(yu_ocean = slice(minlat, maxlat))
    else:
        vararray = vararray.sel(yt_ocean = slice(minlat, maxlat))
    return vararray

########
#Correcting longitude values in a data array so they are between -180 and +180 degrees
def corrlong(array):
    '''
    Inputs:
    array - Data array on which longitude corrections will be applied.
    
    Output:
    Data array with corrected longitude values.
    '''
    
    if array.name in ['u', 'v']:
        long_name = 'xu_ocean'
    else:
        long_name = 'xt_ocean'
    
    #Apply longitude correction
    array[long_name] = ((array[long_name] + 180)%360)-180
    array = array.sortby(array[long_name])
    
    return array


########
#Correcting longitude values in a data array so they are between -180 and +180 degrees
def extract_bottom_layer(da):
    '''
    Inputs:
    da - Data array from which bottom layer needs to be extracted.
    
    Output:
    Data array with a single depth layer.
    '''

    #Give a value of 1 to all cells containing environmental data in a single time step
    mask_2d = xr.where(~np.isnan(da.isel(time = 0)), 1, np.nan)
    #Perform a cumulative sum along depth axis to identify deepest grid cell with environmental data
    mask_2d = mask_2d.cumsum('st_ocean').where(~np.isnan(da.isel(time = 0)))
    #Create a mask identifying deepest cells with a value of 1
    mask_2d = xr.where(mask_2d == mask_2d.max('st_ocean'), 1, np.nan)
    #Apply mask to original data array
    da = (mask_2d*da).sum('st_ocean')
    #Rearrange dimensions to match original dataset
    da = da.transpose('time', 'yt_ocean', 'xt_ocean')
    #Returning bottom layer
    return da


########
#This function calculates distance from each grid cell to its nearest neighbour in a reference data array. Nearest neighbour refers to the search of the point within a predetermined set of points that is located closest (spatially) to a given point.
def nn_dist(target_da, grid_coords_numpy, **kwargs):
    '''
    Inputs:
    target_da (data array) - Reference points to which nearest neighbour distances will be calculated. Maximum values along y axis will be used as reference points.
    grid_coords_numpy (np data array) - Coordinate pairs for each grid cell from which nearest neighbour distance will be calculated
    Optional:
    folder_out (string) - Path to folder where output will be saved
    file_base (string) - Base name to be used to save outputs
    
    Output:
    Data array with distance to nearest neighbour
    '''
    #Getting coordinate pairs for sea ice edge
    ice_coords = np.vstack([target_da.yt_ocean[target_da.argmax(dim = 'yt_ocean')],
                            target_da.xt_ocean]).T
    
    #Set up Ball Tree (nearest neighbour algorithm).
    ball_tree = BallTree(np.deg2rad(ice_coords), metric = 'haversine')
    #The nearest neighbour calculation will give two outputs: distances in radians and indices
    dist_rad, ind = ball_tree.query(grid_coords_numpy, return_distance = True)
    #Transform distances from radians to km and changing data to data array
    earth_radius_km = 6371
    dist_km = xr.DataArray(data = [(dist_rad*earth_radius_km).reshape(target_da.shape)],
                           dims = ['time', 'yt_ocean', 'xt_ocean'],
                           coords = {'time': [target_da.time.values],
                                     'yt_ocean': target_da.yt_ocean.values,
                                     'xt_ocean': target_da.xt_ocean.values},
                           name = 'dist_km')
    dist_km = dist_km.assign_attrs({'units': 'km',
                          'long_name': 'distance to nearest neighbour'})
    
    #If path to folder provided, then save output
    if 'folder_out' in kwargs.keys():
        #Make sure output folder exists
        os.makedirs(kwargs.get('folder_out'), exist_ok = True)
        if 'file_base' in kwargs.keys():
            #Extract year and month to use in filename
            month = str(dist_km.time.dt.month.values[0]).zfill(2)
            year = dist_km.time.dt.year.values[0]
            file_base = kwargs.get('file_base')
            file_out = os.path.join(kwargs.get('folder_out'), 
                                f'{file_base}_{year}-{month}.nc')
            dist_km.to_netcdf(file_out)
        else:
            'File name base is needed to save output. Month and year will be added to this string'
        
    return dist_km


########
#This function creates a colour palette using Crameri's palettes (Crameri, F. (2018), Scientific colour-maps, Zenodo, doi:10.5281/zenodo.1243862)
def colourMaps(colourLibraryPath, palette, rev = True):
    '''
    Inputs:
    colourLibraryPath - the file path where the palettes are currently saved.
    palette - name of the palette to be created.
    rev - Boolean. If True, it will create a normal and reversed version of the palette. If False, it will only return one palette
    
    Outputs:
    One or two palettes based on Crameri (2018) that can be used to colour matplotlib figures
    '''
    #Load relevant libraries to set Scientific Colour Map library
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import ListedColormap

    #Set path where the scientific library is found
    cm_data = np.loadtxt(os.path.join(colourLibraryPath, palette, (palette + '.txt')))
    #Create a colour map based on 'palette' argument
    pal_map_adv = LinearSegmentedColormap.from_list(palette, cm_data)
        
    if rev == True:
        pal_map_ret = ListedColormap(cm_data[::-1])
        return pal_map_adv,pal_map_ret
    else:
        return pal_map_adv


########
#This function creates a single data frame with SDM outputs that can be used to create a data array for plotting
def df_ready(file_path, model, df_coords):
    '''
    Inputs:
    file_path - file path to data location
    model - name of the SDM algorithm used to create outputs
    df_coords - target grid to be used to create data frame
    
    Outputs:
    Data frame containing SDM predictions
    '''
    #Load csv file
    df = pd.read_csv(file_path)
    #Add SDM algorithm to data frame
    df['model'] = model
    #Keep relevant columns 
    df = df[['yt_ocean', 'xt_ocean', 'pred', 'month', 'model']]
    #Add coordinates from target grid
    df = df_coords.merge(df, on = ['xt_ocean', 'yt_ocean', 'month'], how = 'left')
    #Return data frame
    return df


#This function creates a single dataset with SDM outputs from a list of data frames
def ds_sdm(list_df, grid_sample, weights, weights_col):
    '''
    Inputs:
    list_df - a list of data frames to be used in dataset creation
    grid_sample - sample target grid. It must be two dimensional
    weights - weights to be applied to ensemble mean
    weights_chol - name of the column containing weights
    
    Outputs:
    Data frame containing SDM predictions and weighted ensemble mean
    '''
    #Create a single data frame with all predictions
    df = pd.concat(list_df)

    #Initialising empty dictionary to create dataset
    ds = {}

    #Looping through each month
    for m in df.month.unique():
        #Turn data frames into data arrays
        mth = df[df.month == m]
        mods = mth.model.dropna().unique()
        mth_da = xr.DataArray(data = mth.pred.values.reshape((len(mods),*grid_sample.shape)),
                              dims = ['model', 'yt_ocean', 'xt_ocean'],
                              coords = {'model': mods,
                                        'xt_ocean': grid_sample.xt_ocean.values,
                                        'yt_ocean': grid_sample.yt_ocean.values})
        
        #Calculate weighted ensemble
        ensemble = []
        #Loop through all models
        for mod in mods:
            #Extract model data 
            da_mod = mth_da.sel(model = mod)
            #Extract weights for model
            [w] = weights[weights_col][weights.model == mod].to_list()
            #Multiply model data by weight
            da_mod = da_mod*w
            #Attach weighted data to ensemble variable
            ensemble.append(da_mod)
        
        #Concatenate all data frames into single variable, add them up to get ensemble mean
        ensemble = xr.concat(ensemble, dim = 'model').sum('model').expand_dims({'model': ['Ensemble']})
        #Apply land mask from sample grid
        ensemble = xr.where(~np.isnan(grid_sample), ensemble, np.nan)
        
        #Add ensemble to all other models
        mth_da = xr.concat([mth_da, ensemble], dim = 'model')

        #Get name of month based on month digits
        mth_name = calendar.month_name[m]
        #Add to empty list to be used to create dataset
        ds[mth_name] = mth_da

    #Creating datasets
    ds = xr.Dataset(ds)
    
    #Return dataset
    return ds


#This function calculates gross growth potential (GGP) for Antarctic krill
def ggp_krill(log_10_mass):
    '''
    Inputs:
    log_10_mass - data array with krill dry mass at the beginning and end of the month
    
    Outputs:
    ggp - data array with monthly gross growth potential 
    '''
    #Initialise empty variables to save gross growth potential for Antarctic krill
    ggp = []
    #Calculate GGP by mass at the end of the month by mass at the bgeinning of the month
    for m, da in log_10_mass.groupby('time.month'):
        ggp_mth = da.isel(time = -1)/da.isel(time = 0)
        #Get year from dataset
        [yr] = np.unique(log_10_mass.time.dt.year)
        #Add date as middle of the month
        ggp_mth = ggp_mth.expand_dims({'time': [pd.to_datetime(f'{yr}-{m}-16')]})
        #Add results to variable
        ggp.append(ggp_mth)

    #Create data array with GGP
    ggp = xr.concat(ggp, dim = 'time')
    #Add name of variable
    ggp.name = 'krill_ggp'

    #Return data array
    return ggp


#This function calculates dry mass for Antarctic krill based on SST and phytoplankton biomass
def krill_growth(temp_df, phyto_df, **kwargs):
    '''
    Inputs:
    temp_df - data array with daily sea surface values in degrees C
    phyto_df - data array with daily phytoplankton biomass values in mg/m^3

    Optional (kwargs):
    starting_length (int) - starting length for individual krill in mm. Default value is 40 mm.
    
    Outputs:
    mth_gp - data array with monthly growth potential values
    ggp - data array with monthly gross growth potential 
    '''
    # Starting parameters
    alpha = -0.066
    beta = 0.002
    gamma = -0.000061
    delta = 0.385
    epsilon = 0.328
    zeta = 0.0078
    eta = -0.0101

    #Initialise empty variables to save daily growth rate and dry mass for beginning and end of month
    mth_gp = []
    log_10_mass = []

    #Daily calculations
    for t, p in zip(temp_df, phyto_df):
        #Get year, month and day from data array
        year = t.time.dt.year.values.tolist()
        month = t.time.dt.month.values.tolist()
        day = t.time.dt.day.values.tolist()
        #At the beginning of the month reset starting length 
        if (day == 1):
            #If krill starting length provided, then use that as starting length
            if 'starting_length' in kwargs.keys():
                L = kwargs.get('starting_length')
            #Otherwise use 40 mm
            else:
                L = 40
        #Calculate daily growth rates
        dgr = alpha + (beta*L) + (gamma*(L**2)) +  (delta*p/(epsilon+p)) + (zeta*t) + (eta*(t**2))
        #Add results to variable
        mth_gp.append(dgr)
        #Update length before next day calculations
        L = dgr+L
        #If beginning or end of month, then calculate dry mass
        if (day == 1) or (month == 11 and day == 30) or (month == 12 and day == 31):
            dry_mass = (3.89*(np.log10(L)))-4.19
            #Add time dimension
            dry_mass = dry_mass.expand_dims({'time': [t.time.values]})
            #Add results to variable
            log_10_mass.append(dry_mass)

    #Create data array with daily growth rates and calculate monthly growth potential
    mth_gp = xr.concat(mth_gp, dim = 'time').groupby('time.month').sum()
    time = [pd.to_datetime(f'{year}-{m}-16') for m in mth_gp.month.values]
    mth_gp.coords['month'] = time
    mth_gp = mth_gp.rename({'month': 'time'})

    #Create data array with dry mass
    log_10_mass = xr.concat(log_10_mass, dim = 'time')
    #Calculate gross growth potential for Antarctic krill
    ggp = ggp_krill(log_10_mass)

    #Apply land mask to growth rate calculations
    mth_gp = xr.where(np.isnan(ggp), np.nan, mth_gp)
    #Add variable name
    mth_gp.name = 'krill_growth_rate'
    
    #Return monthly growth and GGP
    return mth_gp, ggp

#Calculate the lat-lon coordinates from a dataset in source_crs - Function by Scott Wales
def calculate_latlon_coords(da, source_crs, target_crs):
    '''
    Inputs:
    da - refers to a data array that needs to be reprojected. Dimensions containing spatial data should be labelled as x and y.
    source_crs - original CRS for data array spatial information. Should be provided as a string in the form of 'epsg:4326'.
    target_crs - CRS to which the data array will be transformed. Should be provided as a string in the form of 'epsg:4326'.
    '''
    
    # Convert the 1d coordinates to 2d arrays covering the whole grid
    X, Y = np.meshgrid(da.x, da.y)
    
    # Use proj to create a transformation from the source coordinates to lat/lon
    trans = Transformer.from_crs(source_crs, target_crs)
    
    # Convert the 2d coordinates from the source to the target values
    lat, lon = trans.transform(X, Y)
    
    # Add the coordinates to the dataset
    da.coords['lat'] = (('y','x'), lat)
    da.coords['lon'] = (('y','x'), lon)
    
    # And add standard metadata so the lat and lon get picked up by xesmf
    da.coords['lat'].attrs['units'] = 'degrees_north'
    da.coords['lon'].attrs['units'] = 'degrees_east'
    da.coords['lat'].attrs['axis'] = 'Y'
    da.coords['lon'].attrs['axis'] = 'X'
    da.coords['lat'].attrs['standard_name'] = 'latitude'
    da.coords['lon'].attrs['standard_name'] = 'longitude'
    
    return da


########
def main(inargs):
    '''Run the program.'''

if __name__ == '__main__':
    description = 'This script contains functions used to perform timeseries within different sectors of the Southern Ocean.'
    parser = argparse.ArgumentParser(description = description)

    args = parser.parse_args()
    main(args)