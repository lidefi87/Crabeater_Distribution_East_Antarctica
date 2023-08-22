#Calling libraries
import argparse
import cosima_cookbook as cc
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
#Defining functions

########
#Loads ACCESS-OM2-01 sea ice and ocean data for the Southern Ocean. If ice data is accessed, it corrects the time and coordinate grid to match ocean outputs.
def getACCESSdata_SO(var, start, end, freq, ses, minlat = -90, maxlat = -45, 
                  exp = '01deg_jra55v140_iaf_cycle4', ice_data = False):
    '''
    Defining function that loads data automatically using `cc.querying.getvar()` in a loop. The inputs needed are similar to those for the `cc.querying.getvar()` function, with the addition of inputs to define an area of interest.  
The `getACCESSdata` will achieve the following:  
- Access data for the experiment and variable of interest at the frequency requested and within the time frame specified  
- Apply **time corrections** as midnight (00:00:00) is interpreted differently by the CICE model and the xarray package.
    - CICE reads *2010-01-01 00:00:00* as the start of 2010-01-01, while xarray interprets it as the start of the following day (2010-01-02). To fix this problem, 12 hours are subtracted from the time dimension (also known as *time coordinate*).  
- Latitude and longitude will be corrected in the dataset using the `geolon_t` dataset. The coordinate names are replaced by names that are more intuitive.  
- Minimum and maximum latitudes and longitudes can be specified in the function to access specific areas of the dataset if required.  The **Southern Ocean** is defined as ocean waters south of 45S.

    Inputs:
    var - Short name for the variable of interest
    start - Time from when data has to be returned
    end - Time until when data has to be returned
    freq - Time frequency of the data
    ses - Cookbook session
    minlat - minimum latitude from which to return data. If not set, defaults to -90 to cover the Southern Ocean.
    maxlat - maximum latitude from which to return data. If not set, defaults to -45 to cover the Southern Ocean.
    exp - Experiment name. Default is 01deg_jra55v140_iaf_cycle4.
    ice_data - Boolean, when True the variable being called is related to sea ice, when False is not. Default is set to False (i.e., it assumes variable is related to the ocean).
        
    Output:
    Data array with corrected time and coordinates within the specified time period and spatial bounding box.
    '''
    
    #If data being accessed is an ice related variable, then apply the following steps
    if ice_data == True:
        #Accessing data
        vararray = cc.querying.getvar(exp, var, ses, frequency = freq, start_time = start, end_time = end, decode_coords = False)
        #Accessing corrected coordinate data to update geographical coordinates in the array of interest
        area_t = cc.querying.getvar(exp, 'area_t', ses, n = 1)
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
        vararray = cc.querying.getvar(exp, var, ses, frequency = freq, start_time = start, end_time = end)
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
#This function assigns the same coordinate reference system (CRS) to the data array being clipped as the clipping shapefile and then clips the data array
def clipDataArray(array, shp):
    '''
    Inputs:
    array - Data array to be clipped.
    shp - Shapefile to be used for clipping.
    
    Output:
    Clipped data array.
    '''
        
    #Set the spatial dimensions of the xarray being clipped
    array.rio.set_spatial_dims(x_dim = 'xt_ocean', y_dim = 'yt_ocean', inplace = True) #inplace = True updates the array instead of creating a copy
    #Assign a CRS to the xarray that matches the shapefile used for clipping. CRS included is CF compliant.
    array.rio.write_crs(shp.crs, inplace = True) #inplace = True updates the array instead of creating a copy
    
    #Clipping maintains only those pixels whose center is within the polygon boundaries and drops any data not meeting this requirement.
    clipped = array.rio.clip(shp.geometry, shp.crs, drop = True, invert = False, all_touched = False)
    
    return clipped

########
#Calculates weighted means by season, by month or per timestep
def weightedMeans(array, weights, meanby = 'timestep'):
    '''
    Inputs:
    array - Data array containing variable from which means will be calculated
    weights - Data array containing weights
    meanby - Define how means will be calculate: timestep, month or season. Default set to 'timestep'
    
    Output:
    Data array containing weighted means
    '''
      
    #Calculate weights
    weights = weights/weights.sum()
        
    #Apply weights to variable - Calculate weighted mean over timestep and then calculate the mean by season
    if meanby == 'season':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean')).groupby('time.season').mean()
    elif meanby == 'month':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean')).groupby('time.month').mean()
    elif meanby == 'timestep':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean'))
        
    return weighted_mean

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
#This function performs a linear trend calculation and returns the coefficients as well as p-values for the linear regression
def linearTrends(y, x, rsquared = False):
    '''
    Inputs:
    y - data array with information about dependent variable
    x - data array with information about independent variable
    rsquared - Boolean. If set to True then r squared values will be calculated and returned as outputs
        
    Output:
    Coefficients and p-values of linear regression
    '''
    #To check extra information available in the model use
        #dir(model.fit())
        
    if rsquared == True:
        model = sm.OLS(y, x)
        coef = model.fit().params[1]
        sig = model.fit().pvalues[1]
        rsq_adj = model.fit().rsquared_adj
        return coef, sig, rsq_adj
    else:
        model = sm.OLS(y, x)
        coef = model.fit().params[1]
        sig = model.fit().pvalues[1]
        return coef, sig

    
########
#This function performs a linear trend calculation and returns the coefficients as well as p-values for the linear regression
def lm_yr(y, x):
    '''
    Inputs:
    y - data array with information about dependent variable
    x - data array with information about independent variable
            
    Output:
    Coefficients and p-values of linear regression
    '''
    return ss.linregress(x, y)
    
########
#Defining function that will be applied across dimensions
def lm_lats(arr, lats):
    '''
    Inputs:
    arr - data array containing the dependent and independent variables
    lats - list containing the latitudes for which we will calculate linear regressions
        
    Output:
    Dataset containing the slope, intercept, p and r squared values, std error and predictions
    for the latitudes of interest
    '''
    #Create empty list to store results of linear regression
    slope = []; intercept = []; r_val = []; p_val = []; stderr = []; pred = []
    #Extract values for each value of interest
    for lat in lats:
        sub_lat = arr.sel(yt_ocean = lat, method = 'nearest').dropna('time')
        #Calculate linear regression
        try:
            r_lm = ss.linregress(sub_lat.time.dt.year, sub_lat.values)
        #If a particular latitude only has NA values, create an exception
        #so NA values will be returned for all linear regression results
        except:
            sub_lat = arr.sel(yt_ocean = lat, method = 'nearest')
            r_lm = ss.linregress(sub_lat.time.dt.year, sub_lat.values)
        #Add results to empty lists
        slope.append(r_lm.slope)
        intercept.append(r_lm.intercept)
        r_val.append(r_lm.rvalue)
        p_val.append(r_lm.pvalue)
        stderr.append(r_lm.stderr)
        pred.append(r_lm.intercept+(r_lm.slope*sub_lat.time.dt.year))
    
    #Create a data array with predictions 
    pred = xr.concat(pred, dim = 'yt_ocean')
    #and the rest of variables
    slope = xr.DataArray(data = slope, name = 'slope', dims = ['yt_ocean'], 
                 coords = dict(yt_ocean = arr.yt_ocean.values))
    intercept = xr.DataArray(data = intercept, name = 'intercept', dims = ['yt_ocean'], 
                 coords = dict(yt_ocean = arr.yt_ocean.values))
    p_val = xr.DataArray(data = p_val, name = 'p_val', dims = ['yt_ocean'], 
                 coords = dict(yt_ocean = arr.yt_ocean.values))
    r_val = xr.DataArray(data = r_val, name = 'r_val', dims = ['yt_ocean'], 
                 coords = dict(yt_ocean = arr.yt_ocean.values))
    stderr = xr.DataArray(data = stderr, name = 'stderr', dims = ['yt_ocean'], 
                 coords = dict(yt_ocean = arr.yt_ocean.values))
    
    #Change names prior to creating final dataset
    pred.name = 'predictions'
    arr.name = 'model_data'
    
    #Merge everything into one dataset
    ds = xr.merge([arr, pred, slope, intercept, r_val, p_val, stderr])
    
    return ds    
    
########
#This function calculates anomalies 
def AnomCalc(array, clim_array, std_anom = False):
    '''
    Inputs:
    array - refers to a data array containing information for the period being compared to the baseline. It could include just one year or multiple years (decades)
    clim_array - three dimensional array containing data over the baseline period from which anomalies will be calculated
    std_anom - boolean variable that if set to True will result in standarised anomalies being calculated
      
    Outputs:
    Data array containing anomalies.
    '''
    
    #Calculate long term mean of array
    m_climarray = clim_array.mean('time')
      
    #Calculate anomalies
    #Standarised anomalies
    if std_anom == True:
        s_climarray = clim_array.std('time')
        anom = (array - m_climarray)/s_climarray
    #Non-standarised anomalies
    elif std_anom == False:
        anom = array - m_climarray
    
    #Return anomalies
    return anom


########
#Getting maximum and minimum years included in the analysis
def colbarRange(dict_data, sector, season):
    '''
    Inputs:
    dict_data - refers to a dictionary that contains the datasets being plotted.
    sector - refers to a list of sectors
    season - refers to a list of seasons
    '''
    
    #Define variables that will store the maximum and minimum values
    maxV = []
    minV = []
    
    sector = sector*len(season)
    season = np.concatenate([[i]*len(sector) for i in season], axis = 0)
    for sec, sea in zip(sector, season):
        maxV.append(dict_data[f'{sec}_{sea}'].indexes['time'].year.max())
        minV.append(dict_data[f'{sec}_{sea}'].indexes['time'].year.min())
        
    maxV = max(maxV)
    minV = min(minV)

    return minV, maxV



########
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
#Calculate the lat-lon coordinates from a dataset in source_crs - Function by Scott Wales
def reproject_latlon_coords(da, source_crs, target_crs):
    '''
    Inputs:
    da - refers to a data array that needs to be reprojected. Dimensions containing spatial data should be labelled as x and y.
    source_crs - original CRS for data array spatial information. Should be provided as a string in the form of 'epsg:4326'.
    target_crs - CRS to which the data array will be transformed. Should be provided as a string in the form of 'epsg:4326'.
    '''
    
    # Convert the 1d coordinates to 2d arrays covering the whole grid
    X, Y = np.meshgrid(da.x, da.y)
       
    # Convert the 2d coordinates from the source to the target values
    lon, lat = transform(Proj(init = source_crs), Proj(init = target_crs), X, Y)
    
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
#Calculating climatology
def climCalc(da, clim_period, varname, clim_type = 'overall', **kwargs):
    '''  
    Inputs:
    da - data array, containing information from which a climatology will be calculated
    clim_period - list, time frame to be used in climatology calculation. Only start and end year are needed.
    clim_type - str, what type of climatology needs to be calculated. Default is 'overall', also available
        'seasonal' and 'monthly'.
    Optional:
    folder_out - str, containing file path to folder where results will be stored.
   
    Outputs:
    clim - data array containing the calculated climatology
    '''
    
    #Ensuring time period is a string
    clim_period = [str(yr) for yr in clim_period]
    
    #Select data within climatology period
    clim = da.sel(time = slice(*clim_period))

    #Calculating climatology
    if clim_type == 'overall':
        clim = clim.mean('time')
        fn = varname + f'_Climatology_overall_{clim_period[0]}-{clim_period[1]}.nc'
    elif clim_type == 'seasonal':
        clim = clim.groupby('time.season').mean('time')
        fn = varname + f'_Climatology_seasonal_{clim_period[0]}-{clim_period[1]}.nc'
    elif clim_type == 'monthly':
        clim = clim.groupby('time.month').mean('time')
        fn = varname + f'_Climatology_monthly_{clim_period[0]}-{clim_period[1]}.nc'
 
    if 'folder_out' in kwargs.keys():
        #Define output folder, where climatologies will be stored
        out_folder = kwargs.get('folder_out')
        #Ensure output folder exists
        os.makedirs(out_folder, exist_ok = True)
        #Saving results to disk
        clim.to_netcdf(os.path.join(out_folder, fn))
        
    return clim

########
def main(inargs):
    '''Run the program.'''

if __name__ == '__main__':
    description = 'This script contains functions used to perform timeseries within different sectors of the Southern Ocean.'
    parser = argparse.ArgumentParser(description = description)

    args = parser.parse_args()
    main(args)