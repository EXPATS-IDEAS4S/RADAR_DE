from datetime import datetime
import glob
from scipy.interpolate import griddata
import xarray as xr
from readers.file_dirs import path_radolan_DE
from readers.radar_DWD import read_radar_DWD, read_orography
import pdb
import numpy as np
import scipy
from readers.data_buckets_funcs import download_from_s3
from figures.domain_info import domain_expats, domain_DE_CA
from figures.mpl_style import plot_cities_expats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from figures.one_day_video import plot_radar_map
import matplotlib.pyplot as plt
import os
from time import sleep
from progress.bar import Bar
import time
import pandas as pd
import shutil
import gzip




def generate_regular_grid(lat_min, lat_max, lon_min, lon_max, step_deg, path=None):
    """
    Generate a regular grid for a given bounding box and step size in degrees.

    :param lat_min: Minimum latitude of the bounding box.
    :param lat_max: Maximum latitude of the bounding box.
    :param lon_min: Minimum longitude of the bounding box.
    :param lon_max: Maximum longitude of the bounding box.
    :param step_deg: The step size in degrees for both latitude and longitude.
    :param path: If given save the lat and lon regular points there
    :return: A 1D array of regular lat and lon points.
    """
    lat_points = np.arange(lat_min, lat_max + step_deg, step_deg)
    lon_points = np.arange(lon_min, lon_max + step_deg, step_deg)

    if path:
        np.save(path+'reg_lats.npy', lat_points)
        np.save(path+'reg_lons.npy', lon_points)

    return lat_points, lon_points



def regrid_data(old_lat, old_lon, old_data, new_lat, new_lon, method='linear'):
    """
    Regrid data from an old grid to a new grid.

    :param old_lat: 2D array of latitudes for the old grid.
    :param old_lon: 2D array of longitudes for the old grid.
    :param old_data: 2D array of data corresponding to the old grid.
    :param new_lat: 2D array of latitudes for the new grid.
    :param new_lon: 2D array of longitudes for the new grid.
    :param method: Interpolation method ('linear', 'nearest', 'cubic').
    :return: 2D array of regridded data corresponding to the new grid.
    """
    # Flatten the old grid coordinates and data for interpolation
    old_coords = np.array([old_lat.ravel(), old_lon.ravel()]).T
    old_data_flat = old_data.ravel()

    # Create a mesh of new grid coordinates
    new_coords = np.array([new_lat.ravel(), new_lon.ravel()]).T

    # Interpolate old data to new grid
    new_data_flat = griddata(old_coords, old_data_flat, new_coords, method=method)

    # Reshape the flattened data back into the 2D structure of the new grid
    new_data = new_data_flat.reshape(new_lat.shape)

    return new_data

  
    
def is_valid_date(year, month, day):
    ''' function to find out if the date built with year, month, day exists
    inputs:
    year: string of the year
    month: string of the month
    day: string of the day
    
    returns:
    - true if the date exists
    - false if the date does not exist
    '''
    import datetime
    
    # Convert string inputs to integers if necessary
    if isinstance(year, str):
        year = int(year)
    if isinstance(month, str):
        month = int(month)
    if isinstance(day, str):
        day = int(day)
    
    try:
        datetime.date(year, month, day)
        return True
    except ValueError:
        return False
    
def extract_date_from_filename_ch(filename):
    """
    function to extract date from the filename of swiss radar data 
    with the following rule:
      last 5 characters before .zip are yyddd,
        where yy is the year and 
        ddd is the day of the year  
    input:
        filename: string, name of the file to extract the date from
    output:  
        date: string, date extracted from the filename in the format YYYY-MM-DD
    
    author: Claudia Acquistapace
    date: 16 Febbruary 2026
    email: claudia.acquistapace-at-unipd.it
    
    """
    # extract date from the filename with the following rule: last 5 characters before .zip are yyddd, where yy is the year and ddd is the day of the year
    yy = '20'+ filename[-5:-3]
    # read day of the year and convert it to month and day
    ddd = int(filename[-3:])

    # convert day of the year to month and day including months with 31 days, 30 days and February with 28 or 29 days depending on the year (leap year or not)

    date = datetime.strptime(f"{yy}-{ddd:03d}", "%Y-%j").date()
    month = date.month
    day = date.day
    #print(f"date extracted from filename: {yy}-{month:02d}-{day:02d}")
    return f"{yy}-{month:02d}-{day:02d}"


def select_months_of_interest(files_list, file_path, months):
    """function to select only files from the months of interest, since \n
    the radar data are not good in winter
    """
    dates = []
    for file in files_list:
        filename = file.split("/")[-1].split(".zip")[0]
        date = extract_date_from_filename_ch(filename)
        dates.append(date)
    dates = np.array(dates)
    files_list = np.array(files_list)

    # select only files from the months of interest
    selected_files = []
    selected_dates = []
    for month in months:
        month_files = files_list[np.array([date.split("-")[1] == month for date in dates])]
        selected_files.extend(month_files)
        selected_dates.extend([date for date in dates if date.split("-")[1] == month])

    print(selected_dates)
    return selected_files   