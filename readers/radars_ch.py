



from datetime import date, datetime
import h5py
from matplotlib.pyplot import contourf
import xarray as xr
import pdb
import numpy as np
from pyproj import CRS, Transformer
from readers.radar_DWD import read_orography
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

def read_h5_file(file_path, date):
    """
    Function to read an h5 file and store its content in a xarray dataset. The function reads all 
    h5 files of the day and collects data, time stamp and reconstruct lat lon. 
    All data are then stored in a xarray dataset (time, lat, lon) and returned as output.   
    input: 
    file_path: list of h5files of the day to read
    date: date of the data  
    output: 
    xarray dataset with the content of the h5 file

    author: Claudia Acquistapace
    date: 16 Febbruary 2026
    email: claudia.acquistapace-at-unipd.it
    """

    # loop on files of the day and read data, time stamp and lat lon. store in xarray dataset

    time_arr = []
    time_resolution = []
    RR_arr = []

    # loop on files of the day
    for file_to_read in file_path:

        # open the h5 file and read data, time stamp and lat lon
        with h5py.File(file_to_read, 'r') as f:

            # reading rr in mm/h from the dataset. The path to the data is: f['dataset1']['data1']['data'][:]
            rr = f['dataset1']['data1']['data'][:]

            # reading starting, ending time and time step of the data. The path to the time information is: f['dataset1']['data1']['time'][:]
            time_info = f['dataset1/what']
            time_attrs = dict(time_info.attrs)
            start_time = time_attrs['starttime'][:]
            end_time = time_attrs['endtime'][:]

            # assuming times can be read at hhmmss
            start_time_obj = datetime.strptime(str(start_time)[2:-1], "%H%M%S").time()
            end_time_obj = datetime.strptime(str(end_time)[2:-1], "%H%M%S").time()

            # calculating time resolution 
            time_step = abs(start_time_obj.hour*3600 + start_time_obj.minute*60 + start_time_obj.second - (end_time_obj.hour*3600 + end_time_obj.minute*60 + end_time_obj.second))

            # build timestamp from the file using stard time obj by adding yy mm and day from the date in input
            start_date_obj = datetime.strptime(str(date), "%Y-%m-%d").date()
            start_datetime = datetime.combine(start_date_obj, start_time_obj)

            # read lat and lon from the dataset. The path to the lat and lon information is: f['dataset1']['data1']['lat'][:] and f['dataset1']['data1']['lon'][:]
            space_info = f['where']
            space_attrs = dict(space_info.attrs)

            # reading parameters to construct lat and lon grid
            LL_lat = space_attrs['LL_lat']
            LL_lon = space_attrs['LL_lon']
            UL_lat = space_attrs['UL_lat']
            UL_lon = space_attrs['UL_lon']
            LR_lat = space_attrs['LR_lat']
            LR_lon = space_attrs['LR_lon']
            UR_lat = space_attrs['UR_lat']
            UR_lon = space_attrs['UR_lon']

            LL = (LL_lat, LL_lon)
            UL = (UL_lat, UL_lon)
            LR = (LR_lat, LR_lon)
            UR = (UR_lat, UR_lon)

            xsize = space_attrs['xsize']
            ysize = space_attrs['ysize']

            # Create normalized grid coordinates
            xi = np.linspace(0, 1, xsize)
            yi = np.linspace(0, 1, ysize)
            xx, yy = np.meshgrid(xi, yi)

            # Bilinear interpolation for lat
            lat = (1 - xx) * (1 - yy) * LL[0] + xx * (1 - yy) * LR[0] + (1 - xx) * yy * UL[0] + xx * yy * UR[0]
            # Bilinear interpolation for lon
            lon = (1 - xx) * (1 - yy) * LL[1] + xx * (1 - yy) * LR[1] + (1 - xx) * yy * UL[1] + xx * yy * UR[1]

            # store data in xarray dataset
            ds = xr.Dataset(
                data_vars={
                    'RR': (('time', 'lat', 'lon'), rr[np.newaxis, :, :])

                },
                coords={
                    'time': [start_datetime],
                    'lat': lat[:, 0],
                    'lon': lon[0, :]
                },
                attrs={
                    'time_resolution': time_step,
                    'comment': f"time contains the starting time of the {time_step} interval",
                    'author': 'Claudia Acquistapace',
                    'email': 'claudia.acquistapace-at-unipd.it',
                    'date': '16 Febbruary 2026',
                    'data source': 'Swiss radar data, read from h5 files with the function read_h5_file in readers/radars_ch.py'
                }

            )   

            # append results to the arrays of the day        
            time_arr.append(start_datetime)
            time_resolution.append(time_step)
            RR_arr.append(ds)

    # concatenate data of the day along the time dimension
    ds_day = xr.concat(RR_arr, dim='time')

    # add nan masking to define radar domain
    domain_mask_array = np.where(ds_day.RR.values[0,:,:] >= 0, 1, 0)

    # add domain mask to the dataset as a variable
    ds_day = ds_day.assign(domain_mask=(('lat', 'lon'), domain_mask_array)) 

    # add description of domain mask array as attributesx   

    return ds_day
