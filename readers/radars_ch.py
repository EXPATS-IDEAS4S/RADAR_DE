



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
# import crameri colormaps
import cmcrameri.cm as cmr

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
    note: see bottom of the function fo more info on the files and their naming convention
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
                    'time_resolution': f"{time_step} seconds",
                    'comment': f"time contains the starting time of the {time_step} interval",
                    'author': 'Claudia Acquistapace',
                    'email': 'claudia.acquistapace-at-unipd.it',
                    'date': datetime.now().strftime("%Y-%m-%d"),
                    'data source': 'Swiss radar data, read from h5 files with the function read_h5_file in readers/radars_ch.py'
                }

            )   

            # append results to the arrays of the day        
            time_arr.append(start_datetime)
            time_resolution.append(time_step)
            RR_arr.append(ds)

    # concatenate data of the day along the time dimension
    ds_day = xr.concat(RR_arr, dim='time')
    # sort ds_day by time
    ds_day = ds_day.sortby('time')

    # add nan masking to define radar domain
    domain_mask_array = np.where(ds_day.RR.values[0,:,:] >= 0, 1, 0)

    # calibration factor to add for calibrating RR with rain gauges (summer period, personal comm Francesco Marra)
    ds_day['RR'] = ds_day['RR'] * 0.87

    # set to nan all RR values equal to 0, to avoid plotting them
    ds_day['RR'] = ds_day['RR'].where(ds_day['RR'] != 0., np.nan)

    # add rr units as attribute
    ds_day['RR'].attrs['units'] = 'mm/h'
    ds_day['RR'].attrs['description'] = "istantaneous rain rate in mm/h, with values equal to 0 set to nan."

    # add comment on lat lon resolution
    ds_day['lat'].attrs['comment'] = "latitude values of the radar grid, reconstructed from the lat lon information in the h5 file. "
    ds_day['lon'].attrs['comment'] = "longitude values of the radar grid, reconstructed from the lat lon information in the h5 file. "
 
    # add domain mask to the dataset as a variable
    ds_day = ds_day.assign(domain_mask=(('lat', 'lon'), domain_mask_array)) 

    # add description of domain mask array as attributes
    ds_day.domain_mask.attrs['description'] = "mask to define radar domain, with 1 for grid points \n" \
    "inside the radar domain and 0 for grid points outside the radar domain. The mask is defined based \n" \
    "on the first time step of the day, assuming that the radar domain does not change during the day. \n" \
    "The mask is created by setting to 1 all grid points with RR values greater than or equal to 0,\n" \
    " and to 0 all grid points with RR values less than 0. "

    # add comment with the name of the script used to produce the ds_day dataset
    ds_day.attrs['comment'] = "dataset created with the function read_h5_file in readers/radars_ch.py"
    ds_day.attrs['author'] = 'Claudia Acquistapace'
    ds_day.attrs['email'] = 'claudia.acquistapace-at-unipd.it'
    ds_day.attrs['date'] = datetime.now().strftime("%Y-%m-%d")
    ds_day.attrs['Github repository'] = 'https://github.com/claudia-acquistapace/RADAR_DE'
    ds_day.attrs['data source'] = 'Swiss radar data, read from h5 files with the function read_h5_file in readers/radars_ch.py'
    


    """
    info on files:

    From the size of the file, you can get
      an idea whether it was clear sky 
      (minimum size) or with (widespread) 
      precipitation.
    The data set covers all day;
     The domain of (each) raster file 
     is 710 columns by 640 rows, 
     hich is 454400 pixels. 
     Only ~350400 pixels are 
     under the umbrella of the 
     5 radars, the others are 
     out of range (NaN).
    
    Regarding the conformal grid:
    Bottom-left corner has kilometric
      coordinate  (255 East, -180 North).
    Top-right corner has kilometric 
    coordinate  (965 East, +480 North).
     
    By way of example, the coordinate
      of Bern, the Swiss Capital 
      is 600 East (“Y”), 200 North (“X”).
    Inside Switzerland, such coordinates 
    are always positive; furthermore, Y>X.

    Files with the extension 801 are coded
      on 8 bits, i.e. there are at most 
      256 different values (classes), 
      even though in the HDF5 file 
      they appear as floating point numbers. 
      The 001 files, with the name RZCflt,
        are 32 bits floating point data 
        throughout the processing, i.e.
          they should be slightly more 
          precise as there is no rescaling 
          to 8bit/256 classes. Although
            from a quantitative point
              of view there difference 
              shall be minimal, when you 
              have both options, I would
                use the 001 files.

    files naming convention:
    RZC16001
    yyddd
    year 2016 
    day 1 i.e. jan 1 

    RZC22162
    year 2022 
    day 162 i.e. jun 11 
    """
    return ds_day
