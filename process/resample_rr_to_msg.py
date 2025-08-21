"""
This code has the goal to resample the radar data from DWD ( 5 min rainfall in kgm-2 provided at 5 min temporal resolution and 1kmx1km grid (see link: https://opendata.dwd.de/climate_environment/CDC/help/landing_pages/doi_landingpage_RADKLIM_RW_V2017.002-en.html) to the msg satellite temporal and spatial resolution (15 min and 0.04 degrees )
It then stores the ncdf data on the expats-radar-germany data bucket on EWC from EUMETSAT with filenames of the form yyyymmdd_RR_15min_msg_res.nc.
The code is done to process data from 2023 to 2013

author: Claudia Acquistapace
date: 19 Ago 2025

"""
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

from figures.plotting import plot_msg, plot_radar_mapV2
from readers.data_buckets_funcs import upload_to_bucket, check_file_bucket
from process.utils import generate_regular_grid, regrid_data, is_valid_date

def main():

    # download_from_s3()
    # loop on all days of
    # set the day to process
    
    # list all years from 2024 to 2014
    yy_list = [str(year) for year in range(2023, 2021, -1)]
    # list all months in a year
    mm_list = [f'{month:02d}' for month in range(4, 10)]
    # list all days in a month
    dd_list = [f'{day:02d}' for day in range(1, 32)]
    print(yy_list, mm_list, dd_list)
    
    # loop on all years, months, days to find yy mm dd for processing
    for yy in yy_list:
        for mm in mm_list:
            for dd in dd_list:
                
                # print the date to process
                print(yy, mm, dd)
                
                # process only if the date is valid
                if is_valid_date(yy, mm, dd):
                    
                    # check if the file is already on the bucket
                    file_exist = check_file_bucket(yy+mm+dd+'_RR_15min_msg_res.nc.gz')
                    
                    if file_exist:
                        print('file already exists on bucket - skipping date:', yy, mm, dd)
                        continue
                    else: 
                        
                        print('generating ncdf file for ', yy+mm+dd)
                        
                        # set keyword to decide if plotting or not
                        plotting = False
                        
                        # set path out
                        path_out = '/net/ostro/radolan_5min_composites/'
                    
                        # check that radar file from DWD exists 
                        if os.path.exists('/data/trade_pc/radolan_DE_5min_rain_rate/'+yy+'/'+mm+'/YW_2017.002_'+yy+mm+dd+'.nc'):
                            
                            # run processing to create interpolated radolan radar data on msg time/space res. returns filename created
                            ncdf_radar_gz = process_date(yy, mm, dd, plotting, path_out)

                            # call function to transfer data to bucket, returns true if transfer is complete
                            check = upload_to_bucket(path_out, ncdf_radar_gz)
                            
                        else:
                            # if file radolan does not exist, print message and continue to next day
                            print('radolan file missing - go to next day')
                            continue
                        
                        # check if file was uploaded
                        if check:
                            print(f'{ncdf_radar_gz} file uploaded to bucket')
                            
                            # remove ncdf file produced
                            # cut .gz at the end of ncdf_radar
                            ncdf_radar_name = ncdf_radar_gz[:-3]
                            
                            # if file ncdf_radar+name exists, remove it
                            if os.path.exists(path_out+ncdf_radar_name):
                                os.remove(path_out+ncdf_radar_name)

                        else:
                            print('file not uploaded to bucket - abort')
                            break





def process_date(yy, mm, dd, plotting, path_out):
    '''
    function to read rain rate from dwd data, resample it on the msg temporal (15 min) and spatial (0.04 degrees) resolution, and then store a new netcdf on the small radar domain. It also, if requested by the plotting keyword, produces plots of the rr instantaneous, cumulated over 15 mins, and msg 
    inputs:
    - yy: string of the year
    - mm: string of the month
    - dd: string of the day
    - plotting: boolean, if True produces plots
    - path_out: string, path to output directory where to store ncdf and plots
    returns
    - filename of the file ncdf produced .nc.gz
    Dependencies:
    - read_radar_DWD
    - generate_regular_grid
    - regrid_data
    - read_orography
    - plot_radar_mapV2
    - plot_msg

    '''
    # building date 
    date = yy+mm+dd
    print('processing date to resample RR', date)
    
    # get start time
    st = time.time()

    # read radar data
    ds_radar = read_radar_DWD(path_radolan_DE, date)   

    # read one msg file to get temporal and spatial resolution
    try:
        ds_msg = xr.open_dataset('/data/sat/msg/netcdf/parallax/'+yy+'/'+mm+'/'+date+'-EXPATS-RG.nc')
    except FileNotFoundError:
        print('MSG file not found - opening standard one')
        # abort plotting 
        plotting = False
        
        # reading a random existing file
        ds_msg = xr.open_dataset('/data/sat/msg/netcdf/parallax/2022/08/19/20220819-EXPATS-RG.nc')

    # resample ds_RR_msg to the msg time resolution of 15 minutes. RR values are summer over the time interval to get rain amount over 15 minutes
    ds_radar_msg = ds_radar.resample(time='15T').sum(skipna=True)
    
    # crop msg data to domain_DE_CA
    lon_min, lon_max, lat_min, lat_max = domain_DE_CA

    #print(f'Cropping MSG data to domain: {lat_min}, {lat_max}, {lon_min}, {lon_max}')
    
    # crop msg data to the DE domain (smaller than EXPATS)
    ds_msg_DE = ds_msg.where((ds_msg.lat >= lat_min) & 
                            (ds_msg.lat <= lat_max) & 
                            (ds_msg.lon >= lon_min) & 
                            (ds_msg.lon <= lon_max), 
                            drop=True)

    # if RR file for the day does not exist, read, resample and store
    if not os.path.exists(os.path.join(path_out, date+'_RR_15min_msg_res.nc.gz')):
        
        print('file RR at msg res not existing - resampling now')
        
        # define MSG spatial resolution for the regular grid to be generated
        step_deg = 0.04  # step size in degrees for the regular grid

        # ceate a regular grid for the msg data
        lat_arr, lon_arr = generate_regular_grid(lat_min,
                                                lat_max,  
                                                lon_min, 
                                                lon_max, 
                                                step_deg,
                                                path=None)
        
        # create 2d arrays of lat lon based on the regular lat lon arrays
        lat_reg_grid, lon_reg_grid = np.meshgrid(lat_arr, lon_arr, indexing='ij')

        # define RR_matrix where to store RR resampled
        RR_msg = np.zeros((len(ds_radar_msg.time.values), len(lat_arr), len(lon_arr)))
        
        time1 = time.time()
        
        print(f'Time taken for resampling to msg: {time1 - st:.2f} seconds')
        
        # loop on time stamps to resample 15 min rain amounts to MSG spatial resolution
        with Bar('Processing...') as bar: 
            for i_t, time_val in enumerate(ds_radar_msg.time.values):
                
                print('regrid RR for time ', time_val)

                # setting time obs for loop start
                time_loop_start = time.time()

                # regrid rain rate date from old radar grid to the new msg grid
                RR_msg[i_t, :, :] = regrid_data(ds_radar_msg.lat.values,
                                    ds_radar_msg.lon.values,
                                    ds_radar_msg.RR.values[i_t,:,:],
                                    lat_reg_grid,
                                    lon_reg_grid)
                
                time_loop_end = time.time()
                print(f'Time taken for loop iteration {i_t}: {time_loop_end - time_loop_start:.2f} seconds')
                bar.next()


        # store RR_msg into a xarray dataset
        ds_RR_msg = xr.Dataset(
            {
                "RR": (("time", "lat", "lon"), RR_msg[:, :, :],
                       {
                        "description": "Rain rate data resampled to MSG grid using temporal summation",
                        "long_name": "15 minutes rainfall",
                        "standard_name": "rainfall_amount",
                        "units": "kgm-2",
                        "processing_method": "temporal resampling with nansum aggregation",
                        "temporal_resolution": "15 minutes",
                        "spatial_resolution": "0.04 degrees",
                        "valid_min": 0.0,
                        "valid_max": 1000.0,
                        "_FillValue": np.nan
                       })
            },
            coords={
                "time": ds_radar_msg.time.values,
                "lat": lat_arr, 
                "lon": lon_arr,
            },
            attrs={
                "description": "Rain rate from DWD data resampled to MSG grid as cumulated over 15 minutes",
                "history": "https://opendata.dwd.de/climate_environment/CDC/help/landing_pages/doi_landingpage_RADKLIM_RW_V2017.002-en.html",
                "source": "DWD Radar",
                "reference history":"10.5676/DWD/RADKLIM_YW_V2017.002",
                "created_by": "Claudia Acquistapace",
                "created_on": str(pd.Timestamp.now()),
                "domain": "DE_CA",
                "original_grid": "RADOLAN",
                "target_grid": "MSG-like regular grid"
            }
        )
        # store and compress RR variable to maximum compression
        ds_RR_msg.to_netcdf(os.path.join(path_out, date+'_RR_15min_msg_res.nc'), encoding={'RR': {'zlib': True, 'complevel': 9}})
        
        # gzip the nc file
        with open(os.path.join(path_out, date+'_RR_15min_msg_res.nc'), 'rb') as f_in:
            with gzip.open(os.path.join(path_out, date+'_RR_15min_msg_res.nc.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                
        # remove the local nc file
        os.remove(os.path.join(path_out, date+'_RR_15min_msg_res.nc'))
        

    else:
        
        filename_existing = os.path.join(path_out, date+'_RR_15min_msg_res.nc.gz')
        print(f"{filename_existing} exists, don't do anything")
        # if ncdf file exists
        #print('file RR at msg res exists - reading it ')
        # read stored ncdf
        #ds_RR_msg = xr.open_dataset(os.path.join(path_out, date#+'_RR_15min_msg_res.nc'))    

    # PLOT
    if plotting:
        
        lats = ds_msg_DE.lat.values
        lons = ds_msg_DE.lon.values
        cmap = 'Greys'
        cbar_title = 'IR_10.8 micron'
        domain = domain_DE_CA
        key = 'expats'
        back_transparent = 'True'
        vmin = 200
        vmax = 310
        
        # plot RR maps for comparison at 15 min resolution 
        with Bar('Processing...') as bar: 
            for i_t, time_val in enumerate(ds_RR_msg.time.values):
                
                # generate string of the date with hh and minute
                time_string = pd.Timestamp(time_val).strftime('%Y%m%d_%H%M')
                
                print(f'Plotting radar map for time {time_val} at original and msg resolution')
                label = f'{time_string}_MSG_108'

                plot_radar_mapV2(ds_RR_msg.isel(time=i_t), f'{time_string}_RR_msg', path_out)

                # plot radar data at original resolution on the closest time stamp to compare
                plot_radar_mapV2(ds_radar.sel(time=time_val, method='nearest'), f'{time_string}_RR_orig', path_out)

                print(f'Plotting MSG map for time {str(time_val)}')
                variable = ds_msg_DE.IR_108.values[i_t, :, :]
                
                plot_msg(lons, lats, variable, label, cmap, cbar_title, domain, path_out, key, back_transparent, vmin, vmax) 

                bar.next()

    return(date+'_RR_15min_msg_res.nc.gz')



if __name__ == "__main__":
    main()