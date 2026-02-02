"""
This code reads and plots radar data from ARPAE 

how to run this code:
python3 -m process.main_arpae


"""
import glob
from scipy.interpolate import griddata
import xarray as xr
import pdb
import numpy as np
import scipy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import sleep
from progress.bar import Bar
import time
import pandas as pd
import shutil
import gzip
import os
import xarray as xr
from datetime import datetime
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import boto3
import io
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from readers.file_dirs import path_radolan_DE
from readers.radar_DWD import read_radar_DWD, read_orography
from readers.config import *
from readers.data_buckets_funcs import upload_to_bucket, check_file_bucket
from readers.config import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from figures.domain_info import domain_expats, domain_IT_CA
from figures.one_day_video import plot_radar_map
from process.utils import generate_regular_grid, regrid_data, is_valid_date
from readers.radar_DWD import read_orography
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg
from readers.file_dirs import orography_file
from readers.data_buckets_funcs import read_file_obj, init_s3, upload_file
from readers.data_buckets_funcs import download_from_s3
from figures.one_day_video import plot_radar_map

def ze_to_rr(ze, a=200, b=1.6):
    """
    Convert radar reflectivity Ze (dBZ) to rain rate RR (mm/h).
    ze: float or np.ndarray, reflectivity in dBZ
    a, b: Z-R relationship constants (default: Marshall-Palmer)
    """
    Z = 10 ** (ze / 10)
    rr = (Z / a) ** (1 / b)
    return rr


def print_list_files_in_bucket(s3, S3_BUCKET_NAME):
    """Lists all files in the specified S3 bucket.
    input:
        s3: S3 client
        bucket: Name of the S3 bucket
    output:
        List of file keys in the bucket.
    """
    # Pagination to get all objects
    paginator = s3.get_paginator('list_objects_v2')
    pages = paginator.paginate(Bucket=S3_BUCKET_NAME)
    all_objects = []
    for page in pages:
        if "Contents" in page:
            #for obj in page['Contents']:
                #print(obj['Key'])  # Only print the Key
            all_objects.extend(page["Contents"])
    return all_objects


def read_all_filenames_from_bucket():
    """
    function to read all filenames from the specified S3 bucket and store them in a npy file if the npy file does not exist yet
    otherwise it loads the npy file
    Returns:
        file_names: list of filenames found on the bucket
    Dependencies:
        readers.config
        readers.data_buckets_funcs.init_s3
        
    """
    # start logger and initialize s3 client
    s3 = init_s3()
    
    # check if the npy file already exists
    if os.path.exists('file_names_ARPAE.npy'):
        file_names = np.load('file_names_ARPAE.npy', allow_pickle=True)
        return s3, file_names
    else:
        # Pagination to get all objects
        paginator = s3.get_paginator('list_objects_v2')
        pages = paginator.paginate(Bucket=S3_BUCKET_NAME)
        file_names = []
        for page in pages:
            if "Contents" in page:
                for obj in page['Contents']:
                    #print(obj['Key'])  # Only print the Key
                    file_names.append(obj['Key'])
        
        # store file_names in a npy file
        np.save('file_names_ARPAE.npy', file_names)

        return s3, file_names


def calc_acc_rain_over_deltaT(Ze_5_min, deltaT=5):
    """
    function to calculate accumulated rain over a deltaT intervals from Ze data
    input:
        Ze_5_min: 3D array of Ze data with shape (time, lat, lon)
        deltaT: time interval in minutes over which to calculate accumulated rain (default is 5 min)
    output:
        RR_acc_deltaT: 3D array of accumulated rain over deltaT intervals with shape (time, lat, lon)

    The function first converts Ze from dbZ to mm6/mm3, then applies the Z-R relationship to convert Ze to RR, which is then the rain 
    rate that would occur if conditions observed for the 5 min period would persist for 1 hour.
    Finally, it converts the rain rate in mm/h to accumulated rainfall over the given time step, which is 5 min
    5 min = 5/60 h
    15 min = 15/60 h
    30 min = 30/60 h
    Therefore, RR_acc_5_min = RR_5_min * (5/60)
    1 mm/h = 1 kg/m2/h
    Thus, the output RR_acc_5_min is in mm or kg/m2 over the 5 min interval.
    """

    deltaT_hours = deltaT / 60  # convert deltaT from minutes to hours

    # convert Ze from dBZ to mm6/mm3
    Z_5min_lin = 10 ** (Ze_5_min / 10)

    # apply Z-R relationship to convert Ze to RR in mm/h
    RR_5_min = (Z_5min_lin / 200) ** (1 / 1.6)  # using Marshall-Palmer relationship

    # convert RR in mm/h to accumulated rainfall over 5 min
    RR_acc_deltaT = RR_5_min * deltaT_hours  # converting from mm/h to mm over 5 min

    # set to zero all negative values
    RR_acc_deltaT[RR_acc_deltaT < 0] = np.nan

    return RR_acc_deltaT


def main():

    # output path for processed files
    s3_bucket_new_name = "expats-radar-italy" # name of the destination data bucket where to store daily processed files of RR from italian domain

    # read all filenames from bucket
    s3, file_names = read_all_filenames_from_bucket()

    # list all years from 2024 to 2014
    yy_list = [str(year) for year in range(2017, 2016, -1)]
    # list all months in a year
    mm_list = [f'{month:02d}' for month in range(4, 10)]
    # list all days in a month
    dd_list = [f'{day:02d}' for day in range(1, 32)]
    print(yy_list, mm_list, dd_list)
    
    # loop on all years, months, days to find yy mm dd for processing
    for yyyy in yy_list:
        for mm in mm_list:
            for dd in dd_list:
                
                # print the date to process
                print(f"processing date: {yyyy}-{mm}-{dd}")

                # process only if the date is valid
                if is_valid_date(yyyy, mm, dd):
                    
                    # check if the file is already on the bucket
                    file_exist = check_file_bucket(yyyy+mm+dd+'_RR_IT_15min_msg_res.nc.gz', domain='IT')

                    if file_exist:
                        print('file already exists on bucket - skipping date:', yyyy, mm, dd)
                        continue
                    else: 

                        print('new file - processing')
   
                        # sort filenames of the day and extract dates
                        data = sorted(file_names) # get all filenames from the bucket 
                        dates = [datetime.strptime(f.split("_")[-1].split(".")[0], "%Y%m%d%H%M") for f in data] # extract dates from filenames

                        # select only dates of the specified day
                        data = [f for f, date in zip(data, dates) if date.year == int(yyyy) and date.month == int(mm) and date.day == int(dd)]
                        dates = [date for date in dates if date.year == int(yyyy) and date.month == int(mm) and date.day == int(dd)]

                        print(f"Number of files found for the day {yyyy}-{mm}-{dd}: {len(data)}")

                        if len(data) == 0:
                            print("No files found for the specified day. Exiting.")
                            # add date to a log file: if the file already exists, append to it
                            with open('log_no_data_days_ARPAE.txt', 'a') as log_file:
                                log_file.write(f"{yyyy}-{mm}-{dd}\n")

                            return
                        else:

                            # unzip files, read data and create xr dataset for the entire day
                            ds_list = []

                            for i, file in enumerate(data):

                                # create string of the date with format with 4 digits for year, 2 for month, 2 for day, 2 for hour, 2 for minute
                                date_i = dates[i].strftime("%Y%m%d%H%M")
                                    
                                print(f"Processing file {i+1} of {len(data)}: {file} for date {date_i}")

                                # read file object from S3 bucket
                                file_obj = read_file_obj(s3, file, S3_BUCKET_NAME)

                                if file_obj is not None:

                                    # decompress gzip file object
                                    ds = xr.open_dataset(io.BytesIO(gzip.decompress(file_obj)))

                                    # append to list
                                    ds_list.append(ds)

                            # define date string for the day being processed
                            date_process = f"{yyyy}{mm}{dd}_"

                            # concatenate all datasets along time dimension of the day
                            ds_list = [ds for ds in ds_list if isinstance(ds, xr.Dataset)]

                            if len(ds_list) == 0:
                                print("No valid datasets to concatenate for this day.")
                                continue  # or handle as needed

                            try:
                                ds_day = xr.concat(ds_list, dim="time")
                            except Exception as e:
                                print(f"Error concatenating datasets: {e}")
                                continue  # or handle as needed


                            # add time coordinate
                            ds_day = ds_day.assign_coords(time=("time", dates))

                            # check time resolution of the dataset
                            time_diffs = np.diff(ds_day.time.values)  # differences between consecutive time points
                            unique_diffs = np.unique(time_diffs)
                            print("Unique time differences in the dataset:", unique_diffs)

                            # if only 5 min or only 15 min intervals are present, proceed
                            if len(unique_diffs) == 1 and unique_diffs[0] == np.timedelta64(5, 'm') or len(unique_diffs) == 1 and unique_diffs[0] == np.timedelta64(15, 'm'):   

                                # calculate time resolution in minutes
                                time_resolution_min = unique_diffs[0] / np.timedelta64(1, 'm')
                                print(f"Time resolution of the dataset: {time_resolution_min} minutes") 

                                print(f"Dataset has consistent {time_resolution_min}-minute intervals. Proceeding with processing.")

                                # convert Ze to RR by applying the Z-R relationship of function ze_to_rr
                                ds_day["RR"] = (("time", "lat", "lon"), calc_acc_rain_over_deltaT(ds_day["Z_60"].values, deltaT=time_resolution_min))

                            else:
                                print("Dataset does not have consistent 5-minute intervals. First resampling to 15-minute intervals.")

                                # resample ds_day Ze to the time array from 00 to 24 in 15 min time steps
                                ds_day = ds_day.resample(time='15T').sum(skipna=True)  
                                time_resolution_min = 15
                                # now convert Ze to RR by applying the Z-R relationship of function ze_to_rr
                                ds_day["RR"] = (("time", "lat", "lon"), calc_acc_rain_over_deltaT(ds_day["Z_60"].values, deltaT=time_resolution_min)) 

                                # store date to a log file for further investigation
                                #with open('log_inconsistent_time_intervals_ARPAE.txt', 'a') as log_file:
                                #    log_file.write(f"{yyyy}-{mm}-{dd}\n")   


                            # it time resolution is smaller than 15 min, resample to 15 min by summing up the rain amounts
                            if time_resolution_min < 15:
                                print(f"Resampling data from {time_resolution_min} min to 15 min by summing up rain amounts.")
                                
                                # resample to 15 min by summing up the rain amounts, skipping nans
                                ds_day = ds_day.resample(time='15T').sum(skipna=True)   
                                
                            # create 2D lat and 2d lon arrays for radar data
                            lat_2d_radar, lon_2d_radar = np.meshgrid(ds_day.lat.values, ds_day.lon.values, indexing='ij')
                            
                            # read one msg file to get temporal and spatial resolution
                            ds_msg = xr.open_dataset('/Users/claudia/Documents/Data/20220819-EXPATS-RG.nc')

                            # define domain_IT_CA for cropping msg data
                            lon_min, lon_max, lat_min, lat_max = domain_IT_CA

                            # crop msg data to the IT domain (smaller than EXPATS)
                            ds_msg_IT = ds_msg.where((ds_msg.lat >= lat_min) & 
                                                    (ds_msg.lat <= lat_max) & 
                                                    (ds_msg.lon >= lon_min) & 
                                                    (ds_msg.lon <= lon_max), 
                                                    drop=True)
                            
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
                            RR_msg = np.zeros((len(ds_day.time.values), len(lat_arr), len(lon_arr)))
                                        
                            # loop on time stamps to resample 15 min rain amounts to MSG spatial resolution
                            with Bar('Processing...') as bar: 
                                
                                for i_t, time_val in enumerate(ds_day.time.values):
                                    
                                    print('regrid RR for time ', time_val)

                                    # setting time obs for loop start
                                    time_loop_start = time.time()

                                    # regrid rain rate date from old radar grid to the new msg grid
                                    RR_msg[i_t, :, :] = regrid_data(lat_2d_radar,
                                                        lon_2d_radar,
                                                        ds_day.RR.values[i_t,:,:],
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
                                    "time": ds_day.time.values,
                                    "lat": lat_arr, 
                                    "lon": lon_arr,
                                },
                                attrs={
                                    "description": "Rain rate from ARPAE data resampled to MSG grid as cumulated over 15 minutes",
                                    "history": "Created on " + str(pd.Timestamp.now()) + " by Claudia Acquistapace, ",
                                    "source": "ARPAE Radar",
                                    "reference history":"data provided by Dr. Virginia Politi from ARPAE",
                                    "created_by": "Claudia Acquistapace",
                                    "created_on": str(pd.Timestamp.now()),
                                    "domain": "IT_CA",
                                    "original_grid": "ARPAE radar polar grid",
                                    "target_grid": "MSG-like regular grid"
                                }
                            )

                            # save to ncdf and compress RR variable to maximum compression
                            path_out_nc = "/Users/claudia/Documents/Data/ARPAE_15min_rain_rate/"
                            ds_RR_msg.to_netcdf(path_out_nc+date_process+'RR_IT_15min_msg_res.nc', encoding={'RR': {'zlib': True, 'complevel': 9}})
        
                            # gzip the nc file
                            with open(path_out_nc+date_process+'RR_IT_15min_msg_res.nc', 'rb') as f_in:
                                with gzip.open(path_out_nc+date_process+'RR_IT_15min_msg_res.nc.gz', 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                    
                            # remove the local nc file
                            os.remove(path_out_nc+date_process+'RR_IT_15min_msg_res.nc')
                            
                            # upload the processed file to the new S3 bucket
                            file_to_upload = os.path.join(path_out_nc, date_process+'RR_IT_15min_msg_res.nc.gz')
                            upload_success = upload_file(s3, file_to_upload, s3_bucket_new_name, date_process+'RR_IT_15min_msg_res.nc.gz')
                            

                            if upload_success:
                                print(f"File {date_process+'RR_IT_15min_msg_res.nc.gz'} gzipped and uploaded successfully to bucket {s3_bucket_new_name}.")
                            else:
                                print(f"Failed to upload file {date_process+'RR_IT_15min_msg_res.nc'} to bucket {s3_bucket_new_name}.") 
                            
                            


    
if __name__ == "__main__":
    main()




"""

def plot_radar_timestep_map(ds_day, dates, lon_min, lon_max, lat_min, lat_max):
        
        # plot data for a given time
        time_index = 3

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.spines["top"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.set_extent(domain_expats, crs=ccrs.PlateCarree())


        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'fontsize': 14}
        gl.ylabel_style = {'fontsize': 14}
        

        vmin= -10.
        vmax= 60.
        var_levels = np.arange(vmin, vmax, 10)

        # plot rain rate as filled contours
        lats = ds_day.lat.values
        lons = ds_day.lon.values
        Ze = ds_day.Z_60.values
        
        # reading orography data from raster file
        ds_or = read_orography()
        oro_levels = np.linspace(0, 1500, 20)
        oro = ax.contourf(ds_or.lons.values, 
                            ds_or.lats.values, 
                            ds_or.orography.values, 
                            transform=ccrs.PlateCarree(), 
                            levels=oro_levels, 
                            alpha=.5,
                            cmap='Greys')

        mesh_ZE = ax.contourf(lons, 
                            lats, 
                            Ze[time_index,:,:], 
                            cmap='BuPu_r', 
                            transform=ccrs.PlateCarree(), 
                            vmin=vmin,  
                            vmax=vmax, 
                            levels=var_levels, 
                            extend='max',
                            alpha=0.6)  
        
        
        # define discrete colorbar for reflectivity with interval of 10 dbZ
        cbar = plt.colorbar(mesh_ZE, label='Reflectivity [dBZ]', shrink=0.6, ticks=var_levels)
        cbar.ax.set_yticklabels([f"{int(l)}" for l in var_levels])
        cbar.ax.tick_params(labelsize=12)
        #plot_cities_expats(ax, 'black', 50)
        ax.add_feature(cfeature.BORDERS, linewidth=1., color='black')
        ax.coastlines(linewidth=1., color='black')
        
        # draw a rectangle for the lat max lat min lon max lon min of the radar data
        ax.plot([lon_min, lon_max, lon_max, lon_min, lon_min],
                [lat_min, lat_min, lat_max, lat_max, lat_min],
                color='black', linewidth=2, linestyle=":", transform=ccrs.PlateCarree())
        plt.savefig(
            os.path.join('plots/', 'ZE_'+dates[time_index].strftime("%Y%m%d%H%M")+'.png'),
            dpi=300,
            bbox_inches="tight",
            transparent=True,
            )

        plt.close()

        # drop varibles not needed
        ds_day = ds_day.drop_vars(['Z_60'])


        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.spines["top"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.set_extent(domain_expats, crs=ccrs.PlateCarree())


        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'fontsize': 14}
        gl.ylabel_style = {'fontsize': 14}
        

        vmin= 0.
        vmax= 40.
        var_levels = np.arange(vmin, vmax, 5)

        # plot rain rate as filled contours
        lats = ds_day.lat.values
        lons = ds_day.lon.values
        RR = ds_day.RR.values
        
        # reading orography data from raster file
        ds_or = read_orography()
        oro_levels = np.linspace(0, 1500, 20)
        oro = ax.contourf(ds_or.lons.values, 
                            ds_or.lats.values, 
                            ds_or.orography.values, 
                            transform=ccrs.PlateCarree(), 
                            levels=oro_levels, 
                            alpha=.5,
                            cmap='Greys')

        mesh_ZE = ax.contourf(lons, 
                            lats, 
                            Ze[time_index,:,:], 
                            cmap='BuPu_r', 
                            transform=ccrs.PlateCarree(), 
                            vmin=vmin,  
                            vmax=vmax, 
                            levels=var_levels, 
                            extend='max',
                            alpha=0.6)  
        
        
        # define discrete colorbar for reflectivity with interval of 10 dbZ
        cbar = plt.colorbar(mesh_ZE, label='Rain Rate [mm/h]', shrink=0.6, ticks=var_levels)
        cbar.ax.set_yticklabels([f"{int(l)}" for l in var_levels])
        cbar.ax.tick_params(labelsize=12)
        #plot_cities_expats(ax, 'black', 50)
        ax.add_feature(cfeature.BORDERS, linewidth=1., color='black')
        ax.coastlines(linewidth=1., color='black')
        
        # draw a rectangle for the lat max lat min lon max lon min of the radar data
        ax.plot([lon_min, lon_max, lon_max, lon_min, lon_min],
                [lat_min, lat_min, lat_max, lat_max, lat_min],
                color='black', linewidth=2, linestyle=":", transform=ccrs.PlateCarree())
        plt.savefig(
            os.path.join('plots/', 'RR_'+dates[time_index].strftime("%Y%m%d%H%M")+'.png'),
            dpi=300,
            bbox_inches="tight",
            transparent=True,
            )

        plt.close()


"""