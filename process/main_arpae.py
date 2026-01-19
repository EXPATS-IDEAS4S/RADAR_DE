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
from readers.config import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from figures.domain_info import domain_expats, domain_DE_CA
from figures.one_day_video import plot_radar_map
from process.utils import generate_regular_grid, regrid_data, is_valid_date
from readers.radar_DWD import read_orography
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg
from readers.file_dirs import orography_file
from readers.data_buckets_funcs import read_file_obj, init_s3
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


def main():

    # read all filenames from bucket
    s3, file_names = read_all_filenames_from_bucket()

    # sort filenames in list radar files and extract those of the selected day
    yyyy = "2025"
    mm = "04"
    dd = "01"

    # sort filenames of the day and extract dates
    data = sorted(file_names)
    dates = [datetime.strptime(f.split("_")[-1].split(".")[0], "%Y%m%d%H%M") for f in data]

    # select only dates of the specified day
    data = [f for f, date in zip(data, dates) if date.year == int(yyyy) and date.month == int(mm) and date.day == int(dd)]
    dates = [date for date in dates if date.year == int(yyyy) and date.month == int(mm) and date.day == int(dd)]

    print(f"Number of files found for the day {yyyy}-{mm}-{dd}: {len(data)}")

    if len(data) == 0:
        print("No files found for the specified day. Exiting.")
        return
    else:

        # unzip files, read data and create xr dataset for the entire day
        ds_list = []

        for i, file in enumerate(data):
            #print(f"Reading file {i+1} of {len(data)}: {file}") 
            # read file object from S3 bucket
            file_obj = read_file_obj(s3, file, S3_BUCKET_NAME)

            if file_obj is not None:

                # decompress gzip file object
                ds = xr.open_dataset(io.BytesIO(gzip.decompress(file_obj)))

                # append to list
                ds_list.append(ds)

        #
        ds_day = xr.concat(ds_list, dim="time")

        # add time coordinate
        ds_day = ds_day.assign_coords(time=("time", dates))

        print(np.min(ds_day.Z_60.values), np.max(ds_day.Z_60.values))
        # substitute values of Z_60 < -10 dBZ with NaN
        ds_day["Z_60"] = ds_day["Z_60"].where(ds_day["Z_60"] >= -10., np.nan)

        # convert Ze to RR by applying the Z-R relationship of function ze_to_rr
        ds_day["RR"] = (("time", "lat", "lon"), ze_to_rr(ds_day["Z_60"].values))

        print(np.min(ds_day.Z_60.values), np.max(ds_day.Z_60.values))

        # resample dataset to 15 min intervals
        ds_radar_msg = ds_day.resample(time='15T').sum(skipna=True)


        # read one msg file to get temporal and spatial resolution
        ds_msg = xr.open_dataset('/Users/claudia/Documents/Data/'+yyyy+'/'+mm+'/'+yyyy+mm+dd+'-EXPATS-RG.nc')


        # extract lat max lat mon lon max and lon min from the radar data
        lat_min = ds_day.lat.min().values
        lat_max = ds_day.lat.max().values
        lon_min = ds_day.lon.min().values
        lon_max = ds_day.lon.max().values

        # crop msg data to domain_DE_CA
        lon_min, lon_max, lat_min, lat_max = domain_IT_CA


        # crop msg data to the DE domain (smaller than EXPATS)
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
        RR_msg = np.zeros((len(ds_radar_msg.time.values), len(lat_arr), len(lon_arr)))
        
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

    
if __name__ == "__main__":
    main()

