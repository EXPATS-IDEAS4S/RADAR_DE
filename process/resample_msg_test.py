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
from process.one_day_video import plot_radar_map
import matplotlib.pyplot as plt
import os
from time import sleep
from progress.bar import Bar
import time
import pandas as pd
import shutil

def main():

    # download_from_s3()
    # loop on all days of
    # set the day to process
    
    # list all years from 2024 to 2014
    yy_list = [str(year) for year in range(2023, 2021, -1)]
    # list all months in a year
    mm_list = [f'{month:02d}' for month in range(4, 9)]
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
                    file_exist = check_file_bucket(yy+mm+dd+'_RR_15min_msg_res.nc')
                    
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
                            ncdf_radar = process_date(yy, mm, dd, plotting, path_out)

                            # call function to transfer data to bucket, returns true if transfer is complete
                            check = upload_to_bucket(path_out, ncdf_radar)
                            
                        else:
                            # if file radolan does not exist, print message and continue to next day
                            print('radolan file missing - go to next day')
                            continue
                        
                        # check if file was uploaded
                        if check:
                            print('file uploaded to bucket')
                            
                            # remove ncdf file produced
                            #os.remove(path_out+ncdf_radar)
                            
                        else:
                            print('file not uploaded to bucket - abort')
                            break




def check_file_bucket(file2find):
    """code to check if a file is on the bucket
    input:
    file2find: filename with path to find n bucket
    return:
    filefound (boolean) true if file is found, false if file is not on the bucket
    """
    
    from readers.s3_buckets_credentials import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file
    
    import boto3
    
    s3 = boto3.client(
    's3',
    endpoint_url=S3_ENDPOINT_URL,
    aws_access_key_id=S3_ACCESS_KEY,
    aws_secret_access_key=S3_SECRET_ACCESS_KEY)
    
    # List the objects in our bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    
    # set flag to false if file is not found
    filefound = False
    if "Contents" not in response:
        print(f"No files on bucket")
        
    
    # loop on files found on bucket
    for obj in response["Contents"]:
        key = obj["Key"]
        if key == file2find:
            filefound = True
            return filefound
        else:
            continue

    return filefound


def upload_to_bucket(path_out, filename_ncdf):
    """
    function to upload on a specified S3 bucket the filename ncdf

    Args:
        path_out (string): path to the directory where filename_ncdf is located 
        filename_ncdf (string): name of the file to upload on bucket

    Returns:
        check: boolean, true if file is uploaded, false if not
    """
    
    import time
    from glob import glob

    from readers.s3_buckets_credentials import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file


    start_time = time.time()

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # calling upload function and return boolean for upload status
    # upload_file(s3_client, file_name, bucket, object_name=None):
    file_to_upload = os.path.join(path_out, filename_ncdf)
    check = upload_file(s3, file_to_upload, S3_BUCKET_NAME, filename_ncdf)
    
    # if true, upload done
    if check:
        print("Time taken to upload files: ", time.time() - start_time, flush=True) 
        print('file name of file uploaded on bucket', filename_ncdf)

    return check
    
    
    
    
    
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
    
def process_date(yy,mm,dd, plotting,path_out):
    '''
    function to read rain rate from dwd data, resample it on the msg temporal (15 min) and spatial (0.04 degrees) resolution, and then store a new netcdf on the small radar domain. It also, if requested by the plotting keyword, produces plots of the rr instantaneous, cumulated over 15 mins, and msg 
    inputs:
    - yy: string of the year
    - mm: string of the month
    - dd: string of the day
    - plotting: boolean, if True produces plots
    - path_out: string, path to output directory where to store ncdf and plots
    returns
    - filename of the file ncdf produced.
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
    if not os.path.exists(os.path.join(path_out, date+'_RR_15min_msg_res.nc')):
        
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

    else:
        
        print('file RR at msh res exists - reading it ')
        # read stored ncdf
        ds_RR_msg = xr.open_dataset(os.path.join(path_out, date+'_RR_15min_msg_res.nc'))    

    # gzip the file
    with open(os.path.join(path_out, date+'_RR_15min_msg_res.nc'), 'rb') as f_in:
        with gzip.open(os.path.join(path_out, date+'_RR_15min_msg_res.nc.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

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

    return(date+'_RR_15min_msg_res.nc')


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




def plot_radar_mapV2(data, string_file, path_out):
    """
    function to plot map of radar data from DWD, in rain rate

    Args:
        data (xarray dataset): data from DWD 
    """
    print(np.nanmax(data.RR.values))
    print(np.nanmin(data.RR.values))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.spines["top"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.set_extent(domain_DE_CA)
    
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 14}
    gl.ylabel_style = {'fontsize': 14}
    

    vmin= 0.
    vmax= 5.
    var_levels = np.linspace(vmin, vmax, 15)
    # plot rain rate as filled contours
    lats = data.lat.values
    lons = data.lon.values
    RR = data.RR.values
    
    RR[RR==0.] = np.nan
    # reading orography data from raster file
    ds_or = read_orography()
    oro_levels = np.linspace(0, 1500, 20)
    oro = ax.contourf(ds_or.lons.values, 
                        ds_or.lats.values, 
                        ds_or.orography.values, 
                        transform=ccrs.PlateCarree(), 
                        levels=oro_levels, 
                        alpha=1.,
                        cmap='Greys')

    mesh_rr = ax.contourf(lons, 
                        lats, 
                        RR, 
                        cmap='BuPu_r', 
                        transform=ccrs.PlateCarree(), 
                        vmin=vmin,  
                        vmax=vmax, 
                        levels=var_levels, 
                        extend='max',
                        alpha=0.6)  
        
    cbar = plt.colorbar(mesh_rr, label='Rain rate [mm]', shrink=0.6)
    plot_cities_expats(ax, 'black', 50)
    ax.add_feature(cfeature.BORDERS, linewidth=1., color='black')
       
    plt.savefig(
        os.path.join(path_out, string_file+'.png'),
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        )

    plt.close()
    return()







def plot_msg(lons, lats, variable, label, cmap, cbar_title, domain, path_out, key, back_transparent, vmin, vmax):
    """
    plot map of the variable over the defined domain

    Args:
        lons (array): longitude values
        lats (array): latitudes 
        variable (matrix): variable to plot
        label (string): string for plot filename
        cmap (colormap): color map 
        cbar_title (string): title for color bar
        domain (array): minlon, maxlon, minlat, maxlat     

        path_out (string): path where to save the plot
        key (string): string possible "expats" or "dfg"
        back_transparent (boolean): True or False - true means transparent background
        vmin (float): value min in colorbar
        vmax (float): max value in colorbar
        
    Dependencies:
    plot_cities_expats, 
    plot_local_dfg
    
    """
    from figures.mpl_style import plot_cities_expats

    # Plot the map of MSG channels vs lat/lon
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.spines["top"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.set_extent(domain)
    
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 14}
    gl.ylabel_style = {'fontsize': 14}
    
    

    var_levels = np.linspace(vmin, vmax, 20)
    # plot 10th percentile as filled contours
    mesh = ax.contourf(lons, 
                        lats, 
                        variable, 
                        cmap=cmap, 
                        transform=ccrs.PlateCarree(), 
                        levels=var_levels, 
                        vmin=vmin, # 230
                        vmax=vmax) # 270
    cbar = plt.colorbar(mesh, label=cbar_title, shrink=0.6)

    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    
     #vmin, vmax=220.,
    plt.tick_params(axis='both',which='major',labelsize=14)

    ax.set_xlabel('Latitude [$^{\circ}$]')
    ax.set_ylabel('Longitude [$^{\circ}$]')
    ax.tick_params(which='minor', length=5, width=2)
    ax.tick_params(which='major', length=7, width=3)
    #cbar = plt.colorbar(pc,ax=ax,shrink=0.75)
    #cbar.set_label(label, fontsize=14)
    #cbar.ax.tick_params(labelsize=14)

    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=1.5, color='black')
    #ax.add_feature(cfeature.STATES, linewidth=0.2)
    ax.add_feature(cfeature.BORDERS, linewidth=1.5, color='black')

    # Title and save date_str+'_IR108'
    ax.set_title(label, fontsize=12)
    
    if key =='expats':
        plot_cities_expats(ax, 'grey', 50)

        
        
    plt.savefig(
        os.path.join(path_out, label+"_"+key+".png"),
        dpi=300,
        bbox_inches="tight",
        transparent=back_transparent,
        )

    plt.close()
#f"{i:03d}"

    
    

if __name__ == "__main__":
    main()