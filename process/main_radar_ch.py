"""
Code to read, geolocalize and then process CH radar data to resample them on 
msg spatial resolution and temporal resolution as done in the other radar data

author: Claudia Acquistapace
date: 16 Febbruary 2026
email: claudia.acquistapace@unipd.it

To run this code:
- activate the virtual environment with all the necessary libraries installed
source /Users/claudia/Github/venv/bin/activate 

- run code with python
python3 -m process.main_radar_ch
"""
import zipfile
import os
from glob import glob
from pandas import plotting
import xarray as xr
import numpy as np
import pdb
import subprocess
import gzip
from process.utils import extract_date_from_filename_ch
import h5py
from readers.radars_ch import read_h5_file
from figures.plotting import plot_radar_mapV2, plot_radar_coverage
from figures.mpl_style import plot_cities_expats
from readers.radar_DWD import read_orography
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from figures.domain_info import domain_expats, domain_DE_CA
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm
from progress.bar import Bar
import shutil
from process.utils import select_months_of_interest

def main():

    # set plotting to True to plot the radar maps and coverage, set to False to not plot
    plotting = False

    # reading list of files from folder
    radar_ch_path = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/"
    
    # list all files RZCflt*.zip in the folder (these are the files with 32 bits floating point data)
    files_list = sorted(glob(os.path.join(radar_ch_path, "*RZCflt*.zip")))
    print(f"number of files in the folder: {len(files_list)}")
    files_list = select_months_of_interest(files_list, radar_ch_path, ['04', '05', '06', '07', '08', '09'])  
    print("******************************************************")

    for file in files_list:

        print(f"processing file {file}")
        print("******************************************************")
        
        # content of the loop on the file_list
        zip_path = file # zipped file
        extract_folder_name = zip_path.split("/")[-1].split(".zip")[0]
        # extract date from filename
        date = extract_date_from_filename_ch(extract_folder_name)
        print(f"date of filename: {date}, folder name: {extract_folder_name}")
        print(f"filename: {zip_path}, folder name unzipped: {extract_folder_name}")
        print("******************************************************")

        # if nc.gz file exists, skip the file, otherwise process it
        ncdf_5min_filename = os.path.join("/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/5_minutes/", date+"_radar_swisse_5min_res.nc.gz")
        ncdf_15min_filename = os.path.join("/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/15_minutes/", date+'_RR_CH_15min_msg_res.nc.gz')
        if os.path.exists(ncdf_5min_filename) and os.path.exists(ncdf_15min_filename):
            print(f"file {ncdf_5min_filename} and {ncdf_15min_filename} already exist, skipping file {zip_path}")
            # go to next file in the list and skip the rest of the loop
    
            continue
        
        else:
            print(f"NCDF files do not exist, processing file {zip_path}")
            # set folder path for unzipped files
            folder_path_unzipped = os.path.join(radar_ch_path, extract_folder_name)

            # check if folder with name filename without .zip already exists, if not, unzip the file
            if not os.path.exists(folder_path_unzipped):
                os.makedirs(folder_path_unzipped, exist_ok=True)
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(folder_path_unzipped)

            # check if the unzipped folder contains files with extension .h5, 
            h5_files = sorted(glob(os.path.join(folder_path_unzipped, "*.h5")))
            if len(h5_files) == 0:
                print(f"no .h5 files in the unzipped folder {folder_path_unzipped}")
                return
            else:
                print(f"number of .h5 files in the unzipped folder: {len(h5_files)}")
            

            # read all h5 files in the unzipped folder 
            """ this reading function reads all the data of the day:
                - gets the rr, 
                - adjust RR with cal factor 0.87 (Francesco Marra, personal comm, summer period)
                - constructs the lat lon grid at the original radar resolution, 
                - concatenates all data of the day along the time dimension
                - set to nan all RR values = 0 (these are no data values)
                - creates a domain mask to define the radar domain 
                - stores everything in a xarray dataset and returns it as output."""
            ds_day = read_h5_file(h5_files, date)

            # store ds_day in ncdf file with name RZCflt_yyyymmdd_hhmm.nc applying max data compression
            ncdf_filename = f"{date}_radar_swisse_5min_res.nc"
            ncdf_5min_filename = os.path.join("/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/5_minutes/",
                                    ncdf_filename)
            ds_day.to_netcdf(ncdf_5min_filename, encoding={'RR': {'zlib': True, 'complevel': 9}})

            # convert rr in mm/h to accumulated rain in mm over the time resolution of the dataset
            time_res = int(ds_day.attrs['time_resolution'][0:-8])/ 60 # time res in minutes
            ds_day['RR'] = ds_day['RR'] * time_res / 60 # convert from mm/h to mm over the time resolution of the dataset
            
            # resample to msg time/space grid (15 min time resolution)
            ds_day = ds_day.resample(time='15min').mean(skipna=True)  

            # read one msg file to get temporal and spatial resolution
            ds_msg = xr.open_dataset('/Users/claudia/Documents/Data/20220819-EXPATS-RG.nc')
            
            # interpolate RR to MSG lat lon grid
            lat_msg = ds_msg.lat.values
            lon_msg = ds_msg.lon.values
            ds_day_interp = ds_day.interp(lat=lat_msg, lon=lon_msg, method='nearest')
            # change ds_day_interp attributes to add comment on the interpolation method used
            ds_day_interp.attrs['comment'] = "dataset created by interpolating the original radar dataset \n" \
            "to the MSG lat lon grid with nearest neighbor interpolation, using the interp function of xarray. \n" \
            "The original radar dataset is created with the function read_h5_file in readers/radars_ch.py \n" \
            # store ds_day_interp in ncdf 

            # define destination path 
            path_out_nc = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/15_minutes/"
            ncdf_15min_filename = os.path.join(path_out_nc, date+'_RR_CH_15min_msg_res.nc')
            ds_day_interp.to_netcdf(ncdf_15min_filename, encoding={'RR': {'zlib': True, 'complevel': 9}})

            # plot rain rate
            if plotting:
                plot_radar_mapV2(ds_day, 
                                string_file=f"{date}_radar_ch", 
                                path_out="/Users/claudia/Documents/Data/test_radar/")
                
                plot_radar_mapV2(ds_day_interp, 
                                string_file=f"{date}_radar_ch_15min_msg_res", 
                                path_out="/Users/claudia/Documents/Data/test_radar/")
            
                # if plot of the coverage does not exist plot it
                if not os.path.exists(os.path.join("/Users/claudia/Documents/Data/test_radar/", f"radar_ch_coverage_{date}.png")):
                    plot_radar_coverage(ds_day.sel(time=ds_day.time[0]), # select the first time step of the day to plot the coverage, since the radar domain is not expected to change during the day
                                        string_file=f"radar_ch_coverage_{date}", 
                                        path_out="/Users/claudia/Documents/Data/test_radar/")

            # remove folder and folder zipped
            # If netcdf are in the destination folders then remove h5 files
            if os.path.exists(ncdf_15min_filename) and os.path.exists(ncdf_5min_filename):
                print(f"file {ncdf_15min_filename} exists, removing unzipped folder and h5 files")
                
                try:
                    if os.path.exists(folder_path_unzipped):
                        shutil.rmtree(folder_path_unzipped)
                        print(f"unziped folder {folder_path_unzipped} removed")
                    else:
                        print(f"unziped folder {folder_path_unzipped} does not exist, not removing it")
                except Exception as e:
                    print(f"error while removing unzipped folder {folder_path_unzipped}: {e}")

                try:
                    if os.path.exists(zip_path):
                        os.remove(zip_path)
                        print(f"zipped file {zip_path} removed")
                    else:
                        print(f"zipped file {zip_path} does not exist, not removing it")
                except Exception as e:
                    print(f"error while removing zipped file {zip_path}: {e}")

                # compress file with gzip with open(ncdf_5min_filename, 'rb') as f_in:
                try:
                    with open(ncdf_5min_filename, 'rb') as f_in:
                        with gzip.open(ncdf_5min_filename + '.gz', 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            # remove original ncdf file after compressing it
                            os.remove(ncdf_5min_filename)
                    with open(ncdf_15min_filename, 'rb') as f_in:
                        with gzip.open(ncdf_15min_filename + '.gz', 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            # remove original ncdf file after compressing it
                            os.remove(ncdf_15min_filename)


                except Exception as e:
                    print(f"error while compressing files {ncdf_5min_filename} and {ncdf_15min_filename} with gzip: {e}")

            else:
                print(f"file {ncdf_15min_filename} or {ncdf_5min_filename} does not exist, not removing unzipped folder and h5 files")


if __name__ == "__main__":
    main()