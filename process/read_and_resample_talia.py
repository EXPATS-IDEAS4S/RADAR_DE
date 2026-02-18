"""
Code to read and resample radar data from Talia
 source /Users/claudia/Github/RADAR_DE/venv/bin/activate
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

author: Claudia Acquistapace
date: 16 Febbruary 2026
email: claudia.acquistapace@unipd.it

"""
import zipfile
import os
from glob import glob
import xarray as xr
import numpy as np
import pdb
import subprocess
from compression import gzip
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

def main():

    # reading list of files from folder
    radar_ch_path = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/"
    
    # list all files RZCflt*.zip in the folder (these are the files with 32 bits floating point data)
    files_list = glob(os.path.join(radar_ch_path, "*RZCflt*.zip"))

    print(f"files in the folder: {files_list}")
    print(f"number of files in the folder: {len(files_list)}")
    print(files_list[0])
    print("******************************************************")

    for file in files_list:
        
        print(f"processing file {file}")
        print("******************************************************")
        
        # content of the loop on the file_list
        zip_path = files_list[0] # zipped file
        extract_folder_name = zip_path.split("/")[-1].split(".zip")[0]
        print(f"filename: {zip_path}, folder name unzipped: {extract_folder_name}")
        print("******************************************************")

        folder_path_unzipped = os.path.join(radar_ch_path, extract_folder_name)
        # unzip if not already unzipped
        # check if folder with name filename without .zip already exists, if not, unzip the file

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
            print(f"first 5 .h5 files in the unzipped folder: {h5_files[:5]}")
        
        # extract date from filename
        date = extract_date_from_filename_ch(extract_folder_name)
        print(f"date of filename: {date}")

        # read all h5 files in the unzipped folder 
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

        # store ds_day_interp in ncdf 

        # define destination path 
        path_out_nc = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/15_minutes/"
        ncdf_15min_filename = os.path.join(path_out_nc, date+'_RR_CH_15min_msg_res.nc')
        ds_day_interp.to_netcdf(ncdf_15min_filename, encoding={'RR': {'zlib': True, 'complevel': 9}})

        # plot rain rate
        plot_radar_mapV2(ds_day, 
                        string_file=f"radar_ch_{date}", 
                        path_out="/Users/claudia/Documents/Data/test_radar/")
        
        plot_radar_mapV2(ds_day_interp, 
                        string_file=f"radar_ch_{date}_15min_msg_res", 
                        path_out="/Users/claudia/Documents/Data/test_radar/")
        
        # if plot of the coverage does not exist plot it
        if not os.path.exists(os.path.join("/Users/claudia/Documents/Data/test_radar/", f"radar_ch_coverage_{date}.png")):
          plot_radar_coverage(ds_day,
                            string_file=f"radar_ch_coverage_{date}", 
                              path_out="/Users/claudia/Documents/Data/test_radar/")

        # remove folder and folder zippef
        # If netcdf are in the destination folders then remove h5 files
        if os.path.exists(ncdf_15min_filename) and os.path.exists(ncdf_5min_filename):
            print(f"file {ncdf_15min_filename} exists, removing unzipped folder and h5 files")
            # remove unzipped folder
            shutil.rmtree(folder_path_unzipped)

            # remove zipped file
            shutil.rmtree(zip_path)

            # compress file with gzip with open(ncdf_5min_filename, 'rb') as f_in:
            with gzip.open(ncdf_5min_filename+'.gz', 'wb') as f_out:
                shutil.copyfileobj(ncdf_5min_filename, f_out) 
            with gzip.open(ncdf_15min_filename+'.gz', 'wb') as f_out:
                shutil.copyfileobj(ncdf_15min_filename, f_out)

        else:
            print(f"file {ncdf_15min_filename} or {ncdf_5min_filename} does not exist, not removing unzipped folder and h5 files")


  

if __name__ == "__main__":
    main()
