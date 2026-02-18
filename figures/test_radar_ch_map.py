"""
Code to plot the 21/06/2024
radar map of Switzerland, to check the coverage of the radar network in the area of interest.
"""

from datetime import date, datetime
import glob
import shutil
import zipfile
import h5py
from matplotlib.pyplot import contourf
import xarray as xr
import pdb
import os
import numpy as np
from pyproj import CRS, Transformer
from figures.plotting import plot_radar_mapV2
from process.utils import extract_date_from_filename_ch
from readers.radar_DWD import read_orography
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from readers.radars_ch import read_h5_file
import pandas as pd
import ffmpeg
from figures.make_gif import main as make_gif_main
def main():

    # find out which day of the year is the 21/06 
    date_str = "2024-06-21"
    date_obj = datetime.strptime(date_str, "%Y-%m-%d")
    day_of_year = date_obj.timetuple().tm_yday
    print(f"The day of the year for {date_str} is: {day_of_year}")  

    # construct name as yyddd where ddd is the day of the year and yy is the last two digits of the year
    year = date_obj.year
    yy = str(year)[-2:]
    ddd = str(day_of_year).zfill(3)
    filename = f"RZCflt{yy}{ddd}.zip"
    print(f"filename for the day {date_str} is: {filename}")

    # look for the file in the folder and unzip it, 
    path_h5 = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/"
    file_path = os.path.join(path_h5, filename)

    if os.path.exists(file_path):

        dest_folder = "/Users/claudia/Documents/Data//test_radar/"
        # copy file_path in dest_folder if not already copied
        if not os.path.exists(os.path.join(dest_folder, filename)):
            shutil.copy(file_path, dest_folder)
            print(f"file {filename} copied successfully in {dest_folder}")
        else:
            print(f"file {filename} already copied in {dest_folder}")

        # unzip the file
        extract_folder_name = filename.split(".zip")[0]
        folder_path_unzipped = os.path.join(path_h5, extract_folder_name)
        if not os.path.exists(folder_path_unzipped):
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(folder_path_unzipped)
            print(f"file {filename} unzipped successfully")
        else:
            print(f"file {filename} already unzipped")

    # list the h5 files in the unzipped folder 
    h5_files = sorted(glob.glob(os.path.join(folder_path_unzipped, "*.h5")))
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

    # plot 
    # plot the radar map of Switzerland using the function plot_radar_mapV2 
    # from figures.plotting for each time stamp
    for time in ds_day.time:
        ds_5min = ds_day.sel(time=time)
        # construct filename for the plot as radar_ch_map_yyyymmdd_hhmmss
        time_str = pd.to_datetime(str(time.values)).strftime("%Y%m%d_%H%M%S")
        string_file = f"{time_str}radar_ch_map.png"
        path_out = "/Users/claudia/Documents/Data/test_radar/png/" 
        print(f"plotting radar map for time {time.values}")
        # if plot does not exist already in the folder, plot it, otherwise skip
        if not os.path.exists(os.path.join(path_out, string_file)):
            plot_radar_mapV2(ds_5min, string_file, path_out)
        else:
            print(f"plot {string_file} already exists in {path_out}, skipping plotting")



if __name__ == "__main__":
    main()
    