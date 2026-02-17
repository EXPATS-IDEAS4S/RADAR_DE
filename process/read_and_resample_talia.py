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
from process.utils import extract_date_from_filename_ch
import h5py
from readers.radars_ch import read_h5_file
from figures.plotting import plot_radar_mapV2
from figures.mpl_style import plot_cities_expats
from readers.radar_DWD import read_orography
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from figures.domain_info import domain_expats, domain_DE_CA
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm

def main():

    # reading list of files from folder
    radar_ch_path = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/"
    
    # list all files RZCflt*.zip in the folder (these are the files with 32 bits floating point data)
    files_list = glob(os.path.join(radar_ch_path, "*RZCflt*.zip"))

    print(f"files in the folder: {files_list}")
    print(f"number of files in the folder: {len(files_list)}")
    print(files_list[0])
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
    h5_files = glob(os.path.join(folder_path_unzipped, "*.h5"))
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

    # plot rain rate
    plot_radar_mapV2(ds_day, 
                    string_file=f"radar_ch_{date}", 
                    path_out="/Users/claudia/Documents/Data/test_radar/")
    
    plot_radar_coverage(ds_day,
                        string_file=f"radar_ch_coverage_{date}", 
                        path_out="/Users/claudia/Documents/Data/test_radar/")

    pdb.set_trace()
    # TO DO 



    """
      - resample to msg time/space grid
      - save resampled data on ncdf file with name RZCflt_yyyymmdd_hhmm.nc
      - move ncdf file on bucket

    """
    pdb.set_trace()



if __name__ == "__main__":
    main()
