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

from process.main_arpae import read_all_filenames_from_bucket
from readers.file_dirs import path_radolan_DE
from readers.radar_DWD import read_radar_DWD, read_orography
from readers.config import *
from readers.data_buckets_funcs import delete_file_from_bucket, upload_to_bucket, check_file_bucket
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




def main():

    # output path for processed files
    s3_bucket_new_name = "expats-radar-italy" # name of the destination data bucket where to store daily processed files of RR from italian domain

    # read files to be removed from log_to_delete_and_reprocess.txt
    #with open('log_to_delete_and_reprocess_ARPAE.txt', 'r') as f:
    #    lines = f.readlines()   
    #dates_to_process = [line.strip() for line in lines] 
   # print(dates_to_process)
    
    # initialize s3 connection

    s3 = init_s3()  
    #for date_str in dates_to_process:

        # dates are in the form 2022-08-04
    #    yyyy = date_str[0:4]
    #    mm = date_str[5:7]
        #dd = date_str[8:10]

    #    delete_file_from_bucket(s3, s3_bucket_new_name, yyyy+mm+dd+'_RR_IT_15min_msg_res.nc.gz')    

    #    print(f"file {yyyy+mm+dd+'_RR_IT_15min_msg_res.nc.gz'} removed.")

    # remove all files with yyyy < 2025
    for year in ['2024', '2023', '2022', '2021', '2020', '2019', '2018', '2017', '2016']:
        for month in [f'{m:02d}' for m in range(1, 13)]:
            for day in [f'{d:02d}' for d in range(1, 32)]:
                if is_valid_date(year, month, day):
                    filename = f"{year}{month}{day}_RR_IT_15min_msg_res.nc.gz"
                    delete_file_from_bucket(s3, s3_bucket_new_name, filename)
                    print(f"file {filename} removed.")  
    
if __name__ == "__main__":
    main()
