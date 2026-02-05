
"""
this script identifies all dates in which radar data from ARPAE:
- either have a time resolution that is not constant (e.g., mix of 5 min and 15 min)
- or have a time resolution that is not at 00, 15, 30, 45 minutes
The dates are logged in a text file for further reprocessing with the code main_arpae now adapted to handle these cases

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

from process.main_arpae import read_all_filenames_from_bucket
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



def main():

    # output path for processed files
    s3_bucket_new_name = "expats-radar-italy" # name of the destination data bucket where to store daily processed files of RR from italian domain

    # read all filenames from bucket
    s3, file_names = read_all_filenames_from_bucket()

    # list all years from 2024 to 2014
    yy_list = [str(year) for year in range(2023, 2016, -1)]
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

                        # convert unique time differences in minutes
                        unique_diffs_min = unique_diffs / np.timedelta64(1, 'm')

                        print("Unique time differences in the dataset in minutes:", unique_diffs)

                        # check at which minutes the data is available
                        minutes_available = [pd.to_datetime(str(t)).minute for t in ds_day.time.values]

                        unique_minutes = np.unique(minutes_available)
                        print("Unique minutes available in the dataset:", unique_minutes)

                        # if more than one time res is found or minutes are not 15, 30, 45, 00 , write date to log file for reprocessing
                        if len(unique_diffs) > 1 or ((unique_diffs_min != 5) and any(minute not in [0, 15, 30, 45] for minute in unique_minutes)):
                            print("Inconsistent time resolution detected. Logging date for reprocessing.")
                            with open('log_to_delete_and_reprocess_ARPAE.txt', 'a') as log_file:
                                log_file.write(f"{yyyy}-{mm}-{dd}\n")
                            continue 
                        else:
                            print('Time resolution is consistent. No action needed.')
                            print('*********************************************************')
                            # write date to another log file for info
                            with open('log_consistent_time_res_ARPAE.txt', 'a') as log_file:
                                log_file.write(f"{yyyy}-{mm}-{dd}\n")
                                
if __name__ == "__main__":
    main()  