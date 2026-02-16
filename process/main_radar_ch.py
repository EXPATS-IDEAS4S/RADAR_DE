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
from readers.config import S3_BUCKET_NAME
from readers.data_buckets_funcs import upload_to_bucket
import boto3
import os
import tarfile
import time
from glob import glob

def main():

    # print
    print(f"uploading CH radar data on bucket {S3_BUCKET_NAME}")
    radar_data_file = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/"

    # upload the .tar file on bucket    
    check = upload_to_bucket(radar_data_file, "hdf5_archive_talia_zip.tar")
    if check:
        print("File uploaded on bucket")
    else:
        print("File not uploaded on bucket")


if __name__ == "__main__":
    main()