"""
This code reads and plots radar data from ARPAE 

command to install packages: 
/home/cacquist/Documents/GitHub/EXPATS/.venv/bin/python -m pip install boto3
command to run this script:
/home/cacquist/Documents/GitHub/EXPATS/.venv/bin/python process/main_arpae.py
"""

from readers.s3_buckets_credentials import S3_ACCESS_KEY, S3_ENDPOINT_URL, S3_SECRET_ACCESS_KEY
import xarray as xr
import numpy as np
import glob
from datetime import datetime
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from figures.domain_info import domain_expats
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg
from readers.file_dirs import orography_file
import boto3
import pdb

def read_radar_ARPAE(path_arpae_DE, yyyy, mm, dd):
    """
    function to read the list of radar files of a given day
    Args:
        path_arpae_DE (_type_): _description_
    """

    
    
    s3 = boto3.client(
    's3',
    endpoint_url=S3_ENDPOINT_URL,
    aws_access_key_id=S3_ACCESS_KEY,
    aws_secret_access_key=S3_SECRET_ACCESS_KEY)

    # List the objects in our bucket
    response = s3.list_objects(Bucket="arpae-radar-composite")

    # set flag to false if file is not found
    filefound = False
    if "Contents" not in response:
        print(f"No files on bucket")
        
    
    # loop on files found on bucket to produce file list
    data = []
    for obj in response["Contents"]:
        key = obj["Key"]

        # append only filenames for the specified day
        if (yyyy in  key) and (mm in key) and (dd in key):
            data.append(key)
        else:
            continue
    return(data)


def main():


    yyyy = "2013"
    mm = "08"
    dd = "22"

    path_ARPAE = "/home/vpoli@ARPA.EMR.NET/dati_claudia/composito/"
    data = read_radar_ARPAE(path_ARPAE, yyyy, mm, dd)
    print(data) 
    print(len(data))

if __name__ == "__main__":
    main()