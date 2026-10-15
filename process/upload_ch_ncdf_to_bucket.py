"""
This is a code that uploads the ncdf radar ch files on the corresponding bucket.
it works that it reads the files on the folder and the files on the bucket, and then uploads the files that are in the folder but not in the bucket.
author: Claudia Acquistapace
date: 16 Febbruary 2026
email: claudia.acquistapace-at-unipd.it

how to run the code:
activate the venv and run the script with python3 -m process.upload_ch_ncdf_to_bucket

"""
import os
from glob import glob
from readers.data_buckets_funcs import init_s3, list_objects, upload_to_bucket

import pdb

def main():

    # reading list of files from folder
    five_min_path = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/5_minutes/"
    fifteen_min_path = "/Users/claudia/Documents/Data/RADAR_SVIZZERA/netcdf/15_minutes/"
    data_paths = [fifteen_min_path, five_min_path]
    # data bucket names
    five_min_bucket_name = "radar-ch"
    fifteen_min_bucket_name = "expats-radar-swisse"
    bucket_names = [fifteen_min_bucket_name, five_min_bucket_name]

    # list all files  in the folder (these are the files with 32 bits floating point data)  
    file_strings = [ "_RR_CH_15min_msg_res.nc.gz", "RR_CH_5min_msg_res.nc.gz"]
    print("******************************************************") 

    # loop on file types and bucket of destination
    for data_path, bucket_name, file_string in zip(data_paths, bucket_names, file_strings):

        file_to_upload = sorted(glob(os.path.join(data_path, "*"+file_string)))
        print(f"files in the folder {data_path}: {file_to_upload}")
        print(f"number of files in the folder {data_path}: {len(file_to_upload)}")
        print("******************************************************")
        pdb.set_trace()
        # read list of files on the bucket
        s3 = init_s3()
        files_on_bucket = list_objects(s3, bucket_name)
        print('files on bucket', files_on_bucket)

        for file in file_to_upload:

            # read only filename from the path and check if it is already on bucket. skip if yes    
            filename = file.split("/")[-1]

            file_path = os.path.dirname(file)
            print(f"processing file {filename} to upload on bucket {bucket_name}")

            # check if file is in files_on_bucket. skip if yes
            if filename in files_on_bucket:
                print(f"file {filename} already on bucket {bucket_name}. skipping")
                continue
            else:               
                # reading filename full path
                filename = os.path.basename(file)
                print(f"moving file {filename} on bucket {bucket_name}")
                print("******************************************************")
                check = upload_to_bucket(file_path, filename, bucket_name)
                if check:
                    print(f"file {filename} moved on bucket {bucket_name}")
                else:
                    print(f"file {filename} not moved on bucket {bucket_name}")
                print("******************************************************")
                print("******************************************************")


if __name__ == "__main__":
    main()