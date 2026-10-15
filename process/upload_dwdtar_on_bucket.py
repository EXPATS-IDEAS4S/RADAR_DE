"""
Code to move DWD tar files on dedicated data bucket for further processing

author: Claudia Acquistapace
date: 16 Febbruary 2026
email: claudia.acquistapace@unipd.it

to exectue
 source /Users/claudia/Github/RADAR_DE/venv/bin/activate
python3 -m process.upload_dwdtar_on_bucket

"""


from glob import glob
import os
from readers.config import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from readers.data_buckets_funcs import init_s3, list_objects, Initialize_s3_client, upload_file, check_file_bucket, list_files_bucket, upload_to_bucket
from tqdm import tqdm

import pdb

def main():


    data_bucket_name = "dwd-tar-files" # name of the destination data bucket where to store daily processed files of RR from italian domain
    s3 = init_s3()

    # read filename of files on the bucket
    files_on_bucket = list_objects(s3, data_bucket_name)
    print('files on bucket', files_on_bucket)

    print(f"moving DWD tar files on bucket {data_bucket_name}")
    
    # tar folder path for DWD radar data
    #DWD_archive_path = "/Volumes/Matteo_USB_/radar_germany/"
    DWD_archive_path = "/data/trade_pc/radolan_DE_5min_rain_rate/"

    file_list = glob(os.path.join(DWD_archive_path, "*.tar.gz"))
    print(" files to upload ", file_list)

    # loop to move the files on bucket
    for file in file_list:

        # read only filename from the path and check if it is already on bucket. skip if yes    
        filename = file.split("/")[-1]
        print(filename)

        # check if file is in files_on_bucket. skip if yes
        if filename in files_on_bucket:
            print(f"file {filename} already on bucket {data_bucket_name}. skipping")
            continue
        else:               
            # reading filename full path
            filename = os.path.basename(file)
            print(f"moving file {filename} on bucket {data_bucket_name}")
            print("******************************************************")
            check = upload_to_bucket(DWD_archive_path, filename, data_bucket_name)
            if check:
                print(f"file {filename} moved on bucket {data_bucket_name}")
            else:
                print(f"file {filename} not moved on bucket {data_bucket_name}")
            print("******************************************************")
            print("******************************************************")



if __name__ == "__main__":
    main()
