"""
Docstring for readers.upload_swiss_to_bucket

code to upload swiss radar data .tar located on downloads directory
it then untars the file on the S3 bucket and deletes the .tar file from the bucket to save space,
 since the uncompressed data is much smaller than the original .tar file
 (70 GB) to EWC bucket

"""

import os
import time
from glob import glob

from progressbar import Bar
from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from readers.data_buckets_funcs import Initialize_s3_client, upload_file, check_file_bucket, list_files_bucket
import tarfile
import pdb


def main():


    # define bucket name 
    bucket_name = "radar-ch"

    # define start time to calculate upload time
    start_time = time.time()

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # list files on bucket and local directory
    files_bucket = list_files_bucket(s3, bucket_name)
    files_local = glob("/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/*.zip")

    # calculate the index to extract the filename from the path, since we need to upload only the filename on bucket and not the whole path
    ind_start = len("/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/")
    
    # loop to upload the files on bucket, check if file is already on bucket, if not, upload it
    for file in files_local:

        # check if file is in files_bucket, if not, upload it
        if file not in files_bucket:

            

            filename_ncdf = file[ind_start:]
            print('filename to upload', filename_ncdf)
            file_to_upload = os.path.join("/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/", filename_ncdf)

            # call upload function and return boolean for upload status
            check = upload_file(s3, file_to_upload, bucket_name, filename_ncdf)
            # if true, upload done
            if check:
                print("Time taken to upload files: ", time.time() - start_time, flush=True) 
                print('file name of file uploaded on bucket', filename_ncdf)

if __name__ == "__main__":
    main()

