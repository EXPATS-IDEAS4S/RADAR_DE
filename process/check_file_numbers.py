"""
code to check if the number of .zip files on the bucket is the 
same as the number of .zip files in the local directory, 
to check if all files have been uploaded on the bucket"""


import os
from glob import glob
from readers.data_buckets_funcs import list_files_bucket
from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from readers.data_buckets_funcs import Initialize_s3_client

def main():
    
    # define bucket name 
    bucket_name = "radar-ch"

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # list files on bucket and local directory
    files_bucket = list_files_bucket(s3, bucket_name)
    files_local = glob("/Users/claudia/Documents/Data/RADAR_SVIZZERA/hdf5_archive_talia/*.zip")

    # check if number of files on bucket is the same as number of files in local directory
    if len(files_bucket) == len(files_local):
        print("All files have been uploaded on the bucket")
    else:
        print("Not all files have been uploaded on the bucket")
        print("Number of files on bucket: ", len(files_bucket))
        print("Number of files in local directory: ", len(files_local))
    
if __name__ == "__main__":
    main()
