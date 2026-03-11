"""
code to check if the number of .zip files on the bucket is the 
same as the number of .zip files in the local directory, 
to check if all files have been uploaded on the bucket"""


import os
from glob import glob
from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from readers.data_buckets_funcs import Initialize_s3_client, list_files_bucket
import pdb

def main():
    
    # define bucket name 
    bucket_name = "expats-radar-germany"

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # list files on bucket and local directory
    files_bucket = list_files_bucket(s3, bucket_name)

    # check if 20200401_RR_15min_msg_res.nc.gz is in the list of filenames
    if "20200401_RR_15min_msg_res.nc.gz" in files_bucket:
        print("file 20200401_RR_15min_msg_res.nc.gz is in the list of filenames on bucket")
    else:
        print("file 20200401_RR_15min_msg_res.nc.gz is not in the list of filenames on bucket")
    pdb.set_trace()
    
if __name__ == "__main__":
    main()
