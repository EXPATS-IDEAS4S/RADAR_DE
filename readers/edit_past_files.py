"""
script to read files on the bucket, compress the variable and the gzip. then upload the gzip file and remove all the .nc files left on the bucket
"""

from readers.s3_buckets_credentials import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from readers.data_buckets_funcs import Initialize_s3_client, upload_file

import boto3
import glob
import xarray as xr 
import os
import shutil
import gzip


s3 = boto3.client(
's3',
endpoint_url=S3_ENDPOINT_URL,
aws_access_key_id=S3_ACCESS_KEY,
aws_secret_access_key=S3_SECRET_ACCESS_KEY)
    
from readers.data_buckets_funcs import list_files_bucket, check_file_bucket, upload_to_bucket, download_from_s3

def main():
        
    S3_BUCKET_NAME = 'expats-radar-germany'
    
    outpath = '/net/ostro/radolan_5min_composites/edit_old_files/'
    path_newfiles = '/net/ostro/radolan_5min_composites/edit_old_files/new_ncdf/'
    # download files to outpath
    download_from_s3(outpath, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL, S3_BUCKET_NAME)
    

    # list files on the outpath
    list_files = sorted(glob.glob(outpath+'*.nc'))
    for i_file, filename in enumerate(list_files):
        
        # read ncdf
        ds = xr.open_dataset(filename)
        
        # store ncdf with same name and compressed variable
        ds.to_netcdf(os.path.join(path_newfiles, os.path.basename(filename)), encoding={'RR': {'zlib': True, 'complevel': 9}})

        # gzip the file
        with open(os.path.join(path_newfiles, os.path.basename(filename)), 'rb') as f_in:
            with gzip.open(os.path.join(path_newfiles, os.path.basename(filename)+'.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

                
                # upload gzip file to bucket
                check = upload_to_bucket(path_newfiles, os.path.basename(filename)+'.gz')
    
    # remove all ncdf files from bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    for obj in response["Contents"]:
        if obj["Key"].endswith(".nc"):
            s3.delete_object(Bucket=S3_BUCKET_NAME, Key=obj["Key"])


if __name__ == "__main__":
    main()