
import boto3
import os
import logging
from botocore.exceptions import ClientError
from readers.config import *
import pdb

def list_files_bucket(s3, bucket_name):
    """
    script to list files in the bucket
    - s3: bucket to be checked
    
    """
    # Pagination to get all objects
    paginator = s3.get_paginator('list_objects_v2')
    pages = paginator.paginate(Bucket=bucket_name)

    print(pages)
    all_objects = []
    for page in pages:
        if "Contents" in page:
            for obj in page['Contents']:
                print(obj['Key'])  # Only print the Key
            all_objects.extend(page["Contents"])

    print(f"Total objects in bucket: {len(all_objects)}")

    #files present in the bucket
    all_objects_names = set([obj['Key'].split('/')[0] for obj in all_objects])

    return all_objects_names
    


def check_file_bucket(file2find, domain):
    """code to check if a file is on the destination bucket alredy processed, to avoid 
    reprocessing and reuploading the same file on the bucket.
    input:
    file2find: filename with path to find n bucket
    domain: string, 'DE' or 'IT' for the domain to check the bucket
    
    return:
    filefound (boolean) true if file is found, false if file is not on the bucket
    """
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file
    import boto3

    if domain == 'DE':
        from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
        bucket_name = "expats-radar-germany"    

    elif domain == 'IT':
        from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
        bucket_name = "expats-radar-italy"

    elif domain == 'CH':
        from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
        bucket_name = "radar-ch"
    else:
        print('domain not recognized')
        return(False)
    
    s3 = boto3.client(
        's3',
        endpoint_url=S3_ENDPOINT_URL,
        aws_access_key_id=S3_ACCESS_KEY,
        aws_secret_access_key=S3_SECRET_ACCESS_KEY)
    
    
    # List the objects in our bucket
    response = s3.list_objects(Bucket=bucket_name)

    # set flag to false if file is not found
    filefound = False
    if "Contents" not in response:
        print(f"No files on bucket")
        
    
    # loop on files found on bucket
    for obj in response["Contents"]:
        key = obj["Key"]
        if key == file2find:
            filefound = True
            return filefound
        else:
            continue

    return filefound





def upload_to_bucket(path_out, filename_ncdf, bucket_name):
    """
    function to upload on a specified S3 bucket the filename ncdf

    Args:
        path_out (string): path to the directory where filename_ncdf is located 
        filename_ncdf (string): name of the file to upload on bucket
        bucket_name (string): name of the bucket where to upload the file

    Returns:
        check: boolean, true if file is uploaded, false if not
    """
    
    import time
    from glob import glob

    from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file


    start_time = time.time()

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # calling upload function and return boolean for upload status
    # upload_file(s3_client, file_name, bucket, object_name=None):
    file_to_upload = os.path.join(path_out, filename_ncdf)
    check = upload_file(s3, file_to_upload, bucket_name, filename_ncdf)
    
    # if true, upload done
    if check:
        print("Time taken to upload files: ", time.time() - start_time, flush=True) 
        print('file name of file uploaded on bucket', filename_ncdf)

    return check
    

# method to initialize the S3 client
def Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY):
    """Initialize the S3 client
    :param S3_ENDPOINT_URL: S3 endpoint URL
    :param S3_ACCESS_KEY: S3 access key
    :param S3_SECRET_ACCESS_KEY: S3 secret access key
    :return: S3 client object
    """
    # Initialize the S3 client
    s3 = boto3.client(
        's3',
        endpoint_url=S3_ENDPOINT_URL,
        aws_access_key_id=S3_ACCESS_KEY,
        aws_secret_access_key=S3_SECRET_ACCESS_KEY
    )
    return s3


# methods for reading data
def read_file(s3, file_name, bucket):
    """reading a file from an S3 bucket
    :param s3: Initialized S3 client object
    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :return: object if file was uploaded, else False
    """
    try:
        #with open(file_name, "rb") as f:
        obj = s3.get_object(Bucket=bucket, Key=file_name)
        #print(obj)
        myObject = obj['Body'].read().decode('utf-8')
    except ClientError as e:
        logging.error(e)
        return None
    return myObject


def init_s3():
    """
    Initializes and returns an S3 client using the provided credentials and endpoint.
    """

    return boto3.client(
        's3',
        endpoint_url=S3_ENDPOINT_URL,
        aws_access_key_id=S3_ACCESS_KEY,
        aws_secret_access_key=S3_SECRET_ACCESS_KEY, 
    )

def delete_file_from_bucket(s3, bucket_name, file_key):
    """
    Deletes a file from the specified S3 bucket.

    Args:
        s3: Initialized S3 client object.
        bucket_name (str): Name of the S3 bucket.
        file_key (str): Key (path) of the file to delete.

    Returns:
        bool: True if deletion was successful, False otherwise.
    """
    try:
        s3.delete_object(Bucket=bucket_name, Key=file_key)
        print(f"File {file_key} deleted successfully from bucket {bucket_name}.")
        return True
    except ClientError as e:
        logging.error(e)
        return False
    
    

def read_file_obj(s3, file_path, bucket):

    # print all files in the bucket for debugging
    #all_files = print_list_files_in_bucket(s3, bucket)
    # example: /data/sat/msg/ml_train_crops/IR_108-WV_062-CMA_FULL_EXPATS_DOMAIN/2025/08/merged_MSG_CMSAF_2025-08-22.nc

    try:
        obj = s3.get_object(Bucket=bucket, Key=file_path)
        return obj['Body'].read()
    except ClientError as e:
        logging.warning(f"Failed to read file {file_path}: {e}")
        return None


def list_objects(s3, S3_BUCKET_NAME):
    # List the objects in our bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    items = []
    for item in response['Contents']:
        print(item['Key'])  
        items.append(item['Key'])
    return items
    

# method to upload data to the bucket
def upload_file(s3_client, file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = os.path.basename(file_name)
    try:
        with open(file_name, "rb") as f:    
            s3_client.upload_fileobj(f, bucket, object_name) 

        #response = s3_client.upload_file(file_name, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True

def download_file(s3_client, bucket, object_name, local_file):
    """Download a file from an S3 bucket

    input:
        - s3_client: s3 bucket initialized to download from
        - bucket: string with bucket name to download from
        - object_name: filename string to download
        - local_file: destination full path including filename
    :return: True if file was downloaded, else False
    """
    try:
        with open(local_file, "wb") as f:
            s3_client.download_fileobj(bucket, object_name, f)
        #s3_client.download_file(bucket, object_name, local_file)
    except ClientError as e:
        logging.error(e)
        return False
    return True


# function to download data from bucket
def download_from_s3(outpath, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL, S3_BUCKET_NAME):
    import os
    import boto3
    import logging
    from botocore.exceptions import ClientError
    from datetime import datetime
    import sys

    
    s3 = boto3.client(
        's3',
        endpoint_url=S3_ENDPOINT_URL,
        aws_access_key_id=S3_ACCESS_KEY,
        aws_secret_access_key=S3_SECRET_ACCESS_KEY
    )
    
    os.makedirs(outpath, exist_ok=True)
    
    years = [2023]
    months = range(4, 10)
    days = range(1, 32)
    for year in years:
        for month in months:
            for day in days:
                #prefix = f"output/data/timeseries_crops/{year:04d}/{month:02d}/{day:02d}/MSG_timeseries_{year:04d}-{month:02d}-{day:02d}_"
                prefix = f"{year}{month:02d}{day:02d}_RR_15min_msg_res.nc" #f"/data/sat/msg/ml_train_crops/IR_108-WV_062-CMA_FULL_EXPATS_DOMAIN/{year}/{month:02d}/merged_MSG_CMSAF_{year}-{month:02d}-{day:02d}.nc"
                try:
                    response = s3.list_objects_v2(
                        Bucket=S3_BUCKET_NAME,
                        Prefix=prefix
                    )
                    if "Contents" not in response:
                        print(f"No files for {year}-{month:02d}-{day:02d}")
                        continue
                    for obj in response["Contents"]:
                        key = obj["Key"]
                        if not key.endswith(".nc"):
                            continue
                        filename = os.path.basename(key)
                        local_file = os.path.join(outpath, filename)
                        if os.path.exists(local_file):
                            print(f"Already downloaded: {filename}")
                            continue
                        print(f"Downloading: {key}")
                        with open(local_file, "wb") as f:
                            s3.download_fileobj(S3_BUCKET_NAME, key, f)
                except ClientError as e:
                    print(f"Failed to list/download files for {year}-{month:02d}-{day:02d}: {e}")