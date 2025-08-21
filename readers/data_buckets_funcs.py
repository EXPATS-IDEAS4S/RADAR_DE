
import boto3
import os
import logging
from botocore.exceptions import ClientError

def list_files_bucket(s3, S3_BUCKET_NAME):
    """
    script to list files in the bucket
    - s3: bucket to be checked
    
    """
    # List the objects in our bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    
    # set flag to false if file is not found
    filefound = False
    if "Contents" not in response:
        print(f"No files on bucket")
        return([])
    else:
        return(response["Contents"])
    


def check_file_bucket(file2find):
    """code to check if a file is on the bucket
    input:
    file2find: filename with path to find n bucket
    return:
    filefound (boolean) true if file is found, false if file is not on the bucket
    """
    
    from readers.s3_buckets_credentials import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file
    
    import boto3
    
    s3 = boto3.client(
    's3',
    endpoint_url=S3_ENDPOINT_URL,
    aws_access_key_id=S3_ACCESS_KEY,
    aws_secret_access_key=S3_SECRET_ACCESS_KEY)
    
    # List the objects in our bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    
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





def upload_to_bucket(path_out, filename_ncdf):
    """
    function to upload on a specified S3 bucket the filename ncdf

    Args:
        path_out (string): path to the directory where filename_ncdf is located 
        filename_ncdf (string): name of the file to upload on bucket

    Returns:
        check: boolean, true if file is uploaded, false if not
    """
    
    import time
    from glob import glob

    from readers.s3_buckets_credentials import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
    from readers.data_buckets_funcs import Initialize_s3_client, upload_file


    start_time = time.time()

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)

    # calling upload function and return boolean for upload status
    # upload_file(s3_client, file_name, bucket, object_name=None):
    file_to_upload = os.path.join(path_out, filename_ncdf)
    check = upload_file(s3, file_to_upload, S3_BUCKET_NAME, filename_ncdf)
    
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

def list_objects(s3, S3_BUCKET_NAME):
    # List the objects in our bucket
    response = s3.list_objects(Bucket=S3_BUCKET_NAME)
    for item in response['Contents']:
        print(item['Key'])
    

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