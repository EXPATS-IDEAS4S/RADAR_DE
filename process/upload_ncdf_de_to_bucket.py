"""
code to upload to expats-radar-germany all the ncdf files in the folder that are not 
already on the bucket
"""




from glob import glob
import os
import pdb
import boto3
from botocore.exceptions import ClientError
import logging
# Initializing variables for the client
#S3_BUCKET_NAME = "radar-ch"#"arpae-radar-composite" #"expats-radar-germany"  #Fill this in l
S3_ACCESS_KEY = "8aef2cea93fa46a28ebc43a906a5f2ce"  #Fill this in 
S3_SECRET_ACCESS_KEY = "b84d267a7c0a483985b320d61714d805"  #Fill this in 
S3_ENDPOINT_URL = "https://s3.waw3-1.cloudferro.com"  #Fill this in

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



def list_files_bucket(s3, s3_bucket_name):
    """
    script to list files in the bucket
    - s3: bucket to be checked
    
    """
    # List the objects in our bucket
    response = s3.list_objects(Bucket=s3_bucket_name)
    
    # set flag to false if file is not found
    filefound = False
    if "Contents" not in response:
        print(f"No files on bucket")
        return([])
    else:
        return(response["Contents"])
    
def main():

    # reading list of files from folder
    files_ncdf_folder = sorted(glob("/data/trade_pc/radolan_DE_5min_rain_rate/ncdf/*.nc.gz"))

    # initialize the S3 client to upload the data to bucket
    s3 = Initialize_s3_client(S3_ENDPOINT_URL, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY)
    bucket_name = "expats-radar-germany"

    # list files on bucket and local directory
    files_bucket = list_files_bucket(s3, bucket_name)

    # store filename lists by reading the files on bucket and in local directory
    filenames = [obj['Key'] for obj in files_bucket]
    filenames = sorted(filenames)
    print("total number of files on bucket: ", len(filenames))  

    # find files that are in the local directory but not on the bucket
    files_to_upload = [f for f in files_ncdf_folder if os.path.basename(f) not in filenames]
    print(f"files to upload: {files_to_upload}")
    pdb.set_trace() 

    # loop to upload the files on bucket
    for file in files_to_upload:
        filename = os.path.basename(file)
        print(f"uploading file {filename} to bucket {bucket_name}")
        upload_file(s3, file, bucket_name, filename)

    # check that now files are on bucket
    files_bucket = list_files_bucket(s3, bucket_name)
    filenames = [obj['Key'] for obj in files_bucket]
    filenames = sorted(filenames)
    print("total number of files on bucket after uploading: ", len(filenames))  
    pdb.set_trace()

if __name__ == "__main__":
    main()


