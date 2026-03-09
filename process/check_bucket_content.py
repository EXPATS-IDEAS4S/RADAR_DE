
import boto3
from readers.s3_bucket_credentials import S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL
from botocore.exceptions import ClientError
from readers.data_buckets_funcs import list_files_bucket, check_file_bucket, init_s3

"""
code to list all files in a data bucket, and print the list to a log file
to run the code, activate the venv and run the script with python3 -m process.check_bucket_content

"""


def main():
    
    bucket_name = "expats-radar-swisse"

    # initialize S3 client
    s3 = init_s3()

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
    all_objects_years = set([obj['Key'].split('/')[0] for obj in all_objects])
    print(f"Years present in bucket: {sorted(all_objects_years)}")

    # store list to log file
    with open(f"log_{bucket_name}.txt", "w") as f:
        f.write(f"Total objects in bucket: {len(all_objects_years)}\n")
        f.write("List of all files:\n")
        for obj in all_objects:
            f.write(f"{obj['Key']}\n")
    

if __name__ == "__main__":

    main()
