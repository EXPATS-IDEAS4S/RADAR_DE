"""
Code to access data on the Arcodatahub S3 bucket, 
which is used to store the data of the radar italian network.

activate virtual environment .arco_env before running this code with
source /home/claudia/.arco_env/bin/activate

call the script with 
python3 -m process.access_arcodatahub_data


"""
import xarray as xr
import requests
from readers.arcohub_credentials import username, access_key


def main():


    dataset_url = f"https://{username}:{access_key}@api.arcodatahub.com/S3/italian-radar-dpc-sri.zarr"    
    ds = xr.open_dataset(dataset_url, engine="zarr")

    print(ds.info())


if __name__ == "__main__":
    
    main()
