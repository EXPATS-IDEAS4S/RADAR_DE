import requests
import s3fs
import xarray as xr
from figures.domain_info import domain_expats


def read_arcohub_data(yy, mm, dd):

    """
    function to read data from the Arcodatahub S3 bucket, 
    which is used to store the data of the radar italian network. 
    it selects the specific date to process based on the input
    parameters yy, mm, dd. It crops based on the domain of interest and returns the dataset.

    input: yy, mm, dd as strings

    output: xarray dataset

    """
    from readers.arcohub_credentials import username, access_key    

    # Create an S3 filesystem object using the credentials
    s3 = s3fs.S3FileSystem(key=username, secret=access_key)

    # Example: Access a dataset
    dataset_url = f"https://{username}:{access_key}@api.arcodatahub.com/S3/italian-radar-dpc-sri.zarr"
    ds = xr.open_dataset(dataset_url, engine="zarr")

    # select the specific date to process based on the input parameters yy, mm, dd
    year = int(yy)
    month = int(mm)
    day = int(dd)
    start_date = datetime(year, month, day)
    end_date = datetime(year, month, day, 23, 59, 59)
    ds_day = ds.sel(time=slice(start_date, end_date))
    # delete the original dataset to free up memory
    del ds

    # Remove _FillValue attribute from datetime64 variables to avoid NetCDF/Zarr error
    for var in ds_day.variables:
        if str(ds_day[var].dtype).startswith("datetime64"):
            if '_FillValue' in ds_day[var].encoding:
                del ds_day[var].encoding['_FillValue']

    # crop over domain of interest
    ds_day_cropped = ds_day.where((ds_day.lon >= domain_expats[0]) & (ds_day.lon <= domain_expats[1]) &
                                 (ds_day.lat >= domain_expats[2]) & (ds_day.lat <= domain_expats[3]), drop=True)    
    
    return ds_day_cropped
