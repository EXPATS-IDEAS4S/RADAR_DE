"""
This module contains code to read and 
crop MSG data to the domain of interest.
"""
import xarray as xr
from figures.domain_info import domain_expats, domain_DE_CA

def prepare_msg_data():
    
    """
    code to read and crop msg data to the domain of interest
    returns the cropped dataset to be used for interpolation of 
    radar data on the same grid
    """
    ds_msg = xr.open_dataset('/home/claudia/codes/RADAR_DE/20220819-EXPATS-RG.nc')

    # crop msg data to domain_DE_CA
    lon_min, lon_max, lat_min, lat_max = domain_expats

    #print(f'Cropping MSG data to domain: {lat_min}, {lat_max}, {lon_min}, {lon_max}')
    
    # crop msg data to the DE domain (smaller than EXPATS)
    ds_msg_DE = ds_msg.where((ds_msg.lat >= lat_min) & 
                            (ds_msg.lat <= lat_max) & 
                            (ds_msg.lon >= lon_min) & 
                            (ds_msg.lon <= lon_max), 
                            drop=True)

    return ds_msg_DE