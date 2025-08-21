"""
reading routines for radar data from DWD

"""

import xarray as xr
import numpy as np
import glob
from datetime import datetime
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from figures.domain_info import domain_expats
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg
from readers.file_dirs import orography_file

def read_radar_DWD(path_radolan_DE, day):
    """
    function to read the radar files of a given day
    Args:
        path_radolan_DE (_type_): _description_
    """
    yy = day[0:4]
    mm = day[4:6]
    filename = path_radolan_DE+yy+'/'+mm+'/YW_2017.002_'+day+'.nc'
    data = xr.open_dataset(filename)    
    
    return(data)
    
    

def read_orography():
    
    data = xr.open_dataset(orography_file)
    return data

