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

    

def read_orography():
    
    data = xr.open_dataset(orography_file)
    return data

