
import glob
from scipy.interpolate import griddata
import xarray as xr
from readers.file_dirs import path_radolan_DE
from readers.radar_DWD import read_radar_DWD, read_orography
import pdb
import numpy as np
import scipy
from readers.data_buckets_funcs import download_from_s3
from figures.domain_info import domain_expats, domain_DE_CA
from figures.mpl_style import plot_cities_expats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from figures.one_day_video import plot_radar_map
import matplotlib.pyplot as plt
import os
from time import sleep
from progress.bar import Bar
import time
import pandas as pd
import shutil
import gzip

def plot_msg(lons, lats, variable, label, cmap, cbar_title, domain, path_out, key, back_transparent, vmin, vmax):
    """
    plot map of the variable over the defined domain

    Args:
        lons (array): longitude values
        lats (array): latitudes 
        variable (matrix): variable to plot
        label (string): string for plot filename
        cmap (colormap): color map 
        cbar_title (string): title for color bar
        domain (array): minlon, maxlon, minlat, maxlat     

        path_out (string): path where to save the plot
        key (string): string possible "expats" or "dfg"
        back_transparent (boolean): True or False - true means transparent background
        vmin (float): value min in colorbar
        vmax (float): max value in colorbar
        
    Dependencies:
    plot_cities_expats, 
    plot_local_dfg
    
    """
    from figures.mpl_style import plot_cities_expats

    # Plot the map of MSG channels vs lat/lon
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.spines["top"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.set_extent(domain)
    
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 14}
    gl.ylabel_style = {'fontsize': 14}
    
    

    var_levels = np.linspace(vmin, vmax, 20)
    # plot 10th percentile as filled contours
    mesh = ax.contourf(lons, 
                        lats, 
                        variable, 
                        cmap=cmap, 
                        transform=ccrs.PlateCarree(), 
                        levels=var_levels, 
                        vmin=vmin, # 230
                        vmax=vmax) # 270
    cbar = plt.colorbar(mesh, label=cbar_title, shrink=0.6)

    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    
     #vmin, vmax=220.,
    plt.tick_params(axis='both',which='major',labelsize=14)

    ax.set_xlabel('Latitude [$^{\circ}$]')
    ax.set_ylabel('Longitude [$^{\circ}$]')
    ax.tick_params(which='minor', length=5, width=2)
    ax.tick_params(which='major', length=7, width=3)
    #cbar = plt.colorbar(pc,ax=ax,shrink=0.75)
    #cbar.set_label(label, fontsize=14)
    #cbar.ax.tick_params(labelsize=14)

    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=1.5, color='black')
    #ax.add_feature(cfeature.STATES, linewidth=0.2)
    ax.add_feature(cfeature.BORDERS, linewidth=1.5, color='black')

    # Title and save date_str+'_IR108'
    ax.set_title(label, fontsize=12)
    
    if key =='expats':
        plot_cities_expats(ax, 'grey', 50)

        
        
    plt.savefig(
        os.path.join(path_out, label+"_"+key+".png"),
        dpi=300,
        bbox_inches="tight",
        transparent=back_transparent,
        )

    plt.close()
#f"{i:03d}"

    
    
    

def plot_radar_mapV2(data, string_file, path_out):
    """
    function to plot map of radar data from DWD, in rain rate

    Args:
        data (xarray dataset): data from DWD 
    """
    print(np.nanmax(data.RR.values))
    print(np.nanmin(data.RR.values))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.spines["top"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.set_extent(domain_DE_CA)
    
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 14}
    gl.ylabel_style = {'fontsize': 14}
    

    vmin= 0.
    vmax= 5.
    var_levels = np.linspace(vmin, vmax, 15)
    # plot rain rate as filled contours
    lats = data.lat.values
    lons = data.lon.values
    RR = data.RR.values
    
    RR[RR==0.] = np.nan
    # reading orography data from raster file
    ds_or = read_orography()
    oro_levels = np.linspace(0, 1500, 20)
    oro = ax.contourf(ds_or.lons.values, 
                        ds_or.lats.values, 
                        ds_or.orography.values, 
                        transform=ccrs.PlateCarree(), 
                        levels=oro_levels, 
                        alpha=1.,
                        cmap='Greys')

    mesh_rr = ax.contourf(lons, 
                        lats, 
                        RR, 
                        cmap='BuPu_r', 
                        transform=ccrs.PlateCarree(), 
                        vmin=vmin,  
                        vmax=vmax, 
                        levels=var_levels, 
                        extend='max',
                        alpha=0.6)  
        
    cbar = plt.colorbar(mesh_rr, label='Rain rate [mm]', shrink=0.6)
    plot_cities_expats(ax, 'black', 50)
    ax.add_feature(cfeature.BORDERS, linewidth=1., color='black')
       
    plt.savefig(
        os.path.join(path_out, string_file+'.png'),
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        )

    plt.close()
    return()
