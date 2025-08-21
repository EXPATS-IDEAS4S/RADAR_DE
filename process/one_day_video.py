"""
Code to read and plot the video of 1 dday radar images of rain rates
"""



from readers.file_dirs import path_radolan_DE, path_out
from readers.radar_DWD import read_radar_DWD, read_orography
import xarray as xr
import numpy as np
import glob
from datetime import datetime
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from figures.domain_info import domain_expats, domain_DE_CA
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg
import os
import cartopy.feature as cfeature


def main():
    
    
    # set the day to process
    yy = '2022'
    mm = '06'
    dd = '05'
    date = yy+mm+dd
    
    # read data
    data = read_radar_DWD(path_radolan_DE, date)    

    time_dim = len(data.time.values)
    print(np.nanmax(data.RR.values))
    
    for i_time in range(time_dim):
        
        data_day = data.isel(time=i_time)
        time_string = str(data_day.time.values)
        plot_radar_map(data_day, time_string)
   
        print(time_string)

    # create animated gif video
    gif_maker(path_out, '20220605_radar', path_out, 60, 'rainrate')

def plot_radar_map(data, time_string):
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
    vmax= 10.
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
        os.path.join(path_out, time_string+"_rain_rate.png"),
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        )

    plt.close()
    return()

### TO DO: 1) change color bar to one with white background, 2) change units 3) resample on satellite time resolution 4) plot together satellite and radar data input
    
def gif_maker(image_folder, gif_name, gif_path, gif_duration, channel):
    """
    script to create animated gif from a folder containing images

    Args:
        image_folder (string): folder containing png images
        gif_name (string): string as filename for gif
        gif_path (string): path for gif file
        gif_duration (int): duration for gif (typical 250)
        channel(string): variable string for the gif
        
    """
    from PIL import Image
    import glob
    
    import matplotlib.pyplot as plt
    

    
    # read files into the fil elist
    image_array = []
    print(sorted(glob.glob(path_out+'*.png')))
    for file in sorted(glob.glob(path_out+'*.png')):
                
        image = Image.open(file)
        image_array.append(image)

    im = image_array[0]            
    im.save(gif_path+gif_name+".gif", 
            format='png',
            save_all=True, 
            append_images=image_array, 
            duration=gif_duration, 
            loop=0)
    

    
    return



if __name__ == "__main__":
    main()
    