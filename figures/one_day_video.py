"""
Code to read and plot the video of 1 dday radar images of rain rates
"""



from matplotlib.colors import ListedColormap
from readers.file_dirs import path_radolan_DE, path_out, path_arpae_DE, path_nc
from readers.radar_DWD import read_radar_DWD, read_orography
import xarray as xr
import numpy as np
import glob
from datetime import datetime
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from figures.domain_info import domain_expats, domain_DE_CA, domain_IT_CA
from figures.mpl_style import CMAP, plot_cities_expats, plot_local_dfg, plot_cities_IT
import os
import cartopy.feature as cfeature
import gzip
import shutil
from geopy.distance import geodesic


def main():
    
    
    # set the day to process
    yy = '2018'
    mm = '04'
    dd = '09'
    domain = 'IT'  # 'DE' or 'expats'
    date = yy+mm+dd
    
    if domain == 'DE':
        data = read_radar_DWD(path_radolan_DE, date)    
    elif domain == 'IT':
        # dezip from gz file

        with gzip.open(f'{path_nc}{yy}{mm}{dd}_RR_IT_15min_msg_res.nc.gz', 'rb') as f_in:
            with open(f'{path_nc}{yy}{mm}{dd}_RR_IT_15min_msg_res.nc', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        ncfile = f'{path_nc}{yy}{mm}{dd}_RR_IT_15min_msg_res.nc'
        # unzip and read file 20230401_RR_IT_15min_msg_res.nc.gz
        data = xr.open_dataset(ncfile, 
                              engine='netcdf4', 
                              decode_times=True)    
    # read data
    time_dim = len(data.time.values)
    print(np.nanmax(data.RR.values))
    
    for i_time in range(time_dim):
        
        data_day = data.isel(time=i_time)
        # extract time string of the form yymmdd_HHMM
        time_string = str(data_day.time.values).replace("-", "").replace(":", "").replace("T", "_")[:13]        


        plot_radar_map(data_day, time_string, domain_IT_CA, domain_name='IT')

        #plot_RR_distribution(data_day, time_string, domain_name='IT')
        
        # calculate RR percentiles
        percentiles = np.nanpercentile(data_day.RR.values.flatten(), [95, 96, 97, 98, 99])
        print(f"Rain rate percentiles at {time_string} UTC:")
        for p, val in zip([95, 96, 97, 98, 99], percentiles):
            print(f"  {p}th percentile: {val:.2f} mm/h")    

        print(time_string)

    # create animated gif video
    gif_maker(path_out, f'{yy}{mm}{dd}_radar', path_out, 60, 'rainrate')

def plot_RR_distribution(data, time_string, domain_name='IT'):
    """
    function to plot histogram of rain rate distribution

    Args:
        data (xarray dataset): data from DWD 
    """
    RR = data.RR.values.flatten()

    RR = RR[~np.isnan(RR)]

    
    plt.figure(figsize=(8,6))
    plt.hist(RR, bins=50, range=(0,1), color='blue', alpha=0.7)
    plt.xlabel('Rain rate [mm/h]', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(f'Rain rate distribution at {time_string} UTC', fontsize=16)
    plt.grid()
    plt.savefig(
        os.path.join(path_out, time_string+"_rain_rate_distribution.png"),
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        )
    # save figure
    plt.savefig(f'{path_out}{time_string}_rain_rate_distribution.png', dpi=300) 
    plt.close()
    return()


def plot_radar_map(data, time_string, domain, domain_name='IT'):
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
    ax.set_extent(domain)
    
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontsize': 14}
    gl.ylabel_style = {'fontsize': 14}
    

    # set levels for rain rate from 0 to 20 in steps of 2 mm/h  

    vmin= 0.01
    vmax= 10.
    var_levels = np.arange(0., vmax+0.5, 0.5)
    # define a colormap for rain rate that does not have color for zero rain rate
    cmap_rr = plt.get_cmap('BuPu_r')
    #colors = cmap_rr(np.linspace(0, 1, cmap_rr.N))  
    #colors[0] = [1, 1, 1, 0]  # RGBA, alpha=0 for 0 value
    #custom_cmap = ListedColormap(colors)

    # plot rain rate as filled contours
    lats = data.lat.values
    lons = data.lon.values
    RR = data.RR.values
    
    # set rain rate values < 0.01 to nan
    RR[RR<0.01] = np.nan


    # reading orography data from raster file
    ds_or = read_orography()
    oro_levels = np.linspace(0, 2000, 50)
    oro = ax.contourf(ds_or.lons.values, 
                        ds_or.lats.values, 
                        ds_or.orography.values, 
                        transform=ccrs.PlateCarree(), 
                        levels=oro_levels, 
                        alpha=0.3,
                        cmap='Greys')

    mesh_rr = ax.contourf(lons, 
                        lats, 
                        RR, 
                        cmap=cmap_rr, 
                        transform=ccrs.PlateCarree(), 
                        vmin=vmin,  
                        vmax=vmax, 
                        levels=var_levels, 
                        extend='max',
                        alpha=0.8)  
    
    cbar = plt.colorbar(mesh_rr, label='Rain rate [mm]', shrink=0.6)
    # set colorbar ticks every 2 mm/h
    cbar.set_ticks(np.arange(0., vmax+0.5, 0.5))
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Accumulated Rain rate [mm]', fontsize=14)
    if domain_name == 'DE':
        plot_cities_expats(ax, 'black', 50)
    elif domain_name == 'IT':
        plot_cities_IT(ax, 'black', 50)
        plot_ARPAE_radar_ranges(ax)

    ax.add_feature(cfeature.BORDERS, linewidth=1., color='black')
    ax.add_feature(cfeature.COASTLINE, linewidth=1., color='black')
    

    plt.savefig(
        os.path.join(path_out, time_string+"_rain_rate.png"),
        dpi=300,
        bbox_inches="tight",
        transparent=True,
        )

    plt.close()
    return()

### TO DO: 1) change color bar to one with white background, 2) change units 3) resample on satellite time resolution 4) plot together satellite and radar data input

def plot_ARPAE_radar_ranges(ax):
    """
    Docstring for plot_ARPAE_radar_ranges
    
    :param ax: axis from matplotlib to plot the radar ranges
    :return: None   
    """

    radar_names = ['San Pietro Capofiume', 'Gattatico']
    radar_lats = [44.6547, 44.7914]
    radar_lons = [11.6236, 10.4992]
    radius_km = 212 # km

    # Bearings
    angles = np.linspace(0, 360, 361)


    # calculate and plot circles for each radar
    for i in range(len(radar_names)):
        radar_name = radar_names[i]
        lat0 = radar_lats[i]
        lon0 = radar_lons[i]

        lats = []
        lons = []
        for a in angles:
            destination = geodesic(kilometers=radius_km).destination((lat0, lon0), a)
            lats.append(destination.latitude)
            lons.append(destination.longitude)

        # plot circle lines
        ax.plot(lons, lats, color='white', linewidth=1, linestyle='--', transform=ccrs.PlateCarree())
        ax.scatter(radar_lons[i], radar_lats[i], color='white', s=50, marker='^', transform=ccrs.PlateCarree())
        
    return

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
    