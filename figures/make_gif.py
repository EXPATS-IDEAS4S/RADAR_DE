"""
code to generate gif when a folder containaing images is given
"""

from readers.file_dirs import path_radolan_DE, path_out, path_arpae_DE, path_nc
import os
from figures.one_day_video import gif_maker

def main(mode='radar', date='20150522'):

    # define a dictionary with parameters for thr gif-maker function
    params_radar = {
        "image_folder": path_out,
        "gif_name": f"{date}_radar",
        "gif_path": path_out,
        "gif_duration": 100,  # duration for gif in ms
        "channel": "rainrate"
    }

    params_spacetime_crops = {
        "image_folder": "/Users/claudia/Documents/Data/space_time_crops/",
        "gif_name": f"{date}_spacetime_crops",
        "gif_path": "/Users/claudia/Documents/Data/space_time_crops/",
        "gif_duration": 100,  # duration for gif in ms
        "channel": "rainrate"
    }

    if mode == 'radar':
        dict_plot = params_radar
    elif mode == 'spacetime_crops':
        dict_plot = params_spacetime_crops

    
    path_images = dict_plot["image_folder"]
    gif_name = dict_plot["gif_name"]
    gif_name = dict_plot["gif_name"]
    gif_path = dict_plot["gif_path"]
    gif_duration = dict_plot["gif_duration"]
    channel = dict_plot["channel"]

    gif_maker(image_folder=path_images,
              gif_name=gif_name,
              gif_path=gif_path,
              gif_duration=gif_duration,
              channel=channel)



if __name__ == "__main__":
    main('radar', '20140615')

