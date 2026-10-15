"""
code to generate gif when a folder containaing images is given
"""

from process.create_gifs_mp4 import convert_gif_to_mp4
from readers.file_dirs import path_radolan_DE, path_out, path_arpae_DE, path_nc
import os
from figures.one_day_video import gif_maker
from process.create_gifs_mp4 import convert_gif_to_mp4
import ffmpeg

def main(mode='radar_IT', date='20150522'):


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

    params_radar_CH = {
        "image_folder": "/Users/claudia/Documents/Data/test_radar/png/",
        "gif_name": f"{date}_radar_CH",
        "gif_path": "/Users/claudia/Documents/Data/test_radar/png/",
        "gif_duration": 100,  # duration for gif in ms
        "channel": "rainrate"   
    }


    if mode == 'radar_IT':
        dict_plot = params_radar
    elif mode == 'spacetime_crops':
        dict_plot = params_spacetime_crops
    elif mode == 'radar_CH':
        dict_plot = params_radar_CH
    else:
        raise ValueError(f"mode {mode} not recognized, choose between 'radar_IT', 'radar_CH' and 'spacetime_crops'")

    
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

    # convert gif to mpf using ffmpeg, to reduce the file size and make it easier to visualize
    video_path_mp4 = gif_path.replace('.gif', '.mp4')
    os.system(f'ffmpeg -i {gif_path} -vcodec libx264 -pix_fmt yuv420p {video_path_mp4}')
    print(f'Video quicklook saved at: {video_path_mp4}')

if __name__ == "__main__":
    main('radar_CH', '20240621')

