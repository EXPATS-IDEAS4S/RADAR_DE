"""
code to generate gif when a folder containaing images is given
"""

from readers.file_dirs import path_radolan_DE, path_out, path_arpae_DE, path_nc
import os
from figures.one_day_video import gif_maker

def main():

    date = '20230514'
    gif_name = f"{date}_radar"
    gif_path = path_out
    gif_duration = 100  # duration for gif in ms
    channel = "rainrate"
    gif_maker(image_folder=path_out,
              gif_name=gif_name,
              gif_path=gif_path,
              gif_duration=gif_duration,
              channel=channel)



if __name__ == "__main__":
    main()

