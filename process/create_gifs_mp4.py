"""
script to read lists of .png images, create gifs and then convert to mp4
"""
import pandas as pd


def main():
    

    # set the day to plot
    yy = '2022'
    mm = '08'
    dd = '19'
    date = yy+mm+dd
    path_out = '/net/ostro/radolan_5min_composites/test_1_day_resample_msg/'

    # create a time array of 15 minutes resolution for the selected day
    time_array = pd.date_range(start=f'{yy}-{mm}-{dd} 00:00', end=f'{yy}-{mm}-{dd} 23:59', freq='15T')
    
    list_RR_msg = []
    list_RR_orig = []
    list_msg = []
    # create list of .png for each produced plot
    for i_t, time_val in enumerate(time_array):
        
        # generate string of the date with hh and minute
        time_string = pd.Timestamp(time_val).strftime('%Y%m%d_%H%M')
        
        list_RR_msg.append(path_out+f'{time_string}_RR_msg.png')
        list_RR_orig.append(path_out+f'{time_string}_RR_orig.png')
        list_msg.append(path_out+f'{time_string}_MSG_108_expats.png')

    # create gifs from the lists
    gif_maker(path_out, f'{time_string}_RR_msg', path_out, 250, 'RR_msg', '_RR_msg.png', list_RR_msg)
    print('gif for RR_msg created')
    
    gif_maker(path_out, f'{time_string}_RR_orig', path_out, 250, 'RR_orig', '_RR_orig.png', list_RR_orig)
    print('gif for RR_orig created')

    gif_maker(path_out, f'{time_string}_MSG_108', path_out, 250, 'MSG_108', '_MSG_108_expats.png', list_msg)
    print('gif for MSG_108 created')
    

    

    convert_gif_to_mp4(path_out, path_out)


def gif_maker(image_folder, gif_name, gif_path, gif_duration, channel, png_string, list_png_paths=None):
    """
    script to create animated gif from a folder containing images

    Args:
        image_folder (string): folder containing png images
        gif_name (string): string as filename for gif
        gif_path (string): path for gif file
        gif_duration (int): duration for gif (typical 250)
        channel(string): variable string for the gif
        png_string (string): string to filter images, e.g. "_IR_108_expats.png" or "_OPERA_reflectivity.png"
        list_png_paths (list): optional list of specific png paths to include   

    """
    from PIL import Image
    import glob
    
    import matplotlib.pyplot as plt
    
    # read files into the fil elist
    image_array = []

    # if list_png_paths is provided read from that, otherwise use glob
    if list_png_paths is not None:
        for file in list_png_paths:
            image = Image.open(file)
            image_array.append(image)
    else:
        for file in sorted(glob.glob(image_folder+'*'+png_string)):
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

def convert_gif_to_mp4(gif_path, mp4_path):
    """
    Converts a GIF file to MP4 format using ffmpeg.
    
    Args:
        gif_path (str): Path to the input GIF file.
        mp4_path (str): Path where the output MP4 file will be saved.
    """
    import os
    import subprocess

    # Root directory containing GIFs
    input_root = gif_path

    # Directory where MP4s will be saved
    output_root = mp4_path
    os.makedirs(output_root, exist_ok=True)

    # Output video settings
    output_fps = 10

    for dirpath, _, filenames in os.walk(input_root):
        for filename in filenames:
            if filename.lower().endswith(".gif"):
                gif_path = os.path.join(dirpath, filename)

                # Get relative path from input_root and use it to build output path
                rel_path = os.path.relpath(dirpath, input_root)
                output_dir = os.path.join(output_root, rel_path)
                os.makedirs(output_dir, exist_ok=True)

                mp4_filename = os.path.splitext(filename)[0] + ".mp4"
                mp4_path = os.path.join(output_dir, mp4_filename)

                if os.path.exists(mp4_path):
                    print(f"🔁 Skipping existing MP4: {mp4_path}")
                    continue

                print(f"🎞️ Converting: {gif_path} -> {mp4_path}")

                try:
                    subprocess.run([
                        "ffmpeg",
                        "-y",  # overwrite output if needed
                        "-i", gif_path,
                        "-vf", f"fps={output_fps},pad=ceil(iw/2)*2:ceil(ih/2)*2",  # even dimensions
                        "-c:v", "libx264",
                        "-preset", "slow",         # compression trade-off: slower = smaller file
                        "-crf", "32",
                        "-pix_fmt", "yuv420p",
                        "-movflags", "faststart",
                        mp4_path
                    ], check=True)
                    print(f"✅ Saved MP4: {mp4_path}")
                except subprocess.CalledProcessError as e:
                    print(f"❌ FFmpeg failed for {gif_path}: {e}")
                    
                    
                    
if __name__ == "__main__":
    main()