"""
Code to access data on the Arcodatahub S3 bucket, 
which is used to store the data of the radar italian network.

activate virtual environment .arco_env before running this code with
source /home/claudia/.arco_env/bin/activate

call the script with 
python3 -m process.access_arcodatahub_data


"""
import xarray as xr
from readers.radars_ita import read_arcohub_data
import pdb
from datetime import datetime
from readers.read_msg import prepare_msg_data
import s3fs
from process.utils import generate_regular_grid, regrid_data, is_valid_date
from readers.data_buckets_funcs import upload_to_bucket, check_file_bucket, init_s3, download_file

def main():

    # reading a random existing file from msg to get resolution
    ds_msg = prepare_msg_data()

    # list all years from 2024 to 2014
    yy_list = [str(year) for year in range(2020, 2018, -1)]
    # list all months in a year
    mm_list = [f'{month:02d}' for month in range(4, 10)]
    # list all days in a month
    dd_list = [f'{day:02d}' for day in range(1, 32)]
    print(yy_list, mm_list, dd_list)
    
    # loop on all years, months, days to find yy mm dd for processing
    for yy in yy_list:
        print('processin year ', yy)
        for mm in mm_list:
            print('processin month ', mm)
            for dd in dd_list:
                print('processin day ', dd)
                

                # print the date to process
                print("date", yy, mm, dd)
                
                # process only if the date is valid
                if is_valid_date(yy, mm, dd):
                    
                    # read file list from the bucket
                    s3 = init_s3()
                    file_bucket_list = list_files_bucket(s3, "expats-radar-it-PC")
                    file_to_find = yy+mm+dd+'_RR_15min_msg_res_ita.nc.gz'
                    
                    # check if filename is in file_bucket_list, if yes set file_exist to True, otherwise False
                    if file_to_find in [obj['Key'] for obj in file_bucket_list]:
                        print('file already exists on bucket expats-radar-italy civil protection - skipping date:', yy, mm, dd)
                    else:
                        print('generating ncdf file for ', yy+mm+dd)
                    
                    path_tar = "sat_data/temp_radar_ita/"
                    dest_path_ncdf = 'sat_data/temp_radar_ita/'+yy+'/'
                    
                    # check if there folder exists, otherwise unzip tar.gz file
                    if not os.path.exists(dest_path_ncdf):
                        
                        # creating tar folder for the year if it does not exist
                        os.makedirs(dest_path_ncdf)
                    
                    # read dataset from Arcodatahub S3 bucket for the day cropped on domain
                    ds_radar = read_arcohub_data(yy,mm,dd)
                    pdb.set_trace()

                    # resample radar data to 15 min 
                    ds_radar_msg = ds_radar.resample(time='15T').sum(skipna=True)
    

                    # check if file was uploaded
                    if check:
                        print(f'{ncdf_radar_gz} file uploaded to bucket')
                        
                        # remove ncdf file produced
                        # cut .gz at the end of ncdf_radar
                        ncdf_radar_name = ncdf_radar_gz[:-3]
                        
                        # if file ncdf_radar+name exists, remove it
                        if os.path.exists(path_out+ncdf_radar_name):
                            os.remove(path_out+ncdf_radar_name)

                    else:
                        print('file not uploaded to bucket - abort')
                        break

    # remove the tar.gz file to save space
    if os.path.exists(path_tar+'/'+tar_filename):
        print('removing tar.gz file to save space')
        os.remove(path_tar+'/'+tar_filename)
    


   
    


    # Save the cropped dataset to a local file -zarr
    #ds_year.to_zarr(f"/sat_data/temp_radar_ita/{year}_italian_radar.zarr", mode='w')



if __name__ == "__main__":
    
    main()
