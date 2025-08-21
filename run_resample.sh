#!/bin/bash
# command to launch nohup ./run_resample.sh > bash_script_out.txt
source .venv/bin/activate
python3 process/resample_rr_to_msg.py > radar_de_process_2023_2022.txt
