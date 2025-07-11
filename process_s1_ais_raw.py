#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Created on Thu Jan  4 11:08:52 2024

@author: stefan.graham

Full processing chain to process the Sentinel-1 C&D AIS Instrument Source Packet (ISP) to AIS messages.

Input: SAFE folder containing the AIS ISPs (.dat/.bin)
Outputs include:
- Extracted headers and raw data (.wav)
- .txt file with all detected AIS messages with the naming convention _S1C_AI_L1_YYYYMMDDTHHMM_YYYYMMDDTHHMM.txt
-	Detection statistics for each input stream as a .txt file
-	Text file with potentialially invalid messages.
#############################################
Change log:
    30th Jun 2025
    - Updated file handling to be cross-platform compatible
    - Updated toolbox to be executed from the command line
    20th Jan 2025
    - Added check to ensure that *_detections.txt exist before merge_ais_channels()
"""

import os
import time
import subprocess
from multiprocessing import Pool
import argparse
import glob
import platform

# Add paths to Python's sys.path
import sys
sys.path.append(os.path.join(os.getcwd(), 'src'))


# Import custom packages from src
from parse_ais import decode_binary_ais_messages, merge_ais_channels, generate_file_name, add_quality_flag
from isp_extraction import extract_isp2asc

    
    
    
def process_channel(args):
    
    # process_channel calls the ESA AIS demodulation C-executable with input .wav files
    k, input_dir, output_dir, wav_name = args
    filename_in = os.path.join(input_dir, f"{wav_name}.asc")
    filename_out = os.path.join(output_dir, f"{wav_name}.txt")
    
    if platform.system() == "Windows":
        esa_ais_executable = "AIS_receiver_win.exe"
    elif platform.system() == "Linux":
        esa_ais_executable = './AIS_receiver_linux'
        
    elif platform.system() == "Darwin":
        esa_ais_executable = './AIS_receiver_mac'
    else:
        raise RuntimeError("Unsupported OS")
        
    # Call executable from command line
    # The files from k = 8 to 16 are the sat-AIS channels with message lengths 96 bits.
    #command = ['./AIS_receiver', '168' if k < 8 else '96', filename_out, filename_in]
    command = [esa_ais_executable, '168' if k < 8 else '96', filename_out, filename_in]
    subprocess.run(command)
    

def find_input_dat_file(input_dir):
    # Get all .dat files in the directory
    all_dat_files = glob.glob(os.path.join(input_dir, '*.dat'))

    # Filter out files ending with -annot.dat
    main_dat_files = [f for f in all_dat_files if not f.endswith('-annot.dat')]

    if not main_dat_files:
        raise FileNotFoundError("No main .dat file found in the directory.")

    if len(main_dat_files) > 1:
        raise ValueError(f"Multiple .dat files found (excluding -annot): {main_dat_files}")

    return main_dat_files[0]

def demodulate_all_ais_channels(sInputDir, sOutputDir):
    

    # Find Instrument Source Packet file in SAFE folder
    filename = find_input_dat_file(sInputDir)

    # Input manifest.safe file containing the AIS/SAR Field of view polygon.
    mainfest_l0 =   os.path.join(sInputDir, 'manifest.safe')

    # Extract the Instrument Source Packets (ISPs) into .wav files for further processing
    extract_isp2asc(sInputDir, filename, sOutputDir, isp_idx = [0,10E10])
    
    # ESA AIS Demodulation
    fieldnames = ["ais_ch0_pol_H", "ais_ch0_pol_V", "ais_ch0_pol_HV", "ais_ch0_pol_HmV",
                  "ais_ch1_pol_H", "ais_ch1_pol_V", "ais_ch1_pol_HV", "ais_ch1_pol_HmV",
                  "ais_ch2_pol_H", "ais_ch2_pol_V", "ais_ch2_pol_HV", "ais_ch2_pol_HmV",
                  "ais_ch3_pol_H", "ais_ch3_pol_V", "ais_ch3_pol_HV", "ais_ch3_pol_HmV"]
    
    
    # Define the directory where to load the .asc files to the ESA demodulator.
    # In this case it is the output directory from previous step where the .wav files have been saved.
    sInputDir = sOutputDir
    
    # Paralle processing of AIS channels including the linear polarization combinations
    start_time = time.time()
    with Pool() as pool:
        pool.map(process_channel, [(k, sInputDir, sOutputDir, fieldnames[k]) for k in range(16)])

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # Save run time to file    
    with open(os.path.join(sOutputDir, 'timing_results.txt'), 'w') as file:
        file.write(f'File name: {filename}\n')
        file.write(f'Total time for parallel demodulation: {elapsed_time:.2f} seconds\n')

    
    # Parse AIS - Convert the binary AIS messages to Lat long etc.
    num_channels = 4
    input_timestamps = os.path.join(sInputDir, 'AIFM_timestamps.txt')
    for channel in range(0, num_channels):
        input0 = os.path.join(sInputDir, f'ais_ch{channel}_pol_H.txt')
        input1 = os.path.join(sInputDir, f'ais_ch{channel}_pol_V.txt')
        input2 = os.path.join(sInputDir, f'ais_ch{channel}_pol_HmV.txt')
        input3 = os.path.join(sInputDir, f'ais_ch{channel}_pol_HV.txt')
        output_filename = f'ais_ch{channel}'
        # Parse messages for a single channel over all polarization combinations.
        decode_binary_ais_messages(input0, input1, input2, input3, input_timestamps, output_filename, sOutputDir)
        #########
    
    
    # Generate L1 filename based on the datatake start and end time:
    # S1C_AI_L1_YYYYMMDDTHHMM_YYYYMMDDTHHMM.txt
    l1_filename = generate_file_name(input_timestamps)    
    
    # Merge all channels into a single file
    #filelist = [sInputDir, 'ais_ch0_detections.txt', sInputDir+ '/ais_ch1_detections.txt', sInputDir+ '/ais_ch2_detections.txt', sInputDir+ '/ais_ch3_detections.txt']
    
    filelist = [os.path.join(sInputDir, f'ais_ch{ch}_detections.txt') for ch in range(num_channels)]
    
    existing_files = [file for file in filelist if os.path.exists(file)]
    merge_ais_channels(existing_files, l1_filename, sOutputDir)


    # Add quality flag using the Sentinel-1 AIS/SAR Field of View polygon
    # Only performed if the L0 manifest file is loaded
    try:
        # Output file overides the previous output from 'merge_ais_channels'
        file_out =  os.path.join(sInputDir, l1_filename)
        add_quality_flag(os.path.join(sInputDir, l1_filename), mainfest_l0, file_out)
    except NameError:
        print("Quality Flag NOT added.")

    
    

if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(description="Demodulate and decode Sentinel-1 RAW AIS files")
    parser.add_argument("input_dir", type=str, help="Path to the AIS RAW SAFE folder")
    parser.add_argument("-o", "--output_dir", type=str, default="ais_output", help="Directory to save output (default: ./ais_output)")
    args = parser.parse_args()
    
    input_dir = os.path.abspath(args.input_dir)
    # Default output directory if not provided
    output_dir = args.output_dir or os.path.join(os.getcwd(), "ais_output")
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    demodulate_all_ais_channels(input_dir, output_dir)

