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
=========================================================================
# $Date: 2023-01-17
# $LastChangedBy: graham $
# =========================================================================
# DESCRIPTION: Main algorithm for the AIS message parsing toolbox.
# =========================================================================
# PROJECT       : S1-CPAF
# COMPANY       : ESA
# COMPONENT     : MSG_PARSE
# LANGUAGE      : Python
# -------------------------------------------------------------------------
INPUT: 16 .txt files with detection results; one file for each stream 
    (i.e. outputs from the AIS_DEMOD toolbox)
OUTPUT: For each of the four channels the following output fiels are created 
-	.txt file with ais message detections and a txt file with invalid detections
-	Detection statistics for each input stream as a txt file
==========================================================================
Change log:
24th Jan 2025
	- Added indexing into AIFM timestamp file


20th Jan 2025
    - Updated: In merge_ais_channel: 
        * all_messages = [pd.read_csv(f, low_memory=False) for f in filelist] - low_memory = False
    	to account for columns with different dtypes. 
        * Added check to ensure that the dataframes are not empty before concatenation
        non_empty_dataframes = [df for df in all_messages if not df.empty]
    - Reversed Lat lon value in cleanse_AIS()
    - Added initialisation of: all_detections = []
    - Added additional checks to Parse_AIS() incase the dataframe is empty
    - Replaced exception handling with if statement to only process messages 0-27 (excluding type 24)

"""

import pandas as pd
import sys
import bitarray as bitarray
import os
from datetime import datetime, timedelta
import numpy as np


from pyais.messages import MSG_CLASS
from pyais.util import get_int

# For data cleansing
from shapely.geometry import Polygon, Point
import xml.etree.ElementTree as ET

def GPS_to_UTC(GPS_time):
    GPS_zero_str = "1980-01-06T00:00:00" # UTC time when GPS time was zero (reference)
    GPS_zero_ref =  datetime.strptime(GPS_zero_str,'%Y-%m-%dT%H:%M:%S') # Convert to time object
    leap_sec_str_list =['1981-06-30T23:59:59', # Announced leap seconds in string format
                        '1982-06-30T23:59:59', # The time values are BEFORE the addition of the leap seconds!
                        '1983-06-30T23:59:59',
                        '1985-06-30T23:59:59',
                        '1987-12-31T23:59:59',
                        '1989-12-31T23:59:59',
                        '1990-12-31T23:59:59',
                        '1992-06-30T23:59:59',
                        '1993-06-30T23:59:59',
                        '1994-06-30T23:59:59',
                        '1995-12-31T23:59:59',
                        '1997-06-30T23:59:59',
                        '1998-12-31T23:59:59',
                        '2005-12-31T23:59:59',
                        '2008-12-31T23:59:59',
                        '2012-06-30T23:59:59',
                        '2015-06-30T23:59:59',
                        '2017-12-31T23:59:59']
    # leap_sec_time_list = []  # List with the announced leap seconds in time object format
    num_leap = 0 # Counter containing the added leap seconds for the introduced GPS time
    for ann_leap in leap_sec_str_list: # Convert leap seconds to time object
        # leap_sec_time_list.append(datetime.strptime(ann_leap,'%Y-%m-%dT%H:%M:%S'))  # Convert to time object
        leap_sec_date = datetime.strptime(ann_leap,'%Y-%m-%dT%H:%M:%S')  # Convert to time object
        if (GPS_zero_ref + timedelta(seconds=GPS_time) - leap_sec_date).total_seconds() - num_leap > 0: # GPS is ahead of UTC since 30-06-1981
            num_leap+=1

    UTC_time = GPS_zero_ref + timedelta(seconds = GPS_time) - timedelta(seconds = num_leap)

    return UTC_time.strftime("%Y%m%dT%H%M%S"), UTC_time.microsecond / 1000.0



def decode_binary_ais_messages(InputFile0, InputFile1, InputFile2, InputFile3, InputTimestamps, OutFilename, OutDir):
    #==========================================================================
    # DESCRIPTION:
    # Funciton for Parsing the binary AIS messages from the four streams of each channel
    # INPUTS
    # - (Input0, Input1, Input2, Input3): Detection text files (H, V, HmV, HV) 
    #   for one channel as output from DEMOD_AIS
    # - Out: Filename for outputs
    # - InputTimestamps: .npy file with the AIFM timestamps as extracted from the raw AIS data.
    # - OutDir: Directory for file outputs
    #
    # OUTPUTS
    # -	Pandas data frame with detections for the channel 
    # -	Detection statistics for each input stream
    #==========================================================================

    SaveOutput = os.path.join(OutDir, OutFilename)
    
    # Check if the directory exists if not create it
    if not os.path.exists(OutDir):
      os.makedirs(OutDir)

    # Load files as CSV
    detections_H = pd.read_csv(InputFile0, names = ['Frequency','time','bits'] )
    detections_V = pd.read_csv(InputFile1, names = ['Frequency', 'time','bits'] )
    detections_HmV = pd.read_csv(InputFile2, names = ['Frequency','time','bits'] )
    detections_HV = pd.read_csv(InputFile3, names = ['Frequency','time','bits'] )
    
    # Drop duplcated detections within file
    detections_H = detections_H.drop_duplicates(subset=['bits'])
    detections_V = detections_V.drop_duplicates(subset=['bits'])
    detections_HmV = detections_HmV.drop_duplicates(subset=['bits'])
    detections_HV = detections_HV.drop_duplicates(subset=['bits'])

    # Combine detections from all polarizations
    #detections = pd.concat([detections_H, detections_V, detections_HmV, detections_HV], axis=0)
    
    # Filter out empty DataFrames
    dataframes = [detections_H, detections_V, detections_HmV, detections_HV]
    non_empty_dataframes = [df for df in dataframes if not df.empty]

    # Concatenate only non-empty DataFrames
    if non_empty_dataframes:
        detections = pd.concat(non_empty_dataframes, axis=0)
    else: # Return if all channels are empty
        return
    
    detections = detections.drop_duplicates(subset=['bits'])
    detections = detections.reset_index(drop = True)
    bits = detections['bits']
    msg_list = []

    # Add AIFM timestamp to detected data.
    timestamps = np.loadtxt(InputTimestamps)
    sample_timestamps = np.repeat(timestamps, 168)
    #detections['time'] = sample_timestamps[detections['time']]
    message_timestamps = sample_timestamps[detections['time'].astype(int)]  
    
    
    
    # Parse each AIS message
    invalid_msgs = []
    for i in range(0,len(bits)):
        bit_str = bits[i]
        bits1 = bitarray.bitarray(bit_str)
        msg_type = get_int(bits1,0,6)
        #try:
        
        # Only decodes messages tyep 0-27 excluding type 24 which is not supported.
        if msg_type < 28 and msg_type != 24:
            # Identify message type and decode    
            msg_content = MSG_CLASS[msg_type].from_bitarray(bits1).asdict()

            # Convert GPS timestamp to UTC string and miliseconds
            #UTC_time, mili_seconds = GPS_to_UTC(detections['time'][i])
            UTC_time, mili_seconds = GPS_to_UTC(message_timestamps[i])
            
            # Update dictionary with auxillary values.
            msg_content.update({"time_utc": UTC_time})
            msg_content.update({"mili_seconds": mili_seconds})
            msg_content.update({"doppler_freq": detections['Frequency'][i]})

            # Append to list
            msg_list.append(msg_content.copy())
        else:   
        #except Exception as e:
            invalid_msgs.append(detections.loc[i, ['Frequency','time','bits']])



    # Define Collumn names for output csv file
    column_names = ['mmsi', 'time_utc', 'mili_seconds', 'doppler_freq','msg_type', 'repeat', 'status', 'turn','speed', 'accuracy', 
            'lat', 'lon',  'course', 'heading',
           'second', 'maneuver', 'raim', 'radio', 'gnss',  'assigned', 'draught', 'destination',
           'dte', 'epfd', 'reserved_1', 'reserved_2', 'ship_type', 'shipname',
           'to_bow', 'to_port', 'to_starboard', 'to_stern', 'mmsi1', 'mmsi2',
           'mmsi3', 'mmsi4', 'mmsiseq1', 'mmsiseq2', 'mmsiseq3', 'mmsiseq4',
           'increment1', 'increment2', 'offset1', 'offset2', 'dac', 'data',
           'dest_mmsi', 'fid', 'retransmit', 'seqno', 'day', 'hour', 'minute',
           'month', 'year', 'band', 'cs', 'display', 'dsc', 'msg22', 'interval',
           'ne_lat', 'ne_lon', 'quiet','spare_1', 'spare_2', 'spare_3', 'spare_4', 'station_type',
           'sw_lat', 'sw_lon', 'txrx', 'text', 'callsign', 'model', 'partno',
           'serial', 'vendorid', 'addressed', 'structured', 'ais_version', 'imo', 'alt', 'type1_1',
           'offset1_1', 'type1_2', 'offset1_2', 'type2_1', 'offset2_1', 'number1', 'timeout1',
           'number2', 'timeout2', 'offset3', 'number3', 'timeout3', 'increment3', 'offset4', 'number4', 
           'timeout4', 'increment4', 'aid_type', 'name', 'off_position', 'virtual_aid', 'name_ext', 'channel_a', 
           'channel_b', 'power', 'dest1', 'empty_1', 'dest2', 'empty_2', 'band_a', 'band_b', 'zonesize', 'app_id']
    
    # Convert list to Pandas datafram
    df = pd.DataFrame(msg_list, columns = column_names)
    df = df.sort_values(by='time_utc')  

    df.to_csv(SaveOutput + '_detections.txt', index=False)
    
    invalid_df = pd.DataFrame(invalid_msgs)
    invalid_df.to_csv(SaveOutput + '_invalid.txt', index=False)
    
    
    # Drop duplicated MMSIs
    unique_MMSIs = df.drop_duplicates(subset=['mmsi'])


    # detection report
    with open(SaveOutput + '_summary.txt', "w") as file:
        # Redirect output to the file
        sys.stdout = file

        # Detection Statistics
        print("Detections statistics for: " + str(OutFilename))
        print("Detected messages in H pol: " + str(len(detections_H)))
        print("Total messages in V pol: " + str(len(detections_V)))
        print("Total messages in H + V pol: " + str(len(detections_HV)))
        print("Total messages in H - V pol: " + str(len(detections_HmV)))
        # Total Detection
        print("Total number of unique messages detected: " + str(len(detections)))
        # Unique detections
        print("Total number of unique ship IDs (MMSI): " + str(len(unique_MMSIs)))
        print("Total number of invald message types discarded: " + str(len(invalid_msgs)))
        print("=================================================================")
        sys.stdout = sys.__stdout__



def merge_ais_channels(filelist, output_filename, output_directory):
    #==========================================================================
    # DESCRIPTION:
    # Merges the output txt files to a single file and generates the naming convention.
    # INPUTS
    # - filelist: List of files with the detected AIS messages in each channel to merge
    # - OutFile: Output filename for merged detections
    # - OutDir: Output directory to save file in
    #
    # OUTPUT
    # -- .txt file with containing all detected messages referrenceing the start 
    #    and stop time of aqusitions
    #    The format 'S1C_AI_L1_YYYYMMDDTHHMM_YYYYMMDDTHHMM.txt'
    #==========================================================================
        
    all_messages = [pd.read_csv(f, low_memory=False) for f in filelist]
    
    # Concatenate all non-empty DataFrames into a single DataFrame
    non_empty_dataframes = [df for df in all_messages if not df.empty]
    all_messages = pd.concat(non_empty_dataframes, ignore_index=True)
    all_messages = all_messages.sort_values(by='time_utc')
    
    all_messages.to_csv(os.path.join(output_directory, output_filename), index=False)
    

def generate_file_name(InputTimestamps):
    
    timestamps = np.loadtxt(InputTimestamps)
    start_time ,_ = GPS_to_UTC(timestamps[0])
    end_time ,_ = GPS_to_UTC(timestamps[-1])
    filename = f"S1C_AI_L1_{start_time}_{end_time}.txt"
    return filename



def add_quality_flag(filename_ais, annot_file, output_file):
    
    data_ais = pd.read_csv(filename_ais)
    tree = ET.parse(annot_file)
    root = tree.getroot()
    
    ns = {'gml': 'http://www.opengis.net/gml'}
    coordinates_text = root.find('.//gml:coordinates', ns).text

    # Extract and convert coordinates
    coordinates = [tuple(map(float, coord.split(','))) for coord in coordinates_text.split()]
    lats, lons = zip(*coordinates)


    # Create a Polygon object of the Field of View (FoV)
    polygon = Polygon(coordinates)


    # Function to check if the AIS message is within the FoV
    def is_within_polygon(row):
        point = Point(row['lat'], row['lon'])
        return polygon.contains(point)


    # Update the data frame with a quality indicator:
    # 1= valid AIS data 
    # 2 = invalid AIS message)
    data_ais.insert(4, 'quality', data_ais.apply(lambda row: 1 if is_within_polygon(row) else 2, axis=1))
    
    # Save dataframe to file:
    data_ais.to_csv(output_file, index=False)
    
    # Print numbuer of AIS messages detected
    total_count = data_ais.shape[0]
    print(f"Total number of AIS messages detected: {total_count}")



