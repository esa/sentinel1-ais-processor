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
Created on Tue Mar  5 11:45:03 2024

@author: stefan.graham

Functions to extract raw samples from the Sentinel-1C AIS Instrument Soruce Packets(ISPs)
AIS ISP structure reference:
    S1CD-RD-TAI-PM18-0123
    

Change log:
    20th Jan 2025
    - In extract_isp_data(): Chaned isp_idx = [0, int(1E10)] to start at 0 instead of 1:
    - In read_bin_ais() updated concatenation to end at  isp_cont instead of isp_cont-1
        e.g. user_data_bin = user_data_bin[:isp_cont]
    - In parse_ais_user_data(): Updated the loop to use append instead of indexing into empty array.


"""

from scipy.io import wavfile
import numpy as np
import os
import struct
import time



def s1_parse_binary_data(binary_data, fields_descriptor):
    """
    SS = s1_parse_binary_data(BINARY_DATA,FIELDS_DESCRIPTOR)
    Parses BINARY_DATA into fields described into FIELDS_DESCRIPTOR. Binary
    data is interpreted as unsigned int in Big-endian convention, and it is
    converted to its closer standard numerical representation.
    
    Inputs:
        BINARY_DATA: matrix of class uint8. Each column contains the sequence
        of bytes that is parsed into fields.
        FIELDS_DESCRIPTOR: cell with the name of each field in its first column
        and the corresponding number of bits in its second column.
    
    Outputs:
        SS: struct with parsed fields. Each field contains as many columns as
        BINARY_DATA.
    """
   
    if not isinstance(binary_data, np.ndarray) or binary_data.dtype != np.uint8:
        raise ValueError('Binary data shall be a uint8 array, was {} instead'.format(type(binary_data)))

    bits_v = np.array(fields_descriptor)[:, 1].astype(int)
    if any(bits_v <= 0):
        raise ValueError('Number of bits in fields descriptor shall be positive')
    elif any(bits_v > 64):
        raise ValueError('Number of bits in fields descriptor shall not exceed 64')

    if sum(bits_v) > 8 * binary_data.shape[0]:
        raise ValueError('Total number of bits in fields descriptor exceeds binary data')

    bit_count = 0
    ss = {}
    for k in range(len(fields_descriptor)):
        field = fields_descriptor[k][0]
        bits = bits_v[k]
        bit_indexes = np.arange(bit_count+1, bit_count + bits+1)
        bytes_span = np.unique(np.ceil(bit_indexes / 8).astype(int))


        if bits <= 8:
            cast_fn = np.uint8
        elif 9 <= bits <= 16:
            cast_fn = np.uint16
        elif 17 <= bits <= 32:
            cast_fn = np.uint32
        elif 33 <= bits <= 64:
            cast_fn = np.uint64


        ss[field] = cast_fn(0)
        for bb in bytes_span:
            idx = np.mod(bit_indexes[(bit_indexes >= (bb-1) * 8+1) & (bit_indexes <= bb * 8)]-1, 8)
            bit_shift = (max(idx) - min(idx) - 7)
            if bit_shift < 0:
                ss[field] += cast_fn(np.right_shift(np.left_shift(binary_data[bb-1], min(idx)), -bit_shift)) * 2 ** (8 * (bytes_span[-1] - bb))
            else:        
                ss[field] += cast_fn(np.left_shift(np.left_shift(binary_data[bb-1], min(idx)), bit_shift)) * 2 ** (8 * (bytes_span[-1] - bb))

        bit_count += bits

    return ss

def parse_bin_ais_headers(head_ccsds_bin, head_time_stamp_bin, head_pri_bin, head_sec_bin):
    # Parses binary stream from primary and secondary headers into structs
    
    head_ccsds_fields = [
        ('version_nr', 3),
        ('type', 1),
        ('psh_flag', 1),
        ('apid', 11),
        ('seq_flag', 2),
        ('seq_counter', 14),
        ('packet_len', 16)
    ]
    
    head_time_stamp_fields = [
        ('spare', 3),
        ('sc_pps_time', 45)
    ]
    
    
    head_pri_fields = [
        ('version_nr', 3),
        ('type', 1),
        ('psh_flag', 1),
        ('apid', 11),
        ('seq_flag', 2),
        ('seq_counter', 14),
        ('packet_len', 16)
    ]
    
    head_sec_fields = [
        ('ch', 4),
        ('type', 1),
        ('pol', 1),
        ('spare', 2),
        ('aisr_time', 32),
        ('reserved', 8)
    ]

    # Parse binary data into struct
    head_ccsds = s1_parse_binary_data(head_ccsds_bin, head_ccsds_fields)
    head_time_stamp = s1_parse_binary_data(head_time_stamp_bin, head_time_stamp_fields)
    head_pri = s1_parse_binary_data(head_pri_bin, head_pri_fields)
    head_sec = s1_parse_binary_data(head_sec_bin, head_sec_fields)

    return head_ccsds, head_time_stamp, head_pri, head_sec


def read_user_data(fid, data_size):
    """
    Reads user data from a file and converts it into 12-bit values.
    
    Args:
    fid: File object opened for reading in binary mode.
    data_size: Number of bytes to read (multiple of 6).
    
    Returns:
    user_data: List of lists, where each inner list contains four 12-bit unsigned integers.
    """
    user_data = []
    data_chunk = fid.read(data_size)
    if len(data_chunk) % 3 != 0:
        raise ValueError("Invalid data chunk size. Expected multiple of 6 bytes.")

    for i in range(0, len(data_chunk), 3):
        
        # Extract 3-byte chunk
        data_bytes = data_chunk[i:i+3]
        data_int = int.from_bytes(data_bytes, byteorder='big', signed=False)
        binary_str = format(data_int, '024b')

        in_phase = int(binary_str[0:12],2)
        quadrature = int(binary_str[12::],2)
        
        # Convert to signed integers if necessary
        if binary_str[0] == '1':
            in_phase -= 2 ** 12
        if binary_str[12] == '1':
            quadrature -= 2 ** 12
        
        user_data.append(in_phase)
        user_data.append(quadrature)
        
    return user_data

def read_bin_ais(path_bin, isp_idx):
    # extracts binary data and headers

    # Headers byte sizes (fixed)
    head_ccsds_bytes = 6
    head_time_stamp_bytes = 6
    head_pri_bytes = 6
    head_sec_bytes = 6
    head_bytes = 24  # sum of the above

    # First and last ISPs to read
    isp_start = isp_idx[0]
    isp_end = isp_idx[1]

    # Initialize output variables (they are cropped later)
    headers_bin = np.zeros((head_bytes, int(1E6)), dtype=np.uint8)
    user_data_bin = [None] * int(1E6)
    
    

    with open(path_bin, 'rb') as fid:
        filesize = fid.seek(0, 2)
        fid.seek(0, 0)

        n = 0
        isp_cont = 0
        rem_bytes = filesize

        # Loop starts here to create data matrix
        while rem_bytes > head_bytes and n < isp_end:
        
            # Read header bytes
            heads_bin_aux = np.frombuffer(fid.read(head_bytes), dtype=np.uint8)
        
            # Get length of ISP and read data bytes
            isp_length = struct.unpack('<H', heads_bin_aux[[5,4]])[0]
            user_data_bytes = isp_length - (head_bytes - head_ccsds_bytes) + 1
            
            user_data_aux = read_user_data(fid, int(user_data_bytes))

            user_data_bin[n - isp_start] = user_data_aux
            headers_bin[:, n - isp_start] = heads_bin_aux

            isp_cont += 1
            n += 1
            # Remaining bytes on file
            rem_bytes = filesize - fid.tell()

        # Loop ends here

    # Separate headers. Crop to the actual number of read ISPs
    idx_ccsds = slice(0, head_ccsds_bytes)
    idx_time_stamp = slice(head_ccsds_bytes, head_ccsds_bytes + head_time_stamp_bytes)
    idx_pri = slice(head_ccsds_bytes + head_time_stamp_bytes, head_ccsds_bytes + head_time_stamp_bytes + head_pri_bytes)
    idx_sec = slice(head_ccsds_bytes + head_time_stamp_bytes + head_pri_bytes, head_bytes)


    user_data_bin = user_data_bin[:isp_cont]
    head_ccsds_bin = headers_bin[idx_ccsds, :isp_cont]
    head_time_stamp_bin = headers_bin[idx_time_stamp, :isp_cont]
    head_pri_bin = headers_bin[idx_pri, :isp_cont]
    head_sec_bin = headers_bin[idx_sec, :isp_cont]
    user_data_bin = user_data_bin[:isp_cont]

    return head_ccsds_bin, head_time_stamp_bin, head_pri_bin, head_sec_bin, user_data_bin




def extract_isp_data(sInputDir, filename, sOutputDir, isp_idx = None):
    """
    Main function to extract the AIS Instrument Source Packets (ISPs) into .wav files. 
    Outputs: 
        - Eight .wav files. One for each channel and polarization including the linear pol (H+V) combinations.
        - ISP headers and timestamp file as .npy files. 
    
    Inputs:
        - sInputDir: Input directory containing .dat file.
        - filename: Filename for binary ISP (.dat or .bin)
        - sOutputDir: Output directory to save the files in.
    """
    
    if isp_idx is None:
        isp_idx = [0, int(1E10)]
    elif len(isp_idx) != 2:
        raise ValueError('ISP_IDX length must be 2')
        
    head_ccsds, head_time_stamp, head_pri, head_sec, ais_ch0_pol_H, ais_ch0_pol_V, ais_ch1_pol_H, ais_ch1_pol_V, ais_ch2_pol_H, ais_ch2_pol_V, \
    ais_ch3_pol_H, ais_ch3_pol_V = s1_read_ais_bin(os.path.join(sInputDir, filename),isp_idx)

    # Generate linear polarization combinations for improved detection
    ais_ch0_pol_HV = ais_ch0_pol_H + ais_ch0_pol_V
    ais_ch0_pol_HmV = ais_ch0_pol_H - ais_ch0_pol_V
    ais_ch1_pol_HV = ais_ch1_pol_H + ais_ch1_pol_V
    ais_ch1_pol_HmV = ais_ch1_pol_H - ais_ch1_pol_V
    ais_ch2_pol_HV = ais_ch2_pol_H + ais_ch2_pol_V
    ais_ch2_pol_HmV = ais_ch2_pol_H - ais_ch2_pol_V
    ais_ch3_pol_HV = ais_ch3_pol_H + ais_ch3_pol_V
    ais_ch3_pol_HmV = ais_ch3_pol_H - ais_ch3_pol_V

    fieldnames = ["ais_ch0_pol_H", "ais_ch0_pol_V", "ais_ch0_pol_HV", "ais_ch0_pol_HmV",
                  "ais_ch1_pol_H", "ais_ch1_pol_V", "ais_ch1_pol_HV", "ais_ch1_pol_HmV",
                  "ais_ch2_pol_H", "ais_ch2_pol_V", "ais_ch2_pol_HV", "ais_ch2_pol_HmV",
                  "ais_ch3_pol_H", "ais_ch3_pol_V", "ais_ch3_pol_HV", "ais_ch3_pol_HmV"]

    # Save each output as Wav file
    sampling_freq = 3 * 9600
    for k in range(len(fieldnames)):
        wav_name = fieldnames[k]
        IQ_samples_to_save = eval(wav_name)
        wav_name = os.path.join(sOutputDir, wav_name + '.wav')
        wavfile.write(wav_name, sampling_freq, np.column_stack((np.int16(np.real(IQ_samples_to_save)), np.int16(np.imag(IQ_samples_to_save)))))

    # Save header data and timestamp separately (as .npy files)
    np.save(os.path.join(sOutputDir, 'headers.npy'), {'head_ccsds': head_ccsds, 'head_time_stamp': head_time_stamp, 'head_pri': head_pri, 'head_sec': head_sec})
    #np.save(os.path.join(sOutputDir, 'timestamp.npy'), {'head_time_stamp': head_time_stamp})
    np.savetxt(os.path.join(sOutputDir, 'AIFM_timestamps.txt'), head_time_stamp)

def extract_isp2asc(sInputDir, filename, sOutputDir):
    """
    Main function to extract the AIS Instrument Source Packets (ISPs).
    Outputs: 
        - 16 .asc files. One for each channel and polarization including the linear pol (H +/- V) combinations.
        - ISP headers and timestamp file as .npy files. 
    
    Inputs:
        - sInputDir: Input directory containing .dat file.
        - filename: Filename for binary ISP (.dat or .bin)
        - sOutputDir: Output directory to save the files in.
    """
    head_ccsds, head_time_stamp, head_pri, head_sec, ais_ch0_pol_H, ais_ch0_pol_V, ais_ch1_pol_H, ais_ch1_pol_V, ais_ch2_pol_H, ais_ch2_pol_V, \
    ais_ch3_pol_H, ais_ch3_pol_V = s1_read_ais_bin(os.path.join(sInputDir, filename), [1, 48498])
    
    # Generate linear polarization combinations for improved detection
    ais_ch0_pol_HV = ais_ch0_pol_H + ais_ch0_pol_V
    ais_ch0_pol_HmV = ais_ch0_pol_H - ais_ch0_pol_V
    ais_ch1_pol_HV = ais_ch1_pol_H + ais_ch1_pol_V
    ais_ch1_pol_HmV = ais_ch1_pol_H - ais_ch1_pol_V
    ais_ch2_pol_HV = ais_ch2_pol_H + ais_ch2_pol_V
    ais_ch2_pol_HmV = ais_ch2_pol_H - ais_ch2_pol_V
    ais_ch3_pol_HV = ais_ch3_pol_H + ais_ch3_pol_V
    ais_ch3_pol_HmV = ais_ch3_pol_H - ais_ch3_pol_V

    fieldnames = ["ais_ch0_pol_H", "ais_ch0_pol_V", "ais_ch0_pol_HV", "ais_ch0_pol_HmV",
                  "ais_ch1_pol_H", "ais_ch1_pol_V", "ais_ch1_pol_HV", "ais_ch1_pol_HmV",
                  "ais_ch2_pol_H", "ais_ch2_pol_V", "ais_ch2_pol_HV", "ais_ch2_pol_HmV",
                  "ais_ch3_pol_H", "ais_ch3_pol_V", "ais_ch3_pol_HV", "ais_ch3_pol_HmV"]

    
    # Save each output as Wav file
    sampling_freq = 3 * 9600
    for k in range(len(fieldnames)):
        wav_name = fieldnames[k]
        IQ_samples_to_save = eval(wav_name)
        asc_name = os.path.join(sOutputDir, wav_name + '.asc')
        real_imag_array = np.column_stack((np.int16(np.real(IQ_samples_to_save)), np.int16(np.imag(IQ_samples_to_save))))
        # Save to txt file with comma separation
        np.savetxt(asc_name, real_imag_array, fmt="%d", delimiter=",")
        

    # Save header data and timestamp separately (as .npy files)
    np.save(os.path.join(sOutputDir, 'headers.npy'), {'head_ccsds': head_ccsds, 'head_time_stamp': head_time_stamp, 'head_pri': head_pri, 'head_sec': head_sec})
    #np.save(os.path.join(sOutputDir, 'timestamp.npy'), {'head_time_stamp': head_time_stamp})
    np.savetxt(os.path.join(sOutputDir, 'AIFM_timestamps.txt'), head_time_stamp)
    

def parse_ais_user_data(head_sec, user_data_bin):
    # Reads AISR secondary header and user_data
    # OUTPUT:  complex IQ streams for each AIS channel and polarization
    # NOTES: Assumes all four channels are active.
    # Author: Stefan Graham

    # Check that all four channels are active
    if any(head_sec['ch'] != 15):
        print('READING ERROR: Error in the number of Channels')
        return
    
    ais_ch0_pol_H_list = []
    ais_ch0_pol_V_list = []
    ais_ch1_pol_H_list = []
    ais_ch1_pol_V_list = []
    ais_ch2_pol_H_list = []
    ais_ch2_pol_V_list = []
    ais_ch3_pol_H_list = []
    ais_ch3_pol_V_list = []
    
      

    k = 0
    jj = 0
    nr_samples = 168

    while k < len(user_data_bin):
        user_data = np.array(user_data_bin[k])

        ais_ch0_pol_H_list.append(user_data[0:len(user_data):16] + 1j * user_data[1:len(user_data):16])
        ais_ch0_pol_V_list.append(user_data[2:len(user_data):16] + 1j * user_data[3:len(user_data):16])

        ais_ch1_pol_H_list.append(user_data[4:len(user_data):16] + 1j * user_data[5:len(user_data):16])
        ais_ch1_pol_V_list.append(user_data[6:len(user_data):16] + 1j * user_data[7:len(user_data):16])

        ais_ch2_pol_H_list.append(user_data[8:len(user_data):16] + 1j * user_data[9:len(user_data):16])
        ais_ch2_pol_V_list.append(user_data[10:len(user_data):16] + 1j * user_data[11:len(user_data):16])

        ais_ch3_pol_H_list.append(user_data[12:len(user_data):16] + 1j * user_data[13:len(user_data):16])
        ais_ch3_pol_V_list.append(user_data[14:len(user_data):16] + 1j * user_data[15:len(user_data):16])

        k += 1
        jj += nr_samples

    ais_ch0_pol_H = np.concatenate(ais_ch0_pol_H_list).astype(np.complex64)
    ais_ch0_pol_V = np.concatenate(ais_ch0_pol_V_list).astype(np.complex64)
    ais_ch1_pol_H = np.concatenate(ais_ch1_pol_H_list).astype(np.complex64)
    ais_ch1_pol_V = np.concatenate(ais_ch1_pol_V_list).astype(np.complex64)
    ais_ch2_pol_H = np.concatenate(ais_ch2_pol_H_list).astype(np.complex64)
    ais_ch2_pol_V = np.concatenate(ais_ch2_pol_V_list).astype(np.complex64)
    ais_ch3_pol_H = np.concatenate(ais_ch3_pol_H_list).astype(np.complex64)
    ais_ch3_pol_V = np.concatenate(ais_ch3_pol_V_list).astype(np.complex64)


    return ais_ch0_pol_H, ais_ch0_pol_V, ais_ch1_pol_H, ais_ch1_pol_V, ais_ch2_pol_H, ais_ch2_pol_V, ais_ch3_pol_H, ais_ch3_pol_V



def s1_read_ais_bin(path_bin, isp_idx=None):
    """
    Reads Sentinel-1 AIS ISP file extracting headers and IQ stream on each channel
    INPUTS:
        PATH_BIN: path to S1 .bin file
        ISP_IDX: vector with first and last ISP to read (default=[1,Inf])
    OUTPUT:
        HEAD_CCSDS: CCSDS header decoded into struct for each ISP
        HEAD_TIME_DEC: Spacecraft Pulse Pr Sec timestamp of each ISP in decimal
        HEAD_PRI: Primary header decoded into struct for each ISP
        HEAD_SEC: Secondary header decoded into struct for each ISP
        AIS_CHX_POL_Y: complex IQ samples in on Channel X and polarization Y.
    """
    
    if isp_idx is None:
        isp_idx = [0, int(1E10)]
    elif len(isp_idx) != 2:
        raise ValueError('ISP_IDX length must be 2')

    # read binary streams
    head_ccsds_bin, head_time_stamp_bin, head_pri_bin, head_sec_bin, user_data_bin = read_bin_ais(path_bin, isp_idx)

    # parse headers
    head_ccsds, head_time_stamp, head_pri, head_sec = parse_bin_ais_headers(head_ccsds_bin, head_time_stamp_bin, head_pri_bin, head_sec_bin)

    # Clean useless fieldnames
    head_time_stamp_fnames = list(head_time_stamp.keys())
    for fname in head_time_stamp_fnames:
        if 'spare' in fname:
            del head_time_stamp[fname]

    head_sec_fnames = list(head_sec.keys())
    for fname in head_sec_fnames:
        if 'spare' in fname:
            del head_sec[fname]

    # Create complex IQ streams
    ais_ch0_pol_H, ais_ch0_pol_V, ais_ch1_pol_H, ais_ch1_pol_V, ais_ch2_pol_H, ais_ch2_pol_V, ais_ch3_pol_H, ais_ch3_pol_V = parse_ais_user_data(head_sec, user_data_bin)

    # times stamp in decimal
    time_bits = 2.0 ** np.arange(31, -14, -1)
    array = np.array([list(format(time, '045b')) for time in head_time_stamp['sc_pps_time']])
    array = array.astype(int)
    head_time_dec = np.sum(array * time_bits, axis=1)

    # Check for decoding correctness
    if np.any(head_pri['apid'] != 7):
        print('Possible reading error. Application Process Identifier (APID) is not equal to 0x007')

    return head_ccsds, head_time_dec, head_pri, head_sec, ais_ch0_pol_H, ais_ch0_pol_V, ais_ch1_pol_H, ais_ch1_pol_V, ais_ch2_pol_H, ais_ch2_pol_V, ais_ch3_pol_H, ais_ch3_pol_V
