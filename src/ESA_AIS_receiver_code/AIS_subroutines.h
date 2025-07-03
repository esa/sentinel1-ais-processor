/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/////////////////////////////////////////////////////////////////////////
// File: AIS_subroutines.h
// Author: Tommaso Foggi
/////////////////////////////////////////////////////////////////////////
// DESCRIPTION:
//
// This file contains the definition of all function used to implement
// AIS receiver.
// For almost each function, a floating-point and a fixed-point
// implementation has been defined, named <function_name> and
// <function_name>_q respectively. Detailed comments are present only for
// the fist one, the only difference being the presence of quantization
// parameters for input and output variables.
//
// Related documents:
//   "ITU-R M.1371-3"
//   Patent: PCT/EP2014/051273
//   RECEIVING METHOD AND RECEIVER FOR TIMING AND FREQUENCY OFFSET CORRECTION
//   OF CONTINUOUS PHASE DEMODULATION IN SATELLITE-BASED AUTOMATIC
//   IDENTIFICATION SYSTEMS
/////////////////////////////////////////////////////////////////////////
#ifndef AIS_SUBROUTINES_H
#define AIS_SUBROUTINES_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h> //#CGS_2.4
#include <stdint.h>

#define PI 3.141592653589793238463 // Pi
#define CORDIC_SCALINGFACTOR 1.647 // cordic scaling factor
#define CORDIC_ITER 15			   // number of cordic iteration
#define MINVAL 1e-20			   // #TF_2.3

#define DISABLE_FLAG_ALG 0

#define ZONAL_ADD 1	   // #CGS_2.4
#define ZONAL_REMOVE 1 // #CGS_2.4

/////////////////////////////////////////////////////////////////////////
// This function describes the behavior of the whole AIS receiver
// data_in = received vector; when function ends, this vector contains
// values for received vector - detected frame
// bit_out = detected bit vector
// quant_pre = pre-detection quantization flag: must be set to 1 if
// quantization is required, 0 otherwise
// quant_det = detection quantization flag: must be set to 1 if
// quantization is required, 0 otherwise
// quant_post = post-detection quantization flag: must be set to 1 if
// quantization is required, 0 otherwise
// H1 = low-pass filter impulse response
// H2 = matched filter impulse response
// pulse = CPM pulse
// NcT = number of samples per symbol period (Tx side)
// NcR = number of samples per symbol period (Rx side)
// frame_len = frame length
// Ta = ramp-up duration (# bits)
// Tf = ramp-down duration, part 1 (# bits)
// Tg = ramp-down duration, part 1 (# bits)
// max_delay = number of possible starting positions
// max_stuffing = maximum number of stuffing bits
// frame_len = packet length
// data_len = data field length
// Nsimb =
// lung = matched filter length, expressed as number of symbol periods
// M = constellation order
// num_h = modulation index's numerator
// p = modulation index's denominator
// L = frequency pulse length
// ru_scaling = ramp-up scaling factor
// rd_scaling = ramp-down scaling factor
// The function returns 1 if a valid frame has been detected, 0
// otherwise.
/////////////////////////////////////////////////////////////////////////
// #CGS_2.7 added parameters with error values read from .dat files
int ais_receiver(double complex *data_in, int *bit_out,
				 float *H1, float *H2, float *H3, float *H4, float *pulse, int NcR,
				 int Ta, int Tf, int Tg, int max_delay, int max_stuffing,
				 int pck_len, int frame_len, int data_len,
				 int Nsimb, int lung, int L,
				 int *bit, int *det_pos, double complex *cCancelMessage,
				 int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
				 int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
				 int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin, float *f_out);

// NAME: ZonalFloatingPoint
// PURPOSE:   Add the zonal in floating point
// *INPUT:		double complex * 		data_in				Input data
// *			AISParamets 	sAISParamets		Structure filled with input parameters
// *			double complex *		cFasorStim			Fasor
// *OUTPUT:  	double complex * data_out 					Data with zonal application
void ZonalFloatingPoint(double complex *data_in_pre, double complex *data_out, int pck_len, float zonalValue, int nNcR, int nZonalType);

void remodulateDetection(double complex *data_in, int *bit_out, int *detected, int *diff_dec, int pck_len, int NcR, int detected_frame_len, int detected_position,
						 int L, float *pulse, int STFAS, float tau, float *f_out, double complex *cCancelMessage);

// #CGS_2.4
/////////////////////////////////////////////////////////////////////////
// CPM modulation function.
// state = current state
// symb = current symbol
// M = constellation order
// num_h = modulation index's numerator
// p = modulation index's denominator
// L = frequency pulse length
// n = index
// pulse = pulse vector
// NcT = number of samples per symbol
// s1 = sample portion independent from state
// s2 = sample portion dependent on state (e^(j*phase_state))
// s = calculated sample (s1 * s2)
// See TN1 for details, equations 1.1 and 1.5
/////////////////////////////////////////////////////////////////////////
void mod_CPM(int state, int symb, int M, int num_h, int p, int L, int n,
			 float *pulse, int NcT, double complex *s1, double complex *s2, double complex *s);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for mod_CPM_q function
/////////////////////////////////////////////////////////////////////////
void mod_cpm_q_c_(int *state, int *symb, int *M, int *num_h, int *p,
				  int *L, int *n, float *pulse, int *NcT,
				  float *s1_r, float *s1_i, float *s2_r, float *s2_i,
				  float *s_r, float *s_i,
				  int *p_out, int *q_out, int *p_pulse, int *q_pulse);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of mod_CPM function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void mod_CPM_q(int state, int symb, int M, int num_h, int p, int L, int n,
			   float *pulse, int NcT, double complex *s1, double complex *s2, double complex *s,
			   int p_out, int q_out, int p_pulse, int q_pulse);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for next_state function
/////////////////////////////////////////////////////////////////////////
void next_state_c_(int *state, int *new_state, int *symb, int *M, int *p,
				   int *L);

/////////////////////////////////////////////////////////////////////////
// This function is used during CPM modulation to calculate next CPM
// state
// state = current state
// new_state = new state
// symb = current symbol
// M = constellation order
// p = modulation index's denominator
// L = frequency pulse length
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void next_state(int state, int *new_state, int symb, int M, int p, int L);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for filter function
/////////////////////////////////////////////////////////////////////////
void filter_c_(float *data_in_r, float *data_in_i,
			   float *data_out_r, float *data_out_i,
			   float *H, int *lung, int *pck_len,
			   int *del, int *NcT, int *istcamp);

/////////////////////////////////////////////////////////////////////////
// Downsalmpling filtering function
// data_in = input vector
// data_out = output vector
// H = filter impulse response
// lung = filter length
// pck_len = packet lenght
// del = filter delay
// NcT = number of samples per symbol
// istcamp = sampling time
/////////////////////////////////////////////////////////////////////////
void filter(double complex *data_in, double complex *data_out, float *H, int lung, int pck_len, int del,
			int NcT, int istcamp);

///////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for filter_q function
///////////////////////////////////////////////////////////////////////
void filter_q_c_(float *data_in_r, float *data_in_i,
				 float *data_out_r, float *data_out_i,
				 float *H, int *lung, int *pck_len,
				 int *del, int *NcT, int *istcamp,
				 int *p_h, int *q_h, int *p_in, int *q_in,
				 int *p_out, int *q_out);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of filter function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void filter_q(double complex *data_in, double complex *data_out, float *H, int lung, int pck_len,
			  int del, int NcT, int istcamp, int p_h, int q_h,
			  int p_in, int q_in, int p_out, int q_out);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for filter_nd function
/////////////////////////////////////////////////////////////////////////
void filter_nd_c_(float *data_in_r, float *data_in_i,
				  float *data_out_r, float *data_out_i,
				  float *H, int *lung, int *pck_len, int *del, int *NcT);

/////////////////////////////////////////////////////////////////////////
// Non-downsalmpling filtering function
// data_in = input vector
// data_out = output vector
// H = filter impulse response
// lung = filter length
// pck_len = packet lenght
// del = filter delay
// NcT = number of samples per symbol
/////////////////////////////////////////////////////////////////////////
void filter_nd(double complex *data_in, double complex *data_out, float *H, int lung, int pck_len,
			   int del, int NcT);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for filter_nd_q function
/////////////////////////////////////////////////////////////////////////
void filter_nd_q_c_(float *data_in_r, float *data_in_i,
					float *data_out_r, float *data_out_i,
					float *H, int *lung, int *pck_len,
					int *del, int *NcT, int *p_h, int *q_h,
					int *p_in, int *q_in, int *p_out, int *q_out);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of filter_nd function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void filter_nd_q(double complex *data_i, double complex *data_o, float *H, int lung, int pck_len,
				 int del, int NcT, int p_h, int q_h, int p_in, int q_in,
				 int p_out, int q_out);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for interpolation function
/////////////////////////////////////////////////////////////////////////
void interpolation_c_(float *data_in_r, float *data_in_i,
					  float *data_out_r, float *data_out_i,
					  int *pck_len, int *NcR, float *ist, int *istcamp);

/////////////////////////////////////////////////////////////////////////
// Quadratic interpolation function
// data_in = input vector
// data_out = interpolated vector
// pck_len = packet lenght
// NcR = number of samples per symbol
// ist = resampling time
// istcamp = resampling time (integer value)
/////////////////////////////////////////////////////////////////////////
void interpolation(double complex *data_in, double complex *data_out, int pck_len,
				   int NcR, float ist, int *istcamp);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for interpolation_q function
/////////////////////////////////////////////////////////////////////////
void interpolation_q_c_(float *data_in_r, float *data_in_i,
						float *data_out_r, float *data_out_i,
						int *pck_len, int *NcR, float *ist, int *istcamp,
						int *p_out, int *q_out);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of interpolation function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void interpolation_q(double complex *data_in, double complex *data_out, int pck_len,
					 int NcR, float ist, int *istcamp, int p_out, int q_out);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for interpolation2 function
/////////////////////////////////////////////////////////////////////////
void interpolation2_c_(float *data_in_r, float *data_in_i,
					   float *data_out_r, float *data_out_i,
					   int *pck_len, int *NcR, float *ist);

/////////////////////////////////////////////////////////////////////////
// Quadratic interpolation function, different from interpolation in
// terms of interpolation window width
// data_in = input vector
// data_out = interpolated vector
// pck_len = packet lenght
// ist = resampling time
/////////////////////////////////////////////////////////////////////////
void interpolation2(double complex *data_in, double complex *data_out, int pck_len,
					int NcR, float ist);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for interpolation2_q function
/////////////////////////////////////////////////////////////////////////
void interpolation2_q_c_(float *data_in_r, float *data_in_i,
						 float *data_out_r, float *data_out_i,
						 int *pck_len, int *NcR, float *ist,
						 int *p_out, int *q_out);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of interpolation2 function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void interpolation2_q(double complex *data_in, double complex *data_out, int pck_len,
					  int NcR, float ist, int p_out, int q_out);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for nd_timfreq function
/////////////////////////////////////////////////////////////////////////
void nda_timfreq_c_(float *data_in_r, float *data_in_i, int *NcT,
					int *pck_len, int *offset, int *M, int *L,
					float *ist_re, int *ist_int, float *fstim);

/////////////////////////////////////////////////////////////////////////
// This function is used fot timing and frequency esteem in pre-detection
// procedure
// data_in = received vector after low-pass filtering
// NcT = number of samples per symbol
// pck_len = packet lenght
// offset = position of the first usful bit of data_in
// M = number of correlation terms (30)
// L = signaling intervals (128)
// ist_re = optimum sampling time (real value)
// ist_int = optimum sampling time (integer value)
// fstim = estimated frequency
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void nda_timfreq(double complex *data_in, int NcT, int pck_len, int offset,
				 int M, int L, float *ist_re, int *ist_int,
				 float *fstim);

/////////////////////////////////////////////////////////////////////////
// This function is used fot timing and frequency esteem in pre-detection
// procedure, considering only even elements of correlation
// data_in = received vector after low-pass filtering
// NcT = number of samples per symbol
// pck_len = packet lenght
// offset = position of the first usful bit of data_in
// M = number of correlation terms (30)
// L = signaling intervals (128)
// ist_re = optimum sampling time (real value)
// ist_int = optimum sampling time (integer value)
// fstim = estimated frequency
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void nda_timfreq2T(double complex *data_in, int NcT, int pck_len, int offset,
				   int M, int L, float *ist_re, int *ist_int,
				   float *fstim);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for nd_timfreq2T function
/////////////////////////////////////////////////////////////////////////
void nda_timfreq2T_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *M, int *L,
					  float *ist_re, int *ist_int, float *fstim);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for nd_timfreqMCM function
/////////////////////////////////////////////////////////////////////////
void nda_timfreqMCM_c_(float *data_in_r, float *data_in_i, int *NcT,
					   int *pck_len, int *offset, int *L,
					   int *ist_int, float *fstim);

/////////////////////////////////////////////////////////////////////////
// This function is used for frequency esteem in pre-detection
// procedure, first algorithm
// data_in = received vector after low-pass filtering
// NcT = number of samples per symbol
// pck_len = packet lenght
// offset = position of the first usful bit of data_in
// L = signaling intervals (128)
// ist_int = optimum sampling time (integer value)
// fstim = estimated frequency
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void nda_timfreqMCM(double complex *data_in, int NcT, int pck_len, int offset,
					int L, int *ist_int,
					float *fstim);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for nd_timfreqMM function
/////////////////////////////////////////////////////////////////////////
void nda_timfreqMM_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *L,
					  int *ist_int, float *fstim);

/////////////////////////////////////////////////////////////////////////
// This function is used for frequency esteem in pre-detection
// procedure, second algorithm
// data_in = received vector after low-pass filtering
// NcT = number of samples per symbol
// pck_len = packet lenght
// offset = position of the first usful bit of data_in
// L = signaling intervals (128)
// ist_int = optimum sampling time (integer value)
// fstim = estimated frequency
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void nda_timfreqMM(double complex *data_in, int NcT, int pck_len, int offset,
				   int L, int *ist_int,
				   float *fstim);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for nda_timfreq_q function
/////////////////////////////////////////////////////////////////////////
void nda_timfreq_q_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *M, int *L,
					  float *ist_re, int *ist_int, float *fstim,
					  int *p_in, int *q_in, int *p_f, int *q_f);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of nda_timfreq function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void nda_timfreq_q(double complex *data_in, int NcT, int pck_len, int offset,
				   int M, int L, float *ist_re, int *ist_int,
				   float *fstim, int p_in, int q_in, int p_f, int q_f);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for BCJR function
/////////////////////////////////////////////////////////////////////////
void bcjr_c_(float *data_in_r, float *data_in_i, float *LLR, int *AD,
			 int *pck_len, float *Pd, int *D);

/////////////////////////////////////////////////////////////////////////
// This function is used to perform post-detection frequency estimation
// data_in = input data after pre-detection and matching filtering
// LLR = detected log-likelihood ratios
// AD = detected bits
// pck_len = packet lenght
// Pd = Pdelta in equation 1.23 (TN1)
// D = Q parameter in equation 1.23 (TN1)
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void BCJR(double complex *data_in, float *LLR, int *AD, int pck_len,
		  float Pd, int D);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for BCJR_q function
/////////////////////////////////////////////////////////////////////////
void bcjr_q_c_(float *data_in_r, float *data_in_i, float *LLR, int *AD,
			   int *pck_len, float *Pd, int *D, int *p_in, int *q_in,
			   int *p_llr, int *q_llr);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of BCJR function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void BCJR_q(double complex *data_in, float *LLR, int *AD, int pck_len,
			float Pd, int D, int p_in, int q_in, int p_llr, int q_llr);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimafreqMM function
/////////////////////////////////////////////////////////////////////////
void stimafreqmm_c_(float *data_in1_r, float *data_in1_i,
					float *data_in2_r, float *data_in2_i,
					int *offset, float *freq, int *NcT, int *pck_len);

/////////////////////////////////////////////////////////////////////////
// This function is used to perform post-detection frequency estimation,
// Mengalli-Morelli algorithm
// data_in1 = received vector
// data_in2 = remodulated frame
// NcT = number of samples per symbol
// pck_len = packet lenght
/////////////////////////////////////////////////////////////////////////
void stimafreqMM(double complex *data_in1, double complex *data_in2, float *freq,
				 int offset, int NcT, int pck_len);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimafreqMM_q function
/////////////////////////////////////////////////////////////////////////
void stimafreqmm_q_c_(float *data_in1_r, float *data_in1_i,
					  float *data_in2_r, float *data_in2_i,
					  float *freq, int *offset, int *NcT, int *pck_len,
					  int *p_in1, int *q_in1, int *p_in2,
					  int *q_in2, int *p_f, int *q_f);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of stimafreqMM function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void stimafreqMM_q(double complex *data_in1, double complex *data_in2, float *freq,
				   int offset, int NcT, int pck_len, int p_in1, int q_in1,
				   int p_in2, int q_in2, int p_f, int q_f);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimatiming function
/////////////////////////////////////////////////////////////////////////
void stimatiming_c_(float *data_in1_r, float *data_in1_i,
					float *data_in2_r, float *data_in2_i,
					float *time, int *offset, int *NcT,
					int *pck_len);

/////////////////////////////////////////////////////////////////////////
// This function is used to esteem timing for post-detection operations
// data_in1 = received vector
// data_in2 = received vector after post-processing operation
// time = estimated timing
// NcT = number of samples per symbol
// pck_len = packet lenght
// See report, page 11
/////////////////////////////////////////////////////////////////////////
void stimatiming(double complex *data_in1, double complex *data_in2, float *time,
				 int offset, int NcT, int pck_len);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimatiming_q function
/////////////////////////////////////////////////////////////////////////
void stimatiming_q_c_(float *data_in1_r, float *data_in1_i,
					  float *data_in2_r, float *data_in2_i,
					  float *time, int *offset, int *NcT,
					  int *pck_len, int *p_in1, int *q_in1,
					  int *p_in2, int *q_in2,
					  int *p_t, int *q_t);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of stimatiming function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void stimatiming_q(double complex *data_in1, double complex *data_in2, float *time,
				   int offset, int NcT, int pck_len, int p_in1,
				   int q_in1, int p_in2, int q_in2, int p_t, int q_t);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimafasamp function
/////////////////////////////////////////////////////////////////////////
void stimafasamp_c_(float *data_in1_r, float *data_in1_i,
					float *data_in2_r, float *data_in2_i,
					float *ejfasestim_r, float *ejfasestim_i,
					float *ampstim, int *offset, int *NcT, int *pck_len, int *Ta, int *Tf, int stfas);

/////////////////////////////////////////////////////////////////////////
// This function is used to esteem phase and amplitude for post-detection
// operations
// data_in1 = received vector
// data_in2 = received vector after post-detection operation
// NcT = number of samples per symbol
// pck_len = packet lenght
//
// See TN1 for details
/////////////////////////////////////////////////////////////////////////
void stimafasamp(double complex *data_in1, double complex *data_in2,
				 double complex *ejfasestim, float *ampstim,
				 int offset, int NcT, int pck_len, int Ta, int Tf, int stfas);

/////////////////////////////////////////////////////////////////////////
// Wrapper for fortran simulator interface for stimafasamp_q function
/////////////////////////////////////////////////////////////////////////
void stimafasamp_q_c_(float *data_in1_r, float *data_in1_i,
					  float *data_in2_r, float *data_in2_i,
					  float *ejfasestim_r, float *ejfasestim_i,
					  float *ampstim, int *offset, int *NcT,
					  int *pck_len, int *Ta, int *Tf, int stfas, int *p_data1, int *q_data1,
					  int *p_data2, int *q_data2, int *p_rit,
					  int *q_rit, int *p_fas, int *q_fas,
					  int *p_amp, int *q_amp);

/////////////////////////////////////////////////////////////////////////
// This is the quantized version of stimafasamp function
// p_* and q_* are quantization parameters for input vectors and output
// data
/////////////////////////////////////////////////////////////////////////
void stimafasamp_q(double complex *data_in1, double complex *data_in2,
				   double complex *ejfasestim, float *ampstim, int offset,
				   int NcT, int pck_len, int Ta, int Tf, int stfas, int p_data1, int q_data1,
				   int p_data2, int q_data2, int p_rit, int q_rit,
				   int p_fas, int q_fas, int p_amp, int q_amp);

/////////////////////////////////////////////////////////////////////////
// This function analyzes detected frame to check if it is a valid frame.
// The procedure consists of:
// - remove bit stuffing
// - CRC check
// - check start/stop flag (only if valid CRC)
// - if these fields contain up to 1 errors, these errors are corrected
// and the frame is considered detected
// - otherwise, bit flipping is performed (up to 2 couples flipped)
// - if no frame is detected, the procedure restarts, considering another
// starting position
// frame = detected bits
// LLR = detected bits' LLR
// frame_len = frame length
// data_len = data field length
// max_delay = number of possible starting positions
// max_stuffing = maximum number of stuffing bits
// detected_position = frame starting bit
// detected_frame_len = length of the detected frame
// The function returns 1 in case of successfull detection, 0 otherwise
/////////////////////////////////////////////////////////////////////////
// #CGS_2.7 added parameters with error values read from .dat files
int check_detected(int *frame, float *LLR, int frame_len, int data_len,
				   int max_delay, int max_stuffing, int Ta,
				   int *detected_position, int *detected_frame_len, int T_len, int Icrc,
				   int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
				   int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
				   int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin); // #TF_2.3 //#CGS_2.4

/////////////////////////////////////////////////////////////////////////
// This function removes the differential enconding from LLRs //TF_2.6
/////////////////////////////////////////////////////////////////////////
void DiffLLR(int all_len, float *LLR, float *LLRD);

/////////////////////////////////////////////////////////////////////////
// This function is used by DiffLLR //TF_2.6
/////////////////////////////////////////////////////////////////////////
float transform(float segno, float modulo);

/////////////////////////////////////////////////////////////////////////
// This function computes the CRC of data and crc together //TF_2.6
/////////////////////////////////////////////////////////////////////////
void crc_check_tot(int lt, int *msgDEC, int sum, int *crcout);

/////////////////////////////////////////////////////////////////////////
// This function convertes binary a sequence to integer //TF_2.6
/////////////////////////////////////////////////////////////////////////
int natConv(int a, int *bin);

/////////////////////////////////////////////////////////////////////////
// This function searches error couples in syndrome files //TF_2.6
/////////////////////////////////////////////////////////////////////////
int FindIndex(int a, int n, int *vett);

/////////////////////////////////////////////////////////////////////////
// This function reads possible syndrome errors from files //TF_2.6
// #CGS_2.7 change return value from void to integer (1 OK, 0 errors)
/////////////////////////////////////////////////////////////////////////
int readfile(int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c, int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2, int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin);

/////////////////////////////////////////////////////////////////////////
// This function checks the CRC field of a frame for second postprocessing //TF_2.6
/////////////////////////////////////////////////////////////////////////
void crc_check_3(int lt, int l, int *msgDEC, int *sum); // #CGS_2.7 sum passed by reference

/////////////////////////////////////////////////////////////////////////
// This function analyzes controls CRC, start and stop flag of input
// frame.
// frame = input frame
// frame_len = frame length
// data_len = data field length
// stuffing_bits = number of stuffing bits detected
// The function returns 1 in case of correct frame, 0 otherwise. It also
// writes the number of stuffing bits detected
/////////////////////////////////////////////////////////////////////////
int check_frame(int *frame, int frame_len, int data_len,
				int *stuffing_bits);

/////////////////////////////////////////////////////////////////////////
// This function analyzes controls CRC of input, alternative processing //TF_2.6
// frame = input frame
// frame_len = frame length
// data_len = data field length
// max_stuffing = maximum number of stuffing bits
// T_len = ramps length
// stuffing_bits = number of stuffing bits detected
// The function returns 1 in case of correct frame, 0 otherwise. It also
// writes the number of stuffing bits detected
///////////////////////////////////////////////////////////////////////

// #CGS_2.7 added parameters with error values read from .dat files
// int check_frame_2(int *frame, float *LLRT, int frame_len, int data_len, int max_delay, int max_stuffing, int T_len, int *stuffing_bits);
int check_frame_2(int *frame, float *LLRT, int frame_len, int data_len, int max_delay, int max_stuffing, int T_len, int *stuffing_bits,
				  int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
				  int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
				  int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin);

/////////////////////////////////////////////////////////////////////////
// This function checks the CRC field of a frame.
// CRC calculation is based on generator polinomial
// g(x) = x^16+x^12+x^5+1
// pnFrame = input frame
// L_Frame = frame length
// The function returns 1 in case of correct CRC, 0 otherwise
/////////////////////////////////////////////////////////////////////////
int CRC_check(int pnFrame[], int L_Frame);

/////////////////////////////////////////////////////////////////////////
// This function checks the CRC field of a frame.
// CRC calculation is based on generator polinomial
// g(x) = x^16+x^12+x^5+1
// pnFrame = input frame
// L_Frame = frame length
// The function returns 1 in case of correct CRC, 0 otherwise
/////////////////////////////////////////////////////////////////////////
int CRC_check_2(int pnFrame[], int L_Frame);

/////////////////////////////////////////////////////////////////////////
// This function calculates the CRC field of a frame.
// CRC calculation is based on generator polinomial
// g(x) = x^16+x^12+x^5+1
// pnFrame = input frame
// CRC = calculated CRC
// L_Frame = frame length
// CRC_len = length of CRC field
/////////////////////////////////////////////////////////////////////////
void CRC_calc_2(int pnFrame[], int CRC[], int L_Frame, int CRC_len);

/////////////////////////////////////////////////////////////////////////
// This function calculates the CRC field of a frame.
// CRC calculation is based on generator polinomial
// g(x) = x^16+x^12+x^5+1
// pnFrame = input frame
// L_Frame = frame length
// CRC data appended at the end of pnFrame
/////////////////////////////////////////////////////////////////////////
void CRC_calc(int pnFrame[], int L_Frame);

/////////////////////////////////////////////////////////////////////////
// The following functions are used to implement BCJR algorithm
/////////////////////////////////////////////////////////////////////////
void H_1d(double complex y, int D, float *v);

void H_2d(double complex y, int D, float **v, int i2);

void H_q_1d(double complex y, int D, float *v, int p, int q);

void H_q_2d(double complex y, int D, float **v, int i2, int p, int q);

float integr(float *v, int D);

void normalize(float **v, int D, int i2);

void var_exp_1d(float *v1, float **v2, int D, float *v, int i2);

void var_exp_2d(float *v1, float **v2, int D, float **v, int i2);

void fact_exp1(float **vin, int D, float Pd, float **vout, int iin, int iout);

void fact_exp2(float **vin, int D, float Pd, int alpha, float *vout, int ii);

void normalize_q(float **v, int D, int i2, int p, int q);

void var_exp_1d_q(float *v1, float **v2, int D, float *v, int i2, int p, int q);

void var_exp_2d_q(float *v1, float **v2, int D, float **v, int i2, int p, int q);

void fact_exp1_q(float **vin, int D, float Pd, float **vout, int iin, int iout, int p, int q);

void fact_exp2_q(float **vin, int D, float Pd, int alpha, float *vout, int ii, int p, int q);

float logsumexp(float d1, float d2);

/////////////////////////////////////////////////////////////////////////
// This function is used to quantize a floating point
// number on p+q bits over a quasi-simmetric range
// (-2^(p-1); 2^(p-1)-1/2^q).
// p-1 bits are used for the integer part
// 1 bit for the sign
// q bits for the decimal part
// data_in = data to be quantized
// p = number of bits used to quantize integer part
// q = number of bits used to quantize decimal part
//
// Funcion procedure:
// multiply by 2^q
// truncate to 0 bits after point
// saturate at -2^(p+q-1) and 2^(p+q-1)-1
// divide by 2^q
/////////////////////////////////////////////////////////////////////////
float quantize(float data_in, int p, int q);

/////////////////////////////////////////////////////////////////////////
// This function is used to calculate the phase (in
// radiants) of a double complex number C = re + j*im, using
// CORDIC algorithm
// re = vector's real part
// im = vector's imaginary part
// iterations = number of iterations to be performed
//
// Funcion procedure:
// evaluate if phase of C exceeds limits (-pi/2, pi/2)
// in this case, perform a +/- 90° rotation
// iteratively adjust vector's phase by multiplying
// per R = 1 +/- jK, K = 2**(-i), i = iteration index
//
// Values for atan(2^-i) can be precalculated
/////////////////////////////////////////////////////////////////////////
float cordic_getphase(float re, float im, int iterations);

/////////////////////////////////////////////////////////////////////////
// This function is used to calculate the magnitude
// of a double complex number C = re + j*im, using CORDIC
// algorithm
// re = vector's real part
// im = vector's imaginary part
// iterations = number of iterations to be performed
//
// Funcion procedure:
// evaluate if phase of C exceeds limits (-pi/2, pi/2)
// in this case, perform a +/- 90° rotation
// iteratively adjust vector's phase by multiplying
// per R = 1 +/- jK, K = 2**(-i), i = iteration index
// The final vector will be something like
// re' + j0
// The magnitude of the input vector is re'
// multiplied per a scaling factor (1.647), since
// every rotation modifies the vector's magnitude.
// Algorithm's overall gain converges to 1.647.
//
// Values for atan(2^-i) can be precalculated
/////////////////////////////////////////////////////////////////////////
float cordic_getmag(float re, float im, int iterations);

/////////////////////////////////////////////////////////////////////////
// This function is used to perform polar to
// rectangular conversion, using CORDIC algorithm
// magnitude = vector's magnitude
// phase = vector's phase
// iterations = number of iterations to be performed
//
// Funcion procedure:
// initialize a vector re = magnitude/scaling_factor,
// im = 0 (remember that CORDIC algorithm has a gain
// of approximately 1.647). Rotate this vector until
// the total accumulates phase equals input phase.
//
// Values for atan(2^-i) can be precalculated
/////////////////////////////////////////////////////////////////////////
double complex cordic_polartorect(float magnitude, float phase, int iterations);

/////////////////////////////////////////////////////////////////////////
// This function calculates a-INT(a/b)*b
/////////////////////////////////////////////////////////////////////////
float mod(float a, float b);

/////////////////////////////////////////////////////////////////////////
// This function returns the lower of its inputs
/////////////////////////////////////////////////////////////////////////
float Minf(float a, float b);

/////////////////////////////////////////////////////////////////////////
// This function returns the higher of its inputs
/////////////////////////////////////////////////////////////////////////
float Maxf(float a, float b);

/////////////////////////////////////////////////////////////////////////
// This function returns the lower of its inputs
/////////////////////////////////////////////////////////////////////////
float MinInt(int a, int b);

/////////////////////////////////////////////////////////////////////////
// This function returns the higher of its inputs
/////////////////////////////////////////////////////////////////////////
float MaxInt(int a, int b);

void flip_bytes(int *bitout, const int data_len);

#endif
