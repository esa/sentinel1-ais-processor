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
// File: AIS_subroutines.c
// Author: Tommaso Foggi
/////////////////////////////////////////////////////////////////////////
// DESCRIPTION:
//
// This file contains the implementation of all function used to model
// the AIS FENICE receiver.
// For almost each function, a floating-point and a fixed-point
// implementation has been defined, named <function_name> and
// <function_name>_q respectively. Detailed comments are present only for
// the fist one, the only difference being the presence of quantization
// parameters for input and output variables.
// All quantized variables are represented as floating point variables,
// but are mapped over p+q bits (1 bit = sign, p-1 bits = integer part,
// q bits = decimal part), using "quantize" function.
//
// Related documents:

//   [RD1] Patent: PCT/EP2014/051273
//   	RECEIVING METHOD AND RECEIVER FOR TIMING AND FREQUENCY OFFSET CORRECTION
//   	OF CONTINUOUS PHASE DEMODULATION IN SATELLITE-BASED AUTOMATIC IDENTIFICATION SYSTEMS
//   [RD2] "ITU-R M.1371-3"

#include "AIS_subroutines.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h> // for size_t

/////////////////////////////////////////////////////////////////////////
// #CGS_2.7 added parameters with error values read from .dat files
int ais_receiver(
	double complex *data_in,
	int *bit_out, float *H1, float *H2,
	float *H3, float *H4, float *pulse, int NcR, int Ta, int Tf,
	int Tg, int max_delay, int max_stuffing, int pck_len, int frame_len,
	int data_len, int Nsimb, int lung, int L, int *bit,
	int *det_pos, double complex *cCancelMessage,
	int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
	int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
	int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin, float *f_out)
{
	const int del = 4;	 // delay due to matching filter
	const int norm = 1;	 // normalization enable
	const int STFAS = 1; // Data aided timing and phase estimation in case of remodulation.
	const int Icrc = 2;	 // Syndome error correction on CRC data field

	// Detection algorithm parameters
	const int D = 24;		// Q parameter in equation 3.9 (RD1)
	const float Pd = 0.001; // Pdelta in equation 3.9 (RD1)
	int del1;				// low-pass filter delay

	// Data vectors
	double complex data_in_pre[pck_len * NcR];		 // data after limitator, pre-detection input
	double complex data_in_filtered[pck_len * NcR];	 // data after zonal lowpass filtering
	double complex data_pre_filtered[pck_len * NcR]; // data after first low-pass filtering
	double complex data_in_det[pck_len * NcR];		 // detector input
	double complex data_det_filtered[pck_len];		 // data after matched filtering
	float LLR[pck_len];								 // detected log-likelihood ratios
	float sum;										 // normalization
	int detected[pck_len];							 // detected bits
	int diff_dec[pck_len];							 // differentially decoded bits
	int detection_result;							 // post-processing result
	int detected_position;							 // detected frame starting point
	int detected_frame_len;							 // length of detected frame

	// Frequency estimation results
	float fstim;

	// Time estimation results
	float tau;
	int taucamp;
	int taucamp2;

	// Phase estimation results and relative phase vectors
	float angstim;
	double complex fasorestim;

	// Support variables
	int i, k;
	float t;
	int T_len = Ta + Tg + Tf;

	// Initialize arrays
	for (i = 0; i < pck_len * NcR; i++)
	{
		data_in_pre[i] = 0 + I * 0;
		data_in_filtered[i] = 0 + I * 0;
		data_pre_filtered[i] = 0 + I * 0;
		data_in_det[i] = 0 + I * 0;
		cCancelMessage[i] = 0 + I * 0;
	}

	for (i = 0; i < pck_len; i++)
	{
		data_det_filtered[i] = 0 + I * 0;
		LLR[i] = 0.0;
		detected[i] = 0;
		diff_dec[i] = 0;
	}

	del1 = Nsimb * NcR / 2 + 1;
	// normalization
	if (norm == 1)
	{
		sum = 0;
		for (i = 0; i < pck_len * NcR; i++)
		{
			sum += pow(cabs(data_in[i]), 2) / (pck_len * NcR);
		}
		for (i = 0; i < pck_len * NcR; i++)
		{
			data_in[i] = data_in[i] / sqrt(sum);
		}
	}

	// Pre-detection
	// Low-pass filtering
	filter_nd(data_in, data_in_filtered, H3, Nsimb, pck_len, del1, NcR);

	// Apply limiter
	for (i = 0; i < pck_len * NcR; i++)
	{
		if (cabs(data_in_filtered[i]) > MINVAL)
			data_in_pre[i] = data_in_filtered[i] / cabs(data_in_filtered[i]);

		if (cabs(data_in_filtered[i]) <= MINVAL)
			data_in_pre[i] = MINVAL + I * MINVAL;
	}

	// Low-pass filtering
	filter_nd(data_in_pre, data_pre_filtered, H1, Nsimb, pck_len, del1, NcR);

	// Time and frequency estimation
	nda_timfreq2T(data_pre_filtered, NcR, pck_len, max_delay + Ta, 20,
				  128, &tau, &taucamp, &fstim); // #CGS_2.4

	interpolation(data_in_pre, data_in_det, pck_len, NcR, tau, &taucamp2);

	for (i = 0; i < pck_len * NcR; i++)
	{
		data_in_pre[i] = data_in_det[i];
	}

	// Low-pass filtering
	filter_nd(data_in_pre, data_pre_filtered, H4, Nsimb, pck_len, del1, NcR);

	// Time and frequency estimation
	nda_timfreqMM(data_pre_filtered, NcR, pck_len, max_delay + Ta, 128, &taucamp2, &fstim); // #CGS_2.4

	// Correct input vector
	angstim = 0.0;
	for (i = 0; i < pck_len * NcR; i++)
	{
		angstim = mod(angstim + PI * 2 * fstim / (float)NcR, PI * 2);
		fasorestim = cos(angstim) - I * sin(angstim);
		data_in_det[i] = data_in_pre[i] * fasorestim;
	}

	// Matching filter and decimation
	filter(data_in_det, data_det_filtered, H2, lung, pck_len, del, NcR, taucamp2);

	// Detection
	BCJR(data_det_filtered, LLR, diff_dec, pck_len, Pd, D);

	// Differential decoding
	diff_dec[0] = 0;
	for (i = 1; i < pck_len; i++)
	{
		detected[i] = mod(diff_dec[i] + diff_dec[i - 1] + 1, 2);
	}

	// Post-processing: In case of successfull check, detected_position indicates the position of
	// the correct frame and detected_frame_len indicates the length of detected frame
	detection_result = check_detected(detected, LLR, frame_len, data_len,
									  max_delay, max_stuffing, Ta, &detected_position,
									  &detected_frame_len, T_len, Icrc,
									  err_sing, err_sing_sin, err_doppi, err_doppi_sin, err_sing2c,
									  err_sing2c_sin, err2c1, err2c1_sin, err3, err3_sin, err2c2,
									  err2c2_sin, err2c3, err2c3_sin, err2c4, err2c4_sin);

	f_out[0] = fstim; // Message doppler frequency output to save
	*det_pos = detected_position;

	//  Post-detection; only in case of successfull post-processing
	if (detection_result == 1)
	{
		// Build output bit vector
		for (i = 0; i < detected_frame_len; i++)
		{
			bit_out[i] = detected[detected_position + i];
		}

		// Successive interference cancellation currently disabled as it doesn't improve performance.
		// remodulateDetection(data_in, bit_out, detected, diff_dec, pck_len, NcR, detected_frame_len, detected_position,
		//					M, num_h, p, L, pulse, STFAS, tau, f_out, cCancelMessage);

		return 1; // when a valid frame has been detected
	}
	return 0; // when no frame has been detected
}

void remodulateDetection(double complex *data_in, int *bit_out, int *detected, int *diff_dec, int pck_len, int NcR, int detected_frame_len, int detected_position,
						 int L, float *pulse, int STFAS, float tau, float *f_out, double complex *cCancelMessage)
{

	int state, state_new;
	double complex s1[NcR];
	double complex s2;
	double complex ss[NcR];
	int i, k;
	int cpm_input[pck_len];						 // input vector to CPM re-modulator
	int recod[pck_len];							 // re-encoded data
	double complex remodulated[pck_len * NcR];	 // remodulated vector
	double complex data_in_post[pck_len * NcR];	 // post-detection input (after first interpolation)
	double complex data_in_post2[pck_len * NcR]; // post-detection input (after second interpolation)

	// Amplitude estimation result
	float Astim2[pck_len * NcR];
	float angstim2;
	double complex fasorestim2;
	double complex ejthetastim2[pck_len * NcR];
	int align;
	int ritc;
	float t;
	int sumt;
	const int Ta = 8; // 6; // ramp-up duration (# bits) //#CGS_2.4
	const int Tf = 8; // 8; // ramp-down duration, part 1 (# bits)
	const int Tg = 0; // 17; // ramp-down duration, part 2 (# bits)
	int T_len = Ta + Tg + Tf;
	int max_delay = 86;
	float fstim2;

	const int M = 2;	 // constellation order
	const int num_h = 1; // modulation index's numerator
	const int p = 2;	 // modulation index's denominator

	// Time estimation results
	int taudumb;
	float tau2;

	// Initialize arrays
	for (i = 0; i < pck_len * NcR; i++)
	{
		data_in_post[i] = 0 + I * 0;
		data_in_post2[i] = 0 + I * 0;
		ejthetastim2[i] = 0 + I * 0;
		Astim2[i] = 0.0;
	}

	for (i = 0; i < pck_len; i++)
	{
		cpm_input[i] = 0;
		recod[i] = 0;
	}

	if (detected_position - Ta < 0 || detected_position > max_delay)
	{
		return; // If out of bounds, return.
	}

	recod[detected_position] = mod(detected[detected_position] + 0, 2);
	for (i = detected_position + 1; i < detected_frame_len + 1; i++)
	{
		recod[i] = mod(detected[i] + recod[i - 1] + 1, 2);
	}
	sumt = 0;
	for (i = detected_position + 1; i < detected_frame_len + 1; i++)
	{
		sumt += mod(recod[i] + diff_dec[i], 2);
	}
	// Determine starting bit zero or one.
	int start_bit = 0;
	if (sumt > 10)
		start_bit = 1;

	// Accounts for differential encoding
	cpm_input[0] = mod(bit_out[0] + start_bit, 2);
	for (i = 1; i < detected_frame_len + 1; i++)
	{
		cpm_input[i] = mod(bit_out[i] + cpm_input[i - 1] + 1, 2);
	}

	// clean remodulated vector
	for (i = 0; i < pck_len * NcR; i++)
	{
		remodulated[i] = 0.0f + I * 0.0; // #TF_2.3
	}

	// remodulation
	state = 0;
	for (i = 1; i <= detected_frame_len + 1; i++)
	{
		mod_CPM(state, cpm_input[i - 1], M, num_h, p, L, i, pulse, NcR,
				s1, &s2, ss);
		next_state(state, &state_new, cpm_input[i - 1], M, p, L);
		state = state_new;
		for (k = 0; k < NcR; k++)
		{
			remodulated[(detected_position)*NcR + (i - 1) * NcR + k] = ss[k];
		}
	}

	align = rint(tau);
	t = tau - align;
	ritc = (int)(t * NcR);

	// Translation to correct the integer part of the timing
	for (i = pck_len * NcR - 1; i > 2 * NcR; i--) // #TF_2.3
	{
		if (align == 0)
		{
			remodulated[i] = remodulated[i - ritc];
		}
		else
		{
			remodulated[i] = remodulated[i - NcR - ritc];
		}
	}

	// build ramp-up: linear increment
	angstim2 = atan2(cimag(remodulated[(detected_position + 1) * NcR]),
					 creal(remodulated[(detected_position + 1) * NcR]));
	for (i = Ta * NcR; i > 0; i--)
	{
		if ((detected_position)*NcR - 1 - i >= 0)
		{
			angstim2 = mod(angstim2 + PI * 2 / 12, PI * 2); // mod(atan2(cimag(data_in[detected_position*NcR-Ta*NcR+i]),creal(data_in[detected_position*NcR-Ta*NcR+i])), PI*2);
			fasorestim2 = cos(angstim2) + I * sin(angstim2);
			remodulated[(detected_position + 1) * NcR - 1 - Ta * NcR + i] = fasorestim2; // #CGS_2.4 scommentata la riga
		}
	}

	// build ramp-down: linear decrement
	angstim2 = atan2(cimag(remodulated[(detected_position + detected_frame_len) * NcR - 1]), creal(
																								 remodulated[(detected_position + detected_frame_len) * NcR - 1]));
	for (i = 0; i < Tf * NcR; i++)
	{
		if ((detected_position + detected_frame_len) * NcR + i < pck_len * NcR)
		{
			angstim2 = mod(angstim2 + PI * 2 / 12, PI * 2); // atan2(cimag(data_in[(detected_position+detected_frame_len)*NcR+i]),creal(data_in[(detected_position+detected_frame_len)*NcR+i]));
			fasorestim2 = cos(angstim2) + I * sin(angstim2);
			remodulated[(detected_position + detected_frame_len) * NcR + i] = fasorestim2; // #CGS_2.4 scommentata la riga
		}
	}

	// build last tail: flat
	for (i = 0; i < Tg * NcR; i++)
	{
		if ((detected_position + detected_frame_len + Tf) * NcR + i < pck_len * NcR)
		{
			remodulated[(detected_position + detected_frame_len + Tf) * NcR + i] = remodulated[(detected_position + detected_frame_len + Tf) * NcR - 1];
		}
	}

	// interpolation
	interpolation(remodulated, data_in_post, pck_len, NcR,
				  (float)((int)(NcR * (tau - rint(align)))) / (float)(NcR) - (tau - rint(align)), &taudumb);

	// frequency estimation and correction
	stimafreqMM(data_in, data_in_post, &fstim2, detected_position, NcR,
				detected_frame_len);

	angstim2 = 0.0;
	for (i = NcR; i < pck_len * NcR; i++) // #TF_2.3
	{
		angstim2 = mod(angstim2 + PI * 2 * fstim2 / (float)NcR, PI * 2);
		fasorestim2 = cos(angstim2) + I * sin(angstim2);
		data_in_post[i] = data_in_post[i] * fasorestim2;
	}

	// refine timing estimation
	stimatiming(data_in, data_in_post, &tau2, detected_position,
				NcR, detected_frame_len);

	// interpolation
	interpolation(data_in_post, data_in_post2, pck_len, NcR, tau2,
				  &taudumb);

	// phase and amplitude estimation and correction
	stimafasamp(data_in, data_in_post2, ejthetastim2, Astim2,
				detected_position, NcR, detected_frame_len, Ta, Tf, STFAS);

	for (i = pck_len * NcR - 1; i > (detected_position - Ta) * NcR; i--) // #TF_2.3
	{
		ejthetastim2[i] = ejthetastim2[i - (detected_position - Ta) * NcR];
		Astim2[i] = Astim2[i - (detected_position - Ta) * NcR];
	}

	for (i = 0; i <= (detected_position - Ta) * NcR; i++)
	{
		Astim2[i] = 0;
	}
	for (i = 0; i < pck_len * NcR; i++)
	{
		cCancelMessage[i] = Astim2[i] * data_in_post2[i] * ejthetastim2[i];
	}
}

void ZonalFloatingPoint(double complex *data_in_pre, double complex *data_out, int pck_len, float zonalValue, int nNcR, int nZonalType)
{
	int i;

	float fAngStim = 0.0;
	double dStepAngstim = PI * 2 * zonalValue / nNcR;
	double complex cFasoreStim;

	for (i = 0; i < pck_len * nNcR; i++)
	{

		fAngStim = mod(fAngStim + dStepAngstim, PI * 2);
		cFasoreStim = cos(fAngStim) - I * sin(fAngStim);

		if (nZonalType == ZONAL_ADD)
			data_out[i] = data_in_pre[i] * cFasoreStim;
		else
			data_out[i] = data_in_pre[i] * conj(cFasoreStim);
	}
	return;
}
/// AIS Receiver Sub-functions ////
/////////////////////////////////////////////////////////////////////////
void interpolation_c_(float *data_in_r, float *data_in_i, float *data_out_r,
					  float *data_out_i, int *pck_len, int *NcR, float *ist, int *istcamp)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcR) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*NcR) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcR); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	interpolation(data_in, data_out, *pck_len, *NcR, *ist, istcamp);

	for (i = 0; i < (*pck_len) * (*NcR); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

/////////////////////////////////////////////////////////////////////////
void interpolation(double complex *data_in, double complex *data_out, int pck_len, int NcR,
				   float ist, int *istcamp)
{
	float T;  // bit period
	float Tc; // sampling period
	int i;	  // support variable

	Tc = 1.0 / (float)NcR;
	*istcamp = (int)(ist / Tc);
	T = Tc + mod(ist, Tc);
	T = T / Tc;

	for (i = 2; i < pck_len * NcR; i++)
	{
		data_out[i - 1] = (T - 1.0) * (T - 2.0) * data_in[i - 2] / 2.0 - T * (T - 2.0) * data_in[i - 1] + T * (T - 1.0) * data_in[i] / 2.0;
	}

	// #CGS_2.4 rimosso il calcolo del primo e ultimo elemnto
	data_out[0] = 0.0;
	data_out[pck_len * NcR - 1] = 0.0;
	// T = T - 1.0;
	// data_out[0] = (T-1.0)*(T-2.0)*data_in[0]/2.0-T*(T-2.0)*data_in[1]
	// +T*(T-1.0)*data_in[2]/2.0;

	// T = T + 2.0;
	// data_out[pck_len*NcR-1] = (T-1.0)*(T-2.0)*data_in[pck_len*NcR-3]/2.0
	// -T*(T-2.0)*data_in[pck_len*NcR-2]+T*(T-1.0)*data_in[pck_len*NcR-1]/2.0;

	return;
}

/////////////////////////////////////////////////////////////////////////
void interpolation_q_c_(float *data_in_r, float *data_in_i, float *data_out_r,
						float *data_out_i, int *pck_len, int *NcR, float *ist, int *istcamp,
						int *p_out, int *q_out)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcR) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*NcR) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcR); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	interpolation_q(data_in, data_out, *pck_len, *NcR, *ist, istcamp, *p_out,
					*q_out);

	for (i = 0; i < (*pck_len) * (*NcR); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

/////////////////////////////////////////////////////////////////////////
void interpolation_q(double complex *data_in, double complex *data_out, int pck_len, int NcR,
					 float ist, int *istcamp, int p_out, int q_out)
{
	float T;  // bit period
	float Tc; // sampling period
			  //	float t_tmp; // support variable //#CGS_2.7 unused
	int i; // support variable
	// quantization parameters
	int p_tc, q_tc;
	int p_t, q_t;
	int p_t2, q_t2;

	p_tc = 1; // 15;
	q_tc = 15;
	p_t = 3;   // 15;
	q_t = 13;  // 15;
	p_t2 = 1;  // 15;
	q_t2 = 13; // 15;

	// FILE *interpolation;
	// interpolation = fopen("./interpolation.txt", "w");

	Tc = quantize(1.0 / (float)NcR, p_tc, q_tc);
	*istcamp = (int)(ist / Tc);
	T = Tc + mod(ist, Tc);
	T = T / Tc;
	T = quantize(T, p_t, q_t);

	float c1, c2, c3;

	c1 = (T - 1.0) * (T - 2.0) / 2.0;
	c1 = quantize(c1, p_t2, q_t2);
	c2 = -T * (T - 2.0); // controllare, tolto un diviso 2
	c2 = quantize(c2, p_t2, q_t2);
	c3 = T * (T - 1.0) / 2.0;
	c3 = quantize(c3, p_t2, q_t2);

	// double complex data_out_test[pck_len*NcR];

	for (i = 2; i < pck_len * NcR; i++) // #TF_2.3
	{
		// data_out_test[i-1] = (T-1.0)*(T-2.0)*data_in[i-2]/2.0-T*(T-2.0)*data_in[i-1]
		// +T*(T-1.0)*data_in[i]/2.0;
		data_out[i - 1] = c1 * data_in[i - 2];
		data_out[i - 1] = data_out[i - 1] - c2 * data_in[i - 1];
		data_out[i - 1] = data_out[i - 1] + c3 * data_in[i];
		data_out[i - 1] = quantize(creal(data_out[i - 1]), p_out, q_out) + I * quantize(cimag(data_out[i - 1]), p_out, q_out);
	}

	////#CGS_2.4 rimosso il calcolo del primo e ultimo elemento
	data_out[0] = 0.0;
	data_out[pck_len * NcR - 1] = 0.0;
	// T = T - 1.0;
	// data_out[0] = (T-1.0)*(T-2.0)*data_in[0]/2.0-T*(T-2.0)*data_in[1]
	// +T*(T-1.0)*data_in[2]/2.0;

	// T = T + 2.0;
	// data_out[pck_len*NcR-1] = (T-1.0)*(T-2.0)*data_in[pck_len*NcR-3]/2.0
	// -T*(T-2.0)*data_in[pck_len*NcR-2]+T*(T-1.0)*data_in[pck_len*NcR-1]/2.0;

	// for (i = 0; i < pck_len * NcR; i++)
	// {
	// data_out[i] = quantize(creal(data_out[i]), p_out, q_out) + I*quantize(cimag(data_out[i]), p_out, q_out);
	// }

	return;
}

/////////////////////////////////////////////////////////////////////////
void interpolation2_c_(float *data_in_r, float *data_in_i, float *data_out_r,
					   float *data_out_i, int *pck_len, int *NcR, float *ist)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc(((*NcR) * (*pck_len) + 1) * sizeof(double complex));
	data_out = malloc(((*NcR) * (*pck_len) + 1) * sizeof(double complex));

	for (i = 0; i <= (*pck_len + 2) * (*NcR); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	interpolation2(data_in, data_out, *pck_len + 2, *NcR, *ist);

	for (i = 0; i <= (*pck_len + 2) * (*NcR); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void interpolation2(double complex *data_in, double complex *data_out, int pck_len, int NcR,
					float ist)
{
	float Tc; // sampling period
	float T;  // bit period
	int i;	  // support variable

	Tc = 1.0 / (float)NcR;
	T = Tc + mod(ist, Tc);
	T = T / Tc;

	for (i = 1; i < pck_len * NcR; i++)
	{
		data_out[i] = (T - 1.0) * (T - 2.0) * data_in[i - 1] / 2.0 - T * (T - 2.0) * data_in[i] + T * (T - 1.0) * data_in[i + 1] / 2.0;
	}

	//	T=T - 1.0;
	//	data_out[0] = (T-1.0)*(T-2.0)*data_in[0]/2.0-T*(T-2.0)*data_in[1]
	//		+T*(T-1.0)*data_in[2]/2.0;

	//	T = T + 2.0;
	//	data_out[pck_len*NcR]=(T-1.0)*(T-2.0)*data_in[pck_len*NcR-2]/2.0
	//		-T*(T-2.0)*data_in[pck_len*NcR+NcR-1]
	//		+T*(T-1.0)*data_in[pck_len*NcR]/2.0;

	return;
}

/////////////////////////////////////////////////////////////////////////
void interpolation2_q_c_(float *data_in_r, float *data_in_i, float *data_out_r,
						 float *data_out_i, int *pck_len, int *NcR, float *ist, int *p_out,
						 int *q_out)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc(((*NcR) * (*pck_len + 2) + 1) * sizeof(double complex));
	data_out = malloc(((*NcR) * (*pck_len + 2) + 1) * sizeof(double complex));

	for (i = 0; i <= (*pck_len + 2) * (*NcR); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	interpolation2_q(data_in, data_out, *pck_len + 2, *NcR, *ist, *p_out,
					 *q_out);

	for (i = 0; i <= (*pck_len + 2) * (*NcR); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void interpolation2_q(double complex *data_in, double complex *data_out, int pck_len,
					  int NcR, float ist, int p_out, int q_out)
{
	float Tc; // sampling period
	float T;  // bit period
	int i;	  // support variable
	// quantization parameters
	int p_tc, q_tc;
	int p_t, q_t;

	p_tc = 1; // 15;
	q_tc = 15;
	p_t = 3;  // 15;
	q_t = 13; // 15;

	Tc = quantize(1.0 / (float)NcR, p_tc, q_tc);
	T = Tc + mod(ist, Tc);
	T = T / Tc;
	T = quantize(T, p_t, q_t);

	for (i = 1; i < pck_len * NcR; i++)
	{
		data_out[i] = (T - 1.0) * (T - 2.0) * data_in[i - 1] / 2.0 - T * (T - 2.0) * data_in[i] + T * (T - 1.0) * data_in[i + 1] / 2.0;
	}

	//	T=T - 1.0;
	//	data_out[0] = (T-1.0)*(T-2.0)*data_in[0]/2.0-T*(T-2.0)*data_in[1]
	//		+T*(T-1.0)*data_in[2]/2.0;

	//	T = T + 2.0;
	//	data_out[pck_len*NcR]=(T-1.0)*(T-2.0)*data_in[pck_len*NcR-2]/2.0
	//		-T*(T-2.0)*data_in[pck_len*NcR-1]
	//		+T*(T-1.0)*data_in[pck_len*NcR]/2.0;

	for (i = 0; i < pck_len * NcR + 1; i++)
	{
		data_out[i] = quantize(creal(data_out[i]), p_out, q_out) + I * quantize(cimag(data_out[i]), p_out, q_out);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void mod_CPM(int state, int symb, int M, int num_h, int p, int L, int n,
			 float *pulse, int NcT, double complex *s1, double complex *s2, double complex *s)
{
	int *a;		 // alpha parameter, see [RD1] eq. 1.8
	int phi;	 //
	float phi2;	 // phi parameter, see [RD1] eq. 1.9
	int corr;	 // alpha (signed) parameter, see [RD1] eq. 1.8
	float phase; // signal's phase
	// support variables
	int i;
	int k;
	int j;

	phi = state / pow(M, L - 1);
	corr = mod(state, pow(M, L - 1));

	a = malloc(L * sizeof(int));
	a[0] = 2 * symb - (M - 1);
	for (i = 1; i < L; i++)
	{
		a[i] = 2 * mod(corr, M) - (M - 1);
		corr = corr / M;
	}

	phi2 = 2 * PI * ((float)phi - (float)((M - 1) * (n - L + 1)) / 2.0) * (float)num_h / (float)p;

	*s2 = cos(phi2) + I * sin(phi2);

	for (k = 0; k < NcT; k++)
	{
		phase = 0.0;
		for (j = 0; j < L; j++)
		{
			phase = phase + (float)a[j] * pulse[k + j * NcT];
		}
		phase = phase * 2 * PI * (float)num_h / (float)p;

		s1[k] = cos(phase) + I * sin(phase);

		s[k] = s1[k] * (*s2);
	}

	free(a);

	return;
}

///////////////////////////////////////////////////////////////////////
void mod_cpm_q_c_(int *state, int *symb, int *M, int *num_h, int *p, int *L,
				  int *n, float *pulse, int *NcT, float *s1_r, float *s1_i, float *s2_r,
				  float *s2_i, float *s_r, float *s_i, int *p_out, int *q_out,
				  int *p_pulse, int *q_pulse)
{
	double complex *s1;
	double complex s2;
	double complex *s;
	int i;

	s1 = malloc((*NcT) * sizeof(double complex));
	s = malloc((*NcT) * sizeof(double complex));

	mod_CPM_q(*state, *symb, *M, *num_h, *p, *L, *n, pulse, *NcT, s1, &s2, s,
			  *p_out, *q_out, *p_pulse, *q_pulse);

	for (i = 0; i < (*NcT); i++)
	{
		s1_r[i] = creal(s1[i]);
		s1_i[i] = cimag(s1[i]);
		s_r[i] = creal(s[i]);
		s_i[i] = cimag(s[i]);
	}
	*s2_r = creal(s2);
	*s2_i = cimag(s2);

	free(s1);
	free(s);

	return;
}

///////////////////////////////////////////////////////////////////////
void mod_CPM_q(int state, int symb, int M, int num_h, int p, int L, int n,
			   float *pulse, int NcT, double complex *s1, double complex *s2, double complex *s, int p_out,
			   int q_out, int p_pulse, int q_pulse)
{
	int *a;
	int phi;
	float phi2;
	int corr;
	float phase;
	float pulse_int;
	int i;
	int k;
	int j;
	int p_phi;
	int q_phi;
	int p_phase;
	int q_phase;

	p_phi = 4;
	q_phi = 10;
	p_phase = 4;
	q_phase = 10;

	phi = state / pow(M, L - 1);
	corr = mod(state, pow(M, L - 1));

	a = malloc(L * sizeof(int));
	a[0] = 2 * symb - (M - 1);
	for (i = 1; i < L; i++)
	{
		a[i] = 2 * mod(corr, M) - (M - 1);
		corr = corr / M;
	}

	phi2 = 2 * PI * ((float)phi - (float)((M - 1) * (n - L + 1)) / 2.0) * (float)num_h / (float)p;
	phi2 = mod(phi2, 2 * PI);
	if (phi2 >= PI)
	{
		phi2 = phi2 - 2 * PI;
	}
	if (phi2 < -PI)
	{
		phi2 = phi2 + 2 * PI;
	}
	phi2 = quantize(phi2, p_phi, q_phi);

	*s2 = cordic_polartorect(1.0, phi2, CORDIC_ITER);

	for (k = 0; k < NcT; k++)
	{
		phase = 0.0;
		for (j = 0; j < L; j++)
		{
			pulse_int = quantize(pulse[k + j * NcT], p_pulse, q_pulse);
			phase = phase + (float)a[j] * pulse_int;
		}
		phase = phase * 2 * PI * (float)num_h / (float)p;
		phase = mod(phase, 2 * PI);
		if (phase >= PI)
		{
			phase = phase - 2 * PI;
		}
		if (phase < -PI)
		{
			phase = phase + 2 * PI;
		}
		phase = quantize(phase, p_phase, q_phase);
		s1[k] = cordic_polartorect(1.0, phase, CORDIC_ITER);
		s[k] = s1[k] * (*s2);
		s[k] = quantize(creal(s[k]), p_out, q_out) + I * quantize(cimag(s[k]),
																  p_out, q_out);
	}

	free(a);

	return;
}

///////////////////////////////////////////////////////////////////////
void next_state_c_(int *state, int *new_state, int *symb, int *M, int *p,
				   int *L)
{
	next_state(*state, new_state, *symb, *M, *p, *L);

	return;
}

///////////////////////////////////////////////////////////////////////
void next_state(int state, int *new_state, int symb, int M, int p, int L)
{
	int phi;
	int phi_new;
	int corr;
	int corr_new;
	int symb_old;

	phi = state / pow(M, L - 1);
	corr = mod(state, pow(M, L - 1));
	symb_old = corr / pow(M, L - 2);
	phi_new = mod(phi + p + symb_old, p);
	corr_new = mod(corr, pow(M, L - 2)) * M + symb;
	*new_state = phi_new * pow(M, L - 1) + corr_new;

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_c_(float *data_in_r, float *data_in_i, float *data_out_r,
			   float *data_out_i, float *H, int *lung, int *pck_len, int *del,
			   int *NcT, int *istcamp)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	filter(data_in, data_out, H, *lung, *pck_len, *del, *NcT, *istcamp);

	for (i = 0; i < (*pck_len); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void filter(double complex *data_in, double complex *data_out, float *H, int lung,
			int pck_len, int del, int NcT, int istcamp)
{
	int limi; // integration lower limit
	int lims; // integration upper limit
	int i, l;

	for (i = 1; i <= pck_len; i++)
	{
		data_out[i - 1] = 0;
		// calculate correlation between input data and filter impulse responce
		lims = MinInt((i + del) * NcT + istcamp, lung * NcT - 1);
		limi = MaxInt(0, (i + del) * NcT + istcamp - pck_len * NcT + 1);
		for (l = limi; l <= lims; l++)
		{
			// since this is a downsampling filter, 1 sample out of NcT is considered
			data_out[i - 1] = data_out[i - 1] + H[l] * data_in[(i + del) * NcT + istcamp - l];
		}
		data_out[i - 1] /= (float)NcT;
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_q_c_(float *data_in_r, float *data_in_i, float *data_out_r,
				 float *data_out_i, float *H, int *lung, int *pck_len, int *del,
				 int *NcT, int *istcamp, int *p_h, int *q_h, int *p_in, int *q_in,
				 int *p_out, int *q_out)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	filter_q(data_in, data_out, H, *lung, *pck_len, *del, *NcT, *istcamp, *p_h,
			 *q_h, *p_in, *q_in, *p_out, *q_out);

	for (i = 0; i < (*pck_len); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_q(double complex *data_in, double complex *data_out, float *H, int lung,
			  int pck_len, int del, int NcT, int istcamp, int p_h, int q_h, int p_in,
			  int q_in, int p_out, int q_out)
{
	float H_int;
	int limi;
	int lims;
	int p_int, q_int;
	int i, l;

	p_int = 5;
	q_int = 11;

	for (i = 1; i <= pck_len; i++)
	{
		data_out[i - 1] = 0;
		lims = MinInt((i + del) * NcT + istcamp, lung * NcT - 1);
		limi = MaxInt(0, (i + del) * NcT + istcamp - pck_len * NcT + 1);
		for (l = limi; l <= lims; l++)
		{
			data_in[(i + del) * NcT + istcamp - l] = quantize(creal(data_in[(i + del) * NcT + istcamp - l]), p_in, q_in) + I * quantize(
																																   cimag(data_in[(i + del) * NcT + istcamp - l]), p_in, q_in);
			H_int = quantize(creal(H[l]), p_h, q_h);
			data_out[i - 1] = data_out[i - 1] + H_int * data_in[(i + del) * NcT + istcamp - l];
			data_out[i - 1] = quantize(creal(data_out[i - 1]), p_int, q_int) + I * quantize(cimag(data_out[i - 1]), p_int, q_int);
		}
		data_out[i - 1] = data_out[i - 1] / (float)NcT;
		data_out[i - 1] = quantize(creal(data_out[i - 1]), p_out, q_out) + I * quantize(cimag(data_out[i - 1]), p_out, q_out);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_nd_c_(float *data_in_r, float *data_in_i, float *data_out_r,
				  float *data_out_i, float *H, int *lung, int *pck_len, int *del,
				  int *NcT)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	filter_nd(data_in, data_out, H, *lung, *pck_len, *del, *NcT);

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_nd(double complex *data_in, double complex *data_out, float *H, int lung,
			   int pck_len, int del, int NcT)
{
	int limi; // integration lower limit
	int lims; // integration upper limit
	int Nimp;
	int i, l;

	Nimp = lung * NcT + 1;

	for (i = 1; i <= pck_len * NcT; i++)
	{
		data_out[i - 1] = 0;
		// calculate correlation between input data and filter impulse responce
		lims = MinInt(Nimp, i + del - 1);
		limi = MaxInt(1, i + del - pck_len * NcT);
		for (l = limi; l <= lims; l++)
		{
			data_out[i - 1] = data_out[i - 1] + H[l - 1] * data_in[i + del - l - 1];
		}
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_nd_q_c_(float *data_in_r, float *data_in_i, float *data_out_r,
					float *data_out_i, float *H, int *lung, int *pck_len, int *del,
					int *NcT, int *p_h, int *q_h, int *p_in, int *q_in, int *p_out,
					int *q_out)
{
	int i;
	double complex *data_in;
	double complex *data_out;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_out = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	filter_nd_q(data_in, data_out, H, *lung, *pck_len, *del, *NcT, *p_h, *q_h,
				*p_in, *q_in, *p_out, *q_out);

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_out_r[i] = creal(data_out[i]);
		data_out_i[i] = cimag(data_out[i]);
	}

	free(data_in);
	free(data_out);

	return;
}

///////////////////////////////////////////////////////////////////////
void filter_nd_q(double complex *data_in, double complex *data_out, float *H, int lung,
				 int pck_len, int del, int NcT, int p_h, int q_h, int p_in, int q_in,
				 int p_out, int q_out)
{
	float H_int;
	int limi;
	int lims;
	int Nimp;
	int i, l;

	Nimp = lung * NcT + 1;

	for (i = 1; i <= pck_len * NcT; i++)
	{
		data_out[i - 1] = 0;
		lims = MinInt(Nimp, i + del - 1);
		limi = MaxInt(1, i + del - pck_len * NcT);
		for (l = limi; l <= lims; l++)
		{
			data_in[i + del - l - 1] = quantize(
										   creal(data_in[i + del - l - 1]), p_in, q_in) +
									   I * quantize(cimag(data_in[i + del - l - 1]), p_in, q_in);
			H_int = quantize(creal(H[l - 1]), p_h, q_h);
			data_out[i - 1] = data_out[i - 1] + H_int * data_in[i + del - l - 1];
			data_out[i - 1] = quantize(creal(data_out[i - 1]), p_out, q_out) + I * quantize(cimag(data_out[i - 1]), p_out, q_out);
		}
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq_c_(float *data_in_r, float *data_in_i, int *NcT, int *pck_len,
					int *offset, int *M, int *L, float *ist_re, int *ist_int, float *fstim)
{
	int i;
	double complex *data_in;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	nda_timfreq(data_in, *NcT, *pck_len, *offset, *M, *L, ist_re, ist_int,
				fstim);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq(double complex *data_in, int NcT, int pck_len, int offset, int M,
				 int L, float *ist_re, int *ist_int, float *fstim)
{
	float *A;			   // A1 parameter, see [RD1] eq. 3.3
	double complex *z;	   // z parameter, see [RD1] eq 3.2
	double complex **data; // Rm vector, ssee [RD1] eq 3.2
	// support variables
	double complex xx, xx1;
	double complex ej;
	float ar;
	int i, m, n;

	// initialize A coefficient (see technical note)
	A = malloc(M * sizeof(float));
	A[0] = 1;
	for (i = 1; i < M; i++)
	{
		A[i] = 1.25; // 1.24;
	}

	// initialize z vector
	z = malloc(L * NcT * sizeof(double complex));
	for (i = 0; i < L * NcT; i++)
	{
		z[i] = data_in[i + (offset - 1) * NcT];
	}

	data = malloc(M * sizeof(double complex *));
	for (i = 0; i < M; i++)
	{
		data[i] = malloc(NcT * sizeof(double complex));
	}

	// calculate Rm vector (see technical note)
	for (i = 0; i < NcT; i++)
	{
		for (m = 1; m <= M; m++)
		{
			data[m - 1][i] = 0;
			for (n = m; n < L; n++)
			{
				data[m - 1][i] = data[m - 1][i] + cpow(z[n * NcT + i] * conj(
																			z[(n - m) * NcT + i]),
													   2);
			}
			data[m - 1][i] = data[m - 1][i] / (float)(L - m);
		}
	}

	// esteem timing (see technical note)
	xx = 0;
	for (i = 0; i < NcT; i++)
	{
		ar = 2 * PI * (float)i / (float)NcT;
		ej = cos(ar) - I * sin(ar);
		xx1 = 0;
		for (m = 0; m < M; m++)
		{
			xx1 = xx1 + cabs(data[m][i]) * A[m];
		}

		xx = xx + xx1 * ej;
	}
	*ist_re = mod(-atan2(cimag(xx), creal(xx)) + 2 * PI, 2 * PI) / (2 * PI);
	*ist_int = mod(rintf(*ist_re * (float)NcT), NcT);

	// esteem frequency (see technical note)
	*fstim = 0;
	for (m = 1; m < M; m++)
	{
		xx = -data[m][*ist_int] * conj(data[m - 1][*ist_int]);
		*fstim = *fstim + atan2(cimag(xx), creal(xx));
	}
	*fstim = *fstim / (2 * 2 * PI * (float)M);

	free(A);
	free(z);
	for (i = 0; i < M; i++)
		free(data[i]);
	free(data);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq2T_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *M, int *L, float *ist_re, int *ist_int,
					  float *fstim)
{
	int i;
	double complex *data_in;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	nda_timfreq2T(data_in, *NcT, *pck_len, *offset, *M, *L, ist_re, ist_int,
				  fstim);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq2T(double complex *data_in, int NcT, int pck_len, int offset, int M,
				   int L, float *ist_re, int *ist_int, float *fstim)
{
	float *A;			   // A1 parameter, see RD1 eq. 3.3
	double complex *z;	   // z parameter, see RD1 eq 3.2
	double complex **data; // Rm vector, see RD1 eq 3.2
	// support variables
	double complex xx, xx1;
	double complex ej;
	float ar;
	int i, m, n;

	// initialize A coefficient (see tech\nical note)
	A = malloc(M * sizeof(float));
	A[0] = 1;
	for (i = 1; i < M; i++)
	{
		A[i] = 1.25; // 1.24;
	}

	// initialize z vector
	z = malloc(L * NcT * sizeof(double complex));
	for (i = 0; i < L * NcT; i++)
	{
		z[i] = data_in[i + (offset - 1) * NcT];
	}

	data = malloc(M * sizeof(double complex *));
	for (i = 0; i < M; i++)
	{
		data[i] = malloc(NcT * sizeof(double complex));
	}

	// calculate Rm vector (see technical note)
	for (i = 0; i < NcT; i++)
	{
		for (m = 2; m <= M; m = m + 2)
		{
			data[m - 1][i] = 0;
			for (n = m; n < L; n++)
			{
				data[m - 1][i] = data[m - 1][i] + cpow(z[n * NcT + i] * conj(
																			z[(n - m) * NcT + i]),
													   2);
			}
			data[m - 1][i] = data[m - 1][i] / (float)(L - m);
		}
	}

	// esteem timing (see technical note)
	xx = 0;
	for (i = 0; i < NcT; i++)
	{
		ar = 2 * PI * (float)i / (float)NcT;
		ej = cos(ar) - I * sin(ar);
		xx1 = 0;
		for (m = 1; m < M; m = m + 2)
		{
			xx1 = xx1 + cabs(data[m][i]) * A[m];
		}

		xx = xx + xx1 * ej;
	}
	*ist_re = mod(-atan2(cimag(xx), creal(xx)) + 2 * PI, 2 * PI) / (2 * PI);
	*ist_int = mod(rintf(*ist_re * (float)NcT), NcT);

	// esteem frequency (see technical note)
	*fstim = 0;
	for (m = 3; m < M; m = m + 2)
	{
		xx = data[m][*ist_int] * conj(data[m - 2][*ist_int]);
		*fstim = *fstim + atan2(cimag(xx), creal(xx));
	}
	*fstim = *fstim / (2 * 2 * PI * (float)M);

	free(A);
	free(z);
	for (i = 0; i < M; i++)
		free(data[i]);
	free(data);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreqMCM_c_(float *data_in_r, float *data_in_i, int *NcT,
					   int *pck_len, int *offset, int *L, int *ist_int, float *fstim)
{
	int i;
	double complex *data_in;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	nda_timfreqMCM(data_in, *NcT, *pck_len, *offset, *L, ist_int, fstim);

	free(data_in);

	return;
}

// #CGS_2.7 commentate variabili non utilizzate
///////////////////////////////////////////////////////////////////////
void nda_timfreqMCM(double complex *data_in, int NcT, int pck_len, int offset, int L,
					int *ist_int, float *fstim)
{
	double complex *x2; // z parameter, see RD1 eq 3.2
						//	double complex **data; // Rm vector, see RD1 eq 3.2
	// support variables
	//	double complex xx, xx1;
	//	double complex ej, facc;
	double complex facc;
	//	float ar;
	//  int i, m, n;
	int i;

	// initialize z vector
	x2 = malloc(L * sizeof(double complex));
	facc = 0.f + I * 0.f;
	for (i = 0; i < L; i++)
	{
		x2[i] = cpow(data_in[*ist_int + (i + 1) * NcT + (offset - 1) * NcT], 2);
		if (i > 1)
			facc += x2[i] * conj(x2[i - 1]);
	}

	*fstim = (mod(atan2(cimag(facc), creal(facc)) + 2 * PI, 2 * PI) - PI) / (2 *
																			 2 * PI);

	free(x2);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreqMM_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *L, int *ist_int, float *fstim)
{
	int i;
	double complex *data_in;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	nda_timfreqMM(data_in, *NcT, *pck_len, *offset, *L, ist_int, fstim);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreqMM(double complex *data_in, int NcT, int pck_len, int offset, int L,
				   int *ist_int, float *fstim)
{
	double complex *z2;
	double complex *data;
	float *w;
	// support variables
	double complex cc;
	float mmod;
	int i, m, NN;

	NN = L / 2;
	// initialize z vector
	z2 = malloc(L * sizeof(double complex));
	for (i = 0; i < L; i++)
	{
		z2[i] = pow((-1), i) * cpow(data_in[*ist_int + i * NcT + (offset - 1) * NcT], 2);
	}

	data = malloc((NN + 1) * sizeof(double complex));

	// calculate Rm vector (see technical note)
	data[0] = 1.0 + I * 0.0;
	for (i = 1; i < NN + 1; i++)
	{
		data[i] = 0.0 + I * 0.0;
		for (m = i; m < L; m++)
		{
			data[i] += z2[m] * conj(z2[m - i]);
		}
		data[i] = data[i] / (float)(L - i);
	}
	w = malloc((NN + 1) * sizeof(float));
	for (m = 1; m < NN + 1; m++)
	{
		w[m] = 3. * (float)((L - m) * (L - m + 1) - NN * (L - NN)) / (float)(4 * pow(NN, 2) - 6 * NN * L + 3 * pow(L, 2) - 1);
	}

	// esteem frequency
	*fstim = 0;
	for (m = 1; m < NN + 1; m++)
	{
		cc = data[m] * conj(data[m - 1]);
		mmod = mod(atan2(cimag(cc), creal(cc)) + 2 * PI, 2 * PI);
		if (mmod >= PI)
			mmod = mmod - 2 * PI;
		*fstim += w[m] * mmod;
	}
	*fstim = *fstim / (2 * 2 * PI * (float)NN);

	free(z2);
	free(data);
	free(w);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq_q_c_(float *data_in_r, float *data_in_i, int *NcT,
					  int *pck_len, int *offset, int *M, int *L, float *ist_re, int *ist_int,
					  float *fstim, int *p_in, int *q_in, int *p_f, int *q_f)
{
	int i;
	double complex *data_in;

	data_in = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	nda_timfreq_q(data_in, *NcT, *pck_len, *offset, *M, *L, ist_re, ist_int,
				  fstim, *p_in, *q_in, *p_f, *q_f);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void nda_timfreq_q(double complex *data_in, int NcT, int pck_len, int offset, int M,
				   int L, float *ist_re, int *ist_int, float *fstim, int p_in, int q_in,
				   int p_f, int q_f)
{
	float *A;
	double complex *z;
	double complex z_tmp;
	double complex **data;
	double complex xx; // double complex xx, xx1; //#CGS_2.7 unused
	double complex ej;
	float ar;
	float abs_data;
	float xx1_float;
	float phase;
	float fstim_tmp;
	int p_ztmp, q_ztmp;
	int p_data, q_data, p_data2, q_data2, p_ej, q_ej;
	int p_xx1, q_xx1, p_f1, q_f1;
	int p_xt, q_xt, p_xf, q_xf;
	int p_phase, q_phase;
	int p_ist, q_ist;
	int i, m, n;

	// define quantization parameters for internal variables
	p_ztmp = 1;
	q_ztmp = 15;
	p_data = 8;
	q_data = 15;
	p_ej = 2;
	q_ej = 15;
	p_data2 = 1;
	q_data2 = 15;
	p_xx1 = 5;
	q_xx1 = 15;
	p_xt = 8;
	q_xt = 15;
	p_xf = 8;
	q_xf = 15;
	p_f1 = 7;
	q_f1 = 15;
	p_phase = 3;
	q_phase = 15;
	p_ist = 2;
	q_ist = 15;

	A = malloc(M * sizeof(float));
	A[0] = 1.0;
	for (i = 1; i < M; i++)
	{
		A[i] = 1.25; // 1.24;
	}

	z = malloc(L * NcT * sizeof(double complex));
	for (i = 0; i < L * NcT; i++)
	{
		z[i] = data_in[i + (offset - 1) * NcT];
	}

	data = malloc(M * sizeof(double complex *));
	for (i = 0; i < M; i++)
	{
		data[i] = malloc(NcT * (sizeof(double complex)));
	}

	for (i = 0; i < NcT; i++)
	{
		for (m = 1; m <= M; m++)
		{
			data[m - 1][i] = 0;
			for (n = m; n < L; n++)
			{
				z[n * NcT + i] = quantize(creal(z[n * NcT + i]), p_in, q_in) + I * quantize(cimag(z[n * NcT + i]), p_in, q_in);
				z[(n - m) * NcT + i] = quantize(creal(z[(n - m) * NcT + i]),
												p_in, q_in) +
									   I * quantize(cimag(z[(n - m) * NcT + i]),
													p_in, q_in);
				// 							data[m-1][i] = data[m-1][i] +
				// 								cpow(z[n*NcT+i]*conj(z[(n-m)*NcT+i]), 2);
				z_tmp = z[n * NcT + i] * conj(z[(n - m) * NcT + i]);
				z_tmp = quantize(creal(z_tmp), p_ztmp, q_ztmp) + I * quantize(
																		 cimag(z_tmp), p_ztmp, q_ztmp);
				z_tmp = z_tmp * z_tmp;
				z_tmp = quantize(creal(z_tmp), p_ztmp, q_ztmp) + I * quantize(
																		 cimag(z_tmp), p_ztmp, q_ztmp);
				data[m - 1][i] = data[m - 1][i] + z_tmp;
				data[m - 1][i] = quantize(creal(data[m - 1][i]), p_data, q_data) + I * quantize(cimag(data[m - 1][i]), p_data,
																								q_data);
			}
			data[m - 1][i] = data[m - 1][i] / (float)(L - m);
			data[m - 1][i] = quantize(creal(data[m - 1][i]), p_data2, q_data2) + I * quantize(cimag(data[m - 1][i]), p_data2, q_data2);
		}
	}

	xx = 0;
	for (i = 0; i < NcT; i++)
	{
		ar = 2 * PI * (float)i / (float)NcT;
		if (ar >= PI)
		{
			ar = ar - PI * 2;
		}
		else if (ar <= -PI)
		{
			ar = ar + PI * 2;
		}

		ej = cordic_polartorect(1.0, ar, CORDIC_ITER);
		ej = quantize(creal(ej), p_ej, q_ej) + I * quantize(cimag(ej), p_ej,
															q_ej);
		ej = conj(ej);
		xx1_float = 0;
		for (m = 0; m < M; m++)
		{
			abs_data = quantize(cordic_getmag(creal(data[m][i]), cimag(data[m][i]), CORDIC_ITER), p_data2, q_data2);
			xx1_float = xx1_float + abs_data * A[m];
			xx1_float = quantize(xx1_float, p_xx1, q_xx1);
		}

		//		xx1 = xx1_float + I * 0.0; 	//#CGS_2.7 unused

		xx = xx + xx1_float * ej;
		xx = quantize(creal(xx), p_xt, q_xt) + I * quantize(cimag(xx), p_xt,
															q_xt);
	}
	phase = cordic_getphase(creal(xx), cimag(xx), CORDIC_ITER);
	phase = quantize(phase, p_phase, q_phase);
	*ist_re = mod(2 * PI - phase, 2 * PI) / (2 * PI);

	*ist_re = quantize(*ist_re, p_ist, q_ist);

	*ist_int = mod(rintf(*ist_re * (float)NcT), NcT);

	fstim_tmp = 0;
	for (m = 1; m < M; m++)
	{
		xx = -data[m][*ist_int] * conj(data[m - 1][*ist_int]);
		xx = quantize(creal(xx), p_xf, q_xf) + I * quantize(cimag(xx), p_xf,
															q_xf);
		phase = cordic_getphase(creal(xx), cimag(xx), CORDIC_ITER);
		phase = quantize(phase, p_phase, q_phase);
		fstim_tmp = fstim_tmp + phase;
		fstim_tmp = quantize(fstim_tmp, p_f1, q_f1);
	}

	*fstim = quantize(fstim_tmp / (2 * 2 * PI * (float)M), p_f, q_f);

	free(A);
	free(z);
	for (i = 0; i < M; i++)
		free(data[i]); // #TF_2.3
	free(data);

	return;
}

///////////////////////////////////////////////////////////////////////
void bcjr_c_(float *data_in_r, float *data_in_i, float *LLR, int *AD,
			 int *pck_len, float *Pd, int *D)
{
	int i;
	double complex *data_in;

	data_in = malloc((*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	BCJR(data_in, LLR, AD, *pck_len, *Pd, *D);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void BCJR(double complex *data_in, float *LLR, int *AD, int pck_len, float Pd, int D)
{
	float N0s;	  // N0 parameter in equation 3.10 [RD1]
	float APP[2]; // a posteriori probability
	float **nuf;  // forward recursion variable
	float **nub;  // backward recursion variable

	// support variables
	float ang, somr;
	float *v;
	double complex *y;
	double complex ex1;
	int n, i;

	// vectors generation

	v = malloc(D * sizeof(float));

	nuf = malloc(D * sizeof(float *));
	for (i = 0; i < D; i++)
	{
		nuf[i] = malloc(pck_len * sizeof(float));
	}

	nub = malloc(D * sizeof(float *));
	for (i = 0; i < D; i++)
	{
		nub[i] = malloc(pck_len * sizeof(float));
	}

	N0s = 1.0 / 30.0;

	y = malloc(pck_len * sizeof(double complex));
	for (n = 0; n < pck_len; n++)
	{
		ang = 2 * PI * (n + 1) / 4;
		ex1 = cos(ang) + I * sin(ang);
		y[n] = data_in[n] * ex1 / N0s;
	}

	// forward recursion (see RD, eq 3.10)
	H_2d(y[0], D, nuf, 0);
	for (n = 1; n < pck_len; n++)
	{
		fact_exp1(nuf, D, Pd, nuf, n - 1, n);
		H_1d(y[n], D, v);
		var_exp_2d(v, nuf, D, nuf, n);
		normalize(nuf, D, n);
	}

	// backward recursion
	H_2d(y[pck_len - 1], D, nub, pck_len - 1);
	for (n = pck_len - 1; n > 0; n--)
	{
		fact_exp1(nub, D, Pd, nub, n, n - 1);
		H_1d(y[n - 1], D, v);
		var_exp_2d(v, nub, D, nub, n - 1);
		normalize(nub, D, n - 1);
	}

	// completion (see RD, eq 3.11)
	for (n = 1; n < pck_len; n++)
	{
		for (i = 0; i <= 1; i++)
		{
			fact_exp2(nuf, D, Pd, i, v, n - 1);
			var_exp_1d(v, nub, D, v, n);
			somr = integr(v, D);
			APP[i] = somr;
		}
		// bit estimation
		LLR[n] = APP[0] - APP[1];
		AD[n] = 0;
		if (LLR[n] < 0)
		{
			AD[n] = 1;
		}
	}

	free(v);
	for (i = 0; i < D; i++)
	{
		free(nuf[i]);
		free(nub[i]);
	}
	free(nuf);
	free(nub);
	free(y);

	return;
}

///////////////////////////////////////////////////////////////////////
void bcjr_q_c_(float *data_in_r, float *data_in_i, float *LLR, int *AD,
			   int *pck_len, float *Pd, int *D, int *p_in, int *q_in, int *p_llr,
			   int *q_llr)
{
	int i;
	double complex *data_in;

	data_in = malloc((*pck_len) * sizeof(double complex));

	for (i = 0; i < (*pck_len); i++)
	{
		data_in[i] = data_in_r[i] + I * data_in_i[i];
	}

	BCJR_q(data_in, LLR, AD, *pck_len, *Pd, *D, *p_in, *q_in, *p_llr, *q_llr);

	free(data_in);

	return;
}

///////////////////////////////////////////////////////////////////////
void BCJR_q(double complex *data_in, float *LLR, int *AD, int pck_len, float Pd,
			int D, int p_in, int q_in, int p_llr, int q_llr)
{
	float N0s;
	//	float ang, somr;
	float somr;
	float APP[2];
	float *v;
	float **nuf;
	float **nub;
	double complex *y;
	double complex ex1;
	float cos_ang, sin_ang;
	int n_mod_4;
	int p_y, q_y, p_nuf, q_nuf, p_nub, q_nub, p_v, q_v;
	int n, i;

	p_y = 8;
	q_y = 15;
	p_nuf = 8;
	q_nuf = 15;
	p_nub = 8;
	q_nub = 15;
	p_v = 8;
	q_v = 15;

	v = malloc(D * sizeof(float));

	nuf = malloc(D * sizeof(float *));
	for (i = 0; i < D; i++)
	{
		nuf[i] = malloc(pck_len * sizeof(float));
	}

	nub = malloc(D * sizeof(float *));
	for (i = 0; i < D; i++)
	{
		nub[i] = malloc(pck_len * sizeof(float));
	}

	N0s = 1.0 / 30.0;

	y = malloc(pck_len * sizeof(double complex));
	for (n = 0; n < pck_len; n++)
	{
		n_mod_4 = mod(n, 4);
		if (n_mod_4 == 0)
		{
			cos_ang = 0.0;
			sin_ang = 1.0;
		}
		else if (n_mod_4 == 1)
		{
			cos_ang = -1.0;
			sin_ang = 0.0;
		}
		else if (n_mod_4 == 2)
		{
			cos_ang = 0.0;
			sin_ang = -1.0;
		}
		else if (n_mod_4 == 3)
		{
			cos_ang = 1.0;
			sin_ang = 0.0;
		}
		else
		{
			cos_ang = 0.0;
			sin_ang = 0.0;
		}
		ex1 = cos_ang + I * sin_ang;
		data_in[n] = quantize(creal(data_in[n]), p_in, q_in) + I * quantize(
																	   cimag(data_in[n]), p_in, q_in);
		y[n] = data_in[n] * ex1 / N0s;
		y[n] = quantize(creal(y[n]), p_y, q_y) + I * quantize(cimag(y[n]), p_y,
															  q_y);
	}

	// forward recursion
	H_q_2d(y[0], D, nuf, 0, p_nuf, q_nuf);

	for (n = 1; n < pck_len; n++)
	{

		fact_exp1_q(nuf, D, Pd, nuf, n - 1, n, p_nuf, q_nuf);
		H_q_1d(y[n], D, v, p_v, q_v);
		var_exp_2d_q(v, nuf, D, nuf, n, p_nuf, q_nuf);
		normalize_q(nuf, D, n, p_nuf, q_nuf);
	}

	// backward recursion
	H_q_2d(y[pck_len - 1], D, nub, pck_len - 1, p_nub, q_nub);

	for (n = pck_len - 1; n > 0; n--)
	{

		fact_exp1_q(nub, D, Pd, nub, n, n - 1, p_nub, q_nub);
		H_q_1d(y[n - 1], D, v, p_v, q_v);
		var_exp_2d_q(v, nub, D, nub, n - 1, p_nub, q_nub);
		normalize_q(nub, D, n - 1, p_nub, q_nub);
	}

	// completion
	for (n = 1; n < pck_len; n++)
	{

		for (i = 0; i <= 1; i++)
		{
			fact_exp2_q(nuf, D, Pd, i, v, n - 1, p_v, q_v);
			var_exp_1d_q(v, nub, D, v, n, p_v, q_v);
			somr = quantize(integr(v, D), p_v, q_v);
			APP[i] = somr;
		}

		LLR[n] = quantize(APP[0] - APP[1], p_llr, q_llr);

		AD[n] = 0;
		if (LLR[n] < 0)
		{
			AD[n] = 1;
		}
	}

	free(v);
	for (i = 0; i < D; i++)
	{
		free(nuf[i]);
		free(nub[i]);
	} // #TF_2.3
	free(nuf);
	free(nub);
	free(y);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafreqmm_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					float *data_in2_i, int *offset, float *freq, int *NcT, int *pck_len)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;
	//	float f; //#CGS_2.7

	data_in1 = malloc((*pck_len) * (*NcT) * sizeof(double complex));
	data_in2 = malloc((*pck_len) * (*NcT) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimafreqMM(data_in1, data_in2, freq, *offset, *NcT, *pck_len);

	free(data_in1);
	free(data_in2);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafreqMM(double complex *data_in1, double complex *data_in2, float *freq, int offset,
				 int NcT, int pck_len)
{
	int L;			   // length of the vector employed for frequency estimation
	int NN;			   // design parameter, taken as L/2
	float *w;		   // vector of weights (see Mengali-Morelli algorithm)
	float T;		   // bit period
	double complex *z; // contains the received samples properly derotated by the known useful signal
	double complex *R; // autocorrelation of vector z
	double complex cc; // support variable
	float mmod;		   // support variable
	int m, k;		   // support variables

	// vectors initialization
	w = malloc(NcT * pck_len / 2 * sizeof(float));
	z = malloc(NcT * pck_len * sizeof(double complex));
	R = malloc((NcT * pck_len / 2 + 1) * sizeof(double complex));

	// constants initialization
	L = NcT * pck_len;
	NN = L / 2;
	T = 1.0 / (float)NcT;

	for (k = 0; k < L; k++)
	{
		z[k] = data_in1[k + offset * NcT] * conj(data_in2[k + offset * NcT]);
	}

	// R calculation (see Morelli Mengali, eq. 30)
	R[0] = 1;
	for (m = 1; m <= NN; m++)
	{
		R[m] = 0;
		for (k = m; k < L; k++)
		{
			R[m] = R[m] + z[k] * conj(z[k - m]);
		}
		R[m] = R[m] / (float)(L - m);
	}

	for (m = 1; m <= NN; m++)
	{
		w[m - 1] = 3.0 * (float)((L - m) * (L - m + 1) - NN * (L - NN)) / (float)((4 * pow(NN, 2) - 6 * NN * L + 3 * pow(L, 2) - 1));
	}

	*freq = 0;
	for (m = 1; m <= NN; m++)
	{
		cc = R[m] * conj(R[m - 1]);
		mmod = mod(atan2(cimag(cc), creal(cc)) + 2 * PI, 2 * PI);
		if (mmod >= PI)
		{
			mmod = mmod - 2 * PI;
		}
		*freq = *freq + w[m - 1] * mmod;
	}

	*freq = *freq / (2 * PI * T * NN);

	free(w);
	free(z);
	free(R);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafreqmm_q_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					  float *data_in2_i, float *freq, int *offset, int *NcT, int *pck_len,
					  int *p_in1, int *q_in1, int *p_in2, int *q_in2, int *p_f, int *q_f)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;
	//	float f; //#CGS_2.7

	data_in1 = malloc((*pck_len) * (*NcT) * sizeof(double complex));
	data_in2 = malloc((*pck_len) * (*NcT) * sizeof(double complex));

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*pck_len) * (*NcT); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimafreqMM_q(data_in1, data_in2, freq, *offset, *NcT, *pck_len, *p_in1,
				  *q_in1, *p_in2, *q_in2, *p_f, *q_f);

	free(data_in1);
	free(data_in2);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafreqMM_q(double complex *data_in1, double complex *data_in2, float *freq,
				   int offset, int NcT, int pck_len, int p_in1, int q_in1, int p_in2,
				   int q_in2, int p_f, int q_f)
{
	int L, m, k, NN;
	float *w;
	float mmod, T;
	double complex *z;
	double complex *R;
	double complex cc;
	int p_z, q_z, p_r, q_r, p_cc, q_cc, p_w, q_w;
	//	int p_r2, q_r2, p_mmod, q_mmod, p_fint, q_fint, p_lm, q_lm;  //#CGS_2.7
	int p_r2, q_r2, p_fint, q_fint, p_lm, q_lm;
	float phase, temp;
	int p_phase, q_phase;

	p_z = 4;
	q_z = 4;
	p_r = 15;
	q_r = 0;
	p_cc = 8;
	q_cc = 8;
	p_w = 3;
	q_w = 13;
	p_r2 = 5;
	q_r2 = 11;
	//	p_mmod = 4;  //#CGS_2.7
	//	q_mmod = 12; //#CGS_2.7
	p_fint = 10; // 5
	q_fint = 15;
	p_phase = 1;
	q_phase = 15;
	p_lm = 1;
	q_lm = 13;

	w = malloc(NcT * pck_len / 2 * sizeof(float));
	z = malloc(NcT * pck_len * sizeof(double complex));
	R = malloc((NcT * pck_len / 2 + 1) * sizeof(double complex));

	// L = NcT * pck_len;
	L = 224;  // pck_len; //Ober mettere come define
	NN = 112; // L / 2; //Ober mettere come define
	// T = 1.0 / (float)NcT;
	T = 1.0;

	for (k = 0; k < L; k++)
	{
		data_in1[(k + offset) * NcT] = quantize(creal(data_in1[(k + offset) * NcT]), p_in1, q_in1) + I * quantize(cimag(data_in1[(k + offset) * NcT]), p_in1, q_in1);
		data_in2[(k + offset) * NcT] = quantize(creal(data_in2[(k + offset) * NcT]), p_in2, q_in2) + I * quantize(cimag(data_in2[(k + offset) * NcT]), p_in2, q_in2);
		z[k] = data_in1[(k + offset) * NcT] * conj(data_in2[(k + offset) * NcT]);
		z[k] = quantize(creal(z[k]), p_z, q_z) + I * quantize(cimag(z[k]), p_z,
															  q_z);
	}

	R[0] = 1;
	for (m = 1; m <= NN; m++)
	{
		R[m] = 0;
		temp = quantize(1. / (float)(L - m), p_lm, q_lm);
		for (k = m; k < L; k++)
		{
			R[m] = quantize(creal(R[m]), p_r, q_r) + I * quantize(cimag(R[m]),
																  p_r, q_r);
			R[m] = R[m] + z[k] * conj(z[k - m]);
			R[m] = quantize(creal(R[m]), p_r, q_r) + I * quantize(cimag(R[m]),
																  p_r, q_r);
		}
		R[m] = R[m] * temp;
		R[m] = quantize(creal(R[m]), p_r2, q_r2) + I * quantize(cimag(R[m]),
																p_r2, q_r2);
	}

	for (m = 1; m <= NN; m++)
	{
		w[m - 1] = 3.0 * (float)((L - m) * (L - m + 1) - NN * (L - NN)) / (float)((4 * pow(NN, 2) - 6 * NN * L + 3 * pow(L, 2) - 1));
		w[m - 1] = quantize(w[m - 1], p_w, q_w);
	}

	*freq = 0;
	for (m = 1; m <= NN; m++)
	{
		cc = R[m] * conj(R[m - 1]);
		cc = quantize(creal(cc), p_cc, q_cc) + I * quantize(cimag(cc), p_cc,
															q_cc);
		phase = cordic_getphase(creal(cc), cimag(cc), CORDIC_ITER);
		phase = quantize(phase, p_phase, q_phase);
		//			mmod = mod(phase + 2*PI, 2*PI);
		//			mmod = quantize(mmod, p_mmod, q_mmod);
		//			if (mmod >= PI)
		//				{
		//					mmod = mmod - 2*PI;
		//				}
		//			mmod = quantize(mmod, p_mmod, q_mmod);

		mmod = phase;

		*freq = quantize(*freq + w[m - 1] * mmod, p_fint, q_fint);
	}
	temp = quantize(1. / (T * NN * NcT), p_f, q_f);
	*freq = quantize(*freq * temp, p_f, q_f);

	free(w);
	free(z);
	free(R);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimatiming_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					float *data_in2_i, float *time, int *offset, int *NcT, int *pck_len)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;

	data_in1 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_in2 = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimatiming(data_in1, data_in2, time, *offset, *NcT, *pck_len);

	free(data_in1);
	free(data_in2);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimatiming(double complex *data_in1, double complex *data_in2, float *time, int offset,
				 int NcT, int pck_len)
{
	int L;	   // number of samples per packet length
	float num; // numerator
	float den; // denominator
	double complex cc1, cc2, cc3;
	int k; // support variables

	L = NcT * pck_len;
	cc1 = 0.0 + I * 0.0;
	cc2 = 0.0 + I * 0.0;
	cc3 = 0.0 + I * 0.0;

	for (k = 1; k <= L; k++)
	{
		cc1 = cc1 + data_in1[k - 1 + offset * NcT] * conj(data_in2[offset * NcT + k - 1]);
		cc2 = cc2 + data_in1[k - 1 + offset * NcT] * conj(data_in2[offset * NcT + k - 2]);
		cc3 = cc3 + data_in1[k - 1 + offset * NcT] * conj(data_in2[offset * NcT + k]);
	}

	num = creal(cc1 * conj(-cc2 + cc3)) / 2.0;

	den = pow(cabs(-cc2 / 2.0 + cc3 / 2.0), 2) + creal(cc1 * conj(cc2 + cc3 - 2.0 * cc1));

	*time = -num / (den * (float)NcT);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimatiming_q_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					  float *data_in2_i, float *time, int *offset, int *NcT, int *pck_len,
					  int *p_in1, int *q_in1, int *p_in2, int *q_in2, int *p_t, int *q_t)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;

	data_in1 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_in2 = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimatiming_q(data_in1, data_in2, time, *offset, *NcT, *pck_len, *p_in1,
				  *q_in1, *p_in2, *q_in2, *p_t, *q_t);

	free(data_in1);
	free(data_in2);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimatiming_q(double complex *data_in1, double complex *data_in2, float *time,
				   int offset, int NcT, int pck_len, int p_in1, int q_in1, int p_in2,
				   int q_in2, int p_t, int q_t)
{
	//	int m; //#CGS_2.7
	int k;
	int L;
	float num;
	float den;
	double complex cc1, cc2, cc3;
	double complex cc1_tmp, cc2_tmp, cc3_tmp;
	int p_cc;
	int q_cc;
	int p_n;
	int q_n;
	int p_d;
	int q_d;

	p_cc = 13;
	q_cc = 11;
	p_n = 21;
	q_n = 0;
	p_d = 21;
	q_d = 0;

	L = NcT * pck_len;
	cc1 = 0.0 + I * 0.0;
	cc2 = 0.0 + I * 0.0;
	cc3 = 0.0 + I * 0.0;

	for (k = 1; k <= L; k++)
	{
		cc1 = cc1 + data_in1[offset * NcT + k - 1] * conj(data_in2[offset * NcT + k - 1]);
		cc1 = quantize(creal(cc1), p_cc, q_cc) + I * quantize(cimag(cc1), p_cc,
															  q_cc);

		cc2 = cc2 + data_in1[offset * NcT + k - 1] * conj(data_in2[offset * NcT + k - 2]);
		cc2 = quantize(creal(cc2), p_cc, q_cc) + I * quantize(cimag(cc2), p_cc,
															  q_cc);

		cc3 = cc3 + data_in1[offset * NcT + k - 1] * conj(data_in2[offset * NcT + k]);
		cc3 = quantize(creal(cc3), p_cc, q_cc) + I * quantize(cimag(cc3), p_cc,
															  q_cc);
	}

	//	num = creal(cc1*conj(-cc2+cc3)) / 2.0;
	cc1_tmp = quantize(creal(cc1), p_n, q_n) + I * quantize(cimag(cc1), p_n,
															q_n);
	cc2_tmp = quantize(creal(cc2), p_n, q_n) + I * quantize(cimag(cc2), p_n,
															q_n);
	cc3_tmp = quantize(creal(cc3), p_n, q_n) + I * quantize(cimag(cc3), p_n,
															q_n);
	num = creal(cc1 * conj(-cc2 + cc3)) / 2.0;
	num = quantize(num, p_n, q_n);

	//	den = pow(cabs(-cc2/2.0+cc3/2.0),2) + creal(cc1*conj(cc2+cc3-2.0*cc1));
	cc1_tmp = quantize(creal(cc1), p_d, q_d) + I * quantize(cimag(cc1), p_d,
															q_d);
	cc2_tmp = quantize(creal(cc2), p_d, q_d) + I * quantize(cimag(cc2), p_d,
															q_d);
	cc3_tmp = quantize(creal(cc3), p_d, q_d) + I * quantize(cimag(cc3), p_d,
															q_d);
	den = quantize(pow(cabs(-cc2 / 2.0 + cc3 / 2.0), 2), p_d, q_d);
	den = den + quantize(creal(cc1 * conj(cc2 + cc3 - 2.0 * cc1)), p_d, q_d);
	den = quantize(den, p_d, q_d);

	*time = -num / (den * (float)NcT);
	*time = quantize(*time, p_t, q_t);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafasamp_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					float *data_in2_i, float *ejfasestim_r, float *ejfasestim_i,
					float *ampstim, int *offset, int *NcT, int *pck_len, int *Ta, int *Tf,
					int stfas)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;
	double complex *ejfasestim;

	data_in1 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_in2 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	ejfasestim = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimafasamp(data_in1, data_in2, ejfasestim, ampstim, *offset, *NcT,
				*pck_len, *Ta, *Tf, stfas);

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		ejfasestim_r[i] = creal(ejfasestim[i]);
		ejfasestim_i[i] = cimag(ejfasestim[i]);
	}

	free(data_in1);
	free(data_in2);

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafasamp(double complex *data_in1, double complex *data_in2, double complex *ejfasestim,
				 float *ampstim, int offset, int NcT, int pck_len, int Ta, int Tf,
				 int stfas)
{
	int L;				  // number of samples per packet
	double complex cc0;	  // summation of input vector (data_in1) rotated by detected vector (data_in2)
	int k, i, endt, step; // suppor variable
	// #CGS_2.7 start
	//	float N0, sgd; // thermal and phase noise variances (arbitrary)
	//	double complex ra[(pck_len + Ta + Tf) * NcT]; // support variable
	//	double complex af[(pck_len + Ta + Tf) * NcT]; // forward recursion
	//	double complex ab[(pck_len + Ta + Tf) * NcT]; // backward recursion
	//	double complex ejt[(pck_len + Ta + Tf) * NcT]; // support variable
	float thetf;
	//	float thetb;
	float beta;

	L = NcT * (pck_len + Ta + Tf) + 3;
	//	N0 = 0.05;
	//	sgd = 0.9;
	// #CGS_2.7 end

	cc0 = 0;
	for (k = NcT * Ta; k < L - NcT * Tf; k++)
	{
		cc0 = cc0 + data_in1[k + (offset - Ta) * NcT] * conj(data_in2[k + (offset - Ta) * NcT]);
	}
	for (k = NcT * Ta; k < L - NcT * Tf; k++)
	{
		ampstim[k] = cabs(cc0) / (float)(L - NcT * (Ta + Tf));
	}

	for (k = 0; k < NcT * (Ta + 3); k++)
	{
		endt = 2;
		cc0 = 0;
		for (i = -endt; i < endt + 1; i++)
		{
			cc0 = cc0 + data_in1[k + i + (offset - Ta) * NcT] * conj(data_in2[k + i + (offset - Ta) * NcT]);
		}
		ampstim[k] = cabs(cc0) / (float)(2 * endt + 1);
	}
	for (k = L - NcT * (Tf); k < L; k++)
	{
		endt = 2;
		cc0 = 0;
		for (i = -endt; i < endt + 1; i++)
		{
			cc0 = cc0 + data_in1[k + i + (offset - Ta) * NcT] * conj(data_in2[k + i + (offset - Ta) * NcT]);
		}
		ampstim[k] = cabs(cc0) / (float)(2 * endt + 1);
	}

	//	*ampstim = cabs(cc0) / (float)L;

	if (stfas == 0)
	{
		for (k = 0; k < L; k++)
		{
			endt = 10;
			step = 3;
			if (k < 2 * NcT * Ta || k > L - NcT * Tf)
			{
				endt = 1;
				step = 1;
			}
			cc0 = 0.0;
			for (i = -endt; i < endt + 1; i++)
			{
				if (k + step * i >= 0 && k + step * i <= L - Tf * NcT)
					cc0 = cc0 + data_in1[k + step * i + (offset - Ta) * NcT] * conj(data_in2[k + step * i + (offset - Ta) * NcT]);
			}
			if (cabs(cc0) != 0)
				ejfasestim[k] = cc0 / cabs(cc0);
		}
	}
	else
	{
		thetf = 0.;
		ejfasestim[0] = 0.;
		for (k = 1; k < L; k++)
		{
			cc0 = data_in1[k + (offset - Ta) * NcT] * conj(ejfasestim[k - 1]);
			if (k < 2 * NcT * Ta || k > L - NcT * Tf)
				beta = 0.5;
			else
				beta = 0.01;
			thetf = thetf + beta * cimag(cc0 * conj(data_in2[k + (offset - Ta) * NcT]));
			ejfasestim[k] = cos(thetf) + I * sin(thetf);
		}
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void stimafasamp_q_c_(float *data_in1_r, float *data_in1_i, float *data_in2_r,
					  float *data_in2_i, float *ejfasestim_r, float *ejfasestim_i,
					  float *ampstim, int *offset, int *NcT, int *pck_len, int *Ta, int *Tf,
					  int stfas, int *p_data1, int *q_data1, int *p_data2, int *q_data2,
					  int *p_rit, int *q_rit, int *p_fas, int *q_fas, int *p_amp, int *q_amp)
{
	int i;
	double complex *data_in1;
	double complex *data_in2;
	double complex *ejfasestim;

	data_in1 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	data_in2 = malloc((*NcT) * (*pck_len) * sizeof(double complex));
	ejfasestim = malloc((*NcT) * (*pck_len) * sizeof(double complex));

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in1[i] = data_in1_r[i] + I * data_in1_i[i];
	}

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		data_in2[i] = data_in2_r[i + (*NcT) + 1] + I * data_in2_i[i + (*NcT) + 1];
	}

	stimafasamp_q(data_in1, data_in2, ejfasestim, ampstim, *offset, *NcT,
				  *pck_len, *Ta, *Tf, stfas, *p_data1, *q_data1, *p_data2, *q_data2,
				  *p_rit, *q_rit, *p_fas, *q_fas, *p_amp, *q_amp);

	for (i = 0; i < (*NcT) * (*pck_len); i++)
	{
		ejfasestim_r[i] = creal(ejfasestim[i]);
		ejfasestim_i[i] = cimag(ejfasestim[i]);
	}

	free(data_in1);
	free(data_in2);

	return;
}

// #CGS_2.7 commentate variabili non utilizzate
///////////////////////////////////////////////////////////////////////
void stimafasamp_q(double complex *data_in1, double complex *data_in2, double complex *ejfasestim,
				   float *ampstim, int offset, int NcT, int pck_len, int Ta, int Tf,
				   int stfas, int p_data1, int q_data1, int p_data2, int q_data2,
				   int p_rit, int q_rit, int p_fas, int q_fas, int p_amp, int q_amp)
{
	int L, endt, i;
	//	float maxim, amptmp;
	float amptmp;
	//	double complex cc1, cc_1, cc0;
	double complex cc0;
	//	double complex *z1;
	//	double complex *z_1;
	//	double complex *z0;
	int p_cc1, q_cc1, p_cc2, q_cc2, p_div, q_div;
	//	float z_re, z_im, data1_re, data1_im, data2_re, data2_im, cc_re, cc_im;
	float cc_re, cc_im;
	float abs_cc0, div;
	//	float N0, sgd, ra_re, ra_im, af_re, af_im, ab_re, ab_im, ej_re, ej_im;
	int k, step;
	//	double complex ra[(pck_len + Ta + Tf) * NcT]; // support variable
	//	double complex af[(pck_len + Ta + Tf) * NcT]; // forward recursion
	//	double complex ab[(pck_len + Ta + Tf) * NcT]; // backward recursion
	//	float thetf, thetb, beta;
	float thetf, beta;
	p_cc1 = 14;
	q_cc1 = 9;
	p_cc2 = 7;
	q_cc2 = 9;
	p_div = 1;
	q_div = 15;

	L = NcT * (pck_len + Ta + Tf) + 3;
	//	N0 = 0.05;
	//	sgd = 0.9;

	cc0 = 0;
	for (k = NcT * Ta; k < L - NcT * Tf; k++)
	{
		cc0 = cc0 + data_in1[k + (offset - Ta) * NcT] * conj(data_in2[k + (offset - Ta) * NcT]);
		cc_re = quantize(creal(cc0), p_cc1, q_cc1);
		cc_im = quantize(cimag(cc0), p_cc1, q_cc1);
		cc0 = cc_re + I * cc_im;
	}
	abs_cc0 = cordic_getmag(creal(cc0), cimag(cc0), CORDIC_ITER);
	abs_cc0 = quantize(abs_cc0, p_cc1, q_cc1);
	div = quantize(1. / (float)(L - NcT * (Ta + Tf)), p_div, q_div);
	amptmp = abs_cc0 * div;
	for (k = NcT * Ta; k < L - NcT * Tf; k++)
	{
		ampstim[k] = quantize(amptmp, p_amp, q_amp);
	}

	endt = 2;
	div = quantize(1. / (float)(2 * endt + 1), p_div, q_div);
	for (k = 0; k < NcT * (Ta + 3); k++)
	{
		cc0 = 0;
		for (i = -endt; i < endt + 1; i++)
		{
			cc0 = cc0 + data_in1[k + i + (offset - Ta) * NcT] * conj(data_in2[k + i + (offset - Ta) * NcT]);
			cc_re = quantize(creal(cc0), p_cc1, q_cc1);
			cc_im = quantize(cimag(cc0), p_cc1, q_cc1);
			cc0 = cc_re + I * cc_im;
		}
		abs_cc0 = cordic_getmag(creal(cc0), cimag(cc0), CORDIC_ITER);
		abs_cc0 = quantize(abs_cc0, p_cc1, q_cc1);
		ampstim[k] = abs_cc0 * div;
		ampstim[k] = quantize(ampstim[k], p_amp, q_amp);
	}
	for (k = L - NcT * Tf; k < L; k++)
	{
		cc0 = 0;
		for (i = -endt; i < endt + 1; i++)
		{
			cc0 = cc0 + data_in1[k + i + (offset - Ta) * NcT] * conj(data_in2[k + i + (offset - Ta) * NcT]);
			cc_re = quantize(creal(cc0), p_cc1, q_cc1);
			cc_im = quantize(cimag(cc0), p_cc1, q_cc1);
			cc0 = cc_re + I * cc_im;
		}
		abs_cc0 = cordic_getmag(creal(cc0), cimag(cc0), CORDIC_ITER);
		abs_cc0 = quantize(abs_cc0, p_cc1, q_cc1);
		ampstim[k] = abs_cc0 * div;
		ampstim[k] = quantize(ampstim[k], p_amp, q_amp);
	}

	if (stfas == 0)
	{
		for (k = 0; k < L; k++)
		{
			endt = 10;
			step = 3;
			if (k < 2 * NcT * Ta || k > L - NcT * Tf)
			{
				endt = 1;
				step = 1;
			}
			cc0 = 0.0;
			for (i = -endt; i <= endt; i++)
			{
				cc0 = cc0 + data_in1[k + step * i + (offset - Ta) * NcT] * conj(data_in2[k + step * i + (offset - Ta) * NcT]);
				cc_re = quantize(creal(cc0), p_cc2, q_cc2);
				cc_im = quantize(cimag(cc0), p_cc2, q_cc2);
				cc0 = cc_re + I * cc_im;
			}
			abs_cc0 = cordic_getmag(creal(cc0), cimag(cc0), CORDIC_ITER);
			abs_cc0 = quantize(abs_cc0, p_cc2, q_cc2);
			if (abs_cc0 != 0)
			{
				abs_cc0 = quantize(1. / abs_cc0, 9, 9);
				ejfasestim[k] = cc0 * abs_cc0;
			}
			ejfasestim[k] = quantize(creal(ejfasestim[k]), p_fas, q_fas) + I * quantize(cimag(ejfasestim[k]), p_fas, q_fas);
		}
	}
	else
	{
		thetf = 0.;
		ejfasestim[0] = 0.;
		for (k = 1; k < L; k++)
		{
			if (k < 2 * NcT * Ta || k > L - NcT * Tf)
				beta = 0.5;
			else
				beta = 0.01;
			cc0 = data_in1[k + (offset - Ta) * NcT] * conj(ejfasestim[k - 1]);

			thetf = thetf + beta * cimag(cc0 * conj(data_in2[k + (offset - Ta) * NcT]));

			ejfasestim[k] = cos(thetf) + I * sin(thetf);
		}
	}

	return;
}

// #TF_2.6 Novel postprocessing (check_frame_2 controlled by Icrc)
// #CGS_2.7 added parameters with error values read from .dat files
// #AU_2.8 new start flag search algorithm that includes bits from the training sequence
///////////////////////////////////////////////////////////////////////
int check_detected(int *frame, float *LLR, int frame_len, int data_len,
				   int max_delay, int max_stuffing, int Ta, int *detected_position,
				   int *detected_frame_len, int T_len, int Icrc,
				   int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
				   int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
				   int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin)
{
	const int training_seq_len = 24;												 // training sequence length
	const int flag_len = 8;															 // start/stop flag length
	const int crc_len = 16;															 // CRC field length
																					 //	const int bit_stuffing_len = 5; // number of consecutive '1' before stuffing //#CGS_2.7
	const int data_len_stuf = frame_len - training_seq_len - 2 * flag_len - crc_len; // data field length, considering stuffing bits
	const int flag_err_max = 1;														 // maximum number of mismatching bits in start flag
	int training_seq[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
						  1, 0, 1, 0, 1, 0, 1};	 // training sequence definition
	int start_flag[] = {0, 1, 1, 1, 1, 1, 1, 0}; // start flag definition
												 //	int stop_flag[] = { 0, 1, 1, 1, 1, 1, 1, 0 }; // stop flag definition  //#CGS_2.7
												 //	int offset = training_seq_len + flag_len; // data field starting point //#CGS_2.7
	int stuffing_bits;							 // number of stuffing bits
	int couple1, couple2;						 // position of the couples to be inverted in post-processing
	float min1, min2;							 // minimum reliability (couple1 and couple2 respectively)
	float LLRT[frame_len + max_delay + max_stuffing + T_len];
	int flag_errors; // error counter
	int res;		 // function result
	int sumt;		 // training sequence identification variable
					 //	int i, j, crcn, sum_crc; // support variables //#CGS_2.7
	int i, j;		 // support variables

	int frame2[frame_len + max_delay + max_stuffing + T_len]; // #CGS_2.4

	memset(LLRT, 0.0, sizeof(LLRT)); // #CGS_2.7

	res = 0;

	// #AU_2.8
	//  Cycle on the 2 possible values of the training sequence
	//  The cycle is entered once for the old algorithm (DISABLE_FLAG_ALG=0), twice for the new one (DISABLE_FLAG_ALG=1)
	int ts; // support variable
	for (ts = 0; ts < DISABLE_FLAG_ALG + 1; ts++)
	{

		// scan for all possible starting position, until a good frame is detected
		for (i = 0; i < max_delay + Ta; i++)
		{

			////#CGS_2.4corretto il controllo dello start flag
			int temp;

			if (!DISABLE_FLAG_ALG) // possibilita� di disattivare il controllo preciso
			{
				flag_errors = 0;
				temp = 0;
				for (j = 0; j < flag_len; j++)
					temp = temp + frame[training_seq_len + i + j] * pow(2, j);
				switch (temp)
				{
				case 126:
					break;
				case 254:
					flag_errors++;
					break;
				case 62:
					flag_errors++;
					break;
				case 222:
					flag_errors++;
					break;
				case 46:
					flag_errors++;
					break;
				case 86:
					flag_errors++;
					break;
				case 106:
					flag_errors++;
					break;
				case 116:
					flag_errors++;
					break;
				case 123:
					flag_errors++;
					break;
				case 124: // #AU_2.8 replaced wrong 125 value
					flag_errors++;
					break;
				case 127:
					flag_errors++;
					break;
				default:
					flag_errors = 2;
					break;
				}
			}
			else
			{
				flag_errors = 0;

				// #AU_2.8 new algorithm
				//  training sequence part
				int nts = 6; // number of training sequence bits considered
				for (j = 0; j < nts; j++)
				{
					if (frame[training_seq_len - nts + i + j] != (training_seq[training_seq_len - nts + j] + ts) % 2)
					{
						flag_errors++;
					}
				}

				// start flag part - the same for both training sequences
				for (j = 0; j < flag_len; j++)
				{
					if (frame[training_seq_len + i + j] != start_flag[j])
					{
						flag_errors++;
					}
				}
				int errmax = 3;											 // max number of errors we admit in the sequence
				flag_errors = flag_errors / (errmax + 1) + flag_err_max; // flag_errors==flag_err_max for up to errmax errors, flag_errors>flag_err_max otherwise
																		 // end of new algorithm
			}

			// if there are too many errors,
			// try again

			if (flag_errors > flag_err_max)
			{
				continue;
			}

			// #TF_2.3 #CGS_2.4 J=0
			for (j = 0; j < frame_len + max_delay + max_stuffing + T_len; j++)
			{
				frame2[j] = frame[j];
			}

			// otherwise...
			// #TF_2.6 old postprocessing called with Icrc=1

			if (Icrc == 1)
			{

				// check frame without bit flipping
				if (check_frame(&frame2[i], frame_len, data_len, &stuffing_bits) == 1)
				{
					res = 1;
					*detected_position = i;
					*detected_frame_len = frame_len + stuffing_bits;
					goto loop;
					// break;
				}
				else
				{
					// detect least reliable bits
					// initialization
					min1 = fabs(LLR[training_seq_len + flag_len + i]) + fabs(
																			LLR[training_seq_len + flag_len + i + 1]);
					couple1 = training_seq_len + flag_len + i;
					// scan data and CRC fields
					for (j = training_seq_len + flag_len + i + 1; j < training_seq_len + flag_len + data_len_stuf + crc_len + max_stuffing + i; j++)
					{
						if (min1 > fabs(LLR[j]) + fabs(LLR[j + 1]))
						{
							min1 = fabs(LLR[j]) + fabs(LLR[j + 1]);
							couple1 = j;
						}
					}
					// flip least reliable bits
					frame2[couple1] = (frame2[couple1] + 1) % 2;
					frame2[couple1 + 2] = (frame2[couple1 + 2] + 1) % 2;
					// perform frame check after the first bit flipping
					if (check_frame(&frame2[i], frame_len, data_len, &stuffing_bits) == 1)
					{
						res = 1;
						*detected_position = i;
						*detected_frame_len = frame_len + stuffing_bits;
						goto loop;
						// break;
					}
					else
					{
						// detect least reliable bits
						// initialization
						if (couple1 == i)
						{
							min2 = fabs(LLR[training_seq_len + flag_len + i + 1]) + fabs(
																						LLR[training_seq_len + flag_len + i + 2]);
							couple2 = training_seq_len + flag_len + i + 1;
						}
						else
						{
							min2 = fabs(LLR[training_seq_len + flag_len + i]) + fabs(
																					LLR[training_seq_len + flag_len + i + 1]);
							couple2 = training_seq_len + flag_len + i;
						}
						// scan data and CRC fields
						for (j = training_seq_len + flag_len + i + 1; j < training_seq_len + flag_len + data_len_stuf + crc_len + max_stuffing + i; j++)
						{
							if (couple1 != j)
							{
								if (min2 > fabs(LLR[j]) + fabs(LLR[j + 1]))
								{
									min2 = fabs(LLR[j]) + fabs(LLR[j + 1]);
									couple2 = j;
								}
							}
						}
						// flip least reliable bits
						frame2[couple2] = (frame2[couple2] + 1) % 2;
						frame2[couple2 + 2] = (frame2[couple2 + 2] + 1) % 2;
						// perform frame check after second bit flipping
						if (check_frame(&frame2[i], frame_len, data_len, &stuffing_bits) == 1)
						{
							res = 1;
							*detected_position = i;
							*detected_frame_len = frame_len + stuffing_bits;
							goto loop;
							// break;
						}
					}
				}
			}
			else
			{
				////////////////////////////////////////////////////////////////////////////////
				// error correction through syndrome check (novel postprocessing) #TF_2.6
				////////////////////////////////////////////////////////////////////////////////
				// check frame without bit flipping

				DiffLLR(frame_len + max_delay + max_stuffing + T_len, LLR, LLRT);

				// #CGS_2.7 added parameters with error values read from .dat files

				if (check_frame_2(&frame2[i], &LLRT[i], frame_len, data_len, max_delay, max_stuffing, T_len, &stuffing_bits,
								  err_sing, err_sing_sin, err_doppi, err_doppi_sin, err_sing2c,
								  err_sing2c_sin, err2c1, err2c1_sin, err3, err3_sin, err2c2,
								  err2c2_sin, err2c3, err2c3_sin, err2c4, err2c4_sin) == 1)
				{
					res = 1;
					*detected_position = i;
					*detected_frame_len = frame_len + stuffing_bits;
					goto loop;
					// break;
				}
				else
				{
					// detect least reliable bits
					// initialization
					min1 = fabs(LLR[training_seq_len + flag_len + i]) + fabs(
																			LLR[training_seq_len + flag_len + i + 1]);
					couple1 = training_seq_len + flag_len + i;
					// scan data and CRC fields
					for (j = training_seq_len + flag_len + i + 1; j < training_seq_len + flag_len + data_len_stuf + crc_len + max_stuffing + i; j++)
					{
						if (min1 > fabs(LLR[j]) + fabs(LLR[j + 1]))
						{
							min1 = fabs(LLR[j]) + fabs(LLR[j + 1]);
							couple1 = j;
						}
					}
					// flip least reliable bits
					frame2[couple1] = (frame2[couple1] + 1) % 2;
					frame2[couple1 + 2] = (frame2[couple1 + 2] + 1) % 2;
					// perform frame check after the first bit flipping
					// #CGS_2.7 added parameters with error values read from .dat files
					if (check_frame_2(&frame2[i], &LLRT[i], frame_len, data_len, max_delay, max_stuffing, T_len, &stuffing_bits,
									  err_sing, err_sing_sin, err_doppi, err_doppi_sin, err_sing2c,
									  err_sing2c_sin, err2c1, err2c1_sin, err3, err3_sin, err2c2,
									  err2c2_sin, err2c3, err2c3_sin, err2c4, err2c4_sin) == 1)
					{
						res = 1;
						*detected_position = i;
						*detected_frame_len = frame_len + stuffing_bits;
						goto loop;
						// break;
					}
					else
					{
						// detect least reliable bits
						// initialization
						if (couple1 == i)
						{
							min2 = fabs(LLR[training_seq_len + flag_len + i + 1]) + fabs(
																						LLR[training_seq_len + flag_len + i + 2]);
							couple2 = training_seq_len + flag_len + i + 1;
						}
						else
						{
							min2 = fabs(LLR[training_seq_len + flag_len + i]) + fabs(
																					LLR[training_seq_len + flag_len + i + 1]);
							couple2 = training_seq_len + flag_len + i;
						}
						// scan data and CRC fields
						for (j = training_seq_len + flag_len + i + 1; j < training_seq_len + flag_len + data_len_stuf + crc_len + max_stuffing + i; j++)
						{
							if (couple1 != j)
							{
								if (min2 > fabs(LLR[j]) + fabs(LLR[j + 1]))
								{
									min2 = fabs(LLR[j]) + fabs(LLR[j + 1]);
									couple2 = j;
								}
							}
						}
						// flip least reliable bits
						frame2[couple2] = (frame2[couple2] + 1) % 2;
						frame2[couple2 + 2] = (frame2[couple2 + 2] + 1) % 2;
						// perform frame check after second bit flipping
						// #CGS_2.7 added parameters with error values read from .dat files
						if (check_frame_2(&frame2[i], &LLRT[i], frame_len, data_len, max_delay, max_stuffing, T_len, &stuffing_bits,
										  err_sing, err_sing_sin, err_doppi, err_doppi_sin, err_sing2c,
										  err_sing2c_sin, err2c1, err2c1_sin, err3, err3_sin, err2c2,
										  err2c2_sin, err2c3, err2c3_sin, err2c4, err2c4_sin) == 1)
						{
							res = 1;
							*detected_position = i;
							*detected_frame_len = frame_len + stuffing_bits;
							goto loop;
							// break;
						}
					}
				}
			}
		}
	} // #AU_2.8 end of cycle on ts

// #CGS_2.4 tolto
// for (i = 0; i < frame_len; i++) {
//	frame[i] = frame2[i];
// }
loop:;
	// restore training sequence, start flag and stop flag
	// since they may contain wrong bit
	if (res == 1)
	{
		// #CGS_2.5 aggiunto
		for (i = 0; i < frame_len + max_delay + max_stuffing + T_len; i++)
		{
			frame[i] = frame2[i];
		}

		// training sequence
		sumt = 0; // training sequence identification
		for (i = 0; i < training_seq_len; i++)
		{
			sumt += mod(frame[i + *detected_position] + training_seq[i], 2);
		}
		for (i = 0; i < training_seq_len; i++)
		{
			if (sumt < training_seq_len)
			{
				frame[i + *detected_position] = mod(training_seq[i] + 0, 2);
			}
			else
			{
				frame[i + *detected_position] = mod(training_seq[i] + 1, 2);
			}
		}
		// start and stop flag
		for (i = 0; i < flag_len; i++)
		{
			frame[training_seq_len + *detected_position + i] = start_flag[i];
			frame[training_seq_len + *detected_position + flag_len + data_len // #CGS_2.5
				  + crc_len + stuffing_bits + i] = start_flag[i];
		}
	}
	else
	{
		*detected_position = 0;
		*detected_frame_len = 0;
	}

	// return 0xFFFFFFFF if a valid codeword has not been identified
	// else return the starting bit index
	return res;
}

// #TF_2.6
///////////////////////////////////////////////////////////////////////
// #CGS_2.7 added parameters with error values read from .dat files
// #AU_2.8 new stop flag search algorithm
// int check_frame_2(int *frame, float *LLRT, int frame_len, int data_len, int max_delay, int max_stuffing, int T_len, int *stuffing_bits){
int check_frame_2(int *frame, float *LLRT, int frame_len, int data_len, int max_delay, int max_stuffing, int T_len, int *stuffing_bits,
				  int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c,
				  int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2,
				  int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin)
{

	const int training_seq_len = 24;												 // training sequence length
	const int flag_len = 8;															 // flag length
	const int crc_len = 16;															 // CRC length
	const int bit_stuffing_len = 5;													 // number of consecutive '1' before stuffing
	const int err_max = 1;															 // maximum number of bit mismatch //#TF_2.3
	const int data_len_stuf = frame_len - training_seq_len - 2 * flag_len - crc_len; // data field length, considering bit stuffing
	int bit_stuffing_cntr;															 // stuffing sequence counter
	int start_flag[] = {0, 1, 1, 1, 1, 1, 1, 0};									 // start/stop flag
																					 //	int crcn,nn,ind,e1,e2,sum_crc,index,res,flag_errors;   //#CGS_2.7
	int crcn, nn, ind, sum_crc, res, flag_errors;									 // #CGS_2.7
	int decodOUT[frame_len + max_stuffing + T_len], crcout[16], sum;
	const int exact = 58096; // syndrome without errors
	const int lt = data_len_stuf + crc_len;
	float LLR2[frame_len + max_delay + max_stuffing + T_len];
	// #CGS_2.7 start
	/*	int err_sing_sin[4],err_sing[4];
		int err_doppi_sin[182],err_sing2c_sin[4];
		int err_doppi[182*2],err_sing2c[4];
		int err2c1_sin[8648],err2c1[8648ƒs*4],err3[2*3],err3_sin[2];
		int err2c2_sin[2862],err2c2[2862*8];
		int err2c3_sin[477],err2c3[477*12];
		int err2c4_sin[167],err2c4[167*16];
	*/
	// #CGS_2.7 end
	const float LLRthld = 25.;

	int i, j; // support variable

	memset(decodOUT, 0, sizeof(decodOUT)); // #CGS_2.7
	memset(LLR2, 0.0, sizeof(LLR2));	   // #CGS_2.7

	// remove bit stuffing
	// data and CRC fields only
	bit_stuffing_cntr = 0;
	*stuffing_bits = 0;

	for (i = 0; i < data_len + crc_len; i++)
	{
		decodOUT[i] = frame[training_seq_len + flag_len + *stuffing_bits + i];
		LLR2[i] = -LLRT[training_seq_len + flag_len + *stuffing_bits + i];
		if (frame[training_seq_len + flag_len + *stuffing_bits + i] == 1)
		{
			bit_stuffing_cntr++;
			// after 5 contiguous '1' has been detected
			if (bit_stuffing_cntr == bit_stuffing_len)
			{
				// update stuffing counters
				bit_stuffing_cntr = 0;
				*stuffing_bits = *stuffing_bits + 1;
				// check if stuffing bit is 0 or 1
				// exit in case of 1
				if (frame[training_seq_len + flag_len + *stuffing_bits + i] == 1)
				{
					return 0;
				}
			}
		}
		else
		{
			bit_stuffing_cntr = 0;
		}
	}

	crcn = 0;
	sum_crc = 0;
	res = 0;

	// #CGS_2.7 start
	// error files reading
	//	readfile(err_sing,err_sing_sin,err_doppi,err_doppi_sin,err_sing2c
	//	,err_sing2c_sin,err2c1,err2c1_sin,err3,err3_sin,err2c2,err2c2_sin
	//	,err2c3,err2c3_sin,err2c4,err2c4_sin);
	// #CGS_2.7 end

	// when just one error couple is found, correction is always performed

	// int mcrc[16];
	// res = crc_check_SG(decodOUT, mcrc, data_len);
	// if(res == 1){
	//	goto loop;
	//	};

	sum = 0; // #CGS_2.7
	crc_check_tot(lt, decodOUT, sum, crcout);

	crcn++;
	nn = natConv(16, crcout);
	if (nn == exact)
	{
		res = 1;
		goto loop;
	}

	ind = 0;

	ind = FindIndex(4, nn, err_sing_sin);
	if (ind != 0)
	{
		decodOUT[err_sing[ind]] = (decodOUT[err_sing[ind]] + 1) % 2;
		res = 1;
		goto loop;
	}

	ind = FindIndex(182, nn, err_doppi_sin);
	if (ind != 0)
	{
		decodOUT[err_doppi[ind - 1] - 1] = (decodOUT[err_doppi[ind - 1] - 1] + 1) % 2;
		decodOUT[err_doppi[ind + 182 - 1] - 1] = (decodOUT[err_doppi[ind + 182 - 1] - 1] + 1) % 2;
		res = 1;
		goto loop;
	}

	ind = FindIndex(4, nn, err_sing2c_sin);
	if (ind != 0)
	{
		decodOUT[err_sing2c[ind - 1] - 1] = (decodOUT[err_sing2c[ind - 1] - 1] + 1) % 2;
		res = 1;
		goto loop;
	}
	// when two error couples are found, correction is performed only if LLRs of the couple have low reliability
	// In this way missed detections and performance loss are avoided

	ind = FindIndex(8648, nn, err2c1_sin);
	if (ind != 0)
	{
		if (((fabs(LLR2[err2c1[ind - 1] - 1]) + fabs(LLR2[err2c1[ind + 8648 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c1[ind + 2 * 8648 - 1] - 1]) + fabs(LLR2[err2c1[ind + 3 * 8648 - 1] - 1])) < LLRthld))
		{
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c1[ind + j * 8648 - 1] - 1] = (decodOUT[err2c1[ind + j * 8648 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c1[ind + j * 8648 - 1] - 1] = (decodOUT[err2c1[ind + j * 8648 - 1] - 1]) % 2;
			}
		}
	}

	ind = FindIndex(2, nn, err3_sin);
	if (ind != 0)
	{
		if (((fabs(LLR2[err3[ind - 1] - 1]) + fabs(LLR2[err3[ind + 2 * 2 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err3[ind + 2 - 1] - 1])) < LLRthld))
		{
			for (j = 0; j < 3; j++)
			{
				decodOUT[err3[ind + j * 2 - 1] - 1] = (decodOUT[err3[ind + j * 2 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 0; j < 3; j++)
			{
				decodOUT[err3[ind + j * 2 - 1] - 1] = (decodOUT[err3[ind + j * 2 - 1] - 1]) % 2;
			}
		}
	}

	ind = FindIndex(2862, nn, err2c2_sin);
	if (ind != 0)
	{
		if (((fabs(LLR2[err2c2[ind - 1] - 1]) + fabs(LLR2[err2c2[ind + 2862 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c2[ind + 2 * 2862 - 1] - 1]) + fabs(LLR2[err2c2[ind + 3 * 2862 - 1] - 1])) < LLRthld))
		{
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c2[ind + j * 2862 - 1] - 1] = (decodOUT[err2c2[ind + j * 2862 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c2[ind + j * 2862 - 1] - 1] = (decodOUT[err2c2[ind + j * 2862 - 1] - 1]) % 2;
			}
		}
		if (((fabs(LLR2[err2c2[ind + 4 * 2862 - 1] - 1]) + fabs(LLR2[err2c2[ind + 5 * 2862 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c2[ind + 6 * 2862 - 1] - 1]) + fabs(LLR2[err2c2[ind + 7 * 2862 - 1] - 1])) < LLRthld))
		{
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c2[ind + j * 2862 - 1] - 1] = (decodOUT[err2c2[ind + j * 2862 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c2[ind + j * 2862 - 1] - 1] = (decodOUT[err2c2[ind + j * 2862 - 1] - 1]) % 2;
			}
		}
	}

	ind = FindIndex(477, nn, err2c3_sin);
	if (ind != 0)
	{
		if (((fabs(LLR2[err2c3[ind - 1] - 1]) + fabs(LLR2[err2c3[ind + 477 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c3[ind + 2 * 477 - 1] - 1]) + fabs(LLR2[err2c3[ind + 3 * 477 - 1] - 1])) < LLRthld))
		{
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1]) % 2;
			}
		}
		if (((fabs(LLR2[err2c3[ind + 4 * 477 - 1] - 1]) + fabs(LLR2[err2c3[ind + 5 * 477 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c3[ind + 6 * 477 - 1] - 1]) + fabs(LLR2[err2c3[ind + 7 * 477 - 1] - 1])) < LLRthld))
		{
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1]) % 2;
			}
		}
		if (((fabs(LLR2[err2c3[ind + 8 * 477 - 1] - 1]) + fabs(LLR2[err2c3[ind + 9 * 477 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c3[ind + 10 * 477 - 1] - 1]) + fabs(LLR2[err2c3[ind + 11 * 477 - 1] - 1])) < LLRthld))
		{
			for (j = 8; j < 12; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 8; j < 12; j++)
			{
				decodOUT[err2c3[ind + j * 477 - 1] - 1] = (decodOUT[err2c3[ind + j * 477 - 1] - 1]) % 2;
			}
		}
	}

	ind = FindIndex(167, nn, err2c4_sin);
	if (ind != 0)
	{
		if (((fabs(LLR2[err2c4[ind - 1] - 1]) + fabs(LLR2[err2c4[ind + 167 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c4[ind + 2 * 167 - 1] - 1]) + fabs(LLR2[err2c4[ind + 3 * 167 - 1] - 1])) < LLRthld))
		{
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 0; j < 4; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1]) % 2;
			}
		}
		if (((fabs(LLR2[err2c4[ind + 4 * 167 - 1] - 1]) + fabs(LLR2[err2c4[ind + 5 * 167 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c4[ind + 6 * 167 - 1] - 1]) + fabs(LLR2[err2c4[ind + 7 * 167 - 1] - 1])) < LLRthld))
		{
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 4; j < 8; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1]) % 2;
			}
		}
		if (((fabs(LLR2[err2c4[ind + 12 * 167 - 1] - 1]) + fabs(LLR2[err2c4[ind + 13 * 167 - 1] - 1])) < LLRthld) && ((fabs(LLR2[err2c4[ind + 14 * 167 - 1] - 1]) + fabs(LLR2[err2c4[ind + 15 * 167 - 1] - 1])) < LLRthld))
		{
			for (j = 12; j < 16; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1] + 1) % 2;
			}
			crc_check_3(lt, data_len, decodOUT, &sum_crc); // #CGS_2.7 sum_crc by reference
			crcn++;
			if (sum_crc == 0)
			{
				res = 1;
				goto loop;
			}
			for (j = 12; j < 16; j++)
			{
				decodOUT[err2c4[ind + j * 167 - 1] - 1] = (decodOUT[err2c4[ind + j * 167 - 1] - 1]) % 2;
			}
		}
	}

loop:;

	// stuffing bits restoration
	if (res == 1)
	{
		for (i = 0; i < data_len + crc_len; i++)
		{
			frame[training_seq_len + flag_len + i] = decodOUT[i];
		}
		return 1;
	}
	else
	{
		return 0;
	}

	return 0; // #CGS_2.7
}

// #TF_2.6
///////////////////////////////////////////////////////////////////////
// void crc_check_3(int lt, int l, int *msgDEC, int sum){  //#CGS_2.7
void crc_check_3(int lt, int l, int *msgDEC, int *sum)
{ // #CGS_2.7
	const int crc_len = 16;
	int i, j, temp;
	int crcDEC[crc_len], dati[l], check[crc_len];

	*sum = 0; // #CGS_2.7
	temp = 0;
	for (i = 0; i < crc_len; i++)
	{
		check[i] = 0;
		crcDEC[i] = msgDEC[lt - i];
	}
	for (i = 0; i < l; i++)
	{
		dati[i] = msgDEC[i];
	}
	for (j = 0; j < l; j++)
	{
		temp = (check[crc_len - 1] + dati[j]) % 2;
		for (i = crc_len - 1; i > 0; i--)
		{
			if ((i == 12) || (i == 5))
			{
				check[i] = (check[i - 1] + temp + 1) % 2;
			}
			else
			{
				check[i] = check[i - 1];
			}
		}
		check[0] = temp;
	}
	for (i = 0; i < crc_len; i++)
	{
		//		sum += (check[i] + crcDEC[i]) % 2;    //#CGS_2.7
		*sum += (check[i] + crcDEC[i]) % 2; // #CGS_2.7
	}

	return;
}

// TF_2.6
// #CGS_2.7 change return value from void to integer (1 OK, 0 errors)
///////////////////////////////////////////////////////////////////////
int readfile(int *err_sing, int *err_sing_sin, int *err_doppi, int *err_doppi_sin, int *err_sing2c, int *err_sing2c_sin, int *err2c1, int *err2c1_sin, int *err3, int *err3_sin, int *err2c2, int *err2c2_sin, int *err2c3, int *err2c3_sin, int *err2c4, int *err2c4_sin)
{
	FILE *datafile;
	int e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, j;

	datafile = fopen("error_files/err_singoli_ord.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err_singoli_ord.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 4; j++)
	{
		fscanf(datafile, "%d %d", &e1, &e2);
		err_sing[j] = e1;
		err_sing_sin[j] = e2;
	}
	fclose(datafile);

	datafile = fopen("error_files/err_doppi_ord.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err_doppi_ord.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 182; j++)
	{
		fscanf(datafile, "%d %d %d", &e1, &e2, &e3);
		err_doppi[j] = e1;
		err_doppi[182 + j] = e2;
		err_doppi_sin[j] = e3;
	}
	fclose(datafile);

	datafile = fopen("error_files/err_sing2c_ord.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err_sing2c_ord.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 4; j++)
	{
		fscanf(datafile, "%d %d", &e1, &e2);
		err_sing2c[j] = e1;
		err_sing2c_sin[j] = e2;
	}
	fclose(datafile);

	datafile = fopen("error_files/err2c1.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err2c1.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 8648; j++)
	{
		fscanf(datafile, "%d %d %d %d %d", &e1, &e2, &e3, &e4, &e5);
		err2c1[j] = e1;
		err2c1[j + 8648] = e2;
		err2c1[j + 8648 * 2] = e3;
		err2c1[j + 8648 * 3] = e4;
		err2c1_sin[j] = e5;
	}
	fclose(datafile);

	datafile = fopen("error_files/err3ord.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err3ord.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 2; j++)
	{
		fscanf(datafile, "%d %d %d %d", &e1, &e2, &e3, &e4);
		err3[j] = e1;
		err3[j + 2] = e2;
		err3[j + 2 * 2] = e3;
		err3_sin[j] = e4;
	}
	fclose(datafile);

	datafile = fopen("error_files/err2c2.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err2c2.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 2862; j++)
	{
		fscanf(datafile, "%d %d %d %d %d %d %d %d %d", &e1, &e2, &e3, &e4, &e5, &e6, &e7, &e8, &e9);
		err2c2[j] = e1;
		err2c2[j + 2862] = e2;
		err2c2[j + 2862 * 2] = e3;
		err2c2[j + 2862 * 3] = e4;
		err2c2[j + 2862 * 4] = e5;
		err2c2[j + 2862 * 5] = e6;
		err2c2[j + 2862 * 6] = e7;
		err2c2[j + 2862 * 7] = e8;
		err2c2_sin[j] = e9;
	}
	fclose(datafile);

	datafile = fopen("error_files/err2c3.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err2c3.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 477; j++)
	{
		fscanf(datafile, "%d %d %d %d %d %d %d %d %d %d %d %d %d", &e1, &e2, &e3, &e4, &e5, &e6, &e7, &e8, &e9, &e10, &e11, &e12, &e13);
		err2c3[j] = e1;
		err2c3[j + 477] = e2;
		err2c3[j + 477 * 2] = e3;
		err2c3[j + 477 * 3] = e4;
		err2c3[j + 477 * 4] = e5;
		err2c3[j + 477 * 5] = e6;
		err2c3[j + 477 * 6] = e7;
		err2c3[j + 477 * 7] = e8;
		err2c3[j + 477 * 8] = e9;
		err2c3[j + 477 * 9] = e10;
		err2c3[j + 477 * 10] = e11;
		err2c3[j + 477 * 11] = e12;
		err2c3_sin[j] = e13;
	}
	fclose(datafile);

	datafile = fopen("error_files/err2c4.dat", "r");
	if (datafile == NULL)
	{
		printf("%s doesn't exist\n", "err2c4.dat");
		// exit(0); //#CGS_2.7
		return 0; // #CGS_2.7
	}
	for (j = 0; j < 167; j++)
	{
		fscanf(datafile, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &e1, &e2, &e3, &e4, &e5, &e6, &e7, &e8, &e9, &e10, &e11, &e12, &e13, &e14, &e15, &e16, &e17);
		err2c4[j] = e1;
		err2c4[j + 167] = e2;
		err2c4[j + 167 * 2] = e3;
		err2c4[j + 167 * 3] = e4;
		err2c4[j + 167 * 4] = e5;
		err2c4[j + 167 * 5] = e6;
		err2c4[j + 167 * 6] = e7;
		err2c4[j + 167 * 7] = e8;
		err2c4[j + 167 * 8] = e9;
		err2c4[j + 167 * 9] = e10;
		err2c4[j + 167 * 10] = e11;
		err2c4[j + 167 * 11] = e12;
		err2c4[j + 167 * 12] = e13;
		err2c4[j + 167 * 13] = e14;
		err2c4[j + 167 * 14] = e15;
		err2c4[j + 167 * 15] = e16;
		err2c4_sin[j] = e17;
	}
	fclose(datafile);

	// return;  //#CGS_2.7
	return 1; // #CGS_2.7
}

// TF_2.6
///////////////////////////////////////////////////////////////////////
int FindIndex(int a, int n, int *vett)
{
	int p, u, find, m;
	find = 0;
	u = a;
	p = 1;
	do
	{
		m = (p + u) / 2;
		if (vett[m - 1] == n)
		{
			find = m;
			return find;
		}
		if (vett[m - 1] < n)
		{
			p = m + 1;
		}
		else
		{
			u = m - 1;
		}
	} while (p <= u);
	return find;
}

// TF_2.6
///////////////////////////////////////////////////////////////////////
int natConv(int a, int *bin)
{
	int nat, i;

	nat = 0;
	for (i = 0; i < a; i++)
	{
		nat += bin[a - 1 - i] * pow(2, i);
	}
	return nat;
}

// TF_2.6
///////////////////////////////////////////////////////////////////////
void crc_check_tot(int lt, int *msgDEC, int sum, int *crcout)
{
	const int crc_len = 16;
	int temp, i, j;
	int check[crc_len];

	sum = 0;
	temp = 0;
	for (i = 0; i < crc_len; i++)
	{
		check[i] = 0;
	}
	for (j = 0; j < lt; j++)
	{
		temp = (check[crc_len - 1] + msgDEC[j]) % 2;
		for (i = crc_len - 1; i > 0; i--)
		{
			if ((i == 12) || (i == 5))
			{
				check[i] = (check[i - 1] + temp + 1) % 2;
			}
			else
			{
				check[i] = check[i - 1];
			}
		}
		check[0] = temp;
	}
	for (i = 0; i < crc_len; i++)
	{
		sum += check[i];
	}
	for (i = 0; i < crc_len; i++)
	{
		crcout[i] = check[crc_len - 1 - i];
	}

	return;
}

// TF_2.6
///////////////////////////////////////////////////////////////////////
void DiffLLR(int all_len, float *LLR, float *LLRD)
{
	float modulo, segno;
	float check[all_len];
	int i;

	for (i = 0; i < all_len; i++)
	{
		modulo = fabs(LLR[i]);
		segno = -1.;
		if (LLR[i] >= 0.)
			segno = 1.;
		check[i] = transform(segno, modulo);
	}
	modulo = fabs(check[1]);
	segno = -1.;
	if (check[1] >= 0.)
		segno = 1.;
	LLRD[1] = transform(segno, modulo);
	for (i = 2; i < all_len; i++)
	{
		modulo = fabs(check[i]) + fabs(check[i - 1]);
		if (check[i] * check[i - 1] < 0.)
			segno = -1.;
		if (check[i] * check[i - 1] >= 0.)
			segno = 1.;
		LLRD[i] = transform(segno, modulo);
	}
	return;
}

// TF_2.6
///////////////////////////////////////////////////////////////////////
float transform(float segno, float modulo)
{
	const float b = 15.;
	const float c = 2.;
	const float s = 6.11804637E-07;
	const float d = 0.270670563;
	float a, transf;

	if (modulo >= b)
	{
		transf = s * segno;
	}
	else if (modulo >= c)
	{
		transf = 2 * exp(-modulo) * segno;
	}
	else if (modulo >= d)
	{
		a = exp(-modulo);
		transf = -log((1. - a) / (1. + a)) * segno;
	}
	else if (modulo >= s)
	{
		transf = -segno * log(modulo / 2.);
	}
	else
	{
		transf = b * segno;
	}
	return transf;
}

// #CGS_2.7 commentate variabili non utilizzate
// #AU_2.8 new stop flag search algorithm
///////////////////////////////////////////////////////////////////////
int check_frame(int *frame, int frame_len, int data_len, int *stuffing_bits)
{
	const int training_seq_len = 24; // training sequence length
	const int flag_len = 8;			 // flag length
	const int crc_len = 16;			 // CRC length
	const int bit_stuffing_len = 5;	 // number of consecutive '1' before stuffing
	const int err_max = 1;			 // maximum number of bit mismatch //#TF_2.3
	// const int err_max = 1; // maximum number of bit mismatch //#TF_2.3
	//	const int data_len_stuf = frame_len - training_seq_len - 2* flag_len
	//			- crc_len; // data field length, considering bit stuffing
	int start_flag[] = {0, 1, 1, 1, 1, 1, 1, 0}; // start flag
												 //	int stop_flag[] = { 0, 1, 1, 1, 1, 1, 1, 0 }; // stop flag
	int *clean_frame;							 // frame after destuffing
	int bit_stuffing_cntr;						 // stuffing sequence counter
	int flag_errors, errors2;					 // error counter //#TF_2.3

	int i; // support variable

	clean_frame = malloc((data_len + crc_len) * sizeof(int));
	// remove bit stuffing
	// data and CRC fields only
	bit_stuffing_cntr = 0;
	*stuffing_bits = 0;
	for (i = 0; i < data_len + crc_len; i++)
	{
		clean_frame[i] = frame[training_seq_len + flag_len + *stuffing_bits + i];
		if (frame[training_seq_len + flag_len + *stuffing_bits + i] == 1)
		{
			bit_stuffing_cntr++;
			// after 5 coniguous '1' has been detected
			if (bit_stuffing_cntr == bit_stuffing_len)
			{
				// update stuffing counters
				bit_stuffing_cntr = 0;
				*stuffing_bits = *stuffing_bits + 1;
				// check if stuffing bit is 0 or 1
				// exit in case of 1
				if (frame[training_seq_len + flag_len + *stuffing_bits + i] == 1)
				{
					free(clean_frame); // #CGS_2.4
					return 0;
				}
			}
		}
		else
		{
			bit_stuffing_cntr = 0;
		}
	}

	// perform CRC check
	if (CRC_check_2(clean_frame, data_len) == 1)
	{
		// #TF_2.3
		//  errors1 = 0;
		errors2 = 0;

		// check start flag  //#CGS_2.4 rimosso il doppio controllo dello start flag
		// for (i = 0; i < flag_len-2; i++)
		// {
		// if (frame[training_seq_len+i] != start_flag[i])
		// {
		// if (frame[training_seq_len+i+2] != start_flag[i+2])
		// {
		// errors1++;
		// }
		// }
		// }

		// check stop flag  //#CGS_2.4 corretto il controllo dello stop flag
		// for (i = 0; i < flag_len-2; i++)
		// {
		// if (frame[training_seq_len+flag_len+data_len+crc_len+*stuffing_bits+i] != start_flag[i])
		// {
		// if (frame[training_seq_len+flag_len+data_len+crc_len+*stuffing_bits+i+2] != start_flag[i+2])
		// {
		// errors2++;
		// }
		// }
		// }

		// #CGS_2.4 corretto il controllo dello start flag, fare anche per lo stop
		int j;
		int temp;

		if (!DISABLE_FLAG_ALG) // possibilita� di disattivare il controllo preciso
		{

			flag_errors = 0;
			temp = 0;
			for (j = 0; j < flag_len; j++)
				temp = temp + frame[training_seq_len + flag_len + data_len + crc_len + *stuffing_bits + j] * pow(2, j);
			switch (temp)
			{
			case 126:
				break;
			case 254:
				flag_errors++;
				break;
			case 62:
				flag_errors++;
				break;
			case 222:
				flag_errors++;
				break;
			case 46:
				flag_errors++;
				break;
			case 86:
				flag_errors++;
				break;
			case 106:
				flag_errors++;
				break;
			case 116:
				flag_errors++;
				break;
			case 123:
				flag_errors++;
				break;
			case 124: // #AU_2.8 replaced wrong 125 value
				flag_errors++;
				break;
			case 127:
				flag_errors++;
				break;
			default:
				flag_errors = 2;
				break;
			}
		}
		else
		{
			flag_errors = 0;
			// #AU_2.8 commented and replaced with new algorithm
			// for (j = 0; j < flag_len - 1; j++) {
			//	if (frame[training_seq_len + flag_len + data_len + crc_len
			//			+ *stuffing_bits + j] != start_flag[j]) {
			//		if (frame[training_seq_len + flag_len + data_len + crc_len
			//				+ *stuffing_bits + j + 1] != start_flag[j + 1]) {
			//			flag_errors++;
			//		}
			//	}
			// }

			// #AU_2.8 new algorithm
			for (j = 0; j < flag_len; j++)
			{
				if (frame[training_seq_len + flag_len + data_len + crc_len + *stuffing_bits + j] != start_flag[j])
				{
					flag_errors++;
				}
			}

			int errmax = 3;										// max number of errors we admit in the sequence
			flag_errors = flag_errors / (errmax + 1) + err_max; // flag_errors==err_max for up to errmax errors, flag_errors>err_max otherwise
																// end of new algorithm
		}

		// if the number of wrong bit is too high, check for another option
		// else a valid frame has been identified
		if (flag_errors <= err_max)
		{
			free(clean_frame); // #CGS_2.4
			return 1;
		}
	}

	free(clean_frame);

	return 0;
}

///////////////////////////////////////////////////////////////////////
int CRC_check(int pnFrame[], int L_Frame)
{
	const int CRC_len = 16; // CRC field length
	int FrameLen = L_Frame; // frame length
	int CRC[CRC_len];		// calculated CRC
	int crcFB = 0x11021;	// feedback line
	int stateCRC = 0xFFFF;	// CRC register's state
	int res = 1;			// result variable
	int bit;
	int i; // support variable

	// Calculate CRC
	for (i = 0; i < FrameLen; i++)
	{
		bit = (stateCRC ^ (int)pnFrame[i]) & 1;
		if (bit == 1)
			stateCRC ^= crcFB;
		stateCRC >>= 1;
	}

	for (i = 0; i < CRC_len; i++)
	{
		CRC[i] = stateCRC & 1;
		stateCRC >>= 1;
	}

	// Check CRC
	for (i = 0; i < CRC_len; i++)
	{
		if (CRC[i] != pnFrame[FrameLen + i])
		{
			res = 0;
		}
	}

	return res;
}

///////////////////////////////////////////////////////////////////////
void CRC_calc(int pnFrame[], int L_Frame)
{
	const int CRC_len = 16; // CRC field length
	int FrameLen = L_Frame; // frame length
	int CRC[CRC_len];		// calculated CRC
	int crcFB = 0x11021;	// feedback line
	int stateCRC = 0xFFFF;	// CRC register's state
	int bit;
	int i; // support variable

	// Calculate CRC
	for (i = 0; i < FrameLen; i++)
	{
		bit = (stateCRC ^ (int)pnFrame[i]) & 1;
		if (bit == 1)
			stateCRC ^= crcFB;
		stateCRC >>= 1;
	}

	for (i = 0; i < CRC_len; i++)
	{
		CRC[i] = stateCRC & 1;
		stateCRC >>= 1;
		pnFrame[i + FrameLen] = CRC[i];
	}

	return;
}

///////////////////////////////////////////////////////////////////////
int CRC_check_2(int pnFrame[], int L_Frame)
{
	const int CRC_len = 16; // CRC field length
	int FrameLen = L_Frame; // frame length
	int CRC[CRC_len];		// calculated CRC
	int res = 1;			// result variable
	int i;					// support variable

	CRC_calc_2(pnFrame, CRC, FrameLen, CRC_len);

	// Check CRC
	for (i = 0; i < CRC_len; i++)
	{
		if (CRC[i] != pnFrame[FrameLen + i])
		{
			res = 0;
		}
	}

	return res;
}

///////////////////////////////////////////////////////////////////////
void CRC_calc_2(int pnFrame[], int CRC[], int L_Frame, int CRC_len)
{
	int FrameLen = L_Frame; // frame length
	int bit;
	int i, j; // support variables

	for (i = 0; i < CRC_len; i++)
	{
		CRC[i] = 1;
	}

	// Calculate CRC
	for (i = 0; i < FrameLen; i++)
	{
		bit = CRC[0] ^ pnFrame[i];
		for (j = 0; j < CRC_len - 1; j++)
		{
			if ((j == CRC_len - 1 - 12) || (j == CRC_len - 1 - 5))
			{
				CRC[j] = CRC[j + 1] ^ bit;
			}
			else
			{
				CRC[j] = CRC[j + 1];
			}
		}
		CRC[CRC_len - 1] = bit;
	}
	for (i = 0; i < CRC_len; i++)
	{
		CRC[i] = mod(CRC[i] + 1, 2);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void H_1d(double complex y, int D, float *v)
{
	double complex phase;
	float ang;
	int i;

	// see RD1, eq 3.11
	for (i = 0; i < D; i++)
	{
		ang = 2 * PI * (float)i / (float)D;
		phase = cos(ang) - I * sin(ang);
		v[i] = creal(y * phase);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void H_2d(double complex y, int D, float **v, int i2)
{
	double complex phase;
	float ang;
	int i;

	// see RD1, eq 3.11
	for (i = 0; i < D; i++)
	{
		ang = 2 * PI * (float)i / (float)D;
		phase = cos(ang) - I * sin(ang);
		v[i][i2] = creal(y * phase);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void H_q_1d(double complex y, int D, float *v, int p, int q)
{
	double complex phase;
	float ang;
	int p_ang, q_ang, p_phase, q_phase;
	int i;

	p_ang = 3;
	q_ang = 15;
	p_phase = 7;
	q_phase = 15;

	for (i = 0; i < D; i++)
	{
		ang = 2 * PI * (float)i / (float)D;
		ang = mod(ang, 2 * PI);
		if (ang >= PI)
		{
			ang = ang - 2 * PI;
		}
		ang = quantize(ang, p_ang, q_ang);
		phase = cordic_polartorect(1.0, ang, CORDIC_ITER);
		phase = quantize(creal(phase), p_phase, q_phase) - I * quantize(cimag(
																			phase),
																		p_phase, q_phase);

		v[i] = creal(y * phase);
		v[i] = quantize(v[i], p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void H_q_2d(double complex y, int D, float **v, int i2, int p, int q)
{
	double complex phase;
	float ang;
	int p_ang, q_ang, p_phase, q_phase;
	int i;

	p_ang = 3;
	q_ang = 15;
	p_phase = 7;
	q_phase = 15;

	for (i = 0; i < D; i++)
	{
		ang = 2 * PI * (float)i / (float)D;
		ang = mod(ang, 2 * PI);
		if (ang >= PI)
		{
			ang = ang - 2 * PI;
		}
		ang = quantize(ang, p_ang, q_ang);
		phase = cordic_polartorect(1.0, ang, CORDIC_ITER);
		phase = quantize(creal(phase), p_ang, q_ang) - I * quantize(
															   cimag(phase), p_phase, q_phase);

		v[i][i2] = creal(y * phase);
		v[i][i2] = quantize(v[i][i2], p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
float integr(float *v, int D)
{
	float somr;
	float somr1;
	int i;

	somr = -1e20;
	for (i = 0; i < D; i++)
	{
		somr1 = logsumexp(somr, v[i]);
		somr = somr1;
	}

	return somr;
}

///////////////////////////////////////////////////////////////////////
void normalize(float **v, int D, int i2)
{
	int i;
	float somr, somr1;

	somr = -1e20;

	for (i = 0; i < D; i++)
	{
		somr1 = logsumexp(somr, v[i][i2]);
		somr = somr1;
	}
	for (i = 0; i < D; i++)
	{
		v[i][i2] = v[i][i2] - somr;
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void normalize_q(float **v, int D, int i2, int p, int q)
{
	int i;
	float somr, somr1;

	somr = -1e20;

	for (i = 0; i < D; i++)
	{
		somr1 = logsumexp(somr, v[i][i2]);
		somr = somr1;
	}

	for (i = 0; i < D; i++)
	{
		v[i][i2] = quantize(v[i][i2] - somr, p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void var_exp_1d(float *v1, float **v2, int D, float *v, int i2)
{
	int i;

	for (i = 0; i < D; i++)
	{
		v[i] = v1[i] + v2[i][i2];
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void var_exp_1d_q(float *v1, float **v2, int D, float *v, int i2, int p, int q)
{
	int i;

	for (i = 0; i < D; i++)
	{

		v[i] = quantize(v1[i] + v2[i][i2], p, q);
	}

	return;
}

//////////////////////////////////////////////////////////////////////
void var_exp_2d(float *v1, float **v2, int D, float **v, int i2)
{
	int i;

	for (i = 0; i < D; i++)
	{
		v[i][i2] = v1[i] + v2[i][i2];
	}

	return;
}

//////////////////////////////////////////////////////////////////////
void var_exp_2d_q(float *v1, float **v2, int D, float **v, int i2, int p, int q)
{
	int i;

	for (i = 0; i < D; i++)
	{

		v[i][i2] = quantize(v1[i] + v2[i][i2], p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void fact_exp1(float **vin, int D, float Pd, float **vout, int iin, int iout)
{
	int i;
	float a, b;
	float dd, dd1, dd2, dd3, dd4, dd5, dd6;

	a = log((1.0 - Pd) / 2.0);
	b = log(Pd / 4.0);

	for (i = 0; i < D; i++)
	{
		dd1 = a + vin[i][iin];
		dd2 = a + vin[(int)mod(i + D / 2, D)][iin];
		dd3 = b + vin[(int)mod(i + 1, D)][iin];
		dd4 = b + vin[(int)mod(i + D - 1, D)][iin];
		dd5 = b + vin[(int)mod(i + D / 2 + 1, D)][iin];
		dd6 = b + vin[(int)mod(i + D / 2 - 1, D)][iin];
		dd = logsumexp(dd1, dd2);
		dd = logsumexp(dd, dd3);
		dd = logsumexp(dd, dd4);
		dd = logsumexp(dd, dd5);
		dd = logsumexp(dd, dd6);
		vout[i][iout] = dd;
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void fact_exp1_q(float **vin, int D, float Pd, float **vout, int iin, int iout,
				 int p, int q)
{
	int i;
	float a, b;
	float dd, dd1, dd2, dd3, dd4, dd5, dd6;
	int p_a, q_a;
	int p_b, q_b;

	p_a = 1;
	q_a = 5;
	p_b = 4;
	q_b = 5;

	a = quantize(log((1.0 - Pd) / 2.0), p_a, q_a);
	b = quantize(log(Pd / 4.0), p_b, q_b);

	for (i = 0; i < D; i++)
	{
		dd1 = a + vin[i][iin];
		dd2 = a + vin[(int)mod(i + D / 2, D)][iin];
		dd3 = b + vin[(int)mod(i + 1, D)][iin];
		dd4 = b + vin[(int)mod(i + D - 1, D)][iin];
		dd5 = b + vin[(int)mod(i + D / 2 + 1, D)][iin];
		dd6 = b + vin[(int)mod(i + D / 2 - 1, D)][iin];
		dd = logsumexp(dd1, dd2);
		dd = logsumexp(dd, dd3);
		dd = logsumexp(dd, dd4);
		dd = logsumexp(dd, dd5);
		dd = logsumexp(dd, dd6);
		vout[i][iout] = quantize(dd, p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void fact_exp2(float **vin, int D, float Pd, int alpha, float *vout, int ii)
{
	float a, b;
	float dd, dd1, dd2;
	int i, i1, i2, i3;

	a = log(1.0 - Pd);
	b = log(Pd / 2.0);

	for (i = 0; i < D; i++)
	{
		i1 = mod(i + alpha * D / 2, D);
		i2 = mod(i + alpha * D / 2 + 1, D);
		i3 = mod(i + alpha * D / 2 - 1 + D, D);
		dd1 = a + vin[i1][ii];
		dd2 = b + vin[i2][ii];
		dd = logsumexp(dd1, dd2);
		dd1 = dd;
		dd2 = b + vin[i3][ii];
		dd = logsumexp(dd1, dd2);
		vout[i] = dd;
	}

	return;
}

///////////////////////////////////////////////////////////////////////
void fact_exp2_q(float **vin, int D, float Pd, int alpha, float *vout, int ii,
				 int p, int q)
{
	float a, b;
	float dd, dd1, dd2;
	int i, i1, i2, i3;
	int p_a, q_a;
	int p_b, q_b;

	p_a = 1;
	q_a = 5;
	p_b = 4;
	q_b = 5;

	a = quantize(log(1.0 - Pd), p_a, q_a);
	b = quantize(log(Pd / 2.0), p_b, q_b);

	for (i = 0; i < D; i++)
	{
		i1 = mod(i + alpha * D / 2, D);
		i2 = mod(i + alpha * D / 2 + 1, D);
		i3 = mod(i + alpha * D / 2 - 1 + D, D);
		dd1 = a + vin[i1][ii];
		dd2 = b + vin[i2][ii];
		dd = logsumexp(dd1, dd2);
		dd1 = dd;
		dd2 = b + vin[i3][ii];
		dd = logsumexp(dd1, dd2);
		vout[i] = quantize(dd, p, q);
	}

	return;
}

///////////////////////////////////////////////////////////////////////
float logsumexp(float d1, float d2)
{
	float d;

	d = Maxf(d1, d2); // + log(1.0 + exp(-fabs(d2-d1)));

	return d;
}

///////////////////////////////////////////////////////////////////////
float quantize(float data_in, int p, int q)
{
	float data_out; // quantized output
	float data_tmp; // support variable
	int data_int;	// support variable

	// multiply by 2^q
	data_tmp = data_in * pow(2.0, q);

	// truncate to 0 bits after point
	data_int = (int)data_tmp;
	if (data_tmp > 0) // positive value are rounded up
	// when the error is equal or greater than 0.5
	{
		if (data_tmp - data_int >= 0.5) // >= to be in line with Sysgen rounding
		{
			data_int = data_tmp + 1;
		}
	}
	else // negative value are rounded down
	// when the error is equal or less than -0.5
	{
		if (data_tmp - data_int <= -0.5) // <= to be in line with Sysgen rounding
		{
			data_int = data_tmp - 1;
		}
	}

	// saturate at -2^(p+q-1) and 2^(p+q-1)-1
	if (data_int > (int)pow(2.0, (p + q - 1)) - 1)
	{
		data_int = (int)pow(2.0, (p + q - 1)) - 1;
	}
	if (data_int < -1 * (int)pow(2.0, (p + q - 1)))
	{
		data_int = -1 * (int)pow(2.0, (p + q - 1));
	}

	data_tmp = (float)data_int;

	// divide by 2^q
	data_out = data_tmp / pow(2.0, q);

	return data_out;
}

///////////////////////////////////////////////////////////////////////
float cordic_getphase(float re, float im, int iterations)
{
	float phase;
	float re_tmp;
	int i;

	if (re < 0)
	{
		re_tmp = re;
		if (im > 0)
		{
			// subtract 90°
			re = im;
			im = -re_tmp;
			phase = PI / 2;
		}
		else
		{
			// add 90°
			re = -im;
			im = re_tmp;
			phase = -PI / 2;
		}
	}
	else
	{
		phase = .0;
	}

	// CORDIC iterations
	for (i = 0; i < iterations; i++)
	{
		re_tmp = re;
		if (im > .0)
		{
			// phase is positive, do negative rotation
			re = re + im * pow(2.0, (-i));
			im = im - re_tmp * pow(2.0, (-i));
			phase = phase + atan(pow(2.0, (-i)));
		}
		else if (im < .0)
		{
			// phase is negative, do positive rotation
			re = re - im * pow(2.0, (-i));
			im = im + re_tmp * pow(2.0, (-i));
			phase = phase - atan(pow(2.0, (-i)));
		}
	}

	return phase;
}

///////////////////////////////////////////////////////////////////////
float cordic_getmag(float re, float im, int iterations)
{
	float magnitude;
	float re_tmp;
	float acc_phase;
	int i;

	if (re < 0)
	{
		re_tmp = re;
		if (im > 0)
		{
			// subtract 90°
			re = im;
			im = -re_tmp;
			acc_phase = PI / 2;
		}
		else
		{
			// add 90°
			re = -im;
			im = re_tmp;
			acc_phase = -PI / 2;
		}
	}
	else
	{
		acc_phase = .0;
	}

	// CORDIC iterations
	for (i = 0; i < iterations; i++)
	{
		re_tmp = re;
		if (im >= .0)
		{
			// phase is positive, do negative rotation
			re = re + im * pow(2.0, (-i));
			im = im - re_tmp * pow(2.0, (-i));
			acc_phase = acc_phase + atan(pow(2.0, (-i)));
		}
		else
		{
			// phase is negative, do positive rotation
			re = re - im * pow(2.0, (-i));
			im = im + re_tmp * pow(2.0, (-i));
			acc_phase = acc_phase - atan(pow(2.0, (-i)));
		}
	}

	magnitude = re / CORDIC_SCALINGFACTOR;

	return magnitude;
}

///////////////////////////////////////////////////////////////////////
double complex cordic_polartorect(float magnitude, float phase, int iterations)
{
	double complex data_out;
	float re, re_tmp, re_out;
	float im, im_out;
	float phase_int;
	int i;

	re = magnitude / CORDIC_SCALINGFACTOR;
	im = 0.0;

	// CORDIC algorithm only works if phase is included in
	// (-PI/2; PI/2) interval
	phase_int = phase;
	if (phase >= PI / 2)
	{
		phase_int = phase - PI / 2;
	}
	if (phase <= -PI / 2)
	{
		phase_int = phase + PI / 2;
	}

	// CORDIC iterations
	for (i = 0; i < iterations; i++)
	{
		re_tmp = re;
		if (phase_int < .0)
		{
			re = re + im * pow(2.0, (-i));
			im = im - re_tmp * pow(2.0, (-i));
			phase_int = phase_int + atan(pow(2.0, (-i)));
		}
		else
		{
			re = re - im * pow(2.0, (-i));
			im = im + re_tmp * pow(2.0, (-i));
			phase_int = phase_int - atan(pow(2.0, (-i)));
		}
	}

	// correct result if starting phase wasn't included in
	// (-PI/2; PI/2) interval
	re_out = re;
	im_out = im;
	if (phase >= PI / 2)
	{
		re_out = -im;
		im_out = re;
	}
	if (phase <= -PI / 2)
	{
		re_out = im;
		im_out = -re;
	}

	// compose output vector
	data_out = re_out + I * im_out;

	return data_out;
}

///////////////////////////////////////////////////////////////////////
float mod(float a, float b)
{
	return a - ((int)(a / b)) * b;
}

///////////////////////////////////////////////////////////////////////
float Minf(float a, float b)
{
	if (a < b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

///////////////////////////////////////////////////////////////////////
float Maxf(float a, float b)
{
	if (a > b)
	{
		return a;
	}
	else
	{
		return b;
	}
}
///////////////////////////////////////////////////////////////////////
float MinInt(int a, int b)
{
	if (a < b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

///////////////////////////////////////////////////////////////////////
float MaxInt(int a, int b)
{
	if (a > b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

void flip_bytes(int *bitout, const int data_len)
{
	if (bitout == NULL || data_len < 8)
	{
		return; // Early exit if input is invalid or too short
	}

	for (size_t ll = 0; ll < data_len - 7; ll += 8)
	{
		int tmp;

		tmp = bitout[ll];
		bitout[ll] = bitout[ll + 7];
		bitout[ll + 7] = tmp;

		tmp = bitout[ll + 1];
		bitout[ll + 1] = bitout[ll + 6];
		bitout[ll + 6] = tmp;

		tmp = bitout[ll + 2];
		bitout[ll + 2] = bitout[ll + 5];
		bitout[ll + 5] = tmp;

		tmp = bitout[ll + 3];
		bitout[ll + 3] = bitout[ll + 4];
		bitout[ll + 4] = tmp;
	}
}
