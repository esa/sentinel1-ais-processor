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
// File: AIS_receiver.c
// Author: Stefan Graham
/////////////////////////////////////////////////////////////////////////
// DESCRIPTION:
// Performs the demodulation and detection for a single channel raw AIS data as described in the patent:
// PCT/EP2014/051273
//   RECEIVING METHOD AND RECEIVER FOR TIMING AND FREQUENCY OFFSET CORRECTION
//   OF CONTINUOUS PHASE DEMODULATION IN SATELLITE-BASED AUTOMATIC
//   IDENTIFICATION SYSTEMS
// Use:
//   <executable> <data_len> <outputFile> <inputWavfile>
//   <data_len>  168 for heritage AIS Channels and 96 for SAT AIS Channels
//	 <outputFile> Output directory and txt file to save AIS detections
//	 <inputWavfile> Two channel Wav file containing the RAW complex AIS data sampled at 9600*3
//
// Files loaded by executable:
//   <h1_file> low-pass filter description file
//   <h2_file> matched-filter description file
//   <h3_file> low-pass filter(zonal) description file
//   <h4_file> low-pass filter(mod.Meng.Mor. frequency estim.) description file
//   <pulse_file> pulse description file
//   <error_files> error correction files
// Note that the program must be executed in the source directory (i.e. "./AIS_receiver"") in order to correctly read the error files.
/////////////////////////////////////////////////////////////////////////

#include "AIS_subroutines.h" //#CGS_2.4
#include <time.h>

// Define the WAV file header structure
typedef struct
{
	char chunkId[4];
	uint32_t chunkSize;
	char format[4];
	char subchunk1Id[4];
	uint32_t subchunk1Size;
	uint16_t audioFormat;
	uint16_t numChannels;
	uint32_t sampleRate;
	uint32_t byteRate;
	uint16_t blockAlign;
	uint16_t bitsPerSample;
	char subchunk2Id[4];
	uint32_t subchunk2Size;
} WavHeader;

int main(int argc, char **argv)
{
	const int NcR = 3; // Samples pr symbol
	const int data_len = atoi(argv[1]);
	const int fcs_len = 16;														   // FCS length
	const int training_len = 24;												   // training sequence length
	const int flag_len = 8;														   // start/stop flag length
	const int frame_len = data_len + fcs_len + flag_len + flag_len + training_len; // frame length = training sequence + start + data + FCS + stop
	const int Ta = 8;															   // 6; // ramp-up duration (# bits) //#CGS_2.4
	const int Tf = 8;															   // 8; // ramp-down duration, part 1 (# bits)
	const int Tg = 0;															   // 17; // ramp-down duration, part 2 (# bits)
	const int max_delay = 86;													   // maximum delay //#CGS_2.4
	const int max_stuffing = 4;													   // maximum number of stuffing bits

	int overlap = 240 * NcR;						  // Number of sampels in each sequence to attempt demodulation
	const int nr_samples = 90 * NcR;				  // Number of new samples added after each iteration
	const int pck_len = (overlap + nr_samples) / NcR; // AIS package length
	const int Nsimb = 2;
	const int lung = 6; // matched filter length, expressed as number of symbol periods
	// const int M = 2;	 // constellation order
	// const int num_h = 1; // modulation index's numerator
	// const int p = 2;	 // modulation index's denominator
	const int L = 3; // frequency pulse length

	int nMessagePosition;
	double complex cCancelSingleMsg[pck_len * NcR];
	double complex cCancelSingleMsgNoZonal[pck_len * NcR];

	// Initialize variables for syndrome error correction
	int err_sing_sin[4], err_sing[4];
	int err_doppi_sin[182], err_sing2c_sin[4];
	int err_doppi[182 * 2], err_sing2c[4];
	int err2c1_sin[8648], err2c1[8648 * 4], err3[2 * 3], err3_sin[2];
	int err2c2_sin[2862], err2c2[2862 * 8];
	int err2c3_sin[477], err2c3[477 * 12];
	int err2c4_sin[167], err2c4[167 * 16];
	memset(err_sing_sin, 0, sizeof(err_sing_sin));
	memset(err_sing, 0, sizeof(err_sing));
	memset(err_doppi_sin, 0, sizeof(err_doppi_sin));
	memset(err_sing2c_sin, 0, sizeof(err_sing2c_sin));
	memset(err_doppi, 0, sizeof(err_doppi));
	memset(err_sing2c, 0, sizeof(err_sing2c));
	memset(err2c1_sin, 0, sizeof(err2c1_sin));
	memset(err2c1, 0, sizeof(err2c1));
	memset(err3, 0, sizeof(err3));
	memset(err3_sin, 0, sizeof(err3_sin));
	memset(err2c2_sin, 0, sizeof(err2c2_sin));
	memset(err2c2, 0, sizeof(err2c2));
	memset(err2c3_sin, 0, sizeof(err2c3_sin));
	memset(err2c3, 0, sizeof(err2c3));
	memset(err2c4_sin, 0, sizeof(err2c4_sin));
	memset(err2c4, 0, sizeof(err2c4));

	int training_seq[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}; // training sequence definition
	float h1[Nsimb * NcR + 1];																	   // low-pass filter response
	float h2[lung * NcR];																		   // matched-filter response	#TF_2.3
	float h3[Nsimb * NcR + 1];																	   // low-pass filter response	#TF_2.3
	float h4[Nsimb * NcR + 1];																	   // low-pass filter response	#TF_2.3
	int frame[frame_len + max_stuffing];														   // transmitted frame
	double complex datain_zonal[pck_len * NcR];													   // transmitted vector zonal  //#CGS_2.4
	double complex datato_zonal[pck_len * NcR];													   // transmitted vector zonal  //#CGS_2.4
	double complex datain[pck_len * NcR];														   // transmitted vector //#CGS_2.4
	float pulse[L * NcR + 1];																	   // CPM pulse
	int bitout[pck_len];																		   // received frame
	float fout[0];																				   // Output frequency estimate

	// Support variables
	FILE *fd;
	FILE *databit;
	float f1;
	float f2;
	int ret_eof; // return_end of file
	int k, ll, i, res, detection_time;
	int count = 0;
	int iteration_nr = 0;
	double sample_time = 0;

	// WAV reading
	WavHeader header;
	int16_t buffer1[overlap * 2];
	int nr_samples_wav = nr_samples * 2; // 2*nr_samples since the Wav file is two channel interleaved.
	int16_t buffer[nr_samples_wav];

	// Check inputs
	if (argc < 3)
	{
		printf("Not enough parameters\n");
		printf("Use\n");
		printf("AIS_receiver <data_len> <outputFile> <inputWavfile>\n");
		exit(0);
	}

	if (data_len != 168 && data_len != 96)
	{
		printf("AIS data length must be 168 (Msg Type 1-3) or 96 (Msg type 27)\n");
		exit(0);
	}

	// Frequency shift of each zonal demodulator
	float fZonalValue; // #CGS_2.4
	float ZonalArray[] = {-0.33, 0, 0.33};
	// float ZonalArray[] = {-0.42, 0, 0.42};
	int nrZonals = sizeof(ZonalArray) / sizeof(ZonalArray[0]);

	// Open low-pass filter file (h1.dat)
	fd = fopen("h1.dat", "r");
	if (fd == NULL)
	{
		printf("doesn't exist\n");
		exit(0);
	}
	// Build low-pass filter impulse response
	for (i = 0; i < Nsimb * NcR + 1; i++)
	{
		fscanf(fd, "%f", &f1);
		h1[i] = f1;
	}
	fclose(fd);

	// Open matching filter file
	fd = fopen("h2.dat", "r");
	if (fd == NULL)
	{
		printf("doesn't exist\n");
		exit(0);
	}
	// Build matching filter impulse response
	for (i = 0; i < lung * NcR; i++)
	{
		fscanf(fd, "%f", &f1);
		h2[i] = f1;
	}
	fclose(fd);

	// Open low-pass filter impulse response
	fd = fopen("h3.dat", "r");
	if (fd == NULL)
	{
		printf("doesn't exist\n");
		exit(0);
	}
	// Build low-pass filter impulse response
	for (i = 0; i < Nsimb * NcR + 1; i++)
	{
		fscanf(fd, "%f", &f1);
		h3[i] = f1;
	}
	fclose(fd);

	// Open low-pass filter impulse response
	fd = fopen("h4.dat", "r");
	if (fd == NULL)
	{
		printf("doesn't exist\n");
		exit(0);
	}
	// Build low-pass filter impulse response
	for (i = 0; i < Nsimb * NcR + 1; i++)
	{
		fscanf(fd, "%f", &f1);
		h4[i] = f1;
	}
	fclose(fd);

	// Open modulation pulse file
	fd = fopen("pulse.dat", "r");
	if (fd == NULL)
	{
		printf("doesn't exist\n");
		exit(0);
	}
	// build pulse vector
	for (i = 0; i < L * NcR + 1; i++)
	{
		fscanf(fd, "%f", &f1);
		pulse[i] = f1;
	}
	fclose(fd);

	// open transmitted signal file
	fd = fopen(argv[3], "rb");
	if (fd == NULL)
	{
		printf("%s doesn't exist\n", argv[3]);
		exit(0);
	}

	// Read the WAV file header
	fread(&header, sizeof(WavHeader), 1, fd);
	// Check number of channels
	if (header.numChannels != 2)
	{
		printf("The WAV file must be two channels.\n");
		fclose(fd);
		return 1;
	}
	// Initialize datain
	fread(buffer1, overlap * 2, header.numChannels, fd);
	for (i = 0; i < overlap; i++)
	{
		datain[i] = buffer1[i * header.numChannels] + I * buffer1[i * header.numChannels + 1];
	}

	// File to save correctly decoded bits
	databit = fopen(argv[2], "a");
	if (databit == NULL)
	{
		printf("Failed to open output file\n");
		exit(0);
	}

	clock_t begin = clock();
	while (fread(buffer, nr_samples_wav, header.numChannels, fd) == 2)
	{
		for (i = 0; i < nr_samples; i++)
		{
			datain[i + overlap] = buffer[i * header.numChannels] + I * buffer[i * header.numChannels + 1];
		}

		// Attempt demodulation of AIS message for each zonal section
		for (k = 0; k < nrZonals; k++)
		{
			// Add zonal
			fZonalValue = ZonalArray[k];
			ZonalFloatingPoint(datain, datain_zonal, pck_len, fZonalValue, NcR, ZONAL_ADD);

			// Receiver procedure
			res = ais_receiver(datain_zonal, bitout, h1, h2, h3, h4, pulse, NcR,
							   Ta, Tf, Tg, max_delay, max_stuffing, pck_len, frame_len,
							   data_len, Nsimb, lung, L, frame, &nMessagePosition, cCancelSingleMsg,
							   err_sing, err_sing_sin, err_doppi, err_doppi_sin, err_sing2c,
							   err_sing2c_sin, err2c1, err2c1_sin, err3, err3_sin, err2c2,
							   err2c2_sin, err2c3, err2c3_sin, err2c4, err2c4_sin, fout);

			// #Save Output message detected
			if (res == 1)
			{
				// Discard detected preamble
				for (ll = 0; ll < data_len; ll++)
				{
					bitout[ll] = bitout[training_len + flag_len + ll];
				}

				// Flip bytes
				flip_bytes(bitout, data_len);

				// Save to csv: Doppler frequency, detection time[s], binary message
				detection_time = sample_time + (double)(NcR * nMessagePosition);

				fprintf(databit, "%.1lf, %d, ", (fout[0] + fZonalValue) * 9600, (int)detection_time);
				for (i = 0; i < data_len; i++)
				{
					fputc(bitout[i] + '0', databit);
				}
				fprintf(databit, "\n");
				count = count + 1;

				// Remove Zonal, cancel message and re-detect (Disabled since doesn't improve performance)
				// ZonalFloatingPoint(cCancelSingleMsg, cCancelSingleMsgNoZonal, pck_len, fZonalValue, NcR, ZONAL_REMOVE);
				// Remove message from input data
				// for (i = 0; i < pck_len * NcR; i++)
				//{
				//	datain[i] -= cCancelSingleMsgNoZonal[i];
				//}
			}
		}

		// Initialize next burst with large overlap between previous.
		for (i = 0; i < overlap; i++)
		{
			datain[i] = datain[i + nr_samples];
		}

		sample_time = sample_time + (double)nr_samples;
		iteration_nr += 1;

		// Termination condition
		if (iteration_nr > 10E10)
		{
			break;
		}
	}

	fclose(fd);
	fclose(databit);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%i Detections in %f (sec)\n", count, time_spent);
	return 0;
}