## ESA AIS Receiver - Executable

Performs the demodulation and detection for multiple AIS channels in parallel.
Initially developed under TEC-ESC and adapted to the Sentinel-1 data.
Ref: G. Colavolpe, T. Foggi, A. Ugolini, J. Lizarraga, S. Cioni and A. Ginesi, "A highly efficient receiver for satellite-based Automatic Identification System signal detection," 2014

---

### Stand alone execution:
```
./AIS_receiver <output_directory> <input .asc file1>...<input .asc fileN> <data_len1>....<data_lenN>
```

#### **Command-line Arguments**

* `<output_directory>` — Output directory where the resulting `.txt` files will be saved.
* `<input .asc files>` — Text files containing the **RAW (real, imag)** AIS data sampled at `9600 × 3`.
* `<data_lens>` — Expected bit length for the channel to process:

  * `168` → for heritage AIS channels
  * `96` → for SAT AIS channels

#### **Files Loaded by the Executable**

* `<h1_file>` — Low-pass filter description file.
* `<h2_file>` — Matched-filter description file (oversampled filter matched to the principal pulse of the Laurent decomposition).
* `<h3_file>` — Low-pass filter (zonal) description file.
* `<h4_file>` — Low-pass filter (mod. Meng. Mor. frequency estimation) description file.
* `<pulse_file>` — Pulse description file.
* `<error_files>` — Error correction files.


Note that the program must be executed in the source directory (i.e. `./AIS_receiver`) in order to correctly read the error files.
The source code is located in the folder `/src/ESA_AIS_receiver_code/`

---

### Testcase:

run,
```
./AIS_receiver ./outputdir /ch0_testcase.asc 168
```

Result: Should take around 34 seconds and generate the file ch0_testcase_result.txt with the binary message detections.

---

### Compiling

To compile on MacOS, run:
```
make
```

To compile on Windows use the MSYS Mingw64 terminal:

Install the compiler
```
pacman -S mingw-w64-x86_64-gcc
```
Install make:
```
pacman -S make
```

Ensure that `libwinpthread-1.dll` is installed in the path or just copy it into the directory where it will be compiling.
Then from the ESA_AIS_receiver_code directory run:
make
Using the windows Makefile_win file
Some warnings will occur, but these are expected.

---

