# Generic HLS FFT

## Description
This repository hosts a Generic HLS-ready FFT model. Several parameters are available :
- Input precision format
- Internal precision format
- Twiddle factor precision format
- FFT length
- Butterfly parallelization (streaming width)
- Stage pipelining (by stage groups)
- Radix selection

This repository corresponds to an article : *publication pending*.

## Using the model
Generating LUTs and tools :
```sh
[repo-dir]$ mkdir build
[repo-dir]$ cd build
[repo-dir/build]$ cmake .. #(or run "ccmake .." to configure AP_TYPES_LIB_DIR and other things)
[repo-dir/build]$ make
```
A tool allowing to perform a C/C++ simulation of the FFT model is located at *repo-dir/build/bin/FFT_calc* : 
```sh
usage : FFT_calc --in-signal <input_file.txt> --out-signal <output_file.txt> --constants <twiddle_factor_directory> --fft-model <sbg-radix|fftw-gold> --sw <2|4|8> --radix <2|4|8> --input-quantif <INPUT_WIDTH> --internal-quantif <OUTPUT_WIDTH> --rom-quantif <TWIDDLE_WIDTH>
```
Some input signals are generated when compiling and are available at *repo-dir/signals*

FFT Model sources are available at *repo-dir/src/fft/model*

## Dependencies
> - cmake >= 3.16.3
> - fftw3 (Ubuntu 18.04+ : `apt install libfftw-dev`)
> - tcl-tdom (Ubuntu 18.04+ : `apt install tcl tdom`)
> - arbitrary precision library ([**DrasLorus HLS_arbitrary_Precision_Types**](https://github.com/DrasLorus/HLS_arbitrary_Precision_Types) or [Xilinx HLS_arbitrary_Precision_Types](https://github.com/Xilinx/HLS_arbitrary_Precision_Types) or Xilinx library shipped with Vivado). Library directory path can be set with *ccmake* with **AP_TYPES_LIB_DIR** (defaulting to `/opt/Xilinx/Vivado/2020.2/include/`).

## Known Issues
- Generated architecture is known to use a big amount of BRAM
- Internal format is constant across stages and thus should be log_R(N) bigger than input format. It is also recommanded for now to use the same type for input format (and not to fill all its bits) and internal format.

## To do list
- Fixing the vector format issues in fixed point mode is considered.

## Licensing
This project is released under GNU GPLv3 License. It also contains modified parts of the SPIRAL project under its SPIRAL License. Please refer to the file header for more information.