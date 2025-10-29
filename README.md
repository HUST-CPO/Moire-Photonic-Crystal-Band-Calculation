# MPCBC: Moiré Photonic Crystal Band Calculation
![MoirePhC_schematic](/github_cover.png)

## Overview
The open-source MATLAB package MPCBC is a numerical program for calculating the physical properties of incommensurate Moiré photonic crystals, including the band structures, density of states, and eigenmodes. The program can compute 2D Moiré photonic crystals with arbitrary twist angles based on the improved plane-wave expansion method. 

If you find MPCBC useful for your research, please consider giving a star. Thank you.

## Usage
### Scripts
1. `phc_2D_moire_mynote.m`: The main program to calculate the band structure of 2D Moiré photonic crystals.
2. `phc_2D_moire_DOS.m`: Visualize the band structure. Calculate and visualize the density of states based on bands. Calculate and visualize the eigenmodes based on eigenvectors.
3. `phc_2D_square.m`: Calculate the band structure and DOS of 2D periodic photonic crystals using standard plane-wave expansion method.

### Functions
1. `S_fft.m`: Generate FFT of S (i.e., the distribution of nanodisks) analytically, which is used in Moiré photonic crystal band calculation.
2. `E_z.m`: Calculate the eigenmode $E_z$ for TM mode in plane-wave expansion.

## Run examples
You can follow the following steps to run the program:
1. Run the script `phc_2D_moire_mynote.m` and save the data in `phc_2D_moire_data.mat`.
2. Run the script `phc_2D_moire_DOS.m` to visualize the band structures, density of states, and eigenmodes.
3. Run the script `phc_2D_square.m` to compare with the results from standard plane-wave expansion methods for periodic case.
