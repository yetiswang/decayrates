# Plasmon-enhanced single molecule fluorescence and single-molecule photodynamics calculation based on MNPBEM Matlab toolbox

## Introduction
This numerical calculation is based on [MNPBEM MATLAB toolbox](http://physik.uni-graz.at/mnpbem/). First time users are expected to read the manuals of MNPBEM for debugging and explanation of BEM commands. Surface plasmons are coherently oscillating eletrons in metal nanoparticles. Through coupling to plasmon resonance modes of nanoparticles, the fluorescence process of a single-molecule fluorohore is strongly modified. The current scripts calculate the modified excitation and emission rates as a function of position around a nanoparticle with a given size, shape and material, and also the photodynamics of the modified molecule based on a three-state model. 

## Getting started

Install MNPBEM matlab toolbox first and add all codes to your Matlab path.
Find MNPBEM here: http://physik.uni-graz.at/mnpbem/. Please read the developers' publications and Matlab help file in order to be able to work with the toolbox commands. 

## Steps to follow: 

1. Start by opening `main_FE_1D.m` and modifiy the `height`, `size`, `metal`, excitation `enei_field` and dye emission `enei_dipole` wavelengths in nanometer. The program is designed to calculated at one time the fluorescence enhancement with multiple sizes/shapes of particles, and dyes.  
2. Run `main_FE_1D.m` under desired directory where you want the results to be saved. The folders containing results will be named as `H(height)D(diameter)_(metal)_Exc(enei_field)_Dip_(enei_dipole)`. BEM results will be saved in `decayrates.mat`, and photon count rate will be saved in `PCR.mat`. 
3. Under the same directory, open `Plot_results.m` for reading and plotting the calculated results. 

## More tweaks 

You can find more parameters to vary in the core script `FE_GNR_1D` and `FE_GNS_1D`. More to follow ...

## Data structures 

The results are saved in Matlab structure arrays `struct`, and a list of fields in the results with notes are listed here. More to follow ... 
