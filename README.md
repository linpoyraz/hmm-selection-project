# hmm-selection-project

Final project for CS 4775 Fall 2022. Reimplementation project of: 

Estimating time-varying selection coefficients from time series data of allele frequencies
Iain Mathieson
bioRxiv 2020.11.17.387761; doi: https://doi.org/10.1101/2020.11.17.387761

For single-allele data simulation and analyses run:

` make single-allele`

For multi-allele data simulation and analyses run (takes ~3 hours for 4 cores):

` make multi-allele`

The SLiM simulations were run on SLiM GUI. The Eidos script to run these simulations is  `./slim_scenario3` which outputs segregating site files into `./slim_s2/`


Main reference: [https://github.com/mathii/slattice](https://github.com/mathii/slattice)
