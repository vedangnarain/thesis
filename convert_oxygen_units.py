#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jun 30 13:51:22 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This code converts oxygen concentration units (nM) to partial pressure units (mmHg).

Parameters and methods sourced from: 
    
[1] D. R. Grimes et al., ‘Estimating oxygen distribution from vasculature in 
three-dimensional tumour tissue’, J. R. Soc. Interface., vol. 13, no. 116, p. 
20160070, Mar. 2014, doi: 10.1098/rsif.2016.0070.

[2] T. D. Lewin, P. K. Maini, E. G. Moros, H. Enderling, and H. M. Byrne, ‘The 
Evolution of Tumour Composition During Fractionated Radiotherapy: Implications 
for Outcome’, Bull Math Biol, vol. 80, no. 5, pp. 1207–1235, May 2018, doi: 
10.1007/s11538-018-0391-9.
"""

# Define function to convert nM to mmHg
def convert_nM_to_mmHg(nM, K=22779):
    
    # Convert nM (moles of O2 per cubic decimetre) into kg per cubic metre
    concentration = nM*0.000000001*32
    
    # Calculate partial pressure (mmHg) using formula from [2]
    mmHg = K*concentration
    
    # Return partial pressure
    return mmHg
 
# Define function to mimic mmHg to nM
def convert_mmHg_to_nM(mmHg, K=22779):
    
    # Convert partial pressure (mmHg) to concentration (kg per cubic metre)
    concentration = mmHg/K
    
    # Convert kg per cubic metre into nM
    nM = concentration/(0.000000001*32)
    
    # Return partial pressure
    return nM    
