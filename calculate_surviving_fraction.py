#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:52:45 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This code simulates the surviving fraction from LQ RT functions.
"""

# Import required libraries
import numpy as np
from get_paraview_data import *

# =============================================================================
# PARAMETERS
# =============================================================================

# Common parameters
K_OER = 3.28  # mmHg (lower creates pronounced drop)
# K_OER = 3  # mmHg (lower creates pronounced drop)
d = 2
radioresistant_threshold_mmHg = 3.28
# radioresistant_threshold_nM = 4116
radioresistant_threshold_nM = 4500
# n = 1

# Set the parameters for the Wouters model (Wouters and Brown, 1997) which should give SF of 47% under aerobic conditions
alpha_wb = 0.2 
alpha_beta_ratio_wb = 2.2
beta_wb = alpha_wb/alpha_beta_ratio_wb
OER_am = 2.5
OER_bm = 3.0

# Set the parameters for the Scott model
# alpha_max = 0.2895*3  # Scott code
alpha_max = 0.3  # Grogan
alpha_beta_ratio = 10 
# alpha_beta_ratio = 8.64  
beta_max = alpha_max/alpha_beta_ratio
OER_alpha_max = 1.75 
OER_beta_max = 3.25  
OER_min = 1  # lower creates bigger discrepancy

# Set the parameters for the Lewin model
alpha_lewin = 0.35
beta_lewin = alpha_lewin/alpha_beta_ratio
OER_hyp = 3 

# Set the parameters for a custom function
alpha_my = alpha_lewin
beta_my = alpha_lewin/alpha_beta_ratio_wb

# =============================================================================
# FUNCTIONS
# =============================================================================

# Compute the SF equation based on the basic parameters
def compute_SF(alpha, beta):
    return np.exp(-1*(alpha * d + beta * (d ** 2)))

# Compute the SF equation based on the Lewin parameters
def compute_Lewin_SF(norm_fraction):
    
    # Calculate the radioresistant fraction of cells
    radior_fraction = 1-norm_fraction
    
    # print(norm_fraction, radior_fraction, norm_fraction*compute_SF(alpha_max, beta_max), radior_fraction*compute_SF(alpha_hyp, beta_hyp))
    alpha_hyp, beta_hyp = compute_alpha_beta(alpha_lewin, beta_lewin, OER_hyp, OER_hyp)
    SF = norm_fraction*compute_SF(alpha_lewin, beta_lewin) + radior_fraction*compute_SF(alpha_hyp, beta_hyp)
    
    # Return the SF
    # return SF, OER_hyp
    return SF

# Define a Lewin function to compute the SF for a range of concentrations of c
def calculate_range_Lewin_SF(c):

    # Initialize an array to store the resulting SF values
    SF_array = np.zeros_like(c, dtype=float)

    # Iterate over the array and compute SF for each element
    for i, value in enumerate(c):
        # Set norm_fraction to 1 if value is above threshold, otherwise set it to 0
        if value > radioresistant_threshold_mmHg:
            norm_fraction = 1
        else:
            norm_fraction = 0

        # Calculate the SF using the compute_Lewin_SF function
        # SF_array[i], OER_hyp = compute_Lewin_SF(norm_fraction)
        SF_array[i] = compute_Lewin_SF(norm_fraction)
    
    # Return the SF and the OER in the same dimensions
    # return SF_array, np.full(SF_array.shape, OER_hyp)
    return SF_array

# Define a function to compute the SF for a cell population with an oxygen 
# concentration of c in each cell using the function from Lewin et al. (2018) 
# (in MV Chaste, death probability is calculated for each cell as 1-SF)
def calculate_population_Lewin_SF(c):
    
    # Compute the SFs for the two populations separately and then combine them
    SF = compute_Lewin_SF(compute_cells_NF(c,radioresistant_threshold_nM))
    
    return SF    

# Define a function to compute the SF for an array of concentrations of c (assuming c is specified in mmHg) using 
# the function from Scott et al. (2016) (in MV Chaste, death probability is 
# calculated for each cell as 1-SF)
def calculate_array_Scott_SF(c):
   
    OER_alpha = (((OER_alpha_max - OER_min) * K_OER) / (c + K_OER)) + OER_min
    OER_beta = (((OER_beta_max - OER_min) * K_OER) / (c + K_OER)) + OER_min
    alpha, beta = compute_alpha_beta(alpha_max, beta_max, OER_alpha, OER_beta)
    
    # Compute the surviving fraction
    # SF = np.exp(-n * (alpha * d + beta * (d ** 2)))
    SF = compute_SF(alpha, beta)

    # return SF, alpha, beta, OER_alpha, OER_beta
    # return SF, OER_alpha, OER_beta
    return SF

# Define a function to compute the alpha and beta values for the SF function
def compute_alpha_beta(alpha_num, beta_num, OER_alpha, OER_beta):
    
    alpha = alpha_num/OER_alpha
    beta = beta_num/np.square(OER_beta)
    
    return alpha, beta

# Define a function to compute the SF for an array of concentrations of c (assuming c is specified in mmHg) using 
# the function from Wouters and Brown (1997) (in MV Chaste, death probability is 
# calculated for each cell as 1-SF)
def calculate_array_Wouters_SF(c):

    OER_alpha = (c*OER_am + K_OER)/(c + K_OER)
    OER_beta = (c*OER_bm + K_OER)/(c + K_OER)
    alpha_h, beta_h = compute_alpha_beta(alpha_wb, beta_wb, OER_am, OER_bm)
    alpha_equivalent = alpha_h * OER_alpha
    beta_equivalent = beta_h * (OER_beta**2)     

    # Compute the surviving fraction
    SF = compute_SF(alpha_equivalent, beta_equivalent)
    
    # return SF, OER_alpha, OER_beta
    return SF

# Make a simple Michaelis–Menten equation
def calculate_array_my_SF(c):

    # SF_my = (((OER_hyp-OER_min)*c)/(K_OER+c)) + OER_min
    
    OER_my = OER_hyp-(((OER_hyp-OER_min)*c)/(K_OER+c))
    alpha, beta = compute_alpha_beta(alpha_my, beta_my, OER_my, OER_my)
    SF = compute_SF(alpha, beta)
    # print(OER_my)

    # return SF, OER_my
    return SF

# Make a function to compute the RR relative to the SF if there was no oxygen
def compute_RR(SF_array, anoxic_SF):

    # Compute the RR relative to the minimum sensitivity 
    min_death = 1-anoxic_SF
    actual_death = 1-SF_array
    rr = actual_death / min_death

    return rr

# Plot the functions
# from lqrt_calibration import *