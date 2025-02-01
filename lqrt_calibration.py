#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Sep  2 14:35:30 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This code simulates the surviving fraction from the LQ RT functions.
"""

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Import custom functions
from convert_oxygen_units import *
from calculate_surviving_fraction import *
from pathlib import Path

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 26})
# sns.set(font='STIXGeneral', font_scale=1.7)

# Generate values for c
c_values = np.linspace(0, 100, 10000)  # c specified in mmHg

# Plot the SFs
plt.figure(figsize=(10, 6))
SF_Scott = calculate_array_Scott_SF(c_values)
# SF_Scott, OER_alpha_scott, OER_beta_scott = calculate_array_Scott_SF(c_values)
plt.plot(c_values, SF_Scott, label='Scott')
# print('SF at average level is: ', calculate_range_Scott_SF(convert_nM_to_mmHg(21490)))
SF_Lewin = calculate_range_Lewin_SF(c_values)
# SF_Lewin, OER_lewin = calculate_range_Lewin_SF(c_values)
plt.plot(c_values, SF_Lewin, label='Lewin')
SF_Wouters = calculate_array_Wouters_SF(c_values)
# SF_Wouters, OER_alpha_wouters, OER_beta_wouters = calculate_array_Wouters_SF(c_values)
plt.plot(c_values, SF_Wouters, label='Wouters')
SF_my = calculate_array_my_SF(c_values)
# SF_my, OER_my = calculate_array_my_SF(c_values)
plt.plot(c_values, SF_my, label='my')
# calculate_array_my_SF
# plt.title(r'Surviving Fraction (SF) vs c (mmHg)')
plt.xlabel(r'c (mmHg)')
plt.ylabel(r'Surviving Fraction (SF)')
plt.xlim(0,20)
plt.ylim(0,1.25)
plt.axvline(3, ls='--', label='threshold', c='black')
plt.axvline(10, ls='--', label='mean tumour O2')
plt.grid(True)
plt.legend()
plt.show()

# Create figure of radiosensitivity
plt.figure(figsize=(10, 6))
scott_rr = 1-SF_Scott
plt.plot(c_values, scott_rr/scott_rr[0], label='Scott')
lewin_rr = 1-SF_Lewin
plt.plot(c_values, lewin_rr/lewin_rr[0], label='Lewin')
wouters_rr = 1-SF_Wouters
plt.plot(c_values, wouters_rr/wouters_rr[0], label='Wouters')
my_rr = 1-SF_my
plt.plot(c_values, my_rr/my_rr[0], label='my')
plt.axvline(3, ls='--', label='threshold')
plt.axvline(10, ls='--', label='mean tumour O2')
plt.xlabel(r'c (mmHg)')
plt.ylabel('relative radiosensitivity')
plt.ylim(0,3.25)
plt.grid(True)
plt.legend()
plt.show()
'''
# Create figure of OERs
plt.figure(figsize=(10, 6))
plt.plot(c_values, OER_alpha_scott, label='Scott alpha')
plt.plot(c_values, OER_beta_scott, label='Scott beta')
plt.plot(c_values, OER_alpha_wouters, label='Wouters alpha')
plt.plot(c_values, OER_beta_wouters, label='Wouters beta')
plt.plot(c_values, OER_lewin, label='Lewin')
plt.plot(c_values, OER_my, label='my')
plt.axvline(3, ls='--', label='threshold')
plt.axvline(10, ls='--', label='mean tumour O2')
plt.xlabel(r'c (mmHg)')
plt.ylabel('OER')
plt.ylim(0,3.25)
plt.grid(True)
plt.legend()
plt.show()
'''
# Sweep ranges for K_oER, OER_alpha_max, and OER_beta_max
# K_oER_values = np.linspace(1, 10, 100)  # K_oER ranging from 1 to 5
# OER_alpha_max_values = np.linspace(1, 10, 100)  # OER_alpha_max ranging from 1.5 to 2.5
# OER_beta_max_values = np.linspace(1, 10.0, 100)  # OER_beta_max ranging from 2.5 to 4.0

# Create a meshgrid for K_oER, OER_alpha_max, and OER_beta_max
# K_oER_grid, OER_alpha_max_grid, OER_beta_max_grid = np.meshgrid(K_oER_values, OER_alpha_max_values, OER_beta_max_values)

# Compute the SF function values for each combination of K_oER, OER_alpha_max, and OER_beta_max
# SF_values = calculate_SF(K_oER_grid, OER_alpha_max_grid, OER_beta_max_grid)

# Plotting the parameter sweep (we can choose a slice, e.g., fix one parameter)
# plt.figure(figsize=(12, 8))
# # Example: Fix OER_alpha_max at the midpoint of its range for visualization
# fixed_OER_alpha_max = OER_alpha_max_values[50]
# contour = plt.contourf(K_oER_grid[:, :, 50], OER_beta_max_grid[:, :, 50], SF_values[:, :, 50], 50, cmap='viridis')
# plt.colorbar(contour)
# # plt.title(r'Parameter Sweep of K and $OER_{\beta_{max}}$ for fixed $OER_{\alpha_{max}}$ = {:.2f}'.format(fixed_OER_alpha_max))
# plt.xlabel(r'$K_{oER}$')
# plt.ylabel(r'$OER_{\beta_{max}}$')
# plt.show()

# =============================================================================
# THESIS PLOT
# =============================================================================

# Plot the SFs
fig, ax1 = plt.subplots(figsize=(15, 10))
width = 5

# Primary y-axis for surviving fractions
SF_Scott = calculate_array_Scott_SF(c_values)
SF_Lewin = calculate_range_Lewin_SF(c_values)
SF_Wouters = calculate_array_Wouters_SF(c_values)
SF_my = calculate_array_my_SF(c_values)

# Shade the region below the first vertical line (c=3)
ax1.axvspan(0, 3.28, color='crimson', alpha=0.2, label='hypoxia')

lewin_line, = ax1.plot(c_values, SF_Lewin, label='discrete SF',  ls='--', lw=width, color='C0')
wouters_line, = ax1.plot(c_values, SF_Wouters, label='continuous SF', lw=width, color='C0')

ax1.set_xlabel(r'oxygen partial pressure (mmHg)')
ax1.set_ylabel(r'surviving fraction (SF)')
ax1.set_xlim(0, 15)
# ax1.axvline(3.28, ls=':', c='red', lw=width)
# ax1.axvline(10, ls=':', c='green', lw=width)
ax1.set_ylim(0, 1)
ax1.grid(True, which='both', axis='x')

# Secondary y-axis for relative ratios
ax2 = ax1.twinx()

# Compute and plot the relative radiosensitivity 
lewin_rr = compute_RR(SF_Lewin, anoxic_SF=calculate_range_Lewin_SF([0]))
wouters_rr = compute_RR(SF_Wouters, anoxic_SF=calculate_array_Wouters_SF(0))
ax2.plot(c_values, lewin_rr / lewin_rr[0], label='discrete RR', ls=':', lw=width, color='crimson')
ax2.plot(c_values, wouters_rr / wouters_rr[0], label='continuous RR', ls='dashdot', lw=width, color='crimson')
ax2.set_ylim(1, 3)
ax2.set_ylabel(r'relative radiosensitivity (RR)')

# Combine legends from both axes
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='lower right')

# Save image
main_folder_path = '/Users/vedang/Documents/GitHub/dphil-scripts/demo_datasets/'
figure_folder = main_folder_path
figure_title = 'lqrt_functions'
file_path = Path(figure_folder + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500, bbox_inches = 'tight')
plt.show()
