#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 15:43:38 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Here, we plot the grid refinement tradeoff.
"""

# Initialise libraries
import pandas as pd
import matplotlib.pyplot as plt# import os
import numpy as np
import os
from pathlib import Path

# Import custom functions
from get_paraview_data import *
from convert_oxygen_units import *
from calculate_surviving_fraction import *

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})
# sns.set(font='STIXGeneral', font_scale=1.7)

# Define a function to return the oxygen field
def get_oxygen_field(folder_path):

    field_data, field_spacing = get_vti_data(folder_path + '/oxygen_solution_0.vti')  # import the data for the network
    flat_field = field_data['oxygen'].to_numpy().flatten()
    
    return flat_field

# Define a function to get the normoxic fraction from the oxygen field of a simulation
def get_field_NF(folder_path):
    
    flat_field = get_oxygen_field(folder_path)
    # number_of_points = 1
    # for dim in np.shape(flat_field): number_of_points *= dim
    normoxic_points = (flat_field > 4116).sum()
    field_NF = normoxic_points/len(flat_field)  # calculate the hypoxic fraction

    return field_NF, flat_field   

# Grid sizes in micrometers
grid_sizes = np.array([2, 5, 10, 20, 30, 50, 100])

# Simulation results: each entry in the list represents the simulated values for a different grid size
# Let's say these are 1D arrays of the computed concentrations or diffusion values
# o2_1 = get_oxygen_field('/home/narain/Desktop/Datasets/Grid Spacing/1_um')
suffix = '/TestMantegazzaNetwork/ConstantHaematocrit'
o2_1_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/1_um'+suffix)
o2_2_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/2_um'+suffix)
o2_5_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/5_um'+suffix)
o2_10_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/10_um'+suffix)
o2_20_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/20_um'+suffix)
o2_30_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/30_um'+suffix)
# o2_40 = get_oxygen_field('/home/narain/Desktop/Datasets/Grid Spacing/40_um'+suffix)
o2_50_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/50_um'+suffix)
o2_100_NF, _ = get_field_NF('/home/narain/Desktop/Datasets/Grid Spacing/100_um'+suffix)

simulation_results = [o2_2_NF, o2_5_NF, o2_10_NF, o2_20_NF, o2_30_NF, o2_50_NF, o2_100_NF]
# simulation_results = [o2_10_NF,o2_20_NF,o2_50_NF,o2_100_NF]

# Reference solution
reference_solution = o2_1_NF  

# Computational cost 
computational_cost = np.array([25.25, 5.06, 1.85, 1.48, 1.42, 1.42, 1.03]) 
'''
# Function to calculate accuracy metric (e.g., Mean Squared Error)
def calculate_mse(simulation, reference):
    # Interpolate reference solution to match the size of the simulation result
    interpolated_reference = np.interp(
        np.linspace(0, 1, len(simulation)),
        np.linspace(0, 1, len(reference)),
        reference
    )
    return np.mean((simulation - interpolated_reference) ** 2)
# '''
# Function to calculate accuracy metric 
def calculate_mse(simulation, reference):
    return np.mean((simulation - reference) ** 2)
# '''
# Calculate accuracy (Mean Squared Error) for each grid size
mse_values = [calculate_mse(result, reference_solution) for result in simulation_results]
mse_values = [value * 100000 for value in mse_values]

'''
# Plot Accuracy vs Grid Size
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(grid_sizes, mse_values, marker='o', linestyle='-', color='b')
plt.xlabel('Grid Size (micrometers)')
plt.ylabel('Mean Squared Error')
plt.title('Accuracy vs Grid Size')
plt.xscale('log')
plt.yscale('log')
# plt.gca().invert_xaxis()  # Invert x-axis to show decreasing grid size from left to right
plt.grid(True, which="both", linestyle='--')

# Plot Computational Cost vs Grid Size
plt.subplot(1, 2, 2)
plt.plot(grid_sizes, computational_cost, marker='o', linestyle='-', color='r')
plt.xlabel('Grid Size (micrometers)')
plt.ylabel('Computational Cost (arbitrary units)')
plt.title('Computational Cost vs Grid Size')
plt.xscale('log')
# plt.gca().invert_xaxis()  # Invert x-axis to show decreasing grid size from left to right
plt.grid(True, which="both", linestyle='--')

plt.tight_layout()
plt.show()
'''
# Combined Plot: Accuracy and Computational Cost vs Grid Size
fig, ax1 = plt.subplots(figsize=(10, 5))

color = 'tab:blue'
ax1.set_xlabel('grid size (μm)')
ax1.set_ylabel(r'MSE $(\times 10^{-5})$')
line1 = ax1.plot(grid_sizes, mse_values, marker='o', linestyle='-', color=color, label='MSE')
ax1.set_xscale('log')
# ax1.set_yscale('log')
# ax1.invert_xaxis()  # Invert x-axis to show decreasing grid size from left to right
ax1.tick_params(axis='y')
ax1.grid(True, which="both", linestyle='--')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:red'
ax2.set_ylabel('execution time (seconds)')
line2 = ax2.plot(grid_sizes, computational_cost, marker='o', linestyle='--', color=color, label='execution time')
ax2.set_xscale('log')
ax2.tick_params(axis='y')
# Combine legends
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper center')

fig.tight_layout()  # to prevent overlap
# plt.title('Accuracy and Computational Cost vs Grid Size')

figure_title = 'grid_spacing'
# file_path = Path(figure_folder + figure_title + '.svg').expanduser()
# plt.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('/home/narain/Desktop/Final Figures/' + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500, bbox_inches = 'tight')

plt.show()