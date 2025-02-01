#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:54:28 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Read and plot the data from the Merlo experiments.

Tested in Python 3.11
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import time
from pathlib import Path
# from matplotlib.pyplot import figure
import matplotlib.colors as mcolors
# from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
# import scipy.stats
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error, r2_score

# Import tools for Paraview data
import sys
import os
parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_directory)
from get_paraview_data import *

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 30})

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read the data file
def read_amira_ascii(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Dictionaries to hold various data blocks by their @ index
    data_sections = {
        '1': [],   # VertexCoordinates
        '2': [],   # EdgeConnectivity
        '3': [],   # NumEdgePoints
        '4': [],   # EdgePointCoordinates
        '5': [],   # Thickness
        '6': [],   # Tube_Hematocrit_exp
        '7': [],   # velocityRBC_exp
        '8': [],   # Tube_Hematocrit_sim
        '9': [],   # velocityRBC_sim
        '10': []   # Tube_Hematocrit_sim_corr
    }

    current_section = None

    # Process each line to identify sections and parse data
    for line in lines:
        line = line.strip()
        
        # Skip empty lines and comments
        if line.startswith('#') or not line:
            continue

        # Check if line starts a new section with @ identifier
        if line.startswith('@'):
            current_section = line.strip('@')
            continue

        # Collect data for the current section if it matches one of our keys
        if current_section in data_sections:
            data_sections[current_section].extend(map(float, line.split()))

    # Parse and reshape each section based on the README specifications
    vertices = np.array(data_sections['1']).reshape((-1, 3))          # VERTEX { float[3] }
    edges = np.array(data_sections['2']).reshape((-1, 2)).astype(int)  # EDGE { int[2] }
    num_edge_points = np.array(data_sections['3']).astype(int)         # EDGE { int }
    edge_points = np.array(data_sections['4']).reshape((-1, 3))        # POINT { float[3] }
    thickness = np.array(data_sections['5'])                           # POINT { float }

    # Supplementary fields
    tube_hematocrit_exp = np.array(data_sections['6'])                 # EDGE { float } experimentally measured tube hematocrit
    velocity_rbc_exp = np.array(data_sections['7'])                    # EDGE { float } experimentally measured RBC velocity
    tube_hematocrit_sim = np.array(data_sections['8'])                 # EDGE { float } simulated tube hematocrit (initial uncorrected)
    velocity_rbc_sim = np.array(data_sections['9'])                    # EDGE { float } simulated RBC velocity
    tube_hematocrit_sim_corr = np.array(data_sections['10'])           # EDGE { float } corrected simulated tube hematocrit

    # Summary information for understanding data structure
    summary = {
        'N_v': vertices.shape[0],                   # Number of vertices
        'N_e': edges.shape[0],                      # Number of edges
        'N_p': edge_points.shape[0],                # Total number of discretization points
        'num_edge_points_per_edge': num_edge_points # List of discretization points per edge
    }

    # Store parsed data in a dictionary for easy access
    parsed_data = {
        'vertices': vertices,               # Vertex coordinates (micrometers)
        'edges': edges,                     # Edge connectivity (start and end vertex index)
        'num_edge_points': num_edge_points, # Number of points per edge (for visualization)
        'edge_points': edge_points,         # Coordinates of points along each edge (micrometers)
        'thickness': thickness,             # Diameter at each point (micrometers)
        'tube_hematocrit_exp': tube_hematocrit_exp,           # Experimentally measured tube hematocrit
        'velocity_rbc_exp': velocity_rbc_exp,                 # Experimentally measured RBC velocity
        'tube_hematocrit_sim': tube_hematocrit_sim,           # Simulated tube hematocrit (uncorrected)
        'velocity_rbc_sim': velocity_rbc_sim,                 # Simulated RBC velocity
        'tube_hematocrit_sim_corr': tube_hematocrit_sim_corr  # Corrected simulated tube hematocrit
    }

    return summary, parsed_data

# Define a function to sort the edges and plot the network
def plot_network(vertices, edges, tube_hematocrit_exp):
    # Normalize the hematocrit values for coloring
    norm = mcolors.Normalize(vmin=min(tube_hematocrit_exp), vmax=max(tube_hematocrit_exp))
    cmap = plt.cm.coolwarm  # Use the 'coolwarm' color map

    plt.figure(figsize=(10, 10))  # Set figure to be square

    # Sort edges based on custom logic: minimum x, then y, then secondary x and y
    def edge_key(edge):
        start, end = edge
        v1 = vertices[start]
        v2 = vertices[end]
        
        # Find the lower point in each dimension
        min_x = min(v1[0], v2[0])
        min_y = min(v1[1], v2[1])
        sec_x = max(v1[0], v2[0])
        sec_y = max(v1[1], v2[1])
        
        return (min_x, min_y, sec_x, sec_y)

    sorted_edges = sorted(enumerate(edges), key=lambda e: edge_key(e[1]))
    
    # Create a list to store the hematocrit values in the new order
    sorted_hematocrit_values = []
    
    # Plot each edge with color based on tube_hematocrit_exp and annotate with new label
    for new_label, (original_index, (start, end)) in enumerate(sorted_edges):
        x_values = [vertices[start][0], vertices[end][0]]
        y_values = [vertices[start][1], vertices[end][1]]
        
        # Get color for this edge
        color = cmap(norm(tube_hematocrit_exp[original_index]))
        
        # Plot the edge with the corresponding color
        plt.plot(x_values, y_values, color=color, linewidth=2)
        
        # Calculate midpoint for annotation
        midpoint_x = (x_values[0] + x_values[1]) / 2
        midpoint_y = (y_values[0] + y_values[1]) / 2
        
        # Annotate new edge label and hematocrit value at the midpoint
        # annotation_text = f"{new_label}: {tube_hematocrit_exp[original_index]:.3f}"
        # plt.text(midpoint_x, midpoint_y, annotation_text, color="black", fontsize=8, ha="center", va="center")
        
        # Store the hematocrit value in the new order
        sorted_hematocrit_values.append(tube_hematocrit_exp[original_index])
    
    # Calculate the total plot range to enforce identical x and y axis limits
    x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
    
    # Determine the center and maximum range to ensure the plot is square
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    max_range = max(x_max - x_min, y_max - y_min) / 2
    
    # Set limits for both axes to be the same, centered on the midpoint
    # plt.xlim(-15,600)
    # plt.ylim(-15,600)
    
    # Set aspect to equal for a square plot
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Add a color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(tube_hematocrit_exp)  # Setting the array for color mapping
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('field')
    plt.grid()
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.title('Plot of network using vertices from @1 and edges from @2')
    plt.show()
    
    return sorted_hematocrit_values

# Display summary and parsed data shapes to verify
'''
print("Summary Information:", summary_info)
print("Vertices Shape:", parsed_data['vertices'].shape)
print("Edges Shape:", parsed_data['edges'].shape)
print("NumEdgePoints Shape:", parsed_data['num_edge_points'].shape)
print("EdgePoints Shape:", parsed_data['edge_points'].shape)
print("Thickness Shape:", parsed_data['thickness'].shape)
print("Tube Hematocrit Exp Shape:", parsed_data['tube_hematocrit_exp'].shape)
print("Velocity RBC Exp Shape:", parsed_data['velocity_rbc_exp'].shape)
print("Tube Hematocrit Sim Shape:", parsed_data['tube_hematocrit_sim'].shape)
print("Velocity RBC Sim Shape:", parsed_data['velocity_rbc_sim'].shape)
print("Tube Hematocrit Sim Corr Shape:", parsed_data['tube_hematocrit_sim_corr'].shape)
'''

# Define a function to return the tube haematocrits in the right order
def get_simulation_haematocrit(vtk_path):

    # Extract the data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    default_simulation_haematocrit = cell_data_array['Vessel Haematocrit']
    relative_velocity = cell_data_array['Absolute Vessel Flow Rate m^3/s']/cell_data_array['Absolute Vessel Flow Rate m^3/s'][0]

    return default_simulation_haematocrit, relative_velocity

# Define a function to return the tube haematocrits in the right order (VTK 9+)
# def get_simulation_haematocrit(npz_path):
#     """
#     Reads the haematocrit values from a compressed .npz file that was previously extracted from a VTK file.

#     Parameters:
#         npz_path (str): Path to the .npz file containing the extracted data.

#     Returns:
#         numpy.ndarray: The array of vessel haematocrit values.
#     """
#     # Load the .npz file
#     data = np.load(npz_path)

#     # Access the vessel haematocrit array from the loaded data
#     default_simulation_haematocrit = data['cell_data_Vessel Haematocrit']

#     return default_simulation_haematocrit

# =============================================================================
# READ THE EXPERIMENTAL DATA
# =============================================================================

# Read the experimental data
# file_path = '/Users/vedang/Documents/GitHub/dphil-scripts/microfluidics_validation/exp_data/merlo/hematocrit_velocity_0035_10_mic_expe_honeycomb.am'
file_path = '/home/narain/Desktop/Scripts/microfluidics_validation/exp_data/merlo/hematocrit_velocity_0035_10_mic_expe_honeycomb.am'
summary_info, parsed_data = read_amira_ascii(file_path)

# Rotate the vertices by swapping x and y coordinates
rotated_vertices = np.array([[y, x, z] for x, y, z in parsed_data['vertices']])
# rotated_vertices = np.array([[x, y, z] for x, y, z in parsed_data['vertices']])

# Call the plotting function with rotated vertices
experimental_flow = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['velocity_rbc_exp'])
sim1_f = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['velocity_rbc_sim'])
experimental_hematocrit = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['tube_hematocrit_exp'])
sim1_h = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['tube_hematocrit_sim'])
# sim2_h = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['tube_hematocrit_sim_corr'])

# Convert to numpy arrays if they are not already
experimental_flow = np.array(experimental_flow)
sim1_f = np.array(sim1_f)
experimental_hematocrit = np.array(experimental_hematocrit)
sim1_h = np.array(sim1_h)
# sim2_h = np.array(sim2_h)

# =============================================================================
# READ MY SIMS
# =============================================================================

# Get the solvers (Mac)
# pries_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/PriesHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.npz')
# memory_4_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/MemoryHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.npz')
# memory_28_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_28/TestMerloMicrofluidicsNetwork/MemoryHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.npz')
# fung_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/FungHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.npz')

# Pick the simulation
# prefix= '/scratch/narain/testoutput/Merlo/OldParams'
# prefix= '/scratch/narain/testoutput/Merlo/NewViscParam'
# prefix= '/scratch/narain/testoutput/Merlo/NewABXVisc'
# prefix= '/scratch/narain/testoutput/Merlo/NewViscTermsAndSameParams'
# prefix= '/scratch/narain/testoutput/Merlo/NewHydDiam'
# prefix= '/scratch/narain/testoutput/Merlo/ViscUnchanged'
# prefix= '/scratch/narain/testoutput/Merlo/NewParamsHydDiam220OldViscTerm'
prefix= '/scratch/narain/testoutput'

# Not helpful
# prefix= '/scratch/narain/testoutput/Merlo/LowPress'
# prefix= '/scratch/narain/testoutput/Merlo/NoSquared'

# Read the haematocrit
constant_haematocrit, constant_flow = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/ConstantHaematocrit/FinalHaematocrit.vtp')
pries_haematocrit, pries_flow = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/PriesHaematocrit/FinalHaematocrit.vtp')
# memory_4_haematocrit, _ = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/MemoryHaematocrit/FinalHaematocrit.vtp')
# memory_28_haematocrit, _ = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/Memory28Haematocrit/FinalHaematocrit.vtp')
memory_24_haematocrit, _ = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/Memory24Haematocrit/FinalHaematocrit.vtp')
# memory_24_haematocrit, _ = get_simulation_haematocrit('/tmp/narain/testoutput/TestMerloMicrofluidicsNetwork/MemoryHaematocrit24/FinalHaematocrit.vtp')
fung_haematocrit, _ = get_simulation_haematocrit(prefix+'/TestMerloMicrofluidicsNetwork/FungHaematocrit/FinalHaematocrit.vtp')
# pries_haematocrit, _ = pries_haematocrit[non_zero_indices]

# =============================================================================
# CHECK FLOW
# =============================================================================

# Filter all arrays to only include the non-zero indices
# '''
# non_zero_indices = np.where(experimental_flow != 0)[0]
# experimental_flow = experimental_flow[non_zero_indices]
# sim1_f = sim1_f[non_zero_indices]
# pries_flow = pries_flow[non_zero_indices]
# constant_flow = constant_flow[non_zero_indices]
# non_zero_indices = np.where(experimental_hematocrit != 0)[0]
# experimental_hematocrit = experimental_hematocrit[non_zero_indices]
# sim1_h = sim1_h[non_zero_indices]
# sim2_h = sim2_h[non_zero_indices]
# '''

# Plot the data from the paper
'''
plt.figure(figsize=(10,10))
plt.scatter(experimental_flow, sim1_f, facecolors='none', edgecolors='red', s=75)
plt.scatter(experimental_flow, pries_flow, c='green', marker='x', s=75)
# plt.scatter(experimental_flow, constant_flow, c='green', marker='x', s=75)
min_val = 0
max_val = 1.05
plt.plot([min_val, max_val], [min_val, max_val], '--', color='black')
plt.xlim(min_val,max_val) 
plt.ylim(min_val,max_val)
plt.xlabel('experiment')
plt.ylabel('simulation')
# correlation_l, p_value = pearsonr(experimental_flow,sim1_f)
slope, intercept, r_value, p_value, std_err = linregress(experimental_flow, sim1_f)
plt.title(f'Normalised flow velocity with in vitro viscosity (r = {r_value:.2f}, m = {slope:.2f})')

# Show the plot
figure_title = 'relative_flow'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500)
'''

# =============================================================================
# COMPARE WITH CHASTE NETWORK
# =============================================================================

# Compute the r-values
# correlation_p, p_value = pearsonr(experimental_hematocrit, pries_haematocrit)
correlation_24, p_value = pearsonr(experimental_hematocrit, memory_24_haematocrit)
correlation_f, p_value = pearsonr(experimental_hematocrit, fung_haematocrit)
# correlation_4, p_value = pearsonr(experimental_hematocrit, memory_4_haematocrit)
# correlation_l_h, p_value = pearsonr(experimental_hematocrit,sim1_h)

# Plot the data
plt.figure(figsize=(10,10))
min_val = 0
max_val = 0.1
plt.plot([min_val, max_val], [min_val, max_val], '--' , color='black')
plt.xlim(min_val,max_val) 
plt.ylim(min_val,max_val)
plt.xlabel('experiment')
plt.ylabel('simulation')
slope, intercept, r_value, p_value, std_err = linregress(experimental_hematocrit,sim1_h)
plt.scatter(experimental_hematocrit,sim1_h, label=f'Merlo et al. (r = {r_value:.2f}, r_sq = {r_value**2:.2f}, m = {slope:.2f})', facecolors='none', edgecolors='red', alpha=0.7)
# plt.scatter(experimental_hematocrit,memory_4_haematocrit, label=rf'$\omega=4$ ($r$ = {correlation_4:.2f})', marker='D', c='orange')
plt.scatter(experimental_hematocrit,memory_24_haematocrit, label=rf'$\omega=24$ ($r$ = {correlation_24:.2f})', marker='D', c='b')
slope, intercept, r_value, p_value, std_err = linregress(experimental_hematocrit,pries_haematocrit)
plt.scatter(experimental_hematocrit,pries_haematocrit, label=f'Pries replication (r = {r_value:.2f}, r_sq = {r_value**2:.2f}, m = {slope:.2f})', marker='x', c='g', alpha=0.7)
# plt.scatter(experimental_hematocrit,fung_haematocrit, label=f'Fung ($r$ = {correlation_f:.2f})', marker='+')
# plt.scatter(experimental_hematocrit,memory_4_haematocrit, label=f'Memory ($r$ = {correlation_4:.2f})', marker='x', c='g')
# plt.title(f'Pries: {correlation_p:.2f}, M24: {correlation_24:.2f}, Fung: {correlation_f:.2f}, M4: {correlation_4:.2f})')
plt.title('Tube haematocrit with in vitro viscosity')
plt.legend()

# Show the plot
figure_title = 'haematocrit_comparison'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500)

# =============================================================================
# COMPUTE MSE INSTEAD
# =============================================================================

# Calculate the R-values
r_constant = np.sqrt(abs(r2_score(experimental_hematocrit, constant_haematocrit)))
r_pries = np.sqrt(abs(r2_score(experimental_hematocrit, pries_haematocrit)))
r_memory = np.sqrt(abs(r2_score(experimental_hematocrit, memory_24_haematocrit)))
r_fung = np.sqrt(abs(r2_score(experimental_hematocrit, fung_haematocrit)))
r_values = [r_constant, r_pries, r_memory, r_fung]

# Get the MSEs
constant_error = mean_squared_error(experimental_hematocrit,constant_haematocrit)
pries_error = mean_squared_error(experimental_hematocrit,pries_haematocrit)
memory_error = mean_squared_error(experimental_hematocrit,memory_24_haematocrit)
fung_error = mean_squared_error(experimental_hematocrit,fung_haematocrit)
values_mse = np.array([constant_error, pries_error, memory_error, fung_error])

# Create the plot
fig, ax1 = plt.subplots(figsize=(20, 10))
colors = ['#0072B2', '#009E73','#D55E00']
labels = ['Constant', 'Pries', 'Memory', 'Fung']

# Plot the MSE bars
bars_mse = ax1.bar(range(len(values_mse)), values_mse * 1e3, color='crimson', width=0.25, label='MSE', edgecolor='black', linewidth=2)
ax1.set_ylabel('MSE of network haematocrit (× $10^{-3}$)')
ax1.set_xticks(np.arange(len(labels)) + 0.25/2)
ax1.set_xticklabels(labels)
# ax1.grid(axis='y', linestyle='--', alpha=0.7)

# Create a secondary y-axis for the R-values
ax2 = ax1.twinx()
values_r = np.array(r_values)

# Plot the R-value bars with some offset to avoid overlap
bars_r = ax2.bar(np.arange(len(values_r)) + 0.25, values_r, color=colors[1], width=0.25, label='$r$', edgecolor='black', linewidth=2)
ax2.set_ylabel('correlation coefficient ($r$)', fontsize=30)

# Add a legend to differentiate between MSE and R-value bars
fig.legend(loc='upper right', bbox_to_anchor=(0.895, 0.87))

# Show and save the plot
figure_title = 'solver_bars'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500)
plt.show()








# import matplotlib.pyplot as plt
# import numpy as np

# # Data
# values = [constant_error, pries_error, memory_error, fung_error]
# labels = ['Constant', 'Pries', 'Memory', 'Fung']


# # Create the figure and bar chart
# plt.figure(figsize=(10, 6))  # Adjust the aspect ratio to be wider
# bars = plt.bar(range(len(values)), values, color=colors, edgecolor='black', linewidth=1.5)  # Add colours, edge colour, and line width for more contrast

# # Add value annotations above the bars
# for bar in bars:
#     yval = bar.get_height()
#     plt.text(bar.get_x() + bar.get_width() / 2, yval + 0.02, f'{yval:.2f}', ha='center', va='bottom', fontsize=12, fontweight='bold', color='black')

# # Beautify ticks and labels
# plt.xticks(range(len(values)), labels, rotation=45, ha='right', fontsize=14, fontweight='bold')
# plt.yticks(fontsize=12)
# plt.xlabel('Rule', fontsize=16, fontweight='bold', labelpad=15)  # Increase font size and add padding
# plt.ylabel('MSE (H)', fontsize=16, fontweight='bold', labelpad=15)

# # Add gridlines for better readability of values
# plt.grid(axis='y', linestyle='--', alpha=0.7)

# # Add a title
# plt.title('Comparison of Mean Squared Errors by Rule', fontsize=18, fontweight='bold', pad=20)

# # Add background style to make the plot stand out
# plt.gca().set_facecolor('#f0f0f0')  # Light grey background
# plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2)  # Adjust subplot margins

# # Show the plot
# plt.show()

