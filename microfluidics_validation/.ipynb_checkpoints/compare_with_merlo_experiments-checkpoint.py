#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:54:28 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

<description>

Tested in Python 3.11
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from pathlib import Path
from matplotlib.pyplot import figure
import matplotlib.colors as mcolors
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
import scipy.stats

# Import tools for Paraview data
import sys
import os
parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_directory)
from get_paraview_data import *

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# READ MERLO DATA
# =============================================================================

# import re
# import numpy as np

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

# Example usage
file_path = '/Users/vedang/Documents/GitHub/dphil-scripts/microfluidics_validation/exp_data/merlo/hematocrit_velocity_0035_10_mic_expe_honeycomb.am'
summary_info, parsed_data = read_amira_ascii(file_path)

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

# =============================================================================
# PLOT THE EXPERIMENTAL NETWORK
# =============================================================================

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
        annotation_text = f"{new_label}: {tube_hematocrit_exp[original_index]:.3f}"
        plt.text(midpoint_x, midpoint_y, annotation_text, color="black", fontsize=8, ha="center", va="center")
        
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
    plt.xlim(x_center - max_range, x_center + max_range)
    plt.ylim(y_center - max_range, y_center + max_range)
    
    # Set aspect to equal for a square plot
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Add a color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(tube_hematocrit_exp)  # Setting the array for color mapping
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Tube Hematocrit (Experimental)')
    
    plt.xlabel('X (μm)')
    plt.ylabel('Y (μm)')
    plt.title('Square Network Plot with Tube Hematocrit Experimental Values and Ordered Edge Labels')
    plt.show()
    
    return sorted_hematocrit_values

# Rotate the vertices by swapping x and y coordinates
rotated_vertices = np.array([[y, x, z] for x, y, z in parsed_data['vertices']])

# Call the plotting function with rotated vertices
experimental_hematocrit = plot_network(rotated_vertices, parsed_data['edges'], parsed_data['tube_hematocrit_exp'])

# =============================================================================
# READ CHASTE NETWORK
# =============================================================================

# Specify path
# base_path = 
# sim_name= '/TestMerloMicrofluidicsNetwork/'
# vtk_path = 'Haematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.vtp'
# merlo_path = base_path + sim_name + splitting_rules[2] + vtk_path

# Define a function to return the tube haematocrits in the right order
def get_simulation_haematocrit(vtk_path):

    # Extract the data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    default_simulation_haematocrit = cell_data_array['Vessel Haematocrit']
    
    return default_simulation_haematocrit

# Get the solvers
pries_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/PriesHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.vtp')
memory_4_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/MemoryHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.vtp')
memory_28_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_28/TestMerloMicrofluidicsNetwork/MemoryHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.vtp')
fung_haematocrit = get_simulation_haematocrit('/Users/vedang/Downloads/Microfluidics Network/Merlo Network/old_B/omega_4/TestMerloMicrofluidicsNetwork/FungHaematocrit/Selection1/Sigma8.68/Mu22.76/Kills0/FinalHaematocrit.vtp')

# Plot the data
plt.figure(figsize=(8, 8))
min_val = 0
max_val = 0.1
plt.plot([min_val, max_val], [min_val, max_val], '--' , color='black')
plt.xlim(min_val,max_val) 
plt.ylim(min_val,max_val)
plt.xlabel('Experimental Haematocrit')
plt.ylabel('Model Haematocrit')
correlation_p, p_value = pearsonr(experimental_hematocrit, pries_haematocrit)
correlation_28, p_value = pearsonr(experimental_hematocrit, memory_28_haematocrit)
correlation_f, p_value = pearsonr(experimental_hematocrit, fung_haematocrit)
correlation_4, p_value = pearsonr(experimental_hematocrit, memory_4_haematocrit)
# plt.scatter(merlo_exp_pads,memory_4_haematocrit, label=rf'$\omega=4$ ($r$ = {correlation_4:.2f})', marker='x', c='g')
# plt.scatter(merlo_exp_pads,memory_28_haematocrit, label=rf'$\omega=28$ ($r$ = {correlation_28:.2f})', marker='D', c='r')
# plt.title(f'Pries: {correlation_p:.2f}, M28: {correlation_28:.2f}, Fung: {correlation_f:.2f}, M4: {correlation_4:.2f})')
# plt.title('Average Haematocrit Values for Custom Groups')
plt.legend()

# Show the plot
# figure_title = 'all_solvers'
figure_title = 'memory_solvers'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
plt.savefig(file_path, dpi=500)
plt.show()





