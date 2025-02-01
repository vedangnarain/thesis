#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 21:19:58 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Visualise normalisation.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
#import matplotlib.colors as colors
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import numpy as np
import os
import pandas as pd
# import scipy.io
# import seaborn as sns
import time

# Import tools for Paraview data
from calculate_surviving_fraction import *
from convert_oxygen_units import *
from get_paraview_data import *

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to return the metrics for a simulation
def get_simulation_stats(n_kills, folder_path_pre_kill_number, network_file_name, network_path, flat_field, field_data):    
    
    # Get reference coordinates from the unpruned network
    reference_node_coordinates, reference_network_vessel_nodes, reference_network_coordinates = get_reference_hexagonal_network(os.path.join(folder_path_pre_kill_number + '0', network_file_name))  # Just use the unpruned network as a reference network
    print(network_path)

    # Compute the metrics and get the graph
    metrics = get_bio_predictors(network_path, reference_node_coordinates, reference_network_coordinates, reference_network_vessel_nodes, pq_threshold, field_data, inlet_h)
    cell_NF, field_NF = compare_field_cells_NF(folder_path_pre_kill_number + str(n_kills) + '/results_from_time_0',hypoxic_threshold_list[0])
    # Compute the perfusion metrics
    # total_blood_flow = metrics['total_blood_volume'] * 60000000  # Convert to ml/s
    surface_area = max(field_data['x']) * max(field_data['y'])
    fcd = metrics['fcd_numerator'] / surface_area
    ppv = metrics['n_perfused_vessels'] / metrics['n_vessels'] if metrics['n_vessels'] != 0 else 0
    
    # Get the number of total points
    number_of_points = np.prod(np.shape(flat_field))  # More efficient way to get the product of dimensions

    # Calculate the number of points below the given hypoxic thresholds
    hypoxic_fraction_list = []
    for hypoxic_threshold in hypoxic_threshold_list:
        hypoxic_points = (flat_field < hypoxic_threshold).sum()
        hypoxic_fraction = hypoxic_points / number_of_points  # calculate the hypoxic fraction
        hypoxic_fraction_list.append(hypoxic_fraction)

    # Compute the expected SFs
    lewin_SF_cells = compute_Lewin_SF(cell_NF)
    # print('lewin',lewin_SF)     
    pre_RT = get_pvd_data(folder_path_pre_kill_number + str(n_kills) + '/results_from_time_0/results.pvd')
    wouters_computed_n_surviving = np.sum(calculate_array_Wouters_SF(convert_nM_to_mmHg(pre_RT['oxygen'])))
    wouters_SF_cells = wouters_computed_n_surviving/len(pre_RT)

    # Compute the expected SFs
    lewin_SF_field = compute_Lewin_SF(field_NF)
    # print('lewin',lewin_SF)     
    wouters_computed_n_surviving = np.sum(calculate_array_Wouters_SF(convert_nM_to_mmHg(flat_field)))
    wouters_SF_field = wouters_computed_n_surviving/len(flat_field)

    # Access subdatasets
    h_metrics = metrics['h_metrics']  # Nested dictionary with network metrics
    # network_metrics = metrics['network_metrics']  # Nested dictionary with network metrics
    tda_metrics = metrics['tda_metrics']  # Nested dictionary with network metrics
    
    # Compute average value for the node-based metrics over the entire network
    # pagerank_unweighted_avg = np.mean(list(network_metrics['pagerank_unweighted'].values()))  
    # betweenness_unweighted_avg = np.mean(list(network_metrics['betweenness_unweighted'].values()))    

    # Convert edge-based metrics into a list of tuples (u, v, value) for each metric
    # edge_original_ID = [data['original_ID'] for u, v, data in G_filtered.edges(data=True)]
    # edges = list(G_filtered.edges())  # Extract the edges from the graph

    # Create a DataFrame with descriptive column names
    simulation_df = pd.DataFrame({
        "selection": 0,
        "kills": [n_kills],
        "HF/RF": [hypoxic_fraction_list],
        "cell_NF": [cell_NF],
        "field_NF": [field_NF],
        "total_O2": [np.sum(flat_field)],
        "mean_O2": [np.mean(flat_field)],
        "min_O2": [np.amin(flat_field)],
        "50_O2": [np.percentile(flat_field, 50)],
        "max_O2": [np.amax(flat_field)],
        "SD_O2": [np.std(flat_field)],
        "n_vessels": [metrics['n_vessels']],
        "n_perfused_vessels": [metrics['n_perfused_vessels']],
        "n_unperfused_vessels": [metrics['n_unperfused_vessels']],
        "mean_diameter": [metrics['mean_diameter']],
        "mean_geometric_resistance": [metrics['mean_geometric_resistance']],
        "n_cycles": [metrics['n_cycles']],
        # "TBV": [total_blood_flow],
        "FCD": [fcd],
        "PPV": [ppv],
        "VDi": [metrics['vd_intersections']],
        # "VDu": [metrics['vd_uniques']],
        "n_connected_components": [metrics['n_connected_components']],
        "size_largest_connected_component": [tda_metrics['size_largest_connected_component']],

        # Store haematocrit metrics
        # "h_std_dev": [h_metrics['h_std_dev']],
        # "h_cv": [h_metrics['h_cv']],
        # # "h_gini": [h_metrics['h_gini']],
        # "h_shannon": [h_metrics['h_shannon']],
        "h_inverse_simpson": [h_metrics['h_inverse_simpson']],
        # "h_hoover": [h_metrics['h_hoover']],
        # "h_haemo_potential": [h_metrics['h_haemo_potential']],
        # "h_theil": [h_metrics['h_theil']],
        # "h_atkinson": [h_metrics['h_atkinson']],
        # "h_pietra": [h_metrics['h_pietra']],
        # "h_mld": [h_metrics['h_mld']],

        # Store average of node-based network metrics
        # "pagerank_unweighted_avg": [pagerank_unweighted_avg],
        # "betweenness_unweighted_avg": [betweenness_unweighted_avg], # 

        # Store the radiotherapy outcomes
        "wouters_SF_cells": [wouters_SF_cells],
        "lewin_SF_cells": [lewin_SF_cells],
        "wouters_SF_field": [wouters_SF_field],
        "lewin_SF_field": [lewin_SF_field]
    })

    # Store the metrics for each edge
    # edge_metrics_df = pd.DataFrame({
    #     "edge": [(u, v) for u, v in edges],  # Store each edge as a tuple
    #     "edge_original_ID": edge_original_ID,
    # })

    return simulation_df, reference_node_coordinates, reference_network_vessel_nodes

# Define a function to return statistics for a single network layout with different kills
def get_single_network_stats(main_folder_path):

    layout_folder_path = main_folder_path + solver_name + 'Haematocrit'

    # Initialise empty dfs to store simulations for the whole layout
    layout_df = pd.DataFrame()
    killed_vessels_df = pd.DataFrame()

    # List all directories in the selection folder
    kills_dirs = []
    for entry in os.listdir(layout_folder_path):
        if entry.startswith('Kills') and os.path.isdir(os.path.join(layout_folder_path, entry)):  # Check if the entry is a directory and matches the "Kills" pattern
            kills_dirs.append(entry)
    
    # Sort the Kills directories in ascending order based on their numeric value
    kills_dirs.sort(key=lambda x: int(x.replace('Kills', '')))  # Extract the numeric part and convert to int
    
    # Iterate over all the kills     
    # previous_edge_metrics_df = None  # Keep track of edge metrics of the previous network
    for kills in kills_dirs:

        # Set the file path
        n_kills = int(kills[5:])
        kills_folder_path = os.path.join(layout_folder_path, kills)
        folder_path_pre_kill_number = os.path.join(layout_folder_path, kills[:5])
        network_file_name = 'FinalHaematocrit.vtp'    
        # network_file_name = 'results_from_time_0/VesselNetwork_inc_1.vtp'    
        # field_path =  folder_path + 'results_from_time_0/oxygen_solution_0.vti'
        network_path = os.path.join(kills_folder_path, network_file_name)
        field_path =  os.path.join(kills_folder_path, 'results_from_time_0/oxygen_solution_0.vti')

        # Compute the O2 distribution
        field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
        flat_field = field_data['oxygen'].to_numpy().flatten()
        # print(np.mean(flat_field), np.max(flat_field))
        if np.sum(flat_field)<0.000001:
            print('No oxygen detected!')
        else:
            simulation_df, reference_node_coordinates, reference_network_vessel_nodes = get_simulation_stats(n_kills, folder_path_pre_kill_number, network_file_name, network_path, flat_field, field_data)   
            layout_df = pd.concat([layout_df, simulation_df], ignore_index=True)    

        # Update previous_edge_metrics_df to the current network's metrics for the next iteration
        # previous_edge_metrics_df = edge_metrics_df
    
    return layout_df, killed_vessels_df

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Toggles
seed_selection = 60
# max_layout = 358

# Enter details to allow looping over folders
# sd_list = ['8.68', '13.23', '17.49']
# mean_list = ['22.76', '28.5', '33.64']
solver_list = ['/Constant', '/Pries', '/Memory', '/Fung']
# max_kills = 2
# kills_list = [str(x) for x in range(0, max_kills + 1)]
# hypoxic_threshold_list = [4116, 20578]  # 3 mmHg and 15 mmHg (15-8 is physiological hypoxia, 8-3 is pathological, <3 is radiobiological)
hypoxic_threshold_list = [4500, 10975]  # 3.28 mmHg and 8 mmHg
# diameter_threshold_list = [22, 35, 50]  # in um

solver_name = solver_list[2]
# inlet_radius_m = 7.500000e-05
inlet_h = 0.3
pq_threshold = 1.e-13
# main_folder_path = '/scratch/narain/Voronoi/' + str(seed_selection) + ' Seeds/H ' + str(inlet_h) + '/no_rt_size_pruning/TestVoronoiNetwork'
# save_path = '/home/narain/Desktop/Scripts/demo_datasets/vor_size/raw_voronoi_size.h5'  # Specify the file path
# main_folder_path = '/scratch/narain/Voronoi/' + str(seed_selection) + ' Seeds/H ' + str(inlet_h) + '/no_rt_flow_pruning/TestVoronoiNetwork'
main_folder_path = '/home/narain/Desktop/scratch/narain/Biological/memory_size_pruning/TestMetastaticNetworkTissue'
# main_folder_path = '/scratch/narain/testoutput/memory_flow_pruning/TestMetastaticNetworkTissue'
# save_path = '/home/narain/Desktop/Scripts/demo_datasets/bio_network.h5'  # Specify the file path
hdf_file_name = main_folder_path + solver_name + 'Haematocrit/metrics.h5'

# =============================================================================
# GENERATE AND SAVE DATA FOR PREDICTIONS
# =============================================================================
#'''
# Get the stats for all solvers, varying the value of alpha (change to read=1 to extract from .vti files directly)
layout_df, killed_vessels_df = get_single_network_stats(main_folder_path)
combined_df = layout_df

# Split up HFs
combined_df[['RF', 'HF']] = combined_df['HF/RF'].apply(pd.Series)
combined_df = combined_df.drop('HF/RF', axis=1)

# Create NF column
combined_df['NF'] = 1 - combined_df['RF']

# Create loops/vessel column
# combined_df['loops/vessel'] = combined_df['n_cycles'] / combined_df['n_vessels']

# Create resistance/loop column
# combined_df['resistance/loop'] = np.where(combined_df['loops/vessel'] != 0, combined_df['mean_geometric_resistance'] / combined_df['loops/vessel'], 0)

# Save DataFrames to HDF5

#'''
# =============================================================================
# COMBINE THE DATASETS
# =============================================================================

with pd.HDFStore(hdf_file_name) as store:
    subset = store['combined_df']  # Retrieve the DataFrame using the given key

# # Specify the HDF5 file name and key for the DataFrame
# hdf_file_name = 'metrics.h5'

# # Combine the DataFrames from each selection for 'combined_df'
# imported_raw_df = combine_hdf5_files(main_folder_path + solver_name + 'Haematocrit', max_layout, hdf_file_name, key='combined_df')
# # main_folder_path + solver_name + 'Haematocrit'
# imported_raw_df.to_hdf(save_path, key='df', mode='w')  # Save the DataFrame
# imported_raw_df = pd.read_hdf(save_path, key='df')


# Plotting the best only
fig, ax1 = plt.subplots(figsize=(7, 6))
width = 3

# Ensure LaTeX-style font globally
# sns.reset_orig()
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})

# Compute the percentage of the network pruned
kills_perc = ((subset['n_vessels'].max() - subset['n_vessels']) / subset['n_vessels'].max()) * 100

# Plot NF
NF_line, = ax1.plot(kills_perc, subset['cell_NF'], linestyle='-', lw=width, label='NF')

# Plot SF
woutersSF_line, = ax1.plot(kills_perc, subset['wouters_SF_cells'], label='SF', linestyle='--', lw=width)

# Add vertical lines at x=5 and x=17
# ax1.axvline(5, color='blue', linestyle='--', linewidth=1.5)
# ax1.axvline(17, color='purple', linestyle='--', linewidth=1.5)

# Configure primary y-axis
ax1.set_xlabel('% of network pruned', fontdict={'family': 'STIXGeneral', 'size': 22})
ax1.set_ylabel('NF/SF', fontdict={'family': 'STIXGeneral', 'size': 22})
ax1.set_xlim(0)
ax1.set_ylim(0,0.8)
ax1.set_facecolor('white')
ax1.grid(False)

# Ensure all spines are visible except the top
ax1.spines['left'].set_visible(True)
ax1.spines['bottom'].set_visible(True)
ax1.spines['right'].set_visible(True)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_linewidth(1.5)
ax1.spines['bottom'].set_linewidth(1.5)
ax1.spines['right'].set_linewidth(1.5)

# Secondary y-axis for relative ratios
ax2 = ax1.twinx()

# Compute and plot RR
wouters_rr = compute_RR(subset['wouters_SF_cells'], anoxic_SF=calculate_array_Wouters_SF(0))
ax2.plot(kills_perc, wouters_rr, label='RR', linestyle=':', lw=width, color='red')
ax2.set_ylim(1,2.1)
ax2.set_ylabel('RR', fontdict={'family': 'STIXGeneral', 'size': 22})
ax2.set_facecolor('white')
ax2.grid(False)

# Ensure all spines are visible except the top for secondary y-axis
ax2.spines['left'].set_visible(True)
ax2.spines['bottom'].set_visible(True)
ax2.spines['right'].set_visible(True)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_linewidth(1.5)
ax2.spines['bottom'].set_linewidth(1.5)
ax2.spines['right'].set_linewidth(1.5)

# Highlight intervals where wouters_rr > its initial value at x=0
tolerance = 0.004  # the NF must exceed the original value by this much
initial_wouters_rr = wouters_rr.iloc[0]
above_initial = wouters_rr > initial_wouters_rr+tolerance
in_interval = False
start = None

normalisation_label_added = False  # Ensure the label appears only once

for i in range(len(kills_perc)):
    if above_initial.iloc[i] and not in_interval:
        # Start a new interval
        start = kills_perc.iloc[i]
        in_interval = True
    elif not above_initial.iloc[i] and in_interval:
        # End the interval and shade it
        if not normalisation_label_added:
            ax1.axvspan(start, kills_perc.iloc[i - 1], color='g', alpha=0.3, label='window')
            normalisation_label_added = True
        else:
            ax1.axvspan(start, kills_perc.iloc[i - 1], color='g', alpha=0.3)
        in_interval = False
# Handle the case where the interval goes to the end of the data
if in_interval:
    if not normalisation_label_added:
        ax1.axvspan(start, kills_perc.iloc[-1], color='g', alpha=0.3, label='window')
    else:
        ax1.axvspan(start, kills_perc.iloc[-1], color='g', alpha=0.3)

# Combine legends from both axes
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='lower right')
# ax1.legend(lines + lines2, labels + labels2)