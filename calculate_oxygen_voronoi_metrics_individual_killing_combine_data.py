#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:56:51 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Create a dataframe of all the Voronoi simulations in a folder.

Colour palette from https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7.
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
def get_simulation_stats(which_layout, n_kills, folder_path_pre_kill_number, network_file_name, network_path, flat_field, field_data):    
    
    # Get reference coordinates from the unpruned network
    reference_node_coordinates, reference_network_vessel_nodes, reference_network_coordinates = get_reference_hexagonal_network(os.path.join(folder_path_pre_kill_number + '0', network_file_name))  # Just use the unpruned network as a reference network
    print(network_path)
    
    # Compute the metrics and get the graph
    # n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, n_cycles, total_blood_flow_m3_per_s, fcd_numerator, vd_intersections, vd_uniques, n_connected_components = get_voronoi_predictors(network_path, reference_node_coordinates, reference_network_coordinates, pq_threshold, field_data)
    metrics, G_filtered = get_voronoi_predictors(network_path, reference_node_coordinates, reference_network_coordinates, reference_network_vessel_nodes, pq_threshold, field_data, inlet_h)

    # Compute the perfusion metrics
    total_blood_flow = metrics['total_blood_volume'] * 60000000  # Convert to ml/s
    surface_area = max(field_data['x']) * max(field_data['y'])
    # print(max(field_data['x']),max(field_data['y']))
    # print(fcd_numerator, surface_area)
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
    lewin_SF = compute_Lewin_SF(1 - hypoxic_fraction_list[0])
    # print('lewin',lewin_SF)     
    wouters_computed_n_surviving = np.sum(calculate_array_Wouters_SF(convert_nM_to_mmHg(flat_field)))
    wouters_SF = wouters_computed_n_surviving/len(flat_field)

    # Access subdatasets
    h_metrics = metrics['h_metrics']  # Nested dictionary with network metrics
    network_metrics = metrics['network_metrics']  # Nested dictionary with network metrics
    tda_metrics = metrics['tda_metrics']  # Nested dictionary with network metrics
    
    # Compute average value for the node-based metrics over the entire network
    # closeness_unweighted_avg = np.mean(list(network_metrics['closeness_unweighted'].values()))
    # closeness_weighted_avg = np.mean(list(network_metrics['closeness_weighted'].values()))
    # eigenvector_unweighted_avg = np.mean(list(network_metrics['eigenvector_unweighted'].values()))
    # eigenvector_weighted_avg = np.mean(list(network_metrics['eigenvector_weighted'].values()))
    pagerank_unweighted_avg = np.mean(list(network_metrics['pagerank_unweighted'].values()))
    # pagerank_weighted_avg = np.mean(list(network_metrics['pagerank_weighted'].values()))
    # clustering_unweighted_avg = np.mean(list(network_metrics['clustering_unweighted'].values()))
    # clustering_weighted_avg = np.mean(list(network_metrics['clustering_weighted'].values()))    
    # edge_betweenness_unweighted_avg = np.mean(list(network_metrics['edge_betweenness_unweighted'].values()))    
    # edge_betweenness_weighted_avg = np.mean(list(network_metrics['edge_betweenness_weighted'].values()))    
    betweenness_unweighted_avg = np.mean(list(network_metrics['betweenness_unweighted'].values()))    
    # betweenness_weighted_avg = np.mean(list(network_metrics['betweenness_weighted'].values()))    
    # average_path_length_unweighted = network_metrics['average_path_length_unweighted']    
    # average_path_length_weighted = network_metrics['average_path_length_weighted']    

    # Convert edge-based metrics into a list of tuples (u, v, value) for each metric
    edge_original_ID = [data['original_ID'] for u, v, data in G_filtered.edges(data=True)]
    edges = list(G_filtered.edges())  # Extract the edges from the graph
    # edge_closeness_unweighted = [(u, v, network_metrics['edge_closeness_unweighted'][(u, v)]) for u, v in edges]
    # edge_closeness_weighted = [(u, v, network_metrics['edge_closeness_weighted'][(u, v)]) for u, v in edges]
    # edge_betweenness_unweighted = [(u, v, network_metrics['edge_betweenness_unweighted'][(u, v)]) for u, v in edges]
    # edge_betweenness_weighted = [(u, v, network_metrics['edge_betweenness_weighted'][(u, v)]) for u, v in edges]
    # edge_eigenvector_unweighted = [(u, v, network_metrics['edge_eigenvector_unweighted'][(u, v)]) for u, v in edges]
    # edge_eigenvector_weighted = [(u, v, network_metrics['edge_eigenvector_weighted'][(u, v)]) for u, v in edges]
    # edge_pagerank_unweighted = [(u, v, network_metrics['edge_pagerank_unweighted'][(u, v)]) for u, v in edges]
    # edge_pagerank_weighted = [(u, v, network_metrics['edge_pagerank_weighted'][(u, v)]) for u, v in edges]
    # edge_clustering_unweighted = [(u, v, network_metrics['edge_clustering_unweighted'][(u, v)]) for u, v in edges]
    # edge_clustering_weighted = [(u, v, network_metrics['edge_clustering_weighted'][(u, v)]) for u, v in edges]
    # edge_distance_to_inlet = [(u, v, network_metrics['edge_distance_to_inlet'][(u, v)]) for u, v in edges]
    # edge_distance_to_outlet = [(u, v, network_metrics['edge_distance_to_outlet'][(u, v)]) for u, v in edges]
    
    # Repeat for the max versions
    # edge_closeness_max_unweighted = [(u, v, network_metrics['edge_closeness_max_unweighted'][(u, v)]) for u, v in edges]
    # edge_closeness_max_weighted = [(u, v, network_metrics['edge_closeness_max_weighted'][(u, v)]) for u, v in edges]
    # edge_betweenness_max_unweighted = [(u, v, network_metrics['edge_betweenness_max_unweighted'][(u, v)]) for u, v in edges]
    # edge_betweenness_max_weighted = [(u, v, network_metrics['edge_betweenness_max_weighted'][(u, v)]) for u, v in edges]
    # edge_eigenvector_max_unweighted = [(u, v, network_metrics['edge_eigenvector_max_unweighted'][(u, v)]) for u, v in edges]
    # edge_eigenvector_max_weighted = [(u, v, network_metrics['edge_eigenvector_max_weighted'][(u, v)]) for u, v in edges]
    # edge_pagerank_max_unweighted = [(u, v, network_metrics['edge_pagerank_max_unweighted'][(u, v)]) for u, v in edges]
    # edge_pagerank_max_weighted = [(u, v, network_metrics['edge_pagerank_max_weighted'][(u, v)]) for u, v in edges]
    # edge_clustering_max_unweighted = [(u, v, network_metrics['edge_clustering_max_unweighted'][(u, v)]) for u, v in edges]
    # edge_clustering_max_weighted = [(u, v, network_metrics['edge_clustering_max_weighted'][(u, v)]) for u, v in edges]

    # Create a DataFrame with descriptive column names
    simulation_df = pd.DataFrame({
        "selection": [which_layout],
        "kills": [n_kills],
        "HF/RF": [hypoxic_fraction_list],
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
        "h_std_dev": [h_metrics['h_std_dev']],
        "h_cv": [h_metrics['h_cv']],
        # "h_gini": [h_metrics['h_gini']],
        "h_shannon": [h_metrics['h_shannon']],
        "h_inverse_simpson": [h_metrics['h_inverse_simpson']],
        "h_hoover": [h_metrics['h_hoover']],
        "h_haemo_potential": [h_metrics['h_haemo_potential']],
        "h_theil": [h_metrics['h_theil']],
        "h_atkinson": [h_metrics['h_atkinson']],
        "h_pietra": [h_metrics['h_pietra']],
        # "h_mld": [h_metrics['h_mld']],

        # Store average of node-based network metrics
        # "closeness_unweighted_avg": [closeness_unweighted_avg], #worst 
        # "closeness_weighted_avg": [closeness_weighted_avg],   #worst
        # "eigenvector_unweighted_avg": [eigenvector_unweighted_avg], # 
        # "eigenvector_weighted_avg": [eigenvector_weighted_avg], # 
        "pagerank_unweighted_avg": [pagerank_unweighted_avg],
        # "pagerank_weighted_avg": [pagerank_weighted_avg],#
        # "clustering_unweighted_avg": [clustering_unweighted_avg],# 
        # "clustering_weighted_avg": [clustering_weighted_avg], # 
        "betweenness_unweighted_avg": [betweenness_unweighted_avg], # 
        # "betweenness_weighted_avg": [betweenness_weighted_avg], 
        # "edge_betweenness_weighted_avg": [edge_betweenness_weighted_avg], 
        # "edge_betweenness_unweighted_avg": [edge_betweenness_unweighted_avg],#
        # "average_path_length_unweighted": [average_path_length_unweighted], #   
        # "average_path_length_weighted": [average_path_length_weighted] # 
        
        # Store the radiotherapy outcomes
        "wouters_SF": [wouters_SF],
        "lewin_SF": [lewin_SF]
    })

    # print('n_edges ', len(edges))
    # print('n_edges_act ', len(metrics['segment_haematocrit']))   

    # Store the metrics for each edge
    edge_metrics_df = pd.DataFrame({
        "edge": [(u, v) for u, v in edges],  # Store each edge as a tuple
        "edge_original_ID": edge_original_ID,
        # "edge_closeness_unweighted": [x[2] for x in edge_closeness_unweighted], #worst
        # "edge_closeness_weighted": [x[2] for x in edge_closeness_weighted], #worst 
        # "edge_betweenness_unweighted": [x[2] for x in edge_betweenness_unweighted], # 
        # "edge_betweenness_weighted": [x[2] for x in edge_betweenness_weighted],# 
        # "edge_eigenvector_unweighted": [x[2] for x in edge_eigenvector_unweighted],# 
        # "edge_eigenvector_weighted": [x[2] for x in edge_eigenvector_weighted],# 
        # "edge_pagerank_unweighted": [x[2] for x in edge_pagerank_unweighted],
        # "edge_pagerank_weighted": [x[2] for x in edge_pagerank_weighted], #
        # "edge_clustering_unweighted": [x[2] for x in edge_clustering_unweighted],# 
        # "edge_clustering_weighted": [x[2] for x in edge_clustering_weighted],# 
        # "edge_distance_to_inlet": [x[2] for x in edge_distance_to_inlet], 
        # "edge_distance_to_outlet": [x[2] for x in edge_distance_to_outlet],
        
        # Store the max versions as well
        # "edge_closeness_max_unweighted": [x[2] for x in edge_closeness_max_unweighted], #worst
        # "edge_closeness_max_weighted": [x[2] for x in edge_closeness_max_weighted], #worst 
        # "edge_betweenness_max_unweighted": [x[2] for x in edge_betweenness_max_unweighted], # 
        # "edge_betweenness_max_weighted": [x[2] for x in edge_betweenness_max_weighted],# 
        # "edge_eigenvector_max_unweighted": [x[2] for x in edge_eigenvector_max_unweighted],
        # "edge_eigenvector_max_weighted": [x[2] for x in edge_eigenvector_max_weighted],# 
        # "edge_pagerank_max_unweighted": [x[2] for x in edge_pagerank_max_unweighted] 
        # "edge_pagerank_max_weighted": [x[2] for x in edge_pagerank_max_weighted], #
        # "edge_clustering_max_unweighted": [x[2] for x in edge_clustering_max_unweighted],# 
        # "edge_clustering_max_weighted": [x[2] for x in edge_clustering_max_weighted],# 
    })

    return simulation_df, edge_metrics_df, reference_node_coordinates, reference_network_vessel_nodes

# Define a function to return statistics for a single network layout with different kills
def get_layout_stats(which_layout):

    layout_folder_path = main_folder_path + solver_name + 'Haematocrit/Selection' + str(which_layout)

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
    previous_edge_metrics_df = None  # Keep track of edge metrics of the previous network
    for kills in kills_dirs:

        # Set the file path
        # folder_path_pre_kill = main_folder_path + solver_name + 'Haematocrit/Selection' + str(which_layout) + '/Sigma8.68/Mu22.76/Kills'
        # folder_path_pre_kill = main_folder_path + solver_name + 'Haematocrit/Selection' + str(which_layout) + '/Kills'
        n_kills = int(kills[5:])
        kills_folder_path = os.path.join(layout_folder_path, kills)
        folder_path_pre_kill_number = os.path.join(layout_folder_path, kills[:5])
        network_file_name = 'FinalHaematocrit.vtp'    
        # network_file_name = 'results_from_time_0/VesselNetwork_inc_1.vtp'    
        network_path = os.path.join(kills_folder_path, network_file_name)
        # field_path =  folder_path + 'results_from_time_0/oxygen_solution_0.vti'
        field_path =  os.path.join(kills_folder_path, 'oxygen_solution_0.vti')

        # Compute the O2 distribution
        field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
        # middle_field = get_hex_domain(field_data, field_spacing)  # get the distribution in the middle of the field (replace with designated file)
        # field_data = get_custom_domain(field_data, field_spacing, 0, 100, 0, 100)
        flat_field = field_data['oxygen'].to_numpy().flatten()
        # print(np.mean(flat_field), np.max(flat_field))
        if np.sum(flat_field)<0.000001:
            print('No oxygen detected!')
        else:
            simulation_df, edge_metrics_df, reference_node_coordinates, reference_network_vessel_nodes = get_simulation_stats(which_layout, n_kills, folder_path_pre_kill_number, network_file_name, network_path, flat_field, field_data)   
            layout_df = pd.concat([layout_df, simulation_df], ignore_index=True)    

        # Log the details of the vessel killed
        if int(n_kills)>0:
            previous_network_path = os.path.join(folder_path_pre_kill_number + str(int(n_kills)-1), network_file_name)
            killed_vessel_df = log_killed_vessel(network_path, previous_network_path, previous_edge_metrics_df, reference_node_coordinates, reference_network_vessel_nodes)
            killed_vessel_df.insert(0, 'kills', n_kills)
            killed_vessel_df.insert(0, 'selection', which_layout)
            killed_vessels_df = pd.concat([killed_vessels_df, killed_vessel_df], ignore_index=True)
    
        # Update previous_edge_metrics_df to the current network's metrics for the next iteration
        previous_edge_metrics_df = edge_metrics_df
    
    return layout_df, killed_vessels_df

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Toggles
seed_selection = 60
max_layout = 1001

# Enter details to allow looping over folders
# sd_list = ['8.68', '13.23', '17.49']
# mean_list = ['22.76', '28.5', '33.64']
solver_list = ['/Constant', '/Pries', '/Memory', '/Fung']
# max_kills = 2
# kills_list = [str(x) for x in range(0, max_kills + 1)]
# hypoxic_threshold_list = [4116, 20578]  # 3 mmHg and 15 mmHg (15-8 is physiological hypoxia, 8-3 is pathological, <3 is radiobiological)
hypoxic_threshold_list = [4500, 10975]  # 3.28 mmHg and 8 mmHg
# diameter_threshold_list = [22, 35, 50]  # in um

solver_name = solver_list[1]
# inlet_radius_m = 7.500000e-05
inlet_h = 0.3
pq_threshold = 1.e-13
# main_folder_path = '/scratch/narain/Voronoi/' + str(seed_selection) + ' Seeds/H ' + str(inlet_h) + '/no_rt_size_pruning/TestVoronoiNetwork'
# save_path = '/home/narain/Desktop/Scripts/demo_datasets/vor_size/raw_voronoi_size.h5'  # Specify the file path
main_folder_path = '/scratch/narain/Voronoi/' + str(seed_selection) + ' Seeds/H ' + str(inlet_h) + '/no_rt_flow_pruning/TestVoronoiNetwork'
save_path = '/home/narain/Desktop/Scripts/demo_datasets/vor_flow/raw_voronoi_flow.h5'  # Specify the file path

# =============================================================================
# GENERATE AND SAVE DATA FOR PREDICTIONS
# =============================================================================
# '''
# Repeat for all selections
for which_layout in range(1, max_layout+1):

    # Get the stats for all solvers, varying the value of alpha (change to read=1 to extract from .vti files directly)
    layout_df, killed_vessels_df = get_layout_stats(which_layout)
    combined_df = layout_df
    
    # Split up HFs
    combined_df[['RF', 'HF']] = combined_df['HF/RF'].apply(pd.Series)
    combined_df = combined_df.drop('HF/RF', axis=1)
    
    # Create NF column
    combined_df['NF'] = 1 - combined_df['RF']

    # Create loops/vessel column
    combined_df['loops/vessel'] = combined_df['n_cycles'] / combined_df['n_vessels']
    
    # Create resistance/loop column
    combined_df['resistance/loop'] = np.where(combined_df['loops/vessel'] != 0, combined_df['mean_geometric_resistance'] / combined_df['loops/vessel'], 0)
    
    # # Rectify flow rate error by using values at kills=1 for kills=0
    # for i in range(len(combined_df)-1):
    #     if (combined_df.at[i, 'kills'] == 0):
    #         combined_df.at[i, 'PF'] = combined_df.at[i+1, 'PF']
    #         combined_df.at[i, 'n_perfused_vessels_composite'] = combined_df.at[i+1, 'n_perfused_vessels_composite']
    #         combined_df.at[i, 'n_unperfused_vessels_composite'] = combined_df.at[i+1, 'n_unperfused_vessels_composite']

    # Save DataFrames to HDF5
    hdf_file_name = main_folder_path + solver_name + 'Haematocrit/Selection' + str(which_layout) + '/metrics.h5'
    with pd.HDFStore(hdf_file_name) as store:
        store['combined_df'] = combined_df
        store['killed_vessels_df'] = killed_vessels_df
# '''
# =============================================================================
# COMBINE THE DATASETS
# =============================================================================

# Specify the HDF5 file name and key for the DataFrame
hdf_file_name = 'metrics.h5'

# Combine the DataFrames from each selection for 'combined_df'
imported_raw_df = combine_hdf5_files(main_folder_path + solver_name + 'Haematocrit', max_layout, hdf_file_name, key='combined_df')
# main_folder_path + solver_name + 'Haematocrit'
imported_raw_df.to_hdf(save_path, key='df', mode='w')  # Save the DataFrame
imported_raw_df = pd.read_hdf(save_path, key='df')

# Describe the dataset
print("Total networks: ",imported_raw_df.shape[0])
print("Total number of features (as number of columns): ", imported_raw_df.shape[1])

#Check for null values
null_values = imported_raw_df.isnull().values.any()
if null_values == True:
    print("There are some missing values in imported_raw_df")
else:
    print("There are no missing values in the imported_raw_df dataset")

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))
