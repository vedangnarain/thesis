#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Feb  2 13:48:25 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This script computes different vascular metrics.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
# from collections import deque
#from vtk import *
# from vtk.util.numpy_support import vtk_to_numpy
# from get_perfusion_metrics import *  # import tools for perfusion metrics
# import math
import networkx as nx
import numpy as np
# import os
# import pandas as pd
# import pyvista as pv
# import scipy.interpolate
# import vtk
# import xml.etree.ElementTree as ET
# from scipy.stats import skew, kurtosis

# Set up plotting
import matplotlib
# import matplotlib.pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# PERFUSION METRICS
# =============================================================================

# Define a function to return the following perfusion metrics: Total Blood Vlow (in m3/s); FCD (functional capillary density); Vessel Density; Heterogeneity Index (for the flow velocity)
def get_perfusion_metrics_forking(vessel_data_array, pq_threshold, reference_rank_lengths_um, reference_network_coordinates, segment_xyxy, field_data_um):
    
    # Calculate the TBV (total blood volume) (in m3/s) (NOTE: This is broken. I need to calculate the volume of the vessel but I don't know the definition of TBV really.)
    '''
    first_rank_zero_index = np.where(vessel_data_array['Vessel Owner Rank'] == 0)[0][0]  # Gets the index of the first inlet/outlet vessel in the array (assumes all inlets/outlets are equal)
    total_blood_volume = vessel_data_array['Absolute Vessel Flow Rate m^3/s'][first_rank_zero_index]  # Get the flow rate for the rank 0 vessel
    '''
    total_blood_volume = np.nan
    
    # Calculate the FCD (functional capillary density)
    perfused_vessel_indices = np.where(vessel_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold)[0] 
    perfused_vessel_ranks = vessel_data_array['Vessel Owner Rank'][perfused_vessel_indices]
    perfused_vessel_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in perfused_vessel_ranks]  # get vessel length based on rank    
    fcd_numerator = sum(perfused_vessel_lengths)

    # Create grid based on reference network 
    grid_lines, total_line_length = create_grid_lines(field_data_um)

    # Round the coordinates to allow for some artefacts
    grid_lines = np.round(grid_lines, decimals=0)
    network_lines = np.round(segment_xyxy, decimals=0)
    
    # Compute two possible numerators of vessel density
    # intersection_points = find_intersections(network_lines, grid_lines, unique=0) # number of intersections between the grid and the network
    # intersecting_lines = find_intersections(network_lines, grid_lines, unique=1)  # number of vessels that intersect with at least one line
    intersection_points = find_intersections(network_lines, grid_lines)
    intersecting_lines = {entry['network_line'] for entry in intersection_points}
    
    # Compute two possible definitions of vessel density
    vd_intersections = len(intersection_points)/total_line_length
    vd_uniques = len(intersecting_lines)/total_line_length

    # Calculate another FCD definition
    
    # Calculate the Heterogeneity Index (for the flow velocity)

    return total_blood_volume, fcd_numerator, vd_intersections, vd_uniques

# Define a function to return the following perfusion metrics: Total Blood Vlow (in m3/s); FCD (functional capillary density); Vessel Density; Heterogeneity Index (for the flow velocity)
def get_perfusion_metrics_voronoi(vessel_data_array, pq_threshold, segment_lengths_um, reference_network_coordinates, segment_xyxy, field_data_um):
    
    # Calculate the TBV (total blood volume) (in m3/s) (NOTE: This is broken. I need to calculate the volume of the vessel but I don't know the definition of TBV really.)
    '''
    first_rank_zero_index = np.where(vessel_data_array['Vessel Owner Rank'] == 0)[0][0]  # Gets the index of the first inlet/outlet vessel in the array (assumes all inlets/outlets are equal)
    first_rank_zero_index = 0  # Gets the index of the first inlet/outlet vessel in the array (assumes all inlets/outlets are equal)
    flowing_vessel_indices = np.where(vessel_data_array['Absolute Vessel Flow Rate m^3/s'] > 0.0)[0] 
    vessel_volume = 
    # total_blood_volume = vessel_data_array['Absolute Vessel Flow Rate m^3/s'][first_rank_zero_index]  # Get the flow rate for the rank 0 vessel
    '''
    total_blood_volume = np.nan
    
    # Calculate the FCD (functional capillary density)
    # perfused_vessel_indices = np.where(vessel_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold)[0] 
    # perfused_vessel_ranks = vessel_data_array['Vessel Owner Rank'][perfused_vessel_indices]
    # perfused_vessel_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in perfused_vessel_ranks]  # get vessel length based on rank    
    # fcd_numerator = sum(perfused_vessel_lengths)

    # Calculate the FCD (functional capillary density)
    perfused_vessel_indices = np.where(vessel_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold)[0] 
    
    # perfused_vessel_ranks = vessel_data_array['Vessel Owner Rank'][perfused_vessel_indices]
    
    # perfused_vessel_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in perfused_vessel_ranks]  # get vessel length based on rank    
    perfused_vessel_lengths_um = segment_lengths_um[perfused_vessel_indices]
    fcd_numerator = sum(perfused_vessel_lengths_um)

    # print(segment_lengths_um)

    # Create grid based on reference network 
    grid_lines, total_line_length = create_grid_lines(field_data_um)

    # Round the coordinates to allow for some artefacts
    grid_lines = np.round(grid_lines, decimals=0)
    network_lines = np.round(segment_xyxy, decimals=0)
    
    # Compute two possible numerators of vessel density
    intersection_points = find_intersections(network_lines, grid_lines)
    intersecting_lines = {entry['network_line'] for entry in intersection_points}

    # intersecting_lines = find_intersections(network_lines, grid_lines, unique=1)  # number of vessels that intersect with at least one line
    
    # Compute two possible definitions of vessel density
    vd_intersections = len(intersection_points)/total_line_length
    vd_uniques = len(intersecting_lines)/total_line_length

    # Calculate another FCD definition
    
    # Calculate the Heterogeneity Index (for the flow velocity)

    return total_blood_volume, fcd_numerator, vd_intersections, vd_uniques

# Acquire the grid lines based on the PDE grid under evaluation
def create_grid_lines(field_data_um):
    
    # Get the network extents and note the shorter length to create a square grid
    # print(np.max(reference_network_coordinates))
    # x_extent = np.maximum(np.max(reference_network_coordinates[:, 0]), np.max(reference_network_coordinates[:, 2]))  # 
    # y_extent = np.maximum(np.max(reference_network_coordinates[:, 1]), np.max(reference_network_coordinates[:, 3]))
    x_extent = max(field_data_um['x'])
    y_extent = max(field_data_um['y'])
    shortest_axis = np.minimum(x_extent, y_extent)
    centre_x, centre_y = np.round(x_extent/2,decimals=0), np.round(y_extent/2,decimals=0)
    grid_unit = np.round(shortest_axis/4,decimals=0) # side length of each component square 
    vline_1 = [centre_x-grid_unit, centre_y+2*grid_unit, centre_x-grid_unit, centre_y-2*grid_unit]
    vline_2 = [centre_x, centre_y+2*grid_unit, centre_x, centre_y-2*grid_unit]
    vline_3 = [centre_x+grid_unit, centre_y+2*grid_unit, centre_x+grid_unit, centre_y-2*grid_unit]
    hline_1 = [centre_x-2*grid_unit, centre_y-grid_unit, centre_x+2*grid_unit, centre_y-grid_unit]
    hline_2 = [centre_x-2*grid_unit, centre_y, centre_x+2*grid_unit, centre_y]
    hline_3 = [centre_x-2*grid_unit, centre_y+grid_unit, centre_x+2*grid_unit, centre_y+grid_unit]
    grid_lines = np.vstack([vline_1, vline_2, vline_3, hline_1, hline_2, hline_3])
    total_line_length = grid_unit*4*6  # denominator for VD calculations 
    # print(grid_lines)
    return grid_lines, total_line_length

# Define a function to compute the vessel density numerator 
def find_intersections(network_lines, grid_lines):
    
    # Create an empty array for the number of intersections
    intersections = []
    
    # If we want to count the number of unique lines, create a new variable to store the lines
    # if unique == 1:
    #     processed_lines = set()

    # For each network line i
    for i, network_line in enumerate(network_lines):
        # print(i)
        
        # Skip lines from the network that have already had an intersection (if we only want unique lines)
        # if unique == 1:
        #     if i in processed_lines:
        #         continue
        
        # If we haven't considered the network line alredaya, then for each grid line j        
        for j, grid_line in enumerate(grid_lines):
            # print(j)
            # Get the start and end nodes
            x1, y1, x2, y2 = network_line
            x3, y3, x4, y4 = grid_line

            # Check for intersection
            denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
            
            # If lines are parallel or coincidental, no unique intersection point so continue to the next grid line
            if denominator == 0:
                # print('denom 0')
                continue
            
            ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator
            ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

            if 0 <= ua <= 1 and 0 <= ub <= 1:
                # Intersection point found
                intersection_point = (x1 + ua * (x2 - x1), y1 + ua * (y2 - y1))
                intersections.append({
                    'network_line': i,
                    'grid_line': j,
                    'point': intersection_point
                })
                # print(i)
                # if unique == 1:
                #     processed_lines.add(i)
                # break  # exit the inner loop after the first intersection for the current line
                # if intersections:
                # else:
                    
        # print(unique)
        # print(intersections)
        
    return intersections

# =============================================================================
# TDA METRICS
# =============================================================================

'''
# Define a function to get the network metrics
def get_network_metrics(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes):

    G_filtered, inlet_nodes_graph, outlet_nodes_graph = filter_graph(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes)    

    # Initialize a dictionary to store the metrics
    metrics = {}
  
    # 2. Betweenness Centrality
    # metrics['edge_betweenness_unweighted'] = nx.edge_betweenness_centrality(G_filtered)
    metrics['edge_betweenness_weighted'] = nx.edge_betweenness_centrality(G_filtered, weight='length')

    # Calculate node betweenness centrality 
    # metrics['betweenness_unweighted'] = nx.betweenness_centrality(G_filtered)
    # metrics['edge_betweenness_max_unweighted'] = {(u, v): max(metrics['betweenness_unweighted'][u], metrics['betweenness_unweighted'][v]) for u, v in G_filtered.edges()}    
    metrics['betweenness_weighted'] = nx.betweenness_centrality(G_filtered, weight='length')
    # metrics['edge_betweenness_max_weighted'] = {(u, v): max(metrics['betweenness_weighted'][u], metrics['betweenness_weighted'][v]) for u, v in G_filtered.edges()}    

    # 4. PageRank
    metrics['pagerank_unweighted'] = nx.pagerank(G_filtered)
    metrics['edge_pagerank_unweighted'] = {(u, v): (metrics['pagerank_unweighted'][u] + metrics['pagerank_unweighted'][v]) / 2 for u, v in G_filtered.edges()}
    metrics['edge_pagerank_max_unweighted'] = {(u, v): max(metrics['pagerank_unweighted'][u], metrics['pagerank_unweighted'][v]) for u, v in G_filtered.edges()}
    # metrics['pagerank_weighted'] = nx.pagerank(G_filtered, weight='length')
    # metrics['edge_pagerank_weighted'] = {(u, v): (metrics['pagerank_weighted'][u] + metrics['pagerank_weighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_pagerank_max_weighted'] = {(u, v): max(metrics['pagerank_weighted'][u], metrics['pagerank_weighted'][v]) for u, v in G_filtered.edges()}

    # 5 & 6. Distance to Inlet and Outlet (based on original indices)
    inlet_node = inlet_nodes_graph
    outlet_node = outlet_nodes_graph
    lengths_to_inlet = nx.shortest_path_length(G_filtered, source=inlet_node, weight='length')
    lengths_to_outlet = nx.shortest_path_length(G_filtered, source=outlet_node, weight='length')    
    metrics['edge_distance_to_inlet'] = {(u, v): min(lengths_to_inlet[u], lengths_to_inlet[v]) for u, v in G_filtered.edges()}
    metrics['edge_distance_to_outlet'] = {(u, v): min(lengths_to_outlet[u], lengths_to_outlet[v]) for u, v in G_filtered.edges()}

    return metrics, G_filtered
'''

# Define a function to get the network metrics
def get_network_metrics(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes):

    G_filtered, _, _  = filter_graph(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes)    

    # Initialize a dictionary to store the metrics
    metrics = {}
  
    # 1. Closeness Centrality (underwhelming results)
    # metrics['closeness_unweighted'] = nx.closeness_centrality(G_filtered)
    # metrics['edge_closeness_unweighted'] = {(u, v): (metrics['closeness_unweighted'][u] + metrics['closeness_unweighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_closeness_max_unweighted'] = {(u, v): max(metrics['closeness_unweighted'][u], metrics['closeness_unweighted'][v]) for u, v in G_filtered.edges()}    
    # metrics['closeness_weighted'] = nx.closeness_centrality(G_filtered, distance='length')
    # metrics['edge_closeness_weighted'] = {(u, v): (metrics['closeness_weighted'][u] + metrics['closeness_weighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_closeness_max_weighted'] = {(u, v): max(metrics['closeness_weighted'][u], metrics['closeness_weighted'][v]) for u, v in G_filtered.edges()}
    
    # 2. Betweenness Centrality
    # metrics['edge_betweenness_unweighted'] = nx.edge_betweenness_centrality(G_filtered)
    # metrics['edge_betweenness_weighted'] = nx.edge_betweenness_centrality(G_filtered, weight='length')

    # Calculate node betweenness centrality 
    metrics['betweenness_unweighted'] = nx.betweenness_centrality(G_filtered)
    # metrics['edge_betweenness_max_unweighted'] = {(u, v): max(metrics['betweenness_unweighted'][u], metrics['betweenness_unweighted'][v]) for u, v in G_filtered.edges()}    
    # metrics['betweenness_weighted'] = nx.betweenness_centrality(G_filtered, weight='length')
    # metrics['edge_betweenness_max_weighted'] = {(u, v): max(metrics['betweenness_weighted'][u], metrics['betweenness_weighted'][v]) for u, v in G_filtered.edges()}    
    
    # 3. Eigenvector Centrality (takes too long and not very predictive of anything)
    # metrics['eigenvector_unweighted'] = nx.eigenvector_centrality(G_filtered, max_iter=10000)
    # metrics['edge_eigenvector_unweighted'] = {(u, v): (metrics['eigenvector_unweighted'][u] + metrics['eigenvector_unweighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_eigenvector_max_unweighted'] = {(u, v): max(metrics['eigenvector_unweighted'][u], metrics['eigenvector_unweighted'][v]) for u, v in G_filtered.edges()}
    # metrics['eigenvector_weighted'] = nx.eigenvector_centrality(G_filtered, max_iter=10000, weight='length')
    # metrics['edge_eigenvector_weighted'] = {(u, v): (metrics['eigenvector_weighted'][u] + metrics['eigenvector_weighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_eigenvector_max_weighted'] = {(u, v): max(metrics['eigenvector_weighted'][u], metrics['eigenvector_weighted'][v]) for u, v in G_filtered.edges()}
    
    # 4. PageRank
    metrics['pagerank_unweighted'] = nx.pagerank(G_filtered)
    # metrics['edge_pagerank_unweighted'] = {(u, v): (metrics['pagerank_unweighted'][u] + metrics['pagerank_unweighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_pagerank_max_unweighted'] = {(u, v): max(metrics['pagerank_unweighted'][u], metrics['pagerank_unweighted'][v]) for u, v in G_filtered.edges()}
    # metrics['pagerank_weighted'] = nx.pagerank(G_filtered, weight='length')
    # metrics['edge_pagerank_weighted'] = {(u, v): (metrics['pagerank_weighted'][u] + metrics['pagerank_weighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_pagerank_max_weighted'] = {(u, v): max(metrics['pagerank_weighted'][u], metrics['pagerank_weighted'][v]) for u, v in G_filtered.edges()}
    
    # 5 & 6. Distance to Inlet and Outlet
    # inlet_node = inlet_nodes_graph
    # outlet_node = outlet_nodes_graph
    # lengths_to_inlet = nx.shortest_path_length(G_filtered, source=inlet_node, weight='length')
    # lengths_to_outlet = nx.shortest_path_length(G_filtered, source=outlet_node, weight='length')    
    # metrics['edge_distance_to_inlet'] = {(u, v): min(lengths_to_inlet[u], lengths_to_inlet[v]) for u, v in G_filtered.edges()}
    # metrics['edge_distance_to_outlet'] = {(u, v): min(lengths_to_outlet[u], lengths_to_outlet[v]) for u, v in G_filtered.edges()}

    # 7. Clustering Coefficient (not very predictive)
    # metrics['clustering_unweighted'] = nx.clustering(G_filtered)
    # metrics['edge_clustering_unweighted'] = {(u, v): (metrics['clustering_unweighted'][u] + metrics['clustering_unweighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_clustering_max_unweighted'] = {(u, v): max(metrics['clustering_unweighted'][u], metrics['clustering_unweighted'][v]) for u, v in G_filtered.edges()}
    # metrics['clustering_weighted'] = nx.clustering(G_filtered, weight='length')
    # metrics['edge_clustering_weighted'] = {(u, v): (metrics['clustering_weighted'][u] + metrics['clustering_weighted'][v]) / 2 for u, v in G_filtered.edges()}
    # metrics['edge_clustering_max_weighted'] = {(u, v): max(metrics['clustering_weighted'][u], metrics['clustering_weighted'][v]) for u, v in G_filtered.edges()}
    
    # 8. Average Path Length (underwhelming results)
    # metrics['average_path_length_unweighted'] = nx.average_shortest_path_length(G_filtered)
    # metrics['average_path_length_weighted'] = nx.average_shortest_path_length(G_filtered, weight='length')
    
    return metrics, G_filtered

def filter_graph(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes):
    
    # Initialize adjacency matrices
    adj_matrix = np.array(length_adjacency_matrix)
    flow_rate_matrix = np.zeros_like(adj_matrix)
    haematocrit_matrix = np.zeros_like(adj_matrix)
    length_matrix = np.zeros_like(adj_matrix)
    id_matrix = np.zeros_like(adj_matrix)
    max_possible_nodes = adj_matrix.shape[0]

    # Assume first node is inlet and last node is outlet
    inlet_node_index = 0
    outlet_node_index = max_possible_nodes-1
    
    # Fill the attribute matrices directly using node indices from segment_nodes
    for segment_index, (node_1, node_2) in enumerate(segment_nodes):
        flow_rate_matrix[node_1, node_2] = cell_data_array['Absolute Vessel Flow Rate m^3/s'][segment_index]
        haematocrit_matrix[node_1, node_2] = cell_data_array['Vessel Haematocrit'][segment_index]
        length_matrix[node_1, node_2] = cell_data_array['Vessel Length'][segment_index] * 1e6  # Convert to micrometers

    # Mirror the values in the matrices for undirected graph
    flow_rate_matrix = np.maximum(flow_rate_matrix, flow_rate_matrix.T)
    haematocrit_matrix = np.maximum(haematocrit_matrix, haematocrit_matrix.T)
    length_matrix = np.maximum(length_matrix, length_matrix.T)
    id_matrix = np.maximum(id_matrix, id_matrix.T)

    # Create the graph with the original node indices before filtering disconnected nodes
    G = nx.Graph()

    # Sort the node pairs in reference_network_vessel_nodes for undirected edge handling
    sorted_reference_network_vessel_nodes = [tuple(sorted(edge)) for edge in reference_network_vessel_nodes]

    # Add edges with attributes, keeping the original node indices
    for i in range(max_possible_nodes):
        for j in range(i + 1, max_possible_nodes):
            if adj_matrix[i, j] != 0:  # Only consider existing edges
                
                # Sort the nodes to handle undirected edges
                sorted_edge = tuple(sorted([i, j]))
                
                # Find the index of the edge in the reference network
                if sorted_edge in sorted_reference_network_vessel_nodes:
                    original_ID = sorted_reference_network_vessel_nodes.index(sorted_edge)
                else:
                    raise ValueError(f"Edge {sorted_edge} not found in reference network.")
    
                # Add edge with all attributes, including original_ID
                G.add_edge(i, j, 
                           abs_flow_rate=flow_rate_matrix[i, j], 
                           haematocrit=haematocrit_matrix[i, j],
                           original_ID=original_ID,
                           length=length_matrix[i, j])
    
    # Add node attributes
    for node in G.nodes():
        G.nodes[node]['original_node_index'] = node  # Store the original node index
        G.nodes[node]['is_inlet'] = (node == inlet_node_index)
        G.nodes[node]['is_outlet'] = (node == outlet_node_index)

    # Check for duplicate edges in cell_data_array
    print("Checking for duplicate edges in cell_data_array...")
    seen_edges = set()
    duplicate_edges = []
    for segment_index, edge in enumerate(segment_nodes):
        sorted_edge = tuple(sorted(edge))
        if sorted_edge in seen_edges:
            duplicate_edges.append((sorted_edge, segment_index))
        else:
            seen_edges.add(sorted_edge)
    
    if duplicate_edges:
        print(f"Duplicate edges in cell_data_array ({len(duplicate_edges)}):")
        for edge, index in duplicate_edges:
            print(f"Edge {edge} at segment_nodes index {index}")
    else:
        print("No duplicate edges found in cell_data_array.")


    # # Inspect attributes for edge (16, 17)
    # target_edge = (16, 17)
    # print(f"Inspecting attributes for edge {target_edge}:")
    # for segment_index, edge in enumerate(segment_nodes):
    #     if tuple(sorted(edge)) == target_edge:
    #         print(f"Index {segment_index}:")
    #         for key in cell_data_array.keys():
    #             print(f"  {key}: {cell_data_array[key][segment_index]}")

    # Check that the number of edges in the graph matches the cell data
    num_edges_in_graph = G.number_of_edges()
    num_edges_in_cell_data = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'])
    assert num_edges_in_graph == num_edges_in_cell_data, \
        f"Mismatch: Graph has {num_edges_in_graph} edges, but cell_data_array has {num_edges_in_cell_data} edges."

    # Remove isolated nodes (nodes with no connections)
    G.remove_nodes_from(list(nx.isolates(G)))  # This removes all isolated nodes

    # Now filter out any disconnected components (with no flow) while preserving original node numbering
    '''
    connected_components = list(nx.connected_components(G))
    valid_subgraph = []
    for component in connected_components:
        subgraph = G.subgraph(component).copy()
        total_flow_rate = sum(data.get('abs_flow_rate', 0) for u, v, data in subgraph.edges(data=True))
        if total_flow_rate > 0:
            valid_subgraph.append(subgraph)
    if len(valid_subgraph) > 1:
        raise ValueError("Error: The network has more than one connected component with flow!")
    elif len(valid_subgraph) < 1:
        raise ValueError("Error: The network has no connected components with flow!")        
    else:
        G_filtered = valid_subgraph[0]
    '''
        
    # Now retain only the connected component with the most edges (which should ideally be the one with flow) while preserving original node numbering
    connected_components = list(nx.connected_components(G))
    largest_component = max(connected_components, key=lambda comp: G.subgraph(comp).number_of_edges())  # find the largest connected component based on the number of edges
    G_filtered = G.subgraph(largest_component).copy()

    # Since you know which nodes are the inlet and outlet, you can set them directly
    inlet_nodes_graph = [n for n in G_filtered.nodes() if G_filtered.nodes[n]['is_inlet']]
    outlet_nodes_graph = [n for n in G_filtered.nodes() if G_filtered.nodes[n]['is_outlet']]

    # Ensure that exactly a maximum of one inlet and one outlet exist
    assert len(inlet_nodes_graph) <= 1, f"Expected <=1 inlet node, but found {len(inlet_nodes_graph)}"
    assert len(outlet_nodes_graph) <= 1, f"Expected <=1 outlet node, but found {len(outlet_nodes_graph)}"

    return G_filtered, inlet_nodes_graph, outlet_nodes_graph

# Define a function to return the Betti curves for a network (based on Bernadette Stolz's code)
def PH_betti2(Weighted_adjacency_matrix):
    
#    beta0=np.array([])
#    beta1=np.array([])
#    size_largest_connected_component = np.array([])
#    biggest1 = np.array([])
    
#    for threshold in thresholds:
        
#    print('Threshold = ',step)
    adj = np.array(Weighted_adjacency_matrix)
#    adj[adj <= step] = 0
    
    # we don't want any nodes showing up as separate connected components
    adj = adj[~np.all(adj == 0, axis=1), :] #rows
    adj = adj[:,~np.all(adj == 0, axis=0)] #cols - logically axis should be =1 here, but code only works if I set it to 0

    n_nodes = adj.shape[0]
    
    G = nx.Graph(adj)
    number_of_connected_components = nx.number_connected_components(G)
    
#    beta0=np.append(beta0, number_of_connected_components)
    
    H = max(nx.connected_components(G), key=len)
    H_subgraph = G.subgraph(H).copy()
    size_largest_connected_component = H_subgraph.number_of_nodes()
#        size_largest_connected_component = np.append(size_largest_connected_component, size_largest_cc)

    #computes Betti-1
    n_edges = G.number_of_edges()
    n_cycles = number_of_connected_components - n_nodes + n_edges
#        beta1=np.append(beta1, n_cycles)
    
    
    #computes Beti-1 of largest connected component
    n_cycles_largest_connected_component = 1 - size_largest_connected_component + H_subgraph.number_of_edges()
#        biggest1 = np.append(biggest1,n_cycles_largest_connected_component)
    
    return number_of_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component

# =============================================================================
# HAEMATOCRIT METRICS
# =============================================================================

# Gini Coefficient
def gini_coefficient_optimized(x):
    x = np.array(x).flatten()
    if np.any(x < 0):
        raise ValueError("All values must be non-negative.")
    if np.all(x == 0):
        return 0.0
    x = np.sort(x)
    n = x.size
    index = np.arange(1, n + 1)
    Gini = (2 * np.sum(index * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))
    return Gini

# Calculate Inverse Simpson Index
# def inverse_simpson_index(oxygen_levels):
    
#     # Define bins for oxygen levels
#     # bins = [0, 10, 50, 100]  # Example thresholds
#     n_bins = 1000
    
#     # Categorize data
#     counts, _ = np.histogram(oxygen_levels, bins=n_bins)

#     proportions = counts / np.sum(counts)
#     # Remove zero proportions to avoid division by zero
#     # proportions = proportions[proportions > 0]
#     return 1 / np.sum(proportions**2)

# Inverse Simpson’s Index
def inverse_simpson_index(x):
    x = np.array(x)
    total = np.sum(x)
    proportions = x / total
    return 1 / np.sum(proportions**2)

# Hoover Index
def hoover_index(x):
    mean_x = np.mean(x)
    return 0.5 * np.sum(np.abs(x - mean_x)) / np.sum(x)

# Haemonormalisation Potential (HNP)
def haemonormalisation_potential(x, h_target):
    deviations = np.abs(x - h_target)
    return np.mean(deviations)

# =============================================================================
# ARCHIVE
# =============================================================================

def atkinson_index(x, epsilon=0.5):
    x = np.array(x)
    mean_x = np.mean(x)
    if epsilon == 1:
        A = 1 - np.exp(np.mean(np.log(x[x > 0]))) / mean_x
    else:
        A = 1 - (np.mean(x**(1 - epsilon))**(1 / (1 - epsilon))) / mean_x
    return A

def mean_log_deviation(x):
    x = np.array(x)
    mean_x = np.mean(x)
    small_value = 1e-10  # Replace zeros with a small value
    x = np.where(x == 0, small_value, x)  # Replace zeros
    log_term = np.log(mean_x / x)
    return np.mean(log_term)

def theil_index(x):
    x = np.array(x)
    mean_x = np.mean(x)
    small_value = 1e-10  # Replace zeros with a small value
    x = np.where(x == 0, small_value, x)  # Replace zeros
    log_term = np.log(x / mean_x)
    return np.sum(x * log_term) / np.sum(x)

def shannon_entropy(data, bins=5):
    # Create bins and calculate histogram
    counts, _ = np.histogram(data, bins=bins, density=True)
    
    # Convert counts to probabilities
    probabilities = counts[counts > 0]  # Ignore zero probabilities
    probabilities /= probabilities.sum()  # Normalize to get probabilities
    
    # Calculate Shannon entropy
    entropy = -np.sum(probabilities * np.log(probabilities))
    return entropy

def pietra_index(x):
    x = np.array(x)
    mean_x = np.mean(x)
    return 0.5 * np.sum(np.abs(x - mean_x)) / np.sum(x)