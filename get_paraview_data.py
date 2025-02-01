#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 31 16:44:36 2022

@author: narain

Read VTK, VTU, and VTI files to process them in Python.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
#from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
from get_vascular_metrics import *  # import functions for vascular metrics
# import math
# import networkx as nx
import numpy as np
import os
import pandas as pd
import pyvista as pv
import scipy.interpolate
import vtk
import xml.etree.ElementTree as ET
# from scipy.stats import skew, kurtosis

# Set up plotting
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vtk file and return the point coordinates and data and the cell data
def get_vtk_data(vtk_path):
    
    # print(vtk_path)

    # Set up the file reader
#    if 'vtp' in filename:
#        reader = vtk.vtkXMLPolyDataReader()
#    else:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.Update()
    polydata = reader.GetOutput()
#    polydata.ReleaseDataFlagOn()
    
    # Extract the point coordinates
    point_coordinates_data = polydata.GetPoints().GetData()
    point_coordinates_data_array = vtk_to_numpy(point_coordinates_data)
    
    # Extract the point data 
    point_data_array = {}
    point_data = polydata.GetPointData()
    for column_index in range(point_data.GetNumberOfArrays()):
       column_name =  point_data.GetArrayName(column_index)
       point_data_array[column_name] = vtk_to_numpy(point_data.GetArray(column_index))
    
    # Extract the cell data    
    cell_data_array = {}
    cell_data = polydata.GetCellData()
    for column_index in range(cell_data.GetNumberOfArrays()):
       column_name =  cell_data.GetArrayName(column_index)
       cell_data_array[column_name] = vtk_to_numpy(cell_data.GetArray(column_index))
    
    # Return a dictionary with the point coordinates and data and cell data
    return point_coordinates_data_array, point_data_array, cell_data_array, polydata
           
# Define a function to read a .vti file and return the data
def get_vti_data(vti_path):
    
    # Import the file
    extension = vti_path.split('.').pop()
    reader = None
    if extension == 'vtk':
        reader = vtk.vtkDataSetReader() 
    elif extension == 'vti':
        reader = vtk.vtkXMLImageDataReader() 
    else:
        raise RuntimeError('Unknown File Type: %s ' % vti_path)
    reader.SetFileName( "%s" % (vti_path) ) 
    reader.Update()
    image_data = reader.GetOutput()
    
    # Extract the dimensions, spacing, and origin
    spacing = image_data.GetSpacing()
    
    # Extract the point values 
    field_point_data = image_data.GetPointData() 
    field_values = vtk_to_numpy(field_point_data.GetArray(0)) 
    
    # Get the coordinates of each point
    position_list = deque()
    for index in range(len(field_values)):  # do something 
        position = image_data.GetPoint(index)
        position_list.append(position)
    position_array = np.array(position_list)
    
    # Return the field distribution
    distribution_array = np.column_stack((position_array,field_values))
    distribution = pd.DataFrame(data=distribution_array[:,[0,1,3]], columns=["x", "y", "oxygen"])
    return distribution, spacing#, dimensions, origin       

# Define a function to read a .pvd file and return the time series data using all associated .vtu files in the same folder
def get_pvd_data(pvd_path):
    
    # Load and parse the .pvd file
    pvd_filename = pvd_path
    tree = ET.parse(pvd_filename)
    root = tree.getroot()
    
    # Get the directory of the .pvd file
    pvd_directory = os.path.dirname(pvd_filename)
    
    # Initialize a MultiBlock dataset
    multiblock = pv.MultiBlock()
    
    # Extract .vtu files from the .pvd file and load each one
    for dataset in root.iter('DataSet'):
        vtu_file = dataset.get('file')
        
        # Construct the full path to the .vtu file
        vtu_full_path = os.path.join(pvd_directory, vtu_file)    
        # print(f"Loading {vtu_full_path}...")
        vtu_mesh = pv.read(vtu_full_path)
        multiblock.append(vtu_mesh)
    
    # Check how many blocks were added
    # print(f"Number of blocks loaded: {len(multiblock)}")
    
    '''
    # Now we can access each block as before
    for i, block in enumerate(multiblock):
        print(f"Block {i}:")
        print(block)
    
        if isinstance(block, pv.UnstructuredGrid):
            print(f"Block {i} - Number of Points: {block.n_points}")
            print(f"Block {i} - Number of Cells: {block.n_cells}")
            print(f"Block {i} - Point Data Keys: {block.point_data.keys()}")
            print(f"Block {i} - Cell Data Keys: {block.cell_data.keys()}")
    '''
    
    # Initialize a list to hold DataFrames for each block
    df_list = []
    
    # Iterate over each block in the MultiBlock dataset
    for i, block in enumerate(multiblock):
        if isinstance(block, pv.UnstructuredGrid):
            # Initialize a dictionary to hold the arrays for this block
            data = {}
            
            # Iterate over all arrays in point_data for this block
            for name in block.point_data.keys():
                data[name] = block.point_data[name]
            
            # Convert the dictionary to a DataFrame
            df_block = pd.DataFrame(data)
            
            # Add an identifier for the block (optional)
            df_block['Block_ID'] = i
            
            # Append the DataFrame to the list
            df_list.append(df_block)
    
    # Concatenate all DataFrames into a single DataFrame
    if df_list:
        df = pd.concat(df_list, ignore_index=True)
        # print(df)
    
        # Optionally, save the combined DataFrame to a CSV file
        # df.to_csv('output_all_blocks.csv', index=False)
    else:
        print("No data available in the .pvd file.")

    return df

# Define a function to convert the .vti data to a plottable field with some basic stats
def get_plottable_field(vti_data):
    
    # Change the data type
    vti_data = vti_data.astype('float32')
    
    # Calculate concentration statistics
    O2_stats = vti_data['oxygen'].describe()
    
    # Downsample data if needed
    vti_data = vti_data[::int(1)]
    
    # Convert dataframe into NumPy matrix
    mat = vti_data.to_numpy()
    
    # Get the x and y axes
    x = np.unique(mat[:,0])
    y = np.unique(mat[:,1])
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (X,Y), method='nearest')

    # Return the oxygen stats
    return Z, O2_stats

# Define a function to extract the region of evaluation of the forking network
def get_forking_domain(field_dataset, middle_generation_number, generation_coordinates):
    x_start = generation_coordinates[middle_generation_number-1]
    x_end = generation_coordinates[-middle_generation_number]
    middle_field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
    return middle_field_dataset

# Define a function to extract the region of evaluation of the hexagonal network
def get_hex_domain(field_dataset, field_spacing):
    x_start = 10*field_spacing[0]
    x_end = 195*field_spacing[0]
    y_start = 17*field_spacing[1]
    y_end = 173*field_spacing[1]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to extract a custom region of evaluation
def get_custom_domain(field_dataset, field_spacing, x_start, x_end, y_start, y_end):
    x_start = x_start*field_spacing[0]
    x_end = x_end*field_spacing[0]
    y_start = y_start*field_spacing[1]
    y_end = y_end*field_spacing[1]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x
# '''
# Specify the vessels in terms of node IDs and node coordinates
def get_vessel_nodes(polydata, reference_node_coordinates):
    
    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell IDs store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])  # array to store node IDs
    segment_xyxy = np.array([])  # array to store node coordinates 
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        node_xy = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
            # segment_xxyy = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
            node_xy = np.hstack([node_xy, x, y])
        segment_xyxy = np.vstack([segment_xyxy, node_xy]) if segment_xyxy.size else node_xy
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')

    return segment_nodes, segment_xyxy
# '''

# Define a function to combine all the .csv files containing all the required metrics into one giant .csv file
def combine_csv_files(main_folder_path, max_layout, file_name):
    dfs = []
    for which_layout in range(1, max_layout+1):
        df = pd.read_csv(main_folder_path + '/Selection' + str(which_layout) + '/' + file_name)
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)  # concatenate all DataFrames into a single DataFrame
    combined_df.to_csv(main_folder_path + '/' + file_name, index=False)  # index=False prevents saving the index as a separate column

# Define a function to combine all the HDF5 files containing all the required metrics into one DataFrame
def combine_hdf5_files(main_folder_path, max_layout, hdf_file_name, key):
    dfs = []
    for which_layout in range(1, max_layout+1):
        # Construct the file path for each HDF5 file
        file_path = f'{main_folder_path}/Selection{which_layout}/{hdf_file_name}'
        # Open the HDF5 file and read the specific DataFrame using the key
        with pd.HDFStore(file_path) as store:
            df = store[key]  # Retrieve the DataFrame using the given key
            dfs.append(df)
    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

# Define a function to extract the predictive metrics from the forking network
def get_forking_predictors(vtk_path, reference_rank_lengths, reference_node_coordinates, reference_network_coordinates, pq_threshold, field_data_um):
        
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    reference_rank_lengths_um = reference_rank_lengths*1000000  # convert length (m to um)
    
    # Get the architectural metrics
    n_vessels = len(cell_data_array['Vessel Radius m'])
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    mean_diameter = sum(segment_diameters_um)/n_vessels
    segment_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in cell_data_array['Vessel Owner Rank']]  # get vessel length based on rank    
    mean_geometric_resistance = sum(segment_lengths/(segment_diameters_um**4))/n_vessels

    # Get a list of vessel segments with node IDs for the adjacency matrices
    segment_nodes, segment_xyxy = get_vessel_nodes(polydata, reference_node_coordinates)
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]
        length_adjacency_matrix[col,row] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]

    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(diameter_adjacency_matrix)

    # Get the perfusion metrics
    total_blood_volume, fcd_numerator, vd_intersections, vd_uniques = get_perfusion_metrics_forking(cell_data_array, pq_threshold, reference_rank_lengths_um, reference_network_coordinates, segment_xyxy, field_data_um) 

    # Return the architectural features and the adjacency matrices
    # return n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component
    return n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, n_cycles, total_blood_volume, fcd_numerator, vd_intersections, vd_uniques

'''
# Define a function to return the details of the vessel that was pruned in a certain simulation
def log_killed_vessel(network_path, previous_network_path, previous_edge_metrics_df):
    
    # Get the data for the current network
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(network_path)

    # Round the coordinates to two decimal places
    point_coordinates_data_array = point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    node_coordinates = np.around(point_coordinates_data_array[:,0:2], 2)

    # Express the vessels in terms of node coordinates
    _, vessel_coordinates = get_vessel_nodes(polydata, node_coordinates)
    
    # Repeat for the preceding network
    prev_point_coordinates_data_array, _, prev_cell_data_array, prev_polydata = get_vtk_data(previous_network_path)
    prev_point_coordinates_data_array = prev_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    prev_node_coordinates = np.around(prev_point_coordinates_data_array[:,0:2], 2)
    prev_vessel_nodes, prev_vessel_coordinates = get_vessel_nodes(prev_polydata, prev_node_coordinates)
            
    # Find the edge that was deleted
    prev_edges_set = set(map(tuple, prev_vessel_coordinates))
    current_edges_set = set(map(tuple, vessel_coordinates))

    # Find the missing edge(s)
    missing_edges = list(prev_edges_set - current_edges_set)
    # print(missing_edges)

    # Check if no missing edges were found
    if len(missing_edges) == 0:
        raise ValueError(f"No missing edge detected between the previous network and {network_path}. This is unexpected.")
    
    # Check if more than one missing edge was found
    if len(missing_edges) > 1:
        raise ValueError(f"More than one missing edge detected between the previous network and {network_path}. Multiple edges found: {missing_edges}. This is unexpected.")

    # Process the first (and only) missing edge
    missing_edge = missing_edges[0]
    # missing_edge = list(prev_edges_set - current_edges_set)[0]  # Get the first (and only) missing edge
    index_in_prev = np.where((prev_vessel_coordinates == missing_edge).all(axis=1))[0][0]
    
    # Store the data for the missing edge
    killed_vessel_data = {key: value[index_in_prev] for key, value in prev_cell_data_array.items()}
    df_killed_vessel = pd.DataFrame([killed_vessel_data])

    # Extract the node coordinates of the deleted vessel
    x1, y1, x2, y2 = prev_vessel_coordinates[index_in_prev]
    print(x1, y1, x2, y2)
    # Add the coordinates as new columns
    # df_killed_vessel['x1'] = x1#
    # df_killed_vessel['y1'] = y1#
    # df_killed_vessel['x2'] = x2#
    # df_killed_vessel['y2'] = y2#

    # Convert missing_edge coordinates (x1, y1, x2, y2) to node tuples
    node_1 = np.where((prev_node_coordinates == (x1, y1)).all(axis=1))[0][0]
    node_2 = np.where((prev_node_coordinates == (x2, y2)).all(axis=1))[0][0]
    
    # Create the node tuple representing the edge
    missing_edge_as_nodes = tuple(sorted([node_1, node_2]))  # Sort to handle undirected edges

    # Find the metrics of the missing edge from the previous network's edge_metrics_df
    edge_metrics = previous_edge_metrics_df[previous_edge_metrics_df['edge'] == missing_edge_as_nodes].to_dict('records')
    # if edge_metrics:
    edge_metrics = edge_metrics[0]  # Assuming one row is returned
    print('zazu',edge_metrics['edge_distance_to_outlet'])
    for metric, value in edge_metrics.items():
        if metric != 'edge':  # Skip the 'edge' column itself
            df_killed_vessel[metric] = value
            
    return df_killed_vessel
'''

# Define a function to find the index of the missing edge as numbered in the reference network (i.e., Vessel ID)
def find_missing_edge_in_reference(missing_edge, reference_network_vessel_nodes):
    
    # Sort the nodes in the missing edge (since it's undirected)
    sorted_missing_edge = tuple(sorted(missing_edge))
    
    # Iterate over the reference network vessel nodes and find the matching edge
    for index, edge in enumerate(reference_network_vessel_nodes):
        
        # Sort the nodes in the current reference edge (to ensure undirected comparison)
        sorted_reference_edge = tuple(sorted(edge))
        
        # Check if the sorted missing edge matches the sorted reference edge
        if sorted_missing_edge == sorted_reference_edge:
            return index  # Return the index of the matching edge
    
    # If no match is found, raise an error or return None
    raise ValueError(f"Missing edge {missing_edge} not found in the reference network.")


# Define a function to return the details of the vessel that was pruned in a certain simulation
def log_killed_vessel(network_path, previous_network_path, previous_edge_metrics_df, reference_node_coordinates, reference_network_vessel_nodes):
    
    # Get the data for the current network
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(network_path)

    # Round the coordinates to two decimal places
    point_coordinates_data_array = point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    # node_coordinates = np.around(point_coordinates_data_array[:,0:2], 2)

    # Express the vessels in terms of node coordinates
    vessel_nodes, vessel_coordinates = get_vessel_nodes(polydata, reference_node_coordinates)
    
    # Repeat for the preceding network
    prev_point_coordinates_data_array, _, prev_cell_data_array, prev_polydata = get_vtk_data(previous_network_path)
    prev_point_coordinates_data_array = prev_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    # prev_node_coordinates = np.around(prev_point_coordinates_data_array[:,0:2], 2)
    prev_vessel_nodes, prev_vessel_coordinates = get_vessel_nodes(prev_polydata, reference_node_coordinates)
            
    # Convert vessel_nodes to sets of sorted tuples (to handle undirected edges)
    prev_edges_set = set(map(tuple, map(sorted, prev_vessel_nodes)))
    current_edges_set = set(map(tuple, map(sorted, vessel_nodes)))
    
    # Find the missing edge(s)
    missing_edges = list(prev_edges_set - current_edges_set)
    
    # Check if no missing edges were found
    if len(missing_edges) == 0:
        raise ValueError(f"No missing edge detected between the previous network and {network_path}. This is unexpected.")
    
    # Check if more than one missing edge was found
    if len(missing_edges) > 1:
        raise ValueError(f"More than one missing edge detected between the previous network and {network_path}. Multiple edges found: {missing_edges}. This is unexpected.")
    
    # Process the first (and only) missing edge
    missing_edge = missing_edges[0]
    
    # Use the find_missing_edge_in_reference function to find the index of the missing edge in prev_vessel_nodes
    index_in_prev = find_missing_edge_in_reference(missing_edge, prev_vessel_nodes)
    
    # Store the data for the missing edge
    killed_vessel_data = {key: value[index_in_prev] for key, value in prev_cell_data_array.items()}
    df_killed_vessel = pd.DataFrame([killed_vessel_data])

    # # Extract the node coordinates of the deleted vessel
    # x1, y1, x2, y2 = prev_vessel_coordinates[index_in_prev]
    # print(x1, y1, x2, y2)
    # # Add the coordinates as new columns
    # # df_killed_vessel['x1'] = x1#
    # # df_killed_vessel['y1'] = y1#
    # # df_killed_vessel['x2'] = x2#
    # # df_killed_vessel['y2'] = y2#
    
    # # Convert missing_edge coordinates (x1, y1, x2, y2) to node tuples
    # node_1 = np.where((prev_node_coordinates == (x1, y1)).all(axis=1))[0][0]
    # node_2 = np.where((prev_node_coordinates == (x2, y2)).all(axis=1))[0][0]
    
    # Create the node tuple representing the edge
    # missing_edge_as_nodes = tuple(sorted([node_1, node_2]))  # Sort to handle undirected edges
    index_in_reference = find_missing_edge_in_reference(missing_edge, reference_network_vessel_nodes)
    # print(index_in_reference)
    # print(previous_edge_metrics_df['edge_original_ID'])

    # Check the flow rate for the missing edge
    flow_rate = prev_cell_data_array['Absolute Vessel Flow Rate m^3/s'][index_in_prev]
    
    # Check if edge metrics are available for the missing edge
    edge_metrics = previous_edge_metrics_df[previous_edge_metrics_df['edge_original_ID'] == index_in_reference].to_dict('records')
    
    # If flow rate > 0, raise an error if no metrics are found
    if flow_rate > 0:
        if not edge_metrics:
            raise ValueError(f"No metrics found for the missing edge {index_in_reference} in the previous network.")
        
        # If metrics are found, assign them to df_killed_vessel
        edge_metrics = edge_metrics[0]  # Assuming one row is returned
        for metric, value in edge_metrics.items():
            if metric != 'edge':  # Skip the 'edge' column itself
                df_killed_vessel[metric] = value
    else:
        # If flow rate is 0, assign NaN to all metrics
        metrics_to_assign = previous_edge_metrics_df.columns.tolist()
        for metric in metrics_to_assign:
            if metric != 'edge':  # Skip the 'edge' column itself
                df_killed_vessel[metric] = np.nan  # Assign NaN
    
    # Return the updated dataframe for the killed vessel
    return df_killed_vessel

# Define a function to extract the predictive metrics from the Voronoi network
def get_voronoi_predictors(vtk_path, reference_node_coordinates, reference_network_coordinates, reference_network_vessel_nodes, pq_threshold, field_data_um, inlet_h):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_m = cell_data_array['Vessel Radius m']*2  # convert radii (m) to diameters (um)
    segment_diameters_um = segment_diameters_m*1000000  # convert radii (m) to diameters (um)
    segment_lengths = cell_data_array['Vessel Length']
    segment_lengths_um = segment_lengths*1000000  # convert length (m to um)

    # Get the architectural metrics    
    n_vessels = len(segment_diameters_um)
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    # mean_diameter = sum(non_inlet_segment_diameters_um)/n_vessels
    # mean_geometric_resistance = sum(non_inlet_segment_lengths_um/(non_inlet_segment_diameters_um**4))/n_vessels
    mean_diameter = sum(segment_diameters_um)/n_vessels  
    # mean_diameter = np.mean(segment_lengths_um/segment_diameters_um)  # quick and dirty way to check lambda for network
    mean_geometric_resistance = sum(segment_lengths_um/(segment_diameters_um**4))/n_vessels  # diameter_binned = filter_by_thresholds(segment_diameters_um, diameter_threshold_list)
    # print(vtk_path, segment_lengths_um, segment_diameters_um**4, n_vessels)
    # flow_rate_binned = filter_by_thresholds(cell_data_array['Absolute Vessel Flow Rate m^3/s'], flow_rate_threshold_list)
    segment_nodes, segment_xyxy = get_vessel_nodes(polydata, reference_node_coordinates)

    # Get the haematocrit values for vessels with flow
    h_array = cell_data_array['Vessel Haematocrit'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] > 0]

    # Haematocrit heterogeneity metrics
    h_std_dev = np.std(h_array)
    h_shannon = shannon_entropy(h_array)
    # h_data_kurtosis = kurtosis(h_array)
    # h_gini = gini_coefficient_optimized(h_array)
    h_inverse_simpson = inverse_simpson_index(h_array)
    h_hoover = hoover_index(h_array)
    h_haemo_potential = haemonormalisation_potential(h_array, h_target=inlet_h)
    h_theil = theil_index(h_array)
    h_atkinson = atkinson_index(h_array)
    h_pietra = pietra_index(h_array)
    h_mld = mean_log_deviation(h_array)
    h_cv = np.std(h_array) / np.mean(h_array) * 100

    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    # diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        node_1, node_2 = segment
        # print(node_1,node_2)
        
        # Fill in the adjacency matrices (don't think the diameters themselves are anything beyond identifiers?)
        # diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        # diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[node_1, node_2] = segment_lengths_um[segment_index]
        length_adjacency_matrix[node_2, node_1] = segment_lengths_um[segment_index]
        # print(node_1,node_2,segment_lengths_um[segment_index])
    
    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(length_adjacency_matrix)

    # Get the literature metrics
    total_blood_volume, fcd_numerator, vd_intersections, vd_uniques = get_perfusion_metrics_voronoi(cell_data_array, pq_threshold, segment_lengths_um, reference_network_coordinates, segment_xyxy, field_data_um) 
   
    # Get the network metrics
    network_metrics, G_filtered = get_network_metrics(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes)

    # Create a dictionary to store and return all metrics
    metrics = {
        'n_vessels': n_vessels,
        'n_perfused_vessels': n_perfused_vessels,
        'n_unperfused_vessels': n_unperfused_vessels,
        'mean_diameter': mean_diameter,
        'mean_geometric_resistance': mean_geometric_resistance,
        'n_cycles': n_cycles,
        'total_blood_volume': total_blood_volume,
        'fcd_numerator': fcd_numerator,
        'vd_intersections': vd_intersections,
        'vd_uniques': vd_uniques,
        'n_connected_components': n_connected_components,
        'network_metrics': network_metrics,  # Nested dictionary containing network metrics
        'tda_metrics': {
            'n_cycles_largest_connected_component': n_cycles_largest_connected_component,
            'size_largest_connected_component': size_largest_connected_component
        },
        'h_metrics': {
            'h_std_dev': h_std_dev,
            'h_cv': h_cv,
            # 'h_gini': h_gini,
            'h_shannon': h_shannon ,
            'h_inverse_simpson': h_inverse_simpson,
            'h_hoover': h_hoover,
            'h_haemo_potential': h_haemo_potential,
            'h_theil': h_theil,
            'h_atkinson': h_atkinson,
            'h_pietra': h_pietra,
            'h_mld': h_mld
        }
    }
    # print(segment_haem)
    return metrics, G_filtered

# Define a function to extract the predictive metrics from the biological network (no topological metrics to save time)
def get_bio_predictors(vtk_path, reference_node_coordinates, reference_network_coordinates, reference_network_vessel_nodes, pq_threshold, field_data_um, inlet_h):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_m = cell_data_array['Vessel Radius m']*2  # convert radii (m) to diameters (um)
    segment_diameters_um = segment_diameters_m*1000000  # convert radii (m) to diameters (um)
    segment_lengths = cell_data_array['Vessel Length']
    segment_lengths_um = segment_lengths*1000000  # convert length (m to um)

    # Get the architectural metrics    
    n_vessels = len(segment_diameters_um)
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    # mean_diameter = sum(non_inlet_segment_diameters_um)/n_vessels
    # mean_geometric_resistance = sum(non_inlet_segment_lengths_um/(non_inlet_segment_diameters_um**4))/n_vessels
    mean_diameter = sum(segment_diameters_um)/n_vessels  
    # mean_diameter = np.mean(segment_lengths_um/segment_diameters_um)  # quick and dirty way to check lambda for network
    mean_geometric_resistance = sum(segment_lengths_um/(segment_diameters_um**4))/n_vessels  # diameter_binned = filter_by_thresholds(segment_diameters_um, diameter_threshold_list)
    # print(vtk_path, segment_lengths_um, segment_diameters_um**4, n_vessels)
    # flow_rate_binned = filter_by_thresholds(cell_data_array['Absolute Vessel Flow Rate m^3/s'], flow_rate_threshold_list)
    segment_nodes, segment_xyxy = get_vessel_nodes(polydata, reference_node_coordinates)

    # Get the haematocrit values for vessels with flow
    h_array = cell_data_array['Vessel Haematocrit'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] > 0]

    # Haematocrit heterogeneity metrics
    # h_std_dev = np.std(h_array)
    # h_shannon = shannon_entropy(h_array)
    # h_data_kurtosis = kurtosis(h_array)
    # h_gini = gini_coefficient_optimized(h_array)
    h_inverse_simpson = inverse_simpson_index(h_array)
    # h_hoover = hoover_index(h_array)
    # h_haemo_potential = haemonormalisation_potential(h_array, h_target=inlet_h)
    # h_theil = theil_index(h_array)
    # h_atkinson = atkinson_index(h_array)
    # h_pietra = pietra_index(h_array)
    # h_mld = mean_log_deviation(h_array)
    # h_cv = np.std(h_array) / np.mean(h_array) * 100

    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    # diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        node_1, node_2 = segment
        # print(node_1,node_2)
        
        # Fill in the adjacency matrices (don't think the diameters themselves are anything beyond identifiers?)
        # diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        # diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[node_1, node_2] = segment_lengths_um[segment_index]
        length_adjacency_matrix[node_2, node_1] = segment_lengths_um[segment_index]
        # print(node_1,node_2,segment_lengths_um[segment_index])
    
    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(length_adjacency_matrix)

    # Get the literature metrics
    total_blood_volume, fcd_numerator, vd_intersections, vd_uniques = get_perfusion_metrics_voronoi(cell_data_array, pq_threshold, segment_lengths_um, reference_network_coordinates, segment_xyxy, field_data_um) 
   
    # Get the network metrics
    # network_metrics, G_filtered = get_network_metrics(length_adjacency_matrix, segment_nodes, cell_data_array, point_data_array, reference_network_vessel_nodes)

    # Create a dictionary to store and return all metrics
    metrics = {
        'n_vessels': n_vessels,
        'n_perfused_vessels': n_perfused_vessels,
        'n_unperfused_vessels': n_unperfused_vessels,
        'mean_diameter': mean_diameter,
        'mean_geometric_resistance': mean_geometric_resistance,
        'n_cycles': n_cycles,
        'total_blood_volume': total_blood_volume,
        'fcd_numerator': fcd_numerator,
        'vd_intersections': vd_intersections,
        'vd_uniques': vd_uniques,
        'n_connected_components': n_connected_components,
        # 'network_metrics': network_metrics,  # Nested dictionary containing network metrics
        'tda_metrics': {
            'n_cycles_largest_connected_component': n_cycles_largest_connected_component,
            'size_largest_connected_component': size_largest_connected_component
        },
        'h_metrics': {
            # 'h_std_dev': h_std_dev,
            # 'h_cv': h_cv,
            # 'h_gini': h_gini,
            # 'h_shannon': h_shannon ,
            'h_inverse_simpson': h_inverse_simpson,
            # 'h_hoover': h_hoover,
            # 'h_haemo_potential': h_haemo_potential,
            # 'h_theil': h_theil,
            # 'h_atkinson': h_atkinson,
            # 'h_pietra': h_pietra,
            # 'h_mld': h_mld
        }
    }
    # print(segment_haem)
    return metrics

# Import the generation coordinates, rank lengths, and node coordinates for the forking network 
def get_reference_forking_network(generations, inlet_diameter):
    alpha_value = '1.00'  # only the diameter varies with alpha 
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_inlet_diameter_' + str(inlet_diameter) + '_um/Alpha' + alpha_value + '/FinalHaematocrit.vtp'
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_inlet_diameter_50_um_five_generations/FinalHaematocrit.vtp'
    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_networks/inlet_diameter_' + str(inlet_diameter) + '_um_generations_' + str(generations) + '/Alpha' + alpha_value + '/FinalHaematocrit.vtp'
    reference_point_coordinates_data_array, point_data_array, reference_cell_data_array, polydata = get_vtk_data(reference_network_path)
    generation_coordinates = np.unique(reference_point_coordinates_data_array[:,0])
    rank_diameters = -np.sort(-np.unique(reference_cell_data_array['Vessel Radius m'])*2)
    reference_rank_lengths = rank_diameters*4
    
    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)
  
    # Express the vessels in terms of node coordinates
    reference_network_vessel_nodes, reference_network_vessel_coordinates = get_vessel_nodes(polydata, reference_node_coordinates)
 
    return generation_coordinates, reference_rank_lengths, reference_node_coordinates, reference_network_vessel_nodes, reference_network_vessel_coordinates

# Import the node coordinates for the hexagonal network 
def get_reference_hexagonal_network(reference_network_path):
#def get_reference_hexagonal_network(vessel_length_m=100*(10**-6)):
#    sigma = 1
#    selection = 1
#    radius_threshold = 0
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/Sigma' + str(sigma) + '/Selection' + str(selection) +'/RadiusThreshold' + str(radius_threshold) + '/FinalHaematocrit.vtp'
    
    # reference_point_coordinates_data_array, _, _, _ = get_vtk_data(reference_network_path)
    reference_point_coordinates_data_array, _, _, polydata = get_vtk_data(reference_network_path)

    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)

    # Express the vessels in terms of node coordinates
    reference_network_vessel_nodes, reference_network_vessel_coordinates = get_vessel_nodes(polydata, reference_node_coordinates)

    # return reference_node_coordinates
    return reference_node_coordinates, reference_network_vessel_nodes, reference_network_vessel_coordinates

# Define a function to plot the diameter distribution of a network and return the diameter data
def get_diameter_distribution(reference_network_path):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(reference_network_path)
    
    # Convert the radii in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    mean_diameter = np.mean(segment_diameters_um)
    std_dev = np.std(segment_diameters_um)
    
    # Set the figure layout
    fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)
    fig.subplots_adjust(hspace = 0.75, wspace=.25)

    # Plot the distribution
    n, bins, patches = axs.hist(x=segment_diameters_um, bins='auto', alpha=0.7)
    axs.set_ylabel('number of vessels') 
    axs.set_xlabel('vessel diameter (μm)')    
    axs.title.set_text('${σ}$ = ' + str(round(std_dev,2)) + ' μm')
    axs.tick_params(labelbottom=True)
    axs.axvline(mean_diameter, c='black', ls='--', label='${µ}$ = ' + str(round(mean_diameter,2)) + ' μm')
    axs.set_ylim(0, 75)
    axs.legend()
    
    # Return the diameter stats
    return segment_diameters_um, mean_diameter, std_dev

# Define a function to return the oxygen field
def get_oxygen_field(field_path):

    field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
    flat_field = field_data['oxygen'].to_numpy().flatten()
    
    return flat_field

# Define a function to get the normoxic fraction from the oxygen field of a simulation
def get_field_NF(folder_path, threshold):
    
    flat_field = get_oxygen_field(folder_path+'/oxygen_solution_0.vti')
    number_of_points = 1
    for dim in np.shape(flat_field): number_of_points *= dim
    normoxic_points = (flat_field > threshold).sum()
    field_NF = normoxic_points/number_of_points  # calculate the hypoxic fraction

    return field_NF, flat_field    

# Define a function to return the field and cellular NFs of a simulation
def compare_field_cells_NF(folder_path, threshold):
    
    # Compute the cell oxygen distribution at t=0
    pvd_df = get_pvd_data(folder_path + '/results.pvd')
    pvd_df.rename(columns={'Block_ID': 'time'}, inplace=True)
    pre_RT = pvd_df[pvd_df['time'] == 0]
    cells_NF = compute_cells_NF(pre_RT['oxygen'],threshold)
    
    # Compute the O2 distribution
    field_NF, _ = get_field_NF(folder_path, threshold)

    return cells_NF, field_NF

# Compute the normoxic fraction given the population of cell oxygen concentrations
def compute_cells_NF(c, radioresistant_threshold_nM):
    
    # Count the normoxic fraction of cells
    norm_cells = np.sum(c > radioresistant_threshold_nM)
    norm_fraction = norm_cells/len(c)
    return norm_fraction

# =============================================================================
# ARCHIVE
# =============================================================================

# Define a function to bin an array according to three thresholds
def filter_by_thresholds(array, threshold_list):

    # Extract the thresholds
    threshold_1 = threshold_list[0]
    threshold_2 = threshold_list[1]
    threshold_3 = threshold_list[2]
    
    # Bin the array
    subset_1 = sum(x < threshold_1 for x in array)
    subset_2 = sum(threshold_1 <= x < threshold_2 for x in array)
    subset_3 = sum(threshold_2 <= x < threshold_3 for x in array)
    subset_4 = sum(threshold_3 <= x for x in array)

    # Return the bins
    return subset_1, subset_2, subset_3, subset_4
        
# Define a function to filter the simulation stats based on alpha values and return individual matrices (used in Paper 1)
def filter_by_alpha_old(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,3]
        min_data = alpha_array[:,4]
        half_data = alpha_array[:,5]
        max_data = alpha_array[:,6]
        sd_data = alpha_array[:,7]
        n_vessels_data = alpha_array[:,8]
        n_perfused_vessels_data = alpha_array[:,9]
        n_unperfused_vessels_data = alpha_array[:,10]
        mean_diameter_data = alpha_array[:,11]
        mean_geometric_resistance_data = alpha_array[:,12]  
        diameter_binned_data = alpha_array[:,13]  
        flow_rate_binned_data = alpha_array[:,14]  
        n_connected_components_data = alpha_array[:,15]  
        n_cycles_data = alpha_array[:,16]  
        size_largest_connected_component_data = alpha_array[:,17]  
        n_cycles_largest_connected_component_data = alpha_array[:,18]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 

# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_alpha(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    n_cycles_composite = np.array([])
    total_blood_flow_composite = np.array([])
    fcd_composite = np.array([])
    vd_intersections_composite = np.array([])
    vd_uniques_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,3]
        min_data = alpha_array[:,4]
        half_data = alpha_array[:,5]
        max_data = alpha_array[:,6]
        sd_data = alpha_array[:,7]
        n_vessels_data = alpha_array[:,8]
        n_perfused_vessels_data = alpha_array[:,9]
        n_unperfused_vessels_data = alpha_array[:,10]
        mean_diameter_data = alpha_array[:,11]
        mean_geometric_resistance_data = alpha_array[:,12]  
        n_cycles_data = alpha_array[:,13]  
        total_blood_flow_data = alpha_array[:,14]  
        fcd_data = alpha_array[:,15]  
        vd_intersections_data = alpha_array[:,16]  
        vd_uniques_data = alpha_array[:,17]  
        
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        total_blood_flow_composite = np.vstack([total_blood_flow_composite, total_blood_flow_data]) if total_blood_flow_composite.size else total_blood_flow_data
        fcd_composite = np.vstack([fcd_composite, fcd_data]) if fcd_composite.size else fcd_data
        vd_intersections_composite = np.vstack([vd_intersections_composite, vd_intersections_data]) if vd_intersections_composite.size else vd_intersections_data
        vd_uniques_composite = np.vstack([vd_uniques_composite, vd_uniques_data]) if vd_uniques_composite.size else vd_uniques_data

    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, n_cycles_composite, total_blood_flow_composite, fcd_composite, vd_intersections_composite, vd_uniques_composite

# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_alpha_fork_stoch(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    anoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        anoxic_fraction_data = alpha_array[:,3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,4]
        min_data = alpha_array[:,5]
        half_data = alpha_array[:,6]
        max_data = alpha_array[:,7]
        sd_data = alpha_array[:,8]
        n_vessels_data = alpha_array[:,9]
        n_perfused_vessels_data = alpha_array[:,10]
        n_unperfused_vessels_data = alpha_array[:,11]
        mean_diameter_data = alpha_array[:,12]
        mean_geometric_resistance_data = alpha_array[:,13]  
        diameter_binned_data = alpha_array[:,14]  
        flow_rate_binned_data = alpha_array[:,15]  
        n_connected_components_data = alpha_array[:,16]  
        n_cycles_data = alpha_array[:,17]  
        size_largest_connected_component_data = alpha_array[:,18]  
        n_cycles_largest_connected_component_data = alpha_array[:,19]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        anoxic_fraction_composite = np.vstack([anoxic_fraction_composite, anoxic_fraction_data]) if anoxic_fraction_composite.size else anoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, anoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 

# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_mean_hex(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,3]
        min_data = alpha_array[:,4]
        half_data = alpha_array[:,5]
        max_data = alpha_array[:,6]
        sd_data = alpha_array[:,7]
        n_vessels_data = alpha_array[:,8]
        n_perfused_vessels_data = alpha_array[:,9]
        n_unperfused_vessels_data = alpha_array[:,10]
        mean_diameter_data = alpha_array[:,11]
        mean_geometric_resistance_data = alpha_array[:,12]  
        diameter_binned_data = alpha_array[:,13:17]  
        flow_rate_binned_data = alpha_array[:,17:21]  
        n_connected_components_data = alpha_array[:,21]  
        n_cycles_data = alpha_array[:,22]  
        size_largest_connected_component_data = alpha_array[:,23]  
        n_cycles_largest_connected_component_data = alpha_array[:,24]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 
