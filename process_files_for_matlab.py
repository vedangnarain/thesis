#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:22:29 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This script extracts data from a VTK network and converts it to the format used by the MATLAB solver.
"""

# Initialise libraries
import pathlib
import math
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.io import savemat
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Import custom functions
from get_paraview_data import *

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTeX-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})





# Enter details to allow looping over folders
seed_selection = 100
segment_length_um = 100
max_layout = 101
sd_list = ['8.68', '13.23', '17.49']
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
mean_list = ['22.76', '28.5', '33.64']

# Specify the simulations to examine
solver_name = solver_list[0]
selection = 2
sigma = sd_list[0]
mu = mean_list[1]
kills = 100
# main_folder_path = '/scratch/narain/Hexagonal/With Oxygen/' + str(segment_length_um) + ' VesselLength/TestHexagonalNetwork/'
# vtk_path = main_folder_path + solver_name + 'Haematocrit/Selection' + str(selection) + '/Sigma' + sigma + '/Mu' + mu + '/Kills' + str(kills) + '/FinalHaematocrit.vtp'
# main_folder_path = '/tmp/narain/testoutput/TestDebuggingNetworks/HexagonalNeighbourhood/HeterogeneousDiameters/NoPruning/PriesHaematocrit/Selection1/Sigma8.68/Mu22.76/'
main_folder_path = '/tmp/narain/testoutput/TestDebuggingNetworks/HexagonalNeighbourhood/HomogeneousDiameters/NoPruning/PriesHaematocrit/'
vtk_path = main_folder_path + 'FinalHaematocrit.vtp'
# Get the .vtk data
point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
# reference_node_coordinates, reference_network_vessel_nodes, reference_network_coordinates = get_reference_hexagonal_network(main_folder_path + solver_name + 'Haematocrit/Selection' + str(selection) + '/Sigma' + sigma + '/Mu' + mu + '/Kills0/FinalHaematocrit.vtp')  # Just use the unpruned network as a reference network
# segment_nodes, segment_xyxy = get_vessel_nodes(polydata, reference_node_coordinates)
# segment_nodes, segment_xyxy = get_vessel_nodes(polydata, vtk_path)

# =============================================================================
# seg_nodes
# =============================================================================

# Get a list of vessel segments with node IDs
reference_node_coordinates, reference_network_vessel_nodes, reference_network_coordinates = get_reference_hexagonal_network(vtk_path)  # Just use the unpruned network as a reference network
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
seg_nodes = segment_nodes.T

# =============================================================================
# connections_node
# =============================================================================

# Create a dictionary to store connections for each node
connections_dict = {}

# Iterate through segment_nodes to populate connections_dict
for start, end in segment_nodes:
    if start not in connections_dict:
        connections_dict[start] = []
    if end not in connections_dict:
        connections_dict[end] = []
    
    connections_dict[start].append(end)
    connections_dict[end].append(start)

# Find the maximum node ID to determine the size of the connections_node array
max_node_id = max(connections_dict.keys())

# Create the connections_node array with -1 as default values
connections_node = np.full((max_node_id + 1, 3), -1)

# Populate the connections_node array with the connections
for node, connections in connections_dict.items():
    connections_node[node, :len(connections)] = connections

# Remove rows with no connections
connections_node = connections_node[np.any(connections_node != -1, axis=1)].T

# =============================================================================
# connections_seg
# =============================================================================

# Create a dictionary to store segment IDs for each node
segment_dict = {}

# Iterate through segment_nodes to populate segment_dict
for idx, (start, end) in enumerate(segment_nodes):
    if start not in segment_dict:
        segment_dict[start] = []
    if end not in segment_dict:
        segment_dict[end] = []
    
    segment_dict[start].append(idx)
    segment_dict[end].append(idx)

# Find the maximum node ID to determine the size of the connections_seg array
max_node_id = max(segment_dict.keys())

# Create the connections_seg array with -1 as default values
connections_seg = np.full((max_node_id + 1, 3), -1)

# Populate the connections_seg array with the segment IDs
for node, segments in segment_dict.items():
    connections_seg[node, :len(segments)] = segments

# Remove rows with no segments
connections_seg = connections_seg[np.any(connections_seg != -1, axis=1)].T

# =============================================================================
# seg_diameter
# =============================================================================

# Convert the radii and lengths in metres into the right units
segment_diameters_m = cell_data_array['Vessel Radius m']*2  # convert radii (m) to diameters (um)
segment_diameters_um = segment_diameters_m*1000000  # convert radii (m) to diameters (um)
# segment_viscosities = cell_data_array['Vessel Viscosity Pa.s'] 
# segment_impedances = cell_data_array['Vessel Impedance kg/m^4/s']
# segment_lengths_um = ((segment_impedances*math.pi*((segment_diameters_m/2)**4))/(8.0*segment_viscosities))*1000000  # convert length (m to um)
seg_diameter = np.array([segment_diameters_um])

# =============================================================================
# seg_length
# =============================================================================

segment_lengths_um = [segment_length_um] * len(segment_diameters_um)
seg_length = np.array([segment_lengths_um])
seg_length = seg_length.astype(np.float64)

# =============================================================================
# non_boundary_nodes and boundary_nodes
# =============================================================================

data_dict = point_data_array

# Convert the lists to numpy arrays for easier indexing
node_is_input = np.array(data_dict['Node Is Input'])
node_is_output = np.array(data_dict['Node Is Output'])
node_id = np.array(data_dict['Node Id'])

# Find indices where 'Node Is Input' or 'Node Is Output' is equal to 1
inlets_to_remove = np.where((node_is_input == 1))[0]
outlets_to_remove = np.where((node_is_output == 1))[0]
indices_to_remove = np.concatenate((inlets_to_remove, outlets_to_remove), axis=0)

# Extract the 'Node ID' values at those indices
inlet_node_ids = np.array([node_id[inlets_to_remove]])
outlet_node_ids = np.array([node_id[outlets_to_remove]])
boundary_node_ids = np.array([node_id[indices_to_remove]])

# Remove the 'Node ID' values at those indices from the original 'Node ID' array
remaining_node_ids = np.delete(node_id, indices_to_remove)
non_boundary_nodes = np.array([remaining_node_ids])

# =============================================================================
# WRITE TO FILES
# =============================================================================

# Dictionary to hold arrays
network_architecture = {
    'segment_diameters': seg_diameter,
    'segment_lengths': seg_length,
    'segment_nodes': seg_nodes,
    'connections_node': connections_node,
    'connections_segments': connections_seg,
    'non_boundary_node_IDs': non_boundary_nodes,
    'boundary_node_IDs': boundary_node_ids,
    'inlet_node_ids': inlet_node_ids,
    'outlet_node_ids': outlet_node_ids
    }

# Remove zero indexing for IDs
keys_to_update = ['segment_nodes', 'connections_node', 'connections_segments', 
                  'non_boundary_node_IDs', 'boundary_node_IDs', 'inlet_node_ids', 'outlet_node_ids']

# Add 1 to each element in the specified keys
for key in keys_to_update:
    if key in network_architecture:
        network_architecture[key] = network_architecture[key] + 1


mat_name = main_folder_path + 'network_architecture.mat'
savemat(mat_name, network_architecture)

'''
ids = cell_data_array['Vessel Id']
print(min(diameters))
# print(np.min(cell_data_array['Absolute Vessel Flow Rate m^3/s'][np.nonzero(cell_data_array['Absolute Vessel Flow Rate m^3/s'])]))
print(np.max(cell_data_array['Absolute Vessel Flow Rate m^3/s']))
# print(min(cell_data_array['Absolute Vessel Flow Rate m^3/s']))
print(ids[np.where(diameters == min(diameters))[0]])
'''