#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 18:32:08 2024

@author: narain

A script to copy all the verified edge matrices and diameter distributions into a verified folder.

Use one block at a time. Exercise caution when using operators: make sure you don't keep copying and renaming the same files. Read the comments before using an operator.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise the key libraries
import os
import pandas as pd 
import shutil
import re

# =============================================================================
# SIMULATION DETAILS
# =============================================================================

# Specify the diameter distribution values
# sd_list = ['8.68', '13.23', '17.49']
# mean_list = ['22.76', '28.5', '33.64']

# Toggles 
seeds = 60
# sd = sd_list[2]
# mean = mean_list[2]
# first_selection_number = 102
# last_selection_number = 351

# Folder paths
# main_folder_path = '/scratch/narain/Voronoi/60 Seeds/H 0.3/no_rt_size_pruning/TestVoronoiNetworkUpto500 (535 in original edge names)/'
# main_folder_path = '/scratch/narain/Voronoi/60 Seeds/H 0.3/no_rt_size_pruning/TestVoronoiNetwork/'
main_folder_path = '/scratch/narain/testoutput/TestVoronoiNetwork/'
# main_folder_path = '/scratch/narain/Voronoi/60 Seeds/H 0.3/no_rt_flow_pruning/TestVoronoiNetwork/'
broken_selections_filename = main_folder_path + 'broken_layouts.txt'
simulation_data_path = main_folder_path + 'MemoryHaematocrit'
edges_base_folder = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/voronoi/60SeedPoints/'
edges_source_folder =  edges_base_folder + '/current'
edges_destination_folder = edges_base_folder + '/verified'

# List of good simulations
# good_selection_numbers = [4, 5, 18, 19, 26, 35, 37, 42, 50, 53, 57, 62, 67, 70, 71, 78, 79, 89, 90, 92, 94]  # for up to 100 seeds

# =============================================================================
# READ A LIST OF BROKEN LAYOUTS
# =============================================================================
'''
# Path to the text file
file_path = broken_selections_filename

# List to store the extracted numbers
bad_selection_numbers = []

# Read the file and extract numbers
with open(file_path, 'r') as file:
    for line in file:
        # Use regular expression to find numbers after "Selection" and before "/"
        match = re.search(r'Selection(\d+)/', line)
        if match:
            # Append the found number to the list
            bad_selection_numbers.append(int(match.group(1)))
'''
# =============================================================================
# STEP 1: DELETE THE BAD LAYOUTS
# =============================================================================
'''
# Delete bad selections
# bad_selection_numbers = [3,261,395]

# Iterate over the list and delete matching folders
for number in bad_selection_numbers:
    folder_name = f'Selection{number}'
    folder_path = os.path.join(simulation_data_path, folder_name)
    
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)
        print(f'Deleted folder: {folder_name}')
    else:
        print(f'Folder not found: {folder_name}')
'''
# =============================================================================
# STEP 1.5: RENAME THE LAYOUTS
# =============================================================================
'''
starting_number = 959  # the number that you want the renaming to begin with

# Rename all the selection folders
# Get all folders in the directory
folders = os.listdir(simulation_data_path)

# Filter for folders that match the pattern 'Selection<number>'
selection_folders = [folder for folder in folders if folder.startswith('Selection') and folder[9:].isdigit()]

# Sort folders by the numeric value after 'Selection'
selection_folders.sort(key=lambda x: int(x[9:]))

# Rename folders in ascending order
for i, folder_name in enumerate(selection_folders, start=starting_number):
    old_path = os.path.join(simulation_data_path, folder_name)
    new_folder_name = f'Selection{i}'
    new_path = os.path.join(simulation_data_path, new_folder_name)
    
    if old_path != new_path:
        os.rename(old_path, new_path)
        print(f'Renamed folder: {folder_name} -> {new_folder_name}')
'''
# =============================================================================
# STEP 2: LIST THE GOOD SELECTIONS AND COPY THEIR EDGE MATRICES
# =============================================================================
'''
first_good_simulation_number = 965
last_good_simulation_number = 1250

# Get the good selection numbers
all_numbers = list(range(first_good_simulation_number, last_good_simulation_number + 1))
good_selection_numbers = [number for number in all_numbers if number not in bad_selection_numbers]
# good_selection_numbers = all_numbers

# Copy the good edge matrices into a new verified folder
for selection_number in good_selection_numbers:
    file_pattern = f'EdgesMatrixSampleNumber{selection_number}.txt'
    
    # Check for exact matches in the source folder
    matching_files = [file for file in os.listdir(edges_source_folder) if file == file_pattern]
    for file_name in matching_files:
        source_path = os.path.join(edges_source_folder, file_name)
        destination_path = os.path.join(edges_destination_folder, file_name)

        # Copy the file to the destination folder
        shutil.copy2(source_path, destination_path)
        print(f"Copied: {source_path} to {destination_path}")
'''
# =============================================================================
# STEP 3: RENAME THE GOOD EDGE MATRICES IN THE NEW LOCATION
# =============================================================================
'''
# Rename the good edge matrices
edge_files = [f for f in os.listdir(edges_destination_folder) if os.path.isfile(os.path.join(edges_destination_folder, f)) and f.startswith('EdgesMatrixSampleNumber')]
edge_files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
for index, file_name in enumerate(edge_files, start=1):
    current_path = os.path.join(edges_destination_folder, file_name)
    new_name = f"EdgesMatrixSampleNumber{index}{os.path.splitext(file_name)[1]}"  # Keep the file extension
    new_path = os.path.join(edges_destination_folder, new_name)

    # Rename the file
    try:
        os.rename(current_path, new_path)
        print(f"Renamed: {current_path} to {new_path}")
    except OSError as e:
        print(f"Error renaming file {current_path}: {e}")
'''
# =============================================================================
# FUNCTIONS
# =============================================================================

# For a list of Selections, copy the corresponding Edge Matrices into a new 'verified' folder
def copy_edge_matrices(selection_numbers, edges_source_folder, edges_destination_folder):

    for selection_number in selection_numbers:
        file_pattern = f'EdgesMatrixSampleNumber{selection_number}.txt'
        
        # Check for exact matches in the source folder
        matching_files = [file for file in os.listdir(edges_source_folder) if file == file_pattern]
    
        for file_name in matching_files:
            source_path = os.path.join(edges_source_folder, file_name)
            destination_path = os.path.join(edges_destination_folder, file_name)
    
            # Copy the file to the destination folder
            shutil.copy2(source_path, destination_path)
            print(f"Copied: {source_path} to {destination_path}")

# For a list of Selections, copy the corresponding ID and radii lists into a new 'verified' folder
def copy_radii_and_ids(selection_numbers, radii_and_id_source_folder, radii_and_id_destination_folder):

    # Copy ID lists
    for selection_number in selection_numbers:
        file_pattern = f'id_list_{selection_number}.txt'
        
        # Check for exact matches in the source folder
        matching_files = [file for file in os.listdir(radii_and_id_source_folder) if file == file_pattern]
    
        for file_name in matching_files:
            source_path = os.path.join(radii_and_id_source_folder, file_name)
            destination_path = os.path.join(radii_and_id_destination_folder, file_name)
    
            # Copy the file to the destination folder
            shutil.copy2(source_path, destination_path)
            print(f"Copied: {source_path} to {destination_path}")

    # Copy radii lists
    for selection_number in selection_numbers:
        file_pattern = f'radii_list_{selection_number}.txt'
        
        # Check for exact matches in the source folder
        matching_files = [file for file in os.listdir(radii_and_id_source_folder) if file == file_pattern]
    
        for file_name in matching_files:
            source_path = os.path.join(radii_and_id_source_folder, file_name)
            destination_path = os.path.join(radii_and_id_destination_folder, file_name)
    
            # Copy the file to the destination folder
            shutil.copy2(source_path, destination_path)
            print(f"Copied: {source_path} to {destination_path}")

# Rename 'id_list_' and 'radii_list' files in ascending order
def rename_edge_matrices(edges_destination_folder):

    # Get a list of files in the folder for 'id_list_'
    edge_files = [f for f in os.listdir(edges_destination_folder) if os.path.isfile(os.path.join(edges_destination_folder, f)) and f.startswith('EdgesMatrixSampleNumber')]
    edge_files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
    
    
    for index, file_name in enumerate(edge_files, start=1):
        current_path = os.path.join(edges_destination_folder, file_name)
        new_name = f"EdgesMatrixSampleNumber{index}{os.path.splitext(file_name)[1]}"  # Keep the file extension
        new_path = os.path.join(edges_destination_folder, new_name)
    
        # Rename the file
        try:
            os.rename(current_path, new_path)
            print(f"Renamed: {current_path} to {new_path}")
        except OSError as e:
            print(f"Error renaming file {current_path}: {e}")

# Rename 'id_list_' and 'radii_list' files in ascending order
def rename_radii_and_id_lists(radii_and_id_destination_folder):

    # Get a list of files in the folder for 'id_list_'
    id_list_files = [f for f in os.listdir(radii_and_id_destination_folder) if os.path.isfile(os.path.join(radii_and_id_destination_folder, f)) and f.startswith('id_list_')]
    id_list_files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
    
    
    for index, file_name in enumerate(id_list_files, start=1):
        current_path = os.path.join(radii_and_id_destination_folder, file_name)
        new_name = f"id_list_{index}{os.path.splitext(file_name)[1]}"  # Keep the file extension
        new_path = os.path.join(radii_and_id_destination_folder, new_name)
    
        # Rename the file
        try:
            os.rename(current_path, new_path)
            print(f"Renamed: {current_path} to {new_path}")
        except OSError as e:
            print(f"Error renaming file {current_path}: {e}")

    # Get a list of files in the folder for 'radii_list'
    radii_list_files = [f for f in os.listdir(radii_and_id_destination_folder) if os.path.isfile(os.path.join(radii_and_id_destination_folder, f)) and f.startswith('radii_list')]
    radii_list_files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
        

    for index, file_name in enumerate(radii_list_files, start=1):
        current_path = os.path.join(radii_and_id_destination_folder, file_name)
        new_name = f"radii_list_{index}{os.path.splitext(file_name)[1]}"  # Keep the file extension
        new_path = os.path.join(radii_and_id_destination_folder, new_name)
    
        # Rename the file
        try:
            os.rename(current_path, new_path)
            print(f"Renamed: {current_path} to {new_path}")
        except OSError as e:
            print(f"Error renaming file {current_path}: {e}")

# Rename all the Selection folders in a directory in ascending order
def rename_folders_in_ascending_order(simulation_data_path):

    # Get a list of folder names
    folders = [folder for folder in os.listdir(simulation_data_path) if os.path.isdir(os.path.join(simulation_data_path, folder)) and folder.startswith("Selection")]
    
    # Sort the folders based on the numeric value at the end of each folder name
    folders.sort(key=lambda x: int(x.split("Selection")[-1]))
    
    # Rename the folders in ascending order
    for index, old_folder_name in enumerate(folders, start=1):
        old_folder_path = os.path.join(simulation_data_path, old_folder_name)
        new_folder_name = f"Selection{index}"
        new_folder_path = os.path.join(simulation_data_path, new_folder_name)
    
        # Rename the folder
        try:
            os.rename(old_folder_path, new_folder_path)
            print(f"Renamed: {old_folder_path} to {new_folder_path}")
        except OSError as e:
            print(f"Error renaming folder {old_folder_path}: {e}")
'''
# Copy and rename the radii and ID lists
def copy_and_rename_radii_and_ids(selection_numbers, copy=0):
    for sd in sd_list:
        for mean in mean_list:
            common = 'voronoi_diameter_log_normal_distribution_sigma_' + sd + '/' + 'mu_' + mean
            radii_and_id_source_folder = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/' + str(seeds) + 'SeedPoints/' + common
            radii_and_id_destination_folder = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/' + str(seeds) + 'SeedPoints/verified_diameter_distributions/' + common
            if copy==1:
                copy_radii_and_ids(selection_numbers, radii_and_id_source_folder, radii_and_id_destination_folder)
            rename_radii_and_id_lists(radii_and_id_destination_folder)
'''
# Copy and rename the edge matrices
def copy_and_rename_edge_matrices(selection_numbers, edges_source_folder, edges_destination_folder, copy=0):
    if copy==1:
        copy_edge_matrices(selection_numbers, edges_source_folder, edges_destination_folder)
    rename_edge_matrices(edges_destination_folder)
'''
# Read PF file and remove bad simulations (code assumes the PFs are in a combined file already)
def refine_pf_df(main_folder_path):

    pf_filename = main_folder_path + 'voronoi_lognormal_individual_perfusion_quotients.txt'
    # pf_filename = '/scratch/narain/Voronoi/TestVoronoiNetwork/voronoi_lognormal_individual_perfusion_quotients.txt'
    pf_df = pd.read_csv(pf_filename, index_col=None, delim_whitespace=True, names=["network_name", "solver_name", "selection", "sigma", "mu", "kills", "PPV"])
    filtered_pf_df = pf_df[~pf_df['selection'].isin(bad_selection_numbers['selection'].tolist())]
    # filtered_pf_df = pf_df[pf_df['selection'].isin(selection_numbers)]
    
    # Create a mapping dictionary to map original selection numbers to new numbers
    mapping = {}
    new_number = 1
    for number in sorted(filtered_pf_df['selection'].unique()):
        mapping[number] = new_number
        new_number += 1
    
    # Create a copy of the DataFrame to avoid the warning
    refined_pf_df = filtered_pf_df.copy()
    
    # Apply the mapping to the 'selection' column of the copy
    refined_pf_df['selection'] = refined_pf_df['selection'].map(mapping)

    return refined_pf_df
'''        
# =============================================================================
# OPERATIONS
# =============================================================================
            
# Copy and rename the radii
# copy_and_rename_radii_and_ids(selection_numbers, copy=0)  # only change to copy=1 if there are completely different simulations not already in the directory

# Copy and rename the edge matrices
# copy_and_rename_edge_matrices(selection_numbers, edges_source_folder, edges_destination_folder, copy=0)  # only change to copy=1 if there are completely different simulations not already in the directory

# Rename the simulation folders
# rename_folders_in_ascending_order(simulation_data_path)

# Filter the PFs (code assumes the PFs are in a combined file already)
# refined_pf_df = refine_pf_df(main_folder_path)

# Write the dataframe
# refined_pf_df.to_csv(main_folder_path + 'refined_voronoi_perfusion_quotients.txt', sep=' ', index=False, header=False)
