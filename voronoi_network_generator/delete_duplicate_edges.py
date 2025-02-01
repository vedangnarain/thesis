#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:03:19 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

This code is used to check a folder for duplicate files. Can be useful when 
dealing with lots of data, particularly for checking edge matrices.
"""

# Import necessary libraries
import os
import hashlib
import re

# =============================================================================
# FUNCTIONS
# =============================================================================

def hash_file(filepath):
    """Generate a hash for the file."""
    hasher = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()

def extract_number_from_filename(filename):
    """Extract the numeric part from a filename using regex."""
    match = re.search(r'\d+', filename)  # Look for numbers in the filename
    return int(match.group()) if match else float('inf')  # Return the number or inf if no number is found

def find_duplicates(folder):
    """Find and flag duplicate files in the specified folder."""
    file_hashes = {}
    duplicates = []

    for dirpath, _, filenames in os.walk(folder):
        for filename in filenames:
            if filename.startswith("EdgesMatrixSampleNumber"):  # Check for specific file prefix
                file_path = os.path.join(dirpath, filename)
                file_hash = hash_file(file_path)

                # Extract number from filename
                file_number = extract_number_from_filename(filename)

                if file_hash in file_hashes:
                    # Compare numbers, keep the lower one as original
                    existing_file, existing_number = file_hashes[file_hash]
                    if file_number < existing_number:
                        duplicates.append((existing_file, filename))  # Add previous original as duplicate
                        file_hashes[file_hash] = (filename, file_number)  # Update original to the lower number
                    else:
                        duplicates.append((filename, existing_file))  # Add new file as duplicate
                else:
                    file_hashes[file_hash] = (filename, file_number)

    return duplicates  # Return list of tuples with (duplicate, original)

def main():
    folder = "/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/edges/voronoi/60SeedPoints/verified"  # Your target folder

    if not os.path.isdir(folder):
        print("Invalid directory. Please check the folder path.")
        return

    duplicates = find_duplicates(folder)

    bad_selection_numbers = []  # Initialize the list for duplicate numbers

    if duplicates:
        print("Duplicate files found: ", len(duplicates))
        for dup, original in sorted(duplicates, key=lambda x: extract_number_from_filename(x[0])):
            print(f"{dup} (original: {original})")
            # Extract the number from the duplicate filename and append it
            bad_selection_numbers.append(extract_number_from_filename(dup))
            
            # Delete the duplicate file
            # os.remove(os.path.join(folder, dup))
            # print(f"Deleted duplicate file: {dup}")
    else:
        print("No duplicates found.")

    # Print the list of bad selection numbers
    print(f"bad_selection_numbers = {bad_selection_numbers}")

if __name__ == "__main__":
    main()
