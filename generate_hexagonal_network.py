#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:31:22 2022

@author: narain

Tested in Python 3.7.4.

Generate the edge coordinates for a hexagonal network.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import errno
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import os, os.path

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to calculate the edge coordinates for the network
def calculate_polygons(startx, starty, endx, endy, sl):
    
    # Calculate the short diagonal
    short_diagonal = sl*math.sqrt(3)
    
    # Set the origin
    origx = startx
#    origy = starty

    # Offsets for moving along and up rows
    xoffset = sl+sl/2
    yoffset = short_diagonal/2

    # Initialise counters and arrays
    polygons = []
    row = 1
    counter = 0

    # Calculate the vertices
    while starty < endy:
        if row % 2 == 0:
            startx = origx + xoffset
        else:
            startx = origx
        while startx < endx:
            p1x = startx
            p1y = starty
            p2x = startx + sl
            p2y = starty
            p3x = startx + sl + sl/2
            p3y = starty + short_diagonal/2
            p4x = startx + sl
            p4y = starty + short_diagonal
            p5x = startx
            p5y = starty + short_diagonal
            p6x = startx - sl/2 
            p6y = starty + short_diagonal/2            
            poly_sides = (p1x,p1y,p2x,p2y),(p2x,p2y,p3x,p3y),(p3x,p3y,p4x,p4y),(p4x,p4y,p5x,p5y),(p5x,p5y,p6x,p6y),(p6x,p6y,p1x,p1y)
            polygons.extend(poly_sides)
            counter += 1
            startx += 3*sl
        starty += yoffset
        row += 1
    
    # Return the coordinates
    return polygons

# Define a function to group the edges with the same nodes
def normalize(edge):
    x1, y1, x2, y2 = edge[0],edge[1],edge[2],edge[3]
    if x1 > x2: 
        x1, x2 = x2, x1
        y1, y2 = y2, y1
    return x1, y1, x2, y2

# Define a function to make a directory if needed
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# Define a function to safely open a directory to write a file
def safe_open_w(path):
    mkdir_p(os.path.dirname(path))
    return open(path, 'w')

# =============================================================================
# DIMENSIONS
# =============================================================================

# Set the file name
file_name = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/MicrofluidicsNetworks/MerloEdgesMatrix_Offset.txt'
# file_name = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/MicrofluidicsNetworks/SymQuadMantegazzaEdgesMatrix_Offset.txt'

# Enter the domain specifications
# side_length = 100
side_length = 50
x_start = 0
x_end = 1000
y_start = 0
# y_start = -100
y_end = 1000

# Get a NumPy array with the edges
edge_matrix = calculate_polygons(x_start, y_start, x_end, y_end, side_length)
edge_matrix = np.array(edge_matrix)

# Tailor the network
#x_upper = 2050
#y_upper = 1900
x_lower = 60
x_upper = 590
y_upper = 490
edge_matrix = edge_matrix[np.logical_and(edge_matrix[:,0] >= x_lower, edge_matrix[:,2] >= x_lower)]
edge_matrix = edge_matrix[np.logical_and(edge_matrix[:,0] <= x_upper, edge_matrix[:,2] <= x_upper)]

# Set the y-offset
'''
radius =  side_length*math.sqrt(3)/2    
y_offset = radius
edge_matrix[:,1] += y_offset
edge_matrix[:,3] += y_offset
'''
edge_matrix = edge_matrix[np.logical_and(edge_matrix[:,1] <= y_upper, edge_matrix[:,3] <= y_upper)]

# Remove duplicate edges with the same nodes
for edge in edge_matrix:
    x1, y1, x2, y2 = edge[0],edge[1],edge[2],edge[3]
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1    
        edge[0],edge[1],edge[2],edge[3]=x1, y1, x2, y2
edge_matrix = edge_matrix.round(decimals=10)
edge_matrix = np.unique(edge_matrix, axis=0)

# Manually remove feeders for Mantegazza network
'''
# rows_to_remove = [0, 1, 2, 4, 5, 6, 7, 22, 183, 184, 185, 186, 188, 189, 190, 191]
# edge_matrix = np.delete(edge_matrix, rows_to_remove, axis=0)
'''

# Manually remove feeders for Merlo network
# '''
rows_to_remove = [0, 1, 3, 4, 5, 16, 33, 39, 50, 67, 73, 84, -7, -6, -5, -3, -2, -1]
edge_matrix = np.delete(edge_matrix, rows_to_remove, axis=0)
# '''


# Sort the edges
edge_matrix = edge_matrix[edge_matrix[:, 3].argsort()]
edge_matrix = edge_matrix[edge_matrix[:, 2].argsort(kind='mergesort')]  
edge_matrix = edge_matrix[edge_matrix[:, 1].argsort(kind='mergesort')] 
edge_matrix = edge_matrix[edge_matrix[:, 0].argsort(kind='mergesort')] 

# Set the x-offset to make sure inlet starts at x=0 (not totally sure why this error exists)
# '''
offset =  edge_matrix[0, 0]  
edge_matrix[:,0] -= offset
edge_matrix[:,2] -= offset
offset = side_length  
edge_matrix[:,1] += offset
edge_matrix[:,3] += offset
# '''

# =============================================================================
# PLOTS
# =============================================================================

# Plot the network
figure(figsize=(5, 5), dpi=80)
for x1, y1, x2, y2 in edge_matrix:
    plt.plot([x1, x2], [y1, y2], c='red')
    plt.xlim(-100,2000)
    plt.ylim(-100,2000)
plt.show()

# =============================================================================
# FILE
# =============================================================================
# '''
# Add vessel IDs
vessel_numbers = np.array(range(1,len(edge_matrix)+1))[...,None] #None keeps (n, 1) shape
edge_matrix = np.append(edge_matrix, vessel_numbers, 1)
# '''
# Write the file
with safe_open_w(file_name) as i_f:
    np.savetxt(file_name, edge_matrix, fmt='%1.4f \t %1.4f \t %1.4f \t %1.4f \t %i')
# '''