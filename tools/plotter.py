# Copyright 2023, GITR contributors
# Authors: T. Younkin: some of these scripts were adapted from the original Matlab scripts written by T. Younkin
# Update: A. Diaw: adapted the scripts to Python
"""
This test file is part of GITR (Global Impurity Transport code)

It tests the GITR CPC case (https://github.com/ORNL-Fusion/GITR-CPC-Example/tree/main/input)
for a single ion species, using the multi-ion capability of GITR.

Note: this test the case where all the ions are of the same species.
"""

import periodictable
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


def get_element_from_amu(amu_value):
    closest_element = None
    smallest_difference = float('inf')
    
    # Iterate through all elements to find the closest atomic mass
    for el in periodictable.elements:
        if hasattr(el, 'mass'):  # Ensure the element has the 'mass' attribute
            difference = abs(el.mass - amu_value)
            if difference < smallest_difference:
                smallest_difference = difference
                closest_element = el
                
    return closest_element

def plot_surface(surface_file_path):
    """
    Plot the surface.
    """
    filename = surface_file_path + "surface.nc"
    data = Dataset(filename, "r")
    edist = data.variables["surfEDist"][:]

    Elen = len(data.dimensions["nEnergies"])
    Alen = len(data.dimensions["nAngles"])
    Egrid = np.linspace(0, 1000, Elen)
    Agrid = np.linspace(0, 90, Alen)
    print("edist shape: ", edist.shape)
    edist = np.sum(edist, axis=0)
    
    c = plt.pcolormesh(Egrid, Agrid, edist.T, shading='auto')  # Use LogNorm for log colorscale
    plt.colorbar(c)
    plt.xlim(0,1000)
    plt.ylim(0,90)
    plt.xlabel("$\\theta$")
    plt.ylabel("E (eV)")
    plt.title("Summed Surface W Ion','Energy-Angle Distribution")
    plt.show()



def plot_particles(positions_file_path):
    """
    Plot the particles.
    """
    # Read data from netCDF file
    file_path = positions_file_path + 'positions.nc'
    with Dataset(file_path, 'r') as nc_file:
        x = nc_file.variables['x'][:]
        y = nc_file.variables['y'][:]
        z = nc_file.variables['z'][:]
        hitWall = nc_file.variables['hitWall'][:]
        amu_values = nc_file.variables['amu'][:]

    # Map AMU to element symbol
    elements = [get_element_from_amu(amu).symbol for amu in amu_values]

    # Identify unique elements
    unique_elements = list(set(elements))
    symbols = ['o', 's', '^', 'v', '*', '+', 'x', 'D', '<', '>']  
    colors_hit = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    colors_not_hit = ['pink', 'lightgreen', 'lightblue', 'lightcyan', 'lightcoral', 'lightskyblue', 'lightyellow']
    
    if len(unique_elements) > len(symbols):
        print(f"Warning: More unique elements than defined symbols. Some symbols will be reused.")
    
    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for idx, elem in enumerate(unique_elements):
        element_indices_hit = np.where(np.logical_and(np.array(elements) == elem, hitWall))[0]
        element_indices_not_hit = np.where(np.logical_and(np.array(elements) == elem, np.logical_not(hitWall)))[0]
        
        ax.scatter(x[element_indices_hit], y[element_indices_hit], z[element_indices_hit],
                   c=colors_hit[idx % len(colors_hit)],
                   marker=symbols[idx % len(symbols)],
                   label=f'{elem} Hit')
        
        ax.scatter(x[element_indices_not_hit], y[element_indices_not_hit], z[element_indices_not_hit],
                   c=colors_not_hit[idx % len(colors_not_hit)],
                   marker=symbols[idx % len(symbols)],
                   label=f'{elem} Not Hit')

    ax.set_xlim3d(min(x), max(x))
    ax.set_ylim3d(min(y), max(y))
    ax.set_zlim3d(min(z), max(z))    
    ax.set_title('Final Particle Positions')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()


