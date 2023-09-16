# Copyright 2023, GITR contributors
# Authors: Abdou Diaw
"""
This test file is part of GITR (Global Impurity Transport code)

It solves the CIE equation, parse gitr output data (hisorty.nc) and plot the data.
"""
# -------
# Imports
# -------
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functools import partial

# Define a color map with distinct colors for each charge state
color_map = plt.cm.viridis(np.linspace(0, 1, 6))  # Here, 7 represents the number of charge states

def processGitrData(data_file, steps):
    # Read the data
    data = nc.Dataset(data_file)

    # Extract the charge and material data
    charge = data['charge'][:]
    Z = data['Z'][:]

    # Determine the number of time steps and unique material and  species
    nP = charge.shape[0]
    nT = charge.shape[1]

    unique_species = np.unique(charge)
    unique_materials = np.unique(Z)

    # Create a nested dictionary to store data for each material and species
    material_species_data = {
            material: {species: [] for species in unique_species} for material in unique_materials
        }

     # Calculate normalized number of particles for each species in each material at each time step
    for t in range(nT):
        for material in unique_materials:
            # Calculate the total number of species for this material at this time step
            nspecies = np.sum(Z[:, t] == material)
            for species in unique_species:
                condition = np.logical_and(charge[:, t] == species, Z[:, t] == material)
                count = np.sum(condition)
                normalized_count = count / nspecies
                material_species_data[material][species].append(normalized_count)

    # Create arrays to store the time
    time_steps = np.arange(nT) * steps

    return time_steps, material_species_data

def simulate_ionization_recombination(total_time, numberCharges, initial_conditions, material):
    # Store ionization and recombination rates from GITR interpolation of ADAS data
    ionization_rates = {
        74: [3.08754e+06, 969699, 386312, 206487, 35020.9], #16036.8],
        8: [289691, 29436.8, 5381.04, 759.26, 60.9562] #, 6.44085]  
    }

    recombination_rates = {
        74: [512.104, 1257.22, 1948.76, 2524.34, 3001.53], #, 3311.12],
        8: [33.7523,72.6351, 87.9312, 61.2116, 20.913 ] 

    }
    
    def ionization_rate(material, charge):
        if material in ionization_rates and charge in range(len(ionization_rates[material])):
            return ionization_rates[material][charge]
        else:
            raise ValueError(f"Invalid material {material} or charge {charge} for ionization rate.")

    def recombination_rate(material, charge):
        if material in recombination_rates and charge in range(1, len(recombination_rates[material]) + 1):
            return recombination_rates[material][charge-1]
        else:
            raise ValueError(f"Invalid material {material} or charge {charge} for recombination rate.")

    def rhs(nZ, t, material):
        dnZ_dt = np.zeros_like(nZ)
    
        for j in range(numberCharges):
            # For ionization (Z -> Z+1)
            if j < numberCharges - 1:
                rate = ionization_rate(material, j)
                dnZ_dt[j] -= nZ[j] * rate
                dnZ_dt[j+1] += nZ[j] * rate
                
            # For recombination (Z+1 -> Z)
            if j > 0:
                rate = recombination_rate(material, j)
                dnZ_dt[j] -= nZ[j] * rate
                dnZ_dt[j-1] += nZ[j] * rate
                    
        return dnZ_dt

    # Create an array for time values
    num_points = 1000
    time = np.linspace(0, 1, num_points) * total_time

    # Use odeint to solve the system
    nZ_values_flat = odeint(rhs, initial_conditions, time, args=(material,))

    # Reshape nZ_values_flat to the desired shape
    nZ_values = nZ_values_flat.T

    return time, nZ_values



# Set initial conditions for solver
total_time = 4e-4  # in seconds
dt = 1e-8     # in seconds
numberCharges = 5  # Number of charge states
initial_conditions = [1.0, 0.0, 0.0, 0.0, 0.0]

materials = [8, 74] # tungsten
# materials = ["O", "W"] # tungsten
# read in the data from the gitr output file
data_file = "output/history.nc"
time_steps, material_species_data = processGitrData(data_file, dt)

numberMaterials = len(materials)
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharex=False, sharey=True)  # 1D array of subplots
color_map = plt.get_cmap('tab10')

for idx, material in enumerate(materials):
    time, nZ_values = simulate_ionization_recombination(total_time, numberCharges, initial_conditions, material)

    for j in range(1, numberCharges):
        convert2millisec = 1e3
        
        # Adjust the indexing of axes for 1D
        ax = axes[idx]
        
        # Update the annotation to dynamically reflect the material
        if material == 8:
            materialName = 'O'
        elif material == 74:
            materialName = 'W'
        else:
            print('Invalid material')
        ion_name = f'{materialName}{j-1}+'

        if j == 1:
            ax.plot(time * convert2millisec, nZ_values[j-1, :], label=f'CIE: {ion_name}', color=color_map(j - 1), linewidth=2)
            ax.plot(time_steps * convert2millisec, material_species_data[material][j-1], linestyle='--', label=f'GITR: {ion_name}', color=color_map(j - 1), linewidth=2)
        else:
            ax.plot(time * convert2millisec, nZ_values[j-1, :], label=f'{ion_name}',  color=color_map(j - 1), linewidth=2)
            ax.plot(time_steps * convert2millisec, material_species_data[material][j-1], linestyle='--',color=color_map(j - 1), linewidth=2)

        ax.set_ylabel('$n_Z$ / $n_e$', fontsize=14)
        ax.grid(alpha=0.3)
        ax.legend()

        # Increase the tick font size
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_title(f'{materialName}', fontsize=14)
        if materialName == 'W':
            ax.set_xlim([0, 0.05])
    # Update the x-axis label of the plots
    axes[idx].set_xlabel('Time [ms]', fontsize=14)

plt.grid(alpha=0.3)
# Adjust this for better layout. Especially reduce the h_pad for less space between subplots.
plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.2)
plt.show()
