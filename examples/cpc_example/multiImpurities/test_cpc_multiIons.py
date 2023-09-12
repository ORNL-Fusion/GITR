# Copyright 2023, GITR contributors
# Authors: Abdou Diaw
"""
This test file is part of GITR (Global Impurity Transport code)

It tests the GITR CPC case (https://github.com/ORNL-Fusion/GITR-CPC-Example/tree/main/input)
for a single ion species, using the multi-ion capability of GITR.

Note: this test the case where all the ions are of the same species.
"""
# -------
# Imports
# -------
import libconf
import io
from periodictable import tungsten
import subprocess as sp
import os
import shutil
from os.path import join
from plotter import plot_particles, plot_surface

#
def read_config(filename):
    """Reads file configuration."""
    with io.open(filename, 'r') as f:
        return libconf.load(f)

def write_config(filename, config):
    """Writes a configuration to file."""
    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

def update_species_properties(filename, impurity_species_list, run_info):
    """
    Updates species properties in the input file.
    """
    try:
        config = read_config(filename)
        species = config['impurityParticleSource']['initialConditions']['species']

        total_number_particles = []
        for i, impurity_species in enumerate(impurity_species_list):
                    species[i].update({
                        'impurity_Z': impurity_species.get('number', species[i]['impurity_Z']),
                        'charge': impurity_species.get('charge', species[i]['charge']),
                        'impurity_amu': impurity_species.get('mass', species[i]['impurity_amu']),
                        'energy_eV': impurity_species.get('energy', species[i]['energy_eV']),
                        'phi': impurity_species.get('phi', species[i]['phi']),
                        'theta': impurity_species.get('theta', species[i]['theta']),
                        'numberParticlesPerSpecies': impurity_species.get('numberParticlesPerSpecies', species[i]['numberParticlesPerSpecies']),
                        'x_start': impurity_species.get('x_start', species[i]['x_start']),
                        'y_start': impurity_species.get('y_start', species[i]['y_start']),
                        'z_start': impurity_species.get('z_start', species[i]['z_start'])
                    })

                    total_number_particles.append(impurity_species.get('numberParticlesPerSpecies', species[i]['numberParticlesPerSpecies']))
        # Update the number of particles and the run information
        config['timeStep']['dt'] = run_info['timeStep']
        config['timeStep']['nT'] = run_info['numberTimeSteps']
        if ( config['impurityParticleSource']['nP'] != sum(total_number_particles) ):
            print(f"Warning: the number of particles in the input file is not equal to the sum of the number of particles per species. The number of particles in the input file will be updated to {sum(total_number_particles)}")
            config['impurityParticleSource']['nP'] = sum(total_number_particles)
        write_config(filename, config)

    except (KeyError, IOError) as e:
        print(f"Error updating species properties: {e}")


def run_gitr(base_path):
    """
    Launch GITR using subprocess.Popen.
    """
    gitr_path = join(base_path, "GITR")
    output_path = join(os.curdir, "output")

    # Remove all files in the output directory
    for filename in os.listdir(output_path):
        file_path = join(output_path, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)  # Remove the file
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Remove the directory
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

    print("Running Simulation!")

    # Redirecting the output to the screen
    with sp.Popen([gitr_path], stdout=sp.PIPE, stderr=sp.PIPE, text=True) as proc:
        for line in proc.stdout:
            print(line, end='')  
        for line in proc.stderr:
            print(line, end='')


# Use of update_species_properties
materials = [tungsten, tungsten, tungsten]
numberParticlesPerSpecies = [333, 333, 333]

# Specifiy the impurity species properties: position, angle, energy are the same for all of them. They start all neutral.
base_species_data = {
    'charge': 0,
    'energy': 9.0,
    'phi': 0.0,
    'theta': 0.0,
    'x_start': 0.0,
    'y_start': 0.0,
    'z_start': 0.00001
}

# Create the list of impurity species
impurity_species_list = [
    {**base_species_data, 
     'number': mat.number, 
     'mass': mat.mass,
     'numberParticlesPerSpecies': npps} 
    for mat, npps in zip(materials, numberParticlesPerSpecies)
]
# run information time step and number of time steps
run_info = {
    'timeStep': 1.0e-9,
    'numberTimeSteps': 10000
}
# Update the input file
filename = "input/gitrInput.cfg"
update_species_properties(filename, impurity_species_list, run_info)

# Run GITR
HOME = os.environ['HOME']
base_path = join(HOME, "build")  # Path to the GITR executable
run_gitr(base_path)
print("Simulation finished!")
#EOF

# # Plot the results
path2data = "output/"

plot_particles(path2data)
plot_surface(path2data)
