#from libRustBCA.pybca import *
from libRustBCA import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from formulas import *
import time
import h5py
import netCDF4

def main():

##################
# INPUTS SECTION #
##################

    #scripts/materials.py has a number of potential ions and targets
    ion = tungsten
    target = tungsten

    # Make the ion Z slightly different than the target so that yield doesn't count reflections as well
    ion['Z'] = ion['Z'] + 1.0e-5
    #For smooth distributions and good statistics, you should use at least 10k ions
    number_ions = 10000

    # In this section, parameter ranges for scan and histograms are set
    nE = 5#20; # number of ion energies for building sputtering, reflection, and distributions
    nA = 4#18; # number of ion angles for building sputtering, reflection, and distributions

    E_min = 10.0; # Minimum ion energy for scan
    E_max = 1000.0; # Maximum ion energy for scan

    A_min = 0.0001; # Minimum angle energy for scan - avoid zero for numerical reasons
    A_max = 89.9; # Maximum angle energy for scan - avoid 90 (parallel to surface)
    
    # Number of bins for 3 1-D distributions of sputtered particles
    nE_dist = 100;
    nPhi_dist = 80;
    nTheta_dist = 90;

    E_dist_max = 100.0;

    # Number of bins for 3 1-D distributions of reflected particles
    nE_dist_ref = 100;
    nPhi_dist_ref = 80;
    nTheta_dist_ref = 90;
    
    E_dist_max_ref = E_max;
######################
# END INPUTS SECTION #
######################

    E_vector = np.logspace(np.log10(E_min),np.log10(E_max),nE)
    A_vector = np.linspace(A_min,A_max,nA)

    Edist_grid = np.linspace(0,E_dist_max,nE_dist+1);
    Phidist_grid = np.linspace(0,90.0,nPhi_dist+1);
    Thetadist_grid = np.linspace(-180.0,180.0,nTheta_dist+1);
    

    Edist_grid_ref = np.linspace(0,E_dist_max_ref,nE_dist_ref+1);
    Phidist_grid_ref = np.linspace(0,90.0,nPhi_dist_ref+1);
    Thetadist_grid_ref = np.linspace(-180.0,180.0,nTheta_dist_ref+1);

    Yields = np.zeros((nE,nA))
    Refls = np.zeros((nE,nA))
    
    Edist = np.zeros((nE,nA,nE_dist))
    Phidist = np.zeros((nE,nA,nPhi_dist))
    Thetadist = np.zeros((nE,nA,nTheta_dist))
    
    Edist_ref = np.zeros((nE,nA,nE_dist_ref))
    Phidist_ref = np.zeros((nE,nA,nPhi_dist_ref))
    Thetadist_ref = np.zeros((nE,nA,nTheta_dist_ref))

    for i in range(nE):
        for j in range(nA):

            en = E_vector[i]
            angle = A_vector[j]
            
            print("Simulating energy ", en, " and angle ", angle)
            
            energies_eV = en*np.ones(number_ions)

            #In RustBCA's 0D geometry, +x -> into the surface
            ux = np.cos(angle*np.pi/180.)*np.ones(number_ions)
            uy = np.sin(angle*np.pi/180.)*np.ones(number_ions)
            uz = np.zeros(number_ions)

            start = time.time()

            #Note that simple_bca_list_py expects number densities in 1/Angstrom^3
            output = np.array(simple_bca_list_py(energies_eV, ux, uy, uz, ion['Z'],
                ion['m'], ion['Ec'], ion['Es'], target['Z'], target['m'],
                target['Ec'], target['Es'], target['n']/10**30, target['Eb']))
            
            stop = time.time()
            delta_time = stop - start

            print('Simulation complete. Processing data...')

            Z = output[:, 0]
            m = output[:, 1]
            E = output[:, 2]
            x = output[:, 3]
            y = output[:, 4]
            z = output[:, 5]
            ux = output[:, 6]
            uy = output[:, 7]
            uz = output[:, 8]

            #For the python bindings, these conditionals can be used to distinguish
            #between sputtered, reflected, and implanted particles in the output list
            sputtered = output[np.logical_and(Z == target['Z'], E > 0), :]
            reflected = output[np.logical_and(Z == ion['Z'], x < 0), :]
            implanted = output[np.logical_and(Z == ion['Z'], x > 0), :]

            MI = 1.6737236e-27;
            Q = 1.60217662e-19;

            # Get sputtered velocities
            ux_sputtered = np.sqrt(2*Q*sputtered[:,2]/target['m']/MI)*sputtered[:,7]
            uy_sputtered = np.sqrt(2*Q*sputtered[:,2]/target['m']/MI)*sputtered[:,8]
            uz_sputtered = -np.sqrt(2*Q*sputtered[:,2]/target['m']/MI)*sputtered[:,6]

            ur_sputtered = np.sqrt(np.multiply(ux_sputtered,ux_sputtered) + np.multiply(uy_sputtered,uy_sputtered))
            phi = np.arccos(np.divide(uz_sputtered,ur_sputtered));
            theta = np.arctan(np.divide(uy_sputtered,ux_sputtered));

            # Calculate distributions for sputtered ions
            hist, bin_edges= np.histogram(sputtered[:,2], Edist_grid)
            Edist[i,j,:] = hist
            plt.figure("Edist", figsize=(10, 6), dpi=2000)
            plt.plot(cell_centers(Edist_grid),hist)
            hist, bin_edges = np.histogram(phi, Phidist_grid)
            Phidist[i,j,:] = hist
            hist, bin_edges = np.histogram(theta, Thetadist_grid)
            Thetadist[i,j,:] = hist
            
            # Get reflected velocities 
            ux_reflected = np.sqrt(2*Q*reflected[:,2]/target['m']/MI)*reflected[:,7]
            uy_reflected = np.sqrt(2*Q*reflected[:,2]/target['m']/MI)*reflected[:,8]
            uz_reflected = -np.sqrt(2*Q*reflected[:,2]/target['m']/MI)*reflected[:,6]

            ur_reflected = np.sqrt(np.multiply(ux_reflected,ux_reflected) + np.multiply(uy_reflected,uy_reflected))

            phi_ref = np.arccos(np.divide(uz_reflected,ur_reflected));
            theta_ref = np.arctan(np.divide(uy_reflected,ux_reflected));
            
            # Make histograms for reflected ions
            hist, bin_edges= np.histogram(reflected[:,2], Edist_grid_ref)
            Edist_ref[i,j,:] = hist
            hist, bin_edges = np.histogram(phi_ref, Phidist_grid_ref)
            Phidist_ref[i,j,:] = hist
            hist, bin_edges = np.histogram(theta_ref, Thetadist_grid_ref)
            Thetadist_ref[i,j,:] = hist
            
            hf = h5py.File('velocities_sputtered_'+str(i)+'_'+str(j) +'.h5', 'w')
            hf.create_dataset('vx', data=ux_sputtered)
            hf.create_dataset('vy', data=uy_sputtered)
            hf.create_dataset('vz', data=uz_sputtered)
            hf.close()

            Yields[i,j] = 1.0*len(ux_sputtered)/number_ions
            Refls[i,j] = 1.0*len(ux_reflected)/number_ions

    # Write GITR surface model file
    rootgrp = netCDF4.Dataset("RustBCASelf"+".nc", "w", format="NETCDF4")

    ne = rootgrp.createDimension("nE", nE)
    na = rootgrp.createDimension("nA", nA)
    nedistgrid = rootgrp.createDimension("nEdistBins", nE_dist)
    nedistgridref = rootgrp.createDimension("nEdistBinsRef", nE_dist_ref)
    nphidistgrid =    rootgrp.createDimension("nPhidistBins", nPhi_dist)
    nphidistgridref = rootgrp.createDimension("nPhidistBinsRef", nPhi_dist_ref)
    nthetadistgrid = rootgrp.createDimension("nThetadistBins", nTheta_dist)
    nthetadistgridref = rootgrp.createDimension("nThetadistBinsRef", nTheta_dist_ref)

    ee = rootgrp.createVariable("E","f8",("nE"))
    aa = rootgrp.createVariable("A","f8",("nA"))
    
    edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins"))
    phigrid = rootgrp.createVariable("phiGrid","f8",("nPhidistBins"))
    thetagrid = rootgrp.createVariable("thetaGrid","f8",("nThetadistBins"))
    edistegrid_ref = rootgrp.createVariable("eDistEgrid_ref","f8",("nEdistBins"))
    phigrid_ref = rootgrp.createVariable("phiGrid_ref","f8",("nPhidistBins"))
    thetagrid_ref = rootgrp.createVariable("thetaGrid_ref","f8",("nThetadistBins"))
    
    spyld = rootgrp.createVariable("spyld","f8",("nE","nA"))
    rfyld = rootgrp.createVariable("rfyld","f8",("nE","nA"))
    edist = rootgrp.createVariable("energyDist","f8",("nE","nA","nEdistBins"))
    phidist = rootgrp.createVariable("phiDist","f8",("nE","nA","nPhidistBins"))
    thetadist = rootgrp.createVariable("thetaDist","f8",("nE","nA","nThetadistBins"))
    edist_ref = rootgrp.createVariable("energyDistRef","f8",("nE","nA","nEdistBinsRef"))
    phidist_ref = rootgrp.createVariable("phiDistRef","f8",("nE","nA","nPhidistBinsRef"))
    thetadist_ref = rootgrp.createVariable("thetaDistRef","f8",("nE","nA","nThetadistBinsRef"))
    
    ee[:] = E_vector
    aa[:] = A_vector
    edistegrid[:] = cell_centers(Edist_grid)
    edistegrid_ref[:] = cell_centers(Edist_grid_ref)
    phigrid[:] = cell_centers(Phidist_grid)
    thetagrid[:] = cell_centers(Thetadist_grid)
    phigrid_ref[:] = cell_centers(Phidist_grid_ref)
    thetagrid_ref[:] = cell_centers(Thetadist_grid_ref)

    spyld[:] = Yields
    rfyld[:] = Refls
    edist[:] = Edist
    phidist[:] = Phidist
    thetadist[:] = Thetadist
    edist_ref[:] = Edist_ref
    phidist_ref[:] = Phidist_ref
    thetadist_ref[:] = Thetadist_ref
    rootgrp.close()

    X, Y = np.meshgrid(A_vector, E_vector)

    print(X.shape)
    print(Yields.shape)
    print(Y.shape)
    plt.figure("Yields")
    plt.pcolor(X,Y , Yields,shading='nearest',norm=matplotlib.colors.LogNorm()) #norm=LogNorm(vmin=1e-3, vmax=Yields.max()), cmap=plt.cm.autumn)
    plt.colorbar()
    plt.show()
    plt.savefig('Yields.png')
    
    plt.figure("Edist")
    plt.title("Edist")
    plt.xlabel("E [eV]")
    plt.ylabel("f [Counts]")
    plt.show()
    plt.savefig('Edist.png')
def cell_centers(edges):
    return 0.5*(edges[0:-1] + edges[1:])
if __name__ == '__main__':
    main()
