#!/usr/bin/env python
import sys
import subprocess
import generate_inputs as inputs
import process_output as outputs
import numpy as np
import matplotlib.pyplot as plt
import io
import libconf

def test():
    print("Hyperbolic Secant Python test")
    inputs.generate()

# Run 3 different coulomb collision models
# 1 - callen
    filename = "input/gitrInput.cfg";

    with io.open(filename) as f:
        config = libconf.load(f)
        config['flags']['USECOULOMBCOLLISIONS'] = 1
    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

    print("Running GITR")
    subprocess.run("../../build/GITR", shell=True, check=True)
    
    print("Post-processing")
    z, dens1, flux1, temp1 = outputs.process_out()
# Boris buneman
    with io.open(filename) as f:
        config = libconf.load(f)
        config['flags']['USECOULOMBCOLLISIONS'] = 2
    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

    print("Running GITR")
    subprocess.run("../../build/GITR", shell=True, check=True)
    
    print("Post-processing")
    z, dens2, flux2, temp2 = outputs.process_out()

# Reiser
    with io.open(filename) as f:
        config = libconf.load(f)
        config['flags']['USECOULOMBCOLLISIONS'] = 3
    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

    print("Running GITR")
    subprocess.run("../../build/GITR", shell=True, check=True)
    
    print("Post-processing")
    z, dens3, flux3, temp3 = outputs.process_out()

    tau_s = 1.185*1.151179736678286e-06;
    tau_s = 1.151179736678286e-06;

    ktm = 7.977106753667890e+07;
    
    c1 = 2*np.tanh(1)*tau_s*np.sqrt(ktm) +np.log(np.cosh(1));
    dz = z[1] - z[0];
    x = z;
    
    rho = (-np.log(np.cosh(x)) + c1)/(2*tau_s*ktm);
    a = ktm;
    b = (np.log(np.cosh(x)) - c1)/tau_s;
    c = np.tanh(x)*np.tanh(x);

    rho1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a);
    
    d0 = np.loadtxt("dens_flux_temp0.txt")
    density_difference1 = d0[:,1] - dens1
    flux_difference1 = d0[:,2] - flux1
    temp_difference1 = d0[:,3] - temp1
    density_difference2 = d0[:,1] - dens2
    flux_difference2 = d0[:,2] - flux2
    temp_difference2 = d0[:,3] - temp2
    density_difference3 = d0[:,1] - dens3
    flux_difference3 = d0[:,2] - flux3
    temp_difference3 = d0[:,3] - temp3
    
    plt.plot(z+0.5*dz,dens1,label='GITR1')
    plt.plot(z+0.5*dz,dens2,label='GITR2')
    plt.plot(z+0.5*dz,dens3,label='GITR3')
    plt.plot(z+0.5*dz,d0[:,1],label='GITR_gold')
    plt.plot(x,rho1,'--',label='Analytic')
    plt.title("Case A C Impurity Density")
    plt.xlabel("z [m]")
    plt.ylabel("n [m-3]")
    plt.legend()
    plt.show()
    plt.savefig('density.png')

    rhov = np.tanh(x);
    plt.close()
    plt.plot(z+0.5*dz,flux1,label='GITR1')
    plt.plot(z+0.5*dz,flux2,label='GITR2')
    plt.plot(z+0.5*dz,flux3,label='GITR3')
    plt.plot(z+0.5*dz,d0[:,2],label='GITR_gold')
    plt.plot(x,rhov,'--',label='Analytic')
    plt.title("Case A C Impurity Flux")
    plt.xlabel("z [m]")
    plt.ylabel("nv [m-2s-1]")
    plt.legend()
    plt.show()
    plt.savefig('flux.png')
    
    plt.close()
    plt.plot(z+0.5*dz,temp1/dens1*2/3,label='GITR1')
    plt.plot(z+0.5*dz,temp2/dens2*2/3,label='GITR2')
    plt.plot(z+0.5*dz,temp3/dens3*2/3,label='GITR3')
    plt.plot(z+0.5*dz,d0[:,3]/d0[:,1]*2/3,label='GITR_gold')
    plt.plot(z+0.5*dz,0*temp1 + 10,'--',label='Analytic')
    plt.title("Case A C Impurity Flux")
    plt.xlabel("z [m]")
    plt.ylabel("T [eV]")
    plt.legend()
    plt.show()
    plt.savefig('temp.png')

    pct_diff_dens1 = np.abs(np.sum(density_difference1)/np.sum(d0[:,1]))
    pct_diff_flux1 = np.abs(np.sum(flux_difference1)/np.sum(d0[:,2]))
    pct_diff_temp1 = np.abs(np.sum(temp_difference1)/np.sum(d0[:,3]))
    print(pct_diff_dens1, pct_diff_flux1, pct_diff_temp1)
    pct_diff_dens2 = np.abs(np.sum(density_difference2)/np.sum(d0[:,1]))
    pct_diff_flux2 = np.abs(np.sum(flux_difference2)/np.sum(d0[:,2]))
    pct_diff_temp2 = np.abs(np.sum(temp_difference2)/np.sum(d0[:,3]))
    print(pct_diff_dens2, pct_diff_flux2, pct_diff_temp2)
    pct_diff_dens3 = np.abs(np.sum(density_difference3)/np.sum(d0[:,1]))
    pct_diff_flux3 = np.abs(np.sum(flux_difference3)/np.sum(d0[:,2]))
    pct_diff_temp3 = np.abs(np.sum(temp_difference3)/np.sum(d0[:,3]))
    print(pct_diff_dens3, pct_diff_flux3, pct_diff_temp3)

    # baseline differences 0.011227382556773153 0.0063831482083862336 0.011200203885744164
    if ((pct_diff_dens1 < 0.013) and (pct_diff_flux1 < 0.008) and (pct_diff_temp1 < 0.013)):
        sys.exit(0) #-1 for failed test
    else:
        sys.exit(-1) #-1 for failed test

    return -1 #doesn't affect test outcome

if __name__ == "__main__":
    test()
