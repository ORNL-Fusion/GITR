#!/usr/bin/env python
import sys
import subprocess
import generate_inputs as inputs
import process_output as outputs
import numpy as np
import matplotlib.pyplot as plt

def test():
    print("Hyperbolic Secant Python test")
    inputs.generate()

    print("Running GITR")
    subprocess.run("../../build/GITR", shell=True, check=True)
    
    print("Post-processing")
    z, dens, flux, temp = outputs.process_out()
    tau_s = 1.185*1.151179736678286e-06;

    ktm = 7.977106753667890e+07;
    
    c1 = 2*np.tanh(1)*tau_s*np.sqrt(ktm) +np.log(np.cosh(1));
    dz = z[1] - z[0];
    x = z;
    
    rho = (-np.log(np.cosh(x)) + c1)/(2*tau_s*ktm);
    a = ktm;
    b = (np.log(np.cosh(x)) - c1)/tau_s;
    c = np.tanh(x)*np.tanh(x);

    rho1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a);
    
    plt.plot(z+0.5*dz,dens,label='GITR')
    plt.plot(x,rho1,'--',label='Analytic')
    plt.title("Case A C Impurity Density")
    plt.xlabel("z [m]")
    plt.ylabel("n [m-3]")
    plt.legend()
    plt.show()
    plt.savefig('density.png')

    rhov = np.tanh(x);
    plt.close()
    plt.plot(z+0.5*dz,flux,label='GITR')
    plt.plot(x,rhov,'--',label='Analytic')
    plt.title("Case A C Impurity Flux")
    plt.xlabel("z [m]")
    plt.ylabel("nv [m-2s-1]")
    plt.legend()
    plt.show()
    plt.savefig('flux.png')
    
    plt.close()
    plt.plot(z+0.5*dz,temp/dens*2/3,label='GITR')
    #plt.plot(x,rhov,'--',label='Analytic')
    plt.title("Case A C Impurity Flux")
    plt.xlabel("z [m]")
    plt.ylabel("T [eV]")
    plt.legend()
    plt.show()
    plt.savefig('temp.png')
#    d0 = np.loadtxt("density_cuda.txt")
#    d1 = np.loadtxt("density.txt")
#    density_difference = d0[:,0] - d1[:,1]
#
#    pct_diff = np.sum(density_difference)/np.sum(d0[:,0])
#    print(pct_diff)
#    if pct_diff < 0.01:
#        sys.exit(0) #-1 for failed test
#    else:
#        sys.exit(-1) #-1 for failed test

    return -1 #doesn't affect test outcome

if __name__ == "__main__":
    test()
