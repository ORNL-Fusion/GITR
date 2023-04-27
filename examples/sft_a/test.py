#!/usr/bin/env python
import sys
import subprocess
import generate_inputs as inputs
import process_output as outputs
import numpy as np

def test():
    print("SFT case A Python test")
    inputs.generate()

    print("Running GITR")
    subprocess.run("../../build/GITR", shell=True, check=True)
    
    print("Post-processing")
    outputs.process_out()
    
    d0 = np.loadtxt("density_cuda.txt")
    d1 = np.loadtxt("density.txt")
    density_difference = d0[:,0] - d1[:,1]

    pct_diff = np.sum(density_difference)/np.sum(d0[:,0])
    print(pct_diff)
    if pct_diff < 0.01:
        sys.exit(0) #-1 for failed test
    else:
        sys.exit(-1) #-1 for failed test

    return -1 #doesn't affect test outcome

if __name__ == "__main__":
    test()
