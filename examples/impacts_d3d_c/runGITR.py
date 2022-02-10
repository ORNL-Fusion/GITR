import os
import subprocess
import impacts
import shutil
import numpy as np
target = np.loadtxt('Targets.txt', dtype='float',skiprows=1,delimiter=',')
print(target,target.shape)

for i in range(0,14):
    for j in range(1,7):
        ne = target[i,4]
        te = target[i,5]
        ti = target[i,6]
        vp = target[i,7]
        btot = target[i,8]
        br = target[i,9]
        bt = target[i,10]
        bz = target[i,11]
        impacts.d3d_case_C(charge = j, Tee=te, Tii = ti, n=ne,btot = btot,br = br, bt = bt, bz = bz, sheath_factor = 1.0)
        
        subprocess.run("/home/tqd/Code/GITR/build/GITR", shell=True, check=True)
        shutil.copyfile("output/surface.nc","surface_C"+str(j)+"_loc_"+str(i) +".nc")
        shutil.copyfile("output/positions.nc","positions_C"+str(j)+"_loc_"+str(i)+".nc")
        os.remove("output/surface.nc")
        os.remove("output/positions.nc")
        os.remove("output/positions.m")
        os.remove("output/particleSource.nc")
