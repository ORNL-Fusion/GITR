import gitr
import numpy as np

if __name__ == "__main__":
    GITR_ROOT='/Users/tyounkin/Code/gitr2/'
    x,y,z,r,vx,vy,vz,charge,weight,nP,nT = gitr.nc_plotHist(GITR_ROOT+'examples/operatorTests/ionization/output/history.nc',0)
    nNeutral = np.array(nT)
    for i in range(nT):
        c1 = charge[:][i];
	cond = np.find(c1==0);
	nNeutral[i]=cond.length

        
