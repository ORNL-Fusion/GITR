import os
import subprocess
from multiprocessing import Pool
import numpy as np
import sys
import pickle
import analyze_ftridyn_simulations_gitr
import numpy as np

from mpi4py import MPI
def analyzeFTrun(pathString):
    nEgrid = 100
    maxE = 100.0
    nAgrid = 50
    ffilename = 'W_W_'
    spyl_file = ffilename +'SPYL.DAT'
    print('doing analysis on ' , pathString+"/"+spyl_file)
    WW=analyze_ftridyn_simulations_gitr.ftridyn_output_data(pathString+"/"+ffilename,'0001')
    thisSpyl = WW.calculate_total_sputtering_yield()
    print('sputtering yield', thisSpyl)
    thisRefyl = WW.calculate_total_reflection_yield()
    np.savetxt(pathString+"/"+"YR.out",[thisSpyl,thisRefyl])
    print('reflection yield', thisRefyl)
    nP, nX, binsX,nY,binsY = analyze_ftridyn_simulations_gitr.plot_sputtered_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nX.out",nX)
    nP, nX, binsX,nY,binsY = analyze_ftridyn_simulations_gitr.plot_reflected_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nXref.out",nX)
    nPenergy, nenergy, binsenergy = analyze_ftridyn_simulations_gitr.plot_sputtered_energy_distributions(pathString+"/"+ffilename,nEgrid)
    np.savetxt(pathString+"/"+"energy.out",nenergy)
    nPenergyRef, nenergyRef, binsenergyRef = analyze_ftridyn_simulations_gitr.plot_reflected_energy_distributions(pathString+"/"+ffilename,nEgrid)
    np.savetxt(pathString+"/"+"energyRef.out",nenergyRef)
def runFT(nCall,b):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    name=MPI.Get_processor_name()
    file_namelist = pickle.load( open( "ftridyn_file_namelist.pkl", "rb" ) )
    this_folder = "null"
    fileIndex = int(nCall*b+rank)
    if( fileIndex < len(file_namelist)):
        this_folder = file_namelist[fileIndex]
        analyzeFTrun(this_folder)
    text = 'name and rank and arg '+str(name)+' '+ str(rank) + ' ' + this_folder
    print(text)
    return text
if __name__ == '__main__':
    #parProcessing()
    nCall = int(sys.argv[1])
    b = int(sys.argv[2])
    runFT(nCall,b)
