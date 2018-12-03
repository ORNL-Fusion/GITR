import os
import subprocess
from multiprocessing import Pool
import numpy as np
import sys

from mpi4py import MPI

def parProcessing(arg):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    name=MPI.Get_processor_name()
    text = 'name and rank and arg '+str(name)+' '+ str(rank)+' '+str(arg)
    return text
def runPool(nCall,nCores):
    p = Pool(nCores)
    print(p.map(parProcessing, nCall*np.linspace(0,nCores-1,nCores)))
if __name__ == '__main__':
    #parProcessing()
    nCall = int(sys.argv[1])
    nCores = int(sys.argv[2])
    runPool(nCall,nCores)
