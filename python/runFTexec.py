#!/usr/bin/env python
import os
import subprocess
#from multiprocessing import Pool
import numpy as np
import sys
import pickle
#import analyze_ftridyn_simulations
#import numpy as np

from mpi4py import MPI
def runFTexec(nCall,b,command,inFile):
    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        text1 = 'nCall b command inFile '+str(nCall)+' '+ str(b) + ' ' + command + ' '+inFile
        #print(text1)
        name=MPI.Get_processor_name()
        file_namelist = pickle.load( open( "ftridyn_file_namelist.pkl", "rb" ) )
        this_folder = "null"
        fileIndex = int(nCall*b+rank)
        text2 = 'rank size and fileIndex '+str(rank)+' '+str(size)+' '+ str(fileIndex)
        print(text2)
        sys.stdout.flush()
    except Exception as e:
        sys.stderr.write("Failed at preliminary EXCEPTION %s (check 'errors.txt')\n" % (str(e)))
    
    comm.Barrier()

    try:
        if( fileIndex < len(file_namelist)):
            this_folder = file_namelist[fileIndex]
            input_cmd_sh= command+ " " +inFile
            try:
                fileLog = open(this_folder+'/log.txt','w')
                fileErr = open(this_folder+'/logErr.txt','w')
                sys.stdout.write('input command %s and folder %s\n' %(input_cmd_sh, this_folder))
                sys.stdout.flush()
                ff = subprocess.check_output(input_cmd_sh,cwd=this_folder, shell=True)
                print('ff',ff)
                sys.stdout.flush()
                fileLog.write(ff)
                fileLog.close()
                fileErr.close()
            except subprocess.CalledProcessError as e:
                sys.stderr.write("'ls' failed, returned code %d (check 'errors.txt')\n" % (e.returncode))
    except Exception as e:
        sys.stderr.write("Exited at EXCEPTION %s (check 'errors.txt')\n" % (str(e)))
    text = 'name and rank and arg '+str(name)+' '+ str(rank) + ' ' + this_folder
    #print(text)
    return text
if __name__ == '__main__':
    #parProcessing()
    nCall = int(sys.argv[1])
    b = int(sys.argv[2])
    command = sys.argv[3]
    inFile = sys.argv[4]
    runFTexec(nCall,b,command,inFile)
