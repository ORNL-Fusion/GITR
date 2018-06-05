from mpi4py import MPI
import generate_ftridyn_input
import analyze_ftridyn_simulations
#import ftridyn
import os
import time
import subprocess
import numpy as np
import netCDF4
import getopt, sys
import pickle
WORKTAG=0
DIETAG=1
def func1(path,E,a,r,d):
    print(path)
    cwd=os.getcwd()
    os.mkdir(path)
    os.chdir(path)
    cwd1=os.getcwd()
    print('cwd1 ',cwd1)
    name1 = ''
    name2 = ''
    if(len(d['beam'])>1):
        name1 = d['beam']
    else:
        name1 = d['beam']+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']

    generate_ftridyn_input.beam_and_target(name1+name2,d['beam'],d['target'],sim_number=1,number_histories=d['nH'], incident_energy=E,depth=200.0,incident_angle=a)
    p = subprocess.Popen([d['exe'],name1+name2+'0001.IN'],cwd=cwd+'/'+path)
    p.wait()
    #ftridyn.ftridyn_cpmi('He_W0001.IN')
    ##time.sleep(60)
    #print('finished simulation')
    #sim = analyze_ftridyn_simulations.ftridyn_output_data('He_W','0001')
    #Y = sim.calculate_total_sputtering_yield()
    #print('Energy Yield ',E, Y)
    os.chdir(cwd)
    cwd=os.getcwd()
    analyzeFTrun(path,d)
    print('cwd ',cwd)
    return 1
def analyzeFTrun(pathString,d):
    nEgrid = d['nEdist']
    maxE = d['maxEdist']
    nAgrid = d['nAdist']
    name1 = ''
    name2 = ''
    if(len(d['beam'])>1):
        name1 = d['beam']
    else:
        name1 = d['beam']+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']
    ffilename = name1+name2
    spyl_file = ffilename +'SPYL.DAT'
    print('doing analysis on ' , pathString+"/"+spyl_file)
    WW=analyze_ftridyn_simulations.ftridyn_output_data(pathString+"/"+ffilename,'0001')
    thisSpyl = WW.calculate_total_sputtering_yield()
    print('sputtering yield', thisSpyl)
    thisRefyl = WW.calculate_total_reflection_yield()
    np.savetxt(pathString+"/"+"YR.out",[thisSpyl,thisRefyl])
    print('reflection yield', thisRefyl)
    nP, nX, binsX,nY,binsY,nZ, binsZ = analyze_ftridyn_simulations.plot_sputtered_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nX.out",nX)
    nP, nX, binsX,nY,binsY,nZ, binsZ = analyze_ftridyn_simulations.plot_reflected_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nXref.out",nX)
    nPenergy, nenergy, binsenergy = analyze_ftridyn_simulations.plot_sputtered_energy_distributions(pathString+"/"+ffilename,nEgrid,maxE)
    np.savetxt(pathString+"/"+"energy.out",nenergy)
    nPenergyRef, nenergyRef, binsenergyRef = analyze_ftridyn_simulations.plot_reflected_energy_distributions(pathString+"/"+ffilename,nEgrid)
    np.savetxt(pathString+"/"+"energyRef.out",nenergyRef)

def func2(path,E,a,r):
     print('Path ear ', E, a , r)
class Work():
    def __init__(self,work_items):
        self.work_items = work_items[:] 
 
    def get_next_item(self):
        if len(self.work_items) == 0:
            return None
        return self.work_items.pop()
def master(nList):
    indxList = range(1,nList+1)
    all_data = []
    np = MPI.COMM_WORLD.Get_size()
    current_work = Work(indxList) 
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    size = min(np,nList)
    for i in range(1, size): 
        anext = current_work.get_next_item() 
        if not anext: break
        comm.send(anext-1, dest=i, tag=WORKTAG)
    
    while 1:
        anext = current_work.get_next_item()
        if not anext: break
        data = comm.recv(None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        all_data.append(data)
        comm.send(anext-1, dest=status.Get_source(), tag=WORKTAG)
 
    for i in range(1,size):
        data = comm.recv(None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
        all_data.append(data)
    
    for i in range(1,np):
        comm.send(None, dest=i, tag=DIETAG)
     
    return all_data
        
    
def slave(func,pathList,eList,aList,rList,d):
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(None, source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag(): break
        complete = func(pathList[data],eList[data],aList[data],rList[data],d)
        comm.send(complete, dest=0)

def main(func,argv):
    print 'Argument List:', str(argv)
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    comm = MPI.COMM_WORLD
    exit=0
    print('I am  %s rank %d (total %d)' % (name, rank, size) )
    loadDict = False
    dictName = None
    if(rank ==0):
        try:
            opts, args = getopt.getopt(argv, "d", ["dictionary="])
        except getopt.GetoptError as err:
            # print help information and exit:
            print str(err)  # will print something like "option -a not recognized"
            #usage()
	    exit=1
            #sys.exit(2)
        output = None
        verbose = False
	print("opts ", opts)
	if(exit==0):
            for o, a in opts:
                if o in ("-d", "--dictionary"):
                    #usage()
                    #sys.exit()
		    print("dictionary " , a)
		    loadDict = True
		    dictName = a
                else:
                    assert False, "unhandled option"
                    exit=1

    exit = comm.bcast(exit,root=0)
    if(exit):
        MPI.COMM_WORLD.Abort()
        sys.exit()
    pathList = []
    eList = []
    aList = []
    rList = []
    d = {}
    
    if rank == 0: # Master
        print('master rank',rank)
	if loadDict:
	    f = open(dictName, 'r')   # 'r' for reading; can be omitted
	    d = pickle.load(f)         # load file content as mydict
	    f.close()
	else:
            with open("ftMPI.in") as f:
                for line in f:
                    (key, val) = line.split()
                    if(key == 'beam' or key=='target' or key=='Escale' or key=='exe'):
                        d[key] = val
                    elif(key == 'nE' or key == 'nA' or key=='nR' or key=='nEdist' or key=='nAdist' or key=='nH'):
                        d[key] = int(val)
                    else:
                        d[key] = float(val)
        print d
    
        if(d['Escale'] == 'log'):
            energy = np.logspace(d['energy_start'],d['energy_end'],d['nE'])
        else:
            energy = np.linspace(d['energy_start'],d['energy_end'],d['nE'])
        angle = np.linspace(d['angle_start'],d['angle_end'],d['nA'])
        roughness = np.linspace(d['roughness_start'],d['roughness_end'],d['nR'])
    
        for i in range(d['nE']):
            for j in range(d['nA']):
                for k in range(d['nR']):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    pathList.append(pathString)
                    eList.append(energy[i])
                    aList.append(angle[j])
                    rList.append(roughness[k])
        
    pathList = comm.bcast(pathList,root=0)
    eList = comm.bcast(eList,root=0)
    aList = comm.bcast(aList,root=0)
    rList = comm.bcast(rList,root=0)
    d = comm.bcast(d,root=0)
    
    if rank == 0: # Master
        master(len(pathList))
    else: # Any slave
        print('slave rank',rank)
        slave(func,pathList,eList,aList,rList,d)
    
    comm.Barrier()
    print('Task completed (rank %d)' % (rank) )
    
    if rank == 0:
        nEgrid = d['nEdist']
        maxE = d['maxEdist']
        nAgrid = d['nAdist']
        sputt = np.zeros(shape=(len(energy),len(angle)))
        refl = np.zeros(shape=(len(energy),len(angle)))
        eDistEgrid = np.linspace(0.0,maxE-maxE/nEgrid,nEgrid) 
        cosDistAgrid = np.linspace(0.0,90.0-90.0/nAgrid,nAgrid) 
        cosXDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosZDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosXDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosZDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        eDistribution = np.zeros(shape=(len(energy),len(angle),nEgrid)) 
        eDistributionRef = np.zeros(shape=(len(energy),len(angle),nEgrid))
        totalIndex=0
        for i in range(len(energy)):
            for j in range(len(angle)):
                pathString = pathList[totalIndex]
                yr = np.loadtxt(pathString+"/"+"YR.out", dtype='float')
                sputt[i,j] = yr[0]
                refl[i,j] = yr[1]
                nX = np.loadtxt(pathString+"/"+"nX.out", dtype='float')
                cosXDistribution[i,j,:] = nX
                nXref = np.loadtxt(pathString+"/"+"nXref.out", dtype='float')
                cosXDistributionRef[i,j,:] = nXref
                nenergy = np.loadtxt(pathString+"/"+"energy.out", dtype='float')
                eDistribution[i,j,:] = nenergy
                nenergyRef = np.loadtxt(pathString+"/"+"energyRef.out", dtype='float')
                eDistributionRef[i,j,:] = nenergyRef
                totalIndex = totalIndex+1
        rootgrp = netCDF4.Dataset("ftridyn.nc", "w", format="NETCDF4")
        ne = rootgrp.createDimension("nE", len(energy))
        na = rootgrp.createDimension("nA", len(angle))
        nedistgrid = rootgrp.createDimension("nEdistBins", nEgrid)
        nadistgrid = rootgrp.createDimension("nAdistBins", nAgrid)
        spyld = rootgrp.createVariable("spyld","f8",("nE","nA"))
        rfyld = rootgrp.createVariable("rfyld","f8",("nE","nA"))
        ee = rootgrp.createVariable("E","f8",("nE"))
        aa = rootgrp.createVariable("A","f8",("nA"))
        cosxdist = rootgrp.createVariable("cosXDist","f8",("nE","nA","nAdistBins"))
        cosydist = rootgrp.createVariable("cosYDist","f8",("nE","nA","nAdistBins"))
        coszdist = rootgrp.createVariable("cosZDist","f8",("nE","nA","nAdistBins"))
        cosxdistref = rootgrp.createVariable("cosXDistRef","f8",("nE","nA","nAdistBins"))
        cosydistref = rootgrp.createVariable("cosYDistRef","f8",("nE","nA","nAdistBins"))
        coszdistref = rootgrp.createVariable("cosZDistRef","f8",("nE","nA","nAdistBins"))
        edist = rootgrp.createVariable("energyDist","f8",("nE","nA","nEdistBins"))
        edistref = rootgrp.createVariable("energyDistRef","f8",("nE","nA","nEdistBins"))
        edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins")) 
        cosdistagrid = rootgrp.createVariable("cosDistAgrid","f8",("nAdistBins")) 
        ee[:] = energy
        aa[:] = angle
        edistegrid[:] = eDistEgrid
        cosdistagrid[:] = cosDistAgrid
        spyld[:] = sputt
        rfyld[:] = refl
        cosxdist[:] = cosXDistribution
        cosydist[:] = cosYDistribution
        coszdist[:] = cosZDistribution
        cosxdistref[:] = cosXDistributionRef
        cosydistref[:] = cosYDistributionRef
        coszdistref[:] = cosZDistributionRef
        edist[:] = eDistribution
        edistref[:] = eDistributionRef
        rootgrp.close()

if __name__ == "__main__":
    print 'Argument List:', str(sys.argv[1:])
    main(func1,sys.argv[1:])
    #main(func2)
