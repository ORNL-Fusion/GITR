from mpi4py import MPI
import generate_ftridyn_input_gitr
import analyze_ftridyn_simulations_gitr
import ftridyn
import os
import time
import subprocess
import numpy as np
import netCDF4
import getopt, sys
import pickle
from cStringIO import StringIO
from collections import defaultdict
WORKTAG=0
DIETAG=1
def func1(path,E,a,r,d,specNum):
    print(path)
    cwd=os.getcwd()
    os.mkdir(path)
    os.chdir(path)
    #cwd1=os.getcwd()
    #print('cwd1 ',cwd1)
    beam = d['beam'][specNum]
    name1 = ''
    name2 = ''
    if(len(beam)>1):
        name1 = beam
    else:
        name1 = beam+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']

    generate_ftridyn_input_gitr.beam_and_target(name1+name2,beam,d['target'],sim_number=1,number_histories=d['nH'], incident_energy=E,depth=200.0,incident_angle=a,fluence=1.0E-16)
    #p = subprocess.Popen([d['exe'],name1+name2+'0001.IN'],cwd=cwd+'/'+path,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #stdoutdata, stderrdata = p.communicate()
    #returnCode = p.wait()
    #try:
    #    fileLog = open(path+'/log.txt','w') 
    #    fileErr = open(path+'/logErr.txt','w')
    #    fileIn = open(path+'/'+name1+name2+'0001.IN')
    #    p = subprocess.check_call([d['exe']],stdin=fileIn,cwd=cwd+'/'+path,stdout=fileLog, stderr=fileErr)
    #    fileLog.close()
    #    fileErr.close()
    #except subprocess.CalledProcessError as e:
    #    sys.exit("'ls' failed, returned code %d (check 'errors.txt')" \
    #            % (e.returncode))
    try:
        #original = sys.stdout
        f1 =  open('log.txt', 'w')
        f1.close()
        f2 =  open('logErr.txt', 'w')
        f2.close()
        null_fds = [os.open('log.txt',os.O_RDWR), os.open('logErr.txt',os.O_RDWR)]
        # save the current file descriptors to a tuple
        save = os.dup(1), os.dup(2)
        # put /dev/null fds on 1 and 2
        os.dup2(null_fds[0], 1)
        os.dup2(null_fds[1], 2)
        ftridyn.ftridyn_cpmi(name1+name2+'0001.IN')
        # restore file descriptors so I can print the results
        os.dup2(save[0], 1)
        os.dup2(save[1], 2)
        # close the temporary fds
        os.close(null_fds[0])
        os.close(null_fds[1])
        os.close(save[0])
        os.close(save[1])
        print 'ran tridyn'
        returnCode=0
    #sys.stdout = original
    #print 'ftridyn print', mystdout.getValue()
    except:
        print "FTRIDYN ERROR!"
        returnCode=0
    #print('path, returncode', path, p)
    #file = open(path+'/log.txt','w') 
    #file.write(str(stdoutdata)) 
    #file.close() 
    #file = open(path+'/logErr.txt','w') 
    #file.write(str(stderrdata)) 
    #file.close() 
    #ftridyn.ftridyn_cpmi('He_W0001.IN')
    ##time.sleep(60)
    #print('finished simulation')
    #sim = analyze_ftridyn_simulations.ftridyn_output_data('He_W','0001')
    #Y = sim.calculate_total_sputtering_yield()
    #print('Energy Yield ',E, Y)
    os.chdir(cwd)
    #cwd=os.getcwd()
    analyzeFTrun(path,d,specNum)
    print('cwd ',cwd)
    return returnCode
def analyzeFTrun(pathString,d,specNum):
    nEgrid = d['nEdist']
    maxE = d['maxEdist']
    nAgrid = d['nAdist']
    beam = d['beam'][specNum]
    name1 = ''
    name2 = ''
    if(len(beam)>1):
        name1 = beam
    else:
        name1 = beam+'_'
    if(len(d['target'])>1):
        name2 = d['target']
    else:
        name2 = '_'+d['target']
    ffilename = name1+name2
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
    np.savetxt(pathString+"/"+"nY.out",nY)
    #np.savetxt(pathString+"/"+"nZ.out",nZ)
    nP, nX, binsX,nY,binsY = analyze_ftridyn_simulations_gitr.plot_reflected_angular_distributions(pathString+"/"+ffilename,nAgrid)
    np.savetxt(pathString+"/"+"nXref.out",nX)
    np.savetxt(pathString+"/"+"nYref.out",nY)
    #np.savetxt(pathString+"/"+"nZref.out",nZ)
    nPenergy, nenergy, binsenergy = analyze_ftridyn_simulations_gitr.plot_sputtered_energy_distributions(pathString+"/"+ffilename,nEgrid,maxE)
    np.savetxt(pathString+"/"+"energy.out",nenergy)
    nPenergyRef, nenergyRef, binsenergyRef = analyze_ftridyn_simulations_gitr.plot_reflected_energy_distributions(pathString+"/"+ffilename,nEgrid)
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
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(None, source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag(): break
        else:
	    print('rank %d (total %d) received data %i' % ( rank, size,int(data)) )
	specNum = int(data)%d['nS']
        complete = func(pathList[data],eList[data],aList[data],rList[data],d,specNum)
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
        start_time = time.time()
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
    #d = {}
    d = defaultdict(list)
    
    if rank == 0: # Master
        print('master rank',rank)
        if loadDict:
            f = open(dictName, 'r')   # 'r' for reading; can be omitted
            d = pickle.load(f)         # load file content as mydict
            f.close()
        else:
            with open("ftMPI.in") as f:
                for line in f:
                    split = line.split()
                    if(split[0] == 'beam' or split[0]=='target' or split[0]=='Escale' or split[0]=='exe'):
		        for j in range(1,len(split)):
                            d[split[0]].append(split[j])
                    elif(split[0] == 'nE' or split[0] == 'nA' or split[0]=='nR' or split[0]=='nEdist' or split[0]=='nAdist' or split[0]=='nH'):
                        d[split[0]] = int(split[1])
                    else:
                        d[split[0]] = float(split[1])
        print d
        d['Escale'] = d['Escale'][0]
        d['target'] = d['target'][0]
        if(d['Escale'] == 'log'):
            energy = np.logspace(d['energy_start'],d['energy_end'],d['nE'])
        else:
            energy = np.linspace(d['energy_start'],d['energy_end'],d['nE'])
        angle = np.linspace(d['angle_start'],d['angle_end'],d['nA'])
        roughness = np.linspace(d['roughness_start'],d['roughness_end'],d['nR'])
        d['nS'] = len(d['beam'])
	print 'nS', d['nS'], d['beam']
	for i in range(d['nE']):
            for j in range(d['nA']):
                for k in range(d['nR']):
                    for ii in range(d['nS']):
                        pathString = "FTRIDYN_"+d['beam'][ii]+"_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                        pathList.append(pathString)
                        eList.append(energy[i])
                        aList.append(angle[j])
                        rList.append(roughness[k])
        print pathList
	f = open('pathList.pkl', 'w')
        pickle.dump(pathList,f)
        f.close() 
    comm.Barrier()
    pathList = comm.bcast(pathList,root=0)
    eList = comm.bcast(eList,root=0)
    aList = comm.bcast(aList,root=0)
    rList = comm.bcast(rList,root=0)
    d = comm.bcast(d,root=0)
    comm.Barrier()
    all_data = []
    if rank == 0: # Master
        print('master',rank)
        all_data = master(len(pathList))
    else: # Any slave
        print('slave rank',rank)
        slave(func,pathList,eList,aList,rList,d)
    
    print('Task waiting at barrier (rank %d)' % (rank) )
    comm.Barrier()
    print('Task completed (rank %d)' % (rank) )
    if(rank ==0):
        print('returns ' , all_data)
    
    if rank == 0:
        nEgrid = d['nEdist']
        maxE = d['maxEdist']
        nAgrid = d['nAdist']
        sputt = np.zeros(shape=(d['nS'],len(energy),len(angle)))
        refl = np.zeros(shape=(d['nS'],len(energy),len(angle)))
        eDistEgrid = np.linspace(0.0,maxE-maxE/nEgrid,nEgrid) 
        phiGrid = np.linspace(0.0,90.0-90.0/nAgrid,nAgrid) 
        thetaGrid = np.linspace(0.0,180.0-180.0/nAgrid,nAgrid) 
        cosXDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosYDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        #cosZDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosXDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        cosYDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        #cosZDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nAgrid)) 
        eDistribution = np.zeros(shape=(d['nS'],len(energy),len(angle),nEgrid)) 
        eDistributionRef = np.zeros(shape=(d['nS'],len(energy),len(angle),nEgrid))
        sputtSelf = np.zeros(shape=(len(energy),len(angle)))
        reflSelf = np.zeros(shape=(len(energy),len(angle)))
        cosXDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        #cosZDistributionSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosXDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        #cosZDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        eDistributionSelf = np.zeros(shape=(len(energy),len(angle),nEgrid)) 
        eDistributionRefSelf = np.zeros(shape=(len(energy),len(angle),nEgrid))
        totalIndex=0
        for i in range(len(energy)):
            for j in range(len(angle)):
                for ii in range(d['nS']):
		    specNum = ii%d['nS']
                    pathString = pathList[totalIndex]
                    yr = np.loadtxt(pathString+"/"+"YR.out", dtype='float')
                    sputt[specNum,i,j] = yr[0]
                    refl[specNum,i,j] = yr[1]
                    nX = np.loadtxt(pathString+"/"+"nX.out", dtype='float')
                    cosXDistribution[specNum,i,j,:] = nX
                    nY = np.loadtxt(pathString+"/"+"nY.out", dtype='float')
                    cosYDistribution[specNum,i,j,:] = nY
                    #nZ = np.loadtxt(pathString+"/"+"nZ.out", dtype='float')
                    #cosZDistribution[specNum,i,j,:] = nZ
                    nXref = np.loadtxt(pathString+"/"+"nXref.out", dtype='float')
                    cosXDistributionRef[specNum,i,j,:] = nXref
                    nYref = np.loadtxt(pathString+"/"+"nYref.out", dtype='float')
                    cosYDistributionRef[specNum,i,j,:] = nYref
                    #nZref = np.loadtxt(pathString+"/"+"nZref.out", dtype='float')
                    #cosZDistributionRef[specNum,i,j,:] = nZref
                    nenergy = np.loadtxt(pathString+"/"+"energy.out", dtype='float')
                    eDistribution[specNum,i,j,:] = nenergy
                    nenergyRef = np.loadtxt(pathString+"/"+"energyRef.out", dtype='float')
                    eDistributionRef[specNum,i,j,:] = nenergyRef
                    totalIndex = totalIndex+1
	#for i in range(d['nS']):	    
        rootgrp = netCDF4.Dataset("ftridynBackground"+".nc", "w", format="NETCDF4")
        ne = rootgrp.createDimension("nE", len(energy))
        na = rootgrp.createDimension("nA", len(angle))
        ns = rootgrp.createDimension("nS", d['nS']-1)
        nedistgrid = rootgrp.createDimension("nEdistBins", nEgrid)
        nadistgrid = rootgrp.createDimension("nAdistBins", nAgrid)
        spyld = rootgrp.createVariable("spyld","f8",("nS","nE","nA"))
        rfyld = rootgrp.createVariable("rfyld","f8",("nS","nE","nA"))
        ee = rootgrp.createVariable("E","f8",("nE"))
        aa = rootgrp.createVariable("A","f8",("nA"))
        cosxdist = rootgrp.createVariable("cosXDist","f8",("nS","nE","nA","nAdistBins"))
        cosydist = rootgrp.createVariable("cosYDist","f8",("nS","nE","nA","nAdistBins"))
        #coszdist = rootgrp.createVariable("cosZDist","f8",("nS","nE","nA","nAdistBins"))
        cosxdistref = rootgrp.createVariable("cosXDistRef","f8",("nS","nE","nA","nAdistBins"))
        cosydistref = rootgrp.createVariable("cosYDistRef","f8",("nS","nE","nA","nAdistBins"))
        #coszdistref = rootgrp.createVariable("cosZDistRef","f8",("nS","nE","nA","nAdistBins"))
        edist = rootgrp.createVariable("energyDist","f8",("nS","nE","nA","nEdistBins"))
        edistref = rootgrp.createVariable("energyDistRef","f8",("nS","nE","nA","nEdistBins"))
        edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins")) 
        phigrid = rootgrp.createVariable("phiGrid","f8",("nAdistBins")) 
        thetagrid = rootgrp.createVariable("thetaGrid","f8",("nAdistBins")) 
        ee[:] = energy
        aa[:] = angle
        edistegrid[:] = eDistEgrid
        phigrid[:] = phiGrid
        thetagrid[:] = thetaGrid
        spyld[:] = sputt[0:d['nS']-1,:,:]
        rfyld[:] = refl[0:d['nS']-1,:,:]
        cosxdist[:] = cosXDistribution[0:d['nS']-1,:,:,:]
        cosydist[:] = cosYDistribution[0:d['nS']-1,:,:,:]
        #coszdist[:] = cosZDistribution[0:d['nS']-1,:,:,:]
        cosxdistref[:] = cosXDistributionRef[0:d['nS']-1,:,:,:]
        cosydistref[:] = cosYDistributionRef[0:d['nS']-1,:,:,:]
        #coszdistref[:] = cosZDistributionRef[0:d['nS']-1,:,:,:]
        edist[:] = eDistribution[0:d['nS']-1,:,:,:]
        edistref[:] = eDistributionRef[0:d['nS']-1,:,:,:]
        rootgrp.close()
        rootgrp = netCDF4.Dataset("ftridynSelf"+".nc", "w", format="NETCDF4")
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
        phigrid = rootgrp.createVariable("phiGrid","f8",("nAdistBins")) 
        thetagrid = rootgrp.createVariable("thetaGrid","f8",("nAdistBins")) 
        phigrid[:] = phiGrid
        thetagrid[:] = thetaGrid
        ee[:] = energy
        aa[:] = angle
        edistegrid[:] = eDistEgrid
        spyld[:] = sputt[-1,:,:]
        rfyld[:] = refl[-1,:,:]
        cosxdist[:] = cosXDistribution[-1,:,:,:]
        cosydist[:] = cosYDistribution[-1,:,:,:]
        #coszdist[:] = cosZDistribution[-1,:,:,:]
        cosxdistref[:] = cosXDistributionRef[-1,:,:,:]
        cosydistref[:] = cosYDistributionRef[-1,:,:,:]
        #coszdistref[:] = cosZDistributionRef[-1,:,:,:]
        edist[:] = eDistribution[-1,:,:,:]
        edistref[:] = eDistributionRef[-1,:,:,:]
        rootgrp.close()
        exec_time = time.time()
        print("Execution of FTRIDYN Cases took --- %s seconds ---" % (exec_time - start_time))

if __name__ == "__main__":
    print 'Argument List:', str(sys.argv[1:])
    main(func1,sys.argv[1:])
    #main(func2)
