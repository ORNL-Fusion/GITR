makeCommand = "./makeGITRmac.sh"
import subprocess
from termcolor import colored
import io,libconf
import numpy as np
import os
from subprocess import check_output
import gitr

def buildGITR(examplePath="../examples/operatorTests/straightLine/2Dgeom"):
    env_file = "../env.mac.sh"
    filename = examplePath+"/input/gitrInput.cfg"
    ##Machine specific flags
    USE_CUDA=0
    USE_MPI=1
    USEMPI=0
    USE_OPENMP=0
    USE_BOOST=1
    
    cmake_flags=".. -DTHRUST_INCLUDE_DIR=/Users/tyounkin/code/thrust/ \
        -DCMAKE_C_COMPILER=gcc \
        -DCMAKE_CXX_COMPILER=g++ \
        -DNETCDF_CXX_INCLUDE_DIR=$NETCDFCXX4INCLUDE \
        -DNETCDF_CXX_LIBRARY=$NETCDFLIB_CPP \
        -DNETCDF_DIR=$NETCDFDIR \
        -DNETCDF_INCLUDE_DIR=$NETCDFINCLUDE \
        -DNETCDF_LIBRARY=$NETCDFLIB \
        -DNETCDF_CXX_INCLUDE_DIR=$NETCDFCXX4INCLUDE \
        -DLIBCONFIGPP_INCLUDE_DIR=/Users/tyounkin/Code/libconfigBuild/include \
        -DBoost_DIR=/Users/tyounkin/Code/boostBuild \
        -DBoost_INCLUDE_DIR=/Users/tyounkin/Code/boostBuild/include "
    
    with io.open(filename) as f:
        config = libconf.load(f)
        
    if(USE_CUDA!=config.flags.USE_CUDA):
        print("WARNING: Cuda flag mismatch between binary and input file\n")
    code_flags = "-DUSE_CUDA="+str(USE_CUDA)
    if(USE_MPI!=config.flags.USE_MPI):
        print("WARNING: USE_MPI flag mismatch between binary and input file\n")
    code_flags = code_flags+" -DUSE_MPI="+str(USE_MPI)
    if(USEMPI!=config.flags.USEMPI):
        print("WARNING: USEMPI flag mismatch between binary and input file\n")
    code_flags = code_flags+" -DUSEMPI="+str(USEMPI)
    if(USE_OPENMP!=config.flags.USE_OPENMP):
        print("WARNING: USE_OPENMP flag mismatch between binary and input file\n")
    code_flags = code_flags+" -DUSE_OPENMP="+str(USE_OPENMP)
    if(USE_BOOST!=config.flags.USE_BOOST):
        print("WARNING: USE_BOOST flag mismatch between binary and input file\n")
    code_flags = code_flags+" -DUSE_BOOST="+str(USE_BOOST)
    USEIONIZATION=config.flags.USEIONIZATION
    code_flags = code_flags+" -DUSEIONIZATION="+str(USEIONIZATION)
    USERECOMBINATION=config.flags.USERECOMBINATION
    code_flags = code_flags+" -DUSERECOMBINATION="+str(USERECOMBINATION)
    USEPERPDIFFUSION=config.flags.USEPERPDIFFUSION
    code_flags = code_flags+" -DUSEPERPDIFFUSION="+str(USEPERPDIFFUSION)
    USECOULOMBCOLLISIONS=config.flags.USECOULOMBCOLLISIONS
    code_flags = code_flags+" -DUSECOULOMBCOLLISIONS="+str(USECOULOMBCOLLISIONS)
    USETHERMALFORCE=config.flags.USETHERMALFORCE
    code_flags = code_flags+" -DUSETHERMALFORCE="+str(USETHERMALFORCE)
    USESURFACEMODEL=config.flags.USESURFACEMODEL
    code_flags = code_flags+" -DUSESURFACEMODEL="+str(USESURFACEMODEL)
    USESHEATHEFIELD=config.flags.USESHEATHEFIELD
    code_flags = code_flags+" -DUSESHEATHEFIELD="+str(USESHEATHEFIELD)
    BIASED_SURFACE=config.flags.BIASED_SURFACE
    code_flags = code_flags+" -DBIASED_SURFACE="+str(BIASED_SURFACE)
    USEPRESHEATHEFIELD=config.flags.USEPRESHEATHEFIELD
    code_flags = code_flags+" -DUSEPRESHEATHEFIELD="+str(USEPRESHEATHEFIELD)
    BFIELD_INTERP=config.flags.BFIELD_INTERP
    code_flags = code_flags+" -DBFIELD_INTERP="+str(BFIELD_INTERP)
    LC_INTERP=config.flags.LC_INTERP
    code_flags = code_flags+" -DLC_INTERP="+str(LC_INTERP)
    GENERATE_LC=config.flags.GENERATE_LC
    code_flags = code_flags+" -DGENERATE_LC="+str(GENERATE_LC)
    EFIELD_INTERP=config.flags.EFIELD_INTERP
    code_flags = code_flags+" -DEFIELD_INTERP="+str(EFIELD_INTERP)
    PRESHEATH_INTERP=config.flags.PRESHEATH_INTERP
    code_flags = code_flags+" -DPRESHEATH_INTERP="+str(PRESHEATH_INTERP)
    DENSITY_INTERP=config.flags.DENSITY_INTERP
    code_flags = code_flags+" -DDENSITY_INTERP="+str(DENSITY_INTERP)
    TEMP_INTERP=config.flags.TEMP_INTERP
    code_flags = code_flags+" -DTEMP_INTERP="+str(TEMP_INTERP)
    FLOWV_INTERP=config.flags.FLOWV_INTERP
    code_flags = code_flags+" -DFLOWV_INTERP="+str(FLOWV_INTERP)
    GRADT_INTERP=config.flags.GRADT_INTERP
    code_flags = code_flags+" -DGRADT_INTERP="+str(GRADT_INTERP)
    ODEINT=config.flags.ODEINT
    code_flags = code_flags+" -DODEINT="+str(ODEINT)
    FIXEDSEEDS=config.flags.FIXEDSEEDS
    code_flags = code_flags+" -DFIXEDSEEDS="+str(FIXEDSEEDS)
    PARTICLESEEDS=config.flags.PARTICLESEEDS
    code_flags = code_flags+" -DPARTICLESEEDS="+str(PARTICLESEEDS)
    GEOM_TRACE=config.flags.GEOM_TRACE
    code_flags = code_flags+" -DGEOM_TRACE="+str(GEOM_TRACE)
    GEOM_HASH=config.flags.GEOM_HASH
    code_flags = code_flags+" -DGEOM_HASH="+str(GEOM_HASH)
    GEOM_HASH_SHEATH=config.flags.GEOM_HASH_SHEATH
    code_flags = code_flags+" -DGEOM_HASH_SHEATH="+str(GEOM_HASH_SHEATH)
    PARTICLE_TRACKS=config.flags.PARTICLE_TRACKS
    code_flags = code_flags+" -DPARTICLE_TRACKS="+str(PARTICLE_TRACKS)
    PARTICLE_SOURCE_SPACE=config.flags.PARTICLE_SOURCE_SPACE
    code_flags = code_flags+" -DPARTICLE_SOURCE_SPACE="+str(PARTICLE_SOURCE_SPACE)
    PARTICLE_SOURCE_ENERGY=config.flags.PARTICLE_SOURCE_ENERGY
    code_flags = code_flags+" -DPARTICLE_SOURCE_ENERGY="+str(PARTICLE_SOURCE_ENERGY)
    PARTICLE_SOURCE_ANGLE=config.flags.PARTICLE_SOURCE_ANGLE
    code_flags = code_flags+" -DPARTICLE_SOURCE_ANGLE="+str(PARTICLE_SOURCE_ANGLE)
    PARTICLE_SOURCE_FILE=config.flags.PARTICLE_SOURCE_FILE
    code_flags = code_flags+" -DPARTICLE_SOURCE_FILE="+str(PARTICLE_SOURCE_FILE)
    SPECTROSCOPY=config.flags.SPECTROSCOPY
    code_flags = code_flags+" -DSPECTROSCOPY="+str(SPECTROSCOPY)
    USE3DTETGEOM=config.flags.USE3DTETGEOM
    code_flags = code_flags+" -DUSE3DTETGEOM="+str(USE3DTETGEOM)
    USECYLSYMM=config.flags.USECYLSYMM
    code_flags = code_flags+" -DUSECYLSYMM="+str(USECYLSYMM)
    USEFIELDALIGNEDVALUES=config.flags.USEFIELDALIGNEDVALUES
    code_flags = code_flags+" -DUSEFIELDALIGNEDVALUES="+str(USEFIELDALIGNEDVALUES)
    FLUX_EA=config.flags.FLUX_EA
    code_flags = code_flags+" -DFLUX_EA="+str(FLUX_EA)
    FORCE_EVAL=config.flags.FORCE_EVAL
    code_flags = code_flags+" -DFORCE_EVAL="+str(FORCE_EVAL)
    CHECK_COMPATIBILITY=config.flags.CHECK_COMPATIBILITY
    code_flags = code_flags+" -DCHECK_COMPATIBILITY="+str(CHECK_COMPATIBILITY)
    
    print(USE_CUDA)
    
    #print colored('Beginning clean','green')
    process = subprocess.Popen("./cleanGITRmake.sh", stdout=subprocess.PIPE)
    output, error = process.communicate()
    process.wait()
    p0 = subprocess.Popen("source "+env_file,shell=True, stdout=subprocess.PIPE)
    output, error = p0.communicate()
    p0.wait()
    env = {}
    output = check_output("source ../env.mac.sh", shell=True,executable="/bin/bash")
    env.update(os.environ)
    print(env)
    #print colored('Completed clean','green')
    cmake_command = "cmake " +cmake_flags+code_flags
    p1 = subprocess.Popen(cmake_command,shell=True,env=env, stdout=subprocess.PIPE)
    output, error = p1.communicate()
    p1.wait()
    #p2 = subprocess.Popen(makeCommand, stdout=subprocess.PIPE)
    #output, error = p2.communicate()
    #p2.wait()
    #print colored('Completed make of build files','green')
    p3 = subprocess.Popen("make", stdout=subprocess.PIPE)
    output, error = p3.communicate()
    p3.wait()
    #print colored('GITR successfully built','green')

def runGITR(examplePath="../examples/operatorTests/straightLine/2Dgeom"):
    os.chdir(examplePath)
    os.system('/Users/tyounkin/Code/gitr2/build/runGITR.sh')
def getAnswer(filename,x,y,r,z,charge):
    with io.open(filename) as f:
        config = libconf.load(f)
    relTol=1e-4	
    if(np.isclose(x,float(config.answer.x),rtol=relTol) and np.isclose(y,float(config.answer.y),rtol=relTol)	and np.isclose(z,float(config.answer.z),rtol=relTol) and np.isclose(r,float(config.answer.r),rtol=relTol) and np.isclose(charge,float(config.answer.charge),rtol=relTol)):
        print colored('Test passed','green')
	passed=1
    else:
        print colored('Test failed','red')
        print('x xAns',x,config.answer.x,np.isclose(x,float(config.answer.x),rtol=1e-04))
        print('y yAns',y,config.answer.y,y==float(config.answer.y))
        print('z zAns',z,config.answer.z,z==float(config.answer.z))
        print('r rAns',r,config.answer.r,r==float(config.answer.r))
        print('charge chargeAns',charge,config.answer.charge,charge==float(config.answer.charge))
	passed=0
	
    return passed	
if __name__ == "__main__":
    #cwd = os.getcwd()
    #buildGITR('../examples/operatorTests/straightLine/2Dgeom')
    #runGITR()
    #x,y,r,z,charge = gitr.nc_plotPositions('output/positions.nc')

    #print('x ',x[0])
    #getAnswer('input/answer.cfg',x[0],y[0],r[0],z[0],charge[0])
    #os.chdir(cwd)
    #buildGITR('../examples/operatorTests/straightLine/2DgeomCyl')
    #runGITR('../examples/operatorTests/straightLine/2DgeomCyl')
    #x,y,r,z,charge = gitr.nc_plotPositions('output/positions.nc')

    #print('x ',x[0])
    #getAnswer('input/answer.cfg',x[0],y[0],r[0],z[0],charge[0])
    #os.chdir(cwd)
    #buildGITR('../examples/operatorTests/straightLine/2DgeomPeriodic')
    #runGITR('../examples/operatorTests/straightLine/2DgeomPeriodic')
    #x,y,r,z,charge = gitr.nc_plotPositions('output/positions.nc')

    #print('x ',x[0])
    #getAnswer('input/answer.cfg',x[0],y[0],r[0],z[0],charge[0])
    #os.chdir(cwd)
    #buildGITR('../examples/operatorTests/straightLine/2DgeomCylPeriodic')
    #runGITR('../examples/operatorTests/straightLine/2DgeomCylPeriodic')
    #x,y,r,z,charge = gitr.nc_plotPositions('output/positions.nc')

    #print('x ',x[0])
    #getAnswer('input/answer.cfg',x[0],y[0],r[0],z[0],charge[0])
    #os.chdir(cwd)
    #buildGITR('../examples/operatorTests/straightLine/3Dgeom')
    #runGITR('../examples/operatorTests/straightLine/3Dgeom')
    #x,y,r,z,charge = gitr.nc_plotPositions('output/positions.nc')

    #print('x ',x[0])
    #getAnswer('input/answer.cfg',x[0],y[0],r[0],z[0],charge[0])
    #os.chdir(cwd)
    buildGITR('../examples/operatorTests/singleParticleMotion/gyroMotion')
    runGITR('../examples/operatorTests/singleParticleMotion/gyroMotion')
    x,y,z,r = gitr.nc_plotHist('output/history.nc')

    print('x ',np.min(x))
    print('y ',np.min(y))
    getAnswer('input/answer.cfg',np.min(x),np.min(y),0,0,0)
    os.chdir(cwd)
    #buildGITR('../iter/iter_milestone/3d')
