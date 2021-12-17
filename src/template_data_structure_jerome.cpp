/*
 * template_data_structure_jerome.cpp
 *
 *  Created on: Dec 16, 2021
 *      Author: jguterl
 */


/* Design and templating of data structure are proposed below in symbolic language.
 * Technical implementation technic are let at the discretion of Harry.
 * Class are considered thereby as container with attribute and method in a python fashion like.
 *  I believe that we should avoid inheretance for simplicity and usability by physicists
 *  Class object are therefore equivalent to c structures with methods



/*


/*
 * START TEMPLATE ENTRY POINT: Entry point for gitr */
int main{
GITR gitr(argc, argv) // entry point from inputfile
GITR gitr(libconfig:cfg) // entry point from wrapper with config file already given (can be interface with python easily)

}
/*
class GITR
	attributes:
 	 InputGITR inputgitr (==cfg possibly)
 	 GeomData geomdata
 	 SurfaceData surfacedata
 	 TrackerData trackerdata
 	 PlasmaData plasmadata
 	 SurfaceData surfacedata
 	 ParticleData particledata
 	 NumericsData numericsdata
 	 verbose // handler to print out messages during execution of the code
	methods:
	//constructor overloading
	 gitr(argc,argv){} // start from input file
	 gitr(libconfig::cfg) // start from cfg oject to allow a wrapper to control multiple instances of gitr


	 init()// allocate/initialize/broadcast
	 {
	 numericsdata->init(cfg->inputnumerics)
	 geomdata->init(cfg->inputgeom)
     particledata->init(cfg->inputparticle)
	 geomdata->init(cfg->inputgeom, particledata)
	 surfacedata->init(cfg->inputsurface, particledata)
	 plasmadata->init(cfg->inputplasma,geomdata,surfacedata,particledata)
	 trackerdata->init(cfg->inputtrackers,particledata)
	 }

	 broadcast() // broadcast mpi and sync cuda
	 check() //make sanity check (allocation, size) if needed
	 run()// driver for main loop
	 {

	 }


	 finalize() {}//dealloc +check memory leak + do physics checks (particle conservation)


	 */

/* END TEMPLATE ENTRY POINT */

/* START TEMPLATE INPUT: can be derived from libconfig object. No need to rewrite a class from scratch
 *
 * class InputGITR
 * {
 * attribute
 * InputNumerics inputnumerics // contain all numerics
 * InputPlasma inputplasma //contain plasma conditions
 * InputGeom inputgeom //contain all geom volume info
 * InputSurface inputsurface // contain all surface models (type particles, reflection + sputtering data needed)
 * InputParticles inputparticles // contain all particles info
 * methods
 * read_file() // obvious
 * write_file() //useful to create example and check validity of input and write into simulation folder to keep track of submitted inputs
 * verify() // check that all inputs are conform (Dimension, allowed values)
 *
 * }
 *
  END TEMPLATE INPUT */

/* START TEMPLATE GeomData: contain all geomdata
 * class GeomDATA:
 * 		attr
 * 			vectors with properties of volume elements -> can maybe do 2D and 3D class with class instance as attribute set with constructor
 * 		method
 * 			read_file() //obvious
 * 			write_file() //obvious
 * 			constructor() // read_File and allocate/instanciate vectors of properties
 *
 * END TEMPLATE GeomDATA/
 *
 *
/* START TEMPLATE PlasmaData: contain all plasma data, depends on datageom and datasurface and dataparticle to be sure that we have the right dataset e.g. for collisions between plasma species and impurities particles)
 *
 */
/* END TEMPLATE PlasmaData
 *
 */
/* START TEMPLATE SurfaceData: contain all surface data with model options and data
 *
 */
/* END TEMPLATE SurfaceData*/


/* START TEMPLATE ParticleData: contain all surface data with model options and data
 *
 */
/* END TEMPLATE ParticleData


/* START TEMPLATE TrackerData: contain all data top track particles (spectro+history+loss(not implemented)
 *
 */
/* END TEMPLATE TrackerData
 */

/* START TEMPLATE ParticlesData: contain all particle data: initial mass, type, charge, x,y,z, vx,vy,vz, chemical reactions, collisions options and transport options  if needed
 */
/* END TEMPLATE PlasmaData
