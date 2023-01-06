#include "Boundary.h"
//#include "Fields.h"
#include "Particles.h"
#include "Surfaces.h"
#include "array.h"
#include "boris.h"
#include "boundaryInit.h"
#include "coulombCollisions.h"
#include "crossFieldDiffusion.h"
#include "curandInitialize.h"
#include "fieldLineTrace.h"
#include "geometryCheck.h"
#include "hashGeom.h"
#include "hashGeomSheath.h"
#include "history.h"
#include "interp2d.hpp"
//#include "interpRateCoeff.hpp"
#include "interpolate.h"
#include "ionize.h"
#include <cmath>
//#include "ncFile.h"
#include "recombine.h"
#include "spectroscopy.h"
#include "surfaceModel.h"
//#include "testRoutineCuda.h"
#include "thermalForce.h"
#include "utils.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <netcdf.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "flags.hpp"


#if USE_DOUBLE
typedef double gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncDouble;
#else
typedef float gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncFloat;
#endif

readInput::File()
