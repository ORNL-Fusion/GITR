#pragma once

#include <string>
#include <fstream>

#include "netcdf.h"
#include "ncFile.h"
#include "ncVar.h"
#include "ncDim.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

int readFileDim(const std::string& fileName,const std::string& varName);

int read_ar2Input( std::string fileName, gitr_precision *Bfield[]);
