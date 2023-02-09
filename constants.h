#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "io.h"

const double  h_epsilon = 1.19209e-07f;
const double pi        = 3.14159265358979324;

extern int   h_nx;    // abscissa grid dimension
extern int   h_ny;    // ordinate grid dimension
extern double h_dx;    // abscissa spatial resolution (m)
extern double h_dy;    // ordinate spatial resolution (m)

extern int    b_nx;
extern int    b_ny;
extern double b_xll;
extern double b_yll;

extern int GridSize;

extern double h_xll;
extern double h_yll;
extern double cellsize_original;
extern double cellsize;

extern size_t pitch, pitchBX, pitchBY;

#define DtoH cudaMemcpyDeviceToHost
#define HtoD cudaMemcpyHostToDevice

extern void CheckCUDAError(const char* message);

#endif
