#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "io.h"

const double  h_epsilon = 1.19209e-07f;
const double pi        = 3.14159265358979324;

extern int GridSize;
extern size_t GrowGridSize;

extern size_t pitch, pitchBX, pitchBY;
extern dim3 BlockDim,  GridDim;

#define DtoH cudaMemcpyDeviceToHost
#define HtoD cudaMemcpyHostToDevice

extern void CheckCUDAError(const char* message);

#endif
