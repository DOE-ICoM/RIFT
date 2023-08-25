// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.cu
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-08-25 06:50:26 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>

#include "constants.h"
#include "io.h"
#include "GridSeries.cuh"




// -------------------------------------------------------------
//  class GridSeries
// -------------------------------------------------------------

// -------------------------------------------------------------
// GridSeries:: constructors / destructor
// -------------------------------------------------------------
GridSeries::GridSeries(const std::string& basename,
                       const double& scale,
                       const int& deltat,
                       const struct GridConfig& gc)
  : p_gc(gc),
    p_basename(basename),
    p_scale(scale),
    p_buffer(new double[gc.b_nx*gc.b_ny]),
    p_int_buffer(new double[gc.h_nx*gc.h_ny]),
    p_in_time(-9999.0), p_in_dt(deltat),
    p_current_dev(NULL)
{
  size_t width  = (GridDim.x * BlockDim.x) * sizeof(double);
  size_t height = (GridDim.y * BlockDim.y);

  checkCudaErrors(cudaMallocPitch((void**)&p_current_dev, &pitch, width, height));
  checkCudaErrors(cudaMemset2D(p_current_dev, pitch, 0, width, height));
}

GridSeries::~GridSeries(void)
{
  cudaFree(p_current_dev);
}

// -------------------------------------------------------------
// GridSeries::p_grid_name
// -------------------------------------------------------------
std::string
GridSeries::p_grid_name(const int& index) const
{
  std::stringstream hyetograph_file_ss;

  hyetograph_file_ss
    << p_basename << "-"
    << std::setw(3) << std::setfill('0') << index
    << ".txt";
  return hyetograph_file_ss.str();
}

// -------------------------------------------------------------
// GridSeries::p_read_grid
// -------------------------------------------------------------
void
GridSeries::p_read_grid(void)
{
  int index(trunc(p_in_time/p_in_dt));
  std::string fname(p_grid_name(index));
  
  SetOriginalGrid(p_buffer.get(), fname, p_gc);

  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      int jt = j - 2, it = i - 2;

      int ll = jt     * p_gc.b_nx + it;     // lower-left bathymetry point
      int ul = (jt+1) * p_gc.b_nx + it;     // upper-left bathymetry point
      int ur = (jt+1) * p_gc.b_nx + (it+1); // upper-right bathymetry point
      int lr = jt     * p_gc.b_nx + (it+1); // lower-right bathymetry point

      p_int_buffer[j*p_gc.h_nx+i] = 0.25f*(p_buffer[ur]+p_buffer[lr]+
                                           p_buffer[ll]+p_buffer[ul]);
      p_int_buffer[j*p_gc.h_nx+i] *= p_scale;
    }
  }
}

// -------------------------------------------------------------
// GridSeries::p_update
// -------------------------------------------------------------
void
GridSeries::p_update(const double& t)
{
  if (p_in_time < 0.0 || t > (p_in_time + p_in_dt)) {
    if (p_in_time < 0.0) {
      p_in_time = 0.0;
    } else {
      p_in_time += p_in_dt;
    }

    p_read_grid();
    checkCudaErrors(cudaMemcpy2D(p_current_dev, pitch, p_int_buffer.get(),
                                 p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                 p_gc.h_ny, HtoD));
  }
}

// -------------------------------------------------------------
// GridSeries::p_sum
// -------------------------------------------------------------
double
GridSeries::p_sum(void) const
{
  thrust::device_ptr<double> ptr(p_current_dev);
  double result = thrust::reduce(ptr, ptr + GridSize);
  return result;
}
  
