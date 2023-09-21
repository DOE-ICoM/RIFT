// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.cu
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-09-21 10:37:32 d3g096
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
                       const double& tmax,
                       const struct GridConfig& gc,
                       double *dev_buf)
  : p_gc(gc),
    p_basename(basename),
    p_scale(scale),
    p_buffer(new double[gc.b_nx*gc.b_ny]()),
    p_int_buffer(new double[gc.h_nx*gc.h_ny]()),
    p_in_time(-9999.0), p_in_dt(deltat), p_max_time(tmax),
    p_current_dev(dev_buf),
    p_external(true), 
    p_done(false)
{
  if (p_current_dev == NULL) {
    
    // warning: global variables
    // Call SetDeviceConstants() first
    size_t width  = (GridDim.x * BlockDim.x) * sizeof(double);
    size_t height = (GridDim.y * BlockDim.y);
    
    checkCudaErrors(cudaMallocPitch((void**)&p_current_dev, &pitch, width, height));
    checkCudaErrors(cudaMemset2D(p_current_dev, pitch, 0, width, height));

    p_external = true;
  }
}

GridSeries::~GridSeries(void)
{
  if (!p_external) 
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
GridSeries::p_read_grid(const int& index)
{
  std::string fname(p_grid_name(index));
  
  SetOriginalGrid(p_buffer.get(), fname, p_gc);

  std::cout << "Reading from " << fname << " ..." << std::endl;

  std::uninitialized_fill(p_int_buffer.get(),
                          p_int_buffer.get() + p_gc.h_nx*p_gc.h_ny,
                          0.0);

  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      int jt = j - 2, it = i - 2;

      int ll = jt     * p_gc.b_nx + it;     // lower-left bathymetry point
      int ul = (jt+1) * p_gc.b_nx + it;     // upper-left bathymetry point
      int ur = (jt+1) * p_gc.b_nx + (it+1); // upper-right bathymetry point
      int lr = jt     * p_gc.b_nx + (it+1); // lower-right bathymetry point

      int index(j*p_gc.h_nx+i);
      p_int_buffer[index] = 0.25f*(p_buffer[ur]+p_buffer[lr]+
                                   p_buffer[ll]+p_buffer[ul]);
      p_int_buffer[index] *= p_scale;
      // std::cout << i << ", " << j << ", " << p_int_buffer[index] << std::endl;
    }
  }
}

// -------------------------------------------------------------
// GridSeries::p_update
// -------------------------------------------------------------
void
GridSeries::p_update(const double& t)
{
  bool sendit(false);
  int index;
  
  if (p_in_time < 0.0) {
    p_in_time = t;
    index = trunc(p_in_time/p_in_dt);
    p_read_grid(index);
    sendit = true;
  }

  // after the maximum time is reached, just fill w/ zeroes
  
  if (t >= p_max_time) {
    if (!p_done) {
      std::uninitialized_fill(p_int_buffer.get(),
                              p_int_buffer.get() + p_gc.h_nx*p_gc.h_ny,
                              0.0);
      sendit = true;
    }
    p_done = true;
  } else if (t >= (p_in_time + p_in_dt)) {
    p_in_time += p_in_dt;
    index = trunc(p_in_time/p_in_dt);
    p_read_grid(index);
    sendit = true;
  }

  if (sendit) {
    checkCudaErrors(cudaMemcpy2D(p_current_dev, pitch, &(p_int_buffer[0]),
                                 p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                 p_gc.h_ny, HtoD));
  } 
}

// -------------------------------------------------------------
// GridSeries::p_sum
// -------------------------------------------------------------
extern __global__ void sumReduce(const size_t& pitch, double *x_dev, double *result);

double
GridSeries::p_sum(void) const
{
  double result(0.0);

  // FIXME: do this on the device, somehow
  
  // sumReduce <<< GridDim, BlockDim >>> (pitch, p_current_dev, &result);

  std::unique_ptr<double[]> tmp(new double[p_gc.h_nx*p_gc.h_ny]);
  checkCudaErrors(cudaMemcpy2D(tmp.get(), p_gc.h_nx*sizeof(double), p_current_dev,
                               pitch, p_gc.h_nx*sizeof(double), p_gc.h_ny, DtoH));

  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      result += tmp[j*p_gc.h_nx+i];
      // std::cout << i << ", " << j << ", " << tmp[j*p_gc.h_nx+i] << std::endl;
    }
  }
  
  return result;
}
  
// -------------------------------------------------------------
//  class InterpolatedGridSeries
// -------------------------------------------------------------

// -------------------------------------------------------------
// InterpolatedGridSeries:: constructors / destructor
// -------------------------------------------------------------
InterpolatedGridSeries::InterpolatedGridSeries(const std::string& basename,
                                               const double& scale,
                                               const int& deltat,
                                               const double& tmax,
                                               const struct GridConfig& gc)
  : GridSeries(basename, scale, deltat, tmax, gc),
    p_t0_buffer(new double[gc.h_nx*gc.h_ny]()),
    p_t1_buffer(new double[gc.h_nx*gc.h_ny]())
{}

InterpolatedGridSeries::~InterpolatedGridSeries(void)
{}

// -------------------------------------------------------------
// InterpolatedGridSeries::p_update
// -------------------------------------------------------------
void
InterpolatedGridSeries::p_update(const double& t)
{
  int index;

  if (p_in_time < 0.0) {
    
      p_in_time = t;
      index = trunc(p_in_time/p_in_dt);
      p_read_grid(index);
      std::copy(&(p_int_buffer[0]),
                &(p_int_buffer[0]) + p_gc.b_nx*p_gc.b_ny,
                &(p_t0_buffer[0]));

      p_in_time = index*p_in_dt;
      p_next_time = p_in_time + p_in_dt;
      p_read_grid(index + 1);
      std::copy(&(p_int_buffer[0]),
                &(p_int_buffer[0]) + p_gc.b_nx*p_gc.b_ny,
                &(p_t1_buffer[0]));
      
  } else if (t >= p_max_time) {
    
    p_in_time = p_next_time;
    p_next_time = p_in_time + p_in_dt;
    std::copy(&(p_t1_buffer[0]),
              &(p_t1_buffer[0]) + p_gc.b_nx*p_gc.b_ny,
              &(p_t0_buffer[0]));
    
  } else if (t >= p_next_time) {

    std::copy(&(p_t1_buffer[0]),
              &(p_t1_buffer[0]) + p_gc.b_nx*p_gc.b_ny,
              &(p_t0_buffer[0]));
    
    p_in_time = p_next_time;
    p_next_time = p_in_time + p_in_dt;
    index = trunc(p_next_time/p_in_dt);
    p_read_grid(index);
    std::copy(&(p_int_buffer[0]),
              &(p_int_buffer[0]) + p_gc.b_nx*p_gc.b_ny,
              &(p_t1_buffer[0]));
  }

  std::uninitialized_fill(p_int_buffer.get(),
                          p_int_buffer.get() + p_gc.h_nx*p_gc.h_ny,
                          0.0);

  double factor((t - p_in_time)/(p_next_time - p_in_time));
  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      int index(j*p_gc.h_nx+i);
      p_int_buffer[index] =
        factor*(p_t1_buffer[index] - p_t0_buffer[index]) +
        p_t0_buffer[index];
      // std::cout << i << ", " << j << ", "
      //           << p_t0_buffer[index] << ", "
      //           << p_t1_buffer[index] << ", "
      //           << p_int_buffer[index]
      //           << std::endl;
    }
  }

  checkCudaErrors(cudaMemcpy2D(p_current_dev, pitch, &(p_int_buffer[0]),
                               p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                               p_gc.h_ny, HtoD));

}  
