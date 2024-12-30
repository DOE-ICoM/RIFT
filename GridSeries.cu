// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.cu
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2024-12-30 11:22:55 d3g096
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
#include "grid.h"
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
    p_done(false),
    p_allow_nodata(false),
    p_current_dev_init(true)
{
  if (p_current_dev == NULL) {
    
    // warning: global variables
    // Call SetDeviceConstants() first
    size_t width  = (GridDim.x * BlockDim.x) * sizeof(double);
    size_t height = (GridDim.y * BlockDim.y);
    
    checkCudaErrors(cudaMallocPitch((void**)&p_current_dev, &pitch, width, height));
    p_init_dev();

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
// GridSeries::p_interp
// -------------------------------------------------------------
void
GridSeries::p_interp(void)
{
  std::uninitialized_fill(p_int_buffer.get(),
                          p_int_buffer.get() + p_gc.h_nx*p_gc.h_ny,
                          (p_allow_nodata ? p_gc.nodata : 0.0));

  if (!p_current_dev_init) p_init_dev();

  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      int jt = j - 2, it = i - 2;
      
      int nnd(0);
      double vsum(0.0);
      for (int ni = 0; ni < 2; ++ni) {
        for (int nj = 0; nj < 2; ++nj) {
          int idx = (jt + nj) * p_gc.b_nx + (it + ni);
          if (p_buffer[idx] != p_gc.nodata) {
            nnd++;
            vsum += p_buffer[idx];
          }
        }
      }

      int index(j*p_gc.h_nx+i);
      if (nnd > 2) {
        p_int_buffer[index] = vsum/((double)nnd);
        p_int_buffer[index] *= p_scale;
      } 
    }
  }
}

// -------------------------------------------------------------
// GridSeries::p_read_grid
// -------------------------------------------------------------
void
GridSeries::p_read_grid(const int& index)
{
  std::string fname(p_grid_name(index));
  
  SetOriginalGrid(p_buffer.get(), fname, p_gc, p_allow_nodata);

  std::cout << "Reading from " << fname << " ..." << std::endl;

  this->p_interp();

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
    if (p_in_time < p_max_time) {
      index = trunc(p_in_time/p_in_dt);
      p_read_grid(index);
      sendit = true;
    }
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
    this->p_copy_to_dev();
  } 
}

// -------------------------------------------------------------
// GridSeries::p_sum
// -------------------------------------------------------------
double
GridSeries::p_sum(void) const
{
  // FIXME: do this on the device

  // return (ReduceSumGrid(p_current_dev));

  std::unique_ptr<double[]> tmp(new double[p_gc.h_nx*p_gc.h_ny]);
  checkCudaErrors(cudaMemcpy2D(tmp.get(), p_gc.h_nx*sizeof(double), p_current_dev,
                               pitch, p_gc.h_nx*sizeof(double), p_gc.h_ny, DtoH));

  double result(0.0);

  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      if (tmp[j*p_gc.h_nx+i] != p_gc.nodata) {
        result += tmp[j*p_gc.h_nx+i];
      }
    }
  }

  return result;
}

// -------------------------------------------------------------
// GridSeries::p_init_dev
// -------------------------------------------------------------
void
GridSeries::p_init_dev()
{
  FillGrid(p_current_dev, (p_allow_nodata ? p_gc.nodata : 0.0));
  p_current_dev_init = true;
}

// -------------------------------------------------------------
// GridSeries::p_copy_to_dev
// -------------------------------------------------------------
void
GridSeries::p_copy_to_dev()
{
  checkCudaErrors(cudaMemcpy2D(p_current_dev, pitch, &(p_int_buffer[0]),
                               p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                               p_gc.h_ny, HtoD));
}

// -------------------------------------------------------------
//  class HyetographGridSeries
// -------------------------------------------------------------

const double HyetographGridSeries::scale(1.0/3600.0/1000.0);

// -------------------------------------------------------------
// HyetographGridSeries:: constructors / destructor
// -------------------------------------------------------------
HyetographGridSeries::HyetographGridSeries(const std::string& basename,
                                           const int& deltat,
                                           const double& tmax,
                                           const struct GridConfig& gc,
                                           double *dev_buf)
  : GridSeries(basename,
               HyetographGridSeries::scale, // converts mm/hr to m/s
               deltat, tmax, gc, dev_buf),
    p_sum_cache(0.0)
    
{}

HyetographGridSeries::~HyetographGridSeries(void)
{}

// -------------------------------------------------------------
// HyetographGridSeries::p_update
// -------------------------------------------------------------
void
HyetographGridSeries::p_update(const double& t)
{
  int index;
  
  if (p_in_time < 0.0) {
    GridSeries::p_update(t);
  }

  // rainfall between t and t + p_in_dt should be set from the map at
  // t + p_in_dt

  // Also, ince the hydrograph does not change over the input times,
  // the sum can be computed and cached
  if (t >= p_max_time) {
    GridSeries::p_update(t);
  } else if (t > (p_in_time)) {
    p_in_time += p_in_dt;
    index = trunc(p_in_time/p_in_dt);
    p_read_grid(index);
    this->p_copy_to_dev();
    p_sum_cache = GridSeries::p_sum();
  }
}

// -------------------------------------------------------------
// HyetographGridSeries::p_sum
// -------------------------------------------------------------
double
HyetographGridSeries::p_sum(void) const
{
  return p_sum_cache;
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
                                               const struct GridConfig& gc,
                                               double *dev_buf)
  : GridSeries(basename, scale, deltat, tmax, gc, dev_buf)
{
    // warning: global variables
    // Call SetDeviceConstants() first
    size_t width  = (GridDim.x * BlockDim.x) * sizeof(double);
    size_t height = (GridDim.y * BlockDim.y);
    
    checkCudaErrors(cudaMallocPitch((void**)&p_t0_dev, &pitch, width, height));
    checkCudaErrors(cudaMallocPitch((void**)&p_t1_dev, &pitch, width, height));
}

InterpolatedGridSeries::~InterpolatedGridSeries(void)
{}

// -------------------------------------------------------------
// InterpolatedGridSeries::p_update
// -------------------------------------------------------------
void
InterpolatedGridSeries::p_update(const double& t)
{
  int index;

  if (!p_current_dev_init) p_init_dev();
  
  if (p_in_time < 0.0) {
    
      p_in_time = t;
      index = trunc(p_in_time/p_in_dt);
      p_read_grid(index);

      checkCudaErrors(cudaMemcpy2D(p_t0_dev, pitch, &(p_int_buffer[0]),
                                 p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                 p_gc.h_ny, HtoD));
      
      p_in_time = index*p_in_dt;
      p_next_time = p_in_time + p_in_dt;
      p_read_grid(index + 1);

      checkCudaErrors(cudaMemcpy2D(p_t1_dev, pitch, &(p_int_buffer[0]),
                                 p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                 p_gc.h_ny, HtoD));

  } else if (t >= p_max_time) {

    if (!p_done) {
    
      p_in_time = p_next_time;
      p_next_time = p_in_time + p_in_dt;

      checkCudaErrors(cudaMemcpy2D(p_t0_dev, pitch, p_t1_dev, pitch,
                                   p_gc.h_nx*sizeof(double), p_gc.h_ny,
                                   cudaMemcpyDeviceToDevice));
      
    }
    p_done = true;
      
  } else if (t >= p_next_time) {

    checkCudaErrors(cudaMemcpy2D(p_t0_dev, pitch, p_t1_dev, pitch,
                                 p_gc.h_nx*sizeof(double), p_gc.h_ny,
                                 cudaMemcpyDeviceToDevice));
    
    p_in_time = p_next_time;
    p_next_time = p_in_time + p_in_dt;
    index = trunc(p_next_time/p_in_dt);
    p_read_grid(index);

    checkCudaErrors(cudaMemcpy2D(p_t1_dev, pitch, &(p_int_buffer[0]),
                                 p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                 p_gc.h_ny, HtoD));

  }

  double factor((t - p_in_time)/(p_next_time - p_in_time));

  InterpGrid(factor, p_t0_dev, p_t1_dev, p_current_dev);
  
}  

// -------------------------------------------------------------
// InterpolatedGridSeries::p_init_dev
// -------------------------------------------------------------
void
InterpolatedGridSeries::p_init_dev()
{
  GridSeries::p_init_dev();
  FillGrid(p_t0_dev, (p_allow_nodata ? p_gc.nodata : 0.0));
  FillGrid(p_t1_dev, (p_allow_nodata ? p_gc.nodata : 0.0));
}

// -------------------------------------------------------------
// InterpolatedGridSeries::p_copy_to_dev
// -------------------------------------------------------------
void
InterpolatedGridSeries::p_copy_to_dev(void)
{
}  

