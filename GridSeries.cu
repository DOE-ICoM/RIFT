// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.cu
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-10-12 14:54:04 d3g096
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
    p_done(false),
    p_allow_nodata(false),
    p_result_dev(NULL)
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
  if (p_result_dev != NULL)
    cudaFree(p_result_dev);
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
        // std::cout << nnd << ": "
        //           << index << ", "
        //           << i << ", "
        //           << j << ": "
        //           << p_int_buffer[index]
        //           << std::endl;
      } else {
        p_int_buffer[index] = p_gc.nodata;
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
  
  // if (p_result_dev == NULL) {
  //   checkCudaErrors(cudaMalloc( (void**)&p_result_dev, sizeof(double)));
  // }

  // sumReduce <<< GridDim, BlockDim >>> (pitch, p_current_dev, p_result_dev);

  // cudaMemcpy(&result, p_result_dev, sizeof(double), DtoH);

  // std::cout << "device result = " << result << std::endl;

  // FIXME: do this on the device, somehow
  
  std::unique_ptr<double[]> tmp(new double[p_gc.h_nx*p_gc.h_ny]);
  checkCudaErrors(cudaMemcpy2D(tmp.get(), p_gc.h_nx*sizeof(double), p_current_dev,
                               pitch, p_gc.h_nx*sizeof(double), p_gc.h_ny, DtoH));

  result = 0.0;
  for (int j = 2; j < p_gc.h_ny - 2; j++) {
    for (int i = 2; i < p_gc.h_nx - 2; i++) {
      if (tmp[j*p_gc.h_nx+i] != p_gc.nodata) {
        result += tmp[j*p_gc.h_nx+i];
      }
      // std::cout << i << ", " << j << ", " << tmp[j*p_gc.h_nx+i] << std::endl;
    }
  }

  // std::cout << "host result = " << result << std::endl;
  
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
                                               const struct GridConfig& gc,
                                               double *dev_buf)
  : GridSeries(basename, scale, deltat, tmax, gc, dev_buf),
    p_t0_buffer(new double[gc.h_nx*gc.h_ny]()),
    p_t1_buffer(new double[gc.h_nx*gc.h_ny]())
{}

InterpolatedGridSeries::~InterpolatedGridSeries(void)
{}

// -------------------------------------------------------------
// InterpolatedGridSeries::p_update
// -------------------------------------------------------------
extern __global__ void gridInterp(const double& factor, const size_t& pitch,
                                  double *x0_dev, double *x1_dev, double *x_dev);

void
InterpolatedGridSeries::p_update(const double& t)
{
  int index;

  // FIXME: Do interpolation on device

  
  
  if (p_in_time < 0.0) {
    
      p_in_time = t;
      index = trunc(p_in_time/p_in_dt);
      p_read_grid(index);
      std::copy(&(p_int_buffer[0]),
                &(p_int_buffer[0]) + p_gc.h_nx*p_gc.h_ny,
                &(p_t0_buffer[0]));

      p_in_time = index*p_in_dt;
      p_next_time = p_in_time + p_in_dt;
      p_read_grid(index + 1);
      std::copy(&(p_int_buffer[0]),
                &(p_int_buffer[0]) + p_gc.h_nx*p_gc.h_ny,
                &(p_t1_buffer[0]));
      
  } else if (t >= p_max_time) {
    
    p_in_time = p_next_time;
    p_next_time = p_in_time + p_in_dt;
    std::copy(&(p_t1_buffer[0]),
              &(p_t1_buffer[0]) + p_gc.h_nx*p_gc.h_ny,
              &(p_t0_buffer[0]));
    
  } else if (t >= p_next_time) {

    std::copy(&(p_t1_buffer[0]),
              &(p_t1_buffer[0]) + p_gc.h_nx*p_gc.h_ny,
              &(p_t0_buffer[0]));
    
    p_in_time = p_next_time;
    p_next_time = p_in_time + p_in_dt;
    index = trunc(p_next_time/p_in_dt);
    p_read_grid(index);
    std::copy(&(p_int_buffer[0]),
              &(p_int_buffer[0]) + p_gc.h_nx*p_gc.h_ny,
              &(p_t1_buffer[0]));
  }

  std::uninitialized_fill(p_int_buffer.get(),
                          p_int_buffer.get() + p_gc.h_nx*p_gc.h_ny,
                          (p_allow_nodata ? p_gc.nodata : 0.0));

  double factor((t - p_in_time)/(p_next_time - p_in_time));
  for (int j = 0; j < p_gc.h_ny; j++) {
    for (int i = 0; i < p_gc.h_nx; i++) {
      int index(j*p_gc.h_nx+i);
      if ((p_t1_buffer[index] != p_gc.nodata) &&
          (p_t0_buffer[index] != p_gc.nodata)) {
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
  }

  cudaError_t ierr(cudaMemcpy2D(p_current_dev, pitch, &(p_int_buffer[0]),
                                p_gc.h_nx*sizeof(double), p_gc.h_nx*sizeof(double),
                                p_gc.h_ny, HtoD));

  checkCudaErrors(ierr);
}  
