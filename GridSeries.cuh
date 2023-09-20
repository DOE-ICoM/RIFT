// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-09-20 12:24:50 d3g096
// -------------------------------------------------------------


#ifndef _GridSeries_hpp_
#define _GridSeries_hpp_

#include <string>
#include <memory>
#include "grid_config.hpp"

// -------------------------------------------------------------
// class GridSeries
//
// A thing to read and manipulate a series of grids that are placed on
// the GPU.
// -------------------------------------------------------------
class GridSeries {
public:

  /// Default constructor.
  GridSeries(const std::string& basename,
             const double& scale,
             const int& deltat,
             const double& tmax,
             const struct GridConfig& gc,
             double *dev_buf = NULL);

  /// Destructor
  ~GridSeries(void);

  /// Update grid series to specified time
  void update(const double& t)
  {
    p_update(t);
  }

  /// Sum the entire current grid
  double sum(void) const
  {
    return(p_sum());
  }

  /// Get the on-device buffer
  double *grid_dev(void) const
  {
    return p_current_dev;
  }

protected:

  /// Protected copy constructor to avoid unwanted copies.
  GridSeries(const GridSeries& old);

  /// Construct a grid file name
  std::string p_grid_name(const int& index) const;

  /// Read the next new grid and make ready to use
  virtual void p_read_grid(void);

  /// Update grid series to specified time (specialized)
  virtual void p_update(const double& t);

  /// Sum the entire current grid
  virtual double p_sum(void) const;
  
  /// Local copy
  struct GridConfig p_gc;

  /// Base name for input grid files
  std::string p_basename;

  /// Scale to apply to read series (primarily for units)
  double p_scale;

  /// Input buffer for input grid
  std::unique_ptr<double[]> p_buffer;

  /// Input buffer for spatially interpolated grid
  std::unique_ptr<double[]> p_int_buffer;

  /// Current input time and step, s
  double p_in_time, p_in_dt;

  /// Maximum input time, s
  double p_max_time;

  /// Pointer to on-device grids
  double *p_current_dev;

  /// Was p_current_dev created externally?
  bool p_external;

  /// This is set when time is past maximum
  bool p_done;
};




#endif
