// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: GridSeries.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
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
  virtual ~GridSeries(void);

  /// Update grid series to specified time
  void update(const double& t)
  {
    this->p_update(t);
  }

  /// Sum the entire current grid
  double sum(const bool use_dev = true) const
  {
    return this->p_sum(use_dev);
  }

  /// Get the on-device buffer
  double *grid_dev(void) const
  {
    return p_current_dev;
  }

  /// Recognize no data values in input
  void allow_no_data(const bool& flag)
  {
    p_allow_nodata = flag;
    p_current_dev_init = false; // changes default value
  }

protected:

  /// Protected copy constructor to avoid unwanted copies.
  GridSeries(const GridSeries& old);

  /// Construct a grid file name
  std::string p_grid_name(const int& index) const;

  /// Method used to interpolate input to internal grid
  virtual void p_interp(void); 

  /// Read the next new grid and make ready to use
  virtual void p_read_grid(const int& index);

  /// Update grid series to specified time (specialized)
  virtual void p_update(const double& t);

  /// Sum the entire current grid (on device)
  virtual double p_sum_dev(void) const;

  /// Sum the entire current grid (on host)
  virtual double p_sum_host(void) const;

  /// Sum the entire current grid
  virtual double p_sum(const bool use_dev = true) const
  {
    return (use_dev ? p_sum_dev() : p_sum_host());
  }

  /// Initialize the device buffer(s) with default value
  virtual void p_init_dev(void);

  /// Send the map to the device
  virtual void p_copy_to_dev(void);
  
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

  /// Can the input grid(s) have cells with no data?
  bool p_allow_nodata;
  
  /// An on-device place to put a summation (or other computed value)
  mutable double *p_result_dev;

  /// A flag to indicate if the on-device buffer has been initialized
  bool p_current_dev_init;

};

// -------------------------------------------------------------
// class HyetographGridSeries
//
// A GridSeries specialized for rainfall maps. The rainfall map is
// constant through the time between input maps. The rate map used is
// that from the *next* map input time.  So, map input is handled
// sightly different than GridSeries.  Also, the scale (converts mm/hr
// to m/s) should not be set externally.
// -------------------------------------------------------------
class HyetographGridSeries
  : public GridSeries
{
public:

  /// Scale used for rainfall (converts mm/hr to m/s)
  static const double scale;

  /// Default constructor.
  HyetographGridSeries(const std::string& basename,
                       const int& deltat,
                       const double& tmax,
                       const struct GridConfig& gc,
                       double *dev_buf = NULL);

  /// Destructor
  ~HyetographGridSeries(void);

protected:

  /// Protected copy constructor to avoid unwanted copies.
  HyetographGridSeries(const HyetographGridSeries& old);

  /// Update grid series to specified time (specialized)
  void p_update(const double& t);

  /// Sum the entire current grid
  double p_sum(const bool use_dev) const;

  /// Cache of current map's sum
  double p_sum_cache;

};

// -------------------------------------------------------------
//  class InterpolatedGridSeries
// -------------------------------------------------------------
class InterpolatedGridSeries
  : public GridSeries
{
public:

  /// Default constructor.
  InterpolatedGridSeries(const std::string& basename,
                         const double& scale,
                         const int& deltat,
                         const double& tmax,
                         const struct GridConfig& gc,
                         double *dev_buf = NULL);

  /// Destructor
  ~InterpolatedGridSeries(void);
  
protected:

  /// Protected copy constructor to avoid unwanted copies.
  InterpolatedGridSeries(const InterpolatedGridSeries& old);

  /// Update grid series to specified time (specialized)
  void p_update(const double& t);

  /// Send the map to the device (specialized)
  void p_copy_to_dev(void);
  
  /// Initialize the device buffer(s) with default value
  void p_init_dev(void);

  /// Input time for second grid
  double p_next_time;

  /// Device time plane buffers
  double *p_t0_dev, *p_t1_dev;
};

// -------------------------------------------------------------
// class InterpolatedHyetographGridSeries
//
// An optional GridSeries for rainfall maps that interpolates between
// time steps.
// -------------------------------------------------------------
class InterpolatedHyetographGridSeries
  : public InterpolatedGridSeries
{
protected:

  /// Protected copy constructor to avoid unwanted copies.
  InterpolatedHyetographGridSeries(const InterpolatedHyetographGridSeries& old);

public:

  /// Default constructor.
  InterpolatedHyetographGridSeries(const std::string& basename,
                                   const int& deltat,
                                   const double& tmax,
                                   const struct GridConfig& gc,
                                   double *dev_buf = NULL)
    : InterpolatedGridSeries(basename,
                             HyetographGridSeries::scale,
                             deltat, tmax, gc, dev_buf)
  {}

  /// Destructor
  ~InterpolatedHyetographGridSeries(void)
  {}
};




#endif
