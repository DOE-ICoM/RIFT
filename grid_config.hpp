// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: grid_config.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// -------------------------------------------------------------


#ifndef _grid_config_hpp_
#define _grid_config_hpp_

#include <ostream>

struct GridConfig 
{
  static const double nodata;
  static bool square_cells;
  static bool projected;
  int b_nx, b_ny;
  double b_xll, b_yll;
  double cellsize_original;

  int h_nx, h_ny; 		// abscissa and oridnate grid dimension
  double h_xll, h_yll;
  double cellsize;
  double h_dx, h_dy;		// abscissa and ordinate resolution (m)

  void ComputeCellSize(void);

  void ReportCellSize(std::ostream& out) const;
};

#endif
