// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: grid_config.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-10-04 08:36:33 d3g096
// -------------------------------------------------------------


#ifndef _grid_config_hpp_
#define _grid_config_hpp_

struct GridConfig 
{
  static const double nodata;
  int b_nx, b_ny;
  double b_xll, b_yll;
  double cellsize_original;

  int h_nx, h_ny; 		// abscissa and oridnate grid dimension
  double h_xll, h_yll;
  double cellsize;
  double h_dx, h_dy;		// abscissa and ordinate resolution (m)
};

#endif
