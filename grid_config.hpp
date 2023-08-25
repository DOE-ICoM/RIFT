// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: grid_config.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 24, 2023 by Perkins
// Last Change: 2023-08-24 11:31:43 d3g096
// -------------------------------------------------------------


#ifndef _grid_config_hpp_
#define _grid_config_hpp_

struct GridConfig 
{
  int b_nx, b_ny;
  double b_xll, b_yll;
  double cellsize_original;
  double nodata;

  int h_nx, h_ny; 		// abscissa and oridnate grid dimension
  double h_xll, h_yll;
  double cellsize;
  double h_dx, h_dy;		// abscissa and ordinate resolution (m)
};

#endif
