// -------------------------------------------------------------
// file: GridSeriesTest.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 25, 2023 by Perkins
// Last Change: 2023-08-25 08:02:53 d3g096
// -------------------------------------------------------------

#include <iostream>

#include "constants.h"
#include "io.h"
#include "GridSeries.cuh"


// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  std::string basename;
  std::string g0name;

  GridConfig gc;
  double *b;

  InitBathymetry(b, g0name, gc, true);

  GridSeries g(basename, 1.0, 3600.0, gc);

  return 0;
}
