// -------------------------------------------------------------
// file: GridSeriesTest.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 25, 2023 by Perkins
// Last Change: 2023-10-09 11:46:36 d3g096
// -------------------------------------------------------------

#include <iostream>

#include "constants.h"
#include "io.h"
#include "grid.h"
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

  cudaSetDevice(0);

  if (argc < 3) {
  }

  g0name = argv[1];
  basename = argv[2];
  
  InitBathymetry(b, g0name, gc, true);

  SetDeviceConstants(gc.h_nx, gc.h_ny, gc.h_dx, gc.h_dy, 1.0);

  std::unique_ptr<GridSeries>
    g(new 
#ifdef INTERPOLATE
      InterpolatedGridSeries
#else
      GridSeries
#endif
      (basename, 1.0, 3600.0, 43200.0, gc));
  g->allow_no_data(true);
  double tmax = 86400;
  double dt = 600;

  for (double t = 0; t <= tmax; t += dt) {
    g->update(t);
    std::cout << "Time = " << t << ", "
              << "Sum = " << g->sum()
              << std::endl;
  }    

  return 0;
}
