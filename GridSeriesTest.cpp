// -------------------------------------------------------------
// file: GridSeriesTest.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 25, 2023 by Perkins
// -------------------------------------------------------------

#include <iostream>
#include <libgen.h>
#include <unistd.h>

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
  
  std::string progname(basename(argv[0]));
  std::string usage(progname);
  usage += " [-i] [-H] [-T tmax] [-d i_dt] [-s o_dt] dem basename";

  bool dohyetograph(false);
  bool dointerp(false);

  int opt;
  double instep(600), outstep(instep), tmax(86400);

  while ((opt = getopt(argc, argv, "iHT:d:s:")) != -1) {
    switch (opt) {
    case 'H':
      dohyetograph = true;
      break;
    case 'i':
      dointerp = true;
      break;
    case 'T':
      tmax = atof(optarg);
      break;
    case 'd':
      instep = atof(optarg);
      break;
    case 's':
      outstep = atof(optarg);
      break;
    default: /* '?' */
      std::cerr << "error: unknown option '" << opt << "'" << std::endl;
      std::cerr << "Usage: " << usage << std::endl;
      exit(3);
    }
  }

  if (argc - optind < 2) {
    std::cerr << "Usage: " << usage << std::endl;
    exit(3);
  }

  std::string g0name(argv[optind]);
  std::string basename(argv[optind+1]);
  GridConfig gc;
  double *b;

  cudaSetDevice(0);

  
  InitBathymetry(b, g0name, gc);

  SetDeviceConstants(gc.h_nx, gc.h_ny, gc.h_dx, gc.h_dy, 1.0, gc.nodata);

  std::unique_ptr<GridSeries> g;
  std::string msg("Exercising class "); 
  if (dohyetograph) {
    if (dointerp) {
      msg += "Interpolated";
      g.reset(new InterpolatedHyetographGridSeries(basename, instep, tmax, gc));
    } else {
      g.reset(new HyetographGridSeries(basename, instep, tmax, gc));
    }
    msg += "HyetographGridSeries";
  } else {
    if (dointerp) {
      msg += "Interpolated";
      g.reset(new InterpolatedGridSeries(basename, 1.0, instep, tmax, gc));
    } else {
      g.reset(new GridSeries(basename, 1.0, instep, tmax, gc));
    }
    msg += "GridSeries";
  }
  std::cout << msg << std::endl;
  
  g->allow_no_data(true);

  for (double t = 0; t <= tmax; t += outstep) {
    g->update(t);
    std::cout << "Time = " << t << ", "
              << "Sum = " << g->sum()
              << std::endl;
  }    

  return 0;
}
