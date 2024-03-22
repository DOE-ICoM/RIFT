// -------------------------------------------------------------
// file: grid_config.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October  4, 2023 by Perkins
// -------------------------------------------------------------

#include <math.h>
#include "constants.h"
#include "grid_config.hpp"

/// A global value for internal use
const double GridConfig::nodata(-9999.0);

/// Do not change dx/dy by latitude, assume zero latitude
bool GridConfig::square_cells(true);


/* 
   Blatently plagerized from 
   https://rosettacode.org/wiki/Haversine_formula#C++
*/

inline double DegreeToRadian(double angle)
{
  return M_PI * angle / 180.0;
}

static double HaversineDistance(const double& long1, const double& lat1,
                                const double& long2, const double& lat2)
{
  const static double EarthRadiusKm = 6372.8;
  double latRad1 = DegreeToRadian(lat1);
  double latRad2 = DegreeToRadian(lat2);
  double lonRad1 = DegreeToRadian(long1);
  double lonRad2 = DegreeToRadian(long2);
    
  double diffLa = latRad2 - latRad1;
  double doffLo = lonRad2 - lonRad1;
    
  double computation =
    asin(sqrt(sin(diffLa / 2) * sin(diffLa / 2) +
              cos(latRad1) * cos(latRad2) * sin(doffLo / 2) * sin(doffLo / 2)));
  return 2 * EarthRadiusKm * computation;
}

static void ComputeNonSquare(const int &nx, const int &ny,
                             const double &xmin, const double &ymin,
                             const double &cellsize, double &dx, double &dy) {
  double domain_width_deg(nx*cellsize);
  double domain_height_deg(ny*cellsize);
  double xcenter(xmin + domain_width_deg/2.0), xmax(xmin + domain_width_deg);
  double ycenter(ymin + domain_height_deg/2.0), ymax(ymin + domain_height_deg);
  
  double domain_width_m = HaversineDistance(xmin, ycenter, xmax, ycenter);
  double domain_height_m = HaversineDistance(xcenter, ymin, xcenter, ymax);
  
  dx = domain_width_m/nx*1000.0;
  dy = domain_height_m/ny*1000.0;
}

/// Compute computational cell sizes based on input options
void GridConfig::ComputeCellSize(void)
{
  cellsize = cellsize_original*6378137.0*pi / 180.0;
  if (square_cells) {
    // Convert cellsize from degrees to meters and set dx, dy
    h_dx = (double )cellsize_original*6378137.f*(double )pi / 180.f;
    h_dy = (double )cellsize_original*6378137.f*(double )pi / 180.f;    
  } else {
    ComputeNonSquare(b_nx, b_ny, h_xll, h_yll, cellsize_original, h_dx, h_dy);
  }
  
  h_nx = b_nx + 4 - 1;
  h_ny = b_ny + 4 - 1;
}

/// Make a brief report about cell sizes
void GridConfig::ReportCellSize(std::ostream& out) const
{
  out << "Computed Cell Sizes: " << std::endl;
  out << "\tOriginal (deg): " << cellsize_original << std::endl;
  out << "\tOriginal (m): " << cellsize << std::endl;
  out << "\tInternal (m): "
      << h_dx << ", "
      << h_dy << std::endl;
}
