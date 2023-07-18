#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <math.h>
#include "constants.h"
#include "io.h"

// int    b_nx,  b_ny;
// double b_xll, b_yll;
// double cellsize_original;
// double nodata;

// int h_nx, h_ny;
// double h_xll, h_yll;
// double cellsize;

// double h_dx, h_dy;

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

static void ComputeCellSize(const int &nx, const int &ny,
                     const double &xmin, const double &ymin, const double &cellsize,
                     double &dx, double &dy) {
    double domain_width_deg(nx*cellsize);
    double domain_height_deg(ny*cellsize);
    double xcenter(xmin + domain_width_deg/2.0), xmax(xmin + domain_width_deg);
    double ycenter(ymin + domain_height_deg/2.0), ymax(ymin + domain_height_deg);
    
    double domain_width_m = HaversineDistance(xmin, ycenter, xmax, ycenter);
    double domain_height_m = HaversineDistance(xcenter, ymin, xcenter, ymax);

    dx = domain_width_m/nx*1000.0;
    dy = domain_height_m/ny*1000.0;
}

// depth of water - oceans, seas or lakes
void InitBathymetry(double *&b, std::string filename, GridConfig& grid_config, const bool& square_cells) {

    double minimum = 9999.f;
    std::ifstream ifs;
    int line = 1;
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        ifs.open(filename.c_str());
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": error: cannot open"
                  << std::endl;
        throw e;
    }

    try {
        std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
        ifs >> dummy;
        ifs >> grid_config.b_nx; line++;
        // grid_config.b_nx = b_nx;
	
        ifs >> dummy; 
        ifs >> grid_config.b_ny; line++;
        // grid_config.b_ny = b_ny;

        ifs >> dummy;
        ifs >> grid_config.b_xll; line++;
        // grid_config.b_xll = b_xll;

        ifs >> dummy;
        ifs >> grid_config.b_yll; line++;
        // grid_config.b_yll = b_yll;
    
        ifs >> dummy;
        ifs >> grid_config.cellsize_original; line++;
        // grid_config.cellsize_original = cellsize_original;
    
        double nodata;
        ifs >> dummy;
        ifs >> grid_config.nodata; line++;
        // grid_config.nodata = nodata;

        // This array holds the bathymetry points defined by the input DEM 
        b = (double *)malloc(grid_config.b_nx * grid_config.b_ny * sizeof(double ));

        grid_config.h_xll = grid_config.b_xll + grid_config.cellsize_original/2.0;
        grid_config.h_yll = grid_config.b_yll + grid_config.cellsize_original/2.0;
        // grid_config.h_xll = h_xll;
        // grid_config.h_yll = h_yll;

        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = grid_config.b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < grid_config.b_nx; i++) {
                int id = j*grid_config.b_nx + i;
                ifs >> value;
                if (value == grid_config.nodata) {
                    printf("There are NODATA values present in the DEM. Exiting.\n");
                    exit(1);
                } else {
                    b[id] = value;
                    minimum = (b[id] < minimum) ? b[id] : minimum;
                }
            }
        }

        ifs.close();
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": read error at line " << line
                  << std::endl;
        throw e;
    }
        
#pragma omp parallel for
    for (int j = grid_config.b_ny - 1; j >= 0; j--) {
        for (int i = 0; i < grid_config.b_nx; i++) {
            int id = j*grid_config.b_nx + i;
            b[id] -= minimum;
        }
    }

    grid_config.cellsize = grid_config.cellsize_original*6378137.0*pi / 180.0;
    if (square_cells) {
      // Convert cellsize from degrees to meters and set dx, dy
      grid_config.h_dx = (double )grid_config.cellsize_original*6378137.f*(double )pi / 180.f;
      grid_config.h_dy = (double )grid_config.cellsize_original*6378137.f*(double )pi / 180.f;    
    } else {
      ComputeCellSize(grid_config.b_nx, grid_config.b_ny,
                      grid_config.h_xll, grid_config.h_yll,
                      grid_config.cellsize_original,
                      grid_config.h_dx, grid_config.h_dy);
    }

    grid_config.h_nx = grid_config.b_nx + 4 - 1;
    grid_config.h_ny = grid_config.b_ny + 4 - 1;

    std::cout << "Computed Cell Sizes: " << std::endl;
    std::cout << "\tOriginal (deg): " << grid_config.cellsize_original << std::endl;
    std::cout << "\tOriginal (m): " << grid_config.cellsize << std::endl;
    std::cout << "\tInternal (m): "
              << grid_config.h_dx << ", "
              << grid_config.h_dy << std::endl;

}

void ReadOriginalGrid(double *&G_original, std::string filename, GridConfig& grid_config) {
    std::ifstream ifs;
    std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
	double      junk;
    int ncols, nrows;
    int line = 1;
	
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        ifs.open(filename.c_str());
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": error: cannot open"
                  << std::endl;
        throw e;
    }

    try {

        ifs >> dummy;
        ifs >> ncols; line++;

        ifs >> dummy;
        ifs >> nrows; line++;

        if (grid_config.b_ny != nrows || grid_config.b_nx != ncols) {
            std::ostringstream msg;
            msg << filename << ": incorrect size: "
                << " got (" << nrows << "x" << ncols << ")"
                << " expected (" << grid_config.b_ny << "x" << grid_config.b_nx << ")";
            throw std::runtime_error(msg.str());
        }

        ifs >> dummy;
        ifs >> junk; line++;

        ifs >> dummy;
        ifs >> junk; line++;
    
        ifs >> dummy;
        ifs >> junk; line++;
    
        ifs >> dummy;
        double nodata;
        ifs >> nodata; line++;

        // This array holds the bathymetry points defined by the input DEM 
        G_original = (double *)malloc(grid_config.b_nx * grid_config.b_ny * sizeof(double ));
	
        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = grid_config.b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < grid_config.b_nx; i++) {
                int id = j*grid_config.b_nx + i;
                ifs >> value;
                if (value == nodata) {
                    printf("There are NODATA values present in the DEM. Exiting.\n");
                    exit(1);
                } else {
                    G_original[id] = value;
				
                }
            }
        }
        ifs.close();
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": read error at line " << line
                  << std::endl;
        throw e;
    }
    
}

Grid *ReadGrid(std::string grid_file) {
	struct Grid *G = (struct Grid*)malloc(sizeof(struct Grid));

    std::ifstream ifs;
    int line = 1;
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        ifs.open(grid_file.c_str());
    } catch (std::ifstream::failure& e) {
        std::cerr << grid_file << ": error: cannot open"
                  << std::endl;
        throw e;
    }
    
    try {
        std::string temp_string;

        ifs >> temp_string;
        ifs >> G->num_columns; line++;

        ifs >> temp_string;
        ifs >> G->num_rows; line++;

        ifs >> temp_string;
        ifs >> G->xll; line++;

        ifs >> temp_string;
        ifs >> G->yll; line++;

        ifs >> temp_string;
        ifs >> G->cellsize; line++;

        ifs >> temp_string;
        ifs >> G->nodata; line++;

        G->data = (double *)malloc(G->num_columns*G->num_rows*sizeof(double ));

        for (int j = G->num_rows-1; j != -1; j--, line++) {
            for (int i = 0; i < G->num_columns; i++) {
                ifs >> G->data[j*G->num_columns+i];

                if (G->data[j*G->num_columns+i] == G->nodata) {
                    printf("NODATA value(s) present in %s. Exiting\n",
                           grid_file.c_str());
                    exit(1);
                }
            }
        }
        ifs.close();
    } catch (std::ifstream::failure& e) {
        std::cerr << grid_file << ": read error at line " << line
                  << std::endl;
        throw e;
    }

	return G;
}

Grid *CreateGrid(int num_columns, int num_rows, double xll, double yll,
                 double cellsize, double nodata, double *data) {
	struct Grid *G = (struct Grid*)malloc(sizeof(struct Grid));

	G->num_columns = num_columns;
	G->num_rows    = num_rows;
	G->xll         = xll;
	G->yll         = yll;
	G->cellsize    = cellsize;
	G->nodata      = nodata;

    G->data        = (double *)malloc(G->num_columns*G->num_rows*sizeof(double ));

	return G;
}

void SetOriginalGrid(double *G_original, std::string filename, GridConfig& grid_config) {
    std::ifstream ifs;
	
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        ifs.open(filename.c_str());
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": error: cannot open"
                  << std::endl;
        throw e;
    }
    
    std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
	double      junk;
    int nrows, ncols, line = 1;

    try {
        ifs >> dummy;
        ifs >> ncols; line++;

        ifs >> dummy;
        ifs >> nrows; line++;

        if (grid_config.b_ny != nrows || grid_config.b_nx != ncols) {
            std::ostringstream msg;
            msg << filename << ": incorrect size: "
                << " got (" << nrows << "x" << ncols << ")"
                << " expected (" << grid_config.b_ny << "x" << grid_config.b_nx << ")";
            throw std::runtime_error(msg.str());
        }

        ifs >> dummy;
        ifs >> junk; line++;

        ifs >> dummy;
        ifs >> junk; line++;
    
        ifs >> dummy;
        ifs >> junk; line++;
    
        ifs >> dummy;
        double nodata;
        ifs >> nodata; line++;
        grid_config.nodata = nodata;

        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = grid_config.b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < grid_config.b_nx; i++) {
                int id = j*grid_config.b_nx + i;
                ifs >> value;
                if (value == grid_config.nodata) {
                    printf("There are NODATA values present in the DEM. Exiting.\n");
                    exit(1);
                } else {
                    G_original[id] = value;
                }
            }
        }
        ifs.close();
    } catch (std::ifstream::failure& e) {
        std::cerr << filename << ": read error at line " << line
                  << std::endl;
        throw e;
    }
}

static void writeHeader(std::ofstream &thefile, const int& h_nx, const int& h_ny) {
    thefile.precision(dbl::digits10);
	thefile << "ncols         " << h_nx - 4          << std::endl;
	thefile << "nrows         " << h_ny - 4          << std::endl;
	thefile << "xllcorner     " << h_xll             << std::endl;
	thefile << "yllcorner     " << h_yll             << std::endl;
	thefile << "cellsize      " << cellsize_original << std::endl;
	thefile << "NODATA_value  " << -9999             << std::endl;
    thefile.precision(flt::digits10);
}

void writeGrid(const std::string& fname, double *x, const int& h_nx, const int& h_ny)
{
    std::ofstream out;
    out.open(fname.c_str());
    writeHeader(out, h_nx, h_ny);
    if (out.is_open()) {
        for (int j = h_ny - 3; j >= 2; j--) {
            for (int i = 2; i < h_nx-2; i++) {
                int   id = j*h_nx + i;
                out << x[id] << " ";
            }
            out << std::endl;
        }
    } else {
        std::string msg(fname);
        msg += std::string(": error: cannot open for writing");
        throw std::runtime_error(msg);
    }
    out.close();
}    


void FreeBathymetry(double *&b) {
    free(b);
}

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
