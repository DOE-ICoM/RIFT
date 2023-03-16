#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include "constants.h"
#include "io.h"

int    b_nx,  b_ny;
double b_xll, b_yll;
double cellsize_original;
double nodata;

int h_nx, h_ny;
double h_xll, h_yll;
double cellsize;

double h_dx, h_dy;

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

void InitBathymetry(double *&b, std::string filename, const bool& square_cells) {

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
        ifs >> b_nx; line++;
	
        ifs >> dummy; 
        ifs >> b_ny; line++;

        ifs >> dummy;
        ifs >> b_xll; line++;

        ifs >> dummy;
        ifs >> b_yll; line++;
    
        ifs >> dummy;
        ifs >> cellsize_original; line++;
    
        double nodata;
        ifs >> dummy;
        ifs >> nodata; line++;

        // This array holds the bathymetry points defined by the input DEM 
        b = (double *)malloc(b_nx * b_ny * sizeof(double ));

        h_xll = b_xll + cellsize_original/2.0;
        h_yll = b_yll + cellsize_original/2.0;

        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < b_nx; i++) {
                int id = j*b_nx + i;
                ifs >> value;
                if (value == nodata) {
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
    for (int j = b_ny - 1; j >= 0; j--) {
        for (int i = 0; i < b_nx; i++) {
            int id = j*b_nx + i;
            b[id] -= minimum;
        }
    }

    double tmp_dx, tmp_dy;

    ComputeCellSize(b_nx, b_ny, h_xll, h_yll, cellsize_original, tmp_dx, tmp_dy);

    // Convert cellsize from degrees to meters and set dx, dy
    h_dx     = (double )cellsize_original*6378137.f*(double )pi / 180.f;
    h_dy     = (double )cellsize_original*6378137.f*(double )pi / 180.f;    
    cellsize =        cellsize_original*6378137.0*       pi / 180.0;

    h_dx = tmp_dx;
    h_dy = tmp_dy;

    h_nx = b_nx + 4 - 1;
    h_ny = b_ny + 4 - 1;

    std::cout << "Computed Cell Sizes: " << std::endl;
    std::cout << "\tOriginal (deg): " << cellsize_original << std::endl;
    // std::cout << "\tOriginal (m): " << cellsize << std::endl;
    std::cout << "\tInternal (m): " << h_dx << ", " << h_dy << std::endl;
}

void ReadOriginalGrid(double *&G_original, std::string filename) {
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

        if (b_ny != nrows || b_nx != ncols) {
            std::ostringstream msg;
            msg << filename << ": incorrect size: "
                << " got (" << nrows << "x" << ncols << ")"
                << " expected (" << b_ny << "x" << b_nx << ")";
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
        G_original = (double *)malloc(b_nx * b_ny * sizeof(double ));
	
        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < b_nx; i++) {
                int id = j*b_nx + i;
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

void SetOriginalGrid(double *G_original, std::string filename) {
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

        if (b_ny != nrows || b_nx != ncols) {
            std::ostringstream msg;
            msg << filename << ": incorrect size: "
                << " got (" << nrows << "x" << ncols << ")"
                << " expected (" << b_ny << "x" << b_nx << ")";
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

        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < b_nx; i++) {
                int id = j*b_nx + i;
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

void FreeBathymetry(double *&b) {
    free(b);
}

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
