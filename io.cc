#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <stdexcept>
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

void InitBathymetry(double *&b, std::string filename) {

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

    // Convert cellsize from degrees to meters and set dx, dy
    h_dx     = (double )cellsize_original*6378137.f*(double )pi / 180.f;
    h_dy     = (double )cellsize_original*6378137.f*(double )pi / 180.f;    
    cellsize =        cellsize_original*6378137.0*       pi / 180.0;

    h_nx = b_nx + 4 - 1;
    h_ny = b_ny + 4 - 1;
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
