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

typedef std::numeric_limits<double> dbl;
typedef std::numeric_limits<double>  flt;

// int    b_nx,  b_ny;
// double b_xll, b_yll;
// double cellsize_original;
// double nodata;

// int h_nx, h_ny;
// double h_xll, h_yll;
// double cellsize;

// double h_dx, h_dy;

// depth of water - oceans, seas or lakes
void InitBathymetry(double *&b, std::string filename, GridConfig& grid_config) {

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
        ifs >> nodata; line++;

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
    for (int j = grid_config.b_ny - 1; j >= 0; j--) {
        for (int i = 0; i < grid_config.b_nx; i++) {
            int id = j*grid_config.b_nx + i;
            b[id] -= minimum;
        }
    }

    grid_config.ComputeCellSize();
    grid_config.ReportCellSize(std::cout);

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

void SetOriginalGrid(double *G_original, std::string filename,
                     GridConfig& grid_config, const bool& nodata_ok) {
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

        double value;
        // Gridded data must be read in a "flipped" fashion
        for (int j = grid_config.b_ny - 1; j >= 0; j--, line++) {
            for (int i = 0; i < grid_config.b_nx; i++) {
                int id = j*grid_config.b_nx + i;
                ifs >> value;
                if (value == nodata) {
                    if (nodata_ok) {
                        G_original[id] = grid_config.nodata;
                    } else {
                        printf("There are NODATA values present in the DEM. Exiting.\n");
                        exit(1);
                    }
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

static void writeHeader(std::ofstream &thefile, const GridConfig& gc) {
    thefile.precision(dbl::digits10);
	thefile << "ncols         " << gc.h_nx - 4          << std::endl;
	thefile << "nrows         " << gc.h_ny - 4          << std::endl;
	thefile << "xllcorner     " << gc.h_xll             << std::endl;
	thefile << "yllcorner     " << gc.h_yll             << std::endl;
	thefile << "cellsize      " << gc.cellsize_original << std::endl;
	thefile << "NODATA_value  " << gc.nodata            << std::endl;
    thefile.precision(flt::digits10);
}

void writeGrid(const std::string& fname, double *x, const GridConfig& gc)
{
    std::ofstream out;
    out.open(fname.c_str());
    writeHeader(out, gc);
    if (out.is_open()) {
        for (int j = gc.h_ny - 3; j >= 2; j--) {
            for (int i = 2; i < gc.h_nx-2; i++) {
                int   id = j*gc.h_nx + i;
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
