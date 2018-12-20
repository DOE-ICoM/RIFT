#include <stdlib.h>
#include <fstream>
#include <omp.h>
#include <iostream>
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
    std::ifstream ifs(filename.c_str());

    std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
    ifs >> dummy;
    ifs >> b_nx;
	
    ifs >> dummy;
    ifs >> b_ny;

    ifs >> dummy;
    ifs >> b_xll;

    ifs >> dummy;
    ifs >> b_yll;
    
    ifs >> dummy;
    ifs >> cellsize_original;
    
    double nodata;
    ifs >> dummy;
    ifs >> nodata;

    // This array holds the bathymetry points defined by the input DEM 
    b = (double *)malloc(b_nx * b_ny * sizeof(double ));

	h_xll = b_xll + cellsize_original/2.0;
	h_yll = b_yll + cellsize_original/2.0;

	double minimum = 9999.f;

    double value;
    // Gridded data must be read in a "flipped" fashion
    for (int j = b_ny - 1; j >= 0; j--) {
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
    std::ifstream ifs(filename.c_str());
	
    std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
	double      junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;
    
    ifs >> dummy;
    ifs >> junk;
    
    ifs >> dummy;
    double nodata;
    ifs >> nodata;

    // This array holds the bathymetry points defined by the input DEM 
    G_original = (double *)malloc(b_nx * b_ny * sizeof(double ));
	
    double value;
    // Gridded data must be read in a "flipped" fashion
    for (int j = b_ny - 1; j >= 0; j--) {
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
}

Grid *ReadGrid(std::string grid_file) {
	struct Grid *G = (struct Grid*)malloc(sizeof(struct Grid));

    std::ifstream ifs(grid_file.c_str());

    std::string temp_string;

	ifs >> temp_string;
	ifs >> G->num_columns;

	ifs >> temp_string;
	ifs >> G->num_rows;

	ifs >> temp_string;
	ifs >> G->xll;

	ifs >> temp_string;
	ifs >> G->yll;

	ifs >> temp_string;
	ifs >> G->cellsize;

	ifs >> temp_string;
	ifs >> G->nodata;

    G->data = (double *)malloc(G->num_columns*G->num_rows*sizeof(double ));

    for (int j = G->num_rows-1; j != -1; j--) {
        for (int i = 0; i < G->num_columns; i++) {
			ifs >> G->data[j*G->num_columns+i];

            if (G->data[j*G->num_columns+i] == G->nodata) {
                printf("NODATA value(s) present in %s. Exiting\n",
				       grid_file.c_str());
                exit(1);
			}
		}
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
    std::ifstream ifs(filename.c_str());

    std::string dummy; // used to read in preceding words (e.g., ncols, nrows)
	double      junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;

    ifs >> dummy;
    ifs >> junk;
    
    ifs >> dummy;
    ifs >> junk;
    
    ifs >> dummy;
    double nodata;
    ifs >> nodata;

    double value;
    // Gridded data must be read in a "flipped" fashion
    for (int j = b_ny - 1; j >= 0; j--) {
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
}

void FreeBathymetry(double *&b) {
    free(b);
}
