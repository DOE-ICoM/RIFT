#ifndef IO_H
#define IO_H

#include <string>

typedef struct Grid {
	int num_columns;
	int num_rows;
	double xll;
	double yll;
	double cellsize;
	double nodata;
	double* data;
} Grid;

Grid *ReadGrid  (std::string filename);
Grid *CreateGrid(int num_columns, int num_rows, double xll, double yll,
                 double cellsize, double nodata, double *data);

void InitBathymetry(double *&b, std::string filename, const bool& square_cells);
void ReadOriginalGrid(double *&G_original, std::string filename);
void SetOriginalGrid(double *G_original, std::string filename);
void writeGrid(const std::string& fname, double *x, const int& h_nx, const int& h_ny);
void FreeBathymetry(double *&b);

#endif
