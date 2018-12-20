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

void InitBathymetry(double *&b, std::string filename);
void ReadOriginalGrid(double *&G_original, std::string filename);
void SetOriginalGrid(double *G_original, std::string filename);
void FreeBathymetry(double *&b);

#endif
