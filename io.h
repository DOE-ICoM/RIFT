#ifndef IO_H
#define IO_H

#include <string>
struct GridConfig 
{
	int b_nx, b_ny;
	double b_xll, b_yll;
	double cellsize_original;
	double nodata;

	int h_nx, h_ny; 		// abscissa and oridnate grid dimension
	double h_xll, h_yll;
	double cellsize;
	double h_dx, h_dy;		// abscissa and ordinate resolution (m)
};

typedef struct Grid {
	int num_columns;
	int num_rows;
	double xll;
	double yll;
	double cellsize;
	double nodata;
	double* data;
	int rank;
	int pcount;
	Grid(const int cols, const int rows, const double x_lower_left,
		 const double y_lower_left, const double cellsize_value = 0.0, 
		 const double nodata_value= -9999.0)
		 :num_columns(cols),num_rows(rows),xll(x_lower_left), yll(y_lower_left),
		 cellsize(cellsize_value), nodata(nodata_value)
	{}
	void InitMPI(const int size, const int rank_p) { pcount = size; rank = rank_p;}
	
	void Fill(const double value)
	{
		data = (double*)malloc(sizeof(double)*num_columns*num_rows);
		for(int i = 0; i < num_rows; i++)
			for(int j=0; j<num_columns; j++)
				data[i*num_columns + j] = 0;
	}

	double* get_data() {return data;}
	int get_num_rows() {return num_rows;}
	int get_num_columns() {return num_columns;}
} Grid;

Grid *ReadGrid  (std::string filename);
Grid *CreateGrid(int num_columns, int num_rows, double xll, double yll,
                 double cellsize, double nodata, double *data);


// void InitBathymetry(double *&b, std::string filename);
void InitBathymetry(double *&b, std::string filename, GridConfig& grid_config, const bool& square_cells);
void ReadOriginalGrid(double *&G_original, std::string filename, GridConfig& grid_config);
void SetOriginalGrid(double *G_original, std::string filename, GridConfig& grid_config);
void FreeBathymetry(double *&b);

#endif
