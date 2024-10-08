#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <string>
#include <complex>
#include <sys/time.h>
//#include <time.h>		// for windows, Youcan
#include "source.h"
#include "grid.h"
#include "io.h"
#include "GridSeries.cuh"

class Simulator {

public:
	// Constructor (initialize malloc flags)
	Simulator():f_h_o(0), f_n_o(0), f_K_o(0), f_hyetograph_o(0), f_b(0), f_h_BX(0), f_h_BY(0), 
		f_h_hyetograph(0), f_h_n(0), f_h_F(0), f_h_K(0), f_h_q(0), f_h_h(0), f_h_h_max(0),
		f_h_q_max(0), f_h_t_wet(0), f_source_rate(0), f_source_idx(0) 
		{}

	std::string DEM_file;
	std::string output_file;

	// Simulation time variables
	double t0; // initial time (s) (defined by user)
	double tf; // final time (s) (defined by user)
	double t;  // instantaneous time (s)
	double dt; // timestep (s)
        double courant_max; // Courant limit

	// Integration settings
	bool euler_integration; // specifies if Euler integration is used

	// Manning roughness coefficient variables
	bool  n_gridded;    // specifies if Manning roughness coefficients are gridded
	double n_const;      // averaged Manning's roughness coefficient
	std::string n_file; // file containing gridded Manning roughness coefficients

	// Infiltration variables
	bool infiltration;  // specifies if infiltration occurs
	std::string K_file; // file containing gridded hydraulic conductivities (mm/hr)
	double psi;          // soil suction head (m)
	double dtheta;       // soil moisture deficit (unitless)

	// Nonzero height initialization variables
	bool h_init;
	std::string h_file;

	// Desingularization constant
	double kappa;

	// Simulation output variables
    bool  b_print;           // save internal (interpolated) bathymetry map (for debugging)
    bool  n_print;           // save internal (interpolated) Mannings map (for debugging)
    
	bool  h_print;           // specifies if height data should be printed
    bool ws_print;           // specifies if water surface elevation data should be printed
	bool  q_print;           // specifies if discharge data should be printed
	bool  save_max;          // specifies if maximal data should be printed
	bool  save_arrival_time; // specifies if arrival time data should be printed
	bool steady_state;
	double t_print;           // time used to track output file creation (s)
	double dt_print;          // timestep between output files (s) (defined by user)
	int   count_print;       // number of output files printed

	bool  check_volume; // specifies if volume should be tracked (to verify that volume is conserved)
	double V_added;      // total volume added to the system (m^3)
	double V_computed;   // total volume computed in the system (m^3)

	bool dambreak;               // specifies if a dam breach occurs
	Source hydrograph;           // hydrograph describing the dam breach (s, m^3)
	std::string hydrograph_file; // location of the hydrograph file
    std::complex<double> source_X, source_Y;
	
	/** Added for 1D-2D*/
	int *source_idx, *source_idx_dev;
	double *source_rate, *source_rate_dev;


	bool rainfall_averaged;      // specifies if system-averaged rainfall occurs
	Source hyetograph;           // hyetograph describing the rainfall (s, mm/hr)
	std::string hyetograph_file; // location of the hyetograph file

    bool drain_averaged;        // specifies if system-averaged
                                // drain (reverse of
                                // rainfall_averaged) occurs
    Source drain;
    std::string drain_file;

	bool rainfall_gridded;         // specifies if gridded rainfall occurs
  
    std::unique_ptr<GridSeries> hyetograph_series;
    std::string hyetograph_prefix;  // file location/prefix of gridded rainfall data
    bool hyetograph_interp;         // linearly interpolate hyetograph in time
	double hyetograph_t;            // initial gridded hyetograph time (s)
	double hyetograph_tf;           // final gridded hyetograph time (s)
	double hyetograph_dt;           // timestep between gridded rainfall data (s)

    bool drain_gridded;
    std::unique_ptr<GridSeries> drain_series;
    std::string drain_prefix;
    double drain_t;
    double drain_tf;
    double drain_dt;
    
  bool surge_gridded; 
  std::unique_ptr<GridSeries> surge_series;
  std::string surge_prefix; // file location/prefix of surge data
	double surge_t;            // initial surge time (s)
	double surge_tf;           // final sugre time (s)
	double surge_dt;           // timestep between surgedata (s)

    bool tile_acceleration;     // accelerate computations by only considering wet tiles/blocks

  double flow_rate;
	double volume_old;
	struct timeval start_time;     // start time of simulation
	struct timeval end_time;       // end time of simulation
	//long start_time;			   // for windows only, start time of simulation	//Youcan
	//long end_time;			   // for windows only, end time of simulation	//Youcan

	int device; // specifies NVIDIA device number

	void ReadUserParams(std::string config_file); // reads user-defined parameters
	void InitSimulation(void);                    // initializes simulation arrays and variables
	void ComputeTimestep(void);                   // computes simulation timestep
	void UpdateSource(void);                      // updates water source variables
	double RunSimulation();  // runs the simulation
	void OpenSimulation(std::string config_file, double flowrate);
	void CloseSimulation(void);
	void InitializeSources(void);
	void setSourceValue (double value);
	void updateSources();
	long NumSources;
	int count;

	void PrintData(void);        // prints requested simulation data
	void PrintSummaryData(void); // prints summarized simulation data

	void   StartTimer(void); // starts simulation timer
	double EndTimer(void);   // ends simulation timer

	// Variables (moving here - previous global variables)
	// Define host arrays
	double *h_o, *n_o, *K_o, *hyetograph_o;
	double *b, *h_BC, *h_BX, *h_BY, *h_hyetograph, *h_n, *h_F, *h_K, *h_q, *h_h, *h_h_max,
      *h_q_max, *h_t_wet;
	// malloc flags
	int f_h_o, f_n_o, f_K_o, f_hyetograph_o, f_b, f_h_BC, f_h_BX, f_h_BY, 
		f_h_hyetograph, f_h_n, f_h_F, f_h_K, f_h_q, f_h_h, f_h_h_max,
		f_h_q_max, f_h_t_wet, f_source_idx, f_source_rate;

	Grid *B;

	// Define device arrays
  double *dev_w, *dev_hu, *dev_hv, *dev_w_old, *dev_hu_old, *dev_hv_old, *dev_dw, *dev_dhu, *dev_dhv, *dev_mx, *dev_my, *dev_BC, *dev_BX,
      *dev_BY, *dev_n, *dev_hyetograph_gridded_rate, *dev_drain_gridded_rate, *dev_surge_gridded_depth, *dev_F, *dev_F_old, *dev_dF, *dev_K, *dev_h, *dev_q, *dev_h_max,
      *dev_q_max, *dev_t_wet, *dev_G, *dev_time_peak, *dev_time_dry;
	bool *dev_wet_blocks;
	int *dev_active_blocks;
	int dev_wet_count;

	
	// Grid *B; // grid in device ..
	GridConfig grid_config;
};

#endif

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
