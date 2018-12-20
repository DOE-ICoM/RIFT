#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <string>
#include <sys/time.h>
//#include <time.h>		// for windows, Youcan
#include "source.h"

class Simulator {

public:
	std::string DEM_file;
	std::string output_file;

	// Simulation time variables
	double t0; // initial time (s) (defined by user)
	double tf; // final time (s) (defined by user)
	double t;  // instantaneous time (s)
	double dt; // timestep (s)

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
	bool  h_print;           // specifies if height data should be printed
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
	double source_X;
	double source_Y;
	int source_x;
	int source_y;
	int dambreak_source_idx;     // array index of dam breach source location
	
	/** Added for 1D-2D*/
	int *source_idx, *source_idx_dev;
	double *source_rate, *source_rate_dev;


	bool rainfall_averaged;      // specifies if system-averaged rainfall occurs
	Source hyetograph;           // hyetograph describing the rainfall (s, mm/hr)
	std::string hyetograph_file; // location of the hyetograph file

	bool rainfall_gridded;         // specifies if gridded rainfall occurs 
    std::string hyetograph_prefix; // file location/prefix of gridded rainfall data
	double hyetograph_t;            // initial gridded hyetograph time (s)
	double hyetograph_tf;           // final gridded hyetograph time (s)
	double hyetograph_dt;           // timestep between gridded rainfall data (s)
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
	void InitializeSources(long numSources);
	void SetSourceLocation(int i,double X, double Y);
	void setSourceValue (int i, double value);
	void updateSources();
	long NumSources;
	int count;

	void PrintData(void);        // prints requested simulation data
	void PrintSummaryData(void); // prints summarized simulation data

	void   StartTimer(void); // starts simulation timer
	double EndTimer(void);   // ends simulation timer
};

#endif
