#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include "simulator.h"
#include "constants.h"
#include "config.h"
#include "grid.h"
#include "io.h"

typedef std::numeric_limits<double> dbl;
typedef std::numeric_limits<double>  flt;


// Define device arrays
double *w, *hu, *hv, *w_old, *hu_old, *hv_old, *dw, *dhu, *dhv, *mx, *BC, *BX,
      *BY, *n, *hyetograph_gridded_rate, *F, *F_old, *dF, *K, *h, *q, *h_max,
	*q_max, *t_wet, *time_peak, *time_dry, *G;									//added time_peak and time_dry by Youcan on 20170908
bool *wet_blocks;
int *active_blocks;
int wet_count;

Grid *B;

// Define host arrays
double *h_o, *n_o, *K_o, *hyetograph_o;
double *b, *h_BX, *h_BY, *h_hyetograph, *h_n, *h_F, *h_K, *h_q, *h_h, *h_h_max,
      *h_q_max, *h_t_wet;

void Simulator::ReadUserParams(std::string config_file) {
	ConfigFile cfg(config_file.c_str()); // save the user-defined parameter file
	                                     // into the ConfigFile instance cfg

	// Save the file location of the bathymetry
	DEM_file = cfg.getValueOfKey<std::string>("DEM");

	t0 = cfg.getValueOfKey<double>("t0");
	tf = cfg.getValueOfKey<double>("tf");
	t  = cfg.getValueOfKey<double>("t0");
	dt = 0.0000001f;

	if (cfg.keyExists("euler_integration")) {
		euler_integration = cfg.getValueOfKey<bool>("euler_integration");
	} else {
		euler_integration = false;
	}

	if (cfg.keyExists("output")) {
		output_file = cfg.getValueOfKey<std::string>("output");
	}
	else {
		output_file = "c:/temp/";
	}

    n_gridded = false;
    n_const = 0.f;
	if (cfg.keyExists("n")) {
		if (cfg.keyIsNumber("n")) {
			n_const = cfg.getValueOfKey<double>("n");
		} else {
			n_gridded = true;
			n_file    = cfg.getValueOfKey<std::string>("n");
    		}
    }

    
	if (cfg.keyExists("K")) {
		infiltration = true;
		K_file       = cfg.getValueOfKey<std::string>("K");
		psi          = cfg.getValueOfKey<double>      ("psi");
		dtheta       = cfg.getValueOfKey<double>      ("dtheta");
		std::ifstream test_exist(K_file.c_str());											//added by Youcan on 20170831
		if (!test_exist.good()) std::cout << "Ksat_file: " << K_file << " is not found.";	//added by Youcan on 20170831
	} else {
		infiltration = false;
	}

	if (cfg.keyExists("h0")) {
		h_file = cfg.getValueOfKey<std::string>("h0");
		std::cout<< "going to initialize"<<std::endl;
		h_init = true;
	} else {
		h_init = false;
	}

	if (cfg.keyExists("hydrograph")) {
		dambreak = true;
		hydrograph_file = cfg.getValueOfKey<std::string>("hydrograph");
		hydrograph.ReadSource(hydrograph_file.c_str());
		// Convert the source locations to grid coordinates
		source_X = cfg.getValueOfKey<double>("source_X");
		source_Y = cfg.getValueOfKey<double>("source_Y");
		std::cout << "X location is  " << source_X << " and Y location is " << source_Y << std::endl; 
	} else {
		dambreak = false;
	}

	if (cfg.keyExists("hyetograph")) {
		rainfall_averaged = true;
		hyetograph_file   = cfg.getValueOfKey<std::string>("hyetograph");
		hyetograph.ReadSource(hyetograph_file);
	} else {
		rainfall_averaged = false;
	}

	if (cfg.keyExists("hyetograph_prefix")) {
		rainfall_gridded  = true;
		hyetograph_t      = t0;
		hyetograph_prefix = cfg.getValueOfKey<std::string>("hyetograph_prefix");
		hyetograph_dt     = cfg.getValueOfKey<double>      ("hyetograph_dt");
		hyetograph_tf     = cfg.getValueOfKey<double>      ("hyetograph_tf");

		// added by Youcan on 20170831
		double t_temp = 0.;
		while (t_temp <= hyetograph_tf) {
			std::stringstream hyetograph_file_ss;
			std::string hyetograph_file;
			hyetograph_file_ss << hyetograph_prefix << "-" << std::setw(3)
				<< std::setfill('0') << (int)(t_temp / 3600.0)
				<< ".txt"; //hyetograph_prefix-001.txt	
			hyetograph_file = hyetograph_file_ss.str();
			std::ifstream test_exist(hyetograph_file.c_str());
			if (!test_exist.good()) std::cout << "Hyetograph: " << hyetograph_file << " is not found.";
			t_temp += hyetograph_dt;
		}

	} else {
		rainfall_gridded  = false;
	}

	if (cfg.keyExists("h_print")) {
		h_print = cfg.getValueOfKey<bool>("h_print");
	} else {
		h_print = false;
	}
	if (cfg.keyExists("save_max")) {
		save_max = cfg.getValueOfKey<bool>("save_max");
	}else {
		save_max = false;
	}
       
        if (cfg.keyExists("save_arrival")){

                save_arrival_time = cfg.getValueOfKey<bool>("save_arrival");
        } else{
                save_arrival_time = false;
       }


	if (cfg.keyExists("q_print")) {
		q_print = cfg.getValueOfKey<bool>("q_print");
	} else {
		q_print = false;
	}

	if (cfg.keyExists("check_volume")) {
		check_volume = cfg.getValueOfKey<bool>("check_volume");
	} else {
		check_volume = false;
	}

	if (h_print || q_print || check_volume) {
		t_print     = t0;
		dt_print    = cfg.getValueOfKey<double>("dt_print");
		count_print = 0;
	}

	if (cfg.keyExists("device")) {
		device = cfg.getValueOfKey<int>("device");
	} else {
		device = 0;
	}


	V_added = 0.f;

}

void Simulator::InitSimulation(void) {
	cudaSetDevice(device);
	steady_state = true;
	volume_old = 0;

    // Load the user-defined bathymetry
    InitBathymetry(b, DEM_file);
	B = ReadGrid(DEM_file);

	//if (cfg.keyExists("kappa")) {
	//	h_kappa = cfg.getValueOfKey<double>("kappa");
	//} else {
		// Set the desingularization constant (recommendation of user)
		kappa = sqrtf(0.01f*fmaxf(1.f, fminf(h_dx, h_dy)));
	//}

	if (h_init || h_print) {
		h_h = (double*)malloc(h_nx*h_ny*sizeof(double));
		memset(h_h, 0, h_nx*h_ny*sizeof(double));
		if (h_init) {
			h_o  = (double*)malloc(b_nx*b_ny*sizeof(double));
			ReadOriginalGrid(h_o, h_file);
		}
	}

	if (rainfall_gridded) {
		h_hyetograph = (double*)malloc(h_nx*h_ny*sizeof(double));
		memset(h_hyetograph, 0, h_nx*h_ny*sizeof(double));

		std::stringstream hyetograph_file_ss;
		std::string hyetograph_file;
		//hyetograph_file_ss << hyetograph_prefix << "-" << std::setw(3)					//removed by Youcan on 20170612
		//				   << std::setfill('0') << 6 + (int)(hyetograph_t / 3600.0)			//removed by Youcan on 20170612
		//				   << ".asc"; //hyetograph_prefix-001.txt							//removed by Youcan on 20170612
		hyetograph_file_ss << hyetograph_prefix << "-" << std::setw(3)						//added by Youcan on 20170612
			<< std::setfill('0') << (int)(hyetograph_t / 3600.0)							//added by Youcan on 20170612
			<< ".txt"; //hyetograph_prefix-001.txt											//added by Youcan on 20170612

		hyetograph_file = hyetograph_file_ss.str();

		hyetograph_o = (double*)malloc(b_nx*b_ny*sizeof(double));
		ReadOriginalGrid(hyetograph_o, hyetograph_file);
	}
	
    h_n = (double*)malloc(h_nx*h_ny*sizeof(double));
    memset(h_n, 0, h_nx*h_ny*sizeof(double));
    if (n_gridded) {
        n_o = (double*)malloc(b_nx*b_ny*sizeof(double));
		ReadOriginalGrid(n_o, n_file);
    } 
	
	if (infiltration) {
		K_o = (double*)malloc(b_nx*b_ny*sizeof(double));
		h_K = (double*)malloc(h_nx*h_ny*sizeof(double));
		ReadOriginalGrid(K_o, K_file);
	}

	
	
	// int i;
		if (dambreak) {
					std::cout << "This is a dam break" << std::endl; 

			NumSources = 1;
			InitializeSources(NumSources);
			int j;
			for (j=0; j<=NumSources;j++){
				SetSourceLocation(j,source_X,source_Y);
			}
	}

	

	SetDeviceConstants(B->num_columns, B->num_rows, B->cellsize, kappa);


	AllocateGrid(w, hu, hv, w_old, hu_old, hv_old, dw, dhu, dhv, mx, BC, BX, BY,
                 wet_blocks, active_blocks, n, hyetograph_gridded_rate, F,
                 F_old, dF, K, h, q, h_max, q_max, t_wet, dambreak,
	             rainfall_averaged, rainfall_gridded, infiltration, 
	             euler_integration, check_volume, h_init, h_print, q_print,
	             save_max, save_arrival_time, psi, dtheta, time_peak, time_dry, G);	//added time_peak and time_dry by Youcan on 20170908



    h_BX = (double*)malloc(h_ny*(h_nx+1)*sizeof(double));
    h_BY = (double*)malloc((h_ny+1)*h_nx*sizeof(double));
    memset(h_BX, 0, h_ny*(h_nx+1)*sizeof(double));
    memset(h_BY, 0, (h_ny+1)*h_nx*sizeof(double));

    // Interpolate gridded interfacial points
    for (int j = 2; j < h_ny - 2; j++) {
        for (int i = 2; i < h_nx - 2; i++) {
            int jt = j - 2, it = i - 2;

            int bed = j     * (h_nx+1) + i+1;
            int bwd = j     * (h_nx+1) + i;
            int bnd = (j+1) * (h_nx)   + i;
            int bsd = j     * (h_nx)   + i;

            int ll = jt     * b_nx + it;     // lower-left bathymetry point
            int ul = (jt+1) * b_nx + it;     // upper-left bathymetry point
            int ur = (jt+1) * b_nx + (it+1); // upper-right bathymetry point
            int lr = jt     * b_nx + (it+1); // lower-right bathymetry point

            h_BX[bed] = 0.5f * (B->data[ur]+B->data[lr]);
			h_BX[bwd] = 0.5f * (B->data[ll]+B->data[ul]);
            h_BY[bnd] = 0.5f * (B->data[ul]+B->data[ur]);
			h_BY[bsd] = 0.5f * (B->data[ll]+B->data[lr]);

			if (i == 2) {
				h_BX[j*(h_nx+1)+(i-1)] = h_BX[bed];
				h_BX[j*(h_nx+1)+(i-2)] = h_BX[bwd];
			}  if (i == h_nx-3) {
				h_BX[j*(h_nx+1)+(i+2)] = h_BX[bwd];
			}  if (j == 2) {
				h_BY[(j-1)*h_nx+i] = h_BY[bnd];
				h_BY[(j-2)*h_nx+i] = h_BY[bsd];
			}  if (j == h_ny-3) {
				h_BY[(j+2)*h_nx+i] = h_BY[bsd];
			}

			if (h_init) {
				h_h[j*h_nx+i] = 0.25f * (h_o[ur]+h_o[lr]+h_o[ll]+h_o[ul]);
				
			}

			if (rainfall_gridded) {
				h_hyetograph[j*h_nx+i] = 0.25f*(hyetograph_o[ur]+hyetograph_o[lr]+
												hyetograph_o[ll]+hyetograph_o[ul]);
				h_hyetograph[j*h_nx+i] /= (3600.f*1000.f);
			}

			if (n_gridded) {
				h_n[j*h_nx+i] = 0.25f * (n_o[ur]+n_o[lr]+n_o[ll]+n_o[ul]);
			}  else {
                h_n[j*h_nx+i] = n_const;
            }

			if (infiltration) {
				h_K[j*h_nx+i] = 0.25f * (K_o[ur]+K_o[lr]+K_o[ll]+K_o[ul]);
				h_K[j*h_nx+i] /= (3600.f*1000.f);
			}
        }
    }

    checkCudaErrors(cudaMemcpy2D(BX, pitchBX, h_BX, (h_nx+1)*sizeof(double),
                                 (h_nx+1)*sizeof(double), h_ny, HtoD));
    checkCudaErrors(cudaMemcpy2D(BY, pitchBY, h_BY, h_nx*sizeof(double),
                                 h_nx*sizeof(double), (h_ny+1), HtoD));

	if (h_init) {
		std::cout << "Copying initial Grid to Device" << std::endl;
		checkCudaErrors(cudaMemcpy2D(h, pitch, h_h, h_nx*sizeof(double),
									 h_nx*sizeof(double), h_ny, HtoD));
	}

	if (rainfall_gridded) {
		checkCudaErrors(cudaMemcpy2D(hyetograph_gridded_rate, pitch, h_hyetograph,
									 h_nx*sizeof(double), h_nx*sizeof(double), h_ny,
									 HtoD));
	}

    checkCudaErrors(cudaMemcpy2D(n, pitch, h_n, h_nx*sizeof(double),
                                 h_nx*sizeof(double), h_ny, HtoD));

	if (infiltration) {
		checkCudaErrors(cudaMemcpy2D(K, pitch, h_K, h_nx*sizeof(double),
									 h_nx*sizeof(double), h_ny, HtoD));
	}

	if (q_print) {
		h_q  = (double*)malloc(h_nx*h_ny*sizeof(double));
	}

	if (check_volume && infiltration) {
		h_F = (double*)malloc(h_nx*h_ny*sizeof(double));
	}

	InitGrid(w, hu, hv, w_old, hu_old, hv_old, BC, BX, BY, wet_blocks,
	         active_blocks, h, t_wet, G);
}

void Simulator::StartTimer() {
	gettimeofday(&start_time, NULL);
	//start_time = clock();		//for windows, Youcan
}

void Simulator::UpdateSource(void) {
	
	if (dambreak) {
		// Interpolate source discharge rate
		hydrograph.InterpolateRate(t);

		// Convert hydrograph from CFS to m/s
		//std::cout << "getting discharge rate cfs " << hydrograph.interpolated_rate << std::endl; 
		hydrograph.interpolated_rate /= 35.3147f;
		//std::cout << "getting discharge height m " << hydrograph.interpolated_rate << std::endl; 
                //
		setSourceValue(0,hydrograph.interpolated_rate);
		
		/*if (check_volume) {
			V_added += hydrograph.interpolated_rate * dt;
		}
		int i;*/

		updateSources();
		//checkCudaErrors(cudaMemcpy(source_rate_dev,source_rate,NumSources*sizeof(double),HtoD));
	}

	if (rainfall_averaged) {
		// Interpolate rainfall rate
		hyetograph.InterpolateRate(t);
		// Convert rainfall from mm/hr to m/s
		hyetograph.interpolated_rate /= (3600.f*1000.f);

		if (check_volume) {
			V_added += hyetograph.interpolated_rate * dt *
					   (double)((h_nx-4) * (h_ny-4));
		}
	}

	if (rainfall_gridded) {
		// Interpolate rainfall rate
		if (t > hyetograph_t) {
			hyetograph_t += hyetograph_dt;
		}
	}
}

void Simulator::ComputeTimestep() {
	//std::cout << "compute time step" << std::endl; 
    thrust::device_ptr<double> mxptr(mx);
    thrust::device_ptr<double> mxresptr;
	mxresptr=thrust::max_element(mxptr,mxptr + GridSize);
    //mxresptr = thrust::max_element(mxptr, mxptr + GridSize);

    double max_x = mxresptr[0];

    double c_x = h_dx / fmaxf(4.f*max_x, kappa);
    dt = fminf(c_x, h_dx/10.f);

	dt = (t == 0.f) ? dt = 0.000001f : dt;	
	t += dt;
    // std::cout << "max_x = " << max_x
    //           << ", at mxresptr = " << mxresptr - mxptr
    //           << std::endl;
	// std::cout << "time step is " << dt << std::endl; 
}

void writeHeader(std::ofstream &thefile) {
    thefile.precision(dbl::digits10);
	thefile << "ncols         " << h_nx - 4          << std::endl;
	thefile << "nrows         " << h_ny - 4          << std::endl;
	thefile << "xllcorner     " << h_xll             << std::endl;
	thefile << "yllcorner     " << h_yll             << std::endl;
	thefile << "cellsize      " << cellsize_original << std::endl;
	thefile << "NODATA_value  " << -9999             << std::endl;
    thefile.precision(flt::digits10);
}

void Simulator::PrintData(void) {
	wet_count = 0;
	if (h_print || check_volume) {
		checkCudaErrors(cudaMemcpy2D(h_h, h_nx*sizeof(double), h, pitch,
									 h_nx*sizeof(double), h_ny, DtoH));
	}
	if (save_max) {
		PrintSummaryData();
	}
	if (h_print) {
		std::stringstream filename_h;
		std::cout<<"Interpolated Flow Rate: "<< hydrograph.interpolated_rate << std:: endl;
		filename_h << output_file << "/h" << count_print << ".txt";
		std::ofstream heights;
		heights.open((filename_h.str()).c_str());
        if (heights.is_open()) {
          writeHeader(heights);
          for (int j = h_ny - 3; j >= 2; j--) {
              for (int i = 2; i < h_nx-2; i++) {
                  int   id = j*h_nx + i;
                  heights << h_h[id] << " ";
                  if (h_h[id] > 0) wet_count++;
              }
              heights << std::endl;
          }
        } else {
            std::string msg(filename_h.str());
            msg += std::string(": error: cannot open for writing");
            throw std::runtime_error(msg);
        }
		heights.close();
	}

	if (check_volume) {
		V_computed = 0.f;

		if (infiltration) {
			checkCudaErrors(cudaMemcpy2D(h_F, h_nx*sizeof(double), F, pitch,
										 h_nx*sizeof(double), h_ny, DtoH));
		}
	}

	if (check_volume) {
		
		#pragma omp parallel for reduction(+:V_computed)
		
		for (int j = h_ny - 3; j >= 2; j--) {
			for (int i = 2; i < h_nx-2; i++) {
				int   id       = j*h_nx + i;
				V_computed += h_h[id];
				if (infiltration) {
					V_computed += h_F[id];
				}
			}
		}
		
		double percent = abs((volume_old - V_computed) / (volume_old) * 100);
		//std::cout << "time    " << t << "     Volume   " << V_computed << std::endl;
		
		//if (t > 60 && percent < .02) steady_state = false;
		volume_old = V_computed;
	}

	if (q_print) {
		checkCudaErrors(cudaMemcpy2D(h_q, h_nx*sizeof(double), q, pitch,
									 h_nx*sizeof(double), h_ny, DtoH));
		std::stringstream filename_q;
		filename_q << output_file << "/q" << count_print << ".txt";
		std::ofstream discharge;
		discharge.open((filename_q.str()).c_str());
        if (discharge.is_open()) {
            writeHeader(discharge);
            for (int j = h_ny - 3; j >= 2; j--) {
                for (int i = 2; i < h_nx-2; i++) {
                    int id = j*h_nx + i;
                    discharge << h_q[id] << " ";
                }
                discharge << std::endl;
            }
            discharge.close();
        } else {
            std::string msg(filename_q.str());
            msg += std::string(": error: cannot open for writing");
            throw std::runtime_error(msg);
        }
	}

    count_print++;
}

void Simulator::PrintSummaryData(void) {
	if (save_max) {
		double *h_h_max = (double*)malloc(h_nx*h_ny*sizeof(double));
		double *h_q_max = (double*)malloc(h_nx*h_ny*sizeof(double));
		double *h_t_peak = (double*)malloc(h_nx*h_ny * sizeof(double));						// added by Youcan on 20170424
		double *h_t_dry = (double*)malloc(h_nx*h_ny * sizeof(double));						// added by Youcan on 20170830
		checkCudaErrors(cudaMemcpy2D(h_h_max,  h_nx*sizeof(double), h_max, pitch,
									 h_nx*sizeof(double), h_ny, DtoH));
		checkCudaErrors(cudaMemcpy2D(h_q_max,  h_nx*sizeof(double), q_max,  pitch,
									 h_nx*sizeof(double), h_ny, DtoH));
		checkCudaErrors(cudaMemcpy2D(h_t_peak, h_nx * sizeof(double), time_peak, pitch,		// added by Youcan on 20170424
									 h_nx * sizeof(double), h_ny, DtoH));					// added by Youcan on 20170424
		checkCudaErrors(cudaMemcpy2D(h_t_dry, h_nx * sizeof(double), time_dry, pitch,		// added by Youcan on 20170830
									 h_nx * sizeof(double), h_ny, DtoH));					// added by Youcan on 20170830
		std::ofstream hmax, qmax, tmax, tdry;												// added tmax, tdry by Youcan on 20170908
		std::stringstream hmax_name, qmax_name, tmax_name, tdry_name;						// added tmax_name, tdry_name by Youcan on 20170908
		hmax_name << output_file << "/peak_flood_depth.txt";
		qmax_name << output_file << "/peak_unit_flow.txt";
		tmax_name << output_file << "/time_to_peak.txt";									// added by Youcan on 20170424
		tdry_name << output_file << "/time_to_dry.txt";										// added by Youcan on 20170830
		hmax.open((hmax_name.str()).c_str());
        if (!hmax.is_open()) {
            std::string msg(hmax_name.str());
            msg += ": error: cannot open for writing";
            throw std::runtime_error(msg);
        }

		qmax.open((qmax_name.str()).c_str());
        if (!qmax.is_open()) {
            std::string msg(qmax_name.str());
            msg += ": error: cannot open for writing";
            throw std::runtime_error(msg);
        }

		tmax.open((tmax_name.str()).c_str());									        if (!tmax.is_open()) {
            std::string msg(tmax_name.str());
            msg += ": error: cannot open for writing";
            throw std::runtime_error(msg);
        }
        
			// added by Youcan on 20170424
		tdry.open((tdry_name.str()).c_str());												// added by Youcan on 20170830
        if (!tdry.is_open()) {
            std::string msg(tdry_name.str());
            msg += ": error: cannot open for writing";
            throw std::runtime_error(msg);
        }

		writeHeader(hmax);
		writeHeader(qmax);
		writeHeader(tmax);																	// added by Youcan on 20170424
		writeHeader(tdry);																	// added by Youcan on 20170830


		for (int j = h_ny - 3; j >= 2; j--) {
			for (int i = 2; i < h_nx - 2; i++) {
				int id = j*h_nx + i;
				hmax << h_h_max[id] << " ";
				qmax << h_q_max[id] << " ";
				tmax << h_t_peak[id] << " ";												// added by Youcan on 20170424
				tdry << h_t_dry[id] << " ";													// added by Youcan on 20170830

			}
			hmax << std::endl;
			qmax << std::endl;
			tmax << std::endl;																// added by Youcan on 20170424
			tdry << std::endl;																// added by Youcan on 20170830

		}
	}

	if (save_arrival_time) {
		double *h_t_wet = (double*)malloc(h_nx*h_ny*sizeof(double));
		checkCudaErrors(cudaMemcpy2D(h_t_wet,  h_nx*sizeof(double), t_wet,  pitch,
									 h_nx*sizeof(double), h_ny, DtoH));
		std::ofstream twet;
		std::stringstream arrival_name;
		arrival_name << output_file << "/flood_wave_arrival.txt";
		
		twet.open((arrival_name.str()).c_str());
		writeHeader(twet);

		for (int j = h_ny - 3; j >= 2; j--) {
			for (int i = 2; i < h_nx - 2; i++) {
				int id = j*h_nx + i;
				twet << h_t_wet[id] << " ";
			}
			twet << std::endl;
		}
		twet.close();
	}
}

double Simulator::EndTimer() {
    long seconds, useconds;
    double elapsed;
	//end_time = clock();										// for windows only, Youcan
	//seconds = (end_time - start_time) / CLOCKS_PER_SEC;		// for windows only, Youcan
	//useconds = seconds * 1000;								// for windows only, Youcan
	gettimeofday(&end_time, NULL);
	seconds  = end_time.tv_sec-start_time.tv_sec;
	useconds = end_time.tv_usec-start_time.tv_usec;

/*	end_time = clock();		//for windows, Youcan
	seconds = (end_time - start_time) / CLOCKS_PER_SEC;	//for windows, Youcan
	useconds = seconds * 1000;		//for windows, Youcan
*/

    elapsed  = (seconds*1000 + useconds/1000.0) + 0.5;
	elapsed /= 1000.0;
	return elapsed;
}
/**ADDED for 1d/2d Integration*/

void Simulator::OpenSimulation(std::string filename, double flowrate){
flow_rate = flowrate;
std::cout<<"inside the opend simulation" << std::endl;	
	ReadUserParams(filename);
	std::cout<< "read the user parameters" <<std::endl;
    InitSimulation();
	std::cout<< "initialized the simulation" << std::endl;
    ApplyBoundaries(w, hu, hv, BC);
	count = 0;
 StartTimer();

}
void Simulator::CloseSimulation(){
	if (save_max) {
		PrintSummaryData();
	}

	FreeGrid(w, hu, hv, w_old, hu_old, hv_old, dw, dhu, dhv, mx, BC, BX, BY,
	         wet_blocks, active_blocks, n, hyetograph_gridded_rate, F, F_old,
             dF, K, h, q, h_max, q_max, t_wet, time_peak, time_dry, G);		//added t_peak and t_dry by Youcan on 20170908

    free(b);
    free(h_BX);
    free(h_BY);
	free(h_h);
	free(h_q);
    free(h_hyetograph);
	free(h_n);
	free(h_K);
	free(h_F);
	free(h_o);
    free(hyetograph_o);
	free(n_o);
	free(K_o);
	free(h_h_max);
	free(h_q_max);
	free(h_t_wet);
	free(B->data);

	/** Added for 1D-2D*/
	free(source_idx);
	free(source_rate);
	cudaFree(source_idx_dev);
	cudaFree(source_rate_dev);
}

double Simulator::RunSimulation() {
   
	
	// tf; 
   // while (t < tf) {
	//updateSources();
	checkCudaErrors(cudaMemcpy(source_idx_dev,source_idx,NumSources*sizeof(int),HtoD));

		if (dambreak || rainfall_averaged) {
			//std::cout << "updating source" << std::endl; 
			UpdateSource();
		}
		//std::cout << "Grow Blocks" << std::endl; 
		Grow(wet_blocks, active_blocks, hyetograph_gridded_rate, rainfall_gridded);
		//std::cout << "Compute Fluxes" << std::endl; 
        ComputeFluxes(w, hu, hv, dw, dhu, dhv, mx, BC, BX, BY, G, active_blocks,
		              dt, n, hydrograph.interpolated_rate, dambreak_source_idx,
		              hyetograph.interpolated_rate, hyetograph_gridded_rate, F,
					  F_old, dF, K,source_idx_dev,source_rate_dev,NumSources);
		ComputeTimestep();

		Integrate_1(w, hu, hv, w_old, hu_old, hv_old, dw, dhu, dhv, BC,
					G, wet_blocks, active_blocks, t, dt,
                    hydrograph.interpolated_rate, dambreak_source_idx,
                    hyetograph.interpolated_rate, hyetograph_gridded_rate, F,
                    F_old, dF, K, h, q, h_max, q_max, t_wet,source_idx_dev, source_rate_dev,NumSources, time_peak, time_dry);	//added time_peak and time_dry by Youcan on 20170908
		ApplyBoundaries(w, hu, hv, BC);

		if (!euler_integration) {
			
			ComputeFluxes(w, hu, hv, dw, dhu, dhv, mx, BC, BX, BY, G, active_blocks,
						  dt, n, hydrograph.interpolated_rate, dambreak_source_idx,
						  hyetograph.interpolated_rate, hyetograph_gridded_rate, F,
						  F_old, dF, K,source_idx_dev, source_rate_dev,NumSources);
			Integrate_2(w, hu, hv, w_old, hu_old, hv_old, dw, dhu, dhv, BC,
						G, wet_blocks, active_blocks, t, dt,
						hydrograph.interpolated_rate, dambreak_source_idx,
						hyetograph.interpolated_rate, hyetograph_gridded_rate, F,
						F_old, dF, K, h, q, h_max, q_max, t_wet,source_idx_dev, source_rate_dev,NumSources);
			ApplyBoundaries(w, hu, hv, BC);
		}

        if (t > hyetograph_tf && rainfall_gridded == true) {

            memset(h_hyetograph, 0, h_nx*h_ny*sizeof(double));
            checkCudaErrors(cudaMemcpy2D(hyetograph_gridded_rate, pitch, h_hyetograph,
                                         h_nx*sizeof(double), h_nx*sizeof(double), h_ny,
                                         HtoD));
            rainfall_gridded = false;
        } else if (t > hyetograph_t && rainfall_gridded == true) {
            hyetograph_t += hyetograph_dt;

            std::stringstream hyetograph_file_ss;
            std::string hyetograph_file;
            //hyetograph_file_ss << hyetograph_prefix << "-" << std::setw(3)					//removed by Youcan on 20170612
            //                   << std::setfill('0') << (int)(hyetograph_t / 3600.0)			//removed by Youcan on 20170612
            //                   << ".asc";														//removed by Youcan on 20170612
			hyetograph_file_ss << hyetograph_prefix << "-" << std::setw(3)						//added by Youcan on 20170612
				<< std::setfill('0') << (int)(hyetograph_t / 3600.0)							//added by Youcan on 20170612
				<< ".txt"; //hyetograph_prefix-001.txt											//added by Youcan on 20170612

            hyetograph_file = hyetograph_file_ss.str();

            SetOriginalGrid(hyetograph_o, hyetograph_file);

            for (int j = 2; j < h_ny - 2; j++) {
                for (int i = 2; i < h_nx - 2; i++) {
                    int jt = j - 2, it = i - 2;

                    int ll = jt     * b_nx + it;     // lower-left bathymetry point
                    int ul = (jt+1) * b_nx + it;     // upper-left bathymetry point
                    int ur = (jt+1) * b_nx + (it+1); // upper-right bathymetry point
                    int lr = jt     * b_nx + (it+1); // lower-right bathymetry point

                    h_hyetograph[j*h_nx+i] = 0.25f*(hyetograph_o[ur]+hyetograph_o[lr]+
                                                    hyetograph_o[ll]+hyetograph_o[ul]);
					V_added += h_hyetograph[j*h_nx + i] / 1000.f;					//added by Youcan on 20170612
                    h_hyetograph[j*h_nx+i] /= (3600.f*1000.f);
                }
            }
            checkCudaErrors(cudaMemcpy2D(hyetograph_gridded_rate, pitch, h_hyetograph,
                                         h_nx*sizeof(double), h_nx*sizeof(double), h_ny,
                                         HtoD));
        }

		if (h_print || q_print || check_volume) {
			if (t > t_print) {
				PrintData();
				t_print += dt_print;
				std::cout << t << "\t" << dt;
				double end_time = EndTimer();
				long int computed_cells = (long int)(h_nx-3)*(long int)(h_ny-3)*(long int)count;
				
				if (check_volume) {
					std::cout << "\t" << V_added*h_dx*h_dy << "\t" << V_computed*h_dx*h_dy <<
					"\t" << (V_computed - V_added) / V_added << "\t" << b_nx*b_ny << "\t"  << end_time << "\t" << (double)computed_cells / end_time /1000000.0 << "\t" << wet_count;
				}
				std::cout << std::endl;
			}
			count++;
		}
	//}

	return t;

    
									  
}
	/**
	Changes for 1D-2D

	*/
void Simulator::InitializeSources(long numSources){
	NumSources = numSources; //set the public variable
	source_idx = (int*)malloc(numSources*sizeof(int));
	source_rate = (double*)malloc(numSources*sizeof(double));
	memset(source_idx,0,numSources*sizeof(int));
	memset(source_rate,0,numSources*sizeof(double));
	checkCudaErrors(cudaMalloc((void**)&source_idx_dev,numSources*sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&source_rate_dev,numSources*sizeof(double)));
	
}

void Simulator::SetSourceLocation(int i,double X, double Y){
	source_x = (int)((X-h_xll) / cellsize_original + 0.5);
	std::cout << "source_x is  " << source_x << std::endl; 
	std::cout << "X is  " << X << std::endl;
	std::cout << "h_xll is  " << h_xll << std::endl;
	std::cout << "cellsize_original is  " << cellsize_original << std::endl;
	source_y = (int)((Y-h_yll) / cellsize_original + 0.5);
	std::cout << "Y location is  " << source_y << std::endl; 
	int idx = source_y*h_nx + source_x;
	std::cout << "hnx is  " << h_nx << " idx is " << idx << std::endl;
	source_idx[i] =  idx;
		
}

void Simulator::setSourceValue (int i, double value){
	

	source_rate[i] = value/h_dx/h_dy;
	//std::cout << source_rate[i] <<std::endl;
	if (check_volume) {
		V_added += source_rate[i] * dt;
	}
	
	//printf("\nsource %g",source_rate[i]);
}

void Simulator::updateSources(){
	
	checkCudaErrors(cudaMemcpy(source_rate_dev, source_rate, NumSources*sizeof(double), HtoD));
	//std::cout << "updated the source" << std::endl;
}

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
