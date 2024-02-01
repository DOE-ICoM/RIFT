#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
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
typedef std::numeric_limits<double> flt;

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

    courant_max = 0.25;
    if (cfg.keyExists("courant_limit")) {
        courant_max = cfg.getValueOfKey<double>("courant_limit");
        if (courant_max <=0 || courant_max >= 1.0) {
            std::cerr << "Specifed Courant_limit ("
                      << courant_max << ") "
                      << "is out of range, reverting to default"
                      << std::endl;
            courant_max = 0.25;
        }
    }
    std::cout << "Courant limit = " << courant_max << std::endl;
    
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
        NumSources = 0;
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
		hyetograph_prefix = cfg.getValueOfKey<std::string>("hyetograph_prefix");
		hyetograph_dt     = cfg.getValueOfKey<double>      ("hyetograph_dt");
		hyetograph_tf     = cfg.getValueOfKey<double>      ("hyetograph_tf");
	} else {
		rainfall_gridded  = false;
	} 

	if (cfg.keyExists("drain")) {
        drain_averaged = true;
        drain_file = cfg.getValueOfKey<std::string>("drain");
        drain.ReadSource(drain_file);
    } else {
        drain_averaged = false;
    }

	if (cfg.keyExists("drain_prefix")) {
		drain_gridded  = true;
		drain_prefix = cfg.getValueOfKey<std::string>("drain_prefix");
		drain_dt     = cfg.getValueOfKey<double>      ("drain_dt");
		drain_tf     = cfg.getValueOfKey<double>      ("drain_tf");
	} else {
		drain_gridded  = false;
	} 

    if (cfg.keyExists("surge_prefix")) {
        surge_gridded = true;
        surge_prefix = cfg.getValueOfKey<std::string>("surge_prefix");
        surge_dt = cfg.getValueOfKey<double> ("surge_dt");
        surge_tf = cfg.getValueOfKey<double> ("surge_tf");
    } else {
        surge_gridded = false;
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

    if (cfg.keyExists("b_print")) {
        b_print = true;
    } else {
        b_print = false;
    }

    if (cfg.keyExists("n_print")) {
        n_print = true;
    } else {
        n_print = false;
    }

	if (cfg.keyExists("device")) {
		device = cfg.getValueOfKey<int>("device");
	} else {
		device = 0;
	}

    if (cfg.keyExists("square_cells")) {
        square_cells = cfg.getValueOfKey<bool>("square_cells");
    } else {
        square_cells = true;
    }

	V_added = 0.f;

}

void Simulator::InitSimulation(void) {
	cudaSetDevice(device);
	steady_state = true;
	volume_old = 0;

    // Load the user-defined bathymetry
    InitBathymetry(b, DEM_file, this->grid_config, square_cells);
	f_b = 1; // b is allocated
	B = ReadGrid(DEM_file);

	//if (cfg.keyExists("kappa")) {
	//	h_kappa = cfg.getValueOfKey<double>("kappa");
	//} else {
		// Set the desingularization constant (recommendation of user)
		kappa = sqrtf(0.01f*fmaxf(1.f, fminf(grid_config.h_dx, grid_config.h_dy)));
	//}

	if (h_init || h_print) {
		h_h = (double*)malloc(grid_config.h_nx * grid_config.h_ny * sizeof(double));
		f_h_h = 1;
		memset(h_h, 0, grid_config.h_nx * grid_config.h_ny * sizeof(double));
		if (h_init) {
			h_o  = (double*)malloc(grid_config.b_nx * grid_config.b_ny * sizeof(double));
			f_h_o = 1;
			SetOriginalGrid(h_o, h_file, this->grid_config);
		}
	}

    h_n = (double*)malloc(grid_config.h_nx*grid_config.h_ny*sizeof(double));
	f_h_n = 1;
    memset(h_n, 0, grid_config.h_nx*grid_config.h_ny*sizeof(double));
    if (n_gridded) {
        n_o = (double*)malloc(grid_config.b_nx*grid_config.b_ny*sizeof(double));
		f_n_o = 1;
		SetOriginalGrid(n_o, n_file, this->grid_config);
    } 
	
	if (infiltration) {
		K_o = (double*)malloc(grid_config.b_nx*grid_config.b_ny*sizeof(double));
		f_K_o = 1;
		h_K = (double*)malloc(grid_config.h_nx*grid_config.h_ny*sizeof(double));
		f_h_K = 1;
		SetOriginalGrid(K_o, K_file, this->grid_config);
	}

	
	
    
		if (dambreak) {
					std::cout << "This is a dam break" << std::endl; 

			NumSources = 1;
			InitializeSources(NumSources);
			int j;
			for (j=0; j<NumSources;j++){
				SetSourceLocation(j,source_X,source_Y);
			}
        } else {
            NumSources = 0;
        }

	SetDeviceConstants(grid_config.h_nx, grid_config.h_ny,
                       grid_config.h_dx, grid_config.h_dy,
                       kappa, grid_config.nodata);


	AllocateGrid(dev_w, dev_hu, dev_hv, dev_w_old, dev_hu_old, dev_hv_old, 
				 dev_dw, dev_dhu, dev_dhv, dev_mx, dev_my, dev_BC, dev_BX, dev_BY,
                 dev_wet_blocks, dev_active_blocks, dev_n, 
				 dev_hyetograph_gridded_rate, dev_F, dev_F_old, dev_dF, dev_K,
				 dev_h, dev_q, dev_h_max, dev_q_max, dev_t_wet, dambreak,
	             rainfall_averaged, drain_averaged,
                 rainfall_gridded, drain_gridded, infiltration, 
				 surge_gridded, euler_integration, check_volume, h_init, h_print, q_print,
	             save_max, save_arrival_time, psi, dtheta, dev_time_peak, 
				 dev_time_dry, dev_G, grid_config);	//added time_peak and time_dry by Youcan on 20170908

	if (rainfall_gridded) {
        hyetograph_series.reset(new HyetographGridSeries(hyetograph_prefix,
                                                         hyetograph_dt,
                                                         hyetograph_tf,
                                                         grid_config,
                                                         dev_hyetograph_gridded_rate));
        hyetograph_series->update(t0);
	}

    if (drain_gridded) {
        drain_series.reset(new HyetographGridSeries(drain_prefix,
                                                    drain_dt,
                                                    drain_tf,
                                                    grid_config));
        drain_series->update(t0);
        dev_drain_gridded_rate = drain_series->grid_dev();
    }

    if (surge_gridded) {
        surge_series.reset(new InterpolatedGridSeries(surge_prefix, 1.0,
                                                      surge_dt, surge_tf,
                                                      grid_config));
        // surge data must be surrounded by nodata values
        surge_series->allow_no_data(true);
        surge_series->update(t0);
        dev_surge_gridded_depth = surge_series->grid_dev();
    } else {
        dev_surge_gridded_depth = NULL;
    }
	
    h_BX = (double*)malloc(grid_config.h_ny*(grid_config.h_nx+1)*sizeof(double));
	f_h_BX = 1;
    h_BY = (double*)malloc((grid_config.h_ny+1)*grid_config.h_nx*sizeof(double));
	f_h_BY = 1;
    memset(h_BX, 0, grid_config.h_ny*(grid_config.h_nx+1)*sizeof(double));
    memset(h_BY, 0, (grid_config.h_ny+1)*grid_config.h_nx*sizeof(double));

    // Interpolate gridded interfacial points
    for (int j = 2; j < grid_config.h_ny - 2; j++) {
        for (int i = 2; i < grid_config.h_nx - 2; i++) {
            int jt = j - 2, it = i - 2;

            int bed = j     * (grid_config.h_nx+1) + i+1;
            int bwd = j     * (grid_config.h_nx+1) + i;
            int bnd = (j+1) * (grid_config.h_nx)   + i;
            int bsd = j     * (grid_config.h_nx)   + i;

            int ll = jt     * grid_config.b_nx + it;     // lower-left bathymetry point
            int ul = (jt+1) * grid_config.b_nx + it;     // upper-left bathymetry point
            int ur = (jt+1) * grid_config.b_nx + (it+1); // upper-right bathymetry point
            int lr = jt     * grid_config.b_nx + (it+1); // lower-right bathymetry point

            h_BX[bed] = 0.5f * (B->data[ur]+B->data[lr]);
			h_BX[bwd] = 0.5f * (B->data[ll]+B->data[ul]);
            h_BY[bnd] = 0.5f * (B->data[ul]+B->data[ur]);
			h_BY[bsd] = 0.5f * (B->data[ll]+B->data[lr]);

			if (i == 2) {
				h_BX[j*(grid_config.h_nx+1)+(i-1)] = h_BX[bed];
				h_BX[j*(grid_config.h_nx+1)+(i-2)] = h_BX[bwd];
			}  if (i == grid_config.h_nx-3) {
				h_BX[j*(grid_config.h_nx+1)+(i+2)] = h_BX[bwd];
			}  if (j == 2) {
				h_BY[(j-1)*grid_config.h_nx+i] = h_BY[bnd];
				h_BY[(j-2)*grid_config.h_nx+i] = h_BY[bsd];
			}  if (j == grid_config.h_ny-3) {
				h_BY[(j+2)*grid_config.h_nx+i] = h_BY[bsd];
			}

			if (h_init) {
				h_h[j*grid_config.h_nx+i] = 0.25f * (h_o[ur]+h_o[lr]+h_o[ll]+h_o[ul]);
				
			}

			if (n_gridded) {
				h_n[j*grid_config.h_nx+i] = 0.25f * (n_o[ur]+n_o[lr]+n_o[ll]+n_o[ul]);
			}  else {
                h_n[j*grid_config.h_nx+i] = n_const;
            }

			if (infiltration) {
				h_K[j*grid_config.h_nx+i] = 0.25f * (K_o[ur]+K_o[lr]+K_o[ll]+K_o[ul]);
				h_K[j*grid_config.h_nx+i] /= (3600.f*1000.f);
			}
        }
    }

    checkCudaErrors(cudaMemcpy2D(dev_BX, pitchBX, h_BX, (grid_config.h_nx+1)*sizeof(double),
                                 (grid_config.h_nx+1)*sizeof(double), grid_config.h_ny, HtoD));
    checkCudaErrors(cudaMemcpy2D(dev_BY, pitchBY, h_BY, grid_config.h_nx*sizeof(double),
                                 grid_config.h_nx*sizeof(double), (grid_config.h_ny+1), HtoD));

	if (h_init) {
		std::cout << "Copying initial Grid to Device" << std::endl;
		checkCudaErrors(cudaMemcpy2D(dev_h, pitch, h_h, grid_config.h_nx*sizeof(double),
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, HtoD));
	}

    checkCudaErrors(cudaMemcpy2D(dev_n, pitch, h_n, grid_config.h_nx*sizeof(double),
                                 grid_config.h_nx*sizeof(double), grid_config.h_ny, HtoD));

	if (infiltration) {
		checkCudaErrors(cudaMemcpy2D(dev_K, pitch, h_K, grid_config.h_nx*sizeof(double),
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, HtoD));
	}

	if (q_print) {
		h_q  = (double*)malloc(grid_config.h_nx*grid_config.h_ny*sizeof(double));
		f_h_q = 1;
	}

	if (check_volume && infiltration) {
		h_F = (double*)malloc(grid_config.h_nx*grid_config.h_ny*sizeof(double));
		f_h_F = 1;
	}

	InitGrid(dev_w, dev_hu, dev_hv, dev_w_old, dev_hu_old, dev_hv_old, dev_BC, dev_BX, dev_BY, dev_wet_blocks,
	         dev_active_blocks, dev_h, dev_t_wet, dev_G);

    if (b_print) {
        std::unique_ptr<double[]>
            h_BC(new double[grid_config.h_nx * grid_config.h_ny]);
        checkCudaErrors(cudaMemcpy2D(h_BC.get(), grid_config.h_nx*sizeof(double),
                                     dev_BC, pitch,
									 grid_config.h_nx*sizeof(double),
                                     grid_config.h_ny, DtoH));
		std::stringstream filename_h;
		filename_h << output_file << "/bc.txt";
        writeGrid(filename_h.str(), h_BC.get(), grid_config);
    }

    if (n_print) {
        std::unique_ptr<double[]>
            h_n(new double[grid_config.h_nx * grid_config.h_ny]);
        checkCudaErrors(cudaMemcpy2D(h_n.get(), grid_config.h_nx*sizeof(double),
                                     dev_n, pitch,
									 grid_config.h_nx*sizeof(double),
                                     grid_config.h_ny, DtoH));
		std::stringstream filename_h;
		filename_h << output_file << "/n.txt";
        writeGrid(filename_h.str(), h_n.get(), grid_config);
    }
        
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
					   (double)((grid_config.h_nx-4) * (grid_config.h_ny-4));
		}
	}

    if (drain_averaged) {
        drain.InterpolateRate(t);
        drain.interpolated_rate /= (3600.f*1000.f);
    }        
}

void Simulator::ComputeTimestep() {
	//std::cout << "compute time step" << std::endl; 
    thrust::device_ptr<double> mxptr(dev_mx), myptr(dev_my);
    thrust::device_ptr<double> mxresptr, myresptr;
	mxresptr=thrust::max_element(mxptr,mxptr + GridSize);
	myresptr=thrust::max_element(myptr,myptr + GridSize);

    double max_x = mxresptr[0], max_y = myresptr[0];

    double c_x = grid_config.h_dx / fmaxf(max_x/courant_max, kappa);
    double c_y = grid_config.h_dy / fmaxf(max_y/courant_max, kappa);

    // WAP: This seems like a very arbitrary (and rather small)
    // maximum limit on time step.  I cannot identify a source for
    // it. What does this add over the use of 'kappa' above). Does
    // kappa need to be an input parameter again?
    // dt = fminf(c_x, h_dx/10.f);
    
    dt = min(c_x, c_y);

	dt = (t == 0.f) ? dt = 0.000001f : dt;	
	t += dt;
    // std::cout << "max_x = " << max_x
    //           << ", at mxresptr = " << mxresptr - mxptr
    //           << std::endl;
	// std::cout << "time step is " << dt << std::endl; 
}

void Simulator::PrintData(void) {
	dev_wet_count = 0;
	if (h_print || check_volume) {
		checkCudaErrors(cudaMemcpy2D(h_h, grid_config.h_nx*sizeof(double), dev_h, pitch,
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
	}
	if (save_max) {
		PrintSummaryData();
	}
	if (h_print) {
		std::stringstream filename_h;
		std::cout<<"Interpolated Flow Rate: "<< hydrograph.interpolated_rate << std:: endl;
		filename_h << output_file << "/h" << count_print << ".txt";
        writeGrid(filename_h.str(), h_h, grid_config);
	}

    if (drain_averaged) {
        std::cout << "Drain Rate: " << drain.interpolated_rate << std::endl;
    }
    
	if (check_volume) {
		V_computed = 0.f;

		if (infiltration) {
			checkCudaErrors(cudaMemcpy2D(h_F, grid_config.h_nx*sizeof(double), dev_F, pitch,
										 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
		}
	}

	if (check_volume) {
		
		#pragma omp parallel for reduction(+:V_computed)
		
		for (int j = grid_config.h_ny - 3; j >= 2; j--) {
			for (int i = 2; i < grid_config.h_nx-2; i++) {
				int   id       = j*grid_config.h_nx + i;
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
		checkCudaErrors(cudaMemcpy2D(h_q, grid_config.h_nx*sizeof(double), dev_q, pitch,
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
		std::stringstream filename_q;
		filename_q << output_file << "/q" << count_print << ".txt";
        writeGrid(filename_q.str(), h_q, grid_config);
	}

    count_print++;
}

void Simulator::PrintSummaryData(void) {
	if (save_max) {
        std::unique_ptr< double[] >
            h_h_max(new double[grid_config.h_nx*grid_config.h_ny]),
            h_q_max(new double[grid_config.h_nx*grid_config.h_ny]),
            h_t_peak(new double[grid_config.h_nx*grid_config.h_ny]),
            h_t_dry(new double[grid_config.h_nx*grid_config.h_ny]);
        
		checkCudaErrors(cudaMemcpy2D(&h_h_max[0],  grid_config.h_nx*sizeof(double), dev_h_max, pitch,
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
		checkCudaErrors(cudaMemcpy2D(&h_q_max[0],  grid_config.h_nx*sizeof(double), dev_q_max,  pitch,
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
		checkCudaErrors(cudaMemcpy2D(&h_t_peak[0], grid_config.h_nx * sizeof(double), dev_time_peak, pitch,
									 grid_config.h_nx * sizeof(double), grid_config.h_ny, DtoH));
		checkCudaErrors(cudaMemcpy2D(&h_t_dry[0], grid_config.h_nx * sizeof(double), dev_time_dry, pitch,
									 grid_config.h_nx * sizeof(double), grid_config.h_ny, DtoH));

		std::stringstream out_name;

        out_name.str(std::string());
		out_name << output_file << "/peak_flood_depth.txt";
        writeGrid(out_name.str(), &h_h_max[0], grid_config);
        
        out_name.str(std::string());
		out_name << output_file << "/peak_unit_flow.txt";
        writeGrid(out_name.str(), &h_q_max[0], grid_config);
        
        out_name.str(std::string());
		out_name << output_file << "/time_to_peak.txt";
        writeGrid(out_name.str(), &h_t_peak[0], grid_config);
        
        out_name.str(std::string());
		out_name << output_file << "/time_to_dry.txt";
        writeGrid(out_name.str(), &h_t_dry[0], grid_config);

	}

	if (save_arrival_time) {
        std::unique_ptr< double[] > h_t_wet(new double[grid_config.h_nx*grid_config.h_ny]);
		checkCudaErrors(cudaMemcpy2D(&h_t_wet[0],  grid_config.h_nx*sizeof(double), dev_t_wet,  pitch,
									 grid_config.h_nx*sizeof(double), grid_config.h_ny, DtoH));
		std::ofstream twet;
		std::stringstream arrival_name;
		arrival_name << output_file << "/flood_wave_arrival.txt";

        writeGrid(arrival_name.str(), &h_t_wet[0], grid_config);
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
    ApplyBoundaries(dev_w, dev_hu, dev_hv, dev_BC);
	count = 0;
 StartTimer();

}
void Simulator::CloseSimulation(){
	if (save_max) {
		PrintSummaryData();
	}

	FreeGrid(dev_w, dev_hu, dev_hv, dev_w_old, dev_hu_old, dev_hv_old, dev_dw, dev_dhu, dev_dhv, dev_mx, dev_my, dev_BC, dev_BX, dev_BY,
	         dev_wet_blocks, dev_active_blocks, dev_n, dev_hyetograph_gridded_rate, dev_F, dev_F_old,
             dev_dF, dev_K, dev_h, dev_q, dev_h_max, dev_q_max, dev_t_wet, dev_time_peak, dev_time_dry, dev_G);		//added t_peak and t_dry by Youcan on 20170908
	if(f_b) 			free(b);
	if(f_h_BX) 			free(h_BX);
	if(f_h_BY) 			free(h_BY);
	if(f_h_h) 			free(h_h);
	if(f_h_q)			free(h_q);
    if(f_h_hyetograph) 	free(h_hyetograph);
	if(f_h_n)			free(h_n);
	if(f_h_K)			free(h_K);
	if(f_h_F)			free(h_F);
	if(f_h_o)           free(h_o);
    if(f_hyetograph_o)	free(hyetograph_o);
	if(f_n_o)			free(n_o);
	if(f_K_o)			free(K_o);
	if(f_h_h_max) 		free(h_h_max);
	if(f_h_q_max) 		free(h_q_max);
	if(f_h_t_wet)		free(h_t_wet);
	free(B->data);
	/** Added for 1D-2D*/
	if(f_source_idx) {
        free(source_idx);
        cudaFree(source_idx_dev);
    }
	if(f_source_rate) {
        free(source_rate);
        cudaFree(source_rate_dev);
    }
}

double Simulator::RunSimulation() {
   
	
	// tf; 
   // while (t < tf) {
	//updateSources();

		if (dambreak || rainfall_averaged || drain_averaged) {
            checkCudaErrors(cudaMemcpy(source_idx_dev,source_idx,NumSources*sizeof(int),HtoD));
			//std::cout << "updating source" << std::endl; 
			UpdateSource();
		}
		//std::cout << "Grow Blocks" << std::endl; 
		Grow(dev_wet_blocks, dev_active_blocks, dev_hyetograph_gridded_rate, rainfall_gridded, dev_surge_gridded_depth, surge_gridded);
		//std::cout << "Compute Fluxes" << std::endl; 
        ComputeFluxes(dev_w, dev_hu, dev_hv, dev_dw, dev_dhu, dev_dhv, dev_mx, dev_my, dev_BC, dev_BX, dev_BY, dev_G, dev_active_blocks,
		              dt, dev_n, hydrograph.interpolated_rate, dambreak_source_idx,
		              hyetograph.interpolated_rate, drain.interpolated_rate,
                      dev_hyetograph_gridded_rate,
                      dev_drain_gridded_rate,
                      dev_F,
					  dev_F_old, dev_dF, dev_K,source_idx_dev,source_rate_dev,NumSources);
		ComputeTimestep();

		Integrate_1(dev_w, dev_hu, dev_hv, dev_w_old, dev_hu_old, dev_hv_old, dev_dw, dev_dhu, dev_dhv, dev_BC,
					dev_G, dev_wet_blocks, dev_active_blocks, t, dt,
                    hydrograph.interpolated_rate, dambreak_source_idx,
                    hyetograph.interpolated_rate, dev_hyetograph_gridded_rate,
                    dev_surge_gridded_depth, dev_F,
                    dev_F_old, dev_dF, dev_K, dev_h, dev_q, dev_h_max, dev_q_max, dev_t_wet,source_idx_dev, source_rate_dev,NumSources, dev_time_peak, dev_time_dry);	//added time_peak and time_dry by Youcan on 20170908
		ApplyBoundaries(dev_w, dev_hu, dev_hv, dev_BC);

		if (!euler_integration) {
			
			ComputeFluxes(dev_w, dev_hu, dev_hv, dev_dw, dev_dhu, dev_dhv, dev_mx, dev_my, dev_BC, dev_BX, dev_BY, dev_G, dev_active_blocks,
						  dt, dev_n, hydrograph.interpolated_rate, dambreak_source_idx,
						  hyetograph.interpolated_rate, drain.interpolated_rate,
                          dev_hyetograph_gridded_rate, dev_drain_gridded_rate,
						  dev_F, dev_F_old, dev_dF, dev_K,source_idx_dev, source_rate_dev,NumSources);
			Integrate_2(dev_w, dev_hu, dev_hv, dev_w_old, dev_hu_old, dev_hv_old, dev_dw, dev_dhu, dev_dhv, dev_BC,
						dev_G, dev_wet_blocks, dev_active_blocks, t, dt,
						hydrograph.interpolated_rate, dambreak_source_idx,
						hyetograph.interpolated_rate, dev_hyetograph_gridded_rate, dev_F,
						dev_F_old, dev_dF, dev_K, dev_h, dev_q, dev_h_max, dev_q_max, dev_t_wet,source_idx_dev, source_rate_dev,NumSources);
			ApplyBoundaries(dev_w, dev_hu, dev_hv, dev_BC);
		}

        if (rainfall_gridded) {
            hyetograph_series->update(t);
            V_added += hyetograph_series->sum()*dt;
        }

        if (drain_gridded) {
            drain_series->update(t);
        }

        if (surge_gridded) {
            surge_series->update(t);
        }

		if (h_print || q_print || check_volume) {
			if (t > t_print) {
				PrintData();
				t_print += dt_print;
				std::cout << t << "\t" << dt;
				double end_time = EndTimer();
				long int computed_cells = (long int)(grid_config.h_nx-3)*(long int)(grid_config.h_ny-3)*(long int)count;
				
				if (check_volume) {
					std::cout << "\t" << V_added*grid_config.h_dx*grid_config.h_dy << "\t" << V_computed*grid_config.h_dx*grid_config.h_dy <<
					"\t" << (V_computed - V_added) / V_added << "\t" << grid_config.b_nx*grid_config.b_ny << "\t"  << end_time << "\t" << (double)computed_cells / end_time /1000000.0 << "\t" << dev_wet_count;
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
	f_source_idx = 1;
	source_rate = (double*)malloc(numSources*sizeof(double));
	f_source_rate = 1;
	memset(source_idx,0,numSources*sizeof(int));
	memset(source_rate,0,numSources*sizeof(double));
	checkCudaErrors(cudaMalloc((void**)&source_idx_dev,numSources*sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&source_rate_dev,numSources*sizeof(double)));
	
}

void Simulator::SetSourceLocation(int i,double X, double Y){
	source_x = (int)((X-grid_config.h_xll) / grid_config.cellsize_original + 0.5);
	std::cout << "source_x is  " << source_x << std::endl; 
	std::cout << "X is  " << X << std::endl;
	std::cout << "h_xll is  " << grid_config.h_xll << std::endl;
	std::cout << "cellsize_original is  " << grid_config.cellsize_original << std::endl;
	source_y = (int)((Y-grid_config.h_yll) / grid_config.cellsize_original + 0.5);
	std::cout << "Y location is  " << source_y << std::endl; 
	int idx = source_y*grid_config.h_nx + source_x;
	std::cout << "hnx is  " << grid_config.h_nx << " idx is " << idx << std::endl;
	source_idx[i] =  idx;
		
}

void Simulator::setSourceValue (int i, double value){
	

	source_rate[i] = value/grid_config.h_dx/grid_config.h_dy;
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
