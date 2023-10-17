#ifndef GRID_H
#define GRID_H
#include "io.h"
// #include "simulator.h"

void SetDeviceConstants(int num_columns, int num_rows,
                        double cellxsize, double cellysize,
                        double h_kappa, double h_nodata);

void InterpGrid(double factor, double *x0_dev, double *x1_dev, double *x_dev);
void FillGrid(double *x, double value);

void AllocateGrid(double *&w, double *&hu, double *&hv, double *&w_old,
                  double *&hu_old, double *&hv_old, double *&dw, double *&dhu,
                  double *&dhv, double *&mx, double *&my, double *&BC, double *&BX, double *&BY,
                  bool *&wet_blocks, int *&active_blocks, double *&n,
                  double *&hyetograph_gridded_rate, double *&F, double *&F_old,
                  double *&dF, double *&K, double *&h, double *&q, double *&h_max,
                  double *&q_max, double *&t_wet, bool h_dambreak,
                  bool h_rainfall_averaged, bool h_rainfall_gridded, bool h_infiltration,
                  bool h_surge_gridded, bool h_euler_integration,
                  bool h_check_volume, bool h_h_init, bool h_h_print,
                  bool h_q_print, bool h_save_max, bool h_save_arrival_time,
                  double h_psi, double h_dtheta, double *&t_peak, double *&t_dry, double *&G, 
                  GridConfig& grid_config);	//added time_peak and time_dry by Youcan on 20170908

void InitGrid(double *w, double *hu, double *hv, double *w_old, double *hu_old,
              double *hv_old, double *BC, double *BX, double *BY, bool *wet_blocks,
              int *active_blocks, double *h, double *t_wet, double *G);

void ComputeFluxes(double *w, double *hu, double *hv, double *dw, double *dhu,
                   double *dhv, double *mx, double *my, double *BC, double *BX, double *BY,
                   double *G, int *active_blocks, double dt, double *n,
                   double hydrograph_rate, int hydrograph_source,
                   double hyetograph_rate, double *hyetograph_gridded_rate,
                   double *F, double *F_old, double *dF, double *K,int *source_idx_dev, double *source_rate_dev, long numSources);

void Integrate_1(double *w, double *hu, double *hv, double *w_old, double *hu_old,
                 double *hv_old, double *dw, double *dhu, double *dhv, double *BC,
                 double *G, bool *wet_blocks, int *active_blocks, double t, double dt,
                 double hydrograph_rate, int hydrograph_source,
                 double hyetograph_rate, double *hyetograph_gridded_rate,
                 double *surge_gridded_depth, 
                 double *F, double *F_old, double *dF, double *K, double *h,
                 double *q, double *h_max, double *q_max, double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources, double *t_peak, double *t_dry);	//added time_peak and time_dry by Youcan on 20170908

void Integrate_2(double *w, double *hu, double *hv, double *w_old, double *hu_old,
                 double *hv_old, double *dw, double *dhu, double *dhv, double *BC,
                 double *G, bool *wet_blocks, int *active_blocks, double t, double dt,
                 double hydrograph_rate, int hydrograph_source,
                 double hyetograph_rate, double *hyetograph_gridded_rate,
                 double *F, double *F_old, double *dF, double *K, double *h,
                 double *q, double *h_max, double *q_max, double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources);

void Grow(bool *wet_blocks, int *active_blocks,
          double *hyetograph_gridded_rate, bool h_rainfall_gridded,
          double *surge_gridded_depth, bool h_surge_gridded);

void ApplyBoundaries(double *w, double *hu, double *hv, double *BC);

void FreeGrid(double *&w, double *&hu, double *&hv, double *&w_old, double *&hu_old,
              double *&hv_old, double *&dw, double *&dhu, double *&dhv, double *&mx, double *&my,
              double *&BC, double *&BX, double *&BY, bool *&wet_blocks,
              int *&active_blocks, double *&n, double *&hyetograph_gridded_rate,
              double *&F, double *&F_old, double *&dF, double *&K, double *&h,
              double *&q, double *&h_max, double *&q_max, double *&t_wet, double *&t_peak, 
	      double *&t_dry, double *&G);	//added t_peak and t_dry by Youcan on 20170908

#endif
