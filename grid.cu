#include <stdio.h>
#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <thrust/device_ptr.h>
#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/copy.h>
#include <thrust/remove.h>
#include <thrust/count.h>

#include "constants.h"
#include "grid.h"

#define BLOCK_ROWS 14
#define BLOCK_COLS 16

__constant__ bool  dambreak, rainfall_averaged, rainfall_gridded, infiltration;
__constant__ bool  euler_integration, check_volume;
__constant__ bool  h_init, h_print, q_print, save_max, save_arrival_time;
__constant__ int   nx, ny, nBlocksX, nBlocksY;
__constant__ double dx, dy, kappa, psi, dtheta;
__constant__ double g       = 9.80665f;
__constant__ double epsilon = 1.19209e-07f;
__constant__ double theta   = 1.3f;

dim3 BlockDim,  GridDim;
dim3 BlockDimX, GridDimX;
dim3 BlockDimY, GridDimY;
int nBlocks;
int GridSize;
dim3 GrowGrid;

size_t pitch, pitchBX, pitchBY;

struct is_nonnegative {
    __host__ __device__
    bool operator()(const int x) {
        return (x > -1) == true;
    }
};

struct is_negative {
    __host__ __device__
    bool operator()(const int x) {
        return (x < 0) == true;
    }
};

template <typename T>
__device__ inline T* getElement(T *base, int Pitch, int row,
                                int col) {
	return (T*)((char*)base + row*Pitch) + col;
}

__device__ double signf(double f) {
    if (((int&)f & 0x7FFFFFFF)==0) {
        return 0.f; // test exponent & mantissa bits: is input zero?
    } else {
		// mask sign bit in f, set it in r if necessary
        double r = 1.f;
        (int&)r |= ((int&)f & 0x80000000);
        return r;
    }
}

__device__ double minmod(double a, double b, double c) {
    return 0.25f*signf(a)*(signf(a)+signf(b))*(signf(b)+signf(c))*
           fminf(fminf(fabsf(a), fabsf(b)), fabsf(c));
}

__device__ double limiter(double U, double Up1, double Um1) {
    return 0.5f * minmod(theta*(Up1-U), 0.5f*(Up1-Um1), theta*(U-Um1));
}

__device__ double getVelocity(double h, double hz) {
    return 2.f*h*hz / (h*h + fmaxf(h*h, kappa));
}

__device__ double computeH(double Uij, double Fij, double Ukh, double Fkh, double cm,
                          double cp) {
    return (cp*Fij-cm*Fkh)/(cp-cm) + (cp*cm)/(cp-cm)*(Ukh-Uij);
}

__device__ void Construct(double w,   double hu,   double hv,
                          double wp1, double hup1, double hvp1,
                          double wm1, double hum1, double hvm1,
                          double &wp, double &hup, double &hvp,
                          double &wm, double &hum, double &hvm,
                          double Bp,  double Bm) {
	double wz, huz, hvz;

	wz  = limiter(w,  wp1,  wm1);
    huz = limiter(hu, hup1, hum1);
    hvz = limiter(hv, hvp1, hvm1);

	wp  = w  + wz,  wm  = w  - wz;
	hup = hu + huz, hum = hu - huz;
    hvp = hv + hvz, hvm = hv - hvz;

	double wBar = wp + wm;
	if (wp < Bp) {
		wp = Bp;
		wm = wBar - Bp;
	}
	if (wm < Bm) {
		wp = wBar - Bm;
		wm = Bm;
	}
}

__device__ void getHX(double wp,   double hup,   double hvp,
                      double wp1m, double hup1m, double hvp1m,
                      double &dwp, double &dhup, double &dhvp,
                      double &mx,  double Bp) {
	double hp   = wp   - Bp;
	double hp1m = wp1m - Bp;
    double up   = getVelocity(hp,   hup);
    double up1m = getVelocity(hp1m, hup1m);
    double vp   = getVelocity(hp,   hvp);
    double vp1m = getVelocity(hp1m, hvp1m);

	hup   = hp   * up;
	hup1m = hp1m * up1m;
	hvp   = hp   * vp;
	hvp1m = hp1m * vp1m;

	double ap = fmaxf(fmaxf(up + sqrtf(g*hp), up1m + sqrtf(g*hp1m)), 0.f);
	double am = fminf(fminf(up - sqrtf(g*hp), up1m - sqrtf(g*hp1m)), 0.f);

    double Fp0,   Fp1,   Fp2;
    double Fp1m0, Fp1m1, Fp1m2;

    if (hp > 0.f) {
		Fp0 = hup;
		Fp1 = hup*hup / hp + 0.5f*g*hp*hp;
		Fp2 = hup*hvp / hp;
    } else {
		Fp0 = Fp1 = Fp2 = 0.f;
    }

    if (hp1m > 0.f) {
		Fp1m0 = hup1m;
		Fp1m1 = hup1m*hup1m / hp1m + 0.5f*g*hp1m*hp1m;
		Fp1m2 = hup1m*hvp1m / hp1m;
    } else {
		Fp1m0 = Fp1m1 = Fp1m2 = 0.f;
    }

	if (ap-am > 0.f) {
		dwp  = computeH(wp,  Fp0, wp1m,  Fp1m0, am, ap);
		dhup = computeH(hup, Fp1, hup1m, Fp1m1, am, ap);
		dhvp = computeH(hvp, Fp2, hvp1m, Fp1m2, am, ap);
	} else {
		dwp = dhup = dhvp = 0.f;
    }

    mx = fmaxf(ap, -am);
}

__device__ void getHY(double wp,   double hup,   double hvp,
                      double wp1m, double hup1m, double hvp1m,
                      double &dwp, double &dhup, double &dhvp,
                      double &mx,  double Bp) {
	double hp   = wp   - Bp;
	double hp1m = wp1m - Bp;
    double up   = getVelocity(hp,   hup);
    double up1m = getVelocity(hp1m, hup1m);
    double vp   = getVelocity(hp,   hvp);
    double vp1m = getVelocity(hp1m, hvp1m);

	hup   = hp   * up;
	hup1m = hp1m * up1m;
	hvp   = hp   * vp;
	hvp1m = hp1m * vp1m;

	double bp = fmaxf(fmaxf(vp + sqrtf(g*hp), vp1m + sqrtf(g*hp1m)), 0.f);
	double bm = fminf(fminf(vp - sqrtf(g*hp), vp1m - sqrtf(g*hp1m)), 0.f);

    double Gp0,   Gp1,   Gp2;
    double Gp1m0, Gp1m1, Gp1m2;

    if (hp > 0.f) {
        Gp0 = hvp;
        Gp1 = hup*hvp / hp;
        Gp2 = hvp*hvp / hp + 0.5f*g*hp*hp;
    } else {
        Gp0 = Gp1 = Gp2 = 0.f;
    }

    if (hp1m > 0.f) {
        Gp1m0 = hvp1m;
        Gp1m1 = hup1m*hvp1m / hp1m;
        Gp1m2 = hvp1m*hvp1m / hp1m + 0.5f*g*hp1m*hp1m;
    } else {
        Gp1m0 = Gp1m1 = Gp1m2 = 0.f;
    }

	if (bp-bm > 0.f) {
		dwp  = computeH(wp,  Gp0, wp1m,  Gp1m0, bm, bp);
		dhup = computeH(hup, Gp1, hup1m, Gp1m1, bm, bp);
		dhvp = computeH(hvp, Gp2, hvp1m, Gp1m2, bm, bp);
	} else {
        dwp = dhup = dhvp = 0.f;
    }

    mx = fmaxf(mx, fmaxf(bp, -bm));
}

void SetDeviceConstants(int num_columns, int num_rows, double cellsize,
                        double h_kappa) {
    double cellsize_meters = (double)cellsize*6378137.f*(double)pi / 180.f;
	int   grid_columns    = num_columns + 4 - 1;
	int   grid_rows       = num_rows    + 4 - 1;

    cudaMemcpyToSymbol(nx,    &grid_columns,    sizeof(int),   0, HtoD);
    cudaMemcpyToSymbol(ny,    &grid_rows,       sizeof(int),   0, HtoD);
    cudaMemcpyToSymbol(dx,    &cellsize_meters, sizeof(double), 0, HtoD);
    cudaMemcpyToSymbol(dy,    &cellsize_meters, sizeof(double), 0, HtoD);
    cudaMemcpyToSymbol(kappa, &h_kappa,         sizeof(double), 0, HtoD);

	BlockDim.x = BLOCK_COLS;
	BlockDim.y = BLOCK_ROWS;

	int num_blocks_x = ceil((double)grid_columns / (double)BlockDim.x);
	int num_blocks_y = ceil((double)grid_rows    / (double)BlockDim.y);

    GridDim.x = num_blocks_x;
    GridDim.y = num_blocks_y;

    cudaMemcpyToSymbol(nBlocksX, &num_blocks_x, sizeof(int), 0, HtoD);
    cudaMemcpyToSymbol(nBlocksY, &num_blocks_y, sizeof(int), 0, HtoD);

	GridSize = GridDim.x * GridDim.y;
	nBlocks  = GridDim.x * GridDim.y;
}

void AllocateGrid(double *&w, double *&hu, double *&hv, double *&w_old,
                  double *&hu_old, double *&hv_old, double *&dw, double *&dhu,
                  double *&dhv, double *&mx, double *&BC, double *&BX, double *&BY,
                  bool *&wet_blocks, int *&active_blocks, double *&n,
                  double *&hyetograph_gridded_rate, double *&F, double *&F_old,
                  double *&dF, double *&K, double *&h, double *&q, double *&h_max,
                  double *&q_max, double *&t_wet, bool h_dambreak,
                  bool h_rainfall_averaged, bool h_rainfall_gridded, bool h_infiltration,
                  bool h_euler_integration,
                  bool h_check_volume, bool h_h_init, bool h_h_print,
                  bool h_q_print, bool h_save_max, bool h_save_arrival_time,
                  double h_psi, double h_dtheta, double *&t_peak, double *&t_dry, double *&G) {	//added time_peak and time_dry by Youcan on 20170908
//10000196001000
//100000001000


	std::cout << h_dambreak << h_rainfall_averaged << h_rainfall_gridded <<
h_infiltration << h_euler_integration << h_check_volume <<
h_h_init << h_h_print << h_q_print << h_save_max << h_save_arrival_time <<
std::endl;

    cudaMemcpyToSymbol(dambreak,          &h_dambreak,          sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(rainfall_averaged, &h_rainfall_averaged, sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(rainfall_gridded,  &h_rainfall_gridded,  sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(infiltration,      &h_infiltration,      sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(euler_integration, &h_euler_integration, sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(check_volume,      &h_check_volume,      sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(h_init,            &h_h_init,            sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(h_print,           &h_h_print,           sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(q_print,           &h_q_print,           sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(save_max,          &h_save_max,          sizeof(bool), 0, HtoD);
    cudaMemcpyToSymbol(save_arrival_time, &h_save_arrival_time, sizeof(bool), 0, HtoD);

	size_t width     = (GridDim.x * BLOCK_COLS)     * sizeof(double);
	size_t width_p1  = (GridDim.x * BLOCK_COLS + 1) * sizeof(double);
	size_t height    = (GridDim.y * BLOCK_ROWS);
	size_t height_p1 = (GridDim.y * BLOCK_ROWS + 1);

    checkCudaErrors(cudaMallocPitch((void**)&w,      &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&hu,     &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&hv,     &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&w_old,  &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&hu_old, &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&hv_old, &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&G,      &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&dw,     &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&dhu,    &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&dhv,    &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&BC,     &pitch,   width, height ));
    checkCudaErrors(cudaMallocPitch((void**)&BX, &pitchBX, width_p1, height  ));
    checkCudaErrors(cudaMallocPitch((void**)&BY, &pitchBY, width,   height_p1));

    checkCudaErrors(cudaMemset2D(w,      pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(hu,     pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(hv,     pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(w_old,  pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(hu_old, pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(hv_old, pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(G,      pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(dw,     pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(dhu,    pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(dhv,    pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(BC,     pitch,   0, width,    height   ));
    checkCudaErrors(cudaMemset2D(BX,     pitchBX, 0, width_p1, height   ));
    checkCudaErrors(cudaMemset2D(BY,     pitchBY, 0, width,    height_p1));

	if (h_h_print) {
		checkCudaErrors(cudaMallocPitch((void**)&h, &pitch, width, height));
		checkCudaErrors(cudaMemset2D(h, pitch, 0, width, height));
	}

	if (h_q_print) {
		checkCudaErrors(cudaMallocPitch((void**)&q, &pitch, width, height));
		checkCudaErrors(cudaMemset2D(q, pitch, 0, width, height));
	}

	if (h_rainfall_gridded) {
		checkCudaErrors(cudaMallocPitch((void**)&hyetograph_gridded_rate, &pitch,
										width, height));
		checkCudaErrors(cudaMemset2D(hyetograph_gridded_rate, pitch, 0, width,
									 height));
	}

	if (h_infiltration) {
		checkCudaErrors(cudaMallocPitch((void**)&F,     &pitch, width, height));
		checkCudaErrors(cudaMallocPitch((void**)&F_old, &pitch, width, height));
		checkCudaErrors(cudaMallocPitch((void**)&dF,    &pitch, width, height));
		checkCudaErrors(cudaMallocPitch((void**)&K,     &pitch, width, height));

		checkCudaErrors(cudaMemset2D(F,     pitch, 0, width, height));
		checkCudaErrors(cudaMemset2D(F_old, pitch, 0, width, height));
		checkCudaErrors(cudaMemset2D(dF,    pitch, 0, width, height));
		checkCudaErrors(cudaMemset2D(K,     pitch, 0, width, height));

		checkCudaErrors(cudaMemcpyToSymbol(psi,    &h_psi,    sizeof(double),
										   0, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpyToSymbol(dtheta, &h_dtheta, sizeof(double),
										   0, cudaMemcpyHostToDevice));
	}

    checkCudaErrors(cudaMallocPitch((void**)&n, &pitch, width, height));
    checkCudaErrors(cudaMemset2D(n, pitch, 0, width, height));

	if (h_save_max) {
		checkCudaErrors(cudaMallocPitch((void**)&h_max, &pitch, width, height));
		checkCudaErrors(cudaMallocPitch((void**)&q_max, &pitch, width, height));
		checkCudaErrors(cudaMallocPitch((void**)&t_peak, &pitch, width, height));				// added by Youcan on 20170424
		checkCudaErrors(cudaMallocPitch((void**)&t_dry, &pitch, width, height));				// added by Youcan on 20170830
		checkCudaErrors(cudaMemset2D(h_max, pitch, 0, width, height));
		checkCudaErrors(cudaMemset2D(q_max, pitch, 0, width, height));
		checkCudaErrors(cudaMemset2D(t_peak, pitch, 0, width, height));							// added by Youcan on 20170424
		checkCudaErrors(cudaMemset2D(t_dry, pitch, 0, width, height));							// added by Youcan on 20170830
	}

	if (h_save_arrival_time) {
		checkCudaErrors(cudaMallocPitch((void**)&t_wet, &pitch, width, height));
		checkCudaErrors(cudaMemset2D(t_wet, pitch, 0, width, height));
	}

	BlockDimX.x = 48, BlockDimY.x = 48;
    GridDimX.x  = ceil((double)h_ny / (double)BlockDimX.x);
    GridDimY.x  = ceil((double)h_nx / (double)BlockDimY.x);

    GrowGrid.x = ceil((double)GridDim.x / (double)BlockDim.x);
    GrowGrid.y = ceil((double)GridDim.y / (double)BlockDim.y);

    size_t GrowGridSize = (GrowGrid.x*BLOCK_COLS)*(GrowGrid.y*BLOCK_ROWS);

    checkCudaErrors(cudaMalloc((void**)&wet_blocks,    GrowGridSize*sizeof(bool) ));
    checkCudaErrors(cudaMalloc((void**)&active_blocks, GrowGridSize*sizeof(int)  ));
    checkCudaErrors(cudaMalloc((void**)&mx,            GrowGridSize*sizeof(double)));

    checkCudaErrors(cudaMemset(wet_blocks,    (bool) false, GridSize*sizeof(bool) ));
    checkCudaErrors(cudaMemset(active_blocks, (int)  -1,    GridSize*sizeof(int)  ));
    checkCudaErrors(cudaMemset(mx,            (double) 0,    GridSize*sizeof(double)));
}

__global__ void InitGrid_k(double *w, double *hu, double  *hv, double *w_old,
                           double *hu_old, double *hv_old, double *BC, double *BX,
                           double *BY, bool *wet_blocks, int *active_blocks,
                           size_t pitch, size_t pitchBX, size_t pitchBY,
                           double *h, double *t_wet, double *G) {
    int i = blockIdx.x*BLOCK_COLS + threadIdx.x + 2;
    int j = blockIdx.y*BLOCK_ROWS + threadIdx.y + 2;

    if (i >= nx || j >= ny) return;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
        int blockID = blockIdx.y*nBlocksX + blockIdx.x;
        wet_blocks   [blockID] = true;
        active_blocks[blockID] = blockID;
    }

	double *w_ij      = getElement(w,      pitch, j, i);
	double *hu_ij     = getElement(hu,     pitch, j, i);
	double *hv_ij     = getElement(hv,     pitch, j, i);
	double *w_old_ij  = getElement(w_old,  pitch, j, i);
	double *hu_old_ij = getElement(hu_old, pitch, j, i);
	double *hv_old_ij = getElement(hv_old, pitch, j, i);
	double *G_ij	  = getElement(G,      pitch, j, i);


	double *BC_ij     = getElement(BC, pitch,   j,   i  );
	double *BE_ij     = getElement(BX, pitchBX, j,   i+1);
	double *BW_ij     = getElement(BX, pitchBX, j,   i  );
	double *BN_ij     = getElement(BY, pitchBY, j+1, i  );
	double *BS_ij     = getElement(BY, pitchBY, j,   i  );

	*BC_ij = 0.25f * ((*BE_ij)+(*BW_ij)+(*BN_ij)+(*BS_ij));
	*w_ij  = *w_old_ij  = *BC_ij;
	*hu_ij = *hu_old_ij = 0.f;
	*hv_ij = *hv_old_ij = 0.f;
    *G_ij = 0.0f;

	if (h_init) {
		double *h_ij = getElement(h, pitch, j, i);
		*w_ij     += *h_ij;
		*w_old_ij += *h_ij;
	}

	if (save_arrival_time) {
		double *t_wet_ij = getElement(t_wet, pitch, j, i);
		*t_wet_ij = -9999.f;
	}
}


void InitGrid(double *w, double *hu, double *hv, double *w_old, double *hu_old,
              double *hv_old, double *BC, double *BX, double *BY, bool *wet_blocks,
              int *active_blocks, double *h, double *t_wet, double *G) {
	InitGrid_k <<< GridDim, BlockDim >>> (w, hu, hv, w_old, hu_old, hv_old, BC,
	                                      BX, BY, wet_blocks, active_blocks,
	                                      pitch, pitchBX, pitchBY, h, t_wet, G);
}

__global__ void ComputeFluxes_k(double *w, double *hu, double *hv, double *dw,
                                double *dhu, double *dhv, double *mx, double *BC,
                                double *BX, double *BY, double *G, int *active_blocks,
                                double dt, size_t pitch, size_t pitchBX,
                                size_t pitchBY, double *n, double hydrograph_rate,
                                int hydrograph_source, double hyetograph_rate,
                                double *hyetograph_gridded_rate, double *F,
                                double *F_old, double *dF, double *K,int *source_idx_dev, double *source_rate_dev, long numSources) {
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;

    int blockX = active_blocks[blockIdx.x] % nBlocksX;
    int blockY = active_blocks[blockIdx.x] / nBlocksX;

	int i = blockX*BLOCK_COLS + tidx + 2;
	int j = blockY*BLOCK_ROWS + tidy + 2;

 	__shared__ double sw  [BLOCK_ROWS+4][BLOCK_COLS+4];
 	__shared__ double shu [BLOCK_ROWS+4][BLOCK_COLS+4];
 	__shared__ double shv [BLOCK_ROWS+4][BLOCK_COLS+4];
	__shared__ double sBZ [BLOCK_ROWS+4][BLOCK_COLS+4];
 	__shared__ double swp [BLOCK_ROWS+2][BLOCK_COLS+2];
 	__shared__ double shup[BLOCK_ROWS+2][BLOCK_COLS+2];
 	__shared__ double shvp[BLOCK_ROWS+2][BLOCK_COLS+2];
 	__shared__ double swm [BLOCK_ROWS+2][BLOCK_COLS+2];
 	__shared__ double shum[BLOCK_ROWS+2][BLOCK_COLS+2];
 	__shared__ double shvm[BLOCK_ROWS+2][BLOCK_COLS+2];
	__shared__ double dwp [BLOCK_ROWS+1][BLOCK_COLS+1];
	__shared__ double dhup[BLOCK_ROWS+1][BLOCK_COLS+1];
	__shared__ double dhvp[BLOCK_ROWS+1][BLOCK_COLS+1];
    __shared__ double mxsum[BLOCK_ROWS][BLOCK_COLS];

    __syncthreads();


    // Don't go outside allocated arrays
    if (i >= nBlocksX*BLOCK_COLS || j >= nBlocksY*BLOCK_ROWS) return;

	double *w_ij   = getElement(w,   pitch, j, i);
	double *hu_ij  = getElement(hu,  pitch, j, i);
	double *hv_ij  = getElement(hv,  pitch, j, i);
	double *dw_ij  = getElement(dw,  pitch, j, i);
	double *dhu_ij = getElement(dhu, pitch, j, i);
	double *dhv_ij = getElement(dhv, pitch, j, i);
	double *BC_ij  = getElement(BC,  pitch, j, i);
	double *G_ij   = getElement(G,   pitch, j, i);

	// Fill the inner shared memory cells
	sw [tidy+2][tidx+2] = *w_ij;
	shu[tidy+2][tidx+2] = *hu_ij;
	shv[tidy+2][tidx+2] = *hv_ij;
	sBZ[tidy+2][tidx+2] = *getElement(BX, pitchBX, j, i);

	if (tidy*BLOCK_COLS + tidx < 2*BLOCK_COLS) {
		int jt = blockY*BLOCK_ROWS + tidy;

		sw [tidy  ][tidx+2] = *getElement(w,  pitch,   jt, i);
		shu[tidy  ][tidx+2] = *getElement(hu, pitch,   jt, i);
		shv[tidy  ][tidx+2] = *getElement(hv, pitch,   jt, i);
		sBZ[tidy  ][tidx+2] = *getElement(BY, pitchBY, jt, i);

		jt = blockY*BLOCK_ROWS + tidy + BLOCK_ROWS+2;
        if (i < nx && jt < ny) {
            sw [tidy+BLOCK_ROWS+2][tidx+2] = *getElement(w,  pitch,   jt, i);
            shu[tidy+BLOCK_ROWS+2][tidx+2] = *getElement(hu, pitch,   jt, i);
            shv[tidy+BLOCK_ROWS+2][tidx+2] = *getElement(hv, pitch,   jt, i);
            sBZ[tidy+BLOCK_ROWS+2][tidx+2] = *getElement(BY, pitchBY, jt, i);
        }

		if (tidx < BLOCK_ROWS) {
			int it = blockX*BLOCK_COLS + tidy;
			jt = blockY*BLOCK_ROWS + tidx + 2;
            if (it < nx && jt < ny) {
				sw [tidx+2][tidy  ] = *getElement(w,  pitch,   jt, it);
				shu[tidx+2][tidy  ] = *getElement(hu, pitch,   jt, it);
				shv[tidx+2][tidy  ] = *getElement(hv, pitch,   jt, it);
				sBZ[tidx+2][tidy  ] = *getElement(BX, pitchBX, jt, it);
            }

            it = blockX*BLOCK_COLS + tidy + BLOCK_COLS+2;
            if (it < nx && jt < ny) {
                sw [tidx+2][tidy+BLOCK_COLS+2] = *getElement(w,  pitch,   jt, it);
                shu[tidx+2][tidy+BLOCK_COLS+2] = *getElement(hu, pitch,   jt, it);
                shv[tidx+2][tidy+BLOCK_COLS+2] = *getElement(hv, pitch,   jt, it);
                sBZ[tidx+2][tidy+BLOCK_COLS+2] = *getElement(BX, pitchBX, jt, it);
            }
		}
	}

	__syncthreads();

    double junk = 0.f;

	Construct(sw [tidy+2][tidx+2], shu [tidy+2][tidx+2], shv [tidy+2][tidx+2],
	          sw [tidy+2][tidx+3], shu [tidy+2][tidx+3], shv [tidy+2][tidx+3],
	          sw [tidy+2][tidx+1], shu [tidy+2][tidx+1], shv [tidy+2][tidx+1],
	          swp[tidy+1][tidx+1], shup[tidy+1][tidx+1], shvp[tidy+1][tidx+1],
	          swm[tidy+1][tidx+1], shum[tidy+1][tidx+1], shvm[tidy+1][tidx+1],
	          sBZ[tidy+2][tidx+3], sBZ [tidy+2][tidx+2]);

	if (tidy*BLOCK_COLS + tidx < 2*BLOCK_COLS) {
		if (tidx < BLOCK_ROWS) {
			int tcol = (tidy == 0) ? 0 : BLOCK_COLS+1;

			Construct(sw [tidx+2][tcol+1], shu [tidx+2][tcol+1], shv [tidx+2][tcol+1],
					  sw [tidx+2][tcol+2], shu [tidx+2][tcol+2], shv [tidx+2][tcol+2],
					  sw [tidx+2][tcol  ], shu [tidx+2][tcol  ], shv [tidx+2][tcol  ],
					  swp[tidx+1][tcol  ], shup[tidx+1][tcol  ], shvp[tidx+1][tcol  ],
					  swm[tidx+1][tcol  ], shum[tidx+1][tcol  ], shvm[tidx+1][tcol  ],
					  sBZ[tidx+2][tcol+2], sBZ [tidx+2][tcol+1]);
		}
	}

	__syncthreads();

    double smx;

	getHX(swp[tidy+1][tidx+1], shup[tidy+1][tidx+1], shvp[tidy+1][tidx+1],
		  swm[tidy+1][tidx+2], shum[tidy+1][tidx+2], shvm[tidy+1][tidx+2],
		  dwp[tidy+1][tidx+1], dhup[tidy+1][tidx+1], dhvp[tidy+1][tidx+1],
		  smx,                 sBZ [tidy+2][tidx+3]);

	if (tidx < BLOCK_ROWS) {
		getHX(swp[tidx+1][0], shup[tidx+1][0], shvp[tidx+1][0],
			  swm[tidx+1][1], shum[tidx+1][1], shvm[tidx+1][1],
			  dwp[tidx+1][0], dhup[tidx+1][0], dhvp[tidx+1][0],
			  junk,           sBZ [tidx+2][2]);
	}
			
	__syncthreads();

	*dw_ij  = -(dwp [tidy+1][tidx+1] - dwp [tidy+1][tidx  ]) / dx;
	*dhu_ij = -(dhup[tidy+1][tidx+1] - dhup[tidy+1][tidx  ]) / dx;
	*dhv_ij = -(dhvp[tidy+1][tidx+1] - dhvp[tidy+1][tidx  ]) / dx;

	double S1 = (sBZ[tidy+2][tidx+3] - sBZ[tidy+2][tidx+2]) / dx;

	__syncthreads();

	sBZ[tidy+2][tidx+2] = *getElement(BY, pitchBY, j, i);

    __syncthreads();    

	Construct(sw [tidy+2][tidx+2], shu [tidy+2][tidx+2],  shv [tidy+2][tidx+2],
	          sw [tidy+3][tidx+2], shu [tidy+3][tidx+2],  shv [tidy+3][tidx+2],
	          sw [tidy+1][tidx+2], shu [tidy+1][tidx+2],  shv [tidy+1][tidx+2],
	          swp[tidy+1][tidx+1], shup[tidy+1][tidx+1],  shvp[tidy+1][tidx+1],
	          swm[tidy+1][tidx+1], shum[tidy+1][tidx+1],  shvm[tidy+1][tidx+1],
	          sBZ[tidy+3][tidx+2], sBZ [tidy+2][tidx+2]);

	if (tidy*BLOCK_COLS + tidx < 2*BLOCK_COLS) {
		int trow = (tidy == 0) ? 0 : BLOCK_ROWS+1;
		Construct(sw [trow+1][tidx+2], shu [trow+1][tidx+2], shv [trow+1][tidx+2],
				  sw [trow+2][tidx+2], shu [trow+2][tidx+2], shv [trow+2][tidx+2],
				  sw [trow+0][tidx+2], shu [trow+0][tidx+2], shv [trow+0][tidx+2],
				  swp[trow+0][tidx+1], shup[trow+0][tidx+1], shvp[trow+0][tidx+1],
				  swm[trow+0][tidx+1], shum[trow+0][tidx+1], shvm[trow+0][tidx+1],
				  sBZ[trow+2][tidx+2], sBZ [trow+1][tidx+2]);
	}

	__syncthreads();

	getHY(swp[tidy+1][tidx+1], shup[tidy+1][tidx+1], shvp[tidy+1][tidx+1],
		  swm[tidy+2][tidx+1], shum[tidy+2][tidx+1], shvm[tidy+2][tidx+1],
		  dwp[tidy+1][tidx+1], dhup[tidy+1][tidx+1], dhvp[tidy+1][tidx+1],
		  smx,                 sBZ [tidy+3][tidx+2]);

	if (tidx < BLOCK_COLS) {
		getHY(swp[0][tidx+1], shup[0][tidx+1], shvp[0][tidx+1],
			  swm[1][tidx+1], shum[1][tidx+1], shvm[1][tidx+1],
			  dwp[0][tidx+1], dhup[0][tidx+1], dhvp[0][tidx+1],
			  junk,           sBZ [2][tidx+2]);
	}

	__syncthreads();

	double S2 = (sBZ[tidy+3][tidx+2] - sBZ[tidy+2][tidx+2]) / dy;
	
    *dw_ij  += -(dwp [tidy+1][tidx+1] - dwp [tidy+0][tidx+1]) / dy;
    *dhu_ij += -(dhup[tidy+1][tidx+1] - dhup[tidy+0][tidx+1]) / dy;
    *dhv_ij += -(dhvp[tidy+1][tidx+1] - dhvp[tidy+0][tidx+1]) / dy;

	double hC  = (*w_ij) - (*BC_ij);
	*w_ij = (hC > epsilon) ? *w_ij : *BC_ij;

	S1 *= -g * hC;
	S2 *= -g * hC;

	{
		double *n_ij = getElement(n, pitch, j, i);
		double C = -4.f*g*(*n_ij)*(*n_ij)*sqrtf((*hu_ij)*(*hu_ij) +
											   (*hv_ij)*(*hv_ij)) *
				  powf(hC, 5.f/3.f) / powf((hC*hC + fmaxf(hC*hC, kappa)), 2.f);
		*G_ij = C;

		// S1 += C * (*hu_ij);
		// S2 += C * (*hv_ij);
	}

    double S0 = 0.f;

	if (rainfall_averaged) {
		S0 += hyetograph_rate;
	}

	if (rainfall_gridded) {
		double *hyetograph_gridded_rate_ij =
			getElement(hyetograph_gridded_rate, pitch, j, i);
		S0 += *hyetograph_gridded_rate_ij;
	}
	
	if (dambreak||numSources>0) {
		for (int counter = 0; counter < numSources; counter++){
		S0 = (j*nx + i == source_idx_dev[counter]) ? S0 + source_rate_dev[counter] : S0;
		//if(S0>0.0f) printf("\nSource %g",S0);
		}
		//S0 = (j*nx + i == hydrograph_source) ? S0 + hydrograph_rate : S0;
	}

	if (infiltration) {
		double *F_ij  = getElement(F,  pitch, j, i);
		double *dF_ij = getElement(dF, pitch, j, i);
		double *K_ij  = getElement(K,  pitch, j, i);

		double p1       = (*K_ij) * dt - 2.f * (*F_ij);
		double p2       = (*K_ij) * ((*F_ij) + psi*dtheta);
		double inf_rate = (p1 + sqrtf(p1*p1 + 8.f*p2*dt)) / (2.f*dt);

		inf_rate = (inf_rate > hC / dt) ? hC / dt : inf_rate;
		inf_rate = (inf_rate > 0.f)     ? inf_rate : 0.f;
		*dF_ij   = inf_rate;

		S0 -= inf_rate;
	}

    *dw_ij  += S0;
    *dhu_ij += S1;
    *dhv_ij += S2;

    mxsum[tidy][tidx] = smx;
    if (i > nx - 3 || j > ny - 3) {
        mxsum[tidy][tidx] = 0.f;
    }

    __syncthreads();
			
	if (tidx == 0) {
		for (int k = 0; k < BLOCK_COLS; k++) {
			mxsum[tidy][tidx] = fmaxf(mxsum[tidy][tidx], mxsum[tidy][k]);
		}
	}

	__syncthreads();

	if (tidy == 0 && tidx == 0) {
		for (int l = 0; l < BLOCK_ROWS; l++) {
			mxsum[tidy][tidx] = fmaxf(mxsum[tidy][tidx], mxsum[l][tidx]);
		}
		mx[blockY*nBlocksX+blockX] = mxsum[tidy][tidx];
	}
}

void ComputeFluxes(double *w, double *hu, double *hv, double *dw, double *dhu,
                   double *dhv, double *mx, double *BC, double *BX, double *BY,
                   double *G, int *active_blocks, double dt, double *n,
                   double hydrograph_rate, int hydrograph_source,
                   double hyetograph_rate, double *hyetograph_gridded_rate,
                   double *F, double *F_old, double *dF, double *K, int *source_idx_dev, double *source_rate_dev, long numSources) {
    ComputeFluxes_k <<< nBlocks, BlockDim >>> (w, hu, hv, dw, dhu, dhv, mx, BC,
	                                           BX, BY, G, active_blocks, dt, pitch,
	                                           pitchBX, pitchBY, n,
	                                           hydrograph_rate,
	                                           hydrograph_source,
	                                           hyetograph_rate,
                                               hyetograph_gridded_rate, 
	                                           F, F_old, dF, K, source_idx_dev, source_rate_dev,numSources);
}

__global__ void Integrate_1_k(double *w, double *hu, double *hv, double *w_old,
                              double *hu_old, double *hv_old, double *dw,
                              double *dhu, double *dhv, double *BC,
                              double *G, bool *wet_blocks, int *active_blocks,
                              double t, double dt, size_t pitch,
                              double hydrograph_rate, int hydrograph_source,
                              double hyetograph_rate,
                              double *hyetograph_gridded_rate, double *F,
                              double *F_old, double *dF, double *K, double *h,
                              double *q, double *h_max, double *q_max,
                              double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources, double *t_peak, double *t_dry) {	//added time_peak and time_dry by Youcan on 20170908
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;

    int blockX = active_blocks[blockIdx.x] % nBlocksX;
    int blockY = active_blocks[blockIdx.x] / nBlocksX;

    int i = blockX*BLOCK_COLS + tidx + 2;
    int j = blockY*BLOCK_ROWS + tidy + 2;

	if (i < 2 || j < 2 || i > nx-3 || j > ny-3) return;

	double *w_ij      = getElement(w,      pitch, j, i);
	double *hu_ij     = getElement(hu,     pitch, j, i);
	double *hv_ij     = getElement(hv,     pitch, j, i);
	double *w_old_ij  = getElement(w_old,  pitch, j, i);
	double *hu_old_ij = getElement(hu_old, pitch, j, i);
	double *hv_old_ij = getElement(hv_old, pitch, j, i);
	double *dw_ij     = getElement(dw,     pitch, j, i);
	double *dhu_ij    = getElement(dhu,    pitch, j, i);
	double *dhv_ij    = getElement(dhv,    pitch, j, i);
	double *BC_ij     = getElement(BC,     pitch, j, i);
	double *G_ij      = getElement(G,      pitch, j, i);

	*w_old_ij  = *w_ij;
	*hu_old_ij = *hu_ij;
	*hv_old_ij = *hv_ij;

    *w_ij  += (*dw_ij)  * dt;
    *hu_ij += (*dhu_ij) * dt;
    *hu_ij /= (1.0f - dt*(*G_ij));
    *hv_ij += (*dhv_ij) * dt;
    *hv_ij /= (1.0f - dt*(*G_ij));

	double hC = ((*w_ij)-(*BC_ij) > epsilon) ? *w_ij - *BC_ij : 0.f;
	*w_ij = (hC > epsilon) ? *w_ij : *BC_ij;

	if (infiltration) {
		double *F_ij     = getElement(F,     pitch, j, i);
		double *F_old_ij = getElement(F_old, pitch, j, i);
		double *dF_ij    = getElement(dF,    pitch, j, i);

		*F_old_ij = *F_ij;
		*F_ij += (*dF_ij) * dt;
	}

	if (euler_integration) {
		__shared__ bool wet[BLOCK_ROWS][BLOCK_COLS];
		wet[tidy][tidx] = false;

		if (save_max) {
			double *h_max_ij = getElement(h_max, pitch, j, i);
			double *q_max_ij = getElement(q_max, pitch, j, i);

			double uC = getVelocity(hC, *hu_ij);
			double vC = getVelocity(hC, *hv_ij);

			//added by Youcan on 20170908
			double *t_peak_ij = getElement(t_peak, pitch, j, i);
			*t_peak_ij = (hC > *h_max_ij) ? t : *t_peak_ij;
			double *t_dry_ij = getElement(t_dry, pitch, j, i);
			*t_dry_ij = (hC > 0.0762) ? t : *t_dry_ij;					//3 in ~ 0.0762 m
			

			*h_max_ij = fmaxf(*h_max_ij, hC);			
			*q_max_ij = fmaxf(*q_max_ij, hC*sqrtf(vC*vC + uC*uC));
			}

		if (h_print || check_volume) {
			double *h_ij = getElement(h, pitch, j, i);
			*h_ij = hC;
		}

		if (q_print) {
			double *q_ij = getElement(q, pitch, j, i);
			double uC = getVelocity(hC, *hu_ij);
			double vC = getVelocity(hC, *hv_ij);
			*q_ij = hC*sqrtf(vC * vC + uC * uC);
			// Add per Dave 11/6/2019 by cr
		}

		if (save_arrival_time) {
			double *t_wet_ij = getElement(t_wet, pitch, j, i);
			if (*t_wet_ij < 0.f && hC > epsilon) {
				*t_wet_ij = t;
			}
		}

		wet[tidy][tidx] = ((*w_ij)-(*BC_ij) > epsilon) ? true : false;

		/** Changed for 1D-2D*/
		if (dambreak||numSources>0) {
			for (int counter = 0; counter < numSources; counter++){
					wet[tidy][tidx] = (j*nx + i == source_idx_dev[counter]) ? true : wet[tidy][tidx];
			}

			//if (j*nx + i == hydrograph_source) {
				//wet[tidy][tidx] = true;
			//}
		}

		if (rainfall_averaged) {
			if (hyetograph_rate > 0.f) {
				wet[tidy][tidx] = true;
			}
		}

		if (rainfall_gridded) {
			double *hyetograph_gridded_rate_ij =
				getElement(hyetograph_gridded_rate, pitch, j, i);
			if ((*hyetograph_gridded_rate_ij) > 0.f) {
				wet[tidy][tidx] = true;
			}
		}

		__syncthreads();

		if (tidx == 0) {
			for (int k = 0; k < BLOCK_COLS; k++) {
				if (wet[tidy][k] == true) {
					wet[tidy][0] = true;
					break;
				}
			}
		}

		__syncthreads();

		if (tidx == 0 && tidy == 0) {
			wet_blocks[blockY*nBlocksX+blockX] = false;
			for (int l = 0; l < BLOCK_ROWS; l++) {
				if (wet[l][0]) {
					wet_blocks[blockY*nBlocksX+blockX] = true;
					break;
				}
			}
		}
	}
}

void Integrate_1(double *w, double *hu, double *hv, double *w_old, double *hu_old,
                 double *hv_old, double *dw, double *dhu, double *dhv, double *BC,
                 double *G, bool *wet_blocks, int *active_blocks, double t, double dt,
                 double hydrograph_rate, int hydrograph_source,
                 double hyetograph_rate, double *hyetograph_gridded_rate,
                 double *F, double *F_old, double *dF, double *K, double *h,
                 double *q, double *h_max, double *q_max, double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources, double *t_peak, double *t_dry) {	//added time_peak and time_dry by Youcan on 20170908
    Integrate_1_k <<< nBlocks, BlockDim >>> (w, hu, hv, w_old, hu_old, hv_old,
	                                         dw, dhu, dhv, BC, G, wet_blocks,
	                                         active_blocks, t, dt, pitch,
	                                         hydrograph_rate, hydrograph_source,
	                                         hyetograph_rate,
                                             hyetograph_gridded_rate, F, F_old,
                                             dF, K, h, q, h_max, q_max, t_wet,source_idx_dev, source_rate_dev, numSources, t_peak, t_dry); 	//added time_peak and time_dry by Youcan on 20170908
}

__global__ void Integrate_2_k(double *w, double *hu, double *hv, double *w_old,
                              double *hu_old, double *hv_old, double *dw,
                              double *dhu, double *dhv, double *BC,
                              double *G, bool *wet_blocks, int *active_blocks,
                              double t, double dt, size_t pitch,
                              double hydrograph_rate, int hydrograph_source,
                              double hyetograph_rate,
                              double *hyetograph_gridded_rate, double *F,
                              double *F_old, double *dF, double *K, double *h,
                              double *q, double *h_max, double *q_max,
                              double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources) {
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;

    __shared__ bool wet[BLOCK_ROWS][BLOCK_COLS];
    wet[tidy][tidx] = false;

    int blockX = active_blocks[blockIdx.x] % nBlocksX;
    int blockY = active_blocks[blockIdx.x] / nBlocksX;

    int i = blockX*BLOCK_COLS + tidx + 2;
    int j = blockY*BLOCK_ROWS + tidy + 2;

	if (i > nx-3 || j > ny-3) return;

	double *w_ij      = getElement(w,      pitch, j, i);
	double *hu_ij     = getElement(hu,     pitch, j, i);
	double *hv_ij     = getElement(hv,     pitch, j, i);
	double *w_old_ij  = getElement(w_old,  pitch, j, i);
	double *hu_old_ij = getElement(hu_old, pitch, j, i);
	double *hv_old_ij = getElement(hv_old, pitch, j, i);
	double *dw_ij     = getElement(dw,     pitch, j, i);
	double *dhu_ij    = getElement(dhu,    pitch, j, i);
	double *dhv_ij    = getElement(dhv,    pitch, j, i);
	double *BC_ij     = getElement(BC,     pitch, j, i);
	double *G_ij      = getElement(G,      pitch, j, i);

	*w_ij  = 0.5f*((*w_old_ij)  + (*w_ij)  + (*dw_ij) *dt);
	*hu_ij = 0.5f*((*hu_old_ij) + (*hu_ij) + (*dhu_ij)*dt);
	*hv_ij = 0.5f*((*hv_old_ij) + (*hv_ij) + (*dhv_ij)*dt);

	double hC = ((*w_ij)-(*BC_ij) > epsilon) ? *w_ij - *BC_ij : 0.f;
	*w_ij = (hC > epsilon) ? *w_ij : *BC_ij;

	if (infiltration) {
		double *F_ij     = getElement(F,     pitch, j, i);
		double *F_old_ij = getElement(F_old, pitch, j, i);
		double *dF_ij    = getElement(dF,    pitch, j, i);

		*F_ij = 0.5f*((*F_old_ij) + (*F_ij) + (*dF_ij)*dt);
	}

	if (save_max) {
		double *h_max_ij = getElement(h_max, pitch, j, i);
		double *q_max_ij = getElement(q_max, pitch, j, i);

		double uC = getVelocity(hC, *hu_ij);
		double vC = getVelocity(hC, *hv_ij);

		*h_max_ij = fmaxf(*h_max_ij, hC);
		*q_max_ij = fmaxf(*q_max_ij, hC*sqrtf(vC*vC + uC*uC));
	}

	if (h_print || check_volume) {
		double *h_ij = getElement(h, pitch, j, i);
		*h_ij = hC;
	}
	
	if (q_print) {
		double *q_ij = getElement(q, pitch, j, i);
		double uC = getVelocity(hC, *hu_ij);
		double vC = getVelocity(hC, *hv_ij);
		*q_ij = hC*sqrtf(vC * vC + uC * uC);
		// Add per Dave 11/6/2019 by cr
	}

	if (save_arrival_time) {
		double *t_wet_ij = getElement(t_wet, pitch, j, i);
		if (*t_wet_ij < 0.f && hC > epsilon) {
			*t_wet_ij = t;
		}
	}

    wet[tidy][tidx] = (hC > epsilon) ? true : false;

	if (dambreak||numSources>0) {
			for (int counter = 0; counter < numSources; counter++){
					wet[tidy][tidx] = (j*nx + i == source_idx_dev[counter]) ? true : wet[tidy][tidx];
			}
		
		//if (j*nx + i == hydrograph_source) {
			//wet[tidy][tidx] = true;
		//}
	}

	if (rainfall_averaged) {
		if (hyetograph_rate > 0.f) {
			wet[tidy][tidx] = true;
		}
	}

	if (rainfall_gridded) {
		double *hyetograph_gridded_rate_ij =
			getElement(hyetograph_gridded_rate, pitch, j, i);
		if ((*hyetograph_gridded_rate_ij) > 0.f) {
			wet[tidy][tidx] = true;
		}
	}

    __syncthreads();

	// Do one reduction in parallel
	if (tidx == 0) {
		for (int k = 0; k < BLOCK_COLS; k++) {
			if (wet[tidy][k] == true) {
				wet[tidy][0] = true;
				break;
			}
		}
	}

	__syncthreads();

	// Perform the final reduction
	if (tidx == 0 && tidy == 0) {
		wet_blocks[blockY*nBlocksX+blockX] = false;
		for (int l = 0; l < BLOCK_ROWS; l++) {
			if (wet[l][0]) {
				wet_blocks[blockY*nBlocksX+blockX] = true;
				return;
			}
		}
	}
}

void Integrate_2(double *w, double *hu, double *hv, double *w_old, double *hu_old,
                 double *hv_old, double *dw, double *dhu, double *dhv, double *BC, 
                 double *G, bool *wet_blocks, int *active_blocks, double t, double dt,
                 double hydrograph_rate, int hydrograph_source,
                 double hyetograph_rate, double *hyetograph_gridded_rate,
                 double *F, double *F_old, double *dF, double *K, double *h,
                 double *q, double *h_max, double *q_max, double *t_wet,int *source_idx_dev, double *source_rate_dev, long numSources) {
    Integrate_2_k <<< nBlocks, BlockDim >>> (w, hu, hv, w_old, hu_old, hv_old,
	                                         dw, dhu, dhv, BC, G, wet_blocks,
	                                         active_blocks, t, dt, pitch,
	                                         hydrograph_rate, hydrograph_source,
	                                         hyetograph_rate,
                                             hyetograph_gridded_rate, F, F_old,
                                             dF, K, h, q, h_max, q_max, t_wet,source_idx_dev, source_rate_dev, numSources);
}

__global__ void Grow_k(bool *wet_blocks, int *active_blocks) {
    int i = blockIdx.x*BLOCK_COLS + threadIdx.x;
    int j = blockIdx.y*BLOCK_ROWS + threadIdx.y;
    int col = threadIdx.x + 1, row = threadIdx.y + 1;

 	__shared__ bool swet[BLOCK_ROWS+2][BLOCK_COLS+2];
    swet[threadIdx.y][threadIdx.x] = false;

	if (threadIdx.x == 0 || i == 0) {
		swet[row  ][col-1] = (i > 0 && j < nBlocksY) ? wet_blocks[j*nBlocksX+(i-1)] : false;
	} else if (threadIdx.x == BLOCK_COLS-1 || i == nBlocksX-1) {
		swet[row  ][col+1] = (i < nBlocksX-1 && j < nBlocksY) ? wet_blocks[j*nBlocksX+(i+1)] : false;
	}

	if (threadIdx.y == 0  || j == 0) {
		swet[row-1][col  ] = (j > 0) ? wet_blocks[(j-1)*nBlocksX+i] : false;
	} else if (threadIdx.y == BLOCK_ROWS-1 || j == nBlocksY-1) {
		swet[row+1][col  ] = (j < nBlocksY-1) ? wet_blocks[(j+1)*nBlocksX+i] : false;
	}

    if (i > nBlocksX-1 || j > nBlocksY-1) {
        return;
    } else {
        swet[row][col] = wet_blocks[j*nBlocksX+i];
    }

	__syncthreads();

    if (swet[row  ][col  ] || swet[row  ][col+1] || swet[row  ][col-1] ||
        swet[row+1][col  ] || swet[row-1][col  ]) {
        active_blocks[j*nBlocksX+i] = j*nBlocksX+i;
    } else {
        active_blocks[j*nBlocksX+i] = -(j*nBlocksX + i);
    }
}

__global__ void TrackRainfall_k(bool *wet_blocks, int *active_blocks, double *precip_dist, size_t pitch) {
	
	// already wet set up by other threads
	if (wet_blocks[abs(active_blocks[blockIdx.x])]) return;

	// Goal: find the global block ID# and set the wet_blocks
	// sorted active_blocks give the src locations
	int tidx = threadIdx.x;
	int tidy = threadIdx.y;

	// MJ: Why do these have abs() called? Blocks cant be negative?
	int blockX = abs(active_blocks[blockIdx.x]) % nBlocksX;
	int blockY = abs(active_blocks[blockIdx.x]) / nBlocksX;

	int i = blockX * BLOCK_COLS + tidx + 2;
	int j = blockY * BLOCK_COLS + tidy + 2;

	if (i < 2 || j < 2 || i > nx - 3 || j > ny -3) return;

	// check rain sources
	double *precip_ij = getElement(precip_dist, pitch, j, i);
	if (*precip_ij > epsilon) {
		wet_blocks[abs(active_blocks[blockIdx.x])] = true;
	}
}

void Grow(bool *wet_blocks, int *active_blocks, double *hyetograph_gridded_rate, bool h_rainfall_gridded) {
	if ( h_rainfall_gridded ) {
		thrust::device_ptr<int> dp_active_blocks = thrust::device_pointer_cast(active_blocks);
		thrust::partition(dp_active_blocks, dp_active_blocks + GridSize, is_negative());
		nBlocks = thrust::count_if(dp_active_blocks, dp_active_blocks + GridSize, is_negative());
	
		// Mark Jensen 12/19/19 - Added if statement this was erroring out in CUDA versions g.t. 9.0
		// Thrust now checks thrust::partition below if cudaErrors are present if nBlock = 0 then a CUDA 
		// error is thrown, this just tells the kernel to run when there are valid blocks
		if (nBlocks > 0) {
			TrackRainfall_k << < nBlocks, BlockDim >> > (wet_blocks, active_blocks, hyetograph_gridded_rate, pitch);
		}
	}

    Grow_k <<< GrowGrid, BlockDim >>> (wet_blocks, active_blocks);

    thrust::device_ptr<int> dp_active_blocks = thrust::device_pointer_cast(active_blocks);
    thrust::partition(dp_active_blocks, dp_active_blocks + GridSize, is_nonnegative());

//    thrust::remove(thrust::device, dp_active_blocks, dp_active_blocks + GridSize, -1);
    nBlocks = thrust::count_if(dp_active_blocks, dp_active_blocks + GridSize, is_nonnegative());
//    std::cout << nBlocks << std::endl;
/*    nBlocks = thrust::count_if(dp_active_blocks, dp_active_blocks +
                               GridSize, is_nonnegative());*/
    /*thrust::stable_sort(dp_active_blocks, dp_active_blocks + GridSize,
                        thrust::greater<int>());*/

//    thrust::copy_if(dp_active_blocks, dp_active_blocks + GridSize,
//                    dp_exe_blocks, is_nonnegative());

    thrust::device_ptr<bool> dp_wet_blocks   = thrust::device_pointer_cast(wet_blocks);
    int wetblocks = GridSize - thrust::count(dp_wet_blocks, dp_wet_blocks +
                                       GridSize, (bool)false);
   // std::cout << wetblocks << "\t" << nBlocks << std::endl;

}
    

__global__ void ApplyBoundariesX_k(double *w, double *hu, double *hv,
                                   double *BC, size_t pitch) {
	int j = blockIdx.x*blockDim.x + threadIdx.x;

	if (j > ny-3 || j < 2) return;

	double *w_L    = getElement(w,  pitch, j, 0);
	double *hu_L   = getElement(hu, pitch, j, 0);
	double *hv_L   = getElement(hv, pitch, j, 0);
	double *w_Lp1  = getElement(w,  pitch, j, 1);
	double *hu_Lp1 = getElement(hu, pitch, j, 1);
	double *hv_Lp1 = getElement(hv, pitch, j, 1);
	double *w_Lp2  = getElement(w,  pitch, j, 2);
	double *hu_Lp2 = getElement(hu, pitch, j, 2);
	double *hv_Lp2 = getElement(hv, pitch, j, 2);
	double *BC_Lp2 = getElement(BC, pitch, j, 2);

	double HL = fminf(powf((*hu_Lp2)*(*hu_Lp2) / g, 1.f/3.f),
					 (*w_Lp2)-(*BC_Lp2));
	HL = (HL > epsilon) ? HL : 0.f;
	*w_L  = *w_Lp1  = (*BC_Lp2) + HL;
	*hu_L = *hu_Lp1 = (*hu_Lp2);
	*hv_L = *hv_Lp1 = (*hv_Lp2);

	double *w_R    = getElement(w,  pitch, j, nx-1);
	double *hu_R   = getElement(hu, pitch, j, nx-1);
	double *hv_R   = getElement(hv, pitch, j, nx-1);
	double *w_Rm1  = getElement(w,  pitch, j, nx-2);
	double *hu_Rm1 = getElement(hu, pitch, j, nx-2);
	double *hv_Rm1 = getElement(hv, pitch, j, nx-2);
	double *w_Rm2  = getElement(w,  pitch, j, nx-3);
	double *hu_Rm2 = getElement(hu, pitch, j, nx-3);
	double *hv_Rm2 = getElement(hv, pitch, j, nx-3);
	double *BC_Rm2 = getElement(BC, pitch, j, nx-3);

    double HR = fminf(powf((*hu_Rm2)*(*hu_Rm2) / g, 1.f/3.f),
                     (*w_Rm2)-(*BC_Rm2));
	HR = (HR > epsilon) ? HR : 0.f;
    *w_R  = *w_Rm1  = (*BC_Rm2) + HR;
    *hu_R = *hu_Rm1 = (*hu_Rm2);
    *hv_R = *hv_Rm1 = (*hv_Rm2);
}

__global__ void ApplyBoundariesY_k(double *w, double *hu, double *hv,
                                   double *BC, size_t pitch) {
    int i = blockIdx.x*blockDim.x + threadIdx.x; 
	if (i < 2 || i > nx-3) return;

	double *w_B    = getElement(w,  pitch, 0, i);
	double *hu_B   = getElement(hu, pitch, 0, i);
	double *hv_B   = getElement(hv, pitch, 0, i);
	double *w_Bp1  = getElement(w,  pitch, 1, i);
	double *hu_Bp1 = getElement(hu, pitch, 1, i);
	double *hv_Bp1 = getElement(hv, pitch, 1, i);
	double *w_Bp2  = getElement(w,  pitch, 2, i);
	double *hu_Bp2 = getElement(hu, pitch, 2, i);
	double *hv_Bp2 = getElement(hv, pitch, 2, i);
	double *BC_Bp2 = getElement(BC, pitch, 2, i);

    double HB = fminf(powf((*hv_Bp2)*(*hv_Bp2) / g, 1.f/3.f),
                     (*w_Bp2)-(*BC_Bp2));
    *w_B  = *w_Bp1  = (*BC_Bp2) + HB;
    *hu_B = *hu_Bp1 = (*hu_Bp2);
    *hv_B = *hv_Bp1 = (*hv_Bp2);

	double *w_T    = getElement(w,  pitch, ny-1, i);
	double *hu_T   = getElement(hu, pitch, ny-1, i);
	double *hv_T   = getElement(hv, pitch, ny-1, i);
	double *w_Tm1  = getElement(w,  pitch, ny-2, i);
	double *hu_Tm1 = getElement(hu, pitch, ny-2, i);
	double *hv_Tm1 = getElement(hv, pitch, ny-2, i);
	double *w_Tm2  = getElement(w,  pitch, ny-3, i);
	double *hu_Tm2 = getElement(hu, pitch, ny-3, i);
	double *hv_Tm2 = getElement(hv, pitch, ny-3, i);
	double *BC_Tm2 = getElement(BC, pitch, ny-3, i);

    double HT = fminf(powf((*hv_Tm2)*(*hv_Tm2) / g, 1.f/3.f),
                     (*w_Tm2)-(*BC_Tm2));
    *w_T  = *w_Tm1  = (*BC_Tm2) + HT;
    *hu_T = *hu_Tm1 = (*hu_Tm2);
    *hv_T = *hv_Tm1 = (*hv_Tm2);
}

void ApplyBoundaries(double *w, double *hu, double *hv, double *BC) {
    ApplyBoundariesX_k <<< GridDimX, BlockDimX >>> (w, hu, hv, BC, pitch);
    ApplyBoundariesY_k <<< GridDimY, BlockDimY >>> (w, hu, hv, BC, pitch);
}

void FreeGrid(double *&w, double *&hu, double *&hv, double *&w_old, double *&hu_old,
              double *&hv_old, double *&dw, double *&dhu, double *&dhv, double *&mx,
              double *&BC, double *&BX, double *&BY, bool *&wet_blocks,
              int *&active_blocks, double *&n, double *&hyetograph_gridded_rate,
              double *&F, double *&F_old, double *&dF, double *&K, double *&h,
              double *&q, double *&h_max, double *&q_max, double *&t_wet, double *&t_peak, 
		      double *&t_dry, double *&G) {	//added t_peak and t_dry by Youcan on 20170908
    cudaFree(w);
    cudaFree(hu);
    cudaFree(hv);
    cudaFree(w_old);
    cudaFree(hu_old);
    cudaFree(hv_old);
    cudaFree(dw);
    cudaFree(dhu);
    cudaFree(dhv);
    cudaFree(G);
    cudaFree(BC);
    cudaFree(BX);
    cudaFree(BY);
    cudaFree(mx);
    cudaFree(wet_blocks);
    cudaFree(active_blocks);
	cudaFree(hyetograph_gridded_rate);
	cudaFree(n);
	cudaFree(F);
	cudaFree(F_old);
	cudaFree(dF);
	cudaFree(K);
	cudaFree(h);
	cudaFree(q);
	cudaFree(h_max);
	cudaFree(q_max);
	cudaFree(t_wet);

	cudaFree(t_peak);				//added by Youcan on 20170908
	cudaFree(t_dry);				//added by Youcan on 20170908
}

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
