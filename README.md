# Introduction 

# Build and Test

## PIC

To build on deception or marianas, the following modules are required:
```
module load git
module load gcc/8.1.0
module load cuda/11.4
module load openmpi
module load cmake
```
Check out a copy of the code and make a place to build
```
git checkout [...]
cd RIFT
mkdir build
cd build
cmake \
    -D CUDA_SAMPLES_INC=/share/apps/cuda/11.4/samples/common/inc/ \
    -D CMAKE_VERBOSE_MAKEFILE=TRUE \
    -D CMAKE_CUDA_ARCHITECTURES="60;70;75;80" \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D RIFT_ENABLE_TESTS:BOOL=YES \
    ../src
make -j
```
This should create an executable binary `rift`.  The specifed `CUDA ARCHITECTURES` are those available through the various GPU equipped partitions on PIC.

# RIFT Configuration

Configuration is read from a file containing a list of phrases like
```
keyword = value
```
which define RIFT operation. The available phrases are listed below.  Those marked "req'd" must appear in the file. The default value is shown for those not required.  

## Simulation Time Frame

 - `t0 = float`: simulation start time, s (req'd)
 - `tf = float`: simulation end time, s (req'd)

## Simulation Modes

 - `euler_integration = 0/1`: Use backward Euler integration (0)
 - `courant_limit = float`: Maximum Courant number (0.25)
 - `check_volume = 0/1`: print some messages about mass balance (0)
 - `device = int`: use the specified GPU (0)
 - `square_cells = 0/1`: (1)
 - `projected = 0/1`: (0)
 - `tile_acceleration = 0/1`: (1)

## Static and Initial Conditions

 - `DEM = file`: read the DEM/bathymetry raster map from `file` (req'd)
 - `h0 = file`: read the initial depth map from `file` (0.0 everywhere)
 - `n = float`: use a single Manning's coefficient everywhere (0.0), or
 - `n = file`: read a map of Manning's coefficient from `file`
 - `read_hotstart = file`: read a hotstart file (`hot#.bin`) from another simulation

## Rainfall 

Rainfall is optional. It can be specified as time varying rate over the entire domain using

 - `hyetograph = file`: read rainfall rate (mm/hr) from CSV `file`

Or, as a time varying map of rainfall rate using

 - `hyetograph_prefix = string`: prefix string for rainfall map file
   names, which will be `string-###.txt`, where `###` is the step index
 - `hyetograph_interp = 0/1`: interpolate between rainfall rate maps (0)
 - `hyetograph_dt = float`: time step, s, between rainfall rate maps (req'd)
 - `hyetograph_tf = float`: time, s, of last rainfall map (req'd)


## Infiltration

Simulation of infiltration is optional. To enable, all of the following phrases must appear:

 - `K = file`: read map of infiltration capacity from `file`
 - `psi = `:
 - `dtheta = `:

## Hydrograph

A single discharge hydrograph simulate a source

 - `hydrograph = file`: read a discharge (cfs) from CSV `file`

If `hydrograph` appears, the source location as a point woith

 - `source_X = float`: point longitude/easting
 - `source_Y = float`: point latitude/northing
 
or with

 - `source_pt = (float, float)`:  longitude/easting and latitude/northing of point
 
or as a rectangle

 - `source_rect_X = (float, float)`: minimum and maximum longitude/easting of area
 - `source_rect_Y = (float, float)`: minimum and maximum latitude/northing of area

## Drainage

A drainage rate can be applied to the simulation domain, either as a
time-varying rate over the entire domain, using

 - `drain = file`: read drainage rate (mm/hr) time series from `file` (none)

or as a time-series of maps, using

 - `drain_prefix = string`: prefix string for map file names, like (`hyetograph_prefix`)
 - `drain_dt = float`: time step, s, between drain rate maps (req'd)
 - `drain_tf = float`: time, s, of last drain map (req'd)

## Surge

Surge is a way to represent areas of  known depths in the domain.
This is typically used to represent storm surge in coastal areas.  A
time series of depth (m) maps is specified with

 - `surge_prefix = string`: prefix string for map file names, like (`hyetograph_prefix`)
 - `surge_dt = float`: time step, s, between surge maps (req'd)
 - `surge_tf = float`: time, s, of last surge map (req'd)

## Output

Various raster map outputs can be produced during the simulation.  Output of simulation results is controlled by

 - `output = path`: directory to write output maps (`./`)
 - `h_print = 0/1`: output simulated depth (0)
 - `q_print = 0/1`: output simulated unit discharge (0)
 - `ws_print = 0/1`: output simulated water surface elevation (0)
 - `dt_print = float`: time between output maps, s (req'd if any print enabled)
 - `save_max = 0/1`: save peak depth, and unit discharge (0)
 - `save_arrival = 0/1`: save flood wave arrival time (0)
 - `save_hotstart = 0/1`: save hotstart 

Some static maps can also be output for debugging purposes

 - `b_print = 0/1`: output internal DEM/bathymetry (0)
 - `n_print = 0/1`: output internal Manning's n map (0)


