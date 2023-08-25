cmake_minimum_required(VERSION 3.20)

project(rift LANGUAGES C CXX CUDA)

if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()

# find_package(MPI QUIET)
find_package(MPI REQUIRED)
if(MPI_FOUND)
  set(MPI_CDEFS USE_MPI)
endif()

# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

find_package(CUDAToolkit REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 60)
endif()

# message(STATUS " ** CUDA_INCLUDE_DIR:" ${CUDA_INCLUDE_DIRS})
message(STATUS " ** CUDA_TOOLKIT_INCLUDE_DIR:" ${CUDAToolkit_INCLUDE_DIRS})
# message(STATUS " ** CMAKE_CXX_FLAGS:" ${CMAKE_CXX_FLAGS})
# message(STATUS " ** CXXFLAGS:" ${CXXFLAGS})

if(NOT DEFINED CUDA_SAMPLES_INC)
  message(FATAL_ERROR "cuda samples inc dir not specified")
endif()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")

# set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -m64; -O3;
# 										-use_fast_math;
# 										-ftz=true; -prec-div=false; -prec-sqrt=false;
# 										-Xptxas -dlcm=ca; -lineinfo;
# 										-Wno-deprecated-gpu-targets"
# 										)

add_executable(rift 
				# io.h 
				io.cc 
				# source.h
				source.cc 
				# constants.h 
				# config.h 
				# grid.h 
				# grid.cc
				grid.cu
                                GridSeries.cu
				# simulator.h 
				simulator.cu 
				main.cc
			)


# set_source_files_properties(test.cpp PROPERTIES LANGUAGE CUDA)								
# set_property(TARGET rift PROPERTY CUDA_ARCHITECTURES 60 70)
set_target_properties(rift PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(rift PUBLIC ${CUDAToolkit_INCLUDE_DIRS} ${CUDA_SAMPLES_INC} )

if(MPI_FOUND)
	target_include_directories(rift PUBLIC ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
	target_compile_definitions(rift PUBLIC ${MPI_CDEFS})
	target_link_libraries(rift PUBLIC MPI::MPI_C)
endif()
