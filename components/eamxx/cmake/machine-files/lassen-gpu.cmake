# common setup
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()


message(STATUS "pm-gpu PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")

# set up Kokkos
include(${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# gpu specific
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

# common to lassen-cpu
# set up MPI
include(${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
set(EKAT_MPIRUN_EXE "jsrun -E LD_PRELOAD=/opt/ibm/spectrum_mpi/lib/pami_471/libpami.so" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "--np" CACHE STRING "" FORCE)

# set up TPLs
set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "" FORCE)
set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "" FORCE)

# set up input root
set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata/" CACHE STRING "")
