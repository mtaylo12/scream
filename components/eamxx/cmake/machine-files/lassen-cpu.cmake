include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# set up Kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)


# set up MPI (done in included file?)
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)

set(EKAT_MPIRUN_EXE "jsrun -E LD_PRELOAD=/opt/ibm/spectrum_mpi/lib/pami_471/libpami.so" CACHE STRING "" FORCE)

set(EKAT_MPI_NP_FLAG "--np" CACHE STRING "" FORCE)

# set up TPLs
#set(NetCDF_Fortran_PATH /usr/gdata/climdat/libs/netcdf-fortran/install/lassen/fortran CACHE STRING "")
#set(BLAS_LIBRARIES /usr/gdata/climdat/libs/blas/libblas.a CACHE STRING "")
#set(LAPACK_LIBRARIES /usr/gdata/climdat/libs/lapack/liblapack.a CACHE STRING "")

# set up input root
set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata/" CACHE STRING "")
