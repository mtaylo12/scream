INCLUDE (ScreamUtils)

# Create the test
set (TEST_LABELS "shoc;cld;p3;rrtmgp;physics;driver")
set (NEED_LIBS shoc cld_fraction p3 scream_rrtmgp scream_control scream_share physics_share diagnostics)
CreateUnitTest(shoc_cld_p3_rrtmgp shoc_cld_p3_rrtmgp.cpp "${NEED_LIBS}" LABELS ${TEST_LABELS}
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  PROPERTIES FIXTURES_SETUP shoc_cld_p3_rrtmgp_generate_output_nc_files
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 4h 24h
set (RUN_T0 2021-10-12-45000)

# Determine num subcycles needed to keep shoc dt<=300s
set (SHOC_MAX_DT 300)
math (EXPR MAC_MIC_SUBCYCLES "(${ATM_TIME_STEP} + ${SHOC_MAX_DT} - 1) / ${SHOC_MAX_DT}")

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

## Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_p3_rrtmgp_output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_p3_rrtmgp_output.yaml)

## Finally compare all MPI rank output files against the single rank output as a baseline, using CPRNC
## Only if running with 2+ ranks configurations
# This test requires CPRNC
if (TEST_RANK_END GREATER TEST_RANK_START)
  include (BuildCprnc)
  BuildCprnc()

  set (BASE_TEST_NAME "shoc_cld_p3_rrtmgp")
  math (EXPR CMP_RANK_START ${TEST_RANK_START}+1)
  foreach (MPI_RANKS RANGE ${CMP_RANK_START} ${TEST_RANK_END})

    set (SRC_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${MPI_RANKS}.${RUN_T0}.nc")
    set (TGT_FILE "${BASE_TEST_NAME}_output.INSTANT.nsteps_x${NUM_STEPS}.np${TEST_RANK_START}.${RUN_T0}.nc")
    set (TEST_NAME "${BASE_TEST_NAME}_np${TEST_RANK_START}_vs_np${MPI_RANKS}_bfb")
    add_test (NAME ${TEST_NAME}
              COMMAND cmake -P ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}"
              RESOURCE_LOCK ${BASE_TEST_NAME}
              FIXTURES_REQUIRED shoc_cld_p3_rrtmgp_generate_output_nc_files)
  endforeach()
endif()
