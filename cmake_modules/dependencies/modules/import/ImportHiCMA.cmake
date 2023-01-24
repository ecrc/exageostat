message("")
message(STATUS "Installing Hicma")

execute_process(COMMAND ./InstallHiCMA.sh --prefix ${CMAKE_INSTALL_PREFIX} --setup ${TMP_DIR} --mpi ${MPI_VALUE} --build ${CMAKE_BINARY_DIR}
                WORKING_DIRECTORY ${ECRC_INSTALLATION_SCRIPTS_PATH}
                RESULT_VARIABLE res)
if(${res} EQUAL 0)
    set(HiCMA_INSTALLED TRUE)
endif()
