message("")
message(STATUS "Installing Chameleon")

execute_process(COMMAND ./InstallChameleon.sh --prefix ${CMAKE_INSTALL_PREFIX} --setup ${TMP_DIR} --cuda ${CUDA_VALUE} --mpi ${MPI_VALUE}
        WORKING_DIRECTORY ${ECRC_INSTALLATION_SCRIPTS_PATH}
        RESULT_VARIABLE res)
if(${res} EQUAL 0)
    set(CHAMELEON_INSTALLED TRUE)
endif()
