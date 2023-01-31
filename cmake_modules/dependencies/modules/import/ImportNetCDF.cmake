message("")
message(STATUS "Installing NetCDF")

execute_process(COMMAND ./InstallNetCDF.sh --prefix ${CMAKE_INSTALL_PREFIX} --setup ${TMP_DIR}
        WORKING_DIRECTORY ${ECRC_INSTALLATION_SCRIPTS_PATH}
        RESULT_VARIABLE res)
if(${res} EQUAL 0)
    set(NETCDF_INSTALLED TRUE)
endif()