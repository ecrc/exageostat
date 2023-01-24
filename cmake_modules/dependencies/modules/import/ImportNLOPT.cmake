message("")
message(STATUS "Installing NLOPT")

execute_process(COMMAND ./InstallNLOPT.sh --prefix ${CMAKE_INSTALL_PREFIX} --setup ${TMP_DIR}
                WORKING_DIRECTORY ${ECRC_INSTALLATION_SCRIPTS_PATH}
                RESULT_VARIABLE res)
if(${res} EQUAL 0)
    set(NLOPT_INSTALLED TRUE)
endif()