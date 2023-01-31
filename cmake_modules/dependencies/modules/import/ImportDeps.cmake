message("")
message("---------------------------------------- Dependencies")
message(STATUS "Checking for Exageostat Dependencies")

set(ENV{PKG_CONFIG_PATH} ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})

if(MPI_VALUE)
    message("Trying to find MPI")
    find_package(MPI REQUIRED)
    message("FOUND MPI.....")
endif()

include(ImportGSL)
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL Installation failed")
endif()

include(ImportHWLOC)
if(NOT HWLOC_FOUND)
    message(FATAL_ERROR "HWLOC Installation failed")
endif()

include(ImportStarPu)
if(NOT STARPU_FOUND)
    message(FATAL_ERROR "STARPU Installation failed")
endif()