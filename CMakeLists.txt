# @copyright (c) 2017-2023 King Abdullah University of Science and
#                      Technology (KAUST). All rights reserved.
#
# ExaGeoStat is a software package, provided by King Abdullah
#             University of Science and Technology (KAUST)
#
# @file CMakeLists.txt
# @version 1.2.0
# @author Eduardo Gonzalez Fisher
# @author Sameh Abdulah
# @date 2022-11-09


###############################################################################
#                       THIS IS A TOP-LEVEL CMAKELISTS.txt                    #
#                                                                             #
#        It is intended to find all dependencies (required or optional)       #
#                    and set up corresponding variables                       #
###############################################################################


###############################################################################
##                            PRELIMINARIES                                  ##
###############################################################################

#Minimum CMAKE version.
cmake_minimum_required(VERSION 2.8.12)

# directly make an error if in-source build
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.\n"
            "Please create a build directory first and execute cmake configuration from "
            "this directory. Example: mkdir build && cd build && cmake ..")
endif ()

project(EXAGEOSTAT C)
# set project version number
set(EXAGEOSTAT_VERSION_MAJOR 1)
set(EXAGEOSTAT_VERSION_MINOR 2)
set(EXAGEOSTAT_VERSION_PATCH 0)

#set(CHAMELEON_CMAKE_DIR "" CACHE PATH "Directory of CHAMELEON CMake modules, can be external to the project")
# CMAKE MODULES :: ECRC
# REQUIRED FOR TESTS TO LINK LIBRARIES
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules")
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if (${_res_init} GREATER 0 OR ${_res_update} GREATER 0)
        message(FATAL_ERROR "ECRC CMake modules were not found.\n"
                "We tried: 'git submodule init && git submodule update' and resulted in error")
    endif ()
endif ()
# ECRC INITIALIZATION
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/dependencies/modules/import")
set(ECRC_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules)
set(ECRC_INSTALLATION_SCRIPTS_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/dependencies/scripts)

include(EcrcInit)
include(GenPkgConfig)

set(EXAGEOSTAT_VERSION ${EXAGEOSTAT_VERSION_MAJOR}.${EXAGEOSTAT_VERSION_MINOR}.${EXAGEOSTAT_VERSION_PATCH})
add_definitions("-DVERSION=${EXAGEOSTAT_VERSION}")

include(FindPkgConfig)
find_package(PkgConfig QUIET)

#############################################
#                                           #
#        Compilation of EXAGEOSTAT           #
#                                           #
#############################################

###############################################################################
# Parameters/Options #
######################

# Set the RPATH config
# --------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(EXAGEOSTAT_EXAMPLES "Build ExaGeoStat Examples" ON)

# Define precision supported by EXAGEOSTAT
# -----------------------------------------
set(RP_EXAGEOSTAT_DICTIONNARY ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules/precision_generator/subs.py)
set(RP_EXAGEOSTAT_PRECISIONS "s;d;c;z")
include(RulesPrecisions)

# Install Dependencies for EXAGEOSTAT
# ------------------------------------
if (EXAGEOSTAT_INSTALL_DEPS)
    if (NOT DEFINED TMP_DIR)
        string(RANDOM LENGTH 24 _s)
        set(TMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/${_s})
        execute_process(COMMAND mkdir -p ${TMP_DIR} RESULT_VARIABLE res)
    endif ()
    include(ImportDeps)
endif ()
# Create a list of possible runtime
set(EXAGEOSTAT_SCHED_list "STARPU" "QUARK")

# Initially StarPU runtime is enabled
option(EXAGEOSTAT_SCHED_STARPU
        "Enable StarPU scheduler as the default runtime
    (Conflict with other EXAGEOSTAT_SCHED_* options)" ON)
option(EXAGEOSTAT_SCHED_QUARK
        "Enable Quark scheduler as the default runtime
        (Conflict with other EXAGEOSTAT_SCHED_* options)" OFF)

# For now, we are able to compile only one runtime at a time, so we disable combinations of runtimes
if (EXAGEOSTAT_SCHED_STARPU)
    set(EXAGEOSTAT_SCHED_QUARK OFF)
elseif (EXAGEOSTAT_SCHED_QUARK)
    set(EXAGEOSTAT_SCHED_STARPU OFF)
endif ()
# Set default to StarPU if nothing specific is required by the user
if (NOT EXAGEOSTAT_SCHED_STARPU AND NOT EXAGEOSTAT_SCHED_QUARK)
    set(EXAGEOSTAT_SCHED_STARPU ON)
endif ()
if (EXAGEOSTAT_SCHED_STARPU)
    message("-- ${BoldGreen}EXAGEOSTAT_SCHED_STARPU is set to ON: EXAGEOSTAT uses StarPU runtime\n"
            "   To use EXAGEOSTAT with Quark  runtime: set EXAGEOSTAT_SCHED_QUARK  to ON\n"
            "   (EXAGEOSTAT_SCHED_STARPU will be disabled)${ColourReset}")
elseif (EXAGEOSTAT_SCHED_QUARK)
    message("-- ${BoldGreen}EXAGEOSTAT_SCHED_QUARK is set to ON: EXAGEOSTAT uses Quark runtime\n"
            "   To use EXAGEOSTAT with StarPU runtime: set EXAGEOSTAT_SCHED_STARPU to ON\n"
            "   (EXAGEOSTAT_SCHED_QUARK will be disabled)${ColourReset}")
endif ()


# Check that one, and only one, SCHED option is set to ON
# count number of runtime sets to ON
math(EXPR number_of_active_runtime 0)
foreach (runtime ${EXAGEOSTAT_SCHED_list})
    if (EXAGEOSTAT_SCHED_${runtime})
        math(EXPR number_of_active_runtime "${number_of_active_runtime}+1")
    endif ()
endforeach ()
if (NOT number_of_active_runtime STREQUAL 1)
    message(FATAL_ERROR
            "Number of active runtime is ${number_of_active_runtime}, "
            "the user should activate one (and only one) runtime. ")
endif ()


# Warning if not coherent options
if (EXAGEOSTAT_USE_MPI AND NOT EXAGEOSTAT_SCHED_STARPU)
    message(WARNING "You have activated MPI but EXAGEOSTAT_SCHED_STARPU is OFF.\n"
            "You should use StarPU Runtime system if you want to benefit from MPI.\n"
            "Use -DEXAGEOSTAT_SCHED_QUARK=OFF -DEXAGEOSTAT_SCHED_STARPU=ON at cmake \n"
            "configure to do so.")
endif ()

# Additional options
# ------------------

# Enable the distributed interface (allowed only when StarPU is enabled)
cmake_dependent_option(EXAGEOSTAT_USE_MPI
        "Enable distributed memory through MPI" OFF
        "EXAGEOSTAT_SCHED_STARPU" OFF)
if (NOT EXAGEOSTAT_USE_MPI)
    message("-- ${BoldGreen}EXAGEOSTAT_USE_MPI is set to OFF, turn it ON to use MPI (only with StarPU)${ColourReset}")
endif ()

# Enable NETCDF support
if (NOT EXAGEOSTAT_USE_NETCDF)
    message("-- ${BoldGreen}EXAGEOSTAT_USE_NETCDF is set to OFF, turn it ON to enable NetCDF support${ColourReset}")
endif ()

# Use intermediate variable since cmake_dependent_option doesn't have OR conditions
set(EXAGEOSTAT_ENABLE_MPI OFF CACHE INTERNAL "Tells if MPI might be supported by the runtime")
if (EXAGEOSTAT_SCHED_STARPU)
    set(EXAGEOSTAT_ENABLE_MPI ON FORCE)
endif ()

###############################################################################
# Look for dependencies #
#########################
set(EXAGEOSTAT_DEP "")

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if (THREADS_FOUND)
    list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
set(M_LIBRARIES "")
if (UNIX OR WIN32)
    find_library(
            M_m_LIBRARY
            NAMES m
    )
    mark_as_advanced(M_m_LIBRARY)
    if (M_m_LIBRARY)
        list(APPEND M_LIBRARIES "${M_m_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${M_m_LIBRARY}")
    else ()
        message(FATAL_ERROR "Could NOT find libm on your system."
                " Are you sure to a have a C compiler installed?")
    endif ()
endif ()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
set(RT_LIBRARIES "")
if (UNIX AND NOT APPLE)
    find_library(
            RT_rt_LIBRARY
            NAMES rt
    )
    mark_as_advanced(RT_rt_LIBRARY)
    if (RT_rt_LIBRARY)
        list(APPEND RT_LIBRARIES "${RT_rt_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${RT_rt_LIBRARY}")
    else ()
        message(FATAL_ERROR "Could NOT find librt on your system")
    endif ()
endif ()

# EXAGEOSTAT depends on CBLAS
#---------------------------
find_package(CBLAS COMPONENTS BLASEXT)
if (BLAS_FOUND)
    if (BLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
    endif ()
    if (BLAS_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
    endif ()
else ()
    message(FATAL_ERROR "BLAS library has not been found")
endif ()
if (CBLAS_FOUND)
    include_directories(${CBLAS_INCLUDE_DIRS})
    if (CBLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS}")
    endif ()
    if (CBLAS_LIBRARIES)
        if (CBLAS_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES_DEP})
        else ()
            list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES})
        endif ()
    endif ()
else ()
    if (ECRC_VERBOSE_FIND_PACKAGE)
        if (CBLAS_STANDALONE OR NOT CBLAS_WORKS)
            if (NOT CBLAS_cblas.h_DIRS)
                Print_Find_Header_Status(cblas cblas.h)
            endif ()
            if (NOT CBLAS_cblas_LIBRARY)
                Print_Find_Library_Status(cblas libcblas)
            endif ()
        endif ()
    else ()
        message(WARNING "CBLAS library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
                " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif ()
    message(FATAL_ERROR "A CBLAS library is required but has not been found")
endif ()

# EXAGEOSTAT depends on LAPACKE
#-----------------------------
# standalone version of lapacke seems useless for now
# let the comment in case we meet some problems of non existing lapacke
# functions in lapack library such as mkl, acml, ...
#set(LAPACKE_STANDALONE TRUE)
find_package(LAPACKE COMPONENTS LAPACKEXT)
if (LAPACK_FOUND AND LAPACK_LIBRARY_DIRS)
    # the RPATH to be used when installing
    list(APPEND CMAKE_INSTALL_RPATH "${LAPACK_LIBRARY_DIRS}")
else ()
    message(FATAL_ERROR "A LAPACK library is required but has not been found")
endif ()
if (LAPACKE_FOUND)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    if (LAPACKE_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS}")
    endif ()

    if (LAPACKE_LIBRARIES)
        if (LAPACKE_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${LAPACKE_LIBRARIES_DEP})
        else ()
            list(INSERT EXAGEOSTAT_DEP 0 ${LAPACKE_LIBRARIES})
        endif ()
    endif ()
else ()
    if (ECRC_VERBOSE_FIND_PACKAGE)
        if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
            if (NOT LAPACKE_lapacke.h_DIRS)
                Print_Find_Header_Status(lapacke lapacke.h)
            endif ()
            if (NOT LAPACKE_lapacke_LIBRARY)
                Print_Find_Library_Status(lapacke liblapacke)
            endif ()
        endif ()
    else ()
        message(WARNING "LAPACKE library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
                " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif ()
    message(FATAL_ERROR "A LAPACKE library is required but has not been found")
endif ()

# EXAGEOSTAT depends on GPU
#-------------------------
if (EXAGEOSTAT_USE_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DEXAGEOSTAT_USE_CUDA)
endif ()

# EXAGEOSTAT depends on MPI
#-------------------------
if (EXAGEOSTAT_USE_MPI)
    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if (NOT MPI_C_COMPILER)
        set(MPI_C_COMPILER mpicc)
    endif ()
    find_package(MPI REQUIRED)

    if (MPI_C_FOUND)
        message("-- ${Blue}Add definition EXAGEOSTAT_USE_MPI"
                " - Activate MPI in ExaGeoStat${ColourReset}")
        set(EXAGEOSTAT_USE_MPI ON)
        list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES})
        include_directories(${MPI_C_INCLUDE_PATH})
        # tests for intel mpi
        #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
        #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
        if (MPI_C_LINK_FLAGS)
            if(LAPACKE_LINKER_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${LAPACKE_LINKER_FLAGS}")
            endif()
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS}")
        endif ()
    endif (MPI_C_FOUND)

endif (EXAGEOSTAT_USE_MPI)

# EXAGEOSTAT depends on netcdf
# -------------------------------
if (EXAGEOSTAT_USE_NETCDF)
    message("")
    message("---------------------------------------- NETCDF")
    pkg_search_module(NETCDF netcdf)
    if (NOT NETCDF_FOUND)
        message("")
        message(STATUS "Can't find NetCDF using PkgConfig, Installing it instead.")
        include(ImportNetCDF)
        if(NOT NETCDF_INSTALLED)
            message(FATAL_ERROR "NetCDF Installation failed")
        endif()
        message(STATUS "Now trying to find netCDF")
    endif()

    pkg_search_module(NETCDF netcdf)
    if (NETCDF_FOUND)
        add_definitions(-DEXAGEOSTAT_USE_NETCDF)
        include_directories(${NETCDF_INCLUDE_DIRS})
        link_directories(${NETCDF_LIBRARY_DIRS})
        list(INSERT EXAGEOSTAT_DEP 0 ${NETCDF_LIBRARIES})
        message(STATUS "NETCDF ${NETCDF_VERSION} found")
    endif ()
    message(STATUS "NetCDF Done")
endif ()

# EXAGEOSTAT depends on NLOPT
# -------------------------------
# Check if NLOPT is available to build wrappers
message("")
message("---------------------------------------- NLOPT")
find_package(NLOPT 2.4.2 QUIET)
if (NOT NLOPT_FOUND)

    message("")
    message(STATUS "Can't find NLOPT using PkgConfig, Installing it ..")
    include(ImportNLOPT)
    if(NOT NLOPT_INSTALLED)
        message(FATAL_ERROR "NLOPT Installation failed")
    endif()

    message(STATUS "Now trying to find NLOPT")
endif()

#pkg_search_module(NETCDF netcdf REQUIRED)
find_package(NLOPT 2.4.2 QUIET REQUIRED)
if (NLOPT_FOUND)
    include_directories(${NLOPT_INCLUDE_DIRS})
    link_directories(${NLOPT_LIBRARY_DIRS})
    list(INSERT EXAGEOSTAT_DEP 0 ${NLOPT_LIBRARIES})
    message(STATUS "NLOPT ${NLOPT_VERSION} found")
else ()
    message("NLOPT is NOT found, skipping it")
    set(NLOPT OFF)
endif ()


# EXAGEOSTAT depends on STARSH
# -------------------------------
if (EXAGEOSTAT_USE_STARSH)
    message("")
    message("---------------------------------------- Stars-H")
    find_package(STARSH)
    if (NOT STARSH_FOUND)

        message("")
        message(STATUS "Can't find STARSH using PkgConfig, Installing it through HiCMA")
        include(ImportHiCMA)
        if(NOT HiCMA_INSTALLED)
            message(FATAL_ERROR "HiCMA Installation failed")
        endif()

        message(STATUS "Now trying to find STARSH")
    endif()
    find_package(STARSH REQUIRED)
    if (STARSH_FOUND)
        include_directories(${STARSH_INCLUDE_DIRS_DEP})
        if (STARSH_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARSH_LINKER_FLAGS}")
        endif ()
        if (STARSH_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${STARSH_LIBRARY_DIRS}")
        endif ()
        if (STARSH_LIBRARIES)
            # look for gsl
            find_library(_STARSH_LIB NAME starsh PATHS ${STARSH_LIBRARY_DIRS})
            if (_STARSH_LIB AND NOT "${STARSH_LIBRARIES_DEP}" MATCHES "gsl")
                execute_process(COMMAND nm ${_STARSH_LIB} COMMAND grep gsl RESULT_VARIABLE GSL_IN_STARSH)
                if (${GSL_IN_STARSH} EQUAL 0)
                    message(STATUS "STARSH depends on gsl. Adding it to dependency list")
                    if (STARSH_LIBRARIES_DEP)
                        list(APPEND STARSH_LIBRARIES_DEP "gsl")
                    else ()
                        list(APPEND STARSH_LIBRARIES "gsl")
                    endif ()
                endif ()
            endif ()
            # insert to dependencies
            if (STARSH_LIBRARIES_DEP)
                list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES_DEP})
            else ()
                list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES})
            endif ()
        endif ()
    endif()
    message(STATUS "StarsH Done")
endif ()

# EXAGEOSTAT depends on HICMA
# -------------------------------
if (EXAGEOSTAT_USE_HICMA)
    message("")
    message("---------------------------------------- Hicma")
    find_package(HICMA REQUIRED)
    if (HICMA_FOUND)
        include_directories(${HICMA_INCLUDE_DIRS_DEP})
        add_definitions(-DEXAGEOSTAT_USE_HICMA)
        if (HICMA_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${HICMA_LINKER_FLAGS}")
        endif ()
        if (HICMA_LIBRARY_DIRS)
            # The RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${HICMA_LIBRARY_DIRS}")
            link_directories(${HICMA_LIBRARY_DIRS})
        endif ()
        # Insert to dependencies
        if (HICMA_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${HICMA_LIBRARIES_DEP})
        else ()
            list(INSERT EXAGEOSTAT_DEP 0 ${HICMA_LIBRARIES})
        endif ()
    message(STATUS "Hicma Done")
    endif ()
endif()


# EXAGEOSTAT depends on a runtime
# -------------------------------
set(EXAGEOSTAT_STARPU_VERSION "1.2" CACHE STRING "STARPU version desired")

# create list of components in order to make a single call to find_package(starpu...)
if (NOT EXAGEOSTAT_SIMULATION)
    set(STARPU_COMPONENT_LIST "HWLOC")
else ()
    set(STARPU_COMPONENT_LIST "SIMGRID")
endif ()
if (EXAGEOSTAT_USE_MPI)
    list(APPEND STARPU_COMPONENT_LIST "MPI")
endif ()
if (EXAGEOSTAT_ENABLE_TRACING)
    list(APPEND STARPU_COMPONENT_LIST "FXT")
endif ()

find_package(STARPU ${EXAGEOSTAT_STARPU_VERSION} REQUIRED
        COMPONENTS ${STARPU_COMPONENT_LIST})

# Add definition and include_dir if found
if (STARPU_FOUND)
    message("-- ${Blue}Add definition EXAGEOSTAT_SCHED_STARPU"
            " - Activate StarPU in ExaGeoStat${ColourReset}")
    set(EXAGEOSTAT_SCHED_STARPU 1)
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    if (STARPU_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARPU_LINKER_FLAGS}")
    endif ()
    set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
    foreach (libdir ${STARPU_LIBRARY_DIRS_DEP})
        list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
    endforeach ()
    set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
    if (EXAGEOSTAT_USE_MPI)
        list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
        list(APPEND CMAKE_REQUIRED_FLAGS "${MPI_C_LINK_FLAGS}")
        list(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
    endif ()
    if (EXAGEOSTAT_SIMULATION)
        list(APPEND CMAKE_REQUIRED_FLAGS "-include" "starpu_simgrid_wrap.h")
    endif ()
    string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    check_function_exists(starpu_data_idle_prefetch_on_node STARPU_IDLE_PREFETCH_FOUND)
    if (STARPU_IDLE_PREFETCH_FOUND)
        set(HAVE_STARPU_IDLE_PREFETCH 1)
        message("-- ${Blue}Add definition HAVE_STARPU_IDLE_PREFETCH${ColourReset}")
    endif ()
    check_function_exists(starpu_iteration_push STARPU_ITERATION_PUSH_FOUND)
    if (STARPU_ITERATION_PUSH_FOUND)
        set(HAVE_STARPU_ITERATION_PUSH 1)
        message("-- ${Blue}Add definition HAVE_STARPU_ITERATION_PUSH${ColourReset}")
    endif ()
    check_function_exists(starpu_data_wont_use STARPU_DATA_WONT_USE_FOUND)
    if (STARPU_DATA_WONT_USE_FOUND)
        set(HAVE_STARPU_DATA_WONT_USE 1)
        message("-- ${Blue}Add definition HAVE_STARPU_DATA_WONT_USE${ColourReset}")
    endif ()
    check_function_exists(starpu_data_set_coordinates STARPU_DATA_SET_COORDINATES_FOUND)
    if (STARPU_DATA_SET_COORDINATES_FOUND)
        set(HAVE_STARPU_DATA_SET_COORDINATES 1)
        message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_COORDINATES${ColourReset}")
    endif ()
    check_function_exists(starpu_malloc_on_node_set_default_flags STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
    if (STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
        set(HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS 1)
        message("-- ${Blue}Add definition HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS${ColourReset}")
    endif ()
    if (EXAGEOSTAT_ENABLE_TRACING)
        # check if fxt profiling is accessible in starpu and activate it in exageostat
        unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
        check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
        if (STARPU_FXT_START_PROFILING_FOUND)
            message("-- ${Blue}Add definition HAVE_STARPU_FXT_PROFILING"
                    " - Activate FxT profiling through StarPU${ColourReset}")
            set(HAVE_STARPU_FXT_PROFILING 1)
        else ()
            message("-- ${Red}Looking for starpu with fxt"
                    " - starpu_fxt_start_profiling() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif ()
    endif ()
    if (EXAGEOSTAT_USE_MPI)
        # Check if a specific function exist
        unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
        check_function_exists(starpu_mpi_data_register_comm STARPU_MPI_DATA_REGISTER_FOUND)
        if (STARPU_MPI_DATA_REGISTER_FOUND)
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
                    " use of starpu_mpi_data_register() in ExaGeoStat with StarPU${ColourReset}")
            set(HAVE_STARPU_MPI_DATA_REGISTER 1)
        else ()
            message("-- ${Red}Looking for starpu with starpu_mpi_data_register"
                    " - starpu_mpi_data_register() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif ()
        unset(STARPU_MPI_COMM_RANK_FOUND CACHE)
        check_function_exists(starpu_mpi_comm_rank STARPU_MPI_COMM_RANK_FOUND)
        if (STARPU_MPI_COMM_RANK_FOUND)
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_RANK - Activate"
                    " use of starpu_mpi_comm_rank() in ExaGeoStat with StarPU${ColourReset}")
            set(HAVE_STARPU_MPI_COMM_RANK 1)
        else ()
            message("-- ${Red}Looking for starpu with starpu_mpi_comm_rank"
                    " - starpu_mpi_comm_rank() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif ()
        check_function_exists(starpu_mpi_cached_receive STARPU_MPI_CACHED_RECEIVE)
        if (STARPU_MPI_CACHED_RECEIVE)
            set(HAVE_STARPU_MPI_CACHED_RECEIVE 1)
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_CACHED_RECEIVE${ColourReset}")
        endif ()
    endif ()
    if (HWLOC_FOUND AND HWLOC_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${HWLOC_LIBRARY_DIRS}")
    endif ()
    if (FXT_FOUND AND FXT_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${FXT_LIBRARY_DIRS}")
    endif ()
    if (SIMGRID_FOUND AND SIMGRID_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${SIMGRID_LIBRARY_DIRS}")
    endif ()
    if (STARPU_FOUND AND STARPU_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS}")
    endif ()
    if (STARPU_LIBRARIES)
        if (STARPU_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${STARPU_LIBRARIES_DEP})
        else ()
            list(INSERT EXAGEOSTAT_DEP 0 ${STARPU_LIBRARIES})
        endif ()
    endif ()
endif (STARPU_FOUND)

# EXAGEOSTAT depends on CHAMELEON
# -------------------------------
set(EXAGEOSTAT_USE_CHAMELEON 1)  #Should be enabled
if (EXAGEOSTAT_USE_CHAMELEON)
    message("")
    message("---------------------------------------- Chameleon")
    find_package(CHAMELEON)
    if (NOT CHAMELEON_FOUND)
        message("")
        message(STATUS "Can't find CHAMELEON using PkgConfig, Installing it instead.")
        include(ImportChameleon)
        if(NOT CHAMELEON_INSTALLED)
            message(FATAL_ERROR "CHAMELEON Installation failed")
        endif()

        message(STATUS "Now trying to find Chameleon")
    endif()
    find_package(CHAMELEON REQUIRED)
    if(CHAMELEON_FOUND)
        include_directories(AFTER ${CHAMELEON_INCLUDE_DIRS_DEP})
        include_directories(AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)
        if (CHAMELEON_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${CHAMELEON_LINKER_FLAGS}")
        endif ()
        if (CHAMELEON_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${CHAMELEON_LIBRARY_DIRS}")
        endif ()
        if (CHAMELEON_LIBRARIES)
            if (CHAMELEON_LIBRARIES_DEP)
                list(INSERT EXAGEOSTAT_DEP 0 ${CHAMELEON_LIBRARIES_DEP})
            else ()
                list(INSERT EXAGEOSTAT_DEP 0 ${CHAMELEON_LIBRARIES})
            endif ()
        endif ()
    endif()
    message(STATUS "Chameleon Done")
endif ()

list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# Fix a problem on Mac OS X when building shared libraries
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif ()
#------------------------------------------------------------------------------

if (EXAGEOSTAT_SCHED_STARPU)
    link_directories(${STARPU_LIBRARY_DIRS_DEP})
elseif (EXAGEOSTAT_SCHED_QUARK)
    link_directories(${QUARK_LIBRARY_DIRS})
endif ()

# Save extra dependencies (all required links)
list(APPEND EXAGEOSTAT_DEP ${EXTRA_LIBRARIES})
list(REMOVE_DUPLICATES EXAGEOSTAT_DEP) # WARNING: is it safe, respect order?

# Check for the subdirectories
# ----------------------------
# include headers
# ---------------
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/misc/include)
include_directories(${CMAKE_SOURCE_DIR}/exageostat_exact/src/include)
include_directories(${CMAKE_SOURCE_DIR}/exageostat_approx/src/include)
include_directories(${CMAKE_SOURCE_DIR}/exageostat_exact/cpu_core/include)
include_directories(${CMAKE_SOURCE_DIR}/exageostat_exact/runtime/starpu/include)
include_directories(${CMAKE_SOURCE_DIR}/exageostat_approx/runtime/starpu/include)

set(EXAGEOSTAT_CHAMELEON_SRC "")
set(EXAGEOSTAT_SRC "")

set(EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/compute/MLE.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zmdet.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zzcpy.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zcmg.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zcmg_non_stat.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zcmg_nuggets.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zprint.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_zdotp.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_sdconv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_dsconv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_slag2d.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_dlag2s.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_g_to_ng.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_ng_loglike.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_ng_transform.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/cpu_core/compute/core_ztrace.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zdotp.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zmdet.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_ztrace.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zmse.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zmse_bivariate.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zmse_trivariate.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zcmg.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zzcpy.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zgemv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/runtime/starpu/codelets/codelet_dcmg_diag.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/runtime/starpu/codelets/codelet_dpotrf_diag.c
       	${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_sdconv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_exageostat_sgemm.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_exageostat_strsm.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dsconv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_sdmat_reg.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zmloe_mmom.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_stride_vec.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_ng_transform.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_ng_loglike.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_g_to_ng.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_zcorr.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_tristride_vec.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/control/exact_runtime_descriptor.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/compute/pdpotrf_diag.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/compute/dpotrf_diag.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/compute/MLE_approx.c
        ${CMAKE_CURRENT_SOURCE_DIR}/misc/compute/MLE_misc.c
        ${CMAKE_CURRENT_SOURCE_DIR}/misc/compute/flat_file.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/MLE_dexact.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/MLE_sexact.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/MLE_ng_dexact.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/MLE_sdexact.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/psdportf.c
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/sdpotrf.c
        ${CMAKE_CURRENT_SOURCE_DIR}/r-wrappers/compute/rwrappers.c
        )

if (EXAGEOSTAT_USE_HICMA)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/compute/MLE_lr.c)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/compute/MLE_ng_lr.c)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/runtime/starpu/codelets/codelet_hicma_dmdet.c)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/runtime/starpu/codelets/codelet_hicma_ng_transform.c)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/runtime/starpu/codelets/codelet_hicma_ng_loglike.c)
endif ()

if (EXAGEOSTAT_USE_CUDA)
    message(WARNING "USING CUDA")
    list(APPEND EXAGEOSTAT_SRC exageostat_exact/cuda_core/compute/cuda_conv.cu)
    list(APPEND EXAGEOSTAT_SRC exageostat_exact/cuda_core/compute/cuda_zcmg.cu)
endif ()

if (EXAGEOSTAT_USE_NETCDF)
    list(APPEND EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/misc/compute/nc_file.c)
endif ()

set(INSTALL_HEADERS ${CMAKE_SOURCE_DIR}/src/include/MLE.h
        ${CMAKE_SOURCE_DIR}/include/context.h
        ${CMAKE_SOURCE_DIR}/include/descriptor.h
        ${CMAKE_SOURCE_DIR}/include/auxiliary.h
        ${CMAKE_SOURCE_DIR}/include/chameleon_starpu.h
        ${CMAKE_SOURCE_DIR}/include/common.h
        ${CMAKE_SOURCE_DIR}/include/global.h
        ${CMAKE_SOURCE_DIR}/include/compute_d.h
        ${CMAKE_SOURCE_DIR}/include/async.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include/MLE_exact.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include/MLE_exact_s.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include/MLE_sdexact.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include/MLE_ng_exact.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/runtime/starpu/include/starpu_exageostat.h
        ${CMAKE_SOURCE_DIR}/exageostat_exact/cpu_core/include/exageostatcore.h
        ${CMAKE_SOURCE_DIR}/exageostat_approx/src/include/MLE_approx.h
        ${CMAKE_SOURCE_DIR}/exageostat_approx/runtime/starpu/include/starpu_exageostat_approx.h
        ${CMAKE_SOURCE_DIR}/misc/include/MLE_misc.h
        ${CMAKE_SOURCE_DIR}/misc/include/flat_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/r-wrappers/include/rwrappers.h
        ${CMAKE_SOURCE_DIR}/include/flops.h
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_approx/src/include/diag.h
        ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/include/mixed_prec.h
        )

if (EXAGEOSTAT_USE_HICMA)
    list(APPEND INSTALL_HEADERS ${CMAKE_SOURCE_DIR}/exageostat_approx/src/include/MLE_lr.h)
endif ()
if (EXAGEOSTAT_USE_CUDA)
    list(APPEND INSTALL_HEADERS ${CMAKE_SOURCE_DIR}/exageostat_exact/cuda_core/include/exageostatcudacore.h)
endif ()
if (EXAGEOSTAT_USE_NETCDF)
    list(APPEND INSTALL_HEADERS ${CMAKE_SOURCE_DIR}/misc/include/nc_file.h)
endif ()

link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${STARPU_LIBRARY_DIRS})
link_directories(${STARPU_SHM_STATIC_LIBRARY_DIRS})
link_directories(${CHAMELEON_LIBRARY_DIRS})
link_directories(${HICMA_LIBRARY_DIRS})

if (EXAGEOSTAT_USE_CUDA)
    CUDA_add_library(exageostat ${EXAGEOSTAT_SRC} ${EXAGEOSTAT_CHAMELEON_SRC})
    set_target_properties(exageostat PROPERTIES LINKER_LANGUAGE C)
else ()
    add_library(exageostat ${EXAGEOSTAT_SRC} ${EXAGEOSTAT_CHAMELEON_SRC})
endif ()

set_target_properties(exageostat PROPERTIES PUBLIC_HEADER "${INSTALL_HEADERS}")
target_link_libraries(exageostat
        ${EXAGEOSTAT_DEP}
        )
install(TARGETS exageostat
        DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/
        )

if (EXAGEOSTAT_EXAMPLES)
    enable_testing()
    include(CTest)
    add_subdirectory(examples)
endif ()
###############################################################################
# Documentation #
#################
option(EXAGEOSTAT_ENABLE_DOCS "Build documentation in docs directory" ON)
if (EXAGEOSTAT_ENABLE_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        add_subdirectory("docs")
    else ()
        message(STATUS "Doxygen NOT found, skipping documentation")
    endif ()
endif ()
#------------------------------------------------------------------------------
###############################################################################
# Config files (pkg_config) #
#############################
# Create file chameleon_starpu.pc
# ---------------------------
generate_pkgconfig_file()
#------------------------------------------------------------------------------
option(EXAGEOSTAT_PACKAGE "Enable a packaging system for distribution" OFF)
if (EXAGEOSTAT_PACKAGE)
    ###############################################################################
    # Release source #
    ##################
    set(CPACK_SOURCE_GENERATOR "TGZ")
    set(CPACK_GENERATOR "TGZ")
    set(CPACK_PACKAGE_NAME "exageostat")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ExaGeoStat is a parallel high performance unified framework for geostatistics on manycore systems. Its abbreviation stands for 'Exascale Geostatistics'.")
    set(CPACK_PACKAGE_VERSION "${EXAGEOSTAT_VERSION_MAJOR}.${EXAGEOSTAT_VERSION_MINOR}.${EXAGEOSTAT_VERSION_PATCH}")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-src")
    set(CPACK_PACKAGE_VENDOR "KAUST")
    set(CPACK_PACKAGE_CONTACT "sameh.abdulah@kaust.edu.sa")
    set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_SOURCE_IGNORE_FILES "plasma-conversion;build;.cproject;.settings;.dir-locals.el;.project;.pydevproject;.svn;.git;.gitmodules;.gitlab-ci.yml;Jenkinsfile;jenkins-scripts;chameleon-repo")
    include(CPack)
endif ()
###############################################################################
# Print Options #
#################
get_directory_property(EXAGEOSTAT_DEFINITIONS_LIST DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
include(PrintOpts)

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
            " You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
else ()
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
endif ()
###
### END CMakeLists.txt
###
