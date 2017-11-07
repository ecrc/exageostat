###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2016 Bordemisc INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordemisc. All rights reserved.
# @copyright (c) 2017      King Abdullah University of Science and Technology (KAUST)
#
###
#
#  @file CMakeLists.txt
#
#  @project ExaGeoStat
#  ExaGeoStat is a software package provided by:
#     King Abdullah University of Science and Technology
#
#  @version 0.1.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Eduardo Gonzalez
#  @date 13-07-2017
#
###
cmake_minimum_required(VERSION 3.2.3)

# directly make an error if in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.\n"
   "Please create a build directory first and execute cmake configuration from "
   "this directory. Example: mkdir build && cd build && cmake ..")
endif()

project(EXAGEOSTAT C )

# set project version number
set(EXAGEOSTAT_VERSION_MAJOR 0)
set(EXAGEOSTAT_VERSION_MINOR 1)
set(EXAGEOSTAT_VERSION_PATCH 0)

set(MORSE_CMAKE_DIR "" CACHE PATH "Directory of MORSE CMake modules, can be external to the project")

## CMAKE MODULES :: ECRC
## REQUIRED FOR TESTS TO LINK LIBRARIES
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules" )
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
        message(FATAL_ERROR "ECRC CMake modules were not found.\n"
                            "We tried: 'git submodule init && git submodule update' and resulted in error" )
    endif()
endif()

## ECRC INITIALIZATION
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules")
set(ECRC_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules )

include(EcrcInit)
include(GenPkgConfig)

# # Get Git working directory status and HEAD hash
# include(GetGitRevisionDescription)
# get_git_head_revision(GIT_REFSPEC GIT_SHA1)
# git_local_changes(GIT_DIRTY --ignore-submodules)
# if( ${GIT_DIRTY} STREQUAL "CLEAN" )
#     set( GIT_COMMIT_STATUS "${GIT_SHA1}" )
# else()
#     set( GIT_COMMIT_STATUS "${GIT_SHA1}-${GIT_DIRTY}" )
# endif()
# add_definitions( "-DVERSION='${GIT_COMMIT_STATUS}'" )
set( EXAGEOSTAT_VERSION ${EXAGEOSTAT_VERSION_MAJOR}.${EXAGEOSTAT_VERSION_MINOR}.${EXAGEOSTAT_VERSION_PATCH} )
add_definitions( "-DVERSION=${EXAGEOSTAT_VERSION}" )



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
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define precision supported by EXAGEOSTAT
# -----------------------------------------
set( RP_EXAGEOSTAT_DICTIONNARY ${PROJECT_SOURCE_DIR}/cmake_modules/morse_cmake/modules/precision_generator/subs.py )
set( RP_EXAGEOSTAT_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

# Options to select the runtime
# -----------------------------

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
endif()
# Set default to StarPU if nothing specific is required by the user
if ( NOT EXAGEOSTAT_SCHED_STARPU AND NOT EXAGEOSTAT_SCHED_QUARK )
    set(EXAGEOSTAT_SCHED_STARPU ON)
endif()
if (EXAGEOSTAT_SCHED_STARPU)
    message("-- ${BoldGreen}EXAGEOSTAT_SCHED_STARPU is set to ON: EXAGEOSTAT uses StarPU runtime\n"
            "   To use EXAGEOSTAT with Quark  runtime: set EXAGEOSTAT_SCHED_QUARK  to ON\n"
            "   (EXAGEOSTAT_SCHED_STARPU will be disabled)${ColourReset}")
elseif(EXAGEOSTAT_SCHED_QUARK)
    message("-- ${BoldGreen}EXAGEOSTAT_SCHED_QUARK is set to ON: EXAGEOSTAT uses Quark runtime\n"
            "   To use EXAGEOSTAT with StarPU runtime: set EXAGEOSTAT_SCHED_STARPU to ON\n"
            "   (EXAGEOSTAT_SCHED_QUARK will be disabled)${ColourReset}")
endif()


# Check that one, and only one, SCHED option is set to ON
# count number of runtime sets to ON
math(EXPR number_of_active_runtime 0)
foreach (runtime ${EXAGEOSTAT_SCHED_list})
    if (EXAGEOSTAT_SCHED_${runtime})
        math(EXPR number_of_active_runtime "${number_of_active_runtime}+1")
    endif()
endforeach()
if (NOT number_of_active_runtime STREQUAL 1)
    message(FATAL_ERROR
            "Number of active runtime is ${number_of_active_runtime}, "
            "the user should activate one (and only one) runtime. ")
endif()


# Warning if not coherent options
if (EXAGEOSTAT_USE_MPI AND NOT EXAGEOSTAT_SCHED_STARPU)
    message(WARNING "You have activated MPI but EXAGEOSTAT_SCHED_STARPU is OFF.\n"
    "You should use StarPU Runtime system if you want to benefit from MPI.\n"
    "Use -DEXAGEOSTAT_SCHED_QUARK=OFF -DEXAGEOSTAT_SCHED_STARPU=ON at cmake \n"
    "configure to do so.")
endif()

# Additional options
# ------------------

# Enable the distributed interface (allowed only when StarPU is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(EXAGEOSTAT_USE_MPI
                    "Enable distributed memory through MPI" OFF
                    "EXAGEOSTAT_SCHED_STARPU" OFF)
if (NOT EXAGEOSTAT_USE_MPI)
    message("-- ${BoldGreen}EXAGEOSTAT_USE_MPI is set to OFF, turn it ON to use MPI (only with StarPU)${ColourReset}")
endif()

# Use intermediate variable since cmake_dependent_option doesn't have OR conditions
set(EXAGEOSTAT_ENABLE_MPI OFF CACHE INTERNAL "Tells if MPI might be supported by the runtime")
if ( EXAGEOSTAT_SCHED_STARPU )
     set(EXAGEOSTAT_ENABLE_MPI ON FORCE)
endif()

# Initially we need to generate files for different precisions
# TODO: use this option to avoid generate multiple precisions each time we launch cmake
#option(EXAGEOSTAT_GEN_PREC "Generate source files precisions" ON)
#------------------------------------------------------------------------------


###############################################################################
# Look for dependencies #
#########################
set(EXAGEOSTAT_DEP "")

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if( THREADS_FOUND )
    list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
set(M_LIBRARIES "")
if(UNIX OR WIN32)
    find_library(
        M_m_LIBRARY
        NAMES m
        )
    mark_as_advanced(M_m_LIBRARY)
    if (M_m_LIBRARY)
        list(APPEND M_LIBRARIES "${M_m_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${M_m_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find libm on your system."
            " Are you sure to a have a C compiler installed?")
    endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
set(RT_LIBRARIES "")
if(UNIX AND NOT APPLE)
    find_library(
        RT_rt_LIBRARY
        NAMES rt
        )
    mark_as_advanced(RT_rt_LIBRARY)
    if (RT_rt_LIBRARY)
        list(APPEND RT_LIBRARIES "${RT_rt_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${RT_rt_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
endif()



# EXAGEOSTAT depends on CBLAS
#---------------------------
find_package(CBLAS COMPONENTS BLASEXT)
if(BLAS_FOUND)
    if (BLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
    endif()
    if(BLAS_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
    endif()
else()
    message(FATAL_ERROR "BLAS library has not been found")
endif()
if(CBLAS_FOUND)
    include_directories(${CBLAS_INCLUDE_DIRS})
    if(CBLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS}")
    endif()
    if (CBLAS_LIBRARIES)
        if (CBLAS_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES_DEP})
        else()
            list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES})
        endif()
    endif()
else()
    if(ECRC_VERBOSE_FIND_PACKAGE)
        if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
            if (NOT CBLAS_cblas.h_DIRS)
                Print_Find_Header_Status(cblas cblas.h)
            endif ()
            if (NOT CBLAS_cblas_LIBRARY)
                Print_Find_Library_Status(cblas libcblas)
            endif ()
        endif()
    else()
        message(WARNING "CBLAS library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
        " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif()
    message(FATAL_ERROR "A CBLAS library is required but has not been found")
endif()

# EXAGEOSTAT depends on LAPACKE
#-----------------------------
# standalone version of lapacke seems useless for now
# let the comment in case we meet some problems of non existing lapacke
# functions in lapack library such as mkl, acml, ...
#set(LAPACKE_STANDALONE TRUE)
find_package(LAPACKE COMPONENTS LAPACKEXT)
if(LAPACK_FOUND AND LAPACK_LIBRARY_DIRS)
    # the RPATH to be used when installing
    list(APPEND CMAKE_INSTALL_RPATH "${LAPACK_LIBRARY_DIRS}")
else()
    message(FATAL_ERROR "A LAPACK library is required but has not been found")
endif()
if(LAPACKE_FOUND)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    if(LAPACKE_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS}")
    endif()
    if(LAPACKE_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${LAPACKE_LINKER_FLAGS}")
    endif()
    if (LAPACKE_LIBRARIES)
        if (LAPACKE_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${LAPACKE_LIBRARIES_DEP})
        else()
            list(INSERT EXAGEOSTAT_DEP 0 ${LAPACKE_LIBRARIES})
        endif()
    endif()
else()
    if(ECRC_VERBOSE_FIND_PACKAGE)
        if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
            if (NOT LAPACKE_lapacke.h_DIRS)
                Print_Find_Header_Status(lapacke lapacke.h)
            endif ()
            if (NOT LAPACKE_lapacke_LIBRARY)
                Print_Find_Library_Status(lapacke liblapacke)
            endif ()
        endif()
    else()
        message(WARNING "LAPACKE library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
        " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif()
    message(FATAL_ERROR "A LAPACKE library is required but has not been found")
endif()


# EXAGEOSTAT depends on MPI
#-------------------------
if (EXAGEOSTAT_USE_MPI)

    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if(NOT MPI_C_COMPILER)
        set(MPI_C_COMPILER mpicc)
    endif()
    find_package(MPI REQUIRED)

    if (MPI_C_FOUND)
        message("-- ${Blue}Add definition EXAGEOSTAT_USE_MPI"
        " - Activate MPI in ExaGeoStat${ColourReset}")
        set(EXAGEOSTAT_USE_MPI ON)
        list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
        include_directories( ${MPI_C_INCLUDE_PATH} )
        # tests for intel mpi
        #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
        #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
        if(MPI_C_LINK_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS}")
        endif()
    endif (MPI_C_FOUND)

endif (EXAGEOSTAT_USE_MPI)


# EXAGEOSTAT depends on STARSH
# -------------------------------
if( EXAGEOSTAT_USE_STARSH )
find_package( STARSH REQUIRED)
if( STARSH_FOUND )
    include_directories(${STARSH_INCLUDE_DIRS_DEP})
    if(STARSH_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARSH_LINKER_FLAGS}")
    endif()
    if(STARSH_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${STARSH_LIBRARY_DIRS}")
    endif()
    if (STARSH_LIBRARIES)
        # look for gsl
        find_library( _STARSH_LIB NAME starsh PATHS ${STARSH_LIBRARY_DIRS} )
        if( _STARSH_LIB AND NOT "${STARSH_LIBRARIES_DEP}" MATCHES "gsl" )
            execute_process( COMMAND nm ${_STARSH_LIB} COMMAND grep gsl RESULT_VARIABLE GSL_IN_STARSH )
            if ( ${GSL_IN_STARSH} EQUAL 0)
                message( STATUS "STARSH depends on gsl. Adding it to dependency list")
                if (STARSH_LIBRARIES_DEP)
                    list( APPEND STARSH_LIBRARIES_DEP "gsl")
                else()
                    list( APPEND STARSH_LIBRARIES "gsl")
                endif()
            endif()
    endif()
    # insert to dependencies
        if (STARSH_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES_DEP})
        else()
            list(INSERT EXAGEOSTAT_DEP 0 ${STARSH_LIBRARIES})
        endif()
    endif()
endif()
endif()


# EXAGEOSTAT depends on a runtime
# -------------------------------

#set(EXAGEOSTAT_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")
set(EXAGEOSTAT_STARPU_VERSION "1.2" CACHE STRING "STARPU version desired")

# create list of components in order to make a single call to find_package(starpu...)
if(NOT EXAGEOSTAT_SIMULATION)
    set(STARPU_COMPONENT_LIST "HWLOC")
else()
    set(STARPU_COMPONENT_LIST "SIMGRID")
endif()
if(EXAGEOSTAT_USE_MPI)
    list(APPEND STARPU_COMPONENT_LIST "MPI")
endif()
if(EXAGEOSTAT_ENABLE_TRACING)
    list(APPEND STARPU_COMPONENT_LIST "FXT")
endif()

find_package(STARPU ${EXAGEOSTAT_STARPU_VERSION} REQUIRED
             COMPONENTS ${STARPU_COMPONENT_LIST})

# Add definition and include_dir if found
if ( STARPU_FOUND )
    message("-- ${Blue}Add definition EXAGEOSTAT_SCHED_STARPU"
    " - Activate StarPU in ExaGeoStat${ColourReset}")
    set(EXAGEOSTAT_SCHED_STARPU 1)
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    if(STARPU_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARPU_LINKER_FLAGS}")
    endif()
    set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
    foreach(libdir ${STARPU_LIBRARY_DIRS_DEP})
        list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
    endforeach()
    set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
    if (EXAGEOSTAT_USE_MPI)
        list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
        list(APPEND CMAKE_REQUIRED_FLAGS "${MPI_C_LINK_FLAGS}")
        list(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
    endif()
    if (EXAGEOSTAT_SIMULATION)
        list(APPEND CMAKE_REQUIRED_FLAGS "-include" "starpu_simgrid_wrap.h")
    endif()
    string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    check_function_exists(starpu_data_idle_prefetch_on_node STARPU_IDLE_PREFETCH_FOUND)
    if ( STARPU_IDLE_PREFETCH_FOUND )
        set(HAVE_STARPU_IDLE_PREFETCH 1)
        message("-- ${Blue}Add definition HAVE_STARPU_IDLE_PREFETCH${ColourReset}")
    endif()
    check_function_exists(starpu_iteration_push STARPU_ITERATION_PUSH_FOUND)
    if ( STARPU_ITERATION_PUSH_FOUND )
        set(HAVE_STARPU_ITERATION_PUSH 1)
        message("-- ${Blue}Add definition HAVE_STARPU_ITERATION_PUSH${ColourReset}")
    endif()
    check_function_exists(starpu_data_wont_use STARPU_DATA_WONT_USE_FOUND)
    if ( STARPU_DATA_WONT_USE_FOUND )
        set(HAVE_STARPU_DATA_WONT_USE 1)
        message("-- ${Blue}Add definition HAVE_STARPU_DATA_WONT_USE${ColourReset}")
    endif()
    check_function_exists(starpu_data_set_coordinates STARPU_DATA_SET_COORDINATES_FOUND)
    if ( STARPU_DATA_SET_COORDINATES_FOUND )
        set(HAVE_STARPU_DATA_SET_COORDINATES 1)
        message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_COORDINATES${ColourReset}")
    endif()
    check_function_exists(starpu_malloc_on_node_set_default_flags STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
    if ( STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS )
        set(HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS 1)
        message("-- ${Blue}Add definition HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS${ColourReset}")
    endif()
    if(EXAGEOSTAT_ENABLE_TRACING)
        # check if fxt profiling is accessible in starpu and activate it in exageostat
        unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
        check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
        if ( STARPU_FXT_START_PROFILING_FOUND )
            message("-- ${Blue}Add definition HAVE_STARPU_FXT_PROFILING"
            " - Activate FxT profiling through StarPU${ColourReset}")
            set(HAVE_STARPU_FXT_PROFILING 1)
        else()
            message("-- ${Red}Looking for starpu with fxt"
            " - starpu_fxt_start_profiling() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif()
    endif()
    if (EXAGEOSTAT_USE_MPI)
        # Check if a specific function exist
        unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
        check_function_exists(starpu_mpi_data_register_comm STARPU_MPI_DATA_REGISTER_FOUND)
        if ( STARPU_MPI_DATA_REGISTER_FOUND )
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
            " use of starpu_mpi_data_register() in ExaGeoStat with StarPU${ColourReset}")
            set(HAVE_STARPU_MPI_DATA_REGISTER 1)
        else()
            message("-- ${Red}Looking for starpu with starpu_mpi_data_register"
            " - starpu_mpi_data_register() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif()
        unset(STARPU_MPI_COMM_RANK_FOUND CACHE)
        check_function_exists(starpu_mpi_comm_rank STARPU_MPI_COMM_RANK_FOUND)
        if ( STARPU_MPI_COMM_RANK_FOUND )
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_RANK - Activate"
            " use of starpu_mpi_comm_rank() in ExaGeoStat with StarPU${ColourReset}")
            set(HAVE_STARPU_MPI_COMM_RANK 1)
        else()
            message("-- ${Red}Looking for starpu with starpu_mpi_comm_rank"
            " - starpu_mpi_comm_rank() test fails in StarPU${ColourReset}")
            message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
        endif()
     check_function_exists(starpu_mpi_cached_receive STARPU_MPI_CACHED_RECEIVE)
     if ( STARPU_MPI_CACHED_RECEIVE )
      set(HAVE_STARPU_MPI_CACHED_RECEIVE 1)
      message("-- ${Blue}Add definition HAVE_STARPU_MPI_CACHED_RECEIVE${ColourReset}")
     endif()
    endif()
    if(HWLOC_FOUND AND HWLOC_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${HWLOC_LIBRARY_DIRS}")
    endif()
    if(FXT_FOUND AND FXT_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${FXT_LIBRARY_DIRS}")
    endif()
    if(SIMGRID_FOUND AND SIMGRID_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${SIMGRID_LIBRARY_DIRS}")
    endif()
    if(STARPU_FOUND AND STARPU_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS}")
    endif()
    if (STARPU_LIBRARIES)
        if (STARPU_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${STARPU_LIBRARIES_DEP})
        else()
            list(INSERT EXAGEOSTAT_DEP 0 ${STARPU_LIBRARIES})
        endif()
    endif()
endif ( STARPU_FOUND )

# EXAGEOSTAT depends on CHAMELEON
# -------------------------------
find_package( CHAMELEON REQUIRED )
if( CHAMELEON_FOUND )
    include_directories( AFTER ${CHAMELEON_INCLUDE_DIRS_DEP})
    include_directories( AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)
    if(CHAMELEON_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${CHAMELEON_LINKER_FLAGS}")
    endif()
    if(CHAMELEON_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${CHAMELEON_LIBRARY_DIRS}")
    endif()
    if (CHAMELEON_LIBRARIES)
        if (CHAMELEON_LIBRARIES_DEP)
            list(INSERT EXAGEOSTAT_DEP 0 ${CHAMELEON_LIBRARIES_DEP})
        else()
            list(INSERT EXAGEOSTAT_DEP 0 ${CHAMELEON_LIBRARIES})
        endif()
    endif()
endif()


# Check if NLOPT is available to build wrappers
pkg_search_module(NLOPT nlopt REQUIRED)
if(NLOPT_FOUND)
    include_directories(${NLOPT_INCLUDE_DIRS})
    link_directories(${NLOPT_LIBRARY_DIRS})
    list(INSERT EXAGEOSTAT_DEP 0 ${NLOPT_LIBRARIES})
#     add_definitions("-DNLOPT")
    message(STATUS "NLOPT ${NLOPT_VERSION} found")
else()
    message("NLOPT is NOT found, skipping it")
    set(NLOPT OFF)
endif()

# GSL
find_package(GSL)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
    list(APPEND EXAGEOSTAT_DEP ${GSL_LIBRARIES})
endif()

list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# Fix a problem on Mac OS X when building shared libraries
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif()

# # Add define for Fortran Mangling (should be defined somewhere else)
# # ------------------------------------------------------------------
# message("-- ${Blue}Add definition ADD_"
#         " - For Fortran mangling${ColourReset}")
# set(ADD_ 1)

#------------------------------------------------------------------------------
if(EXAGEOSTAT_SCHED_STARPU)
    link_directories(${STARPU_LIBRARY_DIRS_DEP})
elseif(EXAGEOSTAT_SCHED_QUARK)
    link_directories(${QUARK_LIBRARY_DIRS})
endif()

# Save extra dependencies (all required links)
list(APPEND EXAGEOSTAT_DEP ${EXTRA_LIBRARIES})
list(REMOVE_DUPLICATES EXAGEOSTAT_DEP) # WARNING: is it safe, respect order?

# Check for the subdirectories
# ----------------------------

# include headers
# ---------------
include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${CMAKE_SOURCE_DIR}/misc/include )
include_directories( ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include )
include_directories( ${CMAKE_SOURCE_DIR}/exageostat_exact/runtime/starpu/include )


set( EXAGEOSTAT_CHAMELEON_SRC "" )

set( EXAGEOSTAT_SRC "" )

add_definitions( "-DEXAGEOSTAT_COPY_DIAG -DSTARPU_USE_DEPRECATED_API -DSTARPU_USE_DEPRECATED_ONE_ZERO_API" )

set( EXAGEOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/compute/MLE.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/core/compute/core_dmdet.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/core/compute/core_dzcpy.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/core/compute/core_dcmg.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/core/compute/core_ddotp.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_ddotp.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dmdet.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dmse.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dcmg.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dzcpy.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/runtime/starpu/codelets/codelet_dgemv.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/misc/compute/MLE_misc.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/exageostat_exact/src/compute/MLE_exact.c
  )
set( INSTALL_HEADERS ${CMAKE_SOURCE_DIR}/src/include/MLE.h
                     ${CMAKE_SOURCE_DIR}/include/context.h
                     ${CMAKE_SOURCE_DIR}/include/morse_starpu.h
                     ${CMAKE_SOURCE_DIR}/exageostat_exact/src/include/MLE_exact.h
                     ${CMAKE_SOURCE_DIR}/exageostat_exact/runtime/starpu/include/starpu_exageostat.h
                     ${CMAKE_SOURCE_DIR}/exageostat_exact/core/include/exageostatcore.h
                     ${CMAKE_SOURCE_DIR}/misc/include/MLE_misc.h
                     ${CMAKE_SOURCE_DIR}/include/flops.h
                     )


link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${STARPU_LIBRARY_DIRS})
link_directories(${CHAMELEON_LIBRARY_DIRS})

add_library( exageostat ${EXAGEOSTAT_SRC} ${EXAGEOSTAT_CHAMELEON_SRC} )
set_target_properties( exageostat PROPERTIES PUBLIC_HEADER "${INSTALL_HEADERS}")
target_link_libraries( exageostat
    ${EXAGEOSTAT_DEP}
    )
install( TARGETS exageostat
         DESTINATION lib/
         PUBLIC_HEADER DESTINATION include/ )

add_subdirectory(examples)


###############################################################################
# Documentation #
#################
option(EXAGEOSTAT_ENABLE_DOCS "Build documentation in docs directory" ON)
if(EXAGEOSTAT_ENABLE_DOCS)
    find_package(Doxygen)
    find_package(LATEX COMPONENTS PDFLATEX)
    if(DOXYGEN_FOUND)
        if( LATEX_FOUND )
            add_subdirectory("docs")
        else()
            message(STATUS "Latex NOT found, skipping documentation")
        endif()
    else()
        message(STATUS "Doxygen NOT found, skipping documentation")
    endif()
endif()
#------------------------------------------------------------------------------


###############################################################################
# Config files (pkg_config) #
#############################

# Create file morse_starpu.pc
# ---------------------------
generate_pkgconfig_file()
#------------------------------------------------------------------------------

option( EXAGEOSTAT_PACKAGE "Enable a packaging system for distribution" OFF )
if( EXAGEOSTAT_PACKAGE )
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
#    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
    set(CPACK_SOURCE_IGNORE_FILES "plasma-conversion;build;.cproject;.settings;.dir-locals.el;.project;.pydevproject;.svn;.git;.gitmodules;.gitlab-ci.yml;Jenkinsfile;jenkins-scripts;chameleon-repo")
    include(CPack)
endif()
###############################################################################
# Print Options #
#################
get_directory_property( EXAGEOSTAT_DEFINITIONS_LIST DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
include(PrintOpts)

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
    " You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
else()
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
endif()

###
### END CMakeLists.txt
###
