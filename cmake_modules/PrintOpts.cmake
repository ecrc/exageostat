###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014, 2016 Bordemisc INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordemisc. All rights reserved.
# @copyright (c) 2017, King Abdullah University of Science and Technology. All rights reserved.
#
###
#
#  @file PrintOpts.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordemisc - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @author Florent Pruvost
#  @author Eduardo Gonzalez
#  @date 14-08-2017
#
###
set(dep_message "\nConfiguration of ExaGeoStat:\n"
        "       BUILDNAME ...........: ${BUILDNAME}\n"
        "       SITE ................: ${SITE}\n"
        "\n"
        "       Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
        "       Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n")
if(EXAGEOSTAT_USE_MPI)
  set(dep_message "${dep_message}"
  "       Compiler: MPI .......: ${MPI_C_COMPILER}\n"
  "       compiler flags ......: ${MPI_C_COMPILE_FLAGS}\n")
endif()
set(dep_message "${dep_message}"
"       Linker: .............: ${CMAKE_LINKER}\n"
"\n"
"       Build type ..........: ${CMAKE_BUILD_TYPE}\n"
"       Build shared ........: ${BUILD_SHARED_LIBS}\n"
"       CFlags ..............: ${CMAKE_C_FLAGS}\n"
"       LDFlags .............: ${CMAKE_EXE_LINKER_FLAGS}\n"
"\n"
"       Implementation paradigm\n"
"       CUDA ................: ${EXAGEOSTAT_USE_CUDA}\n"
"       MPI .................: ${EXAGEOSTAT_USE_MPI}\n"
"\n"
"       Runtime specific\n"
"       QUARK ...............: ${EXAGEOSTAT_SCHED_QUARK}\n"
"       PARSEC ..............: ${EXAGEOSTAT_SCHED_PARSEC}\n"
"       StarPU ..............: ${EXAGEOSTAT_SCHED_STARPU}\n"
"\n"
"       Kernels specific\n"
"       BLAS ................: ${BLAS_VENDOR_FOUND}\n"
"       LAPACK...............: ${LAPACK_VENDOR_FOUND}\n"
#"\n"
#"       Trace ...............: ${EXAGEOSTAT_ENABLE_TRACING}\n"
#"       Simulation mode .....: ${EXAGEOSTAT_SIMULATION}\n"
#"\n"
#"       Binaries to build\n"
#"       documentation ........: ${EXAGEOSTAT_ENABLE_DOCS}\n"
#"       example ..............: ${EXAGEOSTAT_ENABLE_EXAMPLE}\n"
#"       testing ..............: ${EXAGEOSTAT_ENABLE_TESTING}\n"
#"       timing ...............: ${EXAGEOSTAT_ENABLE_TIMING}\n"
"\n"
"       EXAGEOSTAT dependencies :\n")
foreach (_dep ${EXAGEOSTAT_DEP})
    set(dep_message "${dep_message}"
    "                                 ${_dep}\n")
endforeach ()
string(REGEX REPLACE ";" " " EXAGEOSTAT_DEFINITIONS_LIST "${EXAGEOSTAT_DEFINITIONS_LIST}")
set(dep_message "${dep_message}"
"\n"
"       Definitions: ${EXAGEOSTAT_DEFINITIONS_LIST}\n")
set(dep_message "${dep_message}"
"\n"
"       INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}\n\n")

string(REPLACE ";" " " dep_message_wsc "${dep_message}")
message(${dep_message})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.log "${dep_message_wsc}")
message(STATUS "Configuration is done - A summary of the current configuration"
"\n   has been written in ${CMAKE_CURRENT_BINARY_DIR}/config.log")
# installation
# ------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.log DESTINATION share/exageostat)
