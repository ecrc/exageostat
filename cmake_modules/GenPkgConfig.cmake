###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2017 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenPkgConfig.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.1
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2020-01-06
#
###

###
#
# CONVERT_LIBSTYLE_TO_PKGCONFIG: convert a libraries list to follow the pkg-config style
#                                used in CLEAN_LIB_LIST
#
###
MACRO(CONVERT_LIBSTYLE_TO_PKGCONFIG _liblist)
    set(${_liblist}_CPY "${${_liblist}}")
    set(${_liblist} "")
    foreach(_dep ${${_liblist}_CPY})
        if (${_dep} MATCHES "^/")
            get_filename_component(dep_libname ${_dep} NAME)
            get_filename_component(dep_libdir  ${_dep} DIRECTORY)
            STRING(REPLACE ".dylib" "" dep_libname "${dep_libname}")
            STRING(REPLACE "lib"    "" dep_libname "${dep_libname}")
            STRING(REPLACE ".so"    "" dep_libname "${dep_libname}")
            STRING(REPLACE ".a"     "" dep_libname "${dep_libname}")
            STRING(REPLACE ".dll"   "" dep_libname "${dep_libname}")
            list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
        elseif(NOT ${_dep} MATCHES "^-")
            list(APPEND ${_liblist} "-l${_dep}")
        else()
            list(APPEND ${_liblist} ${_dep})
        endif()
    endforeach()
ENDMACRO(CONVERT_LIBSTYLE_TO_PKGCONFIG)

###
#
# CLEAN_LIB_LIST: clean libraries lists to follow the pkg-config style
#                 used in GENERATE_PKGCONFIG_FILE
#
###
MACRO(CLEAN_LIB_LIST _package)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS_PRIVATE)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
    CONVERT_LIBSTYLE_TO_PKGCONFIG(${_package}_PKGCONFIG_LIBS)
    CONVERT_LIBSTYLE_TO_PKGCONFIG(${_package}_PKGCONFIG_LIBS_PRIVATE)
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
ENDMACRO(CLEAN_LIB_LIST)

###
#
# GENERATE_PKGCONFIG_FILE: generate files exageostat.pc, coreblas.pc and cudablas.pc
#
###
MACRO(GENERATE_PKGCONFIG_FILE)

    # The definitions that should be given to users (change the API)
    set(EXAGEOSTAT_PKGCONFIG_DEFINITIONS "")

    # The link flags specific to this package and any required libraries
    # that don't support PkgConfig
    set(EXAGEOSTAT_PKGCONFIG_LIBS "-lexageostat")

    # The link flags for private libraries required by this package but not
    # exposed to applications
    set(EXAGEOSTAT_PKGCONFIG_LIBS_PRIVATE "")

    # A list of packages required by this package
    set(EXAGEOSTAT_PKGCONFIG_REQUIRED "")

    # A list of private packages required by this package but not exposed to
    # applications
    set(EXAGEOSTAT_PKGCONFIG_REQUIRED_PRIVATE "")

    if(EXAGEOSTAT_SCHED_STARPU)
        if ( EXAGEOSTAT_USE_MPI )
            list(APPEND EXAGEOSTAT_PKGCONFIG_REQUIRED_PRIVATE libstarpumpi)
        else()
            list(APPEND EXAGEOSTAT_PKGCONFIG_REQUIRED_PRIVATE libstarpu)
        endif()
    elseif(EXAGEOSTAT_SCHED_QUARK)
        list(APPEND EXAGEOSTAT_PKGCONFIG_LIBS_PRIVATE "${QUARK_LIBRARIES_DEP}")
    endif()

    if(NOT EXAGEOSTAT_SIMULATION)

        list(APPEND EXAGEOSTAT_PKGCONFIG_LIBS_PRIVATE
        ${EXTRA_LIBRARIES}
        )
        list(APPEND EXAGEOSTAT_PKGCONFIG_REQUIRED "chameleon")

    else(NOT EXAGEOSTAT_SIMULATION)

        list(APPEND EXAGEOSTAT_PKGCONFIG_LIBS
        -lchameleon
        ${EXTRA_LIBRARIES}
        )

    endif(NOT EXAGEOSTAT_SIMULATION)

    # Define required package
    # -----------------------
    CLEAN_LIB_LIST(EXAGEOSTAT)

    # Create .pc file
    # ---------------
    SET(_output_exageostat_file "${CMAKE_BINARY_DIR}/exageostat.pc")

    # TODO: add url of MORSE releases in .pc file
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/exageostat.pc.in" "${_output_exageostat_file}" @ONLY)

    # installation
    # ------------
    INSTALL(FILES ${_output_exageostat_file} DESTINATION lib/pkgconfig)


ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##
